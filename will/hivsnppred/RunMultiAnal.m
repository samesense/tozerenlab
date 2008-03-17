function RunMultiAnal(PAT_STRUCT,FILENAME,varargin)
%   RunMultiAnal
%       Runs the PredictPatients function on all therapy regimines and
%       records the output into an Excel file.  If the Excel file already
%       exists then it continues from the last iteration recorded.
%
%
%   RunMultiAnal(PAT_STRUCT,EXCEL_FILENAME)   
%
%   PAT_STRUCT      A patient structure as created by PatientStructHelper
%
%   EXCEL_FILENAME  An excel file to save the data to.
%
%                   SHEET1 will have the regimine data
%                   SHEET2 will have the Signum data (OUTPUT1)
%                   SHEET3 will have the scaled CORR_DATA (OUTPUT2)
%                   SHEET4 will have the STD of CORR_DATA (OUTPUT3)
%
%
%   Optional Parameters
%
%   NumPredictReps      The number of independant training/testing
%                       validations to perform in PredictPatients. 
%                       DEFAULT = 100.
%
%   PredictPatientOpt   A cell of options to pass to PredictPatients to
%                       override defaults. DEFAULT = {}
%
%   NumPatientCutoff    The number of patients required for a therapy
%                       regemine. DEFAULT = 200
%
%   CSVWrite            Toggles whether to use csvwrite instead of
%                       xlswrite.  This is important if you don't have a
%                       pre-EXCEL2007 as these are limited in the number of
%                       possible column. DEFAULT = false
%
%
%   SUPER               Option to toggle the use of the BLADECENTER for
%                       parellelization. DEFAULT = false
%
%
%
% 
%

EXCEL_FLAG=true;
SUPER_FLAG=false;
NUM_PRED=100;
NUM_PAT=200;
PRED_PAT_OPT={};

RX_names=PAT_STRUCT{1}.RX_names';


if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case 'numpredictreps'
                if isnumeric(varargin{i+1})&&isscalar(varargin{i+1})&&varargin{i+1}>0
                    NUM_PRED=varargin{i+1};
                else
                    error('RunMultiAnal:BAD_NUMPRED','Arguement for NumPredictReps must be an scalar.')
                end
                
            case 'predictpatientopt'
                if iscell(varargin{i+1})
                    PRED_PAT_OPT=varargin{i+1};
                else
                    error('RunMultiAnal:BAD_PATOPT','Arguement for PredictPatientOpt must be a cell-array.')
                end
                
            case 'numpatientcutoff'
                if isnumeric(varargin{i+1})&&isscalar(varargin{i+1})&&varargin{i+1}>0
                    NUM_PAT=varargin{i+1};
                else
                    error('RunMultiAnal:BAD_NUMPAT','Arguement for NumPatientCutoff must be an scalar.')
                end
                
            case 'csvwrite'
                if islogical(varargin{i+1})
                    EXCEL_FLAG=~varargin{i+1};
                else
                    error('RunMultiAnal:BAD_CSVWRITE','Arguement for CSVWrite must be an logical.')
                end
                
            case 'super'
                if islogical(varargin{i+1})
                    SUPER_FLAG=varargin{i+1};
                else
                    error('RunMultiAnal:BAD_SUPER','Arguement for SUPER must be an logical.')
                end
                
                
            otherwise
                error('RunMultiAnal:BAD_ARGUEMENT','An unknown arguement was provided to RunMultiAnal: %s', varargin{i})
                
        end
        
    end
    
end

[DrugRegimine DR_Inds]=BreakoutPatients(PAT_STRUCT,NUM_PAT);


if any(FILENAME=='.')
    FILENAME=FILENAME(1:find(FILENAME=='.',1)-1);
end

if SUPER_FLAG

    [NUM_PRED_CELL{1,1:length(DR_Inds)}]=deal(NUM_PRED);
    
    PAT_CELL=cell(1,length(DR_Inds));
    
    for i=1:length(DR_Inds)
        PAT_CELL{i}=PAT_STRUCT(DR_Inds{i});
    end
    
    [data1 data2 data3 data4] = dfeval(@SUPER_PredictPatients,...
        PAT_CELL,NUM_PRED_CELL,...
        'jobmanager','WDjobs','lookupURL','biomanageC',...
        'FileDependencies',{'SUPER_PredictPatients','CalculateROC','PatientStructHelper'});
    
    mask=cellfun('isempty',data1);
    
    numcols=size(data1{find(~mask,1)},2);
    
    [data1{mask}]=deal(NaN(1,numcols));
    [data2{mask}]=deal(NaN(1,numcols));
    [data3{mask}]=deal(NaN(1,numcols));
    [data4{mask}]=deal(NaN(1,numcols));
    
    csvwrite([FILENAME '_DRUGS.csv'],[DrugRegimine cellfun('length',DR_Inds)])
    csvwrite([FILENAME '_SIGNUM.csv'],cell2mat(data1));
    csvwrite([FILENAME '_MEANMAG.csv'],cell2mat(data2));
    csvwrite([FILENAME '_STDMAG.csv'],cell2mat(data3));
    csvwrite([FILENAME '_COUNT.csv'],cell2mat(data4));
    
    return
    
end



try
    if EXCEL_FLAG
        [Numeric_2 Text_2 Raw_2]=xlsread(FILENAME,1);
    else
        Numeric_2=csvread([FILENAME '_DRUG.csv']);
    end
    
    LastIND=size(Numeric_2,1);
    
catch
    
    if EXCEL_FLAG
        xlswrite(FILENAME,RX_names,1);
        xlswrite(FILENAME,RX_names,2);
        xlswrite(FILENAME,RX_names,3);
        xlswrite(FILENAME,RX_names,4);
    else
        csvwrite([FILENAME '_DRUG.csv'],RX_names,1);
        csvwrite([FILENAME '_PREDVAL.csv'],0);
        csvwrite([FILENAME '_PREDNORM.csv'],0);
        csvwrite([FILENAME '_PREDSTD.csv'],0);
    end


    LastIND=1;

end




for I=LastIND:length(DR_Inds)

    tic

    [NewPreds NewMeans NewSTDs]=PredictPatients(PAT_STRUCT(DR_Inds{I}),NUM_PRED,PRED_PAT_OPT{:});
    if isempty(NewPreds)
        continue
    end

    if EXCEL_FLAG
        xlswrite(FILENAME,DrugRegimine(I,:),1,['A' num2str(I+1)]);
        xlswrite(FILENAME,NewPreds,2,['A' num2str(I+1)]);
        xlswrite(FILENAME,NewMeans,3,['A' num2str(I+1)]);
        xlswrite(FILENAME,NewSTDs,4,['A' num2str(I+1)]);
    else
        csvwrite([FILENAME '_DRUG.csv'],DrugRegimine(I,:),I,0);
        csvwrite([FILENAME '_PREDVAL.csv'],NewPreds,I,0);
        csvwrite([FILENAME '_PREDNORM.csv'],NewMeans,I,0);
        csvwrite([FILENAME '_PREDSTD.csv'],NewSTDs,I,0);
    end

    time=toc/60;
    
    fprintf('Finished Running Index %d of %d in %s min.\n',I,length(DR_Inds),int2str(time));


end


