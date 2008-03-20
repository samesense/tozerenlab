function GraphMultiAnal(FILENAME,TYPE,varargin)
%   GraphMultiAnal
%       Creates a visualization and summerization of the results from
%       RunMultiAnal.
%
%
%
%
%
%
%
%
%
%
%
%
%

NUM_OBS_CUTOFF=500;
SIG_NUM_CUTOFF=0.7;
MAG_CUTOFF=0.4;



switch lower(TYPE)
    case 'excel'
        [DrugNames DrugRegimines]=xlsread(FILENAME,1);
        NUM_PATS=DrugRegimines(:,end);
        DrugRegimines=DrugRegimines(:,1:end-2);
        
        [junk SIGNUM_DATA]=xlsread(FILENAME,2);
        [junk MEAN_MAG_DATA]=xlsread(FILENAME,3);
        [junk STD_MAG_DATA]=xlsread(FILENAME,4);
        [junk COUNT_DATA]=xlsread(FILENAME,5);
        
    case 'csv'
        DrugRegimines=csvread([FILENAME '_DRUGS.csv']);
        NUM_PATS=DrugRegimines(:,end);
        DrugRegimines=DrugRegimines(:,1:end-2);
        
        SIGNUM_DATA=csvread([FILENAME '_SIGNUM.csv']);
        MEAN_MAG_DATA=csvread([FILENAME '_MEANMAG.csv']);
        STD_MAG_DATA=csvread([FILENAME '_STDMAG.csv']);
        COUNT_DATA=csvread([FILENAME '_COUNT.csv']);
    
    
    otherwise
        error('GraphMultiAnal:UNKNOWN_TYPE','An unknown type was provided: %s',TYPE)
end

[NumRegimines NumDrugs]=size(DrugRegimines);
NumFeatures=size(SIGNUM_DATA,2);



if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case {'rxnames','rx_names','rx names'}
                if iscellstr(varargin{i+1})
                    if length(varargin{i+1})==NumDrug
                        RX_NAMES=varargin{i+1}(:);
                    else
                        error('GraphMultiAnal:BAD_RX_SIZE','Arguement for RX_NAMES did not match the NumDrugs found in the datafile.')
                    end
                else
                    error('GraphMultiAnal:BAD_RX_TYPE','Arguement for RX_NAMES must be a cell-string.')
                end
            case {'feature_names','feature names'}
                if iscellstr(varargin{i+1})
                    if length(varargin{i+1})==NumFeatures
                        FEATURE_NAMES=varargin{i+1}(:);
                    else
                        error('GraphMultiAnal:BAD_FEAT_SIZE','Arguement for FEATURE_NAMES did not match the NumFeatures found in the datafile.')
                    end
                else
                    error('GraphMultiAnal:BAD_FEAT_TYPE','Arguement for FEATURE_NAMES must be a cell-string.')
                end

            case {'elm names','elm_names','elm struct','elm_struct'}
                if iscell(varargin{i+1})
                    ELM_NAMES=varargin{i+1};
                elseif isstruct(varargin{i+1})
                    [ELM_NAMES{1:length(varargin{i+1})}]=deal(varargin{i+1}(:).Name);
                end
            case {'pat_example','pat example'}
                PAT_STRUCT=varargin{i+1};
                RX_NAMES=PAT_STRUCT{1}.RX_names(:);
                
                snp_spots=PAT_STRUCT{1}.SNP_SPOTS;
                %elm_simple=PAT_STRUCT{1}.ELM_simple;
                elm_annot=PAT_STRUCT{1}.ELM_annot;
                
                FEATURE_NAMES=cell(NumFeatures,1);
                FEATURE_NAMES(1:length(RX_NAMES))=RX_NAMES;
                
                counter=1+length(RX_NAMES);
                
                for k=1:length(snp_spots)
                    FEATURE_NAMES{counter}=[int2str(snp_spots(k)) ' A'];
                    FEATURE_NAMES{counter+1}=[int2str(snp_spots(k)) ' C'];
                    FEATURE_NAMES{counter+2}=[int2str(snp_spots(k)) ' G'];
                    FEATURE_NAMES{counter+3}=[int2str(snp_spots(k)) ' T'];
                    FEATURE_NAMES{counter+4}=[int2str(snp_spots(k)) ' GAP'];
                    
                    counter=counter+5;
                end
                
                FEATURE_NAMES(counter:counter+length(ELM_NAMES)-1)=ELM_NAMES;
                
                counter=counter+length(ELM_NAMES);
                
                for k=1:length(elm_annot)
                    FEATURE_NAMES{counter}=[ELM_NAMES{elm_annot(2,k)} ' ' int2str(elm_annot(1,k))];
                    counter=counter+1;
                end
                
            otherwise
                error('GraphMultiAnal:UNKNOWN_ARG','An unknown argument was provided: %s',varargin{i})
        end
    end
    
end



count_mask=COUNT_DATA<NUM_OBS_CUTOFF;
SIGNUM_DATA=SIGNUM_DATA.*count_mask;
MEAN_MAG_DATA=MEAN_MAG_DATA.*count_mask;
STD_MAG_DATA=STD_MAG_DATA.*count_mask;
COUNT_DATA=COUNT_DATA.*count_mask;



sig_cols=(sum(abs(SIGNUM_DATA)>=SIG_NUM_CUTOFF,1)>0)&(sum(abs(MEAN_MAG_DATA)>=MAG_CUTOFF,1)>0);


SIMPLE_ANAL=MEAN_MAG_DATA(:,sig_cols).*sign(SIGNUM_DATA(:,sig_cols));

xlswrite([FILENAME '_ANALRESULT'],SIMPLE_ANAL,2,'B2')

DrugRegimesNames=cell(NumRegimines,1);
for i=1:NumRegimines
    if nnz(DrugRegimines(i,:))==1
        DrugRegimesNames(i)=RX_NAMES(DrugRegimines(i,:)&true);
    else
        for k=find(DrugRegimines(i,:))
            DrugRegimesNames{i}=[DrugRegimesNames{i} ' ' RX_NAMES{k}];
        end
    end
end

xlswrite([FILENAME '_ANALRESULT'],DrugRegimesNames,2,'A2')
xlswrite([FILENAME '_ANALRESULT'],FEATURE_NAMES(sig_cols)',2,'B1')

pcolor(SIMPLE_ANAL')
colorbar
axis([1 NumRegimines 1 nnz(sig_cols)])
set(gca,'ytick',1.5:1:nnz(sig_cols),'yticklabel',FEATURE_NAMES(sig_cols))




