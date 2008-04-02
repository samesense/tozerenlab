function [GeneGOOutput ColHeaders]=DavidParse(INPUT_FILENAME,varargin)
%   DavidParse
%       Parses the output provided by the DAVID gene-enrichment tool.  
%       http://david.abcc.ncifcrf.gov/
%
%   GeneGOOutput = DavidParse(INPUT_FILENAME)
%
%       INPUT_FILENAME      The path to a file downloaded from DAVID.
%
%       GeneGOOutput        A cell-array in which each row corresponds to a
%                           Gene and each column is a different rank in the
%                           GO heiarchy.
%                           {LLID, {GOTerms}}
%
%   [GeneGOOutput ColLabels] = DavidParse(INPUT_FILENAME)
%
%       ColLabels           Returns a 1xN cell-array describing which level
%                           of the GO heiarchy is in each Column.
%
%       DavidParse(INPUT_FILENAME)
%
%                           When called with no output arguements then the
%                           program will create an Excel file of the same
%                           name as INPUT_FILENAME.
%
%   Optional Parameters
%
%       EXCEL_NAME          Overloads the deafult filename.  Also will
%                           force an Excel File even if output arguements
%                           are requested.
%
%       MAKE_FIGURES        [true|FALSE]Produces a set of figures for each 
%                           GO category.
%
%       NUM_INCLUDE         The number of GO categories to include in the
%                           figures. DEFAULT = 10
%
%

XLS_FILENAME=INPUT_FILENAME(1:end-4);
MAKE_EXCEL_FLAG=false;
MAKE_FIGURE_FLAG=false;
P_VAL_CUTOFF=0.05;
NUM_INCLUDE=10;

if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})

            case {'excel_name', 'excel_filename'}
                if ischar(varargin{i+1})
                    XLS_FILENAME=varargin{i+1};
                    MAKE_EXCEL_FLAG=true;
                else
                    error('DavidParse:BAD_EXCELNAME','Arguement to EXCEL_NAME must be a char-array.')
                end
            
            case {'make_figures','make_figure'}
                if islogical(varargin{i+1})
                    MAKE_FIGURE_FLAG=varargin{i+1};
                else
                    error('DavidParse:BAD_MAKEFIGURES','Arguement to MAKE_FIGURES must be a logical.')
                end
                
            case 'num_include'
                if isnumeric(varargin{i+1})&&isscalar(varargin{i+1})
                    NUM_INCLUDE=varargin{i+1};
                else
                    error('DavidParse:BAD_NUMINCLUDE','Arguement to NUM_INCLUDE must be a numeric-scalar.')
                end
            
            otherwise
                error('DavidParse:BAD_ARGUEMENT','An unknown arguement was provided to DavidParse: %s',varargin{i})
        end
    end
end
   

[Category Term Count Percent PValue Genes List_Total Pop_Hits Pop_Total	Fold_Change]=textread(INPUT_FILENAME,'%s%s%s%s%s%s%s%s%s%s','delimiter','\t','bufsize',100000);
    
%remove colheaders and convert to numeric values where needed
InternalCat=Category(2:end);
AllCat=unique(InternalCat);

TermNames=Term(2:end);
NumGenes=str2double(Count(2:end));

ListFraction=str2double(Percent(2:end));

Pvals=str2double(PValue(2:end));

GOGeneList=cellfun(@(x)(str2num(x)),Genes(2:end),'uniformoutput',false);

AllEntrezIDs=unique(cat(2,GOGeneList{:}));

GeneGOOutput=cell(length(AllEntrezIDs),length(AllCat)+1);

GeneGOOutput(:,1)=num2cell(AllEntrezIDs);


for i=1:length(AllCat)
    InternalInds=find(strcmp(AllCat{i},InternalCat));
    
    LogicalCellArray=cellfun(@(x)(ismember(AllEntrezIDs,x)'),GOGeneList(InternalInds),'uniformoutput',false);
    LogicalArray=cat(2,LogicalCellArray{:});
    
    GeneGOOutput(:,i+1)=cellfun(@(x)(TermNames(InternalInds(x))'),num2cell(LogicalArray,2),'uniformoutput',false);
end

[ColHeaders{1:length(AllCat)+1}]=deal('Entrez IDs',AllCat{:});

if nargout==0||MAKE_EXCEL_FLAG
    %nested-cells and empty cells are not allowed in xlswrite
    Excel_Output=[ColHeaders; GeneGOOutput];
    for i=1:numel(Excel_Output)
        if iscell(Excel_Output{i})&&~isempty(Excel_Output{i})
            temp=Excel_Output{i};
            Excel_Output{i}=temp{1};
            for k=2:length(temp)
                Excel_Output{i}=[Excel_Output{i} ', ' temp{k}];
            end
        elseif isempty(Excel_Output{i})
            Excel_Output{i}=' ';
        end
        
    end
    xlswrite(XLS_FILENAME,Excel_Output);
end
    
if MAKE_FIGURE_FLAG
    for i=1:length(AllCat)
       InternalInds=find(strcmp(AllCat{i},InternalCat)&Pvals<=P_VAL_CUTOFF);
        
       [Y I]=sort(NumGenes(InternalInds),'descend');
       
       figure
       barh(1:min(NUM_INCLUDE,length(Y)),Y(1:min(NUM_INCLUDE,length(Y))));
       set(gca,'ytick',1:min(NUM_INCLUDE,length(Y)),'yticklabel',TermNames(InternalInds(I(1:min(NUM_INCLUDE,length(Y))))))
    end
end



















