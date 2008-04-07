function [varargout]=DavidParse(INPUT_FILENAME,varargin)
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
%       DAVID_TABLE         The name of a file output by the DAVID
%                           Annotation table.  This will allow annotation
%                           of ALL GO labels.
%
%       FORCE_ORDER         If provided with a list of EntrezIDs this will
%                           force the output order into that provided.  If
%                           there is no annotation for an ID then the
%                           corresponding cells are filled with
%                           'UNANNOTATED'.  Any IDs which are annotated but
%                           not in the LIST are removed from the output and
%                           a warning is displayed.
%
%
%
%
%

XLS_FILENAME=INPUT_FILENAME(1:end-4);
DAVID_TABLE_FILENAME=[];
MAKE_EXCEL_FLAG=false;
MAKE_FIGURE_FLAG=false;
P_VAL_CUTOFF=0.05;
NUM_INCLUDE=10;
FORCED_ORDER=[];
MISSING_INDS=[];

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

            case 'david_table'
                if ischar(varargin{i+1})
                    DAVID_TABLE_FILENAME=varargin{i+1};
                else
                    error('DavidParse:BAD_DAVIDTABLE','Arguement to DAVID_TABLE must be a character-array.')
                end

                if ~(nargout==3||nargout==0)
                    error('DavidParse:BAD_OUTPUT','If DAVID_TABLE is provided then there must be 0 or 3 output.')
                end
                
            case 'force_order'
                if isnumeric(varargin{i+1})
                    FORCED_ORDER=varargin{i+1};
                else
                    error('DavidParse:BAD_DAVIDTABLE','Arguement to FORCE_ORDER must be a numeric array.')
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

if ~isempty(FORCED_ORDER)
    tempGeneGOOuput=cell(size(FORCED_ORDER,1),size(GeneGOOutput,2));
    tempGeneGOOuput(:,1)=num2cell(FORCED_ORDER);
    [TF LOC]=ismember(cell2mat(GeneGOOutput(:,1)),FORCED_ORDER);
    tempGeneGOOuput(nonzeros(LOC),2:end)=GeneGOOutput(TF,2:end);
    tempGeneGOOuput(~ismember(FORCED_ORDER,cell2mat(GeneGOOutput(:,1))),2:end)={'UNANNOTATED'};
    GeneGOOutput=tempGeneGOOuput;
end


[ColHeaders{1:length(AllCat)+1}]=deal('ENTREZ_GENE_ID',AllCat{:});

if nargout>0
    varargout{1}=GeneGOOutput;
    varargout{2}=ColHeaders;
end

if ~isempty(DAVID_TABLE_FILENAME)
    fid=fopen(DAVID_TABLE_FILENAME,'rt');
    header_row=fgetl(fid);
    fclose(fid);

    header_data=textscan(header_row,'%s','delimiter','\t');
    header_data=header_data{1};

    unwanted_cols={'Gene Name','Species'};

    format=[];
    for i=1:length(header_data)
        if ~any(strcmpi(header_data{i},unwanted_cols))
            format=[format '%s'];
        else
            format=[format '%*s'];
        end
    end

    fid=fopen(DAVID_TABLE_FILENAME,'rt');
    TableData=textscan(fid,format,'delimiter','\t');
    fclose(fid);
    TableData=cat(2,TableData{:});

    TableData(2:end,1)=cellfun(@(x)(str2num(x)),TableData(2:end,1),'uniformoutput',false);

    [TF TableInds]=ismember(cell2mat(GeneGOOutput(:,1)),cell2mat(TableData(2:end,1)));

    [junk HeaderInds]=ismember(TableData(1,:),ColHeaders);

    TableData=cellfun(@CreateNested,TableData(2:end,:),'uniformoutput',false);

    OutputTableData=cell(size(GeneGOOutput,1),length(ColHeaders));

    OutputTableData(TF,:)=TableData(nonzeros(TableInds),nonzeros(HeaderInds));
    
    if nargout==3
        varargout{3}=OutputTableData;
    end

end


if nargout==0||MAKE_EXCEL_FLAG
    OrigTableHeaders={'Category','Term','Num Genes','p-value'};
    OrigTableOutput=cell(size(InternalCat,1),4);

    OrigTableOutput(:,1)=InternalCat;
    OrigTableOutput(:,2)=TermNames;
    OrigTableOutput(:,3)=num2cell(NumGenes);
    OrigTableOutput(:,4)=num2cell(Pvals);

    OrigTableOutput=sortrows(OrigTableOutput,[1 2]);

    xlswrite(XLS_FILENAME,[OrigTableHeaders;OrigTableOutput],1);


    %nested-cells and empty cells are not allowed in xlswrite
    Excel_Output=CleanForExcel([ColHeaders; GeneGOOutput]);

    xlswrite(XLS_FILENAME,Excel_Output,2);

    if ~isempty(DAVID_TABLE_FILENAME)
        Excel_Output=CleanForExcel([ColHeaders; OutputTableData]);
        
        xlswrite(XLS_FILENAME,Excel_Output,3);

    end
end

if MAKE_FIGURE_FLAG
    for i=1:length(AllCat)
        InternalInds=find(strcmp(AllCat{i},InternalCat)&Pvals<=P_VAL_CUTOFF);
        if ~isempty(InternalInds)
            [Y I]=sort(NumGenes(InternalInds),'descend');

            figure
            barh(1:min(NUM_INCLUDE,length(Y)),Y(1:min(NUM_INCLUDE,length(Y))));
            set(gca,'ytick',1:min(NUM_INCLUDE,length(Y)),'yticklabel',TermNames(InternalInds(I(1:min(NUM_INCLUDE,length(Y))))))
            title(AllCat{i})
            xlabel('Number of Genes')
        end
    end
end


    function Output=CreateNested(Input)
        if ischar(Input)&&~isempty(Input)
            Output=textscan(Input,'%s','delimiter',',');
            Output=Output{1};
            Output=strtrim(Output);
        else
            Output=Input;
        end
    end
end










