function OUTPUT=LoadAnnotatedGEData(BASE_FILENAME)
%   LoadAnnotatedGEData
%       Loads Annotated Gene Expression Data
%   
%       OUTPUT=LoadAnnotatedGEData(BASE_FILENAME)
%           BASE_FILENAME       The name of the dataset
%
%           OUTPUT              A structure containing the gene expression
%                               data, the annotated column names, and the
%                               Chip type.
%               OUTPUT.DataName         DataSet Name
%               OUTPUT.ProbeSet         Affymatrix Probes
%               OUTPUT.ChipType         Affymatrix Chip Type
%               OUTPUT.GEData           Gene Expression Data
%               OUTPUT.DiseaseNames     Names of each Disease States
%               OUTPUT.DiseaseStates    A vector indicating which Column
%                                       belongs to which Disease State.
%                                       0:      Normal State
%                                       1-N:    Disease State 1-N
%
%   FILE STRUCTURE
%
%       Annotation File
%               Filename:   BASE_FILENAME'_Annotation.txt'
%               Format:
%
%   D1      Disease State1
%   D2      Disease State2
%   ...
%   Dn      Disease StateN
%   DiseaseHeader1      D1
%   DiseaseHeader2      D1
%   NormalHeader1       N
%   NormalHeader2       N
%
%       Any number of disease states can be combined however only one
%       'Normal' state can be identified.  If you want to perform multiple
%       analysis with multiple definitions of 'Normal' then you must create
%       multiple datasets.
%
%
%       Gene Expression File
%               Filename:   BASE_FILENAME'_RMA_normalized' '_
%               Format:
%
%   First Line:     Header Names
%   First Row:      Probe IDS
%
%   Each Column is one Sample and Each Row is the expression of one Probe





%%%%%%%%%%FILE CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
files=dir([BASE_FILENAME '*']);
catch
    rethrow(lasterror)
    error('Cannot find files')
end

GEfid=fopen(files(2).name);
ANfid=fopen(files(1).name);
    
try
%%%%%%%%%%%READ ANNOTATION FILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    annot_data=textscan(ANfid,'%s %s','delimiter','\t','multipleDelimsAsOne',1);
  
    %%%%%FIND DISEASE STATE NAMES
    defining_rows=strncmpi(annot_data{1},'D',1);
    state_names=annot_data{2}(defining_rows);
    OUTPUT.DiseaseNames=state_names;
    
    temp=fgetl(GEfid);
    header_names=textscan(temp,'%s','delimiter','\t','multipleDelimsAsOne',1);
    
    format=['%s'];
    for i=2:length(header_names{1})
        format=[format '%n'];
    end
        
    temp=textscan(GEfid,format,'delimiter','\t','multipleDelimsAsOne',1);
    OUTPUT.GEData=cell2mat(temp(2:end));
catch
    rethrow(lasterror)
    fclose('all');
end
fclose('all');

%%%%%%%%%%%%%%%FIND THE COLUMNS WHICH CORRESPOND TO DISEASE STATES%%%%%%%%%
temp_disease_states=zeros(1,size(OUTPUT.GEData,2));
for i=1:length(temp_disease_states)
    state_spot=find(strncmp(header_names{1}{i+1},annot_data{1}(~defining_rows),length(header_names{1}{i+1})-4));
    
    if ~strcmp(annot_data{2}(sum(defining_rows)+state_spot),'N')   %leave 'Normal' as zero
        state_val=annot_data{2}(sum(defining_rows)+state_spot);
        temp_disease_states(i)=find(strcmp(state_val,annot_data{1}(defining_rows)));
    end
end
OUTPUT.DiseaseStates=temp_disease_states;
OUTPUT.ProbeSet=temp{1};
OUTPUT.DataName=BASE_FILENAME;

underscore_spots=findstr(files(2).name,'_');
period_spots=findstr(files(2).name,'.');
OUTPUT.ChipType=files(2).name(underscore_spots(end)+1:period_spots-1);
