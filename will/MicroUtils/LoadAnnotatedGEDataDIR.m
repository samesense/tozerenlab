function OUTPUT=LoadAnnotatedGEDataDIR(DIRECTORY_NAME)
%   LoadAnnotatedGEDataDIR
%       Loads All Annotated DataSets in a Directory
%
%       Iteratively loads Datasets using the LoadAnnotatedGEData function
%       into a vector of structs.
%
%   OUTPUT=LoadAnnotatedGEDataDIR
%       Loads all data in current directory.
%
%   OUTPUT=LoadAnnotatedGEDataDIR(DIRECTORY_NAME)
%       Loads all data in DIRECTORY_NAME

if nargin==0
    DIRECTORY_NAME=[];
end

files=struct2cell(dir([DIRECTORY_NAME '*.txt']))';

names=cell(length(files),1);
for i=1:length(names)
    spot=findstr(files{i,1},'_');
    names{i}=files{i,1}(1:spot(1)-1);
end

temp=unique(names);
for i=1:length(temp)
    display(temp{i});
    OUTPUT(i)=LoadAnnotatedGEData(temp{i});
end