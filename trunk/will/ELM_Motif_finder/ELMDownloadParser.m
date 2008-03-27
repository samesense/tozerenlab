function ELM_STRUCT=ELMDownloadParser(DIRECTORY)
%   ELMDownloadParser
%       Parses the HTML files of the ELM database
%       (http://elm.eu.org/browse.html) to create a STRUCT which contains
%       the Name and REGEXP of each ELM motif.
%
%
%   ELM_STRUCT=ELMDownloadParser
%       Parses from the default directory (\'ELM_RAW_DOWNLOAD') and creates
%       an ELM_STRUCT which can be used with ELMFinder.
%
%
%   See also: ELMFinder.

if nargin==0
    DIRECTORY='ELM_RAW_DOWNLOAD\';
end

files=dir(DIRECTORY);

reg_exprs=cell(length(files)-2,1);
for i=3:length(files)
    fid=fopen([DIRECTORY files(i).name]);
    FLAG=false;
    temp=fgetl(fid);
    while(isempty(temp)||~isnumeric(temp))
        if ~isempty(strfind(temp,'Pattern:'))
            FLAG=true;
        elseif FLAG
            start_spot=find(temp=='>',1)+1;
            stop_spot=find(temp=='/',1)-3;
            reg_exprs{i-2}=temp(start_spot:stop_spot);
            break
        end
        temp=fgetl(fid);
    end
    fclose(fid);
end

ELM_STRUCT=struct('Name',arrayfun(@(x)x.name(1:end-5),files(3:end),'uniformoutput',false),'REG_EXPR',reg_exprs);