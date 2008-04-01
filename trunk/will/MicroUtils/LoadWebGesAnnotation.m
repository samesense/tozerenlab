function [varargout]=LoadWebGesAnnotation(FILENAME)
% LoadWebGesAnnotation
%   Loads WebGestault Style Annotation
%
%  Probe_DATA=LoadWebGesAnnotation(FILENAME)
%
%       FILENAME        Filename of the Affymatrix notation
%
%       Probe_DATA     An N X 4 cell-matrix:
%                           Col1:   Probe ID
%                           Col2:   Entrez-Gene ID
%                           Col3:   GO Terms
%                           Col4:   Kegg Pathways
%                           Col5:   Gene Name
%
%  [Probe_DATA CHIP_DATA]=LoadWebGesAnnotation(FILENAME)
%
%       CHIP_DATA       A struct containing information for statistical
%                       analysis on the chip.
%
%                           CHIP_DATA.KeggArray  A boolean array in
%                           which each row coresponds to the probe and each
%                           column refers to one Kegg Pathway
%
%                           CHIP_DATA.KeggNames  A Cell array of the
%                           names indexed in the boolean array
%
%                           CHIP_DATA.GOArray      A boolean array in
%                           which each row coresponds to the probe and each
%                           column referes to one GO ID
%
%                           CHIP_DATA.GONames     A Cell array of the
%                           names indexed in the boolean array


%%%%%%%%%%%%%%%ERROR CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check for file
try
    fid=fopen(FILENAME);
catch
    rethrow(lasterror);
end

wanted_cols={'input_id','locus link id',...
    'gene_name','go_term','kegg_title'};


%%%%%%%%%%%%%%%%%%%%%%%%LOAD FILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=waitbar(0,'Loading File');
try
    temp=fgetl(fid);
    header_names=textscan(temp,'%q','delimiter',',','bufsize',10000);
    junk=fgetl(fid);
    format=[];
    for i=1:length(header_names{1})
        if sum(strcmp(header_names{1}(i),wanted_cols))
            format=[format '%q'];
        else
            format=[format '%*q'];
        end
    end

    temp=textscan(fid,format,'delimiter',',','bufsize',1000000);
    fclose(fid);
catch
    fclose(fid);
    rethrow(lasterror)
end

waitbar(0.1,h,'Collating Data');

OUTPUT_DATA=cell(length(temp{1}),5);
counter=0;
for i=1:length(temp{1})
    waitbar(0.1+0.1*i/length(temp{1}),h)
    
    OUTPUT_DATA{i,1}=temp{1}{i};
    OUTPUT_DATA{i,2}=str2num(temp{2}{i});
    OUTPUT_DATA{i,3}=splitTerms(temp{4}{i});
    OUTPUT_DATA{i,4}=splitTerms(temp{5}{i});
    OUTPUT_DATA{i,5}=temp{3}{i};
end
varargout{1}=OUTPUT_DATA;

%Quit if only getting probe data, otherwise continue for chip statistics
if nargout==1
    return
end

temp=[];
waitbar(0.2,h,'Generating Chip Stats: GO Terms')
cellfun(@collate_IDS,OUTPUT_DATA(:,3));
temp=unique(temp);
CHIP_DATA.GONames=temp;
waitbar(0.5,h,'Generating Chip Stats: GO Terms')
CHIP_DATA.GOArray=cell2mat(cellfun(@MakeArray,OUTPUT_DATA(:,3),'uniformoutput',false));

waitbar(0.7,h,'Generating Chip Stats: KEGG Pathways')
temp=[];
cellfun(@collate_IDS,OUTPUT_DATA(:,4));
temp=unique(temp);
CHIP_DATA.KeggNames=temp';
waitbar(0.9,h,'Generating Chip Stats: KEGG Pathways')
CHIP_DATA.KeggArray=cell2mat(cellfun(@MakeArray,OUTPUT_DATA(:,4),'uniformoutput',false));

close(h)
varargout{2}=CHIP_DATA;







    function output=splitTerms(input)
        if ~isempty(input)
            t_temp=textscan(input,'%s','delimiter','///');
            bad_cells=cellfun(@isempty,t_temp{1});
            output=t_temp{1}(~bad_cells);
            
        else
            output=[];
        end
    end

    function collate_IDS(input)
        temp=union(temp,input);
    end

    function output=MakeArray(input)
        output=false(1,length(temp));
        for p=1:length(input)
            output(strcmp(input(p),temp))=1;
        end
    end
end
































