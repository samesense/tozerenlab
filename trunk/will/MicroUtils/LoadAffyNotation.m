function [varargout]=LoadAffyNotation(FILENAME)
% LoadAffyNotation
%   Loads Affymatrix-style Annotation
%
%  Probe_DATA=LoadAffyNotation(FILENAME)
%
%       FILENAME        Filename of the Affymatrix notation
%
%       Probe_DATA     An N X 4 cell-matrix:
%                           Col1:   Probe ID
%                           Col2:   Entrez-Gene ID
%                           Col3:   GO: Biological Process
%                           Col4:   GO: Cellular Compartment
%                           Col5:   GO: Molecular Function
%                           Col6:   Kegg Pathways
%
%  [Probe_DATA CHIP_DATA]=LoadAffyNotation(FILENAME)
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
%                           CHIP_DATA.CellCompArray      A boolean array in
%                           which each row coresponds to the probe and each
%                           column referes to one GO ID
%
%                           CHIP_DATA.CellCompNames     A Cell array of the
%                           names indexed in the boolean array
%
%                           CHIP_DATA.MolFunArray       A boolean array in
%                           which each row coresponds to the probe and each
%                           column referes to one GO ID
%
%                           CHIP_DATA.MolFunNames       A Cell array of the
%                           names indexed in the boolean array
%
%                           CHIP_DATA.BioProArray       A boolean array in
%                           which each row coresponds to the probe and each
%                           column referes to one GO ID
%
%                           CHIP_DATA.BioProNames       A Cell array of the
%                           names indexed in the boolean array

%%%%%%%%%%%%%%%ERROR CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check for file
try
    fid=fopen(FILENAME);
catch
    rethrow(lasterror);
end

wanted_cols={'Probe Set ID','Entrez Gene',...
    'Gene Ontology Biological Process','Gene Ontology Cellular Component',...
    'Gene Ontology Molecular Function','Pathway','Gene Symbol'};


%%%%%%%%%%%%%%%%%%%%%%%%LOAD FILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=waitbar(0,'Loading File');
try
    temp=fgetl(fid);
    header_names=textscan(temp,'%q','delimiter',',','bufsize',10000);

    format=[];
    for i=1:length(header_names{1})
        if sum(strcmpi(header_names{1}(i),wanted_cols))
            format=[format '%q'];
        else
            format=[format '%*q'];
        end
    end
    
    [junk LOC]=ismember([header_names{:}],wanted_cols);

    [junk INDS]=sort(nonzeros(LOC));
    
    temp=textscan(fid,format,'delimiter',',','bufsize',1000000);
    temp=temp(INDS);
    fclose(fid);
catch
    fclose(fid);
    rethrow(lasterror)
end

waitbar(0.1,h,'Collating Data');
%%%%%%%%%%SEPERATE PROBES WHICH HAVE MULTIPLE ENTREZ IDS%%%%%%%%%%%%%%%%%%
OUTPUT_DATA=cell(10*length(temp{1}),7);
counter=0;
for i=1:length(temp{1})
    waitbar(0.1+0.1*i/length(temp{1}),h)
    junk=findstr('AFFX',temp{1}{i});
    if isempty(junk)
        spots=findstr('///',temp{2}{i});
        bio_proc=splitGO(temp{3}{i});
        cell_comp=splitGO(temp{4}{i});
        mol_fun=splitGO(temp{5}{i});
        kegg_path=splitKegg(temp{6}{i});
        gene_symbol=splitGeneSymbol(temp{7}{i});

        if isempty(spots)
            counter=counter+1;
            OUTPUT_DATA{counter,1}=temp{1}{i};
            OUTPUT_DATA{counter,2}=str2num(temp{2}{i});
            OUTPUT_DATA{counter,3}=bio_proc;
            OUTPUT_DATA{counter,4}=cell_comp;
            OUTPUT_DATA{counter,5}=mol_fun;
            OUTPUT_DATA{counter,6}=kegg_path;
            OUTPUT_DATA{counter,7}=gene_symbol;
        else
            all_spots=[-2 spots length(temp{2}{i})];
            for j=2:length(all_spots)
                counter=counter+1;
                OUTPUT_DATA{counter,1}=temp{1}{i};
                OUTPUT_DATA{counter,2}=str2num(temp{2}{i}(all_spots(j-1)+3:all_spots(j)-1));
                OUTPUT_DATA{counter,3}=bio_proc;
                OUTPUT_DATA{counter,4}=cell_comp;
                OUTPUT_DATA{counter,5}=mol_fun;
                OUTPUT_DATA{counter,6}=kegg_path;
                OUTPUT_DATA{counter,7}=gene_symbol;
            end
        end
    end
end
OUTPUT_DATA=OUTPUT_DATA(1:counter,:);
varargout{1}=OUTPUT_DATA;

%Quit if only getting probe data, otherwise continue for chip statistics
if nargout==1
    close(h);
    return
end

temp=[];
waitbar(0.2,h,'Generating Chip Stats: Biological Process')
cellfun(@collate_IDS,OUTPUT_DATA(:,3));
temp=unique(temp);
CHIP_DATA.BioProNames=temp;
waitbar(0.3,h,'Generating Chip Stats: Biological Process')
CHIP_DATA.BioProArray=cell2mat(cellfun(@MakeArray,OUTPUT_DATA(:,3),'uniformoutput',false));

waitbar(0.4,h,'Generating Chip Stats: Cellular Component')
temp=[];
cellfun(@collate_IDS,OUTPUT_DATA(:,4));
temp=unique(temp);
CHIP_DATA.CellCompNames=temp;
waitbar(0.5,h,'Generating Chip Stats: Cellular Component')
CHIP_DATA.CellCompArray=cell2mat(cellfun(@MakeArray,OUTPUT_DATA(:,4),'uniformoutput',false));

waitbar(0.6,h,'Generating Chip Stats: Molecular Function')
temp=[];
cellfun(@collate_IDS,OUTPUT_DATA(:,5));
temp=unique(temp);
CHIP_DATA.MolFunNames=temp;
waitbar(0.7,h,'Generating Chip Stats: Molecular Function')
CHIP_DATA.MolFunArray=cell2mat(cellfun(@MakeArray,OUTPUT_DATA(:,5),'uniformoutput',false));

waitbar(0.8,h,'Generating Chip Stats: Kegg Pathways')
temp=[];
cellfun(@collate_IDS,OUTPUT_DATA(:,6));
temp=unique(temp);
CHIP_DATA.KeggNames=temp;
waitbar(0.9,h,'Generating Chip Stats: Kegg Pathways')
CHIP_DATA.KeggArray=cell2mat(cellfun(@MakeArray,OUTPUT_DATA(:,6),'uniformoutput',false));
temp=[];

close(h)
varargout{2}=CHIP_DATA;

    function output=splitGO(input)

        triple_spots=findstr('///',input);
        if ~isempty(triple_spots)
            triple_spots=[-2 triple_spots length(input)];
            output=cell(length(triple_spots)-1,1);
            for k=1:length(triple_spots)-1
                double_spots=findstr('//',input(triple_spots(k)+3:triple_spots(k+1)));
                output{k}=input(triple_spots(k)+double_spots(1)+5:triple_spots(k)+double_spots(2));
            end
        else
            output=[];
        end

    end

    function output=splitKegg(input)
        if ~isempty(input)
            triple_spots=findstr('///',input);
            if ~isempty(triple_spots)
                triple_spots=[-2 triple_spots length(input)];
                output=cell(length(triple_spots)-1,1);
                for k=1:length(triple_spots)-2
                    output{k}=input(triple_spots(k)+3:triple_spots(k+1)-1);
                end
                k=length(triple_spots)-1;   %last tag must be treated diferently (no trailing space)
                output{end}=input(triple_spots(k)+3:triple_spots(k+1));
            else
                output={input};
            end
        else
            output=[];
        end
    end

    function output=splitGeneSymbol(input)
        if strcmp(input,'---')
            output='';
        else
            triple_spots=findstr('///',input);
            if isempty(triple_spots)
                output=strtrim(input);
            else
                output=strtrim(input(1:triple_spots-1));
            end
            
            
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