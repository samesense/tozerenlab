function [OUTPUT_DATA]=CheckGeneListOverlap(GENELISTVECTOR,varargin)
%   CheckGeneListOverlap
%       Checks for Overlap of Entrez Gene IDS from a GeneListVector created
%       by GetSigGenesMulti.
%   [OUTPUT_DATA]=CheckGeneListOverlap(GENELISTVECTOR)
%
%   GENELISTVECTOR              A vector of genelists created by
%                               GetSigGenesMulti
%
%   OUTPUT_DATA
%       OUTPUT_DATA.GeneOverlap         A cell array of the overlaping
%                                       Entrez Gene IDS
%
%       OUTPUT_DATA.KeggOverlap         A cell array of the overlaping Kegg
%                                       pathways
%
%       OUTPUT_DATA.GOOverlap           A cell array of the overlaping Gene
%                                       Ontology Categories
%
%       OUTPUT_DATA.P_Overlap           A matrix of the p-values associated
%                                       with the Gene List overlap matrix
%
%       OUTPUT_DATA.DiseaseNameVector   A vector of the disease names as
%                                       they are indexed in the Overlap
%                                       matrix
%
%
%   CheckGeneListOverlap(...,'probe_set1_name',Probe_set1,chip_set1,...)
%       ProbeSetNames
%           'HUGeneFL'
%           'U133A'
%           'U133B'
%           'U133P'
%           'U95A2'
%

if isempty(varargin)
    load Web_Notation_Data
else
    loaded=false(10,1);
    annot_sets={'HU6800_Probes','HU6800_Chip', 'U133A_Probes','U133A_Chip',...
        'U133B_Probes','U133B_Chip', 'U133P2_Probes','U133P2_Chip',...
        'U95Av2_Probes','U95Av2_Chip'};
    for i=1:3:length(varargin)
        switch(varargin{i})
            case 'HUGeneFL'
                if isstruct(varargin{i+1})
                    HU6800_Chip=varargin{i+1};
                    HU6800_Probes=varargin{i+2};
                else
                    HU6800_Chip=varargin{i+2};
                    HU6800_Probes=varargin{i+1};
                end
                display('HUGeneFL Provided')
                loaded(1:2)=[1 1];
            case 'U133A'
                if isstruct(varargin{i+1})
                    U133A_Chip=varargin{i+1};
                    U133A_Probes=varargin{i+2};
                else
                    U133A_Chip=varargin{i+2};
                    U133A_Probes=varargin{i+1};
                end
                display('U133A Provided')
                loaded(3:4)=[1 1];
            case 'U133B'
                if isstruct(varargin{i+1})
                    U133B_Chip=varargin{i+1};
                    U133B_Probes=varargin{i+2};
                else
                    U133B_Chip=varargin{i+2};
                    U133B_Probes=varargin{i+1};
                end
                display('U133B Provided')
                loaded(5:6)=[1 1];
            case 'U133P'
                if isstruct(varargin{i+1})
                    U133P2_Chip=varargin{i+1};
                    U133P2_Probes=varargin{i+2};
                else
                    U133P2_Chip=varargin{i+2};
                    U133P2_Probes=varargin{i+1};
                end
                display('U133P Provided')
                loaded(7:8)=[1 1];
            case 'U95A2'
                if isstruct(varargin{i+1})
                    U95Av2_Chip=varargin{i+1};
                    U95Av2_Probes=varargin{i+2};
                else
                    U95Av2_Chip=varargin{i+2};
                    U95Av2_Probes=varargin{i+1};
                end
                display('U95A2 Provided')
                loaded(9:10)=[1 1];
            otherwise
                warning('Unknown Chip Type Provided')
        end
    end
    if sum(loaded)~=length(loaded)
        load('Web_Notation_data.mat',annot_sets{~loaded})
    end
end


chipnames=unique({GENELISTVECTOR.ChipType});
chip_overlap_genes=cell(length(chipnames));
%%%%%%%%%%%%%%%%MARK ALL OVERLAPING GENES
for i=1:length(chipnames)
    for j=1:length(chipnames)
        switch chipnames{i}
            case 'HUGeneFL'
                invalid_probes=cellfun(@isempty,HU6800_Probes(:,2));
                genes1=cell2mat(HU6800_Probes(~invalid_probes,2));
            case 'U133A'
                invalid_probes=cellfun(@isempty,U133A_Probes(:,2));
                genes1=cell2mat(U133A_Probes(~invalid_probes,2));
            case 'U133B'
                invalid_probes=cellfun(@isempty,U133B_Probes(:,2));
                genes1=cell2mat(U133B_Probes(~invalid_probes,2));
            case 'U133P'
                invalid_probes=cellfun(@isempty,U133P2_Probes(:,2));
                genes1=cell2mat(U133P2_Probes(~invalid_probes,2));
            case {'U95A','U95A2'}
                invalid_probes=cellfun(@isempty,U95Av2_Probes(:,2));
                genes1=cell2mat(U95Av2_Probes(~invalid_probes,2));
        end
        switch chipnames{j}
            case 'HUGeneFL'
                invalid_probes=cellfun(@isempty,HU6800_Probes(:,2));
                genes2=cell2mat(HU6800_Probes(~invalid_probes,2));
            case 'U133A'
                invalid_probes=cellfun(@isempty,U133A_Probes(:,2));
                genes2=cell2mat(U133A_Probes(~invalid_probes,2));
            case 'U133B'
                invalid_probes=cellfun(@isempty,U133B_Probes(:,2));
                genes2=cell2mat(U133B_Probes(~invalid_probes,2));
            case 'U133P'
                invalid_probes=cellfun(@isempty,U133P2_Probes(:,2));
                genes2=cell2mat(U133P2_Probes(~invalid_probes,2));
            case {'U95A','U95A2'}
                invalid_probes=cellfun(@isempty,U95Av2_Probes(:,2));
                genes2=cell2mat(U95Av2_Probes(~invalid_probes,2));
        end
        chip_overlap_genes{i,j} =intersect(genes1,genes2);
    end
end

DiseaseStateVector=[];
DiseaseNameVector=[];
for i=1:length(GENELISTVECTOR)
    for j=1:length(GENELISTVECTOR(i).Entrez)
        DiseaseStateVector=[DiseaseStateVector; GENELISTVECTOR(i).Entrez(j) GENELISTVECTOR(i).ChipType GENELISTVECTOR(i).GONames(j) GENELISTVECTOR(i).KeggPath(j)];
        DiseaseNameVector=[DiseaseNameVector; GENELISTVECTOR(i).DiseaseNames(j)];
    end
end

%%%%%%%%%%%%%%%%%FIND OVERLAP OF GENELISTS
OVERLAP_GENES=cell(length(DiseaseStateVector));
OVERLAP_KEGG=cell(length(DiseaseStateVector));
OVERLAP_GO=cell(length(DiseaseStateVector));
P_OVERLAP_GENES=ones(length(DiseaseStateVector));
for i=1:length(DiseaseStateVector)
    for j=1:length(DiseaseStateVector)
        OVERLAP_GENES{i,j}=intersect(cell2mat(DiseaseStateVector{i,1}),cell2mat(DiseaseStateVector{j,1}));
        OVERLAP_KEGG{i,j}=intersect(DiseaseStateVector{i,4},DiseaseStateVector{j,4});
        OVERLAP_GO{i,j}=intersect(DiseaseStateVector{i,3},DiseaseStateVector{j,3});
        if length(OVERLAP_GENES{i,j})~=0
            chip1=find(strcmp(DiseaseStateVector{i,2},chipnames));
            chip2=find(strcmp(DiseaseStateVector{j,2},chipnames));
            switch DiseaseStateVector{i,2}
                case 'HUGeneFL'
                    genes1=length(HU6800_Probes(:,2));
                case 'U133A'
                    genes1=length(U133A_Probes(:,2));
                case 'U133B'
                    genes1=length(U133B_Probes(:,2));
                case 'U133P'
                    genes1=length(U133P2_Probes(:,2));
                case {'U95A','U95A2'}
                    genes1=length(U95Av2_Probes(:,2));
            end

            defective_picked=length(OVERLAP_GENES{i,j});
            total_lot_size=genes1;
            num_defective=length(chip_overlap_genes{chip1,chip2});
            num_picked=length(DiseaseStateVector{i,1});
            P_OVERLAP_GENES(i,j)=hygecdf(defective_picked,total_lot_size,num_defective,num_picked);
            
            defective_picked=length(OVERLAP_GENES{i,j});
            total_lot_size=genes1;
            num_defective=length(chip_overlap_genes{chip1,chip2});
            num_picked=length(DiseaseStateVector{i,1});
            P_OVERLAP_GENES(i,j)=hygecdf(defective_picked,total_lot_size,num_defective,num_picked);
        else
            P_OVERLAP_GENES(i,j)=-1;
        end
    end
end

OUTPUT_DATA.GeneOverlap=OVERLAP_GENES;
OUTPUT_DATA.KeggOverlap=OVERLAP_KEGG;
OUTPUT_DATA.GOOverlap=OVERLAP_GO;
OUTPUT_DATA.P_Overlap=P_OVERLAP_GENES;
OUTPUT_DATA.DiseaseNameVector=DiseaseNameVector;
