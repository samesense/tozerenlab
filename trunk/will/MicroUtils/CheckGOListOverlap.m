function [OVERLAP_GO P_OVERLAP_GO DiseaseNameVector]=CheckGOListOverlap(GENELISTVECTOR,varargin)
%   CheckGOListOverlap
%       Checks for Overlap of GO Categories from a GeneListVector created
%       by GetSigGenesMulti.
%   [OVERLAP_GENES P_OVERLAP_GENES
%       DiseaseNameVector]=CheckGeneListOverlap(GENELISTVECTOR)
%
%   GENELISTVECTOR              A vector of genelists created by
%                               GetSigGenesMulti
%
%   OVERLAP_GO                  A NxN matrix of all overlapping genes
%   P_OVERLAP_GO                The p-value of the gene-list overlap
%                               (hypergeometric distribution); a value of
%                               -1 implies no overlap
%   DiseaseNameVector           A vector of the disease names as they are
%                               indexed in the OVERLAP_GENES matrix
%
%
%   CheckGOListOverlap(...,'Chip_set1_name',Chip_set1,...)
%       ProbeSetNames
%           'HUGeneFL'
%           'U133A'
%           'U133B'
%           'U133P'
%           'U95A2'
%

if isempty(varargin)
    load Web_Notation_Data HU6800_Chip U133A_Chip U133B_Chip U133P2_Chip U95Av2_Chip
else
    loaded=false(5,1);
    annot_sets={'HUGeneFL','U133A','U133B','U133P','U95A2';...
        'HU6800_Chip' 'U133A_Chip' 'U133B_Chip' 'U133P2_Chip' 'U95Av2_Chip'};
    for i=1:2:length(varargin)
        switch(varargin{i})
            case 'HUGeneFL'
                HU6800_Chip=varargin{i+1};
                display('HUGeneFL Provided')
                loaded(1)=1;
            case 'U133A'
                U133A_Chip=varargin{i+1};
                display('U133A Provided')
                loaded(2)=1;
            case 'U133B'
                U133B_Chip=varargin{i+1};
                display('U133B Provided')
                loaded(3)=1;
            case 'U133P'
                U133P2_Chip=varargin{i+1};
                display('U133P Provided')
                loaded(4)=1;
            case 'U95A2'
                U95Av2_Chip=varargin{i+1};
                display('U95A2 Provided')
                loaded(5)=1;
            otherwise
                warning('Unknown Chip Type Provided')
        end
    end

    if sum(loaded)~=5
        load('Web_Notation_Data.mat',annot_sets{2,~loaded})
    end
end

DiseaseStateVector=[];
DiseaseNameVector=[];
for i=1:length(GENELISTVECTOR)
    for j=1:length(GENELISTVECTOR(i).Entrez)
        DiseaseStateVector=[DiseaseStateVector; GENELISTVECTOR(i).GONames(j) GENELISTVECTOR(i).ChipType];
        DiseaseNameVector=[DiseaseNameVector; GENELISTVECTOR(i).DiseaseNames(j)];
    end
end








































