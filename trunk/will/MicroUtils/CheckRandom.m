function Data=CheckRandom(NUM_CHECKS,ALL_DATA,varargin)

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



for j=1:NUM_CHECKS
display(j)

    display('Shuffling Data')
    for i=1:length(ALL_DATA)
        new_order=randperm(size(ALL_DATA(i).GEData,1));
        ALL_DATA(i).GEData=ALL_DATA(i).GEData(new_order,:);
    end


    display('Analysing Fold-Change')
    FOLD_DATA=GetSigGenesMulti(ALL_DATA,'fold-change','HUGeneFL',HU6800_Probes,HU6800_Chip,'U133A',U133A_Probes,U133A_Chip,'U133B',U133B_Probes,U133B_Chip,'U133P',U133P2_Probes,U133P2_Chip,'U95A2',U95Av2_Probes,U95Av2_Chip);
    FOLD_RESULTS=CheckGeneListOverlap(FOLD_DATA,'HUGeneFL',HU6800_Probes,HU6800_Chip,'U133A',U133A_Probes,U133A_Chip,'U133B',U133B_Probes,U133B_Chip,'U133P',U133P2_Probes,U133P2_Chip,'U95A2',U95Av2_Probes,U95Av2_Chip);


    display('Analysing T-Test')
    T_DATA=GetSigGenesMulti(ALL_DATA,'t-test','HUGeneFL',HU6800_Probes,HU6800_Chip,'U133A',U133A_Probes,U133A_Chip,'U133B',U133B_Probes,U133B_Chip,'U133P',U133P2_Probes,U133P2_Chip,'U95A2',U95Av2_Probes,U95Av2_Chip);
    T_RESULTS=CheckGeneListOverlap(T_DATA,'HUGeneFL',HU6800_Probes,HU6800_Chip,'U133A',U133A_Probes,U133A_Chip,'U133B',U133B_Probes,U133B_Chip,'U133P',U133P2_Probes,U133P2_Chip,'U95A2',U95Av2_Probes,U95Av2_Chip);


    display('Analysing SAM')
    SAM_DATA=GetSigGenesMulti(ALL_DATA,'SAM','HUGeneFL',HU6800_Probes,HU6800_Chip,'U133A',U133A_Probes,U133A_Chip,'U133B',U133B_Probes,U133B_Chip,'U133P',U133P2_Probes,U133P2_Chip,'U95A2',U95Av2_Probes,U95Av2_Chip);
    SAM_RESULTS=CheckGeneListOverlap(SAM_DATA,'HUGeneFL',HU6800_Probes,HU6800_Chip,'U133A',U133A_Probes,U133A_Chip,'U133B',U133B_Probes,U133B_Chip,'U133P',U133P2_Probes,U133P2_Chip,'U95A2',U95Av2_Probes,U95Av2_Chip);


    display('Collating Data')
    Data(j).Fold_P=FOLD_RESULTS.P_Overlap(:);
    Data(j).Fold_L=cellfun(@length,FOLD_RESULTS.GeneOverlap(:));
    Data(j).TTest_P=T_RESULTS.P_Overlap(:);
    Data(j).TTest_L=cellfun(@length,T_RESULTS.GeneOverlap(:));
    Data(j).SAM_P=SAM_RESULTS.P_Overlap(:);
    Data(j).SAM_L=cellfun(@length,SAM_RESULTS.GeneOverlap(:));
end