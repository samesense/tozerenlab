function OUTPUT_DATA=GetSigGenesMulti(INPUT_DATA,ANAL_TYPE,varargin)
%   GetSigGenesMulti
%       Runs GetSigGenes on a vector of values.

h=waitbar(0,'Loading Notation Data');
if isempty(varargin)
    load New_Notation_Data
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
        load('New_Notation_data.mat',annot_sets{~loaded})
    end
end


num_diseases=sum(arrayfun(@(x) length(x.DiseaseNames),INPUT_DATA));
bar_start=0.1;
for i=1:length(INPUT_DATA)
    bar_stop=bar_start+0.9*length(INPUT_DATA(i).DiseaseNames)/num_diseases;
    switch deblank(INPUT_DATA(i).ChipType)
        case 'U133A'
            OUTPUT_DATA(i)=GetSigGenes(INPUT_DATA(i),U133A_Probes,U133A_Chip,ANAL_TYPE,'waitbar',h,[bar_start bar_stop]);
        case 'U133B'
            OUTPUT_DATA(i)=GetSigGenes(INPUT_DATA(i),U133B_Probes,U133B_Chip,ANAL_TYPE,'waitbar',h,[bar_start bar_stop]);
        case 'HUGeneFL'
            OUTPUT_DATA(i)=GetSigGenes(INPUT_DATA(i),HU6800_Probes,HU6800_Chip,ANAL_TYPE,'waitbar',h,[bar_start bar_stop]);
        case 'U133P'
            OUTPUT_DATA(i)=GetSigGenes(INPUT_DATA(i),U133P2_Probes,U133P2_Chip,ANAL_TYPE,'waitbar',h,[bar_start bar_stop]);
        case {'U95A' 'U95A2'}
            OUTPUT_DATA(i)=GetSigGenes(INPUT_DATA(i),U95Av2_Probes,U95Av2_Chip,ANAL_TYPE,'waitbar',h,[bar_start bar_stop]);
        otherwise
            INPUT_DATA.ChipType
    end
    bar_start=bar_stop;
end
close(h)