function GetPatientStatistics(PatPosStruct,EXCEL_FILENAME,SHEET_NUM,JUST_BASIC)
%   GetPatientStatistics
%       Writes the statistics of a patient population to an excel file.
%       
%   GetPatientStatistics(PatPosStruct,EXCEL_FILENAME,SHEET_NUM)
%   
%   PatPosStruct      The patient structure used in the study.
%   
%   EXCEL_FILENAME  The excel filename to append the data.  The default
%                   action is to append the data to the next available
%                   sheet.
%
%
%

if nargin == 3
    JUST_BASIC = true;
end


WantedFields = {'PtID','string';...
                'Study','string';...
                'PR_seqs','length';...
                'RT_seqs','length';...
                'IS_RESPONDER','data';...
                'RX_vals','RXspecial'};

            
GlobalRXData = [];
GlobalData = cell(length(PatPosStruct),size(WantedFields,1));
for i = 1:length(PatPosStruct)
    for k = 1:size(WantedFields,1)
        if isfield(PatPosStruct{i},WantedFields{k,1})
            switch WantedFields{k,2}

                case 'string'
                    GlobalData{i,k} = getfield(PatPosStruct{i},WantedFields{k,1});
                
                case 'length'
                    temp = getfield(PatPosStruct{i},WantedFields{k,1});
                    while iscell(temp)
                        temp=cell2mat(temp);
                    end
                    GlobalData{i,k} = size(temp,1)*(numel(temp)>0);
                    
                case 'data'
                    if iscell(getfield(PatPosStruct{i},WantedFields{k,1}))
                        GlobalData(i,k) = getfield(PatPosStruct{i},WantedFields{k,1});
                    else
                        GlobalData{i,k} = getfield(PatPosStruct{i},WantedFields{k,1});
                    end

                case 'RXspecial'
                    rx_vals=PatPosStruct{i}.RX_vals;
                    initial_spot=find(rx_vals(:,1)>=0&sum(rx_vals(:,2:end),2)>0,1);
                    rx_spots = find(rx_vals(initial_spot,2:end));
                    if isempty(rx_spots)
                        GlobalData{i,k}='Discarded';
                        continue
                    end
                    GlobalRXData = [GlobalRXData;rx_vals(initial_spot,2:end)];
                    rx_temp=PatPosStruct{i}.RX_names(rx_spots(1));
                    for g=2:nnz(rx_spots)
                        rx_temp = [rx_temp ', ' PatPosStruct{i}.RX_names(rx_spots(g))];
                    end
                    GlobalData{i,k}=cat(2,rx_temp{:});
            end
        end
    end
end


StudiesIncluded=unique(GlobalData(:,2));
[tf LOC]=ismember(GlobalData(:,2),StudiesIncluded);
NumPerStudy=histc(LOC,1:length(StudiesIncluded));

NumResponders = sum(cell2mat(GlobalData(:,5)));
NumNonResponders = sum(~cell2mat(GlobalData(:,5)));

NumPRSeqs = sum(cell2mat(GlobalData(:,3)));
NumRTSeqs = sum(cell2mat(GlobalData(:,4)));


TherapyRegimines = unique(GlobalData(:,6));
[tf LOC]=ismember(GlobalData(:,6),TherapyRegimines);
NumPerTherapy = histc(LOC, 1:length(TherapyRegimines));


BasicStats = cell(max(length(StudiesIncluded)+length(TherapyRegimines)+1,4),4);
try
BasicStats(1:length(StudiesIncluded),1)=StudiesIncluded;
BasicStats(1:length(StudiesIncluded),2)=num2cell(NumPerStudy);
BasicStats(length(StudiesIncluded)+1,1)={'Therapies Considered'};
BasicStats(length(StudiesIncluded)+2:end,1)=TherapyRegimines;
BasicStats(length(StudiesIncluded)+2:end,2)=num2cell(NumPerTherapy);
catch
    123
end

BasicStats{1,3}='Responders';
BasicStats{1,4}=NumResponders;
BasicStats{2,3}='NonResponders';
BasicStats{2,4}=NumNonResponders;

BasicStats{3,3}='RT Sequences';
BasicStats{3,4}=NumPRSeqs;
BasicStats{4,3}='PR Sequences';
BasicStats{4,4}=NumRTSeqs;

BasicStats(cellfun('isempty',BasicStats))={' '};

xlswrite(EXCEL_FILENAME,BasicStats,SHEET_NUM)


if ~JUST_BASIC
    xlswrite(EXCEL_FILENAME,{'PatID','Study','PR Seq','RT Seq','Responder','Therapy'},SHEET_NUM,['A' int2str(size(BasicStats,1)+2)])
    xlswrite(EXCEL_FILENAME,GlobalData,SHEET_NUM,['A' int2str(size(BasicStats,1)+3)])
end




