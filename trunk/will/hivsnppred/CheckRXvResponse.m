
[RX_MASK INDS]=BreakoutPatients(NewPatPosStruct,200);


PredictiveFeatures = cell(1,4);
PredictiveFeatures{1}={'BaseCalls'};
PredictiveFeatures{2}={'BaseCalls','SimpleELM'};
PredictiveFeatures{3}={'BaseCalls','SimpleELM','PositionalELM'};
PredictiveFeatures{4}={'BaseCalls','SimpleELM','PositionalELM','PositionalPWM'};


% All_predvals=cell(length(INDS),4);
% All_prednorm=cell(length(INDS),4);
% All_predstd=cell(length(INDS),4);
% All_predcount=cell(length(INDS),4);
% All_ROCs=cell(length(INDS),4);

i=9
k=3

PredictiveFeatures{k}
        [All_predvals{i,k} All_prednorm{i,k} All_predstd{i,k} All_predcount{i,k} All_ROCs{i,k}]=PredictPatients(NewPatPosStruct(INDS{i}),100,'ELM_STRUCT',NewELMPosStruct,'display',false,'predictivefeatures',PredictiveFeatures{k});
        save NewRunData All_predvals All_prednorm All_predstd All_predcount All_ROCs
        
i=9
k=4

PredictiveFeatures{k}
        [All_predvals{i,k} All_prednorm{i,k} All_predstd{i,k} All_predcount{i,k} All_ROCs{i,k}]=PredictPatients(NewPatPosStruct(INDS{i}),100,'ELM_STRUCT',NewELMPosStruct,'display',false,'predictivefeatures',PredictiveFeatures{k});
        save NewRunData All_predvals All_prednorm All_predstd All_predcount All_ROCs

        

for i = 10:length(INDS)
    i
    for k = 1:size(PredictiveFeatures,2)
        PredictiveFeatures{k}
        [All_predvals{i,k} All_prednorm{i,k} All_predstd{i,k} All_predcount{i,k} All_ROCs{i,k}]=PredictPatients(NewPatPosStruct(INDS{i}),100,'ELM_STRUCT',NewELMPosStruct,'display',false,'predictivefeatures',PredictiveFeatures{k});
        save NewRunData All_predvals All_prednorm All_predstd All_predcount All_ROCs
    end
    
end
% 
% PredCount = cat(1,All_predcount{:});
% PredNorm = cat(1,All_prednorm{:});
% PredSTD = cat(1,All_predstd{:});
% PredVals = cat(1,All_predvals{:});
% 
% [FEATURES RESPONSE PAT_INDEX FEATURE_NAMES FEATURE_INDS]=GetPatientFeatures(PatPosStruct,ELM_POS_STRUCT,'DrugRegimine','BaseCalls','SimpleELM','PositionalELM','PositionalPWM');
% 
% 
% 
% InterestingRows = find(any(PredCount>20&abs(PredVals)>.9&PredNorm>0.05,1));
% for i = 1:length(FEATURE_NAMES)
%     FEATURE_NAMES{i}(FEATURE_NAMES{i}=='_') = ' ';
% end
% 
% 
% % xlswrite('PrelimOutput.xlsx',FEATURE_NAMES(InterestingRows))
% % xlswrite('PrelimOutput.xlsx',PredNorm(:,InterestingRows).*sign(PredVals(:,InterestingRows)),1,'A2')
% 
% 
% normFeatures = bsxfun(@rdivide,FEATURES,nanmax(FEATURES));
% 
% 
% spots = zeros(length(InterestingRows),2);
% names = FEATURE_NAMES(InterestingRows);
% 
% for i = 1:length(names)
%     bounder_spots=find(names{i}=='['|names{i}==']');
%     if ~isempty(bounder_spots)
%         spots(i,:) = str2num(names{i}(bounder_spots(1):bounder_spots(2)));
%     end
% end
% 
% goodSpots = spots(:,1)~=0;
% 
% myColor = usercolormap([1 1 1],[0 0 1]);
% 
% HIVFigGen(HIV_STRUCT,spots(goodSpots,:)+rand(nnz(goodSpots),2),nanmean(normFeatures(:,InterestingRows(goodSpots))),names(goodSpots),'stacked_windows','use_aa',true,'colormap',myColor);