
[RX_MASK INDS]=BreakoutPatients(NewPatPosStruct,200);


PredictiveFeatures = cell(1,7);
PredictiveFeatures{1}={'SimpleHumanMiRNA'};
PredictiveFeatures{2}={'SimpleHumanMiRNA','BaseCalls'};
PredictiveFeatures{3}={'BaseCalls'};
PredictiveFeatures{4}={'BaseCalls','SimpleELM'};
PredictiveFeatures{5}={'BaseCalls','SimpleELM','PositionalELM'};
PredictiveFeatures{6}={'BaseCalls','SimpleELM','PositionalELM','PositionalPWM'};
PredictiveFeatures{7}={'BaseCalls','SimpleELM','PositionalELM','PositionalPWM','SimpleHumanMiRNA'};


my_predvals=cell(length(INDS),7);
my_prednorm=cell(length(INDS),7);
my_predstd=cell(length(INDS),7);
my_predcount=cell(length(INDS),7);
my_ROCs=cell(length(INDS),7);

SD_predvals=cell(length(INDS),7);
SD_prednorm=cell(length(INDS),7);
SD_predstd=cell(length(INDS),7);
SD_predcount=cell(length(INDS),7);
SD_ROCs=cell(length(INDS),7);



for i = 1:length(INDS)
    i
    for k = 1:size(PredictiveFeatures,2)
        PredictiveFeatures{k}
        [SD_predvals{i,k} SD_prednorm{i,k} SD_predstd{i,k} SD_predcount{i,k} SD_ROCs{i,k}]=PredictPatients(NewPatPosStruct(INDS{i}),100,'ELM_STRUCT',NewELMPosStruct,'display',false,'predictivefeatures',[{{'ResponderType','SD_method'}}, PredictiveFeatures{k}]);
        [my_predvals{i,k} my_prednorm{i,k} my_predstd{i,k} my_predcount{i,k} my_ROCs{i,k}]=PredictPatients(NewPatPosStruct(INDS{i}),100,'ELM_STRUCT',NewELMPosStruct,'display',false,'predictivefeatures',PredictiveFeatures{k});
        save AnotherRunData my_predvals my_prednorm my_predstd my_predcount my_ROCs SD_predvals SD_prednorm SD_predstd SD_predcount SD_ROCs
    end
    
end
% 
% PredCount = cat(1,my_predcount{:});
% PredNorm = cat(1,my_prednorm{:});
% PredSTD = cat(1,my_predstd{:});
% PredVals = cat(1,my_predvals{:});
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