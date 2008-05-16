% 
% [RX_MASK INDS]=BreakoutPatients(PatPosStruct,200);
% 
% 
% All_predvals=cell(length(INDS),1);
% All_prednorm=cell(length(INDS),1);
% All_predstd=cell(length(INDS),1);
% All_predcount=cell(length(INDS),1);
% for i = 1:length(INDS)
%     i
%     [All_predvals{i} All_prednorm{i} All_predstd{i} All_predcount{i} ]=PredictPatients(PatPosStruct(INDS{i}),100,'ELM_STRUCT',ELM_POS_STRUCT,'display',false);
%     save RunData All_predvals All_prednorm All_predstd All_predcount
%     
% end

% PredCount = cat(1,All_predcount{:});
% PredNorm = cat(1,All_prednorm{:});
% PredSTD = cat(1,All_predstd{:});
% PredVals = cat(1,All_predvals{:});

[FEATURES RESPONSE PAT_INDEX FEATURE_NAMES FEATURE_INDS]=GetPatientFeatures(PatPosStruct,ELM_POS_STRUCT,'DrugRegimine','BaseCalls','SimpleELM','PositionalELM','PositionalPWM');



InterestingRows = find(any(PredCount>20&abs(PredVals)>.9&PredNorm>0.05,1));
for i = 1:length(FEATURE_NAMES)
    FEATURE_NAMES{i}(FEATURE_NAMES{i}=='_') = ' ';
end


% xlswrite('PrelimOutput.xlsx',FEATURE_NAMES(InterestingRows))
% xlswrite('PrelimOutput.xlsx',PredNorm(:,InterestingRows).*sign(PredVals(:,InterestingRows)),1,'A2')


normFeatures = bsxfun(@rdivide,FEATURES,nanmax(FEATURES));


spots = zeros(length(InterestingRows),2);
names = FEATURE_NAMES(InterestingRows);

for i = 1:length(names)
    bounder_spots=find(names{i}=='['|names{i}==']');
    if ~isempty(bounder_spots)
        spots(i,:) = str2num(names{i}(bounder_spots(1):bounder_spots(2)));
    end
end

goodSpots = spots(:,1)~=0;

myColor = usercolormap([1 1 1],[0 0 1]);

HIVFigGen(HIV_STRUCT,spots(goodSpots,:)+rand(nnz(goodSpots),2),nanmean(normFeatures(:,InterestingRows(goodSpots))),names(goodSpots),'stacked_windows','use_aa',true,'colormap',myColor);