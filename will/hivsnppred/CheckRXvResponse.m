
[RX_MASK INDS]=BreakoutPatients(PatPosStruct,100);


% All_predvals=cell(length(INDS),1);
% All_prednorm=cell(length(INDS),1);
% All_predstd=cell(length(INDS),1);
for i = 2:length(INDS)
    i
    [All_predvals{i} All_prednorm{i} All_predstd{i}]=PredictPatients(PatPosStruct(INDS{i}),500,'ELM_STRUCT',ELM_POS_STRUCT,'display',true);
    save RunData All_predvals All_prednorm All_predstd
    
end
