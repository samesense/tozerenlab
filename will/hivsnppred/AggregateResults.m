function [agg_pred_vars agg_resp_vars ]=AggregateResults(INPUT_PRED,INPUT_RESP)
%   AggregateResults
%       A Helper function for DoLogitRegression.  This function takes a
%       matrix of binary vectors of individual patients and then finds each
%       unique observation and collapses them into a format for logistic
%       regression (mnrfit and glmfit).
%
%
%   [AGG_OBV AGG_RESP COL_NUMS]=AggregateResults(INPUT_PRED,INPUT_RESP)
%
%       INPUT_PRED      A matrix of binary observations.
%       INPUT_RESP      A vector of the binary response variable for each
%                       observation in INPUT_PRED.
%
%       AGG_OBV         Aggregated observations in which each row is unique
%                       and each column has variation>0.
%       AGG_RESP        Aggregated response in which:
%                           AGG_RESP(:,1)       # of RESP_VAR==0
%                           AGG_RESP(:,2)       # of RESP_VAR==1
%       COL_NUMS        The columns indices from the orriginal INPUT_PRED
%                       which are mirrored in the output.
%



INPUT_PRED(isnan(INPUT_PRED))=inf;

agg_pred_vars=unique(INPUT_PRED,'rows');
%init_col_nums=find(sum(isnan(agg_pred_vars))==0&sum(agg_pred_vars==1)~=size(agg_pred_vars,1)&sum(agg_pred_vars==0)~=size(agg_pred_vars,1));
%%%remove cols that are all NaN or have no variation

%%%Deal with an odd case where the unique-ed cols have no varaition
% still_bad_nums=~(sum(isnan(agg_pred_vars))==0&sum(agg_pred_vars==1)~=size(agg_pred_vars,1)&sum(agg_pred_vars==0)~=size(agg_pred_vars,1));
% while(any(still_bad_nums)||length(unique(agg_pred_vars))==length(agg_pred_vars))&&~isempty(agg_pred_vars,2)
%     init_col_nums=init_col_nums(~still_bad_nums);
%     agg_pred_vars=unique(INPUT_PRED(:,init_col_nums),'rows');
%     still_bad_nums=~(sum(isnan(agg_pred_vars))==0&sum(agg_pred_vars==1)~=size(agg_pred_vars,1)&sum(agg_pred_vars==0)~=size(agg_pred_vars,1));
% end

agg_resp_vars=zeros(size(agg_pred_vars,1),2);

spots=(1:size(INPUT_PRED,1))';
try
for k=1:length(agg_resp_vars)
    found_spots=sum(repmat(agg_pred_vars(k,:),[nnz(spots) 1])==INPUT_PRED(spots,:),2)==size(agg_pred_vars,2);
    agg_resp_vars(k,:)=[sum(INPUT_RESP(spots(found_spots))==0) sum(INPUT_RESP(spots(found_spots))==1)];
    spots=spots(~found_spots);
end
catch
    display('here')
end
agg_pred_vars(isinf(agg_pred_vars))=NaN;












end