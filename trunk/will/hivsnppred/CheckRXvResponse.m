function CheckRXvResponse(PAT_STRUCT)


[temp_RX_data temp_RX_ind ...
    resp_var RESP_ind]=...
    PatientStructHelper(PAT_STRUCT,...
    {'RX_vals','explodenumeric'},...
    {'IS_RESPONDER','leaveNumeric'});

%%extract the RX data that indicates the start of the study
mask=temp_RX_data(:,1)==0;

temp_RX_data=temp_RX_data(mask,2:end);
temp_RX_ind=temp_RX_ind(mask);

%%%%ensure that the indexes still match up properly
[PAT_IND_MATCHES t_RX_ind t_resp_ind]=intersect(temp_RX_ind,RESP_ind);

FINAL_RX_DATA=temp_RX_data(t_RX_ind,:);
FINAL_RESP_DATA=resp_var(t_resp_ind);









































end