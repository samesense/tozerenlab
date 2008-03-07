function DoLogitRegression(PAT_STRUCT)

NUM_REPS=100;
LEAVE_IN_FRAC=0.66;
NUM_FEATURES=10;
REF_SEQ_LENGTH=9181;

load ANNOT_DATA HIV_ANNOT_CELL


% [temp_RX_data temp_RX_ind ]=...
%     PatientStructHelper(PAT_STRUCT,...
%     {'RX_vals','explodenumeric'});
%
% %%extract the RX data that indicates the start of the study
% mask=(temp_RX_data(:,1)==0);
%
% temp_RX_data=temp_RX_data(mask,2:end);
% temp_RX_ind=temp_RX_ind(mask);
%
% %%%%%Pull out desired triple therapy
% RX_names=PAT_STRUCT{1}.RX_names;
% wanted_drugs={'3TC','AZT','IDV'};
% wanted_drug_mask=ismember(RX_names,wanted_drugs);
%
% mask=sum(temp_RX_data==repmat(wanted_drug_mask',[size(temp_RX_data,1) 1]),2)==length(RX_names);
%
% [PAT_STRUCT]=MakeSNPCalls(HIV_REF,bkg_data,...
%     'alignments',ALL_ALIGNMENTS,...
%     'use_distributed',false,...
%     'pat_struct_out',true,...
%     'pat_struct_input',PAT_STRUCT(temp_RX_ind),...
%     'nested',true,...
%     'no_translate',true,...
%     'hidden',mask);

[temp_RX_data temp_RX_ind ...
    temp_resp_var RESP_ind ...
    temp_BASE_CALLS temp_BASE_CALLS_ind...
    temp_SNP_CALLS temp_SNP_CALLS_ind...
    temp_CLINICAL temp_CLINICAL_ind ...
    temp_SNP_spots temp_SNP_spots_ind]=...
    PatientStructHelper(PAT_STRUCT,...
    {'RX_vals','explodenumeric'},...
    {'IS_RESPONDER','leaveNumeric'},...
    {'BASE_CALLS','explodenumeric'},...
    {'IS_SNPs','explodenumeric'},...
    {'NORM_CLINICAL_DATA','explodecells'},...
    {'SNP_SPOTS','explodeNumeric'});


SNP_SPOTS=temp_SNP_spots(1,:);
%%extract the RX data that indicates the start of the study
mask=(temp_RX_data(:,1)==0);

temp_RX_data=temp_RX_data(mask,2:end);
temp_RX_ind=temp_RX_ind(mask);

%%%%%Pull out desired triple therapy
RX_names=PAT_STRUCT{1}.RX_names;
wanted_drugs={'3TC','AZT','IDV'};
wanted_drug_mask=ismember(RX_names,wanted_drugs);

mask=sum(temp_RX_data==repmat(wanted_drug_mask',[size(temp_RX_data,1) 1]),2)==length(RX_names);

temp_RX_data=temp_RX_data(mask,:);
temp_RX_ind=temp_RX_ind(mask);


%%%%ensure that the indexes still match up properly
PAT_IND_MATCHES=intersect(intersect(intersect(intersect(temp_RX_ind,RESP_ind),temp_BASE_CALLS_ind),temp_SNP_CALLS_ind),temp_CLINICAL_ind);

BASE_CALLS=zeros(length(PAT_IND_MATCHES),size(temp_BASE_CALLS,2)*5);
SNP_CALLS=zeros(length(PAT_IND_MATCHES),size(temp_SNP_CALLS,2));





INITIAL_CLINICAL=cell2mat(cellfun(@(x)x(1,2:end),temp_CLINICAL(PAT_IND_MATCHES),'uniformoutput',false));
FINAL_CLINICAL=cell2mat(cellfun(@(x)x(end,2:end),temp_CLINICAL(PAT_IND_MATCHES),'uniformoutput',false));

ACTUAL_IMPROVED=(INITIAL_CLINICAL(:,1)-FINAL_CLINICAL(:,1))>0;



m=size(BASE_CALLS,2);

translate_vars=[1 2 3 4 16];
translate_mat=[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; NaN NaN NaN NaN NaN];

RESP_VAR=temp_resp_var(PAT_IND_MATCHES)&1;
SNP_CALLS=temp_SNP_CALLS(PAT_IND_MATCHES,:);


unex_BASE_CALLS=temp_BASE_CALLS(PAT_IND_MATCHES,:);

for i=1:length(PAT_IND_MATCHES)
    [TF LOC]=ismember(unex_BASE_CALLS(i,:),translate_vars);
    LOC(~TF)=6;
    BASE_CALLS(i,:)=reshape(translate_mat(LOC,:),[],m);
end

WANTED_SENS=0:0.01:1;

AGG_train_AUC=zeros(NUM_REPS,1);
AGG_train_SPEC=zeros(NUM_REPS,length(WANTED_SENS));

AGG_test_AUC=zeros(NUM_REPS,1);
AGG_test_SPEC=zeros(NUM_REPS,length(WANTED_SENS));

CORR_SPOTS=NaN(NUM_REPS,2+size(BASE_CALLS,2));



num_resp_include=ceil(nnz(RESP_VAR)*LEAVE_IN_FRAC);
num_non_resp_include=ceil(nnz(~RESP_VAR)*LEAVE_IN_FRAC);

num_training=ceil(LEAVE_IN_FRAC*size(BASE_CALLS,1));


for i=1:NUM_REPS

    %%%%%%%%%%Pick random training/testing data
    shuff_inds=randperm(size(BASE_CALLS,1));

    this_TRAIN_inds=shuff_inds(1:num_training);
    not_TRAIN_inds=shuff_inds(num_training+1:end);
    
    
    this_COL_NUMS=SelectFeatures(BASE_CALLS(this_TRAIN_inds,:),RESP_VAR(this_TRAIN_inds));
    
    TRAINING_DATA=[INITIAL_CLINICAL(this_TRAIN_inds,:) BASE_CALLS(this_TRAIN_inds,this_COL_NUMS)];

%    TRAINING_RESP=RESP_VAR(this_TRAIN_inds);

    B_values=glmfit(TRAINING_DATA,FINAL_CLINICAL(this_TRAIN_inds,1));
    
    CORR_SPOTS(i,[1 2 this_COL_NUMS+2])=B_values(2:end);
    
    PRED_TRAIN=glmval(B_values,TRAINING_DATA,'identity');
    PRED_IMPROV_TRAIN=PRED_TRAIN-INITIAL_CLINICAL(this_TRAIN_inds,1);
    [AGG_train_AUC(i) AGG_train_SPEC(i,:)]=CalculateROC(PRED_IMPROV_TRAIN,RESP_VAR(this_TRAIN_inds),WANTED_SENS);
    
    PRED_TEST=glmval(B_values,[INITIAL_CLINICAL(not_TRAIN_inds,:) BASE_CALLS(not_TRAIN_inds,this_COL_NUMS)],'identity');
    PRED_IMPROV_TEST=PRED_TEST-INITIAL_CLINICAL(not_TRAIN_inds,1);
    [AGG_test_AUC(i) AGG_test_SPEC(i,:)]=CalculateROC(PRED_IMPROV_TEST,RESP_VAR(not_TRAIN_inds),WANTED_SENS);
    
end


train_AUC=mean(AGG_train_AUC);
test_AUC=mean(AGG_test_AUC);

train_SPEC=mean(AGG_train_SPEC);
test_SPEC=mean(AGG_test_SPEC);  

figure
plot(1-test_SPEC,WANTED_SENS,'r',1-train_SPEC,WANTED_SENS,'b',0:.1:1,0:.1:1,'g')
legend(['Testing ' num2str(test_AUC)],['Training ' num2str(train_AUC)],'Random')


CORR_SPOTS=CORR_SPOTS.*repmat(1./nanmax(abs(CORR_SPOTS),[],2),[1 size(CORR_SPOTS,2)]);
POSS_VALS=CORR_SPOTS(:,3:end);
POSS_VALS(POSS_VALS<0)=NaN;

NEG_VALS=CORR_SPOTS(:,3:end);
NEG_VALS(NEG_VALS>0)=NaN;

POSS_CORR_VALS=reshape(nanmean(POSS_VALS),5,[]);
NEG_CORR_VALS=reshape(nanmean(NEG_VALS),5,[]);

figure

subplot(2,1,1),bar(1:length(SNP_SPOTS),POSS_CORR_VALS');
hold;
bar(1:length(SNP_SPOTS),NEG_CORR_VALS');
title('Correlation With Drug Response')
legend('A','C','G','T','Gap','location','northwest')
axis([0 nnz(SNP_SPOTS) -1.1 1.1])

subplot(2,1,2), plot(1:length(SNP_SPOTS),ones(1,length(SNP_SPOTS)))
AnnotFig(SNP_SPOTS,HIV_ANNOT_CELL)
%title('Negatively Correlated with Drug Response')
%legend('A','C','G','T','Gap','location','northwest')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HELPER FUNCTIONS

    function [AGG_PRED_VARS AGG_RESP_VARS COL_NUM]=CollectData(INPUT_MAT,INPUT_RESP)
        %find the top NUM_FEATURES ranked columns
        ratios=RatioBSStoWSS(INPUT_MAT,INPUT_RESP);
        ratios(isnan(ratios))=-inf;
        [sorted_ratios ranks]=sort(ratios,'descend');
        
        COL_NUM=ranks(1:min([find(isinf(sorted_ratios),1) NUM_FEATURES]));

        %%%%%%Conform the input results for mnrfit
        [AGG_PRED_VARS AGG_RESP_VARS]=AggregateResults(INPUT_MAT(:,COL_NUM),INPUT_RESP);

    end

    function [COL_NUMS]=SelectFeatures(INPUT_MAT,INPUT_RESP)
        ratios=RatioBSStoWSS(INPUT_MAT,INPUT_RESP);
        ratios(isnan(ratios))=-inf;
        [sorted_ratios ranks]=sort(ratios,'descend');
        
        COL_NUMS=ranks(1:min([find(isinf(sorted_ratios),1) NUM_FEATURES]));
    end




end