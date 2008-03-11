function PRED_VALS=PredictPatients(PAT_STRUCT,NUM_REPS,varargin)
%   PredictPatients
%       Uses (leave-1/3)-out cross validation to determine the predictive
%       power of the BASE_CALLS and ELM vectors in predicting patient
%       response.  This method uses Nearest-Neighbor as defined by the
%       'Hamming' distance.
%
%
%       MEAN_AUC=PredictPatients(PAT_STRUCT,NUM_REPS)
%
%
%

DISPLAY_FLAG=false;

LEAVE_IN_FRAC=0.66;
MAX_FEATURES=10;


[temp_RX_data temp_RX_inds ...
    temp_resp_var RESP_inds ...
    temp_BASE_CALLS temp_BASE_CALLS_inds ...
    temp_CLINICAL temp_CLINICAL_inds ...
    temp_SNP_spots temp_SNP_spots_inds ...
    temp_ELM_simple temp_ELM_simple_inds...
    temp_ELM_vec temp_ELM_vec_inds...
    temp_ELM_annot temp_ELM_annot_inds]=...
    PatientStructHelper(PAT_STRUCT,...
    {'RX_vals','explodenumeric'},...
    {'IS_RESPONDER','leaveNumeric'},...
    {'BASE_CALLS','explodenumeric'},...
    {'NORM_CLINICAL_DATA','explodecells'},...
    {'SNP_SPOTS','explodeNumeric'},...
    {'ELM_simple','explodeNumeric'},...
    {'ELM_vec','explodeNumeric'},...
    {'ELM_annot','leaveCell'});


%%extract the RX data that indicates the start of the study
% mask=(temp_RX_data(:,1)==0);
% 
% temp_RX_data=temp_RX_data(mask,2:end);
% temp_RX_inds=temp_RX_inds(mask);
% 
% %%%%%Pull out desired triple therapy
% RX_names=PAT_STRUCT{1}.RX_names;
% %wanted_drugs={'3TC','AZT','IDV'};
% wanted_drugs={'3TC'};
% wanted_drug_mask=ismember(RX_names,wanted_drugs);
% 
% %mask=sum(temp_RX_data==repmat(wanted_drug_mask',[size(temp_RX_data,1) 1]),2)==length(RX_names);
% mask=temp_RX_data(:,find(wanted_drug_mask))&true;

mask=true(size(temp_RX_inds));
temp_RX_data=temp_RX_data(mask,:);
temp_RX_inds=temp_RX_inds(mask);

%%%%ensure that the indexes still match up properly
PAT_IND_MATCHES=intersect(intersect(intersect(intersect(temp_RX_inds,RESP_inds),temp_BASE_CALLS_inds),temp_CLINICAL_inds),temp_ELM_simple_inds);

BASE_CALLS=zeros(length(PAT_IND_MATCHES),size(temp_BASE_CALLS,2)*5);

INITIAL_CLINICAL=cell2mat(cellfun(@(x)x(1,2:end),temp_CLINICAL(PAT_IND_MATCHES),'uniformoutput',false));
FINAL_CLINICAL=cell2mat(cellfun(@(x)x(end,2:end),temp_CLINICAL(PAT_IND_MATCHES),'uniformoutput',false));

m=size(BASE_CALLS,2);

translate_vars=[1 2 3 4 16];
translate_mat=[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; NaN NaN NaN NaN NaN];

RESP_VAR=temp_resp_var(PAT_IND_MATCHES)&1;
unex_BASE_CALLS=temp_BASE_CALLS(PAT_IND_MATCHES,:);

ELM_features=[temp_ELM_simple(PAT_IND_MATCHES,:) temp_ELM_vec(PAT_IND_MATCHES,:)];

for i=1:length(PAT_IND_MATCHES)
    [TF LOC]=ismember(unex_BASE_CALLS(i,:),translate_vars);
    LOC(~TF)=6;
    BASE_CALLS(i,:)=reshape(translate_mat(LOC,:),[],m);
end

INITIAL_CLINICAL=cell2mat(cellfun(@(x)x(1,2:end),temp_CLINICAL(PAT_IND_MATCHES),'uniformoutput',false));
FINAL_CLINICAL=cell2mat(cellfun(@(x)x(end,2:end),temp_CLINICAL(PAT_IND_MATCHES),'uniformoutput',false));

ACTUAL_IMPROVED=(INITIAL_CLINICAL(:,1)-FINAL_CLINICAL(:,1))>0;

%RESP_VAR=ACTUAL_IMPROVED;

WANTED_SENS=0:0.01:1;


% RESP_VAR=[true(30,1); false(30,1)];
% test_data=double(rand(60,100)>0.5);
% test_data(:,10)=[rand(30,1)+.4>0.5; rand(30,1)-.4>0.5];
% test_data(:,15)=[rand(30,1)+.4>0.5; rand(30,1)-.4>0.5];
% test_data(:,5)=[rand(30,1)-.2>0.5; rand(30,1)+.2>0.5];
% CORR_SPOTS=zeros(NUM_REPS,size(test_data,2));



resp_inds=find(RESP_VAR);
non_resp_inds=find(~RESP_VAR);

num_resp_include=ceil(nnz(RESP_VAR)*LEAVE_IN_FRAC);
num_non_resp_include=ceil(nnz(~RESP_VAR)*LEAVE_IN_FRAC);

num_training=ceil(LEAVE_IN_FRAC*size(BASE_CALLS,1));


CORR_LABELS=cumsum([1 size(temp_RX_data,2) size(BASE_CALLS,2) size(temp_ELM_simple,2) size(temp_ELM_vec,2)]);

CLASS_PERF=classperf(RESP_VAR);

if DISPLAY_FLAG
    %kNN_fig_handle=figure;
    %SVM_fig_handle=figure;
    SWR_fig_handle=figure;
end

kNN_CORR_SPOTS=zeros(NUM_REPS,size(temp_RX_data,2)+size(BASE_CALLS,2)+size(ELM_features,2));
kNN_TRAIN_CLASS_CORRECT=zeros(NUM_REPS,1);

SWR_CORR_SPOTS=zeros(NUM_REPS,size(temp_RX_data,2)+size(BASE_CALLS,2)+size(ELM_features,2));
SWR_AUC_VALS=zeros(NUM_REPS,1);
SWR_SPEC_VALS=zeros(NUM_REPS,length(WANTED_SENS));
SWR_REG_VALS=zeros(NUM_REPS,size(temp_RX_data,2)+size(BASE_CALLS,2)+size(ELM_features,2));


for i=1:NUM_REPS
    shuffle_resp=resp_inds(randperm(length(resp_inds)));
    shuffle_non_resp=non_resp_inds(randperm(length(non_resp_inds)));


    TRAINING_FEATURES=[temp_RX_data(shuffle_resp(1:num_resp_include),:) BASE_CALLS(shuffle_resp(1:num_resp_include),:) ELM_features(shuffle_resp(1:num_resp_include),:); ...
        temp_RX_data(shuffle_non_resp(1:num_non_resp_include),:) BASE_CALLS(shuffle_non_resp(1:num_non_resp_include),:) ELM_features(shuffle_non_resp(1:num_non_resp_include),:)];

    TRAINING_RESP=[RESP_VAR(shuffle_resp(1:num_resp_include)); RESP_VAR(shuffle_non_resp(1:num_non_resp_include))];

    TESTING_FEATURES=[temp_RX_data(shuffle_resp(num_resp_include+1:end),:) BASE_CALLS(shuffle_resp(num_resp_include+1:end),:) ELM_features(shuffle_resp(num_resp_include+1:end),:); ...
        temp_RX_data(shuffle_non_resp(num_non_resp_include+1:end),:) BASE_CALLS(shuffle_non_resp(num_non_resp_include+1:end),:) ELM_features(shuffle_non_resp(num_non_resp_include+1:end),:)];

    TESTING_RESP=[RESP_VAR(shuffle_resp(num_resp_include+1:end)); RESP_VAR(shuffle_non_resp(num_non_resp_include+1:end))];

    %     TRAINING_FEATURES=[test_data(shuffle_resp(1:num_resp_include),:); test_data(shuffle_non_resp(1:num_non_resp_include),:)];
    %     TRAINING_RESP=[RESP_VAR(shuffle_resp(1:num_resp_include)); RESP_VAR(shuffle_non_resp(1:num_non_resp_include))];
    %
    %     TESTING_FEATURES=[test_data(shuffle_resp(num_resp_include+1:end),:); test_data(shuffle_non_resp(num_non_resp_include+1:end),:)];
    %     TESTING_RESP=[RESP_VAR(shuffle_resp(num_resp_include+1:end)); RESP_VAR(shuffle_non_resp(num_non_resp_include+1:end))];


    %    [I Z]=rankfeatures(TRAINING_FEATURES',TRAINING_RESP);
    %     RATIOS=RatioBSStoWSS(TRAINING_FEATURES,TRAINING_RESP,[0 1]);

    %     [Y I]=sort(RATIOS,'descend');
    %     spot=find(~isnan(Y),1);

    %    figure

    %
% 
%    good_features=find(sum(isnan(TRAINING_FEATURES),1)==0);

%    [I Z]=rankfeatures(TRAINING_FEATURES(:,good_features)',TRAINING_RESP,'criterion','roc');
%    I=good_features(I);
%     [K NUM_FEATURES]=FindOptkNN(TRAINING_FEATURES(:,I),TRAINING_RESP);
%     kNN_CORR_SPOTS(i,I(1:NUM_FEATURES))=1;
% 
%     Nearest_Testing_vals=knnclassify(TESTING_FEATURES(:,I(1:NUM_FEATURES)),TRAINING_FEATURES(:,I(1:NUM_FEATURES)),TRAINING_RESP,K,'hamming','consensus');
% 
%     CLASS_PERF=classperf(CLASS_PERF,Nearest_Testing_vals,[shuffle_resp(num_resp_include+1:end); shuffle_non_resp(num_non_resp_include+1:end)]);
%     kNN_TRAIN_CLASS_CORRECT(i)=CLASS_PERF.LastCorrectRate;
%     
%     
%     figure(kNN_fig_handle)
%     subplot(2,2,4), barh(mean(kNN_CORR_SPOTS(1:i,:),1))
%     axis([0 1 0 size(kNN_CORR_SPOTS,2)])
%     set(gca,'ytick',CORR_LABELS,'Yticklabel',{'RX Data','SNP Data','Simple ELM','Positional ELM'});
%     xlabel('Frequency Picked')
%     subplot(2,2,3), hist(kNN_TRAIN_CLASS_CORRECT(1:i),0:0.01:1)
%     ylabel('Frequency')
%     xlabel('CorrectRate')
      
    
NAN_MASK=find(sum(isnan(TRAINING_FEATURES))==0);

intial_inds=FindOptSWR(TRAINING_FEATURES(:,NAN_MASK),TRAINING_RESP);

[B_values,se,pval,inmodel,stats,nextstep,history]=stepwisefit(TRAINING_FEATURES(:,NAN_MASK),TRAINING_RESP,'inmodel',ismember(NAN_MASK,intial_inds),'display','off');

%[AGG_OBV_TEST AGG_RESP_TEST]=AggregateResults(TESTING_FEATURES(:,I(inmodel)),TESTING_RESP);

VALS=glmval([B_values(inmodel);0],TESTING_FEATURES(:,NAN_MASK(inmodel)),'identity');
           
[temp_AUC SWR_SPEC_VALS(i,:)]=CalculateROC(VALS,TESTING_RESP,WANTED_SENS);

SWR_AUC_VALS(i)=abs(temp_AUC-0.5);


SWR_REG_VALS(i,NAN_MASK(inmodel))=sign(B_values(inmodel));

if DISPLAY_FLAG
    figure(SWR_fig_handle)
    subplot(2,2,2), bar(mean(SWR_REG_VALS(1:i,:),1))
    axis([0 size(SWR_REG_VALS,2) -1.1 1.1])
    % title(['Mean AUC: ' num2str(mean(SWR_AUC_VALS(1:i)))])

    AUC_VALS_HIST=hist(SWR_AUC_VALS(1:i),0:0.05:1);
    subplot(2,2,4), bar(0:0.05:1,AUC_VALS_HIST/i)
    axis([0 0.5 0 1])
    drawnow
end










    %
    %    figure
    %    last_spot=find(diff(Y(spot:end))==0,1);
    %    subplot(2,2,1), plot(1:(length(Y)-spot+1),Y(spot:end),'b',last_spot,Y(spot+last_spot),'kx')
    %
    %
    %     [AGG_OBV AGG_RESP]=AggregateResults(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),TRAINING_RESP);
    %
    %     B_values=glmfit(AGG_OBV,[AGG_RESP(:,2) sum(AGG_RESP,2)],'binomial','link','logit');
    %
    %     VALS=glmval(B_values,TESTING_FEATURES(:,I(spot:spot+MAX_FEATURES)),'logit');
    %
    %     t_vals=glmval(B_values,TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),'logit');
    %
    % %     [auc spec]=CalculateROC(t_vals,TRAINING_RESP,WANTED_SENS);
    % %    subplot(2,2,4), plot(1-spec,WANTED_SENS);
    % %    title('Training Data Logistic Regression')
    %
    %    [auc spec]=CalculateROC(VALS,TESTING_RESP,WANTED_SENS);
    %    subplot(2,1,2), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    %    title(['Testing Data Logistic Regression, AUC:' num2str(auc)])
    %
    %
    %
    %
    %
    % figure
    %
    % [xsup,w,b,pos,timeps,alpha,obj]=svmclass(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),TRAINING_RESP*2-1,100,10,'gaussian',.5);
    % [y]=svmval(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),xsup,w,b,'gaussian',.5);
    %
    % [auc spec]=CalculateROC(y,TRAINING_RESP,WANTED_SENS);
    % subplot(4,2,1), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    % title('SVM Gaussian 0.5 training')
    % legend(['AUC: ' num2str(auc)],'Random')
    %
    % [y]=svmval(TESTING_FEATURES(:,I(spot:spot+MAX_FEATURES)),xsup,w,b,'gaussian',.5);
    %
    % [auc spec]=CalculateROC(y,TESTING_RESP,WANTED_SENS);
    % subplot(4,2,2), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    % title('SVM Gaussian 0.5 testing')
    % legend(['AUC: ' num2str(auc)],'Random')
    %
    %
    % [xsup,w,b,pos,timeps,alpha,obj]=svmclass(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),TRAINING_RESP*2-1,100,10,'gaussian',1);
    % [y]=svmval(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),xsup,w,b,'gaussian',1);
    %
    % [auc spec]=CalculateROC(y,TRAINING_RESP,WANTED_SENS);
    % subplot(4,2,3), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    % title('SVM Gaussian 1 training')
    % legend(['AUC: ' num2str(auc)],'Random')
    %
    % [y]=svmval(TESTING_FEATURES(:,I(spot:spot+MAX_FEATURES)),xsup,w,b,'gaussian',1);
    %
    % [auc spec]=CalculateROC(y,TESTING_RESP,WANTED_SENS);
    % subplot(4,2,4), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    % title('SVM Gaussian 1 testing')
    % legend(['AUC: ' num2str(auc)],'Random')
    %
    %
    % [xsup,w,b,pos,timeps,alpha,obj]=svmclass(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),TRAINING_RESP*2-1,100,10,'poly',1);
    % [y]=svmval(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),xsup,w,b,'poly',1);
    %
    % [auc spec]=CalculateROC(y,TRAINING_RESP,WANTED_SENS);
    % subplot(4,2,5), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    % title('SVM poly 1 training')
    % legend(['AUC: ' num2str(auc)],'Random')
    %
    % [y]=svmval(TESTING_FEATURES(:,I(spot:spot+MAX_FEATURES)),xsup,w,b,'poly',1);
    %
    % [auc spec]=CalculateROC(y,TESTING_RESP,WANTED_SENS);
    % subplot(4,2,6), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    % title('SVM poly 1 testing')
    % legend(['AUC: ' num2str(auc)],'Random')
    %
    %
    % [xsup,w,b,pos,timeps,alpha,obj]=svmclass(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),TRAINING_RESP*2-1,100,10,'poly',5);
    % [y]=svmval(TRAINING_FEATURES(:,I(spot:spot+MAX_FEATURES)),xsup,w,b,'poly',5);
    %
    % [auc spec]=CalculateROC(y,TRAINING_RESP,WANTED_SENS);
    % subplot(4,2,7), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    % title('SVM poly 5 training')
    % legend(['AUC: ' num2str(auc)],'Random')
    %
    % [y]=svmval(TESTING_FEATURES(:,I(spot:spot+MAX_FEATURES)),xsup,w,b,'poly',5);
    %
    % [auc spec]=CalculateROC(y,TESTING_RESP,WANTED_SENS);
    % subplot(4,2,8), plot(1-spec,WANTED_SENS,WANTED_SENS,WANTED_SENS);
    % title('SVM poly 5 testing')
    % legend(['AUC: ' num2str(auc)],'Random')
    %
    %

    %     [AGG_test_AUC(i) AGG_test_SPEC(i,:)]=CalculateROC(VALS,TESTING_RESP,WANTED_SENS);
    %     [t_AUC t_SPEC]=CalculateROC(t_vals,TRAINING_RESP,WANTED_SENS);
    %
    %     subplot(2,2,2), plot(1-AGG_test_SPEC(i,:),WANTED_SENS,WANTED_SENS,WANTED_SENS);
    %     title('Current Plot')
    %     legend(['AUC: ' num2str(AGG_test_AUC(i))],'Random')
    %
    %     subplot(2,2,3), plot(1-mean(AGG_test_SPEC(1:i,:),1),WANTED_SENS,WANTED_SENS,WANTED_SENS);
    %     title('Accumulated Plot')
    %     legend(['AUC: ' num2str(mean(AGG_test_AUC(1:i)))],'Random')
    %
    %     subplot(2,2,4),plot(1-t_SPEC(:),WANTED_SENS,WANTED_SENS,WANTED_SENS);
    %     legend(['AUC: ' num2str(t_AUC)],'Random')
    %    subplot(2,2,4), bar(mean(CORR_SPOTS(1:i,:),1))
    %   set(gca,'xtick',CORR_LABELS)



end

PRED_VALS=mean(SWR_REG_VALS,1);

    function [K_val F_num]=FindOptkNN(FEATURES,RESP)

        MAX_K=30;
        [NUM_OBSERVATIONS MAX_FEATURES]=size(FEATURES);

        OUTPUT_MAT=zeros(MAX_FEATURES,MAX_K);


        parfor (K_IND = 1:MAX_K)
            temp_CLASS_PERF=classperf(RESP);

            F_IND=1;
            OUTPUT_SLICE=OUTPUT_MAT(:,K_IND);

            break_var=true;
            while(F_IND<MAX_FEATURES&&break_var)
                if numel(unique(reshape(FEATURES(:,1:F_IND),1,[])))==1
                    F_IND=F_IND+1;    
                    continue
                end
                    
                temp_mat=zeros(20,1);
                for inside=1:20
                    [this_train this_test]=crossvalind('leaveMout',NUM_OBSERVATIONS,round(.1*NUM_OBSERVATIONS));
                    this_classified=knnclassify(FEATURES(this_test,1:F_IND),FEATURES(this_train,1:F_IND),RESP(this_train),K_IND,'hamming','consensus');
                    temp_CLASS_PERF=classperf(temp_CLASS_PERF,this_classified,this_test);
                end

                OUTPUT_SLICE(F_IND)=temp_CLASS_PERF.CorrectRate;

                if F_IND>10&&(mean(diff(OUTPUT_SLICE(F_IND-4:F_IND)))<0.00001)%&&(rand*F_IND)/MAX_FEATURES>0.25)
                    break_var=false;
                end
                F_IND=F_IND+1;
            end
            
            OUTPUT_MAT(:,K_IND)=OUTPUT_SLICE;
            
        end

        figure(kNN_fig_handle);
        subplot(2,2,1), pcolor(OUTPUT_MAT(1:max(sum(OUTPUT_MAT>0)),:))
        shading flat
        ylabel('# Features')
        xlabel('K-NearestNeighbors')
        caxis([0 1])
        colorbar

        %        AUC=AUC(:,2:end);

        [max_val max_ind]=max(OUTPUT_MAT(:));

        [F_num K_val]=ind2sub(size(OUTPUT_MAT),max_ind);
        %        F_num=F_num+2;

    end


function [FINAL_INDS]=FindOptSWR(FEATURES,RESP)

        NUM_OPT_REPS=50;
        [NUM_OBSERVATIONS MAX_FEATURES]=size(FEATURES);

        OUTPUT_MAT=zeros(NUM_OPT_REPS,size(SWR_CORR_SPOTS,2));
        AUC_MAP=zeros(NUM_OPT_REPS,1);


        parfor (IND = 1:NUM_OPT_REPS)
           
            [this_train this_test]=crossvalind('leaveMout',NUM_OBSERVATIONS,round(.1*NUM_OBSERVATIONS));
            
%            [AGG_OBS AGG_RESP]=AggregateResults(FEATURES(this_train,:),RESP(this_train));
            

%            [B_values,se,pval,inmodel,stats,nextstep,history]=stepwisefit(AGG_OBS,AGG_RESP(:,1)-AGG_RESP(:,2),'inmodel',rand(MAX_FEATURES,1)<0.1,'display','off');
            [B_values,se,pval,inmodel,stats,nextstep,history]=stepwisefit(FEATURES(this_train,:),RESP(this_train),'inmodel',rand(MAX_FEATURES,1)<0.1,'display','off');
            VALS=glmval([B_values(inmodel);0],FEATURES(this_test,inmodel),'identity');
            
            AUC_MAP(IND)=abs(CalculateROC(VALS,RESP(this_test))-0.5);
            
            OUTPUT_SLICE=zeros(1,size(SWR_CORR_SPOTS,2));
            OUTPUT_SLICE(NAN_MASK(inmodel))=1;
            
            OUTPUT_MAT(IND,:)=OUTPUT_SLICE;
            
        end

        
        auc=hist(AUC_MAP,0:0.05:1);
        NORM_AUC=auc/NUM_OPT_REPS;

        SPOT_FREQ=mean(OUTPUT_MAT,1);
        AUC_NORM_SPOT_FREQ=mean(OUTPUT_MAT.*repmat(AUC_MAP,[1 size(OUTPUT_MAT,2)]),1);
        if DISPLAY_FLAG
            figure(SWR_fig_handle);

            subplot(2,2,1), bar(0:0.05:1,NORM_AUC)
            axis([0 0.5 0 1])
            ylabel('Frequency')
            xlabel('AUC Val')

            subplot(2,2,3), barh([SPOT_FREQ',AUC_NORM_SPOT_FREQ']);
            set(gca,'ytick',CORR_LABELS,'Yticklabel',{'RX Data','SNP Data','Simple ELM','Positional ELM'});
            legend('SPOT OCCURANCE','AUC Normalized')
        end
        [Y_val I_val]=sort(AUC_NORM_SPOT_FREQ,'descend');
        
        FINAL_INDS=I_val(1:10);
    end


end