function [varargout]=PredictPatients(PAT_STRUCT,NUM_REPS,varargin)
%   PredictPatients
%       Uses (leave-1/3)-out cross validation to determine the predictive
%       power of the BASE_CALLS and ELM vectors in predicting patient
%       response.
%
%       [PRED_VALS PRED_NORM PRED_STD]=PredictPatients(PAT_STRUCT,NUM_REPS)
%
%       PAT_STRUCT          A patient structure for testing.
%
%       NUM_REPS            The number of independant replications of the
%                           trianing/testing.
%
%       PRED_VALS           The average of the SIGNUM of each feature
%                           variable.
%
%       PRED_NORM           The average magnitude of each feature variable.
%       PRED_STD            The standard-deviation.
%
%       Optional Parameters
%
%       LeaveInFrac         The fraction of samples to use a training data
%                           in each iteration. DEFAULT = 0.66
%
%       Display             Toggles a figure displaying the prediction
%                           process as it happens. DEFAULT = false
%
%       NumCrossValid       The number of Cross-Validations to perform on
%                           the training data to determine the seeds for
%                           the training set. DEFAULT = 50
%
%       DistComputing       Toggles the use of distributed computing
%                           toolbox. DEFAULT = true
%
%       Method              Which method to for classification:
%                           'NearestNeighbor', 'kNN'
%                           'StepWise Linear Regression','SWR'    DEFAULT
%
%                           'NearestNeighbor' only returns PRED_VAR
%
%
%       See also: RunMultiAnal, PatientStructHelper, knnclassify,
%       stepwisefit, CalculateROC.
%
%



%%%%%%%%%%%%%%%%%SET DEFAULT VALUES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DISPLAY_FLAG=false;
LEAVE_IN_FRAC=0.66;
MAX_FEATURES=10;
NUM_OPT_REPS=50;
USE_DISTCOMP_FLAG=true;
kNN_FLAG=false;
SWR_FLAG=true;
varargout=cell(3,1);


%%%%%%%%%%%%%%%%%PARSE INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            
            case 'leaveinfrac'
                if isnumeric(varargin{i+1})&&isscalar(varargin{i+1})&&varargin{i+1}>0&&varargin{i+1}<1
                    LEAVE_IN_FRAC=varargin{i+1};
                else
                    error('PredictPatients:BAD_LEAVEINFRAC','Arguement to LeaveInFrac must be a numeric-scalar between 0 and 1.')
                end
            
            case 'numcrossvalid'
                if isscalar(varargin{i+1})&&varargin{i+1}>0
                    NUM_OPT_REPS=round(varargin{i+1});
                else
                    error('PredictPatients:BAD_NUMCROSSVALID','Arguement to NumCrossValis must be a scalar integer.')
                end
            
            case 'display'
                if islogical(varargin{i+1})&&isscalar(varargin{i+1})
                    DISPLAY_FLAG=varargin{i+1};
                else
                    error('PredictPatients:BAD_DIPLAY','Arguement to Display must be a logical.')
                end
                
            case 'distcomputing'
                if islogical(varargin{i+1})&&isscalar(varargin{i+1})
                    USE_DISTCOMP_FLAG=varargin{i+1};
                else
                    error('PredictPatients:BAD_DISTCOMP','Arguement to DistComputing must be a logical.')
                end
                
            case 'method'
                switch lower(varargin{i+1})
                    case {'nearestneighbor', 'knn'}
                        kNN_FLAG=true;
                        SWR_FLAG=false;
                        varargout=cell(1);
                        if nargout>1
                            warning('PredictPatients:kNN_OUTPUT','NearestNeighbor only returns PRED_VAR')
                        end
                    case {'stepwise linear legression','swr'}
                        kNN_FLAG=false;
                        SWR_FLAG=true;
                        varargout=cell(3,1);
                    otherwise
                        error('PredictPatients:UNKNOWN_METHOD','An unknown Method was provided: %s',varargin{i+1})
                end
            
            otherwise
                error('PredictPatients:UNKNOWN_ARG','An unknown argument was provided: %s',varargin{i})
        end
    end
end


%%%%%%%%%%%%%%%%%%EXTRACT DATA FROM PATIENT STRUCTURE%%%%%%%%%%%%%%%%%%%%%%
[temp_RX_data temp_RX_inds ...
    temp_resp_var RESP_inds ...
    temp_BASE_CALLS temp_BASE_CALLS_inds ...
    temp_CLINICAL temp_CLINICAL_inds ...
    temp_SNP_spots temp_SNP_spots_inds ...
    temp_ELM_simple temp_ELM_simple_inds...
    temp_ELM_vec temp_ELM_vec_inds]=...
    PatientStructHelper(PAT_STRUCT,...
    {'RX_vals','explodenumeric'},...
    {'IS_RESPONDER','leaveNumeric'},...
    {'BASE_CALLS','explodenumeric'},...
    {'NORM_CLINICAL_DATA','explodecells'},...
    {'SNP_SPOTS','explodeNumeric'},...
    {'ELM_simple','explodeNumeric'},...
    {'ELM_vec','explodeNumeric'});


if all(temp_resp_var)||all(~temp_resp_var)
    warning('PredictPatients:SINGLE_CLASS','This dataset only contained One Class, Cannot classify.')
    
    return
end



%%%%ensure that the indexes still match up properly
PAT_IND_MATCHES=intersect(intersect(intersect(intersect(temp_RX_inds,RESP_inds),temp_BASE_CALLS_inds),temp_CLINICAL_inds),temp_ELM_simple_inds);

BASE_CALLS=zeros(length(PAT_IND_MATCHES),size(temp_BASE_CALLS,2)*5);

m=size(BASE_CALLS,2);

translate_vars=[1 2 3 4 16];
translate_mat=[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; NaN NaN NaN NaN NaN];

RESP_VAR=temp_resp_var(PAT_IND_MATCHES)&1;
unex_BASE_CALLS=temp_BASE_CALLS(PAT_IND_MATCHES,:);

ELM_features=[temp_ELM_simple(PAT_IND_MATCHES,:) temp_ELM_vec(PAT_IND_MATCHES,:)];

%Explode BASE_CALLS into a feature vector
for i=1:length(PAT_IND_MATCHES)
    [TF LOC]=ismember(unex_BASE_CALLS(i,:),translate_vars);
    LOC(~TF)=6;
    BASE_CALLS(i,:)=reshape(translate_mat(LOC,:),[],m);
end

WANTED_SENS=0:0.01:1;

CORR_LABELS=cumsum([1 size(temp_RX_data,2) size(BASE_CALLS,2) size(temp_ELM_simple,2) size(temp_ELM_vec,2)]);

CLASS_PERF=classperf(RESP_VAR);

if DISPLAY_FLAG
    kNN_fig_handle=figure;
    %SVM_fig_handle=figure;
    SWR_fig_handle=figure;
end

kNN_CORR_SPOTS=zeros(NUM_REPS,size(temp_RX_data,2)+size(BASE_CALLS,2)+size(ELM_features,2));
kNN_TRAIN_CLASS_CORRECT=zeros(NUM_REPS,1);


SWR_CORR_SPOTS=zeros(NUM_REPS,size(temp_RX_data,2)+size(BASE_CALLS,2)+size(ELM_features,2));
SWR_AUC_VALS=zeros(NUM_REPS,1);
SWR_SPEC_VALS=zeros(NUM_REPS,length(WANTED_SENS));
SWR_REG_VALS=zeros(NUM_REPS,size(temp_RX_data,2)+size(BASE_CALLS,2)+size(ELM_features,2));
SWR_NORM_VALS=zeros(NUM_REPS,size(temp_RX_data,2)+size(BASE_CALLS,2)+size(ELM_features,2));


ALL_groups={'NR','R'};

CLASSES=ALL_groups(RESP_VAR+1);

for i=1:NUM_REPS

    [this_train_var this_test_var]=crossvalind('holdout',CLASSES,1-LEAVE_IN_FRAC,'classes',ALL_groups);


    TRAINING_FEATURES=[temp_RX_data(this_train_var,:) BASE_CALLS(this_train_var,:) ELM_features(this_train_var,:)];
    TRAINING_RESP=RESP_VAR(this_train_var);

    TESTING_FEATURES=[temp_RX_data(this_test_var,:) BASE_CALLS(this_test_var,:) ELM_features(this_test_var,:)];
    TESTING_RESP=RESP_VAR(this_test_var);

    if kNN_FLAG

        good_features=find(sum(isnan(TRAINING_FEATURES),1)==0);

        [I Z]=rankfeatures(TRAINING_FEATURES(:,good_features)',TRAINING_RESP,'criterion','roc');
        I=good_features(I);
        [K NUM_FEATURES]=FindOptkNN(TRAINING_FEATURES(:,I),TRAINING_RESP);
        kNN_CORR_SPOTS(i,I(1:NUM_FEATURES))=1;

        Nearest_Testing_vals=knnclassify(TESTING_FEATURES(:,I(1:NUM_FEATURES)),TRAINING_FEATURES(:,I(1:NUM_FEATURES)),TRAINING_RESP,K,'hamming','consensus');

        CLASS_PERF=classperf(CLASS_PERF,Nearest_Testing_vals,[shuffle_resp(num_resp_include+1:end); shuffle_non_resp(num_non_resp_include+1:end)]);
        kNN_TRAIN_CLASS_CORRECT(i)=CLASS_PERF.LastCorrectRate;

        if DISPLAY_FLAG
            figure(kNN_fig_handle)
            subplot(2,2,4), barh(mean(kNN_CORR_SPOTS(1:i,:),1))
            axis([0 1 0 size(kNN_CORR_SPOTS,2)])
            set(gca,'ytick',CORR_LABELS,'Yticklabel',{'RX Data','SNP Data','Simple ELM','Positional ELM'});
            xlabel('Frequency Picked')
            subplot(2,2,3), hist(kNN_TRAIN_CLASS_CORRECT(1:i),0:0.01:1)
            ylabel('Frequency')
            xlabel('CorrectRate')
        end
    end
    
    if SWR_FLAG
        NAN_MASK=find(sum(isnan(TRAINING_FEATURES))==0);

        intial_inds=FindOptSWR(TRAINING_FEATURES(:,NAN_MASK),TRAINING_RESP);

        [B_values,se,pval,inmodel,stats,nextstep,history]=stepwisefit(TRAINING_FEATURES(:,NAN_MASK),TRAINING_RESP,'inmodel',ismember(NAN_MASK,intial_inds),'display','off');

        VALS=glmval([B_values(inmodel);0],TESTING_FEATURES(:,NAN_MASK(inmodel)),'identity');

        [temp_AUC SWR_SPEC_VALS(i,:)]=CalculateROC(VALS,TESTING_RESP,WANTED_SENS);

        SWR_AUC_VALS(i)=abs(temp_AUC-0.5);
        SWR_REG_VALS(i,NAN_MASK(inmodel))=sign(B_values(inmodel));
        SWR_NORM_VALS(i,NAN_MASK(inmodel))=abs(B_values(inmodel))/sum(abs(B_values(inmodel)));

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
    end
end

if SWR_FLAG
    PRED_VALS=mean(SWR_REG_VALS,1);
    PRED_NORM=mean(SWR_NORM_VALS,1);
    PRED_STD=std(SWR_NORM_VALS,1);
    
    varargout{1}=PRED_VALS;
    varargout{2}=PRED_NORM;
    varargout{3}=PRED_STD;
end

if kNN_FLAG

    varargout{1}=mean(kNN_CORR_SPOTS,1);
end



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
        % Finds the optimal features for step-wise regression using
        % cross validation on the training data.

        
        [NUM_OBSERVATIONS MAX_FEATURES]=size(FEATURES);

        OUTPUT_MAT=zeros(NUM_OPT_REPS,size(SWR_CORR_SPOTS,2));
        AUC_MAP=zeros(NUM_OPT_REPS,1);

        
        if USE_DISTCOMP_FLAG
            parfor (IND = 1:NUM_OPT_REPS)

                groups={'NR','R'};

                [this_train this_test]=crossvalind('holdout',groups(RESP+1),0.1,'classes',groups);

                [B_values,se,pval,inmodel,stats,nextstep,history]=stepwisefit(FEATURES(this_train,:),RESP(this_train),'inmodel',rand(MAX_FEATURES,1)<0.1,'display','off');
                VALS=glmval([B_values(inmodel);0],FEATURES(this_test,inmodel),'identity');

                AUC_MAP(IND)=abs(CalculateROC(VALS,RESP(this_test))-0.5);

                OUTPUT_SLICE=zeros(1,size(SWR_CORR_SPOTS,2));
                OUTPUT_SLICE(NAN_MASK(inmodel))=1;

                OUTPUT_MAT(IND,:)=OUTPUT_SLICE;

            end
        else
            groups={'NR','R'};
            for IND=1:NUM_OPT_REPS

                [this_train this_test]=crossvalind('holdout',groups(RESP+1),0.1,'classes',groups);

                [B_values,se,pval,inmodel,stats,nextstep,history]=stepwisefit(FEATURES(this_train,:),RESP(this_train),'inmodel',rand(MAX_FEATURES,1)<0.1,'display','off');
                VALS=glmval([B_values(inmodel);0],FEATURES(this_test,inmodel),'identity');

                AUC_MAP(IND)=abs(CalculateROC(VALS,RESP(this_test))-0.5);

                OUTPUT_SLICE=zeros(1,size(SWR_CORR_SPOTS,2));
                OUTPUT_SLICE(NAN_MASK(inmodel))=1;

                OUTPUT_MAT(IND,:)=OUTPUT_SLICE;
            end
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