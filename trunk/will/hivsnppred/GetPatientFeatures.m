function varargout=GetPatientFeatures(PAT_STRUCT,ELM_STRUCT,varargin)
%   GetPatientFeatures
%       A helper function to retrieve the features for classification of
%       patients.
%
%   [FEATURES RESPONSE PAT_INDEX FEATURE_NAMES FEATURE_INDS]= ...
%               GetPatientFeatures(PAT_STRUCT,feature1,feature2, featureN)
%
%   PAT_STRUCT      The patient structure to extract features from.
%
%   ELM_STRUCT      A structure describing the ELMs (as created by
%                   ELMDownloadParser).  Can be left as [] if you don't
%                   care about FEATURE_NAMES.
%
%   FEATURES        The feature-vectors for each patient.
%
%   RESPONSE        The response variable for each patient.
%
%   PAT_INDEX       An index vector associating each row in FEATURES to a
%                   patient in PAT_STRUCT.
%
%   FEATURE_NAMES   A cell-array giving a descriptive name for each
%                   feature.
%
%   FEATURE_INDS    The indexes to the first feature for each type
%                   provided.
%
%
%   Possible Feature Choices
%
%   DrugRegimine    A binary indicator vector indicating the presence or
%                   abscence of each Drug during the therapy.
%
%   BaseCalls       A binary indicator vector for the SNP identity at spots
%                   determined by MakeSNPCalls
%
%   SimpleELM       A binary indicatory vector for the presence or absence
%                   of each ELM in the sequence.
%
%   PositionalELM   A binary indicator vector for the presence or absence
%                   of ELMs at specific locations.
%
%   PositionalPWM   A continious valued vector indicating the conservation
%                   of the ELM relative to the reference genomic sequences.
%
%


if isempty(ELM_STRUCT)
    %create a dummy variable incase needed
    ELM_STRUCT=struct('Name',arrayfun(@(x)(['ELM ' int2str(x)]),1:500,'uniformoutput',false));
end


feature_cell=cell(length(varargin),1);
name_cell=cell(length(varargin),1);
index_cell=cell(length(varargin),1);

for i=1:length(varargin)
    switch(lower(varargin{i}))
    
        case 'drugregimine'
            [temp_rx temp_rx_inds]=PatientStructHelper(PAT_STRUCT,{'RX_vals','explodenumeric'});
            
            mask=temp_rx(:,1)==0;
            
            feature_cell{i}=temp_rx(mask,2:end);
            index_cell{i}=temp_rx_inds(mask);
            
            name_cell{i}=PAT_STRUCT{1}.RX_names';
            
        case 'basecalls'
            [temp_mat index_cell{i} ...
                 snp_spots junk_ind]=PatientStructHelper(PAT_STRUCT,...
                 {'BASE_CALLS','explodenumeric'},{'SNP_SPOTS','explodeNumeric'});
        
            BASE_CALLS=zeros(size(temp_mat,1),size(temp_mat,2)*5);
            m=size(BASE_CALLS,2);
            translate_vars=[1 2 3 4 16];
            translate_mat=[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; NaN NaN NaN NaN NaN];
            
            %Explode BASE_CALLS into a feature vector
            for k=1:size(temp_mat,1)
                [TF LOC]=ismember(temp_mat(k,:),translate_vars);
                LOC(~TF)=6;
                BASE_CALLS(k,:)=reshape(translate_mat(LOC,:),[],m);
            end
            
            feature_cell{i}=BASE_CALLS;
            
            snp_spots=snp_spots(1,:);
            temp_names=cell(1,length(snp_spots)*5);
            counter=1;
            
            for k=1:length(snp_spots)
                temp_names{counter}=[int2str(snp_spots(k)) ' A'];
                temp_names{counter+1}=[int2str(snp_spots(k)) ' C'];
                temp_names{counter+2}=[int2str(snp_spots(k)) ' G'];
                temp_names{counter+3}=[int2str(snp_spots(k)) ' T'];
                temp_names{counter+4}=[int2str(snp_spots(k)) ' GAP'];

                counter=counter+5;
            end
            
            name_cell{i}=temp_names;
            
            clear temp_names;
            
        case 'simpleelm'
            [feature_cell{i} index_cell{i}]=PatientStructHelper(PAT_STRUCT,{'ELM_simple','explodeNumeric'});
            [temp_names{1:length(ELM_STRUCT)}]=deal(ELM_STRUCT.Name);

            name_cell{i}=temp_names;
            clear temp_names;
            
        case 'positionalelm'
            [feature_cell{i} index_cell{i}]=PatientStructHelper(PAT_STRUCT,{'ELM_vec','explodeNumeric'});
            elm_annot=PAT_STRUCT{1}.ELM_annot;
            temp_names=cell(1,size(elm_annot,2));

            for k=1:size(temp_names,2)
                temp_names{k}=[ELM_STRUCT(elm_annot(1,k)).Name ' at [' int2str(elm_annot(2,k)) ',' int2str(elm_annot(3,k)) ']'];
            end

            name_cell{i}=temp_names;
            clear temp_names;

        case 'positionalpwm'
            [feature_cell{i} index_cell{i}]=PatientStructHelper(PAT_STRUCT,{'ELM_PWM','explodeNumeric'});
            elm_annot=PAT_STRUCT{1}.ELM_PWMannot;
            temp_names=cell(1,size(elm_annot,2));

            for k=1:size(temp_names,2)
                temp_names{k}=[ELM_STRUCT(elm_annot(1,k)).Name ' PWM at [' int2str(elm_annot(2,k)) ',' int2str(elm_annot(3,k)) ']'];
            end

            name_cell{i}=temp_names;
            clear temp_names;
            
            
            
        case 'responder'
            %do nothing, responder is put in the last cell automatically.
            
        otherwise
            error('GetPatientFeatures:BAD_ARGUEMENT', 'An unknown arguement was provided: %s',varargin{i});
    end
end

[resp_var resp_ind]=PatientStructHelper(PAT_STRUCT,{'IS_RESPONDER','leaveNumeric'});

PAT_INDEX=resp_ind;
%find the patient indexes that are in all of the features
for i=1:length(index_cell)
    PAT_INDEX=intersect(PAT_INDEX,index_cell{i});
end

num_feat=cumsum([0; cellfun('size',feature_cell,2)]);

FEATURES=zeros(length(PAT_INDEX),num_feat(end));

for i=1:length(num_feat)-1
    [TF LOC]=ismember(PAT_INDEX,index_cell{i});
    
    FEATURES(:,num_feat(i)+1:num_feat(i+1))=feature_cell{i}(LOC,:);
end

[TF LOC]=ismember(PAT_INDEX,resp_ind);

RESPONSE=resp_var(LOC);

varargout=cell(nargout,1);
varargout{1}=FEATURES;
varargout{2}=RESPONSE;

if nargout>2
    varargout{3}=PAT_INDEX;
end

if nargout>3
    varargout{4}=cat(2,name_cell{:});
end

if nargout>4
    varargout{5}=num_feat;
end
    
























