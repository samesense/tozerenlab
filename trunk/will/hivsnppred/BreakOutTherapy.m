function BreakOutTherapy(PAT_STRUCT,vargin)


if isempty(vargin)  %want histogram
    [RX_vals RX_inds]=PatientStructHelper(PAT_STRUCT,{'RX_vals','explodeNumeric'});
    RX_names=PAT_STRUCT{1}.RX_names;
    uni_therapies=unique(RX_vals(:,2:20),'rows');

    [num_entries junk]=size(RX_vals);
    [num_therapies num_drugs]=size(uni_therapies);
    
    before_study_mask=RX_vals(:,1)<0;
    during_study_mask=RX_vals(:,1)>=0&RX_vals(:,1)<=24;
    
    study_hist=zeros(num_therapies,2);
    
    for i=1:num_therapies
        mask=sum(repmat(uni_therapies(i,:),[num_entries 1])==RX_vals(:,2:20),2)==num_drugs;
        study_hist(i,:)=[nnz(mask&before_study_mask) nnz(mask&during_study_mask)];
    end
    
    common_therapies=find(study_hist(:,2)>100)';
    barh(study_hist(common_therapies,:))
    common_therapy_names=cell(nnz(common_therapies),1);
    
    for i=1:length(common_therapies)
        common_therapy_names{i}={cell2mat(RX_names([uni_therapies(common_therapies(i),:) zeros(1,num_drugs-size(common_therapies,2))]&true))};
    end
    
    
    
    
    set(gca,'ytick',1:nnz(common_therapies))
    
    
    
    



end







































end