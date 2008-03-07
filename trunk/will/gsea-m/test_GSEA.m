function test_GSEA(INPUT_DATA,annot_data,SIGNATURES)
%create TEST_DATA ... the rows of the first column of gene_sigs have a TRUE
%difference in means
%
% test_data=rand(1000,100);
% N_genes=1000;
% class_vector=[ones(1,50) zeros(1,50)];
% junk=randperm(1000);
% gene_sigs=[{junk(1:100)}; {junk(51:150)}; {junk(76:175)}; {junk(111:210)}; {junk([1:50 75:100 500:550])}];
% N_genes_sig=cellfun('length',gene_sigs);
% test_data(junk(1:100),:)=repmat(5*(class_vector+1),[100 1]).*rand(100,100);

sig_cell=struct2cell(SIGNATURES)';


sizes=cellfun('length',sig_cell(:,2));

mask=(~cellfun('isempty',regexpi(sig_cell(:,3),'Curated'))&sizes>5&sizes<500);

valid_signatures=sig_cell(mask,:);


uni_gene_symbols=unique(annot_data(:,7));

[tf loc]=ismember(annot_data(:,7),uni_gene_symbols);


tic
COLAPSED_EXP_MAT=zeros(length(uni_gene_symbols),size(INPUT_DATA.GEData,2));

[tf probe_loc]=ismember(annot_data(:,1),INPUT_DATA.ProbeSet);
parfor (i=1:length(uni_gene_symbols))
    COLAPSED_EXP_MAT(i,:)=median(INPUT_DATA.GEData(loc==i,:),1);
end
toc




N_genes=size(test_data,1);
N_genes_sig=cellfun('length',gene_sigs);
p_weight=5.0;

%Rank genes by Singal-Noise ratio.

SNR=Signal_Noise(test_data,class_vector);

[S_SNR I_SNR]=sort(SNR,'descend');

HITS=cellfun(@(x)(ismember(I_SNR,x)),gene_sigs,'uniformoutput',false);

HITS_mat=cat(2,HITS{:});

N_r=sum(abs(S_SNR).^p_weight);


P_val=(abs(repmat(S_SNR,[1 size(HITS_mat,2)]).*HITS_mat).^p_weight)./N_r-(~HITS_mat./repmat(N_genes-N_genes_sig',[size(HITS_mat,1) 1]));
P_stat=cumsum(P_val,1);

[ES leading_val]=max(P_stat,[],1);

NUM_PERMS=100;

ES_null=zeros(NUM_PERMS,length(gene_sigs));

for i=1:NUM_PERMS
    temp=randperm(length(class_vector));
    perm_SNR=Signal_Noise(test_data,class_vector(temp));
    
    [perm_S_SNR perm_I_SNR]=sort(perm_SNR,'descend');
    
    perm_HITS=cellfun(@(x)(ismember(perm_I_SNR,x)),gene_sigs,'uniformoutput',false);
    perm_HITS_mat=cat(2,perm_HITS{:});
    
    perm_P_val=(abs(repmat(perm_S_SNR,[1 size(perm_HITS_mat,2)]).*perm_HITS_mat).^p_weight)./N_r-(~perm_HITS_mat./repmat(N_genes-N_genes_sig',[size(perm_HITS_mat,1) 1]));
    perm_P_stat=cumsum(perm_P_val,1);
    
    ES_null(i,:)=max(perm_P_stat,[],1);
end


pos_mask=ES_null>0;
neg_mask=ES_null<=0;
pos_mean=(sum(ES_null.*pos_mask)./sum(pos_mask));
neg_mean=(sum(ES_null.*neg_mask)./sum(neg_mask));

pos_mean(isnan(pos_mean))=0;
neg_mean(isnan(neg_mean))=0;

NES=(ES.*ES>0)./pos_mean-(ES.*ES<=0)./neg_mean;
NES_null=(ES_null.*pos_mask)./repmat(pos_mean,[NUM_PERMS 1])-(ES_null.*neg_mask)./repmat(neg_mean,[NUM_PERMS 1]);

p=zeros(1,length(gene_sigs));
FDR=zeros(1,length(gene_sigs));
for i=1:length(gene_sigs)
    if ES(i)>0     %consider the distribution with the same sign
        p(i)=sum(ES_null(:,i)>=ES(i))/NUM_PERMS;
        FDR(i)=sum(NES_null(:,i)>=NES(i))/sum(NES_null(:,i)>=0);
    else
        p(i)=sum(ES_null(:,i)<=ES(i))/NUM_PERMS;
        FDR(i)=sum(NES_null(:,i)<=NES(i))/sum(NES_null(:,i)<=0);
    end
end








    function SN_val=Signal_Noise(data,class)
        SN_val=(mean(data(:,true&class),2)-mean(data(:,~class),2))./(std(data(:,true&class),[],2)+std(data(:,~class),[],2));
    end
end



