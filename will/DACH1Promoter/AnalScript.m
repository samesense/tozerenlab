KNOWN_SEQ='AA[AT]ACAA[AT]TAA[AT]';


[headers PromoterSeqs]=fastaread('unique_IDS.fa');
LLIDS=cell(size(headers));
GeneNames=cell(size(headers));
for i=1:length(headers)
    temp=textscan(headers{i},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
    GeneNames(i)=temp{1};
    LLIDS(i)=temp{2};
end

display('Finding Matches')
[GeneSymbol LLID POS STRAND ORIENT SEQS]=DNAPromoterMatcher(KNOWN_SEQ,'unique_IDS.fa',2);

SeqProf=seqprofile(SEQS,'alphabet','nt');

%[normal, reverse, complement, reverse-complement]
SeqPos=zeros(length(PromoterSeqs),4);
SeqMaxVal=zeros(length(PromoterSeqs),4);

display('Calculating normal pos')
SeqVals=PWMEvaluator(SeqProf,PromoterSeqs);
SeqVals_mat=cell2mat(SeqVals);

[SeqMaxVal(:,1) SeqPos(:,1)]=max(SeqVals_mat,[],2);

display('Calculating seqreverse')
SeqVals=PWMEvaluator(SeqProf,seqreverse(PromoterSeqs));
SeqVals_mat=cell2mat(SeqVals);

[SeqMaxVal(:,2) SeqPos(:,2)]=max(SeqVals_mat,[],2);

display('Calculating complement')
PromoterComplement=cellfun(@(x)(seqcomplement(x)),PromoterSeqs,'uniformoutput',false);
SeqVals=PWMEvaluator(SeqProf,PromoterComplement);
SeqVals_mat=cell2mat(SeqVals);

[SeqMaxVal(:,3) SeqPos(:,3)]=max(SeqVals_mat,[],2);
clear PromoterComplement

display('Calculating reverse-complement')
PromoterReverseComplement=cellfun(@(x)(seqrcomplement(x)),PromoterSeqs,'uniformoutput',false);
SeqVals=PWMEvaluator(SeqProf,PromoterReverseComplement);
SeqVals_mat=cell2mat(SeqVals);

[SeqMaxVal(:,4) SeqPos(:,4)]=max(SeqVals_mat,[],2);

display('loading RNA Null-data')
[junk_headers rna_seq]=fastaread('rna.fa');
shuffle=randperm(length(rna_seq));

display('Calculating Null-distribution')
RNA_NullVals=PWMEvaluator(SeqProf,rna_seq(shuffle(1:10000)));
RNA_NullDist=[RNA_NullVals{:}];


RNA_NullDist=sort(RNA_NullDist(~isnan(RNA_NullDist)&RNA_NullDist~=0));

clear rna_seq junk_headers RNA_NullVals
% % 
% display('loading pseudo_without_product Null-data')
% [junk_headers pseudo_seq]=fastaread('pseudo_without_product.fa');
% 
% display('Calculating pseudo_seq Null-distribution')
% pseudo_NullVals=PWMEvaluator(SeqProf,pseudo_seq);
% pseudo_NullDist=[pseudo_NullVals{:}];
% 
% pseudo_NullDist=sort(pseudo_NullDist(~isnan(pseudo_NullDist)&pseudo_NullDist~=0));
% 
% clear pseudo_seq pseudo_NullVals pseudo_NullVals
% 
% 
% display('loading CHR1 Null-data')
% [junk_headers ref_chr1_seq]=fastaread('ref_chr1.fa');
% 
% display('Calculating ref_chr1 Null-distribution')
% ref_chr1_NullVals=PWMEvaluator(SeqProf,ref_chr1_seq);
% ref_chr1_NullDist=[ref_chr1_NullVals{:}];
% 
% ref_chr1_NullDist=sort(ref_chr1_NullDist(~isnan(ref_chr1_NullDist)&ref_chr1_NullDist~=0));
% 
% clear ref_chr1_seq ref_chr1_NullVals ref_chr1_NullVals



display('Calculating P-vals')
RNA_Pvals=zeros(size(SeqPos));
psuedo_Pvals=zeros(size(SeqPos));
ref_chr1_Pvals=zeros(size(SeqPos));
tic
for pos=1:numel(SeqPos)
    spot=find(RNA_NullDist>SeqMaxVal(pos),1);
    if ~isempty(spot)
        RNA_Pvals(pos)=1-spot/length(RNA_NullDist);
    else
        RNA_Pvals(pos)=1/length(RNA_NullDist);
    end
    
%     spot=find(pseudo_NullDist>SeqMaxVal(pos),1);
%     if ~isempty(spot)
%         psuedo_Pvals(pos)=1-spot/length(pseudo_NullDist);
%     else
%         psuedo_Pvals(pos)=0;
%     end
%     
%     spot=find(ref_chr1_NullDist>SeqMaxVal(pos),1);
%     if ~isempty(spot)
%         ref_chr1_Pvals(pos)=1-spot/length(ref_chr1_NullDist);
%     else
%         ref_chr1_Pvals(pos)=1/length(ref_chr1_NullDist);
%     end

    
    fprintf('Numel: %d of %d done in: %f min\n',pos,numel(RNA_Pvals),(toc/pos)*(numel(RNA_Pvals)-pos)/60)
end


% display('Setting-up Output')
% %{Name, LLID, POS, STRAND, Orient, Seq, p_val}
% [I J]=find(Pvals<P_val_cutoff);
% excel_output=cell(nnz(I),7);
% 
% 
% for i=1:length(I)
%     excel_output(i,1)=GeneNames(I(i));
%     excel_output(i,2)=LLIDS(I(i));
%     excel_output{i,7}=Pvals(I(i),J(i));
%     switch I(i)
%         case 1
%             excel_output{i,3}=SeqPos(I(i),J(i))-4000;
%             excel_output{i,4}='3/prime to 5/prime';
%             excel_output{i,5}='Sense';
%             temp=PromoterSeqs{J(i)};
%             excel_output{i,6}=temp(SeqPos(I(i),J(i))+length(KNOWN_SEQ));
%             
%         case 2
%             excel_output{i,3}=-SeqPos(I(i),J(i));
%             excel_output{i,4}='5/prime to 3/prime';
%             excel_output{i,4}='Sense';
%             temp=seqreverse(PromoterSeqs{J(i)});
%             excel_output{i,6}=temp(SeqPos(I(i),J(i))+length(KNOWN_SEQ));
%             
%         case 3
%             excel_output{i,3}=SeqPos(I(i),J(i))-4000;
%             excel_output{i,4}='3/prime to 5/prime';
%             excel_output{i,4}='Anti-Sense';
%             temp=seqcomplement(PromoterSeqs{J(i)});
%             excel_output{i,6}=temp(SeqPos(I(i),J(i))+length(KNOWN_SEQ));
%             
%         case 4
%             excel_output{i,3}=-SeqPos(I(i),J(i));
%             excel_output{i,4}='5/prime to 3/prime';
%             excel_output{i,4}='Anti-Sense';
%             temp=seqrcomplement(PromoterSeqs{J(i)});
%             excel_output{i,6}=temp(SeqPos(I(i),J(i))+length(KNOWN_SEQ));
%     end
% end
% 
% xlswrite('PvalOutput',excel_output);

