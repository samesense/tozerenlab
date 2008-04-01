% KNOWN_SEQ='AATACAATTAAAT';
% 
% 
% [headers PromoterSeqs]=fastaread('unique_IDS.fa');
% LLIDS=cell(size(headers));
% GeneNames=cell(size(headers));
% for i=1:length(headers)
%     temp=textscan(headers{i},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
%     GeneNames(i)=temp{1};
%     LLIDS(i)=temp{2};
% end
% 
% display('Finding Matches')
% [GeneSymbol LLID POS STRAND ORIENT SEQS]=DNAPromoterMatcher(KNOWN_SEQ,'unique_IDS.fa',2);
% 
% SeqProf=seqprofile(SEQS,'alphabet','nt');
% 
% %[normal, reverse, complement, reverse-complement]
% SeqPos=zeros(length(PromoterSeqs),4);
% SeqMaxVal=zeros(length(PromoterSeqs),4);
% 
% display('Calculating normal pos')
% SeqVals=PWMEvaluator(SeqProf,PromoterSeqs);
% SeqVals_mat=cell2mat(SeqVals);
% 
% [SeqMaxVal(:,1) SeqPos(:,1)]=max(SeqVals_mat,[],2);
% 
% display('Calculating seqreverse')
% SeqVals=PWMEvaluator(SeqProf,seqreverse(PromoterSeqs));
% SeqVals_mat=cell2mat(SeqVals);
% 
% [SeqMaxVal(:,2) SeqPos(:,2)]=max(SeqVals_mat,[],2);
% 
% display('Calculating complement')
% PromoterComplement=cellfun(@(x)(seqcomplement(x)),PromoterSeqs,'uniformoutput',false);
% SeqVals=PWMEvaluator(SeqProf,PromoterComplement);
% SeqVals_mat=cell2mat(SeqVals);
% 
% [SeqMaxVal(:,3) SeqPos(:,3)]=max(SeqVals_mat,[],2);
% clear PromoterComplement
% 
% display('Calculating reverse-complement')
% PromoterReverseComplement=cellfun(@(x)(seqrcomplement(x)),PromoterSeqs,'uniformoutput',false);
% SeqVals=PWMEvaluator(SeqProf,PromoterReverseComplement);
% SeqVals_mat=cell2mat(SeqVals);
% 
% [SeqMaxVal(:,4) SeqPos(:,4)]=max(SeqVals_mat,[],2);
% 
% display('loading Null-data')
% [junk_headers rna_seq]=fastaread('rna.fa');
% 
% shuffle=randperm(length(rna_seq));
% bkg_seq=rna_seq(shuffle(1:10000));
% 
% clear rna_seq junk_headers
% 
% display('Calculating Null-distribution')
% NullVals=PWMEvaluator(SeqProf,bkg_seq);
% NullDist=sort([NullVals{:}]);
% 
% P_val_cutoff=0.05;
% cutoff_val=NullDist(round(length(NullDist)*(1-P_val_cutoff)));

display('Calculating P-vals')
Pvals=zeros(size(SeqPos));
for pos=1:4
    for i=1:length(Pvals)
        if SeqMaxVal(i,pos)>cutoff_val
            spot=find(NullDist>SeqMaxVal(i,pos),1);
            if ~isempty(spot)
                Pvals(i,pos)=1-spot/length(NullDist);
            else
                Pvals(i,pos)=NaN;
            end
        else
            Pvals(i,pos)=1;
        end
    end
end


display('Setting-up Output')
%{Name, LLID, POS, STRAND, Orient, Seq, p_val}
[I J]=find(Pvals<P_val_cutoff);
excel_output=cell(nnz(I),7);


for i=1:length(I)
    excel_output(i,1)=GeneNames(I(i));
    excel_output(i,2)=LLIDS(I(i));
    excel_output{i,7}=Pvals(I(i),J(i));
    switch I(i)
        case 1
            excel_output{i,3}=SeqPos(I(i),J(i))-4000;
            excel_output{i,4}='3/prime to 5/prime';
            excel_output{i,5}='Sense';
            temp=PromoterSeqs{J(i)};
            excel_output{i,6}=temp(SeqPos(I(i),J(i))+length(KNOWN_SEQ));
            
        case 2
            excel_output{i,3}=-SeqPos(I(i),J(i));
            excel_output{i,4}='5/prime to 3/prime';
            excel_output{i,4}='Sense';
            temp=seqreverse(PromoterSeqs{J(i)});
            excel_output{i,6}=temp(SeqPos(I(i),J(i))+length(KNOWN_SEQ));
            
        case 3
            excel_output{i,3}=SeqPos(I(i),J(i))-4000;
            excel_output{i,4}='3/prime to 5/prime';
            excel_output{i,4}='Anti-Sense';
            temp=seqcomplement(PromoterSeqs{J(i)});
            excel_output{i,6}=temp(SeqPos(I(i),J(i))+length(KNOWN_SEQ));
            
        case 4
            excel_output{i,3}=-SeqPos(I(i),J(i));
            excel_output{i,4}='5/prime to 3/prime';
            excel_output{i,4}='Anti-Sense';
            temp=seqrcomplement(PromoterSeqs{J(i)});
            excel_output{i,6}=temp(SeqPos(I(i),J(i))+length(KNOWN_SEQ));
    end
end

xlswrite('PvalOutput',excel_output);

