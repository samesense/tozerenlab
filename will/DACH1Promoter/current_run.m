search_seq ='AATACAATTAAAT';

M=length(search_seq);

[header seqs]=fastaread('unique_IDS.fa');
shorter_seqs=cellfun(@(x)(x(2000:end)),seqs,'uniformoutput',false);
[shorter_seqs I J]=unique(shorter_seqs);

header=header(I);


trans_mat=[eye(M);eye(M)*2;eye(M)*3;eye(M)*4];

final_mat=[];
for i=1:M*4
    final_mat=[final_mat;repmat(trans_mat(i,:),[M*4 1])+[diag(~trans_mat(i,:),0);diag(~trans_mat(i,:),0)*2;diag(~trans_mat(i,:),0)*3;diag(~trans_mat(i,:),0)*4]];
end

input_seq=nt2int(search_seq);
for i=1:M
    final_mat=final_mat(final_mat(:,i)~=input_seq(i),:);
end

match_mat=zeros(length(shorter_seqs),length(final_mat),4);
for i=1:length(final_mat)
    this_seq=input_seq;
    this_seq(final_mat(i,:)~=0)=nonzeros(final_mat(i,:));
    this_seq=int2nt(this_seq);
    match_mat(:,i,1)=FASTBoyer_Moore({this_seq},shorter_seqs');
    match_mat(:,i,2)=FASTBoyer_Moore({seqreverse(this_seq)},shorter_seqs');
    match_mat(:,i,3)=FASTBoyer_Moore({seqcomplement(this_seq)},shorter_seqs');
    match_mat(:,i,4)=FASTBoyer_Moore({seqrcomplement(this_seq)},shorter_seqs');
    display(i/length(final_mat))
end

geneinfo=textscan([header{10}],'%*s%*s%s%s%*s%*s%*s','delimiter','|');
headers=geneinfo{1};
LLIDS=geneinfo{2};

mask=sum(final_mat~=0,2)==1;
spots=find(match_mat);
[I J K]=ind2sub(size(match_mat),spots);

single_mask=ismember(J,find(mask));

I=I(single_mask);
J=J(single_mask);
K=K(single_mask);


excel_output=cell(length(J),5);
for i=1:length(J)
    temp=textscan(header{I(i)},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
    excel_output(i,1)=temp{1};
    excel_output(i,2)=temp{2};
    excel_output{i,3}=2000-match_mat(I(i),J(i),K(i));
    switch (K(i))
        case 1
            excel_output{i,4}='Sense';
            excel_output{i,5}='3prime - 5prime';
        case 2
            excel_output{i,4}='Sense';
            excel_output{i,5}='5prime - 3prime';
        case 3
            excel_output{i,4}='Anti-Sense';
            excel_output{i,5}='3prime - 5prime';
        case 4
            excel_output{i,4}='Anti-Sense';
            excel_output{i,5}='5prime - 3prime';
    end
end

xlswrite('TFtargets.xls',[{'Gene Symbol','EntrezID','BP Upstream','Strand','Orientation'};excel_output],1)

exact_mat=zeros(length(shorter_seqs),4);
this_seq=search_seq;
exact_mat(:,1)=FASTBoyer_Moore({this_seq},shorter_seqs');
exact_mat(:,2)=FASTBoyer_Moore({seqreverse(this_seq)},shorter_seqs');
exact_mat(:,3)=FASTBoyer_Moore({seqcomplement(this_seq)},shorter_seqs');
exact_mat(:,4)=FASTBoyer_Moore({seqrcomplement(this_seq)},shorter_seqs');

[I J]=find(exact_mat);


excel_output=cell(length(I),5);
for i=1:length(I)
    temp=textscan(header{I(i)},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
    excel_output(i,1)=temp{1};
    excel_output(i,2)=temp{2};
    excel_output{i,3}=2000-exact_mat(I(i),J(i));
    switch (J(i))
        case 1
            excel_output{i,4}='Sense';
            excel_output{i,5}='3prime - 5prime';
        case 2
            excel_output{i,4}='Sense';
            excel_output{i,5}='5prime - 3prime';
        case 3
            excel_output{i,4}='Anti-Sense';
            excel_output{i,5}='3prime - 5prime';
        case 4
            excel_output{i,4}='Anti-Sense';
            excel_output{i,5}='5prime - 3prime';
    end
end

xlswrite('TFtargets.xls',[{'Gene Symbol','EntrezID','BP Upstream','Strand','Orientation'};excel_output],2)












