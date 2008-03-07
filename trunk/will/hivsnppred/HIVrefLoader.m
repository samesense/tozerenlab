function HIV_REF=HIVrefLoader(FILENAME)
%   HIVrefLoader
%       HIVrefLoader reads a FASTA formatted file which has the complete
%       HIV genomic region followed by the the NT sequence of each protien
%       region.  This produces and HIV_REF which is compatible with all of
%       the GenerateRead functions.
%
%   HIV_REF=HIVrefLoader(FILENAME)
%
%

[headers seqs]=fastaread(FILENAME);

gene_names=cell(length(headers)-1,1);
locs=zeros(length(headers)-1,2);
for i=2:length(headers)
    temp=textscan(headers{i},'gi|4558520:%n-%n %s HIV-1, complete genome');
    locs(i-1,:)=[temp{1:2}];
    gene_names(i-1)=temp{3};
end

trans_aa=nt2aa(seqs(2:end));

HIV_REF=struct('Sequence',seqs{1},'GeneNames',{gene_names},'AAseqs',{trans_aa'},'GenePos',locs,'AppendedGenome',[trans_aa{:}],'AAPos',cumsum([1 cellfun('length',trans_aa(1:end-1))]));