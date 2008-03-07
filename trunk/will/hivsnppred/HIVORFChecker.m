function [OUTPUT_SEQS MAPPING]=HIVORFChecker(HIV_REF,INPUT_SEQ,varargin)
%   HIVORFChecker
%       Determines which of the input sequences (assuming each is from a
%       different ORF) corresponds to which HIV protein.  It does this by
%       aligning each INPUT sequence to each HIV protein sequence and then
%       returning the index of maximum homology.
%
%
%   [TRANS_SEQS MAPPING]=HIVORFChecker(HIV_REF,INPUT_SEQ)
%
%       TRANS_SEQS  The translated Amino-Acid sequences which match HIV
%                   proteins.
%
%       MAPPING     The mapping to the REFERENCE AMINO ACID
%
%   HIVORFChecker(...,'WHOLE_GENOME',true)
%
%       Use this when you know that the sequence is an ENTIRE genome and
%       this function will return cell-arrays which correspond with each
%       protein.
%
%   HIVORFChecker(...,'need_translate',false)
%
%       Use this when you are providing translated sequence.
%
%
%
%



WHOLE_GENOME_FLAG=false;
NEED_TRANSLATE=true;
if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case 'whole_genome'
                WHOLE_GENOME_FLAG=true&varargin{i+1};
            case 'need_translate'
                NEED_TRANSLATE=true&varargin{i+1};
            otherwise
                error('HIVORFChecker:BAD_ARG','An unknown arguement was provided: %s',lower(varargin))
        end

    end
end

if NEED_TRANSLATE
    %%%%If all of the sequences are ACGT then we can use built-in matlab
    %%%%functions ... otherwise we have to translate ourselves.
    if iscell(INPUT_SEQ)&&length(INPUT_SEQ)==1
        EASY_FLAG=all(ismember(upper(INPUT_SEQ{1}),'ACGT'));
        SEQ=INPUT_SEQ{1};
    elseif ischar(INPUT_SEQ)
        EASY_FLAG=all(ismember(upper(INPUT_SEQ),'ACGT'));
        SEQ=INPUT_SEQ;
    end

    map=geneticcode;
    trans_table=[[fieldnames(map) struct2cell(map)];{'NNN'  'X'}];

    start_spots=regexpi(SEQ,map.Starts{1},'start')';
    if EASY_FLAG    %USE built-in function if possible
        TRANS_SEQS=arrayfun(@(x)(nt2aa(SEQ(x:end))),start_spots,'uniformoutput',false);
        TRANS_SEQS=cellfun(@(x)(x(1:min(find(x=='*',1)-1,length(x)))),TRANS_SEQS,'uniformoutput',false);
    else
        TRANS_SEQS=cell(length(start_spots),1);
        for k=1:length(start_spots)
            start_spot=start_spots(k);
            num_seqs=zeros(1,length(SEQ)-start_spot);
            counter=0;
            while start_spot+2<length(SEQ)&&(counter==0||trans_table{num_seqs(counter),2}~='*')
                counter=counter+1;
                this_spot=find(strcmpi(SEQ(start_spot:start_spot+2),trans_table(:,1)));
                if ~isempty(this_spot)
                    num_seqs(counter)=this_spot;
                else
                    num_seqs(counter)=length(trans_table);
                end
                start_spot=start_spot+3;
            end
            TRANS_SEQS{k}=char(trans_table(num_seqs(1:counter-1),2))';
        end
    end
else
    TRANS_SEQS=INPUT_SEQ;
end

scores=zeros(length(HIV_REF.AAseqs),length(TRANS_SEQS));
alignments=cell(length(HIV_REF.AAseqs),length(TRANS_SEQS));

for i=1:length(HIV_REF.AAseqs)
    for j=1:length(TRANS_SEQS)
        if WHOLE_GENOME_FLAG||(2*length(HIV_REF.AAseqs{i})>length(TRANS_SEQS{j})&&0.1*length(HIV_REF.AAseqs{i})<length(TRANS_SEQS{j}))
            [scores(i,j) alignments{i,j}]=nwalign_mod(HIV_REF.AAseqs{i},TRANS_SEQS{j});
        end
    end
end

%deal with case where none of the input sequences match
if all(scores(:)==0)
    OUTPUT_SEQS=[];
    MAPPING=[];
    warning('HIVORFChecker:NO_SEQ','None of the ORFs corresponded to a HIV protein.')
    return
end


if ~WHOLE_GENOME_FLAG
    [Y I]=sort(scores(:),'descend');
    [HIV_spot ORF_spot]=ind2sub(size(scores),I(find(Y~=0,1)));

    OUTPUT_SEQS=TRANS_SEQS{ORF_spot};
    temp_ref=cumsum([HIV_REF.AAPos(HIV_spot) isletter(alignments{HIV_spot,ORF_spot}(1,:))]);
    temp_clinical=cumsum([1 isletter(alignments{HIV_spot,ORF_spot}(3,:))]);
    MAPPING=temp_ref(temp_clinical);
else

    OUTPUT_SEQS=cell(length(HIV_REF.AAseqs),1);
    MAPPING=cell(length(HIV_REF.AAseqs),1);

    for i=1:length(HIV_REF.AAseqs)
        [vals orf_ind]=max(scores,[],2);
        [val protein_ind]=max(vals);

        OUTPUT_SEQS(protein_ind)=TRANS_SEQS(orf_ind(protein_ind));

        temp_ref=cumsum([HIV_REF.AAPos(protein_ind) isletter(alignments{protein_ind,orf_ind(protein_ind)}(1,:))]);
        temp_clinical=cumsum([1 isletter(alignments{protein_ind,orf_ind(protein_ind)}(3,:))]);

        MAPPING{protein_ind}=temp_ref(temp_clinical);

        scores(:,orf_ind(protein_ind))=NaN;
        scores(protein_ind,:)=NaN;

    end

end

end

