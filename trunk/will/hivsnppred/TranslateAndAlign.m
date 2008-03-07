function TranslateAndAlign(REF_SEQ,TEST_SEQS,PROVIDED_ALIGNMENTS)

load TESTING_DATA HIV_ANNOT_CELL
[INT_NT_MAP ALLOWABLE_CHARS]=AmbigiousNTMap;

STRICT_ALLOWABLE_CHARS='ACGT';


FIG_TOTAL=5;

DISPLAY_FIG_FLAG=false;
NEED_ALIGNMENT_FLAG=true;
NEED_SNP_SPOTS=true;
NEED_TRANSLATE=false;
NEED_BASE_CALLS=true;
OUTPUT_ALL_ALIGMENT=true;
DISTRIBUTED_FLAG=false;
DISTRIBUTED_POOL_FLAG=false;
NEED_PAT_STRUCT=false;
FIGURES_ONLY=false;
FIGURES_TO_GENERATE=1:FIG_TOTAL;

%PROVIDED_ALIGNMENTS={};

WINDOW_CUTOFF=25;
SNP_CUTOFF_MEAN=.4;

SNP_SPOTS=[];
INPUT_PAT_ARRAY=[];

NEED_EXCEL_OUT_FLAG=false;
EXCEL_OUT_FILENAME=[];

%%%%%%%CHECK REF_SEQ%%%%%%%%%
if isstruct(REF_SEQ)&&isfield(REF_SEQ,'Sequence')&&numel(REF_SEQ)==1
    REF_SEQ_CELL={SeqPREP(REF_SEQ.Sequence)};
elseif iscell(REF_SEQ)&&numel(REF_SEQ)==1&&ischar(REF_SEQ{1})
    REF_SEQ_CELL={SeqPREP(REF_SEQ{1})};
elseif ischar(REF_SEQ)&&size(REF_SEQ,1)==1&&size(REF_SEQ,2)>1
    REF_SEQ_CELL={SeqPREP(REF_SEQ)};
else
    error('GenerateSimFigures:REF_SEQ','REF_SEQ must be a cell, a char-array or a struct with a Sequence field');
end

REF_SEQ_LENGTH=length(REF_SEQ_CELL{1});

%%%%%%%CHECK TEST_SEQS%%%%%%%%%
if isstruct(TEST_SEQS)&&isfield(TEST_SEQS,'Sequence')
    TEST_SEQS_CELL=arrayfun(@(x)SeqPREP(x.Sequence),TEST_SEQS,'uniformoutput',false);
elseif iscell(TEST_SEQS)
    TEST_SEQS_CELL=cellfun(@SeqPREP,TEST_SEQS,'uniformoutput',false);
else
    error('GenerateSimFigures:TEST_SEQS','TEST_SEQS must be a cell array or a struct-array with a Sequence field');
end
TEST_SEQS_CELL=TEST_SEQS_CELL(:);

bad_seqs=cellfun('isempty',TEST_SEQS_CELL);

if any(bad_seqs)
    display([int2str(nnz(bad_seqs)) ' Sequences skipped']);
    TEST_SEQS_CELL=TEST_SEQS_CELL(~bad_seqs);
end

NUM_TEST_SEQS=length(TEST_SEQS_CELL);
FRAGMENT_ARRAY=false(NUM_TEST_SEQS,1);
PAT_IND_ARRAY=zeros(NUM_TEST_SEQS,1);



%%%%%%%%%%%%%%%%%%%%%DONE ARGUEMENT CHECKING


%%%%%%%%%%%DO alignments if needed
if NEED_ALIGNMENT_FLAG
    ALIGNMENTS=cell(NUM_TEST_SEQS,1);
    %%%%%use provided alignments if possible
    if ~isempty(PROVIDED_ALIGNMENTS)
        provided_indexes=MatchAlignments(REF_SEQ_CELL,TEST_SEQS_CELL,PROVIDED_ALIGNMENTS,false);

        ALIGNMENTS(provided_indexes~=0)=PROVIDED_ALIGNMENTS(nonzeros(provided_indexes));
    end

    still_needed=cellfun('isempty',ALIGNMENTS);

    if any(still_needed)
        needed_alignments=find(still_needed)';
        if DISTRIBUTED_FLAG
            parfor (i=needed_alignments)
                [junk ALIGNMENTS{i}]=nwalign_mod(REF_SEQ_CELL{1},TEST_SEQS_CELL{i},ALIGNMENT_PROP{:});
            end
        else
            tic
            for i=1:length(needed_alignments)
                [junk ALIGNMENTS{needed_alignments(i)}]=nwalign_mod(REF_SEQ_CELL{1},TEST_SEQS_CELL{needed_alignments(i)},ALIGNMENT_PROP{:});
                time=(toc/i)*(length(needed_alignments)-i);
            end
        end

        %%%%add the alignments preformed so they can be returned easily if
        %%%%desired
        PROVIDED_ALIGNMENTS=[PROVIDED_ALIGNMENTS;ALIGNMENTS(still_needed)];
    end
else
    ALIGNMENTS=PROVIDED_ALIGNMENTS;
end


map=geneticcode;
trans_table=[fieldnames(map) struct2cell(map)];
stop_indexes=find(strcmp('*',trans_table(:,2)));


REF_PROT=cell(1,size(HIV_ANNOT_CELL,1)-2);

guess_spots=[HIV_ANNOT_CELL{2:end-1,1}]';
guess_starts=guess_spots(1:2:end);

guess_lengths=round((guess_spots(2:2:end)-guess_starts)/3);

tic
this_genome=REF_SEQ_CELL{1};
for i=1:length(REF_PROT)
    REF_PROT{i}=Translate(this_genome,guess_starts(i),guess_lengths(i));
end
toc




















%%%%%%%%%%%%%%%%%INPUT PARSING FUNCTIONS



    function AA_seqs=Translate(this_seq,guess_start,guess_length)
        num_seqs=zeros(1,2*guess_length);
        AA_seqs=[];
        test_AA_seqs=[];
        best_length_offset=guess_length;
        if isempty(this_seq)
            warning('GetAlignedSequences:NoSeq',['No sequence '...
                'information found in entry %d'])
        else
            stripped_seq=this_seq(isletter(this_seq));
            offset=-100:10:100;
            offset_counter=1;

            start_spot_log=zeros(size(offset));

            while(offset_counter<length(offset)&&(best_length_offset>=0.1*guess_length))

                start_spot=regexpi(stripped_seq(guess_start+offset(offset_counter):end),map.Starts{1},'once');
                if isempty(start_spot)
                    continue
                end

                start_spot=start_spot+guess_start+offset(offset_counter)-1;
                if ~any(start_spot==start_spot_log)
                    offset_counter=find(staoffset_counter);
                else
                    offset_counter=offset_counter+1;
                    continue
                end

                counter=0;
                while counter==0||trans_table{num_seqs(counter),2}~='*'
                    counter=counter+1;
                    num_seqs(counter)=find(strcmpi(stripped_seq(start_spot:start_spot+2),trans_table(:,1)));
                    start_spot=start_spot+3;
                end
                test_AA_seqs=char(trans_table(num_seqs(1:counter-1),2));
                offset_counter=offset_counter+1;

                length_offset=abs(length(test_AA_seqs)-guess_length);


                if length_offset<best_length_offset
                    AA_seqs=test_AA_seqs;
                    best_length_offset=length_offset;
                end

            end
        end
    end


    function [intmap charorder]=AmbigiousNTMap
        %The rows are in ACGT order and the columns are in charorder
        intmap=zeros(4,16);
        charorder='ACGTRYKMSWBDHVX-';
        %A      %C      %G      %T
        intmap(:,1)=    [1;     0;      0;      0];     %A
        intmap(:,2)=    [0;     1;      0;      0];     %C
        intmap(:,3)=    [0;     0;      1;      0];     %G
        intmap(:,4)=    [0;     0;      0;      1];     %T
        intmap(:,5)=    [0.5;   0;      0.5;    0];     %R
        intmap(:,6)=    [0;     0.5;    0;      0.5];   %Y
        intmap(:,7)=    [0;     0;      0.5;    0.5];   %K
        intmap(:,8)=    [0.5;   0.5;    0;      0];     %M
        intmap(:,9)=    [0;     0.5;    0.5;    0];     %S
        intmap(:,10)=   [0.5;   0;      0;      0.5];   %W
        intmap(:,11)=   [0;     1/3;    1/3;    1/3];   %B
        intmap(:,12)=   [1/3;   0;      1/3;    1/3];   %D
        intmap(:,13)=   [1/3;   1/3;    0;      1/3];   %H
        intmap(:,14)=   [1/3;   1/3;    1/3;    0];     %V
        intmap(:,15)=   [0.25;  0.25;   0.25;   0.25];  %N,X
        intmap(:,16)=   [0;     0;      0;      0];     %-


    end

    function stripped_seq=SeqPREP(input)
        input=upper(input);
        if any(~ismember(input,ALLOWABLE_CHARS))
            warning('MakeSNPCalls:BAD_CHAR','All sequences must be composed of %s. Skipping the incorrect sequences and proceding',ALLOWABLE_CHARS);
            warning('OFF','MakeSNPCalls:BAD_CHAR');
            stripped_seq=[];
        else
            stripped_seq=input(isletter(input));
        end

    end

end
