function [SEQ_ARRAY HASH_VALUES MISMATCH_POS]=ReadAndHashSeqs(INPUT,NUM)
%   ReadAndHashSeqs
%       Reads Sequences from a file and returns the nt2int value and the
%       hash value.  If given an nt2int sequence it will produce a desired
%       number of mismatches and return their nt2int sequences, hash-values
%       and thier mis-match positions
%
%
%       [SEQ_ARRAY HASH_VALUES]=ReadAndHashSeqs(FILE_HANDLE);
%           Reads the next sequence from the file identified by FILE_HANDLE
%           and then returns the nt2int sequence and the hash value
%           associated with the sequence.
%
%
%       [SEQ_ARRAY HASH_VALUES MISMATCH_POS]=ReadAndHashSeqs(INPUT,NUM)
%           Takes the seqeunce given by INPUT and finds the 1:NUM possbile
%           and returns the SEQ_ARRAY, HASH_values and the
%           MISMATCH_POSitions where 0 indicates no-other mismatches.



HASH_FUN=@(input)(prod(cumsum(input.*fliplr(1:length(input)))));


if nargin==2    %Create Mismatches
    MISMATCH_MAT=[2 1 1 1; 3 3 2 2; 4 4 4 3];

    [SEQ_ARRAY MISMATCH_POS]=GenMisMatch(INPUT,0);
    HASH_VALUES=cell2mat(cellfun(@(x)HASH_FUN(x),num2cell(SEQ_ARRAY,2),'uniformoutput',false));

    if NUM>1
        for k=2:NUM


            new_pos=find(sum(MISMATCH_POS==0,2)==0);

            
            seq_array_cell=num2cell(SEQ_ARRAY(new_pos,:),2);
            mismatch_cell=num2cell(MISMATCH_POS(new_pos,:),2);

            [SEQ_ARRAY2 MISMATCH_POS2]=cellfun(@(x,y)(GenMisMatch(x,y)),seq_array_cell,mismatch_cell,'uniformoutput',false);
            SEQ_ARRAY2=cell2mat(SEQ_ARRAY2);
            MISMATCH_POS2=cell2mat(MISMATCH_POS2);

            [SEQS indexes]=unique(SEQ_ARRAY2,'rows','first');
            SEQ_ARRAY=[SEQ_ARRAY; SEQS];

            MISMATCH_POS=[MISMATCH_POS zeros(length(MISMATCH_POS),1); MISMATCH_POS2(indexes,:)];
        end
        [SEQ_ARRAY indexes]=unique(SEQ_ARRAY,'rows','first');
        HASH_VALUES=cell2mat(cellfun(@(x)HASH_FUN(x),num2cell(SEQ_ARRAY,2),'uniformoutput',false));

        MISMATCH_POS=MISMATCH_POS(indexes,:);
    end


else    %Get original sequence from a file
    
    HASH_VALUES=HASH_FUN(SEQ_ARRAY);
    MISMATCH_POS=[];
end




    function [ARRAY POS]=GenMisMatch(input_vec,already_mismatched)

        mis_match_pos=1:length(input_vec);
        mis_match_pos=mis_match_pos(mis_match_pos>already_mismatched(end));

        corr_spots=[1 3*mis_match_pos(1:end-1)+1];

        ARRAY=zeros(3*length(mis_match_pos),length(input_vec));
        POS=zeros(3*length(mis_match_pos),1);

        for i=1:length(mis_match_pos)
            ARRAY(corr_spots(i):corr_spots(i)+2,:)=...
                [ones(3,1)*input_vec(1:mis_match_pos(i)-1)...
                MISMATCH_MAT(:,input_vec(mis_match_pos(i)))...
                ones(3,1)*input_vec(mis_match_pos(i)+1:end)];
            POS(corr_spots(i):corr_spots(i)+2)=ones(3,1)*mis_match_pos(i);
        end
        if nnz(already_mismatched)~=0
            POS=[ones(length(POS),1)*already_mismatched POS];
        end
    end
end
