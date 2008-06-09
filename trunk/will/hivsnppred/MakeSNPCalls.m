function [varargout]=MakeSNPCalls(REF_SEQ,TEST_SEQS,varargin)
%   MakeSNPCalls
%       Aligns all TEST sequences to the REFERENCE sequence to find
%       conserved islands of sequence and SNPs within those islands.
%
%   [SNP_SPOTS BASE_CALLS]=MakeSNPCalls(REF_SEQ,TEST_SEQS)
%
%       REF_SEQ         A defined reference sequence to compare against
%       TEST_SEQS       A set of test sequences (Assumed to be whole genome
%                       sequences)
%                           See FRAGMENTS
%
%       BASE_CALLS      A matrix with the base calls at each spot defined
%                       as a SNP by the provided criteria
%
%    [SNP_SPOTS BASE_CALLS ALIGNMENTS]=MakeSNPCalls(REF_SEQ,TEST_SEQS)
%
%       ALIGNMENTS      A Nx1 cell of the alignments between each test
%                       sample and the reference.  These alignments can be
%                       out of order and the program will determine which
%                       (if any) need to be aligned (see MatchAlignments).
%
%
%    [SNP_SPOTS BASE_CALLS ALIGNMENTS MUT_CALLS]=...
%                                       MakeSNPCalls(REF_SEQ,TEST_SEQS)
%
%
%       MUT_CALLS       A matrix the same size as BASE_CALLS which
%                       describes whether sample creates an Amino Acid
%                       mutation.
%                           MUT_CALLS(I,J)==0
%                               The Jth SNP in the Ith sample is a 'Silent
%                               Mutation'.
%                           MUT_CALLS(I,J)==1
%                               The Jth SNP in the Ith sample created a
%                               'Mis-sense Mutation' (Amino Acid Change)
%                           MUT_CALLS(I,J)==2
%                               The Jth SNP in the Ith sample created a
%                               'Non-sense Mutation' (Frame-shift with a
%                               stop codon)
%
%
%
%
%   Creating A Patient Structure
%
%   [PAT_STRUCT ...]=MakeSNPCalls(REF_SEQ,TEST_SEQS,'PAT_STRUCT_OUT',true,...)
%
%       PAT_STRUCT      A structure with the information determined from
%                       MakeSNPCalls and associated with a PATIENT.  This
%                       can be used by other programs to display the data.
%                       This program will add or update the following
%                       fields.
%
%                           Sequence_CELL
%                           Alignment_CELL
%                           SNP_POS
%                           BASE_CALLS
%                           MUT_CALLS
%                           IS_SNP_VEC
%                           IS_FRAGMENT
%                       Any subsequent output arguements are appended as
%                       desbribed above.
%
%   Optional Properties
%
%       ALIGNMENTS      ALIGNMENT_CELL
%                           A Nx1 cell of alignments to the reference
%                           sample
%       ALIGNMENT_PROP  A cell array of properties to pass to NWALIGN
%
%       TRUST_ORDER     true|[false]
%                           Will blindly assume that ALIGNMENT_CELL is in
%                           the proper order and will not preform any error
%                           checking (USE WITH CARE).
%
%       ORDERED_ALIGN_OUTPUT    true|[false]
%                                   Will output ALIGNMENTS such that it
%                                   will satisfy TRUST_ORDER in any
%                                   subsequent calls.
%
%
%       FIG_ONLY        [true]|false        [1:5]
%                       This arguement can be used to control which (if any
%                       figures) to generate.  If the arguement is TRUE
%                       then all figures will be created while FALSE
%                       creates no figures.  You can also enter a vector
%                       indicating which figures to create.
%
%       NESTED          true|[false]
%                       Set this arguement to TRUE if the function is being
%                       used as a nested function.  This will prevent
%                       WAITBARS and FIGURES which makes it easier to
%                       integrate properly.
%
%       SNP_SPOTS       Nx1 array of spots to use as the SNP locations
%
%       FRAGMENTS       Either a 1x1 logical indicating that all sequences
%                       are fragments [TRUE] or all are whole genome
%                       [FALSE].  Or a Nx1 logical array indicating which
%                       sequences are fragments [TRUE] or whole genome
%                       [FALSE].
%
%                       The BASE_CALLS of any spot outside of the
%                       alignment of a fragment is returned as NaN.  They
%                       are also ingnored in determining the SNP frequency
%                       and all other calculations.
%
%       PATDB_ARRAY     Uses the output of LoadPATDB as the source of
%                       patient information.
%
%       PAT_STRUCT_INPUT    Uses a previously created PAT_STRUCT as input.
%
%       USE_DISTRIBUTED     true|[false]
%                               Will use an existing MATLABPOOL and PARFOR
%                               to distribute the work across multiple
%                               cores.
%       CREATE_POOL         true|[false]
%                               Will create its own MATLABPOOL.
%
%
%   See also PatDBbaseCalls, MatchAlignments.
%



%%%%%%%%%%%%%%%%%%%%%%%%%ARGUEMENT CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[INT_NT_MAP ALLOWABLE_CHARS]=AmbigiousNTMap;

ALIGNMENT_PROP={'alphabet', 'nt'};
load ANNOT_DATA HIV_ANNOT_CELL


FIG_TOTAL=5;

DISPLAY_FIG_FLAG=false;
NEED_ALIGNMENT_FLAG=true;
NEED_SNP_SPOTS=true;
NEED_TRANSLATE=false;
NEED_BASE_CALLS=true;
OUTPUT_ALL_ALIGMENT=true;
DISTRIBUTED_FLAG=true;
DISTRIBUTED_POOL_FLAG=false;
NEED_PAT_STRUCT=false;
FIGURES_ONLY=false;
FIGURES_TO_GENERATE=1:FIG_TOTAL;
NESTED_FLAG=false;
FORCE_NO_TRANLATE=false;
FORCE_SNP_SPOT_FLAG=false;


PROVIDED_ALIGNMENTS={};

% WINDOW_CUTOFF=10;
% SNP_CUTOFF_MEAN=.4;

WINDOW_CUTOFF=30;
SNP_CUTOFF_MEAN=.05;

SNP_SPOTS=[];
INPUT_PAT_ARRAY=[];




%%%%%%%CHECK REF_SEQ%%%%%%%%%
if isstruct(REF_SEQ)&&isfield(REF_SEQ,'Sequence')&&numel(REF_SEQ)==1
    REF_SEQ_CELL={SeqPREP(REF_SEQ.Sequence)};
elseif iscell(REF_SEQ)&&numel(REF_SEQ)==1&&ischar(REF_SEQ{1})
    REF_SEQ_CELL={SeqPREP(REF_SEQ{1})};
elseif ischar(REF_SEQ)&&size(REF_SEQ,1)==1&&size(REF_SEQ,2)>1
    REF_SEQ_CELL={SeqPREP(REF_SEQ)};
else
    error('MakeSNPCalls:REF_SEQ','REF_SEQ must be a cell, a char-array or a struct with a Sequence field');
end

REF_SEQ_LENGTH=length(REF_SEQ_CELL{1});


%%%%%%%CHECK TEST_SEQS%%%%%%%%%
if isstruct(TEST_SEQS)&&isfield(TEST_SEQS,'Sequence')
    TEST_SEQS_CELL=arrayfun(@(x)SeqPREP(x.Sequence),TEST_SEQS,'uniformoutput',false);
elseif iscell(TEST_SEQS)
    TEST_SEQS_CELL=cellfun(@SeqPREP,TEST_SEQS,'uniformoutput',false);
else
    error('MakeSNPCalls:TEST_SEQS','TEST_SEQS must be a cell array or a struct-array with a Sequence field');
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
HIDDEN_MASK=false(NUM_TEST_SEQS,1);


%%%%%%%PARSE input arguements
if ~isempty(varargin)
    ParseInputArg
end

NEED_TRANSLATE=NEED_TRANSLATE&~FORCE_NO_TRANLATE;
NEED_SNP_SPOTS=NEED_SNP_SPOTS||FORCE_SNP_SPOT_FLAG;


%%%%%%%DONE ARGUEMENT CHECKING


if ~NESTED_FLAG
    WAITBAR_HANDLE=waitbar(0,'Running');
end

%%%%%%%%%%%DO alignments if needed
if NEED_ALIGNMENT_FLAG
    if ~NESTED_FLAG
        waitbar(0,WAITBAR_HANDLE,'Processing Alignments')
    end
    ALIGNMENTS=cell(NUM_TEST_SEQS,1);

    %%%%%use provided alignments if possible
    if ~isempty(PROVIDED_ALIGNMENTS)
        if ~NESTED_FLAG
            waitbar(0,WAITBAR_HANDLE,'Checking Provided Alignments')
        end
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
                if ~NESTED_FLAG
                    waitbar(0.8*i/NUM_TEST_SEQS,WAITBAR_HANDLE,['Processing Alignments: ' int2str(time/60) ':' int2str(mod(time,60))])
                end
            end
        end

        %%%%add the alignments preformed so they can be returned easily if
        %%%%desired
        PROVIDED_ALIGNMENTS=[PROVIDED_ALIGNMENTS;ALIGNMENTS(still_needed)];


    end
else
    ALIGNMENTS=PROVIDED_ALIGNMENTS;
end


%%%%%%%%%Deal with Fragments%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FRAGMENT_POS=zeros(NUM_TEST_SEQS,2);

if any(FRAGMENT_ARRAY)
    spots=find(FRAGMENT_ARRAY)';
    ind_array=1:REF_SEQ_LENGTH;
    for i=spots
        ref_mask=isletter(ALIGNMENTS{i}(1,:));
        frag_mask=isletter(ALIGNMENTS{i}(3,:));

        %find the last letter in REF seq that occurs before the fragment
        %starts
        FRAGMENT_POS(i,1)=max([0 find(ref_mask(1:find(frag_mask,1,'first')),1,'last')]);


        %find the first letter in REF seq that occurs after the end of the
        %fragment
        FRAGMENT_POS(i,2)=min([find(ref_mask&(ind_array>=find(frag_mask,1,'last')),1,'first') REF_SEQ_LENGTH]);

    end
end


%%%%%%%%%%%%%%Determine SNP_SPOTS if needed
if any(FIGURES_TO_GENERATE==3)||NEED_SNP_SPOTS
    if ~NESTED_FLAG
        waitbar(0.8,WAITBAR_HANDLE,'Determining SNP SPOTS')
    end

    snp_spot=false(NUM_TEST_SEQS,REF_SEQ_LENGTH);
    for i=1:NUM_TEST_SEQS
        mask=isletter(ALIGNMENTS{i}(1,:));
        snp_spot(i,:)=ALIGNMENTS{i}(1,mask)~=ALIGNMENTS{i}(3,mask);
    end
    [forward_looking reverse_looking]=VFindRuns(~snp_spot);
    snp_spot=double(snp_spot);

    if any(FRAGMENT_ARRAY)
        spots=find(FRAGMENT_ARRAY)';
        ind_array=1:REF_SEQ_LENGTH;
        for i=spots
            outside_spots=ind_array<=FRAGMENT_POS(i,1)|ind_array>=FRAGMENT_POS(i,2);
            snp_spot(i,outside_spots)=NaN(1,nnz(outside_spots));
            forward_looking(i,outside_spots)=NaN(1,nnz(outside_spots));
            reverse_looking(i,outside_spots)=NaN(1,nnz(outside_spots));
        end
    end

    if ~NESTED_FLAG
        waitbar(0.9,WAITBAR_HANDLE,'Determining SNP SPOTS')
    end

    values=nanmean(forward_looking(~HIDDEN_MASK,:)+reverse_looking(~HIDDEN_MASK,:)>WINDOW_CUTOFF)>SNP_CUTOFF_MEAN;
    snps=nanmean(snp_spot(~HIDDEN_MASK,:))>SNP_CUTOFF_MEAN;

    if NEED_SNP_SPOTS
        SNP_SPOTS=find(snps&values);
    end
    IS_SNPs=snp_spot(:,SNP_SPOTS);

end

%%%%%%%%%%%%%%%%%%NO SNPs found so exit and display message
if isempty(SNP_SPOTS)

    varargout=SetupOutput(nargout);
    warning('MakeSNPCalls:NO_SNP_SPOTS','The Inputs did not have any recognizable SNPs.  Try lowering the SNP_CUTOFF_MEAN.');
    return

end

BASE_CALLS=cell2mat(arrayfun(@GetLetters,ALIGNMENTS(:),FRAGMENT_ARRAY(:),FRAGMENT_POS(:,1),FRAGMENT_POS(:,2),'uniformoutput',false));


if NEED_TRANSLATE||(any(FIGURES_TO_GENERATE==4)&&~NESTED_FLAG)
    if ~NESTED_FLAG
        waitbar(0.9,WAITBAR_HANDLE,'Determining MUTATIONS')
    end

    map=geneticcode;
    trans_table=[fieldnames(map) struct2cell(map)];
    stop_indexes=find(strcmp('*',trans_table(:,2)));

    MUT_MAT=zeros(size(BASE_CALLS));

    start_locs=[HIV_ANNOT_CELL{2:end-1,1}];
    start_locs=start_locs(1:2:end);

    ref_seq=REF_SEQ_CELL{1};

    [junk correct_locs]=arrayfun(@(x)(Translate(ref_seq(x:end))),start_locs,'uniformoutput',false);

    tic
    for i=1:NUM_TEST_SEQS

        test_seq=ref_seq;
        test_seq(SNP_SPOTS(~isnan(BASE_CALLS(i,:))))=int2nt(BASE_CALLS(i,~isnan(BASE_CALLS(i,:))));

        for k=1:length(start_locs)
            [test_spots test_locs]=Translate(test_seq(start_locs(k):end));

            test_spots=test_spots+start_locs(k);

            %Check to make sure lengths are equal before checking for spot
            %equality
            if ~(length(test_locs)==length(correct_locs{k})&&all(test_locs==correct_locs{k}))
                if length(test_locs)==length(correct_locs{k})
                    miss_sense_rows=test_spots(test_locs~=correct_locs{k},:);
                    [junk locs]=ismember(miss_sense_rows(:),SNP_SPOTS);
                    MUT_MAT(i,nonzeros(locs))=ones(size(nonzeros(locs)));

                else
                    min_val=min([length(test_locs) length(correct_locs{k})]);
                    miss_sense_rows=test_spots(test_locs(1:min_val)~=correct_locs{k}(1:min_val),:);
                    [junk locs]=ismember(miss_sense_rows(1:end-3),SNP_SPOTS);
                    MUT_MAT(i,nonzeros(locs))=ones(size(nonzeros(locs)));

                    [junk locs]=ismember(miss_sense_rows(end-2:end),SNP_SPOTS);
                    MUT_MAT(i,nonzeros(locs))=2*ones(size(nonzeros(locs)));
                end
            end
        end
        time=(toc/i)*(NUM_TEST_SEQS-i);
        if ~NESTED_FLAG
            waitbar(0.9+0.1*i/NUM_TEST_SEQS,WAITBAR_HANDLE,['Processing Mutations: ' int2str(time/60) ':' int2str(mod(time,60))])
        end
    end
end


if any(FIGURES_TO_GENERATE)&~NESTED_FLAG
    GenerateFigures
end


if ~NESTED_FLAG
    waitbar(1,WAITBAR_HANDLE,'Finishing')
end
%%%%%%%%%%%%%%%%Setup outputs
varargout=SetupOutput(nargout);
if ~NESTED_FLAG
    close(WAITBAR_HANDLE)
end

%%%%%%%%%%%%%%%%%%SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function calls=GetLetters(input,fragment_flag,start,stop)
        spots=find(isletter(input{1}(1,:)));
        if ~fragment_flag
            calls=double(nt2int(input{1}(3,spots(SNP_SPOTS))));
        else
            valid_spots=(SNP_SPOTS>start&SNP_SPOTS<stop);
            calls=NaN(size(SNP_SPOTS));
            calls(valid_spots)=double(nt2int(input{1}(3,spots(SNP_SPOTS(valid_spots)))));
        end
    end




    function [SPOTS LOCS]=Translate(this_seq)
        AA_seq=[];
        if isempty(this_seq)
            warning('GetAlignedSequences:NoSeq',['No sequence '...
                'information found in entry %d'])
        else
            stripped_seq=this_seq(isletter(this_seq));
            start_spot=regexpi(stripped_seq,map.Starts{1},'once');
            if isempty(start_spot)
                warning('GetAlignedSequences:NoStart',['No start codon '...
                    'found in Seq %d'])
            else

                num_rows=floor((length(stripped_seq)-start_spot)/3);

                codon_array=num2cell(reshape(stripped_seq(start_spot:start_spot+num_rows*3-1),3,num_rows)',2);

                SPOTS=reshape(start_spot:start_spot+num_rows*3-1,3,num_rows)';

                [junk LOCS]=ismember(codon_array,trans_table(:,1));
                stop_spot=find(ismember(LOCS,stop_indexes),1);
                if ~isempty(stop_spot)
                    LOCS=LOCS(1:stop_spot);
                    SPOTS=SPOTS(1:stop_spot,:);
                end

            end
        end
    end


    function GenerateFigures()


        %%%%%%Figure displaying the alignment of fragments if needed

        if any(FRAGMENT_ARRAY)&&any(FIGURES_TO_GENERATE==1)
            figure
            coverage=zeros(1,REF_SEQ_LENGTH);
            spots=find(FRAGMENT_ARRAY)';
            for k=spots
                coverage(FRAGMENT_POS(k,1):FRAGMENT_POS(k,2))=coverage(FRAGMENT_POS(k,1):FRAGMENT_POS(k,2))+1;
            end

            bar(1:REF_SEQ_LENGTH,coverage)
            axis([0 REF_SEQ_LENGTH 0 1.1*nnz(FRAGMENT_ARRAY)])
            title('Fragment Coverage Map')
            ylabel('#Aligned Fragments')
            xlabel('HIV Sequence')
            AnnotFig(1:REF_SEQ_LENGTH,HIV_ANNOT_CELL)
        end


        %%%%%%%%%%%%%%%%%%Figure displaying the BASES at the SNP LOCATIONS
        if any(FIGURES_TO_GENERATE==2)
            replace_mat=zeros(5,size(BASE_CALLS,2));
            for k=0:4
                replace_mat(k+1,:)=nansum(BASE_CALLS==k);
            end

            figure
            title_color={'Gap','r';'A','g';'C','b';'G','k';'T','c'};
            for fig=1:5
                subplot(5,1,fig), bar(SNP_SPOTS,replace_mat(fig,:)',title_color{fig,2});
                title(title_color{fig,1})
                axis([0 length(ref_seq) 0 1.1*size(BASE_CALLS,1)]);
                xlabel('Location')
                ylabel('Number')
                AnnotFig(1:REF_SEQ_LENGTH,HIV_ANNOT_CELL)
            end
        end

        %%%%%%%%%%%%%%%Figure displaying the frequency of SNPs


        run_vals=log2(nanmean(double(forward_looking+reverse_looking),2));
        snps_per_seq=nansum(snp_spot,2);
        seqs_with_snps=nansum(snp_spot,1);

        if any(FIGURES_TO_GENERATE==3)
            figure
            subplot(3,1,1), hist(snps_per_seq,50)
            title('SNPs per Sequence')

            subplot(3,1,2), hist(seqs_with_snps,50)
            title('Seqs with SNP')

            subplot(3,1,3), hist(run_vals,100)
            title('Average Island Length')
            xlabel('Log_2(Island Length)')
            ylabel('Frequency')
        end
        %%%%%%%%%%%%%%%Figure displaying the Frequeny of MUTATIONS at a SNP LOCATION

        if any(FIGURES_TO_GENERATE==4)
            figure
            title_color={'Silent','g';'Mis-sense','k';'Non-sense','r'};
            for fig=1:3
                subplot(3,1,fig), bar(SNP_SPOTS,mean(MUT_MAT==(fig-1)),title_color{fig,2})
                title(title_color{fig,1})
                axis([0 REF_SEQ_LENGTH 0 1.1]);
                xlabel('Location')
                ylabel('Frequency')
                AnnotFig(1:REF_SEQ_LENGTH,HIV_ANNOT_CELL)
            end
        end

        if any(FIGURES_TO_GENERATE==5)
            [ELM_spots ELM_name ELM_match]=ParseDatabaseResults('ELM_RESULTS.elm',HIV_ANNOT_CELL,'ELM');
            ELM_annot=[num2cell(ELM_spots,2) ELM_name];

            good_elms=cellfun(@(x)(any(ismember(x(1):x(2),SNP_SPOTS))),ELM_annot(:,1));

            %             figure
            %             subplot(3,1,1), bar(SNP_SPOTS,mean(MUT_MAT==2),title_color{fig,2})
            %             axis([0 REF_SEQ_LENGTH 0 1.1]);
            %             AnnotFig(1:REF_SEQ_LENGTH,HIV_ANNOT_CELL)
            %             AnnotFig(1:REF_SEQ_LENGTH,ELM_annot)
            %             xlabel('Location')
            %             ylabel('Frequency')
            %             title('ALL ELM locations')
            %
            %             subplot(3,1,2), bar(SNP_SPOTS,mean(MUT_MAT==2),title_color{fig,2})
            %             axis([0 REF_SEQ_LENGTH 0 1.1]);
            %             AnnotFig(1:REF_SEQ_LENGTH,HIV_ANNOT_CELL)
            %             AnnotFig(1:REF_SEQ_LENGTH,ELM_annot(good_elms,:))
            %             xlabel('Location')
            %             ylabel('Frequency')
            %             title('SNP-ed ELM locations')
            %
            %
            %             [ELM_snps ELM_counts]=cellfun(@CheckSingleELM,ELM_annot(:,1));
            %
            %
            %             if any(ELM_counts)
            %
            %                 ELM_snp_sums=sum(MUT_MAT(:,unique(nonzeros(ELM_snps)))==2);
            %
            %                 subplot(3,1,3), bar(SNP_SPOTS(unique(nonzeros(ELM_snps))),ELM_snp_sums,title_color{fig,2})
            %                 axis([0 REF_SEQ_LENGTH 0 1.1*max(ELM_snp_sums)]);
            %                 AnnotFig(1:REF_SEQ_LENGTH,HIV_ANNOT_CELL)
            %                 AnnotFig(1:REF_SEQ_LENGTH,ELM_annot(good_elms,:))
            %                 xlabel('Location')
            %                 ylabel('Count')
            %                 title('Mutated ELM locations')
            %
            %
            %                 figure
            %                 barh(ELM_counts(ELM_counts>0));
            %                 set(gca,'ytick',1:nnz(ELM_counts),'yticklabel',ELM_annot(ELM_counts>0,2));
            %             end

            WIN_SIZE=100;
            win_sim=zeros(1,REF_SEQ_LENGTH);
            for p=1:REF_SEQ_LENGTH-WIN_SIZE
                temp=snp_spot(:,p:p+25);
                win_sim(p)=1-mean(nanmean(temp(:)));
            end

            figure
            plot(win_sim)
            AnnotFig(1:REF_SEQ_LENGTH,HIV_ANNOT_CELL)
            AnnotFig(1:REF_SEQ_LENGTH,ELM_annot)
            xlabel('Location')
            ylabel('% Similarity with REF')

            ELM_sim=cellfun(@(x)(mean(win_sim(x(1):x(2)))),ELM_annot(:,1));

            figure
            barh(ELM_sim)
            set(gca,'ytick',1:length(ELM_sim),'yticklabel',ELM_annot(:,2));
            xlabel('% Similarity with REF')


        end

        %%%%%%%%%%%%%%%%%Display legend for the HIV proteins

        AnnotFig([],HIV_ANNOT_CELL,1)

    end

    function [SNPs NUM_OCCURS]=CheckSingleELM(start_stop)
        [TF LOC]=ismember(start_stop(1):start_stop(2),SNP_SPOTS);

        if any(TF)
            NUM_OCCURS=sum(sum(MUT_MAT(:,LOC(TF))==2));
            SNPs=max(LOC);
        else
            NUM_OCCURS=0;
            SNPs=0;
        end
    end




%%%%%%%%%%%%%%ARGUEMENT CHECKING FUNCTIONS
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

    function ValidateAlign(align_input,test_seq)
        ref_align=SeqPREP(align_input(1,:));
        test_align=SeqPREP(align_input(3,:));

        if ~strcmp(ref_align,REF_SEQ_CELL{1})||~strcmp(test_align,test_seq)
            error('MakeSNPCalls:BAD_ALIGN','ALIGNMENT must a cell array of the outputs of NWALIGN(REF_SEQ,TEST_SEQ)');
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

    function output_cell=SetupOutput(input)
        output_cell=cell(input,1);
        start=1;
        temp=[];
        if NEED_PAT_STRUCT
            %%%%%ADD PAT_STRUCT IF DESIRED
            if NEED_TRANSLATE
                output_cell{start}=PatientStructHelper(INPUT_PAT_ARRAY,...
                    {'Alignment_CELL',PAT_IND_ARRAY,ALIGNMENTS},...
                    {'BASE_CALLS',PAT_IND_ARRAY,BASE_CALLS},...
                    {'SNP_SPOTS',PAT_IND_ARRAY,ones(length(PAT_IND_ARRAY),1)*SNP_SPOTS},...     %SNP_SPOTS should be copied identically to each element
                    {'MUT_CALLS',PAT_IND_ARRAY,MUT_MAT},...
                    {'IS_FRAGMENT',PAT_IND_ARRAY,FRAGMENT_ARRAY},...
                    {'Sequence_CELL',PAT_IND_ARRAY,TEST_SEQS_CELL},...
                    {'IS_SNPs',PAT_IND_ARRAY,IS_SNPs});
            else
                output_cell{start}=PatientStructHelper(INPUT_PAT_ARRAY,...
                    {'Alignment_CELL',PAT_IND_ARRAY,ALIGNMENTS},...
                    {'BASE_CALLS',PAT_IND_ARRAY,BASE_CALLS},...
                    {'SNP_SPOTS',PAT_IND_ARRAY,ones(length(PAT_IND_ARRAY),1)*SNP_SPOTS},...     %SNP_SPOTS should be copied identically to each element
                    {'IS_FRAGMENT',PAT_IND_ARRAY,FRAGMENT_ARRAY},...
                    {'Sequence_CELL',PAT_IND_ARRAY,TEST_SEQS_CELL},...
                    {'IS_SNPs',PAT_IND_ARRAY,IS_SNPs});
            end
            start=2;
        end

        switch input-(start-1)
            case 1
                output_cell{start}=SNP_SPOTS;
            case 2

                output_cell{start}=SNP_SPOTS;
                output_cell{start+1}=IS_SNPs;
            case 3

                output_cell{start}=SNP_SPOTS;
                output_cell{start+1}=BASE_CALLS;
                if OUTPUT_ALL_ALIGMENT
                    output_cell{start+2}=PROVIDED_ALIGNMENTS;
                else
                    output_cell{start+2}=ALIGNMENTS;
                end
            case 4
                output_cell{start}=SNP_SPOTS;
                output_cell{start+1}=BASE_CALLS;
                if OUTPUT_ALL_ALIGMENT
                    output_cell{start+2}=PROVIDED_ALIGNMENTS;
                else
                    output_cell{start+2}=ALIGNMENTS;
                end
                output_cell{start+2}=MUT_CALLS;

        end
    end

    function ParseInputArg
        for i=1:2:length(varargin)
            switch lower(varargin{i})
                case {'force snp_spots'}
                    FORCE_SNP_SPOT_FLAG=varargin{i+1};
                case {'alignments'}
                    if iscell(varargin{i+1})
                        PROVIDED_ALIGNMENTS=[PROVIDED_ALIGNMENTS;varargin{i+1}];

                    else
                        error('MakeSNPCalls:BAD_ALIGNMENTS','Alignments must be ');
                    end
                case {'alignment_prop', 'alignment prop'}
                    if iscell(varargin{i+1})
                        ALIGNMENT_PROP=varargin{i+1};
                        try     %test with trivial example to make sure the input is valid
                            nwalign_mod('ACGT','AAGTC',ALIGNMENT_PROP{:});
                        catch
                            warning('MakeSNPCalls:BAD_ALIGNMENT_PROP','Invalid nwalign arguements provided.');
                            rethrow(lasterror)
                        end
                    else
                        error('MakeSNPCalls:BAD_ALIGNMENT_PROP','Invalid nwalign arguements provided.');
                    end
                case {'snp_spots', 'snp spots'}
                    if isnumeric(varargin{i+1})
                        SNP_SPOTS=sort(varargin{i+1}(:))';
                        if SNP_SPOTS(end)>REF_SEQ_LENGTH
                            error('MakeSNPCalls:BAD_SNP_SPOT','A SNP_SPOT was provided which is outside of the range of REF_SEQ [0 %d] \n found: %d',REF_SEQ_LENGTH,SNP_SPOTS(end));
                        end
                        NEED_SNP_SPOTS=false;
                    else
                        error('MakeSNPCalls:BAD_SNP_SPOT_ARG','Arguement to SNP_SPOT must be a numeric array');
                    end
                case {'fragments'}
                    if islogical(varargin{i+1})
                        if numel(varargin{i+1})==1
                            FRAGMENT_ARRAY=true(NUM_TEST_SEQS,1)&varargin{i+1};
                        elseif numel(varargin{i+1})==NUM_TEST_SEQS
                            FRAGMENT_ARRAY=varargin{i+1}(:);
                        else
                            error('MakeSNPCalls:BAD_FRAGMENTS_ARG','Arguement to FRAGMENTS must be a logical array of either 1x1 or the same length as TEST_SEQS');
                        end
                    else
                        error('MakeSNPCalls:BAD_FRAGMENTS_ARG','Arguement to FRAGMENTS must be a logical array');
                    end
                case {'patdb_array','patbd array'}
                    INPUT_PAT_ARRAY=varargin{i+1};
                    [temp_seqs temp_fragment temp_pat_ind_array]=PatDBbaseCalls(REF_SEQ,INPUT_PAT_ARRAY,1,'seq_only');

                    PAT_IND_ARRAY=[PAT_IND_ARRAY;temp_pat_ind_array];

                    TEST_SEQS_CELL=[TEST_SEQS_CELL; temp_seqs];
                    FRAGMENT_ARRAY=[FRAGMENT_ARRAY;temp_fragment];
                    NUM_TEST_SEQS=length(TEST_SEQS_CELL);

                case {'trust order','trust_order'}
                    if islogical(varargin{i+1})&&varargin{i+1}
                        NEED_ALIGNMENT_FLAG=false;
                    else
                        error('MakeSNPCalls:TRUST_ORDER_ARG','Arguement to TRUST_ORDER must be logical');
                    end

                case {'ordered alignment output','ordered_alignment_output'}
                    if islogical(varargin{i+1})
                        OUTPUT_ALL_ALIGMENT=~varargin{i+1};
                    else
                        error('MakeSNPCalls:ORDERED_ALIGN_ARG','Arguement to ORDERED_ALIGN_OUTPUT must be logical');
                    end
                case {'use_distributed'}
                    if islogical(varargin{i+1})
                        DISTRIBUTED_FLAG=varargin{i+1};
                    else
                        error('MakeSNPCalls:USE_DISTRIBUTED_ARG','Arguement to USE_DISTRIBUTED must be logical');
                    end
                case {'create_pool'}
                    if islogical(varargin{i+1})
                        DISTRIBUTED_POOL_FLAG=varargin{i+1};
                    else
                        error('MakeSNPCalls:USE_DISTRIBUTED_ARG','Arguement to DISTRIBUTED_POOL must be logical');
                    end
                case {'pat_struct_out'}
                    if islogical(varargin{i+1})
                        NEED_PAT_STRUCT=varargin{i+1};
                    else
                        error('MakeSNPCalls:NEED_PAT_STRUCT_ARG','Arguement to NEED_PAT_STRUCT must be logical');
                    end
                case {'pat_struct_input'}
                    %%%%%Retrieve the needed information from the Patient_Structure
                    INPUT_PAT_ARRAY=varargin{i+1};

                    [temp_SNP_SPOTS SNP_SPOTS_IND ...
                        temp_ALIGNMENT ALIGNMENTS_IND ...
                        temp_seqs temp_seqs_ind ...
                        temp_fragment temp_fragment_ind] =...
                        PatientStructHelper(INPUT_PAT_ARRAY,...
                        {'SNP_SPOTS','explodeNumeric'},...
                        {'Alignment_CELL','explodeCell'},...
                        {'Sequence_CELL','explodeCell'},...
                        {'IS_FRAGMENT','explodeNumeric'});
                    %%%%Check to make sure everything came back properly

                    PAT_IND_ARRAY=[PAT_IND_ARRAY ; temp_seqs_ind];
                    TEST_SEQS_CELL=[TEST_SEQS_CELL; temp_seqs];
                    NUM_TEST_SEQS=length(TEST_SEQS_CELL);
                    HIDDEN_MASK=false(NUM_TEST_SEQS,1);

                    %%SNP_SPOTS should be the same on every row
                    if isempty(temp_SNP_SPOTS)||any(sum(temp_SNP_SPOTS)/size(temp_SNP_SPOTS,1)~=temp_SNP_SPOTS(1,:))
                        warning('MakeSNPCalls:BAD_PAT_STRUCT_SNP_SPOTS','The information from the PAT_STRUCT on SNP_SPOTS was missing or inconsitant.  Will calculate again.');
                    else
                        NEED_SNP_SPOTS=false;
                        SNP_SPOTS=temp_SNP_SPOTS(1,:);
                    end

                    if isempty(temp_ALIGNMENT)
                        warning('MakeSNPCalls:BAD_PAT_STRUCT_ALIGNEMNTS','The information from the PAT_STRUCT on ALIGNEMNTS was missing or inconsitant.  Will calculate again.');
                        [temp_seqs temp_fragment PAT_IND_ARRAY]=PatDBbaseCalls(REF_SEQ,INPUT_PAT_ARRAY,1,'seq_only');

                        PAT_IND_ARRAY=[zeros(NUM_TEST_SEQS,1) ; PAT_IND_ARRAY];

                        TEST_SEQS_CELL=[TEST_SEQS_CELL; temp_seqs];
                        FRAGMENT_ARRAY=[FRAGMENT_ARRAY;temp_fragment];
                        NUM_TEST_SEQS=length(TEST_SEQS_CELL);
                        HIDDEN_MASK=false(NUM_TEST_SEQS,1);

                    else
                        PROVIDED_ALIGNMENTS=[PROVIDED_ALIGNMENTS;temp_ALIGNMENT];
                    end

                    if isempty(temp_fragment)
                        warning('MakeSNPCalls:BAD_PAT_STRUCT_FRAGMENT','The information from the PAT_STRUCT on FRAGMENT was missing or inconsitant.  Assuming all patient data is a fragment.');
                        FRAGMENT_ARRAY=[FRAGMENT_ARRAY;ones(length(temp_seqs),1)];
                    else
                        FRAGMENT_ARRAY=[FRAGMENT_ARRAY;temp_fragment];
                    end
                case {'fig_only','fig only','figure_only','figure only'}
                    if islogical(varargin{i+1})
                        FIGURES_TO_GENERATE=1:FIG_TOTAL*varargin{i+1};
                    elseif isvector(varargin{i+1})&&all(varargin{i+1}>0&&varargin{i+1}<=FIG_TOTAL)
                        FIGURES_TO_GENERATE=varargin{i+1};
                    else
                        error('MakeSNPCalls:BAD_FIG_ONLY_ARG','Arguement to FIG_ONLY must either be a logical scalar or a numerical vector [1 %d].',FIG_TOTAL);
                    end
                case {'nested'}
                    if islogical(varargin{i+1})
                        NESTED_FLAG=varargin{i+1};
                    else
                        error('MakeSNPCalls:BAD_NESTED_ARG','Arguement to NESTED must be a logical scalar.',FIG_TOTAL);
                    end
                case {'hidden_patients','hidden_pats','hidden'}
                    if islogical(varargin{i+1})&&length(varargin{i+1})==length(INPUT_PAT_ARRAY)
                        temp_hidden=varargin{i+1}(:);
                        hidden_patients=find(temp_hidden);

                        TF=ismember(PAT_IND_ARRAY,hidden_patients);

                        HIDDEN_MASK=false(NUM_TEST_SEQS,1);
                        HIDDEN_MASK(TF)=true(nnz(TF),1);
                        %%%%%%%%%%%%%%%%%%STILL NEED TO CORRECT HERE


                    else
                        error('MakeSNPCalls:BAD_HIDDEN_ARG','Arguement to HIDDEN must be a logical vector the same size as INPUT_PAT_ARRAY.',FIG_TOTAL);

                    end
                case {'no_translate'}
                    FORCE_NO_TRANLATE=varargin{i+1};



                otherwise
                    error('MakeSNPCalls:BAD_ARG','An unknown input arguements was provided: %s',varargin{i});
            end
        end

    end


end