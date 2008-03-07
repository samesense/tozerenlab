function GenerateSimFigures(REF_SEQ,TEST_SEQS,varargin)

load TESTING_DATA HIV_ANNOT_CELL
[INT_NT_MAP ALLOWABLE_CHARS]=AmbigiousNTMap;

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

PROVIDED_ALIGNMENTS={};

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

%%%%%%%PARSE input arguements
if ~isempty(varargin)
    ParseInputArg
end


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


%%%%%%%%%%%%%%Determine the SNP at each position

snp_spot=zeros(NUM_TEST_SEQS,REF_SEQ_LENGTH);

for i=1:NUM_TEST_SEQS
    spots=find(isletter(ALIGNMENTS{i}(1,:)));
    snp_spot(i,:)=ALIGNMENTS{i}(1,spots)~=ALIGNMENTS{i}(3,spots);
end

if any(FRAGMENT_ARRAY)
    spots=find(FRAGMENT_ARRAY)';
    ind_array=1:REF_SEQ_LENGTH;
    for i=spots
        outside_spots=ind_array<=FRAGMENT_POS(i,1)|ind_array>=FRAGMENT_POS(i,2);
        snp_spot(i,outside_spots)=NaN(1,nnz(outside_spots));
    end
end


[ELM_spots ELM_name ELM_match]=ParseDatabaseResults('ELM_RESULTS.elm',HIV_ANNOT_CELL,'ELM');
ELM_annot=[num2cell(ELM_spots,2) ELM_name];

good_elms=cellfun(@(x)(any(ismember(x(1):x(2),SNP_SPOTS))),ELM_annot(:,1));

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
title('All ELM')


[MATCH_spots MATCH_name MATCH_match]=ParseDatabaseResults('MATCH_results.txt',HIV_ANNOT_CELL,'MATCH');
MATCH_annot=[num2cell(MATCH_spots,2) MATCH_name];

figure
plot(win_sim)
AnnotFig(1:REF_SEQ_LENGTH,HIV_ANNOT_CELL)
AnnotFig(1:REF_SEQ_LENGTH,MATCH_annot)
xlabel('Location')
ylabel('% Similarity with REF')
title('All Transcription Factor BS')

NEED_EXCEL_OUT_FLAG



%%%%%%%%%%%%%%%%%%%%%%%HELPER FUNCTIONS
    function GenerateExcelReport
        xlswrite(EXCEL_OUT_FILENAME,[])
        
        
        
        
    end





%%%%%%%%%%%%%%%ARGUEMENT CHEKCING FUNCTIONS
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

    function ParseInputArg
        for i=1:2:length(varargin)
            switch lower(varargin{i})
                case {'alignments'}
                    if iscell(varargin{i+1})
                        PROVIDED_ALIGNMENTS=[PROVIDED_ALIGNMENTS;varargin{i+1}];

                    else
                        error('GenerateSimFigures:BAD_ALIGNMENTS','Alignments must be ');
                    end
                case {'alignment_prop', 'alignment prop'}
                    if iscell(varargin{i+1})
                        ALIGNMENT_PROP=varargin{i+1};
                        try     %test with trivial example to make sure the input is valid
                            nwalign_mod('ACGT','AAGTC',ALIGNMENT_PROP{:});
                        catch
                            warning('GenerateSimFigures:BAD_ALIGNMENT_PROP','Invalid nwalign arguements provided.');
                            rethrow(lasterror)
                        end
                    else
                        error('GenerateSimFigures:BAD_ALIGNMENT_PROP','Invalid nwalign arguements provided.');
                    end
                case {'fragments'}
                    if islogical(varargin{i+1})
                        if numel(varargin{i+1})==1
                            FRAGMENT_ARRAY=true(NUM_TEST_SEQS,1)&varargin{i+1};
                        elseif numel(varargin{i+1})==NUM_TEST_SEQS
                            FRAGMENT_ARRAY=varargin{i+1}(:);
                        else
                            error('GenerateSimFigures:BAD_FRAGMENTS_ARG','Arguement to FRAGMENTS must be a logical array of either 1x1 or the same length as TEST_SEQS');
                        end
                    else
                        error('GenerateSimFigures:BAD_FRAGMENTS_ARG','Arguement to FRAGMENTS must be a logical array');
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
                        error('GenerateSimFigures:TRUST_ORDER_ARG','Arguement to TRUST_ORDER must be logical');
                    end

                case {'ordered alignment output','ordered_alignment_output'}
                    if islogical(varargin{i+1})
                        OUTPUT_ALL_ALIGMENT=~varargin{i+1};
                    else
                        error('GenerateSimFigures:ORDERED_ALIGN_ARG','Arguement to ORDERED_ALIGN_OUTPUT must be logical');
                    end
                case {'use_distributed'}
                    if islogical(varargin{i+1})
                        DISTRIBUTED_FLAG=varargin{i+1};
                    else
                        error('GenerateSimFigures:USE_DISTRIBUTED_ARG','Arguement to USE_DISTRIBUTED must be logical');
                    end
                case {'create_pool'}
                    if islogical(varargin{i+1})
                        DISTRIBUTED_POOL_FLAG=varargin{i+1};
                    else
                        error('GenerateSimFigures:USE_DISTRIBUTED_ARG','Arguement to DISTRIBUTED_POOL must be logical');
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

                    %%SNP_SPOTS should be the same on every row
                    if isempty(temp_SNP_SPOTS)||any(sum(temp_SNP_SPOTS)/size(temp_SNP_SPOTS,1)~=temp_SNP_SPOTS(1,:))
                        warning('GenerateSimFigures:BAD_PAT_STRUCT_SNP_SPOTS','The information from the PAT_STRUCT on SNP_SPOTS was missing or inconsitant.  Will calculate again.');
                    else
                        NEED_SNP_SPOTS=false;
                        SNP_SPOTS=temp_SNP_SPOTS(1,:);
                    end

                    if isempty(temp_ALIGNMENT)
                        warning('GenerateSimFigures:BAD_PAT_STRUCT_ALIGNEMNTS','The information from the PAT_STRUCT on ALIGNEMNTS was missing or inconsitant.  Will calculate again.');
                        [temp_seqs temp_fragment PAT_IND_ARRAY]=PatDBbaseCalls(REF_SEQ,INPUT_PAT_ARRAY,1,'seq_only');

                        PAT_IND_ARRAY=[zeros(NUM_TEST_SEQS,1) ; PAT_IND_ARRAY];

                        TEST_SEQS_CELL=[TEST_SEQS_CELL; temp_seqs];
                        FRAGMENT_ARRAY=[FRAGMENT_ARRAY;temp_fragment];
                        NUM_TEST_SEQS=length(TEST_SEQS_CELL);
                    else
                        PROVIDED_ALIGNMENTS=[PROVIDED_ALIGNMENTS;temp_ALIGNMENT];
                    end

                    if isempty(temp_fragment)
                        warning('GenerateSimFigures:BAD_PAT_STRUCT_FRAGMENT','The information from the PAT_STRUCT on FRAGMENT was missing or inconsitant.  Assuming all patient data is a fragment.');
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
                        error('GenerateSimFigures:BAD_FIG_ONLY_ARG','Arguement to FIG_ONLY must either be a logical scalar or a numerical vector [1 %d].',FIG_TOTAL);
                    end
                    
                case {'excel output','excel_output','excel'}
                    if ischar(varargin{i+1})
                        NEED_EXCEL_OUT_FLAG=true;
                        EXCEL_OUT_FILENAME=varargin{i+1};
                        if strcmp(EXCEL_OUT_FILENAME(end-4:end),'.xls')
                            EXCEL_OUT_FILENAME=[EXCEL_OUT_FILENAME '.xls'];
                        end
                    end
                   
                otherwise
                    error('GenerateSimFigures:BAD_ARG','An unknown input arguements was provided: %s',varargin{i});
            end
        end

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