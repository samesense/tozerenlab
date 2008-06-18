function [outputData, varargout]=FindmiRNAislands(REF_SEQ,TEST_SEQS,MIN_LENGTH,CONS_CUTOFF,varargin)
%   FindmiRNAislands
%       Finds islands of homology within the highly varying sequence using
%       an alignment to a reference sample.
%
%       OUTPUT_DATA = FindmiRNAislands(REF_SEQ,TEST_SEQS,MIN_LENGTH,CONS_CUTOFF)
%
%       REF_SEQ         The reference sequence for aligning.
%
%       TEST_SEQS       The whole genome sequences to test.
%
%       MIN_LENGTH      The minumum island length required.
%
%       CONS_CUTOFF     A vector of conservation cutoffs to check the
%                       islands.
%
%
%   Optional Parameters:
%       ALIGNMENTS      ALIGNMENT_CELL
%                           A Nx1 cell of alignments to the reference
%                           sample.  They do not need to be in any order,
%                           they are re-arranged as needed.
%
%
%       ALIGNMENT_PROP  A cell array of properties to pass to NWALIGN
%
%
%       EXCEL_OUTPUT    A filename to write the sequence islands to.
%
%       AXES_HANDLE     A handle to an AXES for plotting into.
%
%
%
%






GLOBAL_CUTOFF=0;
PARFOR_BLOCK_SIZE=20;
ALIGN_ONLY_FLAG=false;
[INT_NT_MAP ALLOWABLE_CHARS]=AmbigiousNTMap;
FIG_TOTAL=1;
AX_HANDLE=[];
FIGURES_TO_GENERATE=1:FIG_TOTAL;
EXCEL_FILENAME=[];
ALIGNMENT_PROP={'alphabet', 'nt'};

%%%%%%%CHECK INPUT NUMS%%%%%%
if ~isnumeric(MIN_LENGTH)||~isscalar(MIN_LENGTH)
    error('FindmiRNAislands:MIN_LENGTH','MIN_LENGTH must be a scalar double!');
end

if ~isnumeric(CONS_CUTOFF)
    error('FindmiRNAislands:CONS_CUTOFF','CONS_CUTOFF must be a double!');
end


%%%%%%%%%%%%%%DEAL WITH MULTI REF_SEQS

if isstruct(REF_SEQ)&&length(REF_SEQ)>1

    outputData = cell(size(REF_SEQ));

    for ind = 1:length(REF_SEQ)
        outputData{ind} = FindmiRNAislands(REF_SEQ(ind),TEST_SEQS,MIN_LENGTH,CONS_CUTOFF,varargin{:},'fig_only',false);
    end

    
    return

end






%%%%%%%%%%%%%%%DEAL WITH SINGLE REF_SEQS

%%%%%%%CHECK REF_SEQ%%%%%%%%%
if isstruct(REF_SEQ)&&isfield(REF_SEQ,'Sequence')&&numel(REF_SEQ)==1
    REF_SEQ_CELL={SeqPREP(REF_SEQ.Sequence)};
elseif iscell(REF_SEQ)&&numel(REF_SEQ)==1&&ischar(REF_SEQ{1})
    REF_SEQ_CELL={SeqPREP(REF_SEQ{1})};
elseif ischar(REF_SEQ)&&size(REF_SEQ,1)==1&&size(REF_SEQ,2)>1
    REF_SEQ_CELL={SeqPREP(REF_SEQ)};
else
    error('FindmiRNAislands:REF_SEQ','REF_SEQ must be a cell, a char-array or a struct with a Sequence field');
end

REF_SEQ_LENGTH=length(REF_SEQ_CELL{1});


%%%%%%%CHECK TEST_SEQS%%%%%%%%%
if isstruct(TEST_SEQS)&&isfield(TEST_SEQS,'Sequence')
    TEST_SEQS_CELL=arrayfun(@(x)SeqPREP(x.Sequence),TEST_SEQS,'uniformoutput',false);
elseif iscell(TEST_SEQS)
    TEST_SEQS_CELL=cellfun(@SeqPREP,TEST_SEQS,'uniformoutput',false);
else
    error('FindmiRNAislands:TEST_SEQS','TEST_SEQS must be a cell array or a struct-array with a Sequence field');
end
TEST_SEQS_CELL=TEST_SEQS_CELL(:);


bad_seqs=cellfun('isempty',TEST_SEQS_CELL);

if any(bad_seqs)
    display([int2str(nnz(bad_seqs)) ' Sequences skipped']);
    TEST_SEQS_CELL=TEST_SEQS_CELL(~bad_seqs);
end

NUM_TEST_SEQS=length(TEST_SEQS_CELL);

NEED_ALIGNMENT_FLAG=true;
NEED_ALIGNMENT_OUTPUT=false;
PROVIDED_ALIGNMENTS=[];

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case {'alignments'}
            if iscell(varargin{i+1})
                PROVIDED_ALIGNMENTS=[PROVIDED_ALIGNMENTS;varargin{i+1}];

            else
                error('FindmiRNAislands:BAD_ALIGNMENTS','Alignments must be a cell-array');
            end
        case {'alignment_prop', 'alignment prop'}
            if iscell(varargin{i+1})
                ALIGNMENT_PROP=varargin{i+1};
                try     %test with trivial example to make sure the input is valid
                    nwalign_mod('ACGT','AAGTC',ALIGNMENT_PROP{:});
                catch
                    warning('FindmiRNAislands:BAD_ALIGNMENT_PROP','Invalid nwalign arguements provided.');
                    rethrow(lasterror)
                end
            else
                error('FindmiRNAislands:BAD_ALIGNMENT_PROP','Invalid nwalign arguements provided.');
            end
        case {'alignment_output'}
            if islogical(varargin{i+1})
                NEED_ALIGNMENT_OUTPUT=varargin{i+1}&true;
            else
                error('FindmiRNAislands:BAD_ALIGNMENT_OUTPUT','Arguement to ALIGNMENT_OUTPUT must be a logical.');
            end



        case {'trust order','trust_order'}
            if islogical(varargin{i+1})&&varargin{i+1}
                NEED_ALIGNMENT_FLAG=false;
            else
                error('FindmiRNAislands:TRUST_ORDER_ARG','Arguement to TRUST_ORDER must be logical');
            end

        case {'fig_only','fig only','figure_only','figure only'}
            if islogical(varargin{i+1})
                FIGURES_TO_GENERATE=1:FIG_TOTAL*varargin{i+1};
            elseif isvector(varargin{i+1})&&all(varargin{i+1}>0&&varargin{i+1}<=FIG_TOTAL)
                FIGURES_TO_GENERATE=varargin{i+1};
            else
                error('FindmiRNAislands:BAD_FIG_ONLY_ARG','Arguement to FIG_ONLY must either be a logical scalar or a numerical vector [1 %d].',FIG_TOTAL);
            end

        case {'axes_handle','axes handle'}
            if ishandle(varargin{i+1})
                AX_HANDLE = varargin{i+1};
            else
                error('FindmiRNAislands:BAD_AXES_HANDLE','Arguement to AXES_HANDLE must be an axes handle.');
            end

        case {'excel_output'}
            if ischar(varargin{i+1})
                EXCEL_FILENAME=varargin{i+1};
            else
                error('FindmiRNAislands:BAD_EXCEL_FILENAME','Arguement to EXCEL_FILENAME must be a char-array.');
            end

        case {'alignments_only'}
            if islogical(varargin{i+1})
                ALIGN_ONLY_FLAG = varargin{i+1}&true;
            else
                error('FindmiRNAislands:BAD_ALIGNMENTS_ONLY','Arguement to ALIGNMENTS_ONLY must be a logical.');
            end
            
        case {'global_cutoff'}
            if isnumeric(varargin{i+1})&&varargin{i+1}<1
                GLOBAL_CUTOFF = varargin{i+1};
            else
                error('FindmiRNAislands:BAD_GLOBAL_CUTOFF','Arguement to GLOBAL_CUTOFF must be a numeric between [0,1].');
            end

        otherwise
            error('FindmiRNAislands:BAD_ARG','An unknown input arguements was provided: %s',varargin{i});
    end
end

if nargout == 2
    NEED_ALIGNMENT_OUTPUT=true;
end

%%%%%%%DONE ARGUEMENT CHECKING



%%%%%%%%%%%DO alignments if needed
if NEED_ALIGNMENT_FLAG
    ALIGNMENTS=cell(NUM_TEST_SEQS,1);

    %%%%%use provided alignments if possible
    if ~isempty(PROVIDED_ALIGNMENTS)
        provided_indexes=MatchAlignments(REF_SEQ_CELL,TEST_SEQS_CELL,PROVIDED_ALIGNMENTS);
        ALIGNMENTS(provided_indexes~=0)=PROVIDED_ALIGNMENTS(nonzeros(provided_indexes));
    end


    still_needed=cellfun('isempty',ALIGNMENTS);
    keepTrack = still_needed;
    DONE_COUNTER=0;
    tic
    while any(still_needed)
        display('Doing Alignments')

        needed_alignments = find(still_needed,PARFOR_BLOCK_SIZE)';
        tempAlignments = cell(min(PARFOR_BLOCK_SIZE,length(needed_alignments)),1);
        tempSeqs = TEST_SEQS_CELL(needed_alignments);

        parfor (LOOP_IND=1:min(PARFOR_BLOCK_SIZE,length(tempAlignments)))
            [junk tempAlignments{LOOP_IND}]=nwalign(REF_SEQ_CELL{1},tempSeqs{LOOP_IND},ALIGNMENT_PROP{:});
        end

        DONE_COUNTER = DONE_COUNTER+PARFOR_BLOCK_SIZE;

        still_needed(needed_alignments)=false;
        ALIGNMENTS(needed_alignments)=tempAlignments;


        time=(toc/DONE_COUNTER)*(nnz(still_needed))/60

        %%%%add the alignments preformed so they can be returned easily if
        %%%%desired

    end
    PROVIDED_ALIGNMENTS=[PROVIDED_ALIGNMENTS;ALIGNMENTS(keepTrack)];
else
    ALIGNMENTS=PROVIDED_ALIGNMENTS;
end

if ALIGN_ONLY_FLAG
    outputData = PROVIDED_ALIGNMENTS;
    return
end

snp_spot=false(NUM_TEST_SEQS,REF_SEQ_LENGTH);
for i=1:NUM_TEST_SEQS
    mask=isletter(ALIGNMENTS{i}(1,:));
    snp_spot(i,:)=ALIGNMENTS{i}(1,mask)~=ALIGNMENTS{i}(3,mask);
end

[exactforward exactreverse]=VFindRuns(~snp_spot);

windowCons=mean((exactforward>MIN_LENGTH));

validSpots=bsxfun(@ge,windowCons,CONS_CUTOFF(:));

[validForward validReverse]=VFindRuns(validSpots);



if any(FIGURES_TO_GENERATE==1)

    if isempty(AX_HANDLE)
        AX_HANDLE = axes;
    end

    axes(AX_HANDLE);

    image(validForward+validReverse)
    labels = arrayfun(@(x)(num2str(x)),CONS_CUTOFF,'uniformoutput',false);
    set(gca,'ytick',10*CONS_CUTOFF,'yticklabel',labels)
    ylabel('Conservation')
    xlabel('Genomic Position')

    colormap(flipud(gray))
end

headMask = diff([zeros(size(validForward,1),1) validForward],1,2)>0;
headSpots = ((validForward.*headMask)>MIN_LENGTH).*validForward;



%[text genome_pos length tested_cons true_cons]
outputData = cell(nnz(headSpots),5);

spots=find(headSpots);

for i = 1:length(spots)
    [I J]=ind2sub(size(headSpots),spots(i));

    outputData{i,1}=REF_SEQ_CELL{1}(J:J+headSpots(spots(i)));
    outputData{i,2}=J;
    outputData{i,3}=headSpots(spots(i));
    outputData{i,4}=CONS_CUTOFF(I);

    temp = strfind(TEST_SEQS_CELL,REF_SEQ_CELL{1}(J:J+headSpots(spots(i))));


    outputData{i,5}=nnz(~cellfun('isempty',temp))/length(TEST_SEQS_CELL);
end

mask = cellfun(@(x)(x>GLOBAL_CUTOFF),outputData(:,5));

outputData=outputData(mask,:);


if ~isempty(EXCEL_FILENAME)
    xlswrite(EXCEL_FILENAME,outputData);
end


if NEED_ALIGNMENT_OUTPUT
    varargout{1}=PROVIDED_ALIGNMENTS;
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






end


















