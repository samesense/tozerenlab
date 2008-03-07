function varargout=PatDBbaseCalls(REF_SEQ,PATDB_ARRAY,SNP_SPOTS,varargin)
%   PatDBbaseCalls
%       Aligns all TEST sequences to the REFERENCE sequence to find
%       conserved islands of sequence and SNPs within those islands.
%
%   [SNP_SPOTS BASE_CALLS]=PatDBbaseCalls(REF_SEQ,PATDB_ARRAY,SNP_SPOTS)
%
%       REF_SEQ         A defined reference sequence to compare against
%       TEST_SEQS       A set of test sequences
%
%       BASE_CALLS      A matrix with the base calls at each spot defined
%                       as a SNP by the provided criteria
%       SNP_SPOTS       A vector of the spots at which the BASE_CALLS were
%                       made: relative to the reference sequence
%
%    [SNP_SPOTS BASE_CALLS ALIGNMENTS]= ...
%                           PatDBbaseCalls(REF_SEQ,TEST_SEQS,SNP_SPOTS)
%
%       ALIGNMENTS      A Nx1 cell of the alignments between each test
%                       sample and the reference
%
%
%    [SNP_SPOTS BASE_CALLS ALIGNMENTS MUT_CALLS]=...
%                               PatDBbaseCalls(REF_SEQ,TEST_SEQS,SNP_SPOTS)
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
%   Optional Properties
%
%       ALIGNMENTS      alignment_cell
%                           A Nx1 cell of alignments to the reference
%                           sample
%       ALIGNMENT_PROP  A cell array of properties to pass to NWALIGN
%
%       GRAPH_DISPLAY   true|[false]
%                           Toggle whether to display a set of graphs
%                           describing the output
%
%       SNP_SPOTS       Nx1 array of spots to use as the SNP locations
%
%
%
%   [VALID_MASK]=PatDBbaseCalls(TEST_SEQS)
%       Since this algorithm requires that all sequences be composed of
%       only 'ACGT-' then any input sequences which do not conform are
%       removed in pre-processing steps.  This call returns a boolean mask
%       indicating which sequences are in the final output.
%
%
%
%   See also MakeSNPCalls.
%


%%%%%%%%%%%%%%%%%%%%%Arguement parsing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[INT_NT_MAP ALLOWABLE_CHARS]=AmbigiousNTMap;
warning('ON','PatDBbaseCalls:BAD_CHAR');

ALIGNMENT_PROP={'alphabet', 'nt'};
MakeSNPCall_PROPS={'fragments',true};

%DEAL WITH SPECIAL CASE
if nargin==1&&(nargout==1||nargout==0)
    CheckSeq=@(input)(all(ismember(input,ALLOWABLE_CHARS)));

    if isstruct(REF_SEQ)&&isfield(REF_SEQ,'Sequence')
        varargout{1}=arrayfun(@(x)CheckSeq(x.Sequence),REF_SEQ);
    elseif iscell(REF_SEQ)
        varargout{1}=cellfun(CheckSeq,REF_SEQ);
    else
        error('PatDBbaseCalls:BAD_MASK','TEST_SEQS must be a cell-array or a struct-array with a Sequence field');
    end
    return
end

%%%%%%%CHECK REF_SEQ%%%%%%%%%
if isstruct(REF_SEQ)&&isfield(REF_SEQ,'Sequence')&&numel(REF_SEQ)==1
    REF_SEQ_CELL={SeqPREP(REF_SEQ.Sequence)};
elseif iscell(REF_SEQ)&&numel(REF_SEQ)==1&&ischar(REF_SEQ{1})
    REF_SEQ_CELL={SeqPREP(REF_SEQ{1})};
elseif ischar(REF_SEQ)&&size(REF_SEQ,1)==1&&size(REF_SEQ,2)>1
    REF_SEQ_CELL={SeqPREP(REF_SEQ)};
else
    error('PatDBbaseCalls:REF_SEQ','REF_SEQ must be a cell, a char-array or a struct with a Sequence field');
end




%%%%%%%CHECK PATDB_ARRAY%%%%%%%%%
if iscell(PATDB_ARRAY)
    T_TEST_SEQS_CELL=cellfun(@GetSeqFromPatDB,PATDB_ARRAY(:),'uniformoutput',false);
else
    error('PatDBbaseCalls:PATDB_ARRAY','TEST_SEQS must be a cell array or a struct-array with a Sequence field');
end

counter=1;
TEST_SEQS_CELL=cell(10*length(T_TEST_SEQS_CELL),1);
PAT_IND_ARRAY=zeros(10*length(T_TEST_SEQS_CELL),1);
for i=1:length(T_TEST_SEQS_CELL);
    temp=T_TEST_SEQS_CELL{i};
    if ~isempty(temp)
        num_seqs=length(temp);
        TEST_SEQS_CELL(counter:counter+num_seqs-1)=temp(:);
        PAT_IND_ARRAY(counter:counter+num_seqs-1)=i;
        counter=counter+num_seqs;
    end
end

TEST_SEQS_CELL=TEST_SEQS_CELL(1:counter-1);
PAT_IND_ARRAY=PAT_IND_ARRAY(1:counter-1);

valid_patients=ismember(1:length(PATDB_ARRAY),PAT_IND_ARRAY);
if any(~valid_patients)
    display([int2str(nnz(~valid_patients)) ' patients did not have sequences.'])
end




NUM_TEST_SEQS=length(TEST_SEQS_CELL);
FRAGMENT_ARRAY=true(NUM_TEST_SEQS,1);

%%%%%CHECK SNP_SPOTS
if isnumeric(SNP_SPOTS)&&max(SNP_SPOTS)<=length(REF_SEQ_CELL{1})&&min(SNP_SPOTS)>=1
    SNP_SPOTS=sort(SNP_SPOTS(:));
else
    error('PatDBbaseCalls:SNP_SPOTS','SNP_SPOTS must be a numeric array of values between 1 and length(REF_SEQ) [%d].',length(REF_SEQ_CELL{1}));
end



if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch(lower(varargin{i}))
            case {'seq_only'}
                %deal with special case and return the sequences
                varargout{1}=TEST_SEQS_CELL;
                varargout{2}=FRAGMENT_ARRAY;
                varargout{3}=PAT_IND_ARRAY;
                return                
            
            case{'alignments'}
                if iscell(varargin{i+1})&&size(varargin{i+1},1)==NUM_TEST_SEQS
                    ALIGNMENTS=varargin{i+1};
                    cellfun(@ValidateAlign,ALIGNMENTS(:),TEST_SEQS_CELL(:));
                    NEED_ALIGNMENT_FLAG=false;
                    MakeSNPCall_PROPS=[MakeSNPCall_PROPS,'alignments',ALIGNMENTS];
                else
                    error('PatDBbaseCalls:BAD_ALIGNMENTS','Alignments must be a cell array.');
                end
            case {'alignment_prop', 'alignment prop'}
                if iscell(varargin{i+1})
                    ALIGNMENT_PROP=varargin{i+1};
                    try     %test with trivial example to make sure the input is valid
                        nwalign_mod('ACGT','AAGTC',ALIGNMENT_PROP{:})
                    catch
                        warning('MakeSNPCalls:BAD_ALIGNMENT_PROP','Invalid nwalign arguements provided.');
                        rethrow(lasterror)
                    end
                else
                    error('PatDBbaseCalls:BAD_ALIGNMENT_PROP','Invalid nwalign arguements provided.');
                end
            case {'graph_display', 'graph display'}
                if islogical(varargin{i+1})
                    DISPLAY_FIG_FLAG=varargin{i+1};
                    MakeSNPCall_PROPS=[MakeSNPCall_PROPS,'graph_display',DISPLAY_FIG_FLAG];
                else
                    error('PatDBbaseCalls:BAD_GRAPH_DISPLAY_ARG','Arguement to GRAPH_DISPLAY must be logical');
                end
            otherwise
                error('PatDBbaseCalls:BAD_ARG','A unknown arguement was provided: %s',varargin{i})
        end
    end
end

%%%%Determine desired output
switch nargout
    case 0
        SNP_SPOTS_OUT_FLAG=false;
        BASE_CALLS_OUT_FLAG=false;
        ALIGNMENT_OUT_FLAG=false;
    case 2
        SNP_SPOTS_OUT_FLAG=true;
        BASE_CALLS_OUT_FLAG=true;
        ALIGNMENT_OUT_FLAG=false;
    case 3
        SNP_SPOTS_OUT_FLAG=true;
        BASE_CALLS_OUT_FLAG=true;
        ALIGNMENT_OUT_FLAG=true;
    otherwise
        warning('MakeSNPCalls:BAD_OUT','Incorrect number of output arguements provided: %d',nargout);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%Now that arguements are set, throw everything into MakeSNPCalls

if ALIGNMENT_OUT_FLAG
    [SNP_SPOTS BASE_CALLS ALIGNMENTS]=MakeSNPCalls(REF_SEQ_CELL,TEST_SEQS_CELL,'fragments',FRAGMENT_ARRAY,'snp_spots',SNP_SPOTS,'alignment_prop',ALIGNMENT_PROP,MakeSNPCall_PROPS{:});
else
    [SNP_SPOTS BASE_CALLS]=MakeSNPCalls(REF_SEQ_CELL,TEST_SEQS_CELL,'fragments',FRAGMENT_ARRAY,'snp_spots',SNP_SPOTS,'alignment_prop',ALIGNMENT_PROP,MakeSNPCall_PROPS{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%SETUP OUTPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargout
    case 0
        varargout=[];
    case 2
        varargout{1}=SNP_SPOTS;
        varargout{2}=BASE_CALLS;
    case 3
        varargout{1}=SNP_SPOTS;
        varargout{2}=BASE_CALLS;
        varargout{3}=ALIGNMENTS;
    otherwise
        varargout=[];
end


%%%%%%%%%%%%%%%%%%%ARGUEMENT CHECKING FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%

    function stripped_seq=SeqPREP(input)
        input=upper(input);
        if any(~ismember(input,ALLOWABLE_CHARS))
            warning('PatDBbaseCalls:BAD_CHAR','All sequences must be composed of %s. Skipping the incorrect sequences and proceding',ALLOWABLE_CHARS);
            warning('OFF','PatDBbaseCalls:BAD_CHAR');
            stripped_seq=[];
        else
            stripped_seq=input(isletter(input));
        end
    end

    function ValidateAlign(align_input,test_seq)
        ref_align=SeqPREP(align_input(1,:));
        test_align=SeqPREP(align_input(3,:));

        if ~strcmp(ref_align,REF_SEQ_CELL{1})||~strcmp(test_align,test_seq)
            error('PatDBbaseCalls:BAD_ALIGN','ALIGNMENT must a cell array of the outputs of NWALIGN(REF_SEQ,TEST_SEQ)');
        end
    end

    function seq=GetSeqFromPatDB(input)
        if isstruct(input)&&isfield(input,'PR_seqs')&&isfield(input,'RT_seqs')
            t_seqs=[input.PR_seqs(:);input.RT_seqs(:)];
            
            
            ind=~cellfun('isempty',t_seqs);
            if any(ind)
                seq=cellfun(@SeqPREP,t_seqs(ind),'uniformoutput',false);
                seq=seq(~cellfun('isempty',seq));
            else
                seq={};
            end
        else
            error('PatDBbaseCalls:BAD_PAT_DB','The PAT_DB is in an invalid format.  See HIVPatDBLoad for specifications.')
        end
    end

    function [intmap charorder]=AmbigiousNTMap
        %The rows are in ACGT order and the columns are in charorder
        intmap=zeros(4,17);
        charorder=['ACGTRYKMSWBDHVX-.'];
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
        intmap(:,17)=   [0;     0;      0;      0];     %.
        
        
    end



end
