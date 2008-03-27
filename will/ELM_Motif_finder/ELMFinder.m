function varargout=ELMFinder(SEQS,ELM_STRUCT,varargin)
%   ELMFinder
%       Searches through a set of sequences for the ELM Motiffs provided.
%       Use ELMDownloadParser to create a STRUCT suitable for input.
%
%
%   MATCH_SPOTS=ELMFinder(SEQS,ELM_STRUCT)
%
%       SEQS        An M x 1 set of search sequences.  These can be 
%                   char-arrays, cells or structures with a Sequence-field.
%
%       ELM_STRUCT  An N x 1 structure of ELM regular expressions as 
%                   created by ELMDownloadParser.
%
%       MATCH_SPOTS An M x N cell-array where each cell contains the
%                   matched locations of the ELM motif N in the Sequence M.
%
%       [MATCH_SPOTS MATCH_SEQ]=ELMFinder(SEQS,ELM_STRUCT)
%
%       MATCH_SEQ   The matched sequence which corresponds to a particular
%                   ELM motif.
%
%
%
%       Optional Properties
%
%       NESTED      [true|FALSE].  If set to true then the waitbar is
%                   disabled so it will play nice with as a nested
%                   function.
%

%
%   See also: ELMDownloadParser.
%

NESTED_FLAG=false;

if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            
            case 'nested'
                if islogical(varargin{i+1})
                    NESTED_FLAG=true&varargin{i+1};
                else
                    error('ELMFinder:BAD_NESTED','Arguement to NESTED must be a logical.')
                end
        
            otherwise
                error('ELMFinder:BAD_ARG','Unrecognized arguement provided: %s.',varargin{i})
        end
        
    end
end




%%%%%%%%%%%%%%%%ERROR CHECKING
if ~NESTED_FLAG
    WAITBAR_HANDLE=waitbar(0,'Processing Inputs');
end

if isstruct(ELM_STRUCT)&&isfield(ELM_STRUCT,'REG_EXPR')
    [REG_EXPRS{1:length(ELM_STRUCT)}]=deal(ELM_STRUCT.REG_EXPR);
else
    error('ELMFinder:BAD_STRUCT','ELM_STRUCT must be a STRUCT as created by ELMDownloadParser.')
end

if isstruct(SEQS)&&isfield(SEQS,'Sequence')
    [SEQUENCES{1:length(SEQS)}]=deal(SEQS.Sequence);
elseif iscell(SEQS)
    LOCS=~cellfun('isempty',SEQS);
    SEQUENCES=cell(size(SEQS));
    SEQUENCES(LOCS)=SEQS(LOCS); %#ok<FNDSB>
end

if nargout==2
    BOTH_FLAG=true;
else
    BOTH_FLAG=false;
end

%%%%%%%%%%%%%%%%%DONE ERROR CHEKCING

MATCH_SPOTS=cell(length(SEQUENCES),length(REG_EXPRS));

if BOTH_FLAG
    MATCH_SEQ=cell(length(SEQUENCES),length(REG_EXPRS));
end

for i=1:length(REG_EXPRS)
    if ~NESTED_FLAG
        waitbar(i/length(REG_EXPRS),WAITBAR_HANDLE,['Processing: ' ELM_STRUCT(i).Name])
    end
    if BOTH_FLAG
        [MATCH_SPOTS(LOCS,i) MATCH_SEQ(LOCS,i)]=regexpi(SEQUENCES(LOCS),REG_EXPRS{i},'start','match');
    else
        MATCH_SPOTS(LOCS,i)=regexpi(SEQUENCES(LOCS),REG_EXPRS{i});
    end
end

   
    
 if BOTH_FLAG
    varargout=cell(2,1);
    varargout{1}=MATCH_SPOTS;
    varargout{2}=MATCH_SEQ;
else
    varargout{1}=MATCH_SPOTS;
end


if ~NESTED_FLAG
    close(WAITBAR_HANDLE)
end











end
