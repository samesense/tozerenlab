function varargout=PWMEvaluator(PWM_mat, INPUT_SEQS,varargin)
%   PWMEvaluator
%       Applies the Position Weight Matrix to every location in the
%       provided sequences.  This function can use either an optimized MEX
%       function PWMEvaluatorFAST or use a slower matlab implemented loop.
%
%   SeqVals=PWMEvaluator(PWM_mat, INPUT_SEQS)
%
%       PWM_mat         A position weight matrix for evaluation.  If
%                       size(PWM_mat,1)==4|5 then NT is assumed, if
%                       size(PWM_mat,1)==20|21 then AA is assumed.  This 
%                       can be over-ridden with optional parameters
%
%       INPUT_SEQS      A cell-array of sequences for evaluation.
%
%       SeqVals         A cell-array of doubles where X(i) is the value of
%                       the PWM on the window X(i:i+size(PWM,2))
%
%   SeqVals=PWMEvaluator(PWM_mat, INPUT_SEQS, 'TRUST_INPUT')
%
%       When called with this option then no input checking is performed
%       and the data is immediately sent to the MEX function.  INPUT_SEQS
%       must be a cell-array of DOUBLES (double(nt2int) or double(aa2int))
%       and that there are no Zero's or >size(PWM_mat,1).
%
%   Optional Properties
%
%       MISSING_VAL     Defines how to treat Gaps and Ambigious Characters.
%                       If set to NaN then any window that has a missing
%                       value will return NaN.  If set to Zero then the
%                       missing value will be ommitted.  DEFAULT=NaN.
%                       MISSING_VAL can also be a vector of window-position
%                       specific values.
%
%       SEQ_TYPE        If the default choice of 'AA' or 'NT' is not 
%                       possible then override with this option.  
%
%       MAPPING_FUN     A handle to mapping fun which can map abitrary the
%                       values in INPUT_SEQ to indexes in the ROWs of
%                       PWM_mat.
%
%
%       [newPWM_MAT newINPUT_SEQS]=PWMEvaluator(PWM_mat, INPUT_SEQS, ...
%           'SETUP_TRUST',true)
%
%       SETUP_TRUST     [true|FALSE] If set to TRUE then the function will
%                       return newPWM_MAT newINPUT_SEQS which can be
%                       re-used with the 'TRUST_INPUT' setting.
%
%


%deal with the special case without extra overhead
if length(varargin)==1&&strcmpi(varargin{1},'TRUST_INPUT')
    varargout{1}=PWMEvaluatorFAST(PWM_mat,INPUT_SEQS);
    return
end

SETUP_TRUST_FLAG=false;
MISSING_VAL=NaN;



if size(PWM_mat,1)==20||size(PWM_mat,1)==20
    MapFun=@(x)(double(aa2int(x)));
elseif size(PWM_mat,1)==4||size(PWM_mat,1)==5
    MapFun=@(x)(double(nt2int(x)));
else
    MapFun=@(x)(double(aa2int(x)));
end

if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case 'missing_val'
                if isnumeric(varargin{i+1})&&(numel(varargin{i+1})==1||size(varargin{i+1},2)==size(PWM_mat,2))
                    MISSING_VAL=varargin{i+1};
                else
                    error('PWMEvaluator:BAD_MISSING_VAL','Arguement for MISSING_VAL must either be a scalar or a vector the same size as PWM_mat');
                end
                
            case 'seq_type'
                switch lower(varargin{i+1})
                    case 'aa'
                        MapFun=@(x)(double(aa2int(x)));
                    case 'nt'
                        MapFun=@(x)(double(nt2int(x)));
                    otherwise
                        error('PWMEvaluator:BAD_SEQTYPE','An unkown arguement was provided to SEQ_TYPE: %s',varargin{i+1});
                end
            case 'mapping_fun'
                if ishandle(varargin{i+1})&&nargin(varargin{i+1})==1
                    MapFun=varargin{i+1};
                else
                    error('PWMEvaluator:BAD_MAPPINGFUN','Arguement for MAPPING_FUN be a function handle with 1 input and 1 output.');
                end
                
            case 'setup_trust'
                SETUP_TRUST_FLAG=true;

            otherwise
                error('PWMEvaluator:BAD_ARG','An unkown arguement was provided to PWMEvaluator: %s',varargin{i});
        end
    end
end

%Map sequences to Indexes into PWM_mat
NormalizedSeqs=cellfun(MapFun,INPUT_SEQS,'uniformoutput',false);


%check for values that are either Zero or >size(PWM_mat,2) because they
%will disrupt indexing within the MEX function.
MaxVal=max(cellfun(@max,NormalizedSeqs));
AnyZeros=any(cellfun(@(x)(any(x==0)),NormalizedSeqs));
if MaxVal>size(PWM_mat,1)||AnyZeros
    %add an extra column to the PWM_mat to give a value to the missing
    %values
    PWM_mat=[PWM_mat; zeros(1,size(PWM_mat,2))];
    PWM_mat(end,:)=MISSING_VAL;
    NormalizedSeqs=cellfun(@ReviseSeq,NormalizedSeqs,'uniformoutput',false);
end

if SETUP_TRUST_FLAG
    varargout{1}=PWM_mat;
    varargout{2}=NormalizedSeqs;
    return
end

varargout{1}=PWMEvaluatorFAST(PWM_mat,NormalizedSeqs);



    function input_seq=ReviseSeq(input_seq)
        %remove the bad values and replace them with an index into the last
        %column
        input_seq(input_seq==0|input_seq>size(PWM_mat,2))=size(PWM_mat,2);
    end
end
















