function ELM_POS_STRUCT=CreatePosAnnot(ELM_STRUCT,INPUT_SEQS,varargin)
%   CreatePosAnnot
%       This function modifies ELM_STRUCT to have positional information
%       for later use.  This includes the binning of each ELM and the PWMs
%       for each bin.
%
%
%   ELM_POS_STRUCT = CreatePosAnnot(ELM_STRUCT,INPUT_SEQS,...)
%
%       ELM_STRUCT          An ELM struct as created by ELMDownloadParser
%
%       INPUT_SEQS          A cell-array of sequences that will be fed to
%                           the Positional functions:
%                               ELMPositional
%                               ELMTrainPWM
%
%       ELM_POS_STRUCT      An ELM_STRUCT with the following feilds added
%                           to each element.
%
%       Optional Parameters:
%
%       ELMPositional_RESULTS   The results from a previous run of
%                               ELMPositional.  Including this will skip
%                               the step and add the results directly to
%                               the output.
%
%       ELMTrainPWM_RESUTLS     The results from a previous run of
%                               ELMTrainPWM.  Including this will skip
%                               the step and add the results directly to
%                               the output.
%
%
%
%
%


ELM_POS = [];
ELM_PWM = [];
MAPPING_CELL = [];


if ~isempty(varargin)
    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'elmpositional_results'
                if iscell(varargin{i+1})&&numel(varargin{i+1})==2
                    ELM_POS = varargin{i+1};
                else
                    error('CreatePosAnnot:BAD_POSRESULTS','The arguement of ELMPositional_RESULTS must be a 2x1 cell-array of POSTIONAL_CALL and ANNOT_VEC')
                end
            case 'elmtrainpwm_results'
                ELM_PWM = varargin{i+1};
                
            otherwise
                error('CreatePosAnnot:BAD_ARG','An unknown arguement was provided: %s',varargin{i})
        end
    end
end

if isempty(ELM_POS)
    ELM_POS = cell(2,1);
    
    [ELM_POS{1} ELM_POS{2}] = ELMPositional(INPUT_SEQS,MAPPING_CELL,ELM_STRUCT);

end

if isempty(ELM_PWM)
   
    ELM_PWM = ELMTrainPWM(INPUT_SEQS,ELM_STRUCT,...
        'POSITIONAL_CALL',ELM_POS{1},...
        'ANNOT_VECTOR',ELM_POS{2});

end

ELM_PWM_INDS = cell2mat(ELM_PWM(:,1));
ELM_PWM_INDS(2:2:end,:)=[];



ELM_POS_STRUCT = struct('Name',[],'REG_EXPR',[],'PosBins',[],'PosPWMs',[]);

for i = 1:length(ELM_STRUCT)
    thisELM = ELM_STRUCT(i);
    
    AnnotMask = find(ELM_POS{2}(1,:)==i);
    
    if any(AnnotMask)
        thisELM.PosBins = ELM_POS{2}(2:end,AnnotMask);
        try
        thisELM.PosPWMs = ELM_PWM(ELM_PWM_INDS(:,1)==i,2);
        catch
            123
        end

    else
        thisELM.PosBins=[];
        thisELM.PosPWMs=[];
        
    end
    
    ELM_POS_STRUCT(i) = thisELM;
end













