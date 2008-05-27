function PWM_CELL=ELMTrainPWM(INPUT_SEQS,ELM_STRUCT,varargin)
%   ELMTrainPWM
%       This creates a set of Position-Weight-Matrices (PWMs) for each ELM
%       such that the difference in score between the known Match Spots and
%       the background is maximized. As described:
%           iff I = [all spots where the RegExp matches]
%           iff J = [all other locations]
%           iff PWM(X) returns the value at spot X
%       min(mean(PWM(J))-mean(PWM(I)))
%
%
%   PWM_CELL=ELMTrainPWM(INPUT_SEQS,ELM_STRUCT)
%
%       INPUT_SEQS      The training AA sequences.
%
%       ELM_STRUCT      A Struct describing the ELMs as created by
%                       ELMDownloadParser.
%
%       PWM_CELL        A cell-array of Position Weight Matrices.  It will
%                       be in the same order as returned by AnnotVector
%                       from ELMPositional.
%
%
%
%   Optional Inputs
%   
%       MAPPING_DATA    A cell-array of mapping data from HIVORFChecker.
%                       Beware, only minimal error-checking is performed,
%                       so be SURE that it is the right data.
%       
%       MATCH_SPOTS     The output of ELMFinder.  This can drastically
%                       speed up computation.
%
%       POSITION_CALL   The first output of ELMPositional.
%
%       ANNOT_VECTOR    The second output of ELMPositional.
%
%
%   Optional Properities
%
%       WHOLE_SEQ       [true|FALSE] if set to true then the algorithm uses
%                       the entire sequence for calculating the PWM.
%                       Multiple RegExp matches will be maximized together.
%
%       DISPLAY         [TRUE|false] When set to true the algorithm shows
%                       the optimization at each iteration.
%
%
%

GetMatched=@(x,y)(x(y));
GetNotMatched=@(x,y)(x(~y));



MAPPING_DATA=[];
WHOLE_SEQ_FLAG=false;
DISPLAY_FLAG=true;
MatchedSpots=[];
PositCall=[];
AnnotVec=[];
MaxIter=500;
FILENAME=[];

if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            
            case 'whole_seq'
                if islogical(varargin{i+1})
                    WHOLE_SEQ_FLAG=varargin{i+1};
                else
                    error('ELMTrainPWM:BAD_WHOLESEQ','The arguement provided to WHOLE_SEQ must be a logical')
                end
                
            case 'display'
                if islogical(varargin{i+1})
                    DISPLAY_FLAG=varargin{i+1};
                else
                    error('ELMTrainPWM:BAD_DISPLAY','The arguement provided to DISPLAY must be a logical')
                end
                    
            case 'mapping_data'
                if iscell(varargin{i+1})&&numel(varargin{i+1})==numel(INPUT_SEQS)
                    MAPPING_DATA=varargin{i+1};
                    
                else
                    error('ELMTrainPWM:BAD_MAPPINGDATA','The arguement provided to MAPPING_DATA must be a cell-array the same size as INPUT_SEQS')
                end
            
            case 'match_spots'
                if iscell(varargin{i+1})&&length(INPUT_SEQS)==size(varargin{i+1},1)&&length(ELM_STRUCT)==size(varargin{i+1},2)
                    MatchedSpots=varargin{i+1};
                else
                    error('ELMTrainPWM:BAD_MATCHSPOTS','The arguement provided to MATCH_SPOTS must be a cell-array with size of: [%d %d], was provided: [%d %d]',...
                        length(INPUT_SEQS),length(ELM_STRUCT),size(varargin{i+1},1),size(varargin{i+1},2));
                end
            
            case 'position_call'
                if (isnumeric(varargin{i+1})||islogical(varargin{i+1}))&&size(varargin{i+1},1)==length(INPUT_SEQS)
                    PositCall=varargin{i+1};
                else
                    error('ELMTrainPWM:BAD_POSITIONCALL','The arguement provided to POSITION_CALL must be a logical-array of which has length(INPUT_SEQS) rows');
                end
                
            case 'annot_vector'
                if isnumeric(varargin{i+1})&&size(varargin{i+1},1)==3
                    AnnotVec=varargin{i+1};
                else
                    error('ELMTrainPWM:BAD_ANNOTVECTOR','The arguement provided to ANNOT_VECTOR must be a numeric-array of which has 3 rows');
                end
                
            case 'from_file'
                if ischar(varargin{i+1})
                    FILENAME=varargin{i+1};
                    
                else
                    error('ELMTrainPWM:BAD_FROMFILE','The arguement provided to FROM_FILE must be a char-array.');
                    
                end

               
            otherwise
                error('ELMTrainPWM:BAD_ARG','An unknown arguement was provided: %s',varargin{i})
            
        end
        
    end
    if (~isempty(PositCall)&&~isempty(AnnotVec))&&size(PositCall,2)~=size(AnnotVec,2)
        error('ELMTrainPWM:BAD_ANNOTVEC_POSITCALL_SIZE','The number of columns in POSITION_CALL must match the columns in ANNOT_VECTOR');
    end
    
    
    
end

if DISPLAY_FLAG
    options=optimset('OutputFcn', @PWMDisplay,'MaxIter',MaxIter);
else
    options=optimset('MaxIter',MaxIter);
end


[reg_exps{1:length(ELM_STRUCT)}]=deal(ELM_STRUCT(:).REG_EXPR);

INITIAL_PWMS=cellfun(@TransELMtoPWM,reg_exps,'uniformoutput',false);

skipped_mask=(cellfun('isempty',INITIAL_PWMS));
if nnz(skipped_mask)~=0
    display(['There are ' int2str(nnz(skipped_mask)) ' ELMs that will be skipped.'])
end

if isempty(MatchedSpots)
    MatchedSpots=ELMFinder(INPUT_SEQS,ELM_STRUCT);
end

if isempty(FILENAME)
    if isempty(AnnotVec)
        [PositCall AnnotVec]=ELMPositional(INPUT_SEQS,MAPPING_DATA,ELM_STRUCT,'Matched_Spots',MatchedSpots);
    end
    
else
    fid=fopen(FILENAME,'r');
    binning=cell(0,2);
    temp1=fgetl(fid);
    temp2=fgetl(fid);
    while temp1~=-1
        temp_bin = cell(1,2);
        temp_bin{2} = str2num(temp2);
        temp_bin{1} = temp1;
        
        if length(temp_bin{2}) > 2
            binning=[binning; temp_bin];
        end

        temp1=fgetl(fid);
        temp2=fgetl(fid);
    end
    fclose(fid);
    
    [correctOrdering{1:length(ELM_STRUCT)}] = deal(ELM_STRUCT(:).Name);
    
    [TF ordering]=ismember(binning(:,1),correctOrdering);
    
    AnnotVec = zeros(3,0);
    
    for i=1:length(ordering)
        for k=1:length(binning{i,2})-1
            thisSlice=zeros(3,1);
            thisSlice(1)=ordering(i);
            thisSlice(2)=binning{i,2}(k)+1;
            thisSlice(3)=binning{i,2}(k+1)+1;
            AnnotVec = [AnnotVec thisSlice];
        end
    end

end

INT_SEQS_MASTER=cellfun(@(x)(double(aa2int(x))),INPUT_SEQS,'uniformoutput',false);
SeqLengths=cellfun('length',INT_SEQS_MASTER);
SeqLengths_CELL=cellfun(@length,INT_SEQS_MASTER,'uniformoutput',false);


if WHOLE_SEQ_FLAG
    INT_SEQS=INT_SEQS_MASTER;

    for i=1:length(AnnotVec)
        CurrentFval=zeros(1,1000);
        SeqVals=cell(length(INT_SEQS),1);
        OriginalPWM=INITIAL_PWMS{AnnotVec(2,i)};
        if isempty(OriginalPWM)
            continue
        end

        PWMLength=size(OriginalPWM,2);
        OriginalPWM=repmat(1./sum(OriginalPWM,1),[20 1]).*OriginalPWM;
        ChangeSpots=find(OriginalPWM>0&OriginalPWM<1);
        [find_I find_J]=find(OriginalPWM>0&OriginalPWM<1);

        this_MatchedSpots=cellfun(@(x,y)(ismember(1:y,x)),MatchedSpots(:,AnnotVec(2,i)),SeqLengths_CELL,'uniformoutput',false);

        Aeq=zeros(length(ChangeSpots));
        for j=1:length(ChangeSpots)
            Aeq(j,find_J(j)==find_J)=1;
        end
        Beq=ones(size(ChangeSpots));

        figure
        pwm_figure=subplot(2,2,1);
        fval_figure=subplot(2,2,2);
        score_figure=subplot(2,1,2);
        drawnow

        
        [x fval]=fmincon(@PWMObjective,OriginalPWM(ChangeSpots),[],[],[],[],zeros(length(ChangeSpots),1),[],[],options);
        
    end

else
    FoundELMs=unique(AnnotVec(1,:));
    counter=0;
    PWM_CELL = cell(0,2);
    for i=1:length(FoundELMs)
        this_ELM=FoundELMs(i);

        if isempty(INITIAL_PWMS{this_ELM})
            %could not convert the reg-exp to PWM, nothing to do
            continue
        end
        
        TheseBins=AnnotVec(2:3,AnnotVec(1,:)==this_ELM);
        
        for g=1:size(TheseBins,2)
            if TheseBins(1,g)>TheseBins(2,g)
                continue
            end
            CurrentFval=zeros(1,1000);
            INT_SEQS=cellfun(@SeqRetrieve,INT_SEQS_MASTER,'uniformoutput',false);
            

            
            NormMatchedSpots=cellfun(@NormalizeSpots,MatchedSpots(:,this_ELM),'uniformoutput',false);
            SeqVals=cell(length(INT_SEQS),1);
            SeqLengths=cellfun('length',INT_SEQS);
            
            OriginalPWM=INITIAL_PWMS{this_ELM};
            OriginalPWM=repmat(1./sum(OriginalPWM,1),[20 1]).*OriginalPWM;

            [OriginalPWM INT_SEQS]=PWMEvaluator(OriginalPWM, INT_SEQS,'SETUP_TRUST',true);
            
            ChangeSpots=find(OriginalPWM>0);

            [find_I find_J]=find(OriginalPWM>0);

            this_MatchedSpots=cellfun(@(x)ismember(1:(TheseBins(2,g)-TheseBins(1,g)+1),x),NormMatchedSpots,'uniformoutput',false);
            this_MatchedSpots=cellfun(@(x,y)(x(1:length(y))),this_MatchedSpots,INT_SEQS,'uniformoutput',false);
            
            pwm_figure=subplot(2,2,1);
            fval_figure=subplot(2,2,2);
            score_figure=subplot(2,1,2);
            
            drawnow

            [x fval]=fminsearchbnd(@PWMObjective,OriginalPWM(ChangeSpots),zeros(length(ChangeSpots),1),10*ones(length(ChangeSpots),1),options);
            counter=counter+1;
            ThisPWM = OriginalPWM;
            ThisPWM(ChangeSpots) = x
            PWM_CELL{counter,2}=ThisPWM;
            PWM_CELL{counter,1}=[this_ELM*ones(2,1) TheseBins(:,g)];
            %[x fval]=fminsearch(@PWMObjective,OriginalPWM(ChangeSpots),options);
        end
        
    end

end



    function ObjVal=PWMObjective(INPUT_VALS)
        
        CurrentPWM=OriginalPWM;
        CurrentPWM(ChangeSpots)=INPUT_VALS;
        
        SeqVals=PWMEvaluator(CurrentPWM,INT_SEQS,'TRUST_INPUT');

        MatchedVals=cellfun(GetMatched,SeqVals,this_MatchedSpots,'uniformoutput',false);
        NonMatchedVals=cellfun(GetNotMatched,SeqVals,this_MatchedSpots,'uniformoutput',false);

        
        ObjVal=nanmean([NonMatchedVals{:}])/nanmean([MatchedVals{:}]);
        
    end

    function stop=PWMDisplay(x,OptimStruct,state)
        stop=false;
        switch state
            case 'done'
                if OptimStruct.fval<0.01
                    stop=true;
                end
                %CurrentFval(OptimStruct.iteration+1)=OptimStruct.fval;
                this_PWM=OriginalPWM;
                this_PWM(ChangeSpots)=x;

                subplot(pwm_figure), pcolor(this_PWM);
                subplot(fval_figure), plot(1:OptimStruct.iteration+1,CurrentFval(1:OptimStruct.iteration+1));              
                
                noise_display_mat=zeros(length(INT_SEQS),max(SeqLengths));
                signal_display_mat=NaN(length(INT_SEQS),max(SeqLengths));
                
                for p=1:length(INT_SEQS)
                    noise_display_mat(p,1:SeqLengths(p))=SeqVals{p}(:);
                    if ~isempty(NormMatchedSpots{p})
                        noise_display_mat(p,NormMatchedSpots{p})=NaN;
                        signal_display_mat(p,NormMatchedSpots{p})=SeqVals{p}(NormMatchedSpots{p});
                    end
                end
                
                noise_data=nanmean(noise_display_mat,1);
                noise_data(isnan(noise_data))=0;
                
                signal_data=nanmean(signal_display_mat,1);
                signal_data(isnan(signal_data))=0;

                subplot(score_figure), plot(1:length(noise_data),noise_data,'k',1:length(signal_data),signal_data,'b');
                
                drawnow
            case 'interupt'
            case  'init'
                
            case 'iter'
                CurrentFval(OptimStruct.iteration+1)=OptimStruct.fval;
        end
    end

    function PWM_mat=TransELMtoPWM(reg_exp)
        initial_mat=NaN(20,length(reg_exp)+1);
        
        col_counter=1;
        not_flag=false;
        inside_bracket_flag=false;
        try
            for k=1:length(reg_exp)
                switch reg_exp(k)
                    case '.'
                        initial_mat(:,col_counter)=1;
                        col_counter=col_counter+1;
                    case {'[','('}
                        inside_bracket_flag=true;
                    case {']',')'}
                        inside_bracket_flag=false;
                        col_counter=col_counter+1;
                        not_flag=false;
                    case '|'
                        col_counter=col_counter-1;
                    case '^'
                        not_flag=true;
                        initial_mat(:,col_counter)=1;

                    case '$'
                    otherwise
                        if isletter(reg_exp(k))
                            initial_mat(aa2int(reg_exp(k)),col_counter)=1&~not_flag;
                            col_counter=col_counter+double(1&~inside_bracket_flag);
                        else
                            PWM_mat=[];
                            return
                        end
                end
            end

            PWM_mat=initial_mat(:,1:find(mean(isnan(initial_mat),1)==1,1)-1);
            PWM_mat(isnan(PWM_mat))=0;
        catch
            PWM_mat=[];
        end

    end

    function SEQ=SeqRetrieve(input_seq)
        SEQ=input_seq(TheseBins(1,g):min(TheseBins(2,g),length(input_seq)));
    end

    function NormedSpot=NormalizeSpots(match_spots)
        if isempty(match_spots)
            NormedSpot=[];
            return
        else
            NormedSpot=match_spots(match_spots>=TheseBins(1,g)&match_spots<TheseBins(2,g))-TheseBins(1,g)+1;
        end
    end

end





