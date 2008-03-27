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
%   Optional Properities
%
%       WHOLE_SEQ       [true|FALSE] if set to true then the algorithm uses
%                       the entire sequence for calculating the PWM.
%                       Multiple RegExp matches will be maximized together.
%
%       DISPLAY         [TRUE|false] When set to true the algorithm shows
%                       the optimization at each iteration.
%
%       MAPPING_DATA    A cell-array of mapping data from HIVORFChecker.
%                       Beware, only minimal error-checking is performed,
%                       so be SURE that it is the right data.
%
%

MAPPING_DATA=[];
WHOLE_SEQ_FLAG=false;
DISPLAY_FLAG=true;


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
            
            
            otherwise
                error('ELMTrainPWM:BAD_ARG','An unknown arguement was provided: %s',varargin{i})
            
        end
        
    end
    
end

if DISPLAY_FLAG
    options=optimset('OutputFcn', @PWMDisplay);
end


[reg_exps{1:length(ELM_STRUCT)}]=deal(ELM_STRUCT(:).REG_EXPR);

INITIAL_PWMS=cellfun(@TransELMtoPWM,reg_exps,'uniformoutput',false);

skipped_mask=(cellfun('isempty',INITIAL_PWMS));
if nnz(skipped_mask)~=0
    display(['There are ' int2str(nnz(skipped_mask)) ' ELMs that will be skipped.'])
end


MatchedSpots=ELMFinder(INPUT_SEQS,ELM_STRUCT);
[PositCall AnnotVec Centers]=ELMPositional(INPUT_SEQS,MAPPING_DATA,ELM_STRUCT,'Matched_Spots',MatchedSpots);


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

        

        [x fval]=fmincon(@PWMObjective,OriginalPWM(ChangeSpots),[],[],Aeq,Beq,zeros(length(ChangeSpots),1),ones(length(ChangeSpots),1),[],options);
    end

else
    CenterInd=find(~cellfun('isempty',Centers));
    
    for i=1:length(CenterInd)
        this_ELM=CenterInd(i);

        if isempty(INITIAL_PWMS{this_ELM})
            continue
        end
        
        TheseCenters=Centers{CenterInd(i)}(1,:);
        
        CentroidWidth=round(2*min(TheseCenters));
        CentroidBound=1:CentroidWidth:max(SeqLengths);
        
        for g=1:length(CentroidBound)-1
            CurrentFval=zeros(1,1000);
            INT_SEQS=cellfun(@SeqRetrieve,INT_SEQS_MASTER,'uniformoutput',false);
            NormMatchedSpots=cellfun(@NormalizeSpots,MatchedSpots(:,this_ELM),'uniformoutput',false);
            SeqVals=cell(length(INT_SEQS),1);

            OriginalPWM=INITIAL_PWMS{this_ELM};
            OriginalPWM=repmat(1./sum(OriginalPWM,1),[20 1]).*OriginalPWM;
            ChangeSpots=find(OriginalPWM>0);
            

            [find_I find_J]=find(OriginalPWM>0);


            this_MatchedSpots=cellfun(@(x)ismember(1:CentroidWidth,x),NormMatchedSpots,'uniformoutput',false);
            

%             Aeq=zeros(length(ChangeSpots));
%             for j=1:length(ChangeSpots)
%                 Aeq(j,find_J(j)==find_J)=1;
%             end
%             Beq=ones(size(ChangeSpots));

            figure
            pwm_figure=subplot(3,2,1);
            fval_figure=subplot(3,2,2);
            score_figure=subplot(3,1,2);
            hist_figure=subplot(3,1,3);
            subplot(hist_figure), hist([NormMatchedSpots{:}],CentroidWidth);
            drawnow

            [x fval]=fmincon(@PWMObjective,OriginalPWM(ChangeSpots),[],[],[],[],zeros(length(ChangeSpots),1),[],[],options);
            %[x fval]=fminsearch(@PWMObjective,OriginalPWM(ChangeSpots),options);
        end
        
    end

end



    function ObjVal=PWMObjective(INPUT_VALS)
        
        CurrentPWM=OriginalPWM;
        CurrentPWM(ChangeSpots)=INPUT_VALS;
        %loop implemented in mexFunction PWMEvaluator
%         tic
%         for CELL_IND=1:length(INT_SEQS)
%             temp_seq=INT_SEQS{CELL_IND};
%             temp_mat=zeros(1,length(temp_seq));
%             for SEQ_IND=1:length(temp_mat)-PWMLength
%                 temp_mat(SEQ_IND)=sum(CurrentPWM(sub2ind(size(CurrentPWM),temp_seq(SEQ_IND:(SEQ_IND+PWMLength-1)),1:PWMLength)));
%             end
%             SeqVals{CELL_IND}=temp_mat;
%         end
%         toc
        
        SeqVals=PWMEvaluator(CurrentPWM,INT_SEQS);
        
        MatchedVals=cellfun(@(x,y)(x(y)),SeqVals,this_MatchedSpots,'uniformoutput',false);
        NonMatchedVals=cellfun(@(x,y)(x(~y)),SeqVals,this_MatchedSpots,'uniformoutput',false);
        
        ObjVal=mean([NonMatchedVals{:}])-mean([MatchedVals{:}]);
        
    end

    function stop=PWMDisplay(x,OptimStruct,state)
        stop=false;
        switch state
            case 'iter'
                this_PWM=OriginalPWM;
                this_PWM(ChangeSpots)=x;

                CurrentFval(OptimStruct.iteration+1)=OptimStruct.fval;
                
                subplot(pwm_figure), pcolor(this_PWM);
                subplot(fval_figure), plot(1:OptimStruct.iteration+1,CurrentFval(1:OptimStruct.iteration+1));
                
                
                noise_display_mat=zeros(length(INT_SEQS),CentroidWidth);
                signal_display_mat=NaN(length(INT_SEQS),CentroidWidth);
                
                for p=1:length(INT_SEQS)
                    noise_display_mat(p,1:CentroidWidth)=SeqVals{p}(1:CentroidWidth);
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
                
            case 'init'
                

                
               

            case 'done'

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
        SEQ=input_seq(CentroidBound(g):min(CentroidBound(g+1),length(input_seq)));
    end

    function NormedSpot=NormalizeSpots(match_spots)
        if isempty(match_spots)
            NormedSpot=[];
            return
        else
            NormedSpot=match_spots(match_spots>=CentroidBound(g)&match_spots<CentroidBound(g+1))-CentroidBound(g)+1;
        end
    end

end





