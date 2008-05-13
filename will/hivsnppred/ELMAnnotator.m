function PAT_STRUCT=ELMAnnotator(HIV_STRUCT,PAT_STRUCT,ELM_STRUCT)
%   ELMAnnotator
%       This function searches through all sequences in the Patient_Struct
%       to find ELM motiffs.  It then maps these onto positions on the
%       REF_SEQ.  Each patient is then annotated with a binary vector
%       indicating the presence or absence of each unique instance of each
%       motiff.
%
%       PAT_STRUCT=ELMAnnotator(PAT_STRUCT,ELM_STRUCT)
%
%       ELM_STRUCT      A structure containing the ELM regular expressions.
%
%
%
%       See also: ELMDownloadParser
%

%WAITBAR_HANDLE=waitbar(0,'Processing Inputs');

if ~isstruct(ELM_STRUCT)||~isfield(ELM_STRUCT,'REG_EXPR')
    error('ELMAnnotator:BAD_ELM_STRUCT','ELM_STRUCT must contain a REG_EXPR field as created by ELMDownloadParser.')
end


map=geneticcode;
trans_table=[[fieldnames(map) struct2cell(map)];{'NNN'  'X'}];

stop_indexes=find(strcmp('*',trans_table(:,2)));



[ALIGNMENT_CELLS ALIGNMENT_INDS]=PatientStructHelper(PAT_STRUCT,{'Alignment_CELL','explodeCells'});

%waitbar(0,WAITBAR_HANDLE,'Translating')
SEQS=cell(size(ALIGNMENT_CELLS));
MAPPING=cell(size(ALIGNMENT_CELLS));
%tic
for i=1:length(ALIGNMENT_CELLS)
    %time=(toc/i)*(length(ALIGNMENT_CELLS)-i);
    %waitbar(i/length(ALIGNMENT_CELLS),WAITBAR_HANDLE, ['Translating: ' num2str(time/60)])
    temp=ALIGNMENT_CELLS{i};
    [SEQS{i} MAPPING{i}]=HIVORFChecker(HIV_STRUCT,temp(3,isletter(temp(3,:))));
    %    [SEQS{i} MAPPING{i}]=Translate(temp(3,:),1,nnz(isletter(temp(3,:)))/3);
end


MATCH_SPOTS=ELMFinder(SEQS,ELM_STRUCT);

%%%%normalize spots to the reference sequence

for i=1:length(ELM_STRUCT)
    MATCH_SPOTS(:,i)=cellfun(@(x,y)(y(x)),MATCH_SPOTS(:,i),MAPPING,'uniformoutput',false);
end

ELMSimple = num2cell(~cellfun('isempty',MATCH_SPOTS),2);
ELMBinned = cell(length(SEQS),length(ELM_STRUCT));
ELMPWMScore = cell(length(SEQS),length(ELM_STRUCT));

%tic
for i = 1:length(ELM_STRUCT)
    %time=(toc/i)*(length(length(ELM_STRUCT))-i);
    %waitbar(i/length(length(ELM_STRUCT)),WAITBAR_HANDLE, ['Applying Bins: ' num2str(time/60)])
    if ~isempty(ELM_STRUCT(i).PosBins)
        ELMBinned(:,i) = cellfun(@(x)(CheckBins(ELM_STRUCT(i).PosBins,x)),MATCH_SPOTS(:,i),'uniformoutput',false);
    end
    if ~isempty(ELM_STRUCT(i).PosPWMs)
        PWMMat = NaN(size(SEQS,1),length(ELM_STRUCT(i).PosPWMs));
        
        for k = 1:length(ELM_STRUCT(i).PosPWMs)
            WindowedSeqs=cellfun(@(x,y)(GetSeqWindow(ELM_STRUCT(i).PosBins(:,k),x,y)),MAPPING,SEQS,'uniformoutput',false);
            
            PWMVals = NaN(size(SEQS,1),1);
            
            windowMask = ~cellfun('isempty',WindowedSeqs);
   
            if any(windowMask)
                determinedVals=PWMEvaluator(ELM_STRUCT(i).PosPWMs{k},WindowedSeqs(windowMask));
                PWMVals(windowMask)=cellfun(@max,determinedVals);
                PWMMat(:,k)=PWMVals;
            end
           
        end
        ELMPWMScore(:,i)=num2cell(PWMMat,2);

    end
end

ELMBinAnnot=cell(length(ELM_STRUCT),1);
ELMPWMAnnot=cell(length(ELM_STRUCT),1);

[ELMBins{1:length(ELM_STRUCT)}]=deal(ELM_STRUCT(:).PosBins);
for i=1:length(ELMBinAnnot)
    if ~isempty(ELM_STRUCT(i).PosBins)
        ELMBinAnnot{i} = [i*ones(1,size(ELMBins{i},2)); ELMBins{i}];
    end
    if ~isempty(ELM_STRUCT(i).PosPWMs)
        ELMPWMAnnot{i} = [i*ones(1,size(ELMBins{i},2)); ELMBins{i}];
    end
        
end

AllBins = cat(2,ELMBinAnnot{:});
PWMBins = cat(2,ELMPWMAnnot{:});

for i = 1:length(PAT_STRUCT)
    SeqInds = find(ALIGNMENT_INDS==i);
    
    if ~isempty(SeqInds)
        tempSimple = sum(cell2mat(ELMSimple(SeqInds)),1)>0;

        tempBinned = ELMBinned(SeqInds,:);
        finalBinned=[];
        for k=1:size(tempBinned,2)
            if ~isempty(tempBinned{1,k})
                finalBinned = [finalBinned sum(cat(1,tempBinned{:,k}),1)>0];
            end
        end
        
        tempPWMScore = ELMPWMScore(SeqInds,:);
        finalScore=[];
        for k=1:size(tempPWMScore,2)
            if ~isempty(tempPWMScore{1,k})
                finalScore = [finalScore max(cat(1,tempPWMScore{:,k}),[],1)];
            end
        end

        PAT_STRUCT{i}.ELM_simple = tempSimple;
        PAT_STRUCT{i}.ELM_vec = finalBinned;
        PAT_STRUCT{i}.ELM_PWM = finalScore;
        PAT_STRUCT{i}.ELM_annot = AllBins;
        PAT_STRUCT{i}.ELM_PWMannot = PWMBins;
        
    else
        PAT_STRUCT{i}=rmfield(PAT_STRUCT{i},{'ELM_simple','ELM_vec','ELM_annot'});
    end

end




%close(WAITBAR_HANDLE)

    function Vals=CheckBins(Bins,MatchedSpots)
        if isempty(MatchedSpots)
            Vals = false(1,size(Bins,2));
        else
            Vals = histc(MatchedSpots,[Bins(1,:) Bins(2,end)])>0;
            Vals = Vals(1:end-1);
        end
    end

    function thisSeq=GetSeqWindow(Bin,Mapping,Seq)
        TF = Mapping>Bin(1) & Mapping<Bin(2);
        thisSeq=Seq(TF);

    end

end