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

WAITBAR_HANDLE=waitbar(0,'Processing Inputs');

if ~isstruct(ELM_STRUCT)||~isfield(ELM_STRUCT,'REG_EXPR')
    error('ELMAnnotator:BAD_ELM_STRUCT','ELM_STRUCT must contain a REG_EXPR field as created by ELMDownloadParser.')
end


map=geneticcode;
trans_table=[[fieldnames(map) struct2cell(map)];{'NNN'  'X'}];

stop_indexes=find(strcmp('*',trans_table(:,2)));



[ALIGNMENT_CELLS ALIGNMENT_INDS]=PatientStructHelper(PAT_STRUCT,{'Alignment_CELL','explodeCells'});

waitbar(0,WAITBAR_HANDLE,'Translating')
SEQS=cell(size(ALIGNMENT_CELLS));
MAPPING=cell(size(ALIGNMENT_CELLS));
tic
for i=1:length(ALIGNMENT_CELLS)
    time=(toc/i)*(length(ALIGNMENT_CELLS)-i);
    waitbar(i/length(ALIGNMENT_CELLS),WAITBAR_HANDLE, ['Translating: ' num2str(time/60)])
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
ELMAnnot=cell(length(ELM_STRUCT),1);

tic
for i = 1:length(ELM_STRUCT)
    time=(toc/i)*(length(length(ELM_STRUCT))-i);
    waitbar(i/length(length(ELM_STRUCT)),WAITBAR_HANDLE, ['Applying Bins: ' num2str(time/60)])
    if ~isempty(ELM_STRUCT(i).PosBins)
        ELMBinned(:,i) = cellfun(@(x)(CheckBins(ELM_STRUCT(i).PosBins,x)),MATCH_SPOTS(:,i),'uniformoutput',false);
        
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

for i = 1:length(PAT_STRUCT)
    SeqInds = find(ALIGNMENT_INDS==i);
    
    if ~isempty(SeqInds)
        tempSimple = sum(cell2mat(ELMSimple(SeqInds)),1)>0;

        tempBinned = ELMBinned(SeqInds,:);
        tempBinned = cat(2,tempBinned{:});
        tempBinned = reshape(tempBinned,nnz(SeqInds),[]);
        tempBinned = sum(tempBinned,1)>0;

        tempPWMScore = ELMPWMScore(SeqInds,:);
        tempPWMScore = cat(2,tempPWMScore{:});
        tempPWMScore = reshape(tempPWMScore,nnz(SeqInds),[]);
        tempPWMScore = max(tempPWMScore,[],1)>0;

        PAT_STRUCT{i}.ELM_simple = tempSimple;
        PAT_STRUCT{i}.ELM_vec = tempBinned;
        PAT_STRUCT{i}.ELM_PWM = tempPWMScore;
        
    else
        PAT_STRUCT{i}=rmfield(PAT_STRUCT{i},{'ELM_simple','ELM_vec','ELM_annot'});
    end
    
    
    
    
end




close(WAITBAR_HANDLE)

    function Vals=CheckBins(Bins,MatchedSpots)
        if isempty(MatchedSpots)
            Vals = false(1,size(Bins,2));
        else
            Vals = histc(MatchedSpots,[Bins(1,:) Bins(2,end)])>0;
        end
    end

    function thisSeq=GetSeqWindow(Bin,Mapping,Seq)
        TF = Mapping>Bin(1) & Mapping<Bin(2);
        thisSeq=Seq(TF);

    end

end