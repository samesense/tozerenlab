function varargout=DNAPromoterMatcher(SEARCH_SEQ,SEQ_DB_FILENAME,NUM_MISMATCH,EXCEL_FILENAME,varargin)
%   DNAPromoterMatcher
%       Searches through a database of sequences to find the SEARCH_SEQ.
%       This is optimized to find exact (or near exact) matches throughout
%       the promoter-eome.  If you are looking for the closest match of
%       SEARCH_SEQ in a small subset of promoters use ClosestDNAMatch.
%
%   [ GeneSymbol LLID POS STRAND ORIENT SEQ] =...
%               DNAPromoterMatcher(SEARCH_SEQ,SEQ_DB_FILENAME,NUM_MISMATCH)
%
%   SEARCH_SEQ          The querry sequence (can be a regular-expression).
%
%   SEQ_DB_FILENAME     The filename of the sequences.  These must be in
%                       the format provided by PAINT's Upstreamer Function.
%                       http://www.dbi.tju.edu/dbi/tools/paint/
%
%   NUM_MISMATCH        The number of allowable mismatches.
%
%   LLIDS               The EntrezIDS of the matches.
%
%   NUM_MIS             The number of missed base-pairs.
%
%   POS                 The number of base-pairs upstream of the
%                       translation start site.
%
%   STRAND              Whether the match is on the Sense or Anti-sense
%                       strand.
%
%   ORIENTATION         Whether the match is from 3'-5' or 5'-3'.
%
%   SEQ                 The matched GENOMIC sequence.
%
%
%   ...=DNAPromoterMatcher(...,EXCEL_FILENAME)
%
%   EXCEL_FILENAME      A filename to output the data in EXCEL format.
%
%   
%
%
%   See also: ClosestDNAMatch.
%
%

% SEARCH_SEQ ='AA[AT]ACAA[AT]TAA[AT]';
% SEQ_DB_FILENAME='unique_IDS.fa';

WAITBAR_HANDLE=waitbar(0,'Loading Sequences');

%%%%%%%%%%%%%%%INPUT CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GLOBAL_LENGTH=2000;

if nargin==3
    EXCEL_FILENAME=[];
end

try
    [headers seqs]=fastaread(SEQ_DB_FILENAME);
    LLIDS=cell(size(headers));
    GeneNames=cell(size(headers));
    for i=1:length(headers)
        temp=textscan(headers{i},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
        GeneNames(i)=temp{1};
        LLIDS(i)=temp{2};
    end
    
    SeqLength=cellfun('length',seqs);
    seqs=arrayfun(@(x,y)(x{1}(max(y-GLOBAL_LENGTH,1):end)),seqs,SeqLength,'uniformoutput',false);
    
    SeqLength=cellfun('length',seqs);
    
catch
    warning('DNAPromoterMatcher:BAD_SEQ_DB','Cannot read Sequence Database.')
    close(WAITBAR_HANDLE);
    rethrow(lasterror)
end

%%shorten database to 2000bp.

%%%Parse the Regular-Expression into sets of cells

counter=1;
cell_counter=1;
RegCell=cell(nnz(isletter(SEARCH_SEQ)),1);
INSIDE_TYPE_FLAG=false;

while(counter<=length(SEARCH_SEQ))
    if SEARCH_SEQ(counter)=='['
        INSIDE_TYPE_FLAG=true;
        counter=counter+1;
    elseif SEARCH_SEQ(counter)==']'
        INSIDE_TYPE_FLAG=false;
        counter=counter+1;
        cell_counter=cell_counter+1;
    elseif isletter(SEARCH_SEQ(counter))
        RegCell{cell_counter}=[RegCell{cell_counter} SEARCH_SEQ(counter)];
        cell_counter=cell_counter+(1*~INSIDE_TYPE_FLAG);
        counter=counter+1;
    end
end

RegCell=RegCell(1:cell_counter-1)';
RegCell=cellfun(@(x)(['[' x ']']),RegCell,'uniformoutput',false);

if NUM_MISMATCH==2
    MisMatchInds=zeros(1,length(RegCell));
    for i=1:length(RegCell)
        for j=i:length(RegCell)
            temp=zeros(1,length(RegCell));
            temp([i j])=1;
            MisMatchInds=[MisMatchInds; temp];
        end
    end
    [Y I]=sort(sum(MisMatchInds,2),'descend');
    MisMatchInds=MisMatchInds(I,:);
elseif NUM_MISMATCH==1
    MisMatchInds=[zeros(1,length(RegCell)); eye(length(RegCell))];
elseif NUM_MISMATCH==0
    MisMatchInds=zeros(1,length(RegCell));
end




waitbar(0,WAITBAR_HANDLE,'Creating Sequence-Cells')
PromoterSeqs=[seqs'; seqreverse(seqs)'; ...
    cellfun(@(x)(seqcomplement(x)),seqs','uniformoutput',false); ...
    cellfun(@(x)(seqrcomplement(x)),seqs','uniformoutput',false)];
PromoterType=[(1:length(seqs))' ones(length(seqs),1); (1:length(seqs))' 2*ones(length(seqs),1); ...
    (1:length(seqs))' 3*ones(length(seqs),1); (1:length(seqs))' 4*ones(length(seqs),1)];

waitbar(0,WAITBAR_HANDLE,'Searching Seqs')
MatchCell=cell(length(PromoterSeqs),length(MisMatchInds));
PosCell=cell(length(PromoterSeqs),length(MisMatchInds));
for i=1:size(MisMatchInds,1)
    waitbar(i/size(MisMatchInds,1))
    fprintf('%d of %d\n',i,size(MisMatchInds,1))
    temp=RegCell;
    change_inds=find(MisMatchInds(i,:));
    if ~isempty(change_inds)
        for j=change_inds
            temp{j}=[temp{j}(1) '^' temp{j}(2:end)];
        end
        RegSearch=cat(2,temp{:});
    else
        RegSearch=cat(2,temp{:});
    end
   
    [PosCell(:,i) MatchCell(:,i)]=regexpi(PromoterSeqs,RegSearch,'start','match');
end

%%%%%Explodes cells in which there are multiple matches for each
%%%%%MisMatchInds
counter=1;
while counter<=size(PosCell,2)
    NeedExpand=cellfun('length',PosCell(:,counter));
    if any(NeedExpand>1)
        NewTempPos=cellfun(@(x)(x(end)),PosCell(NeedExpand>1,counter),'uniformoutput',false);
        ShortenedPos=cellfun(@(x)(x(1:end-1)),PosCell(NeedExpand>1,counter),'uniformoutput',false);
        
        NewTempMatch=cellfun(@(x)(x(end)),MatchCell(NeedExpand>1,counter),'uniformoutput',false);
        ShortenedMatch=cellfun(@(x)(x(1:end-1)),MatchCell(NeedExpand>1,counter),'uniformoutput',false);
        
        
        PosCell(NeedExpand>1,counter)=ShortenedPos;
        MatchCell(NeedExpand>1,counter)=ShortenedMatch;
        
        if counter~=size(PosCell,2)
            PosCell=[PosCell(:,1:counter) cell(size(PosCell,1),1) PosCell(:,counter+1:end)];
            PosCell(NeedExpand>1,counter+1)=NewTempPos;
            
            MatchCell=[MatchCell(:,1:counter) cell(size(MatchCell,1),1) MatchCell(:,counter+1:end)];
            MatchCell(NeedExpand>1,counter+1)=NewTempMatch;
            
            MisMatchInds=[MisMatchInds(1:counter,:); MisMatchInds(counter,:); MisMatchInds(counter+1:end,:)];

        else
            PosCell=[PosCell(:,1:counter) cell(size(PosCell,1),1)];
            PosCell(NeedExpand>1,counter+1)=NewTempPos;
            
            MatchCell=[MatchCell(:,1:counter) cell(size(MatchCell,1),1)];
            MatchCell(NeedExpand>1,counter+1)=NewTempMatch;
            MisMatchInds=[MisMatchInds(1:counter,:); MisMatchInds(counter,:)];        
        end
        
    else
        counter=counter+1;
    end
    
end


FoundMask=~cellfun('isempty',PosCell(:,sum(MisMatchInds,2)<3));
[I J]=find(FoundMask);

waitbar(0,WAITBAR_HANDLE,'Generating Output')
UniGene=PromoterType(I);

UniCols=J;
UniRows=I;

%%%%%%%[GeneSymbol LLID POS STRAND ORIENT SEQ NUM_MISS]
ExcelOutput=cell(length(UniGene),7);
for i=1:length(UniGene)
    ExcelOutput(i,1)=GeneNames(UniGene(i));
    ExcelOutput(i,2)=LLIDS(UniGene(i));
    
    switch ceil(UniRows(i)/length(seqs))
        case 1
            ExcelOutput{i,3}=SeqLength(UniGene(i))-PosCell{UniRows(i),UniCols(i)};
            ExcelOutput{i,4}='Sense';
            ExcelOutput{i,5}='3prime - 5prime';
        case 2
            ExcelOutput{i,3}=PosCell{UniRows(i),UniCols(i)};
            ExcelOutput{i,4}='Sense';
            ExcelOutput{i,5}='5prime - 3prime';
            
        case 3
            ExcelOutput{i,3}=SeqLength(UniGene(i))-PosCell{UniRows(i),UniCols(i)};
            ExcelOutput{i,4}='AntiSense';
            ExcelOutput{i,5}='3prime - 5prime';
            
        case 4
            ExcelOutput{i,3}=PosCell{UniRows(i),UniCols(i)};
            ExcelOutput{i,4}='AntiSense';
            ExcelOutput{i,5}='5prime - 3prime';
    end
    temp=upper(MatchCell{UniRows(i),UniCols(i)});
    if iscell(temp)
        temp=temp{1};
    end
    
    temp(MisMatchInds(UniCols(i),:)&true)=lower(temp(MisMatchInds(UniCols(i),:)&true));

    ExcelOutput{i,6}=temp;
    ExcelOutput{i,7}=sum(MisMatchInds(UniCols(i),:));
    
end


if nargout>0
    varargout=cell(1,nargout);
    for i=1:nargout
        varargout{i}=ExcelOutput(:,i);
    end
end

if ~isempty(EXCEL_FILENAME)
    xlswrite(EXCEL_FILENAME,ExcelOutput);
end





