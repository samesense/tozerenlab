function ClosestDNAMatch(SEARCH_SEQ,SEQ_DB_FILENAME,DESIRED_LLIDS,EXCEL_FILENAME)
%   ClosestDNAMatch
%       Searches through a subset of the database of sequences to find the
%       closest match to SEARCH_SEQ.  This is optimized to find the closest
%       match in a small number of sequences.  If you are looking for exact
%       (or near exact) matches use DNAPromoterMatcher.
%
%   ClosestDNAMatch(SEARCH_SEQ,SEQ_DB_FILENAME,LLIDS,EXCEL_FILENAME)
%
%   SEARCH_SEQ          The querry sequence.
%
%   SEQ_DB_FILENAME     The filename of the sequences.  These must be in
%                       the format provided by PAINT's Upstreamer Function.
%                       http://www.dbi.tju.edu/dbi/tools/paint/
%
%   EXCEL_FILENAME      A filename to output the data in EXCEL format.
%
%
%
%
%   See also: DNAPromoterMatcher.
%
%

WAITBAR_HANDLE=waitbar(0,'Loading Sequences');

PROMOTER_DB_LENGTH=4000;
WANTED_PROMOTER_LENGTH=2000;


%search_seq ='AATACAATTAAAT';

%%%%%%%%%%%%%%%INPUT CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==3
    EXCEL_FILENAME=[];
end

try
    [header seqs]=fastaread(SEQ_DB_FILENAME);
    shorter_seqs=cellfun(@(x)(x(PROMOTER_DB_LENGTH-WANTED_PROMOTER_LENGTH:end)),seqs,'uniformoutput',false);
    LLIDS=zeros(length(header),1);
    for i=1:length(header)
        LLIDS(i)=cell2mat(textscan(header{i},'%*s%*s%*s%d%*s%*s%*s','delimiter','|'));
    end
catch
    warning('DNAPromoterMatcher:BAD_SEQ_DB','Cannot read Sequence Database.')
    close(WAITBAR_HANDLE)
    rethrow(lasterror)
end

%%%shorten database to 2000bp.


M=length(SEARCH_SEQ);



available_IDS=ismember(DESIRED_LLIDS,LLIDS);
if all(~available_IDS)
    error('ClosestDNAMatch:ALL_MISSING_IDS','All of the LLIDS [%d] were not in the Database.  Cannot Continue!',DESIRED_LLIDS(~available_IDS))
elseif any(~available_IDS)
    warning('ClosestDNAMatch:SOME_MISSING_IDS','Some of the LLIDS [%d] were not in the Database',DESIRED_LLIDS(~available_IDS))
end

%%%%%%%%%%%%%%DONE INPUT CHECKING

spots=find(ismember(LLIDS,DESIRED_LLIDS));
this_search=[SEARCH_SEQ; seqreverse(SEARCH_SEQ); seqcomplement(SEARCH_SEQ); seqrcomplement(SEARCH_SEQ)];

output=cell(1,length(spots));



for k=1:length(spots)
    waitbar(k/length(spots),WAITBAR_HANDLE,['Checking LLID: ' num2str(LLIDS(spots(k)))])
    desired_MATCH=zeros(4,WANTED_PROMOTER_LENGTH);
    seq=repmat(cell2mat(shorter_seqs(spots(k))),[4 1]);

    for i=1:(WANTED_PROMOTER_LENGTH-length(SEARCH_SEQ)-1)
        desired_MATCH(:,i)=sum(seq(:,i:i+M-1)==this_search,2);
    end
    num_matches=histc(desired_MATCH(:),1:M);
    best_match_score=find(num_matches>0,1,'last');

    [I J]=find(desired_MATCH==best_match_score);
    temp_output=cell(length(J),6);
    for i=1:length(J)
        this_seq=seq(1,J(i):J(i)+M-1);
        switch(I(i))
            case 1
                temp_output{i,4}='Sense';
                temp_output{i,5}='3prime - 5prime';
            case 2
                temp_output{i,4}='Sense';
                temp_output{i,5}='5prime - 3prime';
                this_seq=seqreverse(this_seq);
            case 3
                temp_output{i,4}='Anti-Sense';
                temp_output{i,5}='3prime - 5prime';
                this_seq=seqcomplement(this_seq);
            case 4
                temp_output{i,4}='Anti-Sense';
                temp_output{i,5}='5prime - 3prime';
                this_seq=seqrcomplement(this_seq);
        end

        miss_mask=this_seq~=SEARCH_SEQ;
        temp_output{i,1}=LLIDS(spots(k));
        temp_output{i,2}=nnz(miss_mask);
        temp_output{i,3}=J(i);
        this_seq(miss_mask)=lower(this_seq(miss_mask));
        temp_output{i,6}=this_seq;
    end
    output{k}=temp_output;

end

if ~isempty(EXCEL_FILENAME)
    waitbar(1,WAITBAR_HANDLE,'Writting Excel File')
    for i=1:length(spots)
        xlswrite(EXCEL_FILENAME,[{'EntrezID','Mismatches','BP Upstream','Strand','Orientation','Sequence'};output{i}],i)
    end
end

end