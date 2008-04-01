function varargout=DNAPromoterMatcher(SEARCH_SEQ,SEQ_DB_FILENAME,NUM_MISMATCH,EXCEL_FILENAME)
%   DNAPromoterMatcher
%       Searches through a database of sequences to find the SEARCH_SEQ.
%       This is optimized to find exact (or near exact) matches throughout
%       the promoter-eome.  If you are looking for the closest match of
%       SEARCH_SEQ in a small subset of promoters use ClosestDNAMatch.
%
%   [ GeneSymbol LLID POS STRAND ORIENT SEQ] =...
%               DNAPromoterMatcher(SEARCH_SEQ,SEQ_DB_FILENAME,NUM_MISMATCH)
%
%   SEARCH_SEQ          The querry sequence.
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

%search_seq ='AATACAATTAAAT';


WAITBAR_HANDLE=waitbar(0,'Loading Sequences');

%%%%%%%%%%%%%%%INPUT CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==3
    EXCEL_FILENAME=[];
end

try
    [header seqs]=fastaread(SEQ_DB_FILENAME);
    shorter_seqs=cellfun(@(x)(x(2000:end)),seqs,'uniformoutput',false);
catch
    warning('DNAPromoterMatcher:BAD_SEQ_DB','Cannot read Sequence Database.')
    close(WAITBAR_HANDLE);
    rethrow(lasterror)
end

%%%shorten database to 2000bp.



%%%%%%%%%%%%%%DONE INPUT CHECKING

M=length(SEARCH_SEQ);

%%%SINGLE MUT
if NUM_MISMATCH>0
    trans_mat=[eye(M);eye(M)*2;eye(M)*3;eye(M)*4];
else
    final_mat=[]; %#ok<NASGU>
end

if NUM_MISMATCH>1
    final_mat=[];

    for k=1:NUM_MISMATCH-1
        for i=1:M*4
            final_mat=[final_mat;repmat(trans_mat(i,:),[M*4 1])+[diag(~trans_mat(i,:),0);diag(~trans_mat(i,:),0)*2;diag(~trans_mat(i,:),0)*3;diag(~trans_mat(i,:),0)*4]]; %#ok<AGROW>
        end
        trans_mat=final_mat;
    end
else
    final_mat=trans_mat;
end

%%%add the exact match to the search
final_mat=[zeros(1,M);final_mat];


%%%%remove mis-matches don't change the search sequence
input_seq=nt2int(SEARCH_SEQ);
for i=1:M
    final_mat=final_mat(final_mat(:,i)~=input_seq(i),:);
end



match_mat=zeros(length(shorter_seqs),length(final_mat),4);
waitbar(0,WAITBAR_HANDLE,'Checking');
for i=1:length(final_mat)
    this_seq=input_seq;
    this_seq(final_mat(i,:)~=0)=nonzeros(final_mat(i,:));
    this_seq=int2nt(this_seq);
    match_mat(:,i,1)=FASTBoyer_Moore({this_seq},shorter_seqs');
    match_mat(:,i,2)=FASTBoyer_Moore({seqreverse(this_seq)},shorter_seqs');
    match_mat(:,i,3)=FASTBoyer_Moore({seqcomplement(this_seq)},shorter_seqs');
    match_mat(:,i,4)=FASTBoyer_Moore({seqrcomplement(this_seq)},shorter_seqs');
    waitbar(i/length(final_mat),WAITBAR_HANDLE)
end



J_nums=sum(final_mat~=0,2);
spots=find(match_mat);
[I J K]=ind2sub(size(match_mat),spots);

excel_output=cell(length(J),7);
temp2=cell(length(J),1);
for i=1:length(J)
    temp=textscan(header{I(i)},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
    excel_output(i,1)=temp{1};
    excel_output(i,2)=temp{2};
    excel_output{i,3}=2000-match_mat(I(i),J(i),K(i));
    switch (K(i))
        case 1
            excel_output{i,4}='Sense';
            excel_output{i,5}='3prime - 5prime';
        case 2
            excel_output{i,4}='Sense';
            excel_output{i,5}='5prime - 3prime';
        case 3
            excel_output{i,4}='Anti-Sense';
            excel_output{i,5}='3prime - 5prime';
        case 4
            excel_output{i,4}='Anti-Sense';
            excel_output{i,5}='5prime - 3prime';
    end
    this_seq=input_seq;
    this_seq(final_mat(J(i),:)~=0)=nonzeros(final_mat(J(i),:));
    this_seq=int2nt(this_seq);
    this_seq(final_mat(J(i),:)~=0)=lower(this_seq(final_mat(J(i),:)~=0));
    excel_output{i,6}=this_seq;
    excel_output{i,7}=J_nums(J(i));
    temp2{i}=[excel_output{i,:}];
end

[junk I J]=unique(temp2); %#ok<NASGU,ASGLU>

excel_output=excel_output(I,:);



if ~isempty(EXCEL_FILENAME)
    waitbar(1,WAITBAR_HANDLE,'Writing Excel File')
    for i=0:NUM_MISMATCH
        mismatch_num=cellfun(@(x)(x==i),excel_output(:,7));
        xlswrite(EXCEL_FILENAME,[{'Gene Symbol','EntrezID','BP Upstream','Strand','Orientation','Sequence','Mismatches'};excel_output(mismatch_num,:)],i+1)
    end
end

if nargout>0
    for i=1:nargout
        varargout{i}=excel_output(:,i);
    end
else
    varargout=[];
end

close(WAITBAR_HANDLE);







% exact_mat=zeros(length(shorter_seqs),4);
% this_seq=SEARCH_SEQ;
% exact_mat(:,1)=FASTBoyer_Moore({this_seq},shorter_seqs');
% exact_mat(:,2)=FASTBoyer_Moore({seqreverse(this_seq)},shorter_seqs');
% exact_mat(:,3)=FASTBoyer_Moore({seqcomplement(this_seq)},shorter_seqs');
% exact_mat(:,4)=FASTBoyer_Moore({seqrcomplement(this_seq)},shorter_seqs');
% 
% [I J]=find(exact_mat);
% 
% 
% excel_output=cell(length(I),5);
% for i=1:length(I)
%     temp=textscan(header{I(i)},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
%     excel_output(i,1)=temp{1};
%     excel_output(i,2)=temp{2};
%     excel_output{i,3}=2000-exact_mat(I(i),J(i));
%     switch (J(i))
%         case 1
%             excel_output{i,4}='Sense';
%             excel_output{i,5}='3prime - 5prime';
%         case 2
%             excel_output{i,4}='Sense';
%             excel_output{i,5}='5prime - 3prime';
%         case 3
%             excel_output{i,4}='Anti-Sense';
%             excel_output{i,5}='3prime - 5prime';
%         case 4
%             excel_output{i,4}='Anti-Sense';
%             excel_output{i,5}='5prime - 3prime';
%     end
% end
% 
% xlswrite('TFtargets.xls',[{'Gene Symbol','EntrezID','BP Upstream','Strand','Orientation'};excel_output],2)













