function varargout=DNAPromoterMatcher(SEARCH_SEQ,SEQ_DB,NUM_MISMATCH,varargin)
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
%
%   Optional Properties
%
%   EXCEL_OUTPUT        A filename to output the data in EXCEL format.
%
%   HTML_OUTPUT         A filename to output marked-up seqeunces.
%
%   COLUMN_WIDTH        A property to adjust the fixed-width of columns in
%                       the HTML_OUTPUT. DEFAULT=50;
%
%   MAX_PRO_LENGTH      The maximum length of allowed promoters.  All
%                       promoter regions are shortened to this length.
%                       Anything shorter than the given length is left
%                       alone. DEFAULT = 2000.
%
%
%
%   See also: ClosestDNAMatch.
%
%

% SEARCH_SEQ ='AA[AT]ACAA[AT]TAA[AT]';
% SEQ_DB_FILENAME='unique_IDS.fa';

%WAITBAR_HANDLE=waitbar(0,'Loading Sequences');

%%%%%%%%%%%%%%%INPUT CHECKING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GLOBAL_LENGTH=2000;
EXCEL_FLAG=false;
EXCEL_FILENAME=[];
HTML_FLAG=false;
HTML_FILENAME=[];
COLUMN_WIDTH=50;

if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case 'excel_output'
                if ischar(varargin{i+1})
                    EXCEL_FLAG=true;
                    EXCEL_FILENAME=varargin{i+1};
                else
                    error('DNAPromoterMatcher:BAD_EXCELFILENAME','The arguement to EXCEL_FILENAME must be a char-array')
                end

            case 'html_output'
                if ischar(varargin{i+1})
                    HTML_FLAG=true;
                    HTML_FILENAME=varargin{i+1};
                else
                    error('DNAPromoterMatcher:BAD_HTMLFILENAME','The arguement to HTML_FILENAME must be a char-array')
                end
                
            case 'column_width'
                if isnumeric(varargin{i+1})&&isscalar(varargin{i+1})&&varargin{i+1}>1
                    COLUMN_WIDTH=varargin{i+1};
                else
                    error('DNAPromoterMatcher:BAD_COLUMNWIDTH','The arguement to COLUMN_WIDTH must be a numeric scalar')
                end
                
            case 'max_pro_length'
                if isnumeric(varargin{i+1})&&isscalar(varargin{i+1})&&varargin{i+1}>1
                    GLOBAL_LENGTH=varargin{i+1};
                else
                    error('DNAPromoterMatcher:BAD_MAXPROLENGTH','The arguement to MAX_PRO_LENGTH must be a numeric scalar')
                end
                
            otherwise
                error('DNAPromoterMatcher:BAD_ARG','An unknown arguement was provided: %s', varargin{i})

        end
    end
end


%deal with loading FASTA-file
if ischar(SEQ_DB)
    try
        [headers seqs]=fastaread(SEQ_DB);
        LLIDS=cell(size(headers));
        GeneNames=cell(size(headers));
        for i=1:length(headers)
            temp=textscan(headers{i},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
            GeneNames(i)=temp{1};
            LLIDS(i)=temp{2};
        end

    catch
        warning('DNAPromoterMatcher:BAD_SEQ_DB','Cannot read Sequence Database.')
        %close(WAITBAR_HANDLE);
        rethrow(lasterror)
    end
elseif isstruct(SEQ_DB)
    [seqs{1:length(SEQ_DB)}]=deal(SEQ_DB.Sequence);

    %if data is easily available, then use it
    if isfield(SEQ_DB,'LLID')&&isfield(SEQ_DB,'GeneName')
        [LLIDS{1:length(SEQ_DB)}]=deal(SEQ_DB.LLID);
        [GeneNames{1:length(SEQ_DB)}]=deal(SEQ_DB.GeneName);
    else
        %otherwise take from header information
        try
            [headers{1:length(SEQ_DB)}]=deal(SEQ_DB.Headers);

            LLIDS=cell(size(headers));
            GeneNames=cell(size(headers));
            for i=1:length(headers)
                temp=textscan(headers{i},'%*s%*s%s%s%*s%*s%*s','delimiter','|');
                GeneNames(i)=temp{1};
                LLIDS(i)=temp{2};
            end
        catch
            warning('DNAPromoterMatcher:BAD_NAMEINFO','Cannot understand header information, continuing anyway.')
            LLIDS=num2cell(1:length(SEQ_DB));
            [GeneNames{1:length(SEQ_DB)}]=deal('UNKNOWN');
        end
    end
end

%%shorten database to 2000bp.
SeqLength=cellfun('length',seqs);
seqs=arrayfun(@(x,y)(x{1}(max(y-GLOBAL_LENGTH,1):end)),seqs,SeqLength,'uniformoutput',false);

SeqLength=cellfun('length',seqs);

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




%waitbar(0,WAITBAR_HANDLE,'Creating Sequence-Cells')
PromoterSeqs=[seqs'; ...
    cellfun(@(x)(seqreverse(x)),seqs','uniformoutput',false); ...
    cellfun(@(x)(seqcomplement(x)),seqs','uniformoutput',false); ...
    cellfun(@(x)(seqrcomplement(x)),seqs','uniformoutput',false)];


PromoterType=[(1:length(seqs))' ones(length(seqs),1); (1:length(seqs))' 2*ones(length(seqs),1); ...
    (1:length(seqs))' 3*ones(length(seqs),1); (1:length(seqs))' 4*ones(length(seqs),1)];

%waitbar(0,WAITBAR_HANDLE,'Searching Seqs')
MatchCell=cell(length(PromoterSeqs),length(MisMatchInds));
PosCell=cell(length(PromoterSeqs),length(MisMatchInds));
EndCell=cell(length(PromoterSeqs),length(MisMatchInds));
for i=1:size(MisMatchInds,1)
    %waitbar(i/size(MisMatchInds,1))
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

    [PosCell(:,i) MatchCell(:,i) EndCell(:,i)]=regexpi(PromoterSeqs,RegSearch,'start','match','end');
end

%%%%%Explodes cells in which there are multiple matches for each
%%%%%MisMatchInds
counter=1;
while counter<=size(PosCell,2)
    NeedExpand=cellfun('length',PosCell(:,counter));
    if any(NeedExpand>1)
        NewTempPos=cellfun(@(x)(x(end)),PosCell(NeedExpand>1,counter),'uniformoutput',false);
        ShortenedPos=cellfun(@(x)(x(1:end-1)),PosCell(NeedExpand>1,counter),'uniformoutput',false);

        NewTempEnd=cellfun(@(x)(x(end)),EndCell(NeedExpand>1,counter),'uniformoutput',false);
        ShortenedEnd=cellfun(@(x)(x(1:end-1)),EndCell(NeedExpand>1,counter),'uniformoutput',false);

        NewTempMatch=cellfun(@(x)(x(end)),MatchCell(NeedExpand>1,counter),'uniformoutput',false);
        ShortenedMatch=cellfun(@(x)(x(1:end-1)),MatchCell(NeedExpand>1,counter),'uniformoutput',false);


        PosCell(NeedExpand>1,counter)=ShortenedPos;
        MatchCell(NeedExpand>1,counter)=ShortenedMatch;
        EndCell(NeedExpand>1,counter)=ShortenedEnd;

        if counter~=size(PosCell,2)
            PosCell=[PosCell(:,1:counter) cell(size(PosCell,1),1) PosCell(:,counter+1:end)];
            PosCell(NeedExpand>1,counter+1)=NewTempPos;

            EndCell=[EndCell(:,1:counter) cell(size(EndCell,1),1) EndCell(:,counter+1:end)];
            EndCell(NeedExpand>1,counter+1)=NewTempEnd;

            MatchCell=[MatchCell(:,1:counter) cell(size(MatchCell,1),1) MatchCell(:,counter+1:end)];
            MatchCell(NeedExpand>1,counter+1)=NewTempMatch;

            MisMatchInds=[MisMatchInds(1:counter,:); MisMatchInds(counter,:); MisMatchInds(counter+1:end,:)];

        else
            PosCell=[PosCell(:,1:counter) cell(size(PosCell,1),1)];
            PosCell(NeedExpand>1,counter+1)=NewTempPos;

            EndCell=[EndCell(:,1:counter) cell(size(EndCell,1),1)];
            EndCell(NeedExpand>1,counter+1)=NewTempEnd;

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

%waitbar(0,WAITBAR_HANDLE,'Generating Output')
UniGene=PromoterType(I);

UniCols=J;
UniRows=I;

%%%%%%%[GeneSymbol LLID POS STRAND ORIENT SEQ NUM_MISS]
ExcelOutput=cell(length(UniGene),7);
for i=1:length(UniGene)
    ExcelOutput{i,1}=['''' GeneNames{UniGene(i)}];
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

if HTML_FLAG
    fid=fopen(HTML_FILENAME,'wt');
    fprintf(fid,'<html> \n <body> \n <FONT face="Courier"> \n');
    for i=1:length(seqs)
        tempMarkup=false(8,SeqLength(i));
        for j=1:4
            StartPos=[PosCell{i+(j-1)*length(seqs),:}];
            EndPos=[EndCell{i+(j-1)*length(seqs),:}];

            if ~isempty(StartPos)
                tempMarkup(1+2*(j-1),StartPos)=true;
                tempMarkup(2+2*(j-1),EndPos)=true;
            end
        end
        %only print those with markup
        if any(tempMarkup(:))
            tempMarkup(3:4,:)=fliplr(tempMarkup([4 3],:));
            tempMarkup(7:8,:)=fliplr(tempMarkup([8 7],:));

            fprintf(fid,'<br> %s \t %s <br>\n <br>\n',GeneNames{i},LLIDS{i});
            if any(tempMarkup(1,:))
                fprintf(fid,'<FONT COLOR="#FF0000"> 3-5 Sense Pos: %d </FONT> <br> \n',find(tempMarkup(1,:)));
            end
            if any(tempMarkup(3,:))
                fprintf(fid,'<FONT COLOR="#00FF00"> 5-3 Sense Pos: %d </FONT> <br> \n',find(tempMarkup(3,:)));
            end
            if any(tempMarkup(5,:))
                fprintf(fid,'<FONT COLOR="#0000FF"> 3-5 Anti-Sense Pos: %d </FONT> <br> \n',find(tempMarkup(5,:)));
            end
            if any(tempMarkup(7,:))
                fprintf(fid,'<FONT COLOR="#00FFFF"> 5-3 Anti-Sense Pos: %d </FONT> <br> \n',find(tempMarkup(7,:)));
            end

            counter=1;
            k=1;
            while k<=size(tempMarkup,2)
                %print an entire line if possible
                if k+COLUMN_WIDTH<size(tempMarkup,2)&&~any(sum(tempMarkup(:,k:k+COLUMN_WIDTH-1)))&&counter==1;
                    fprintf(fid,'%s <br> \n',seqs{i}(k:k+COLUMN_WIDTH-1));
                    k=k+COLUMN_WIDTH;
                else    %otherwise go char-by-char
                    spot=find(tempMarkup(:,k));
                    if ~isempty(spot)
                        switch spot(1)
                            case 1
                                %red
                                fprintf(fid,'%s<FONT COLOR="#FF0000">',seqs{i}(k));
                            case 2
                                fprintf(fid,'%s</FONT>',seqs{i}(k));
                            case 3
                                %green
                                fprintf(fid,'%s<FONT COLOR="#00FF00">',seqs{i}(k));
                            case 4
                                fprintf(fid,'%s</FONT>',seqs{i}(k));
                            case 5
                                %blue
                                fprintf(fid,'%s<FONT COLOR="#0000FF">',seqs{i}(k));
                            case 6
                                fprintf(fid,'%s</FONT>',seqs{i}(k));
                            case 7
                                %cyan
                                fprintf(fid,'%s<FONT COLOR="#00FFFF">',seqs{i}(k));
                            case 8
                                fprintf(fid,'%s</FONT>',seqs{i}(k));
                        end
                    else
                        fprintf(fid,'%s',seqs{i}(k));
                    end
                    k=k+1;
                    counter=counter+1;
                    if counter>COLUMN_WIDTH
                        fprintf(fid,'<br> \n');
                        counter=1;
                    end
                end
            end
            fprintf(fid,'<br> \n');
        end

    end
    fprintf(fid,'</FONT> \n </body> \n </html> \n');
    fclose(fid);
end


if nargout>0
    varargout=cell(1,nargout);
    for i=1:nargout
        varargout{i}=ExcelOutput(:,i);
    end
end

if EXCEL_FLAG
    xlswrite(EXCEL_FILENAME,ExcelOutput);
end





