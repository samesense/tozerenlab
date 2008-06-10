function INDEXES = MatchAlignments(TOP_SEQ,BOTTOM_SEQS,PROVIDED_ALIGNMENTS)
%   MatchAlignments
%       Matches a cell-array of pairwise alignments to a set of alignment
%       sequences.  This allows the submission of a batch of alignments and
%       sequences and returning the indexes of the alignments which match
%       the sequences.
%
%   INDEX = MatchAlignments(TOP_SEQ,BOTTOM_SEQS,PROVIDED_ALIGNMENTS)
%
%   TOP_SEQ                 The sequence that should be on the first-row of
%                           all alignments.
%
%   BOTTOM_SEQS             A cell-array of sequences that should match the
%                           bottom-row of the alignment.
%
%   PROVIDED_ALIGNMENTS     A cell-array of pair-wise alignments to
%                           retrieve.
%
%
%   INDEX                   An index vector where:
%       PROVIDED_ALIGNMENTS(INDEX(i)) is the alignment for TOP_SEQ X
%                                                          BOTTOM_SEQS(i)
%
%                           and 0 where no such alignment is in
%                           PROVIDED_ALIGNMENTS
%
%
%
%

if iscell(TOP_SEQ)
    TOP_SEQ = TOP_SEQ{1};
end

TopSeqs = cellfun(@(x)(isequal(x(1,isletter(x(1,:))),TOP_SEQ(isletter(TOP_SEQ)))),PROVIDED_ALIGNMENTS);
if all(~TopSeqs)
    INDEXES = zeros(size(BOTTOM_SEQS));
    return
end

validSeqs = find(TopSeqs);


BottomSeqHashes = cellfun(@(x)(hash(x(isletter(x)),'MD5')),BOTTOM_SEQS,'uniformoutput',false);
ProvidedHashes = cellfun(@(x)(hash(x(end,isletter(x(end,:))),'MD5')),PROVIDED_ALIGNMENTS,'uniformoutput',false);

[junk INDEXES] = ismember(BottomSeqHashes,ProvidedHashes(TopSeqs));

INDEXES(INDEXES~=0) = validSeqs(INDEXES(INDEXES~=0));


