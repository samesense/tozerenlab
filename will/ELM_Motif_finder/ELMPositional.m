function [POSITIONAL_CALL ELM_vec ELM_annot]=ELMPositional(SEQS,MAPPING_CELL,ELM_STRUCT,varargin)
%   ELMPositional
%       Finds positional based ELM presence/absence calls by using kmeans
%       clustering on the location of ELM calls.
%
%
%   [POSITIONAL_CALL ANNOT_VEC]=ELMPositional(INPUT_SEQS,MAPPING_CELL,ELM_STRUCT)
%
%       SEQS        An M x 1 set of search sequences.  These can be 
%                   char-arrays, cells or structures with a Sequence-field.
%       MAPPING_CELL    A cell-array of the mapping of each position in
%                       SEQS to a reference sequence.  If this is left
%                       empty then it is assumed that there are no gaps
%                       relative to a reference sequence.
%
%       ELM_STRUCT  An N x 1 structure of ELM regular expressions as 
%                   created by ELMDownloadParser.
%
%       POSITIONAL_CALL     A cell-array of the presence/absence calls 
%                           for each feature in each sequence.
%
%       ANNOT_VEC   A 2xN vector where the first row represents the postion
%                   in the sequence and the second row is the index of the
%                   ELM_STRUCT.
%
%
%   [CALLS ANNOT_VEC CENTERS]=ELMPositional(...)
%
%       CENTERS     An Nx1 cell-array of the centroids found by the kmeans
%                   algorithm for each ELM_STRUCT.
%
%   Optional Parameters
%
%       MATCHED_SPOTS   If you have already run ELMFinder then you can
%                       provided the input to avoid re-running the
%                       time-consuming process.
%                       
%
%
%
%
%
%
%

MATCH_SPOTS=[];

if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case 'matched_spots'
                MATCH_SPOTS=varargin{i+1};
        end
    end
end


WAIT_HANDLE=waitbar(0,'Processing');

if isstruct(SEQS)&&isfield(SEQS,'Sequence')
    [SEQUENCES{1:length(SEQS)}]=deal(SEQS.Sequence);
elseif iscell(SEQS)
    LOCS=~cellfun('isempty',SEQS);
    SEQUENCES=cell(size(SEQS));
    SEQUENCES(LOCS)=SEQS(LOCS); %#ok<FNDSB>
end

if isempty(MATCH_SPOTS)
    MATCH_SPOTS=ELMFinder(SEQUENCES,ELM_STRUCT,'nested',true);
end

if ~isempty(MAPPING_CELL)
    for i=1:length(ELM_STRUCT)
        MATCH_SPOTS(:,i)=cellfun(@(x,y)(y(x)),MATCH_SPOTS(:,i),MAPPING_CELL,'uniformoutput',false);
    end
end

warning('off','stats:kmeans:EmptyCluster')
N=size(MATCH_SPOTS,1);
ELM_annot=cell(1,length(ELM_STRUCT));
current_center=0;
for (i=1:length(ELM_STRUCT))
    sizes=cellfun('length',MATCH_SPOTS(:,i));
    
    waitbar(i/length(ELM_STRUCT),WAIT_HANDLE)
    
    if nnz(sizes)>length(SEQS)*0.1  %ONLY consider ELMs that are in at least 10% of sequences

        max_length=max(sizes);

        found_inds=sizes~=0;
        ELM_locs=cell2mat(cellfun(@(x)([x NaN(1,max_length-length(x))]),MATCH_SPOTS(found_inds,i),'uniformoutput',false));

%                 figure
%                 hist(ELM_locs(:),1:600)
                %hold

        [ident centers]=kmeans(ELM_locs(~isnan(ELM_locs(:))),round(1.5*max_length),'emptyaction','drop','replicates',50,'start','cluster');

        %%%%%Remove and rename clusters accounting for empty clusters
        new_centers=centers(~isnan(centers));
        
        new_index=cumsum([current_center; ~isnan(centers)]);
        new_ident=new_index(ident);

        current_center=new_index(end);


        ELM_locs(~isnan(ELM_locs(:)))=new_ident;
        

        MATCH_SPOTS(found_inds,i)=num2cell(ELM_locs,2);
        if ~all(found_inds) %avoid a special case
            MATCH_SPOTS(~found_inds,i)=num2cell(NaN(nnz(~found_inds),size(ELM_locs,2)),2);
        end
        
        ELM_annot{i}=[new_centers'; i*ones(1,length(new_centers))];

%                 bar(centers(~isnan(centers)),nonzeros(histc(ident,1:length(centers))),'r')
    end
end

ELM_vec=[ELM_annot{:}];
close(WAIT_HANDLE)


POSITIONAL_CALL=false(size(SEQUENCES,1),size(ELM_vec,2));
for i=1:size(POSITIONAL_CALL,1)
    POSITIONAL_CALL(i,:)=ismember(1:current_center,[MATCH_SPOTS{i,:}]);
end





