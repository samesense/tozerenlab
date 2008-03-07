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


warning('off','stats:kmeans:EmptyCluster')
N=size(MATCH_SPOTS,1);
ELM_annot=cell(1,length(ELM_STRUCT));
current_center=0;
for i=1:length(ELM_STRUCT)
    sizes=cellfun('length',MATCH_SPOTS(:,i));

    if nnz(sizes)>length(SEQS)*0.1  %ONLY consider ELMs that are in at least 10% of sequences

        max_length=max(sizes);

        found_inds=sizes~=0;
        ELM_locs=cell2mat(cellfun(@(x)([x NaN(1,max_length-length(x))]),MATCH_SPOTS(found_inds,i),'uniformoutput',false));

%                 figure
%                 hist(ELM_locs(:),1:3600)
%                 hold

        [ident centers]=kmeans(ELM_locs(~isnan(ELM_locs(:))),round(1.5*max_length),'emptyaction','drop','replicates',100,'start','cluster');

        %%%%%Remove and rename clusters accounting for empty clusters
        new_centers=centers(~isnan(centers));
        new_index=cumsum([current_center; ~isnan(centers)]);
        new_ident=new_index(ident);

        current_center=new_index(end);


        ELM_locs(~isnan(ELM_locs(:)))=new_ident;

        MATCH_SPOTS(found_inds,i)=num2cell(ELM_locs,2);
        MATCH_SPOTS(~found_inds,i)=num2cell(NaN(nnz(~found_inds),size(ELM_locs,2)),2);
        ELM_annot{i}=[new_centers'; i*ones(1,length(new_centers))];

%                 bar(centers(~isnan(centers)),nonzeros(histc(ident,1:length(centers))),'r')
    end
end

ELM_vec=[ELM_annot{:}];



for i=1:length(PAT_STRUCT)

    bool_vec=ismember(1:current_center,unique([MATCH_SPOTS{ALIGNMENT_INDS==i,:}]));

    PAT_STRUCT{i}.ELM_simple=ismember(1:length(ELM_STRUCT),unique(ELM_vec(2,bool_vec)));
    PAT_STRUCT{i}.ELM_vec=bool_vec;
    PAT_STRUCT{i}.ELM_annot=ELM_vec;
end

close(WAITBAR_HANDLE)

end