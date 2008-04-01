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
%       ITER_DISPAY     [true|FALSE] If set to TRUE then the optimization
%                       technique will display the progress at each
%                       iterateration.  When set to FALSE then only the
%                       final window display is shown.
%
%       NO_DISPLAY      [true|FALSE] If set to TRUE then all displays are
%                       removed.  Good for batch processing and nested
%                       functions.
%
%       MaxIter         The Maximum Number of iterations to perform on each
%                       optimization round. DEFAULT=300
%
%       MaxTime         The Maximim Time (in seconds) to spend on each 
%                       optimization. DEFAULT=inf;
%
%
%
%
%
%

ITER_DISP_FLAG=false;
NO_DISPLAY_FLAG=false;
MaxIter=300;
MaxTime=inf;

MATCH_SPOTS=[];
OPT_METHOD_FLAG=true;


if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case {'matched_spots','match_spots'}
                if iscell(varargin{i+1})
                    MATCH_SPOTS=varargin{i+1};
                    if size(MATCH_SPOTS,1)~=size(SEQS,1)
                        error('ELMPositional:BAD_MATCHEDSPOTS_ROWS','MATCHED_SPOTS have the same number of rows as SEQS.')
                    elseif size(MATCH_SPOTS,2)~=length(ELM_STRUCT)
                        error('ELMPositional:BAD_MATCHEDSPOTS_COLS','MATCHED_SPOTS have the same number of columns as ELM_STRUCT.')
                    end
                else
                    error('ELMPositional:BAD_MATCHEDSPOTS','Arguement to MATCHED_SPOTS must be a cell-array.')
                end
            case 'opt_method'
                if islogical(varargin{i+1})
                    OPT_METHOD_FLAG=varargin{i+1};
                else
                    error('ELMPositional:BAD_OPTMETHOD','Arguement to OPT_METHOD must be a logical.')
                end
            case 'iter_display'
                if islogical(varargin{i+1})
                    ITER_DISP_FLAG=varargin{i+1};
                else
                    error('ELMPositional:BAD_ITERDISPLAY','Arguement to ITER_DISPLAY must be a logical.')
                end
            case 'no_display'
                if islogical(varargin{i+1})
                    NO_DISPLAY_FLAG=varargin{i+1};
                else
                    error('ELMPositional:BAD_NODISPLAY','Arguement to NO_DISPLAY must be a logical.')
                end
            case 'maxiter'
            otherwise
                error('ELMPositional:BAD_ARG','An unknown arguement was provided: %s',varargin{i})
        end
    end
end


%WAIT_HANDLE=waitbar(0,'Processing');

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

if ~OPT_METHOD_FLAG
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

        end
    end

    ELM_vec=[ELM_annot{:}];
    close(WAIT_HANDLE)


    POSITIONAL_CALL=false(size(SEQUENCES,1),size(ELM_vec,2));
    for i=1:size(POSITIONAL_CALL,1)
        POSITIONAL_CALL(i,:)=ismember(1:current_center,[MATCH_SPOTS{i,:}]);
    end

else
    if NO_DISPLAY_FLAG
        options=optimset('diffminchange',1,'diffmaxchange',10,'MaxIter',MaxIter);
    else
        options=optimset('OutputFcn', @EvalDisplay,'diffminchange',1,'diffmaxchange',10,'MaxIter',MaxIter);
    end
    
    ELM_Annot_Cell=cell(length(ELM_STRUCT),1);
    ELM_PosCall_Cell=cell(length(ELM_STRUCT),1);
    
    for i=1:length(ELM_STRUCT)
        CurrentFval=zeros(1,MaxIter);
        sizes=cellfun('length',MATCH_SPOTS(:,i));
        %waitbar(i/length(ELM_STRUCT),WAIT_HANDLE)

        LogicalData=false(length(SEQUENCES),max(sizes));

        for p=1:size(LogicalData,1)
            if ~isempty(MATCH_SPOTS{p,i})
                LogicalData(p,MATCH_SPOTS{p,i})=true;
            end
        end

        subplot(2,2,3), bar(1:size(LogicalData,2),sum(LogicalData,1))
        FvalLog=zeros(10,1);
        FposLog=unique(ceil(linspace(.5*max(sizes),2*max(sizes),10)));
        WindowLog=cell(length(FposLog),1);
        
        if ~ITER_DISP_FLAG&&~NO_DISPLAY_FLAG
            WAIT_HANDLE=waitbar(0,ELM_STRUCT(i).Name);
        end
        

        for testSlice=FposLog
            
            NumSlices=ceil(testSlice);
            x0=linspace(1,size(LogicalData,2),NumSlices+2);
            x0=x0(2:end-1);
            StartTime=clock;
            [WindowLog{nnz(FvalLog)+1} fval]=fminsearchbnd(@EvalWindows,x0,ones(size(x0)),size(LogicalData,2)*ones(size(x0)),options);
            FvalLog(nnz(FvalLog)+1)=fval;
            if ~NO_DISPLAY_FLAG
                subplot(2,2,2), plot(FposLog(1:nnz(FvalLog)),FvalLog(1:nnz(FvalLog)))
            end

        end

        [MinVal MinSpot]=min(FvalLog);
        tempWindows=WindowLog{MinSpot};
        
        ELM_Annot_Cell{i}=[repmat(i,[1 length(WindowLog{MinSpot})+1]);1 tempWindows; tempWindows size(LogicalData,2)];
        
        tempPosCall=cellfun(@CheckWindows,MATCH_SPOTS(:,i),'uniformoutput',false);
        
        ELM_PosCall_Cell{i}=cell2mat(tempPosCall);
        
        if ~ITER_DISP_FLAG&&~NO_DISPLAY_FLAG
            close(WAIT_HANDLE);
        end
    end
    
    POSITIONAL_CALL=cat(2,ELM_PosCall_Cell{:});
    ELM_vec=cat(2,ELM_Annot_Cell{:});   
    

end


    function ObjVal=EvalWindows(INPUT)
        %count the number of times that more than 2 matches occur in a window
        ObjVal=0;
        INPUT=sort([1 ceil(INPUT) size(LogicalData,2)]);

        for k=1:length(INPUT)-1
            ObjVal=ObjVal+sum(sum(LogicalData(:,INPUT(k):INPUT(k+1)-1),2)>1);
        end

    end

    function stop=EvalDisplay(x,OptimStruct,state)
        stop=etime(clock,StartTime)>MaxTime;
        
        if ITER_DISP_FLAG
            switch state
                case 'iter'
                    CurrentFval(OptimStruct.iteration+1)=OptimStruct.fval;

                    subplot(2,2,4), plot(1:OptimStruct.iteration+1,CurrentFval(1:OptimStruct.iteration+1))

                    display_mat=double(LogicalData);
                    display_mat(LogicalData)=30;
                    display_mat(:,ceil(x))=200;
                    display_mat(:,ceil(x)+1)=200;
                    display_mat(:,ceil(x)-1)=200;

                    subplot(2,2,1), image(display_mat)
                    
                    subplot(2,2,3), bar(1:size(LogicalData,2),sum(LogicalData,1))
                    for lineind=1:length(x)
                        line('xdata',[x(lineind) x(lineind)],[0 100],'color','r');
                    end

                    drawnow


                case 'interupt'

                case 'init'

                case 'done'

            end
        else
            switch state
                case 'iter'
                    CurrentFval(OptimStruct.iteration+1)=OptimStruct.fval;
                    waitbar(OptimStruct.iteration/300,WAIT_HANDLE)
                case 'interupt'
                case 'init'
                    waitbar(0,WAIT_HANDLE,[ELM_STRUCT(i).Name 'Iter ' int2str(nnz(FvalLog)+1) '/' int2str(length(FposLog))])
                case 'done'

                    subplot(2,2,4), plot(1:OptimStruct.iteration+1,CurrentFval(1:OptimStruct.iteration+1))

                    display_mat=double(LogicalData);
                    display_mat(LogicalData)=30;
                    display_mat(:,ceil(x))=200;
                    display_mat(:,ceil(x)+1)=200;
                    display_mat(:,ceil(x)-1)=200;

                    subplot(2,2,1), image(display_mat)
                    
                    subplot(2,2,3), bar(1:size(LogicalData,2),sum(LogicalData,1))
                    for lineind=1:length(x)
                        line([x(lineind) x(lineind)],[0 max(sum(LogicalData,1))],'color','r');
                    end
                    
                    
                    drawnow

            end

        end
    end

    function OutCalls=CheckWindows(TempMatchSpot)
        if ~isempty(TempMatchSpot)
        tempRight=repmat([1 tempWindows],[length(TempMatchSpot) 1]);
        tempLeft=repmat([tempWindows size(LogicalData,2)],[length(TempMatchSpot) 1]);
        tempMatch=repmat(TempMatchSpot',[1 size(tempRight,2)]);
        
        OutCalls=sum(tempMatch>tempRight&tempMatch<tempLeft,1)>0;
        else
            OutCalls=false(1,size(tempWindows,2)+1);
        end
        
    end


end


