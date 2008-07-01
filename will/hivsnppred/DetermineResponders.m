function [varargout]=DetermineResponders(PAT_STRUCT, varargin)
%   DetermineResponders
%       Determines which of a set of patients have a statistically
%       improved prognosis.  This is determined by looking for a decreased
%       viral load and increased CD4 cell count over many time-points.
%
%   RESPONDERS=DetermineResponders(PAT_STRUCT)
%
%       RESPONDERS      An Nx1 boolean vector where true indicates GOOD
%                       PROGNOSIS
%
%   [RESPONDERS TIME_POINTS]=DetermineResponders(PAT_STRUCT)
%
%
%   [PAT_STRUCT ...]=DetermineResponders(PAT_STRUCT,'PAT_STRUCT_OUT',true,...)
%
%       PAT_STRUCT      A structure with the information determined from
%                       MakeSNPCalls and associated with a PATIENT.  This
%                       can be used by other programs to display the data.
%                       This program will add or update the following
%                       fields.
%
%                           NORM_CLINICAL_DATA
%                               [TIME(WEEKS) CD4_count RNA_count]
%
%                       Any subsequent output arguements are appended as
%                       desbribed above.%
%
%
%


%%%%%%%%%%%%%%%PARSE INPUT ARGUEMENTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIG_PROVIDED=false;
NUM_BINS=6;
TIMEPOINT_FLAG=false;
NEED_PAT_STRUCT=false;
NUM_TIME_POINTS_NEEDED=3;
SD_METHOD_FLAG = false;

if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case {'timepoints','time_points'}
                if isnumeric(varargin{i+1})&&isvector(varargin{i+1})
                    time_points=sort(varargin{i+1}(:));
                    TIMEPOINT_FLAG=true;
                    PLOT_NEED_CORR=true;

                    NUM_BINS=length(time_points);
                else
                    error('DetermineResponders:BADTIMEPOINTS','The arguement to TIMEPOINTS must be a numerical vector')
                end
            case {'numbins','num_bins'}
                if isnumeric(varargin{i+1})&&isscalar(varargin{i+1})&&varargin{i+1}>0
                    NUM_BINS=varargin{i+1};
                else
                    error('DetermineResponders:BADNUMBINS','The arguement to NUM_BINS must be a scalar')
                end
            case {'pat_struct_out'}
                if islogical(varargin{i+1})
                    NEED_PAT_STRUCT=varargin{i+1};
                else
                    error('DetermineResponders:NEED_PAT_STRUCT_ARG','Arguement to NEED_PAT_STRUCT must be logical');
                end
            case 'sdmethod'
                if islogical(varargin{i+1})
                    SD_METHOD_FLAG=varargin{i+1}&true;
                else
                    error('DetermineResponders:SD_METHOD','Arguement to SD_METHOD must be logical');
                end
            otherwise
                error('DetermineResponders:BADPARAM','An unknown parameter was provided: %s',varargin{i})
        end
    end

end


%%%%%%%%%%Get Data from the Patient Structure

[CD_vecs CD_ind_array ...
    RNA_vecs RNA_ind_array]=...
    PatientStructHelper(PAT_STRUCT,...
    {'CD4_vec','leaveCells'},...
    {'RNA_vec','leaveCells'});

CD_max=max(cell2mat(cellfun(@(x)(max(x(:,2),[],1)),CD_vecs,'uniformoutput',false)));
CD_min=min(cell2mat(cellfun(@(x)(min(x(:,2),[],1)),CD_vecs,'uniformoutput',false)));

RNA_max=max(cell2mat(cellfun(@(x)(max(x(:,2),[],1)),RNA_vecs,'uniformoutput',false)));
RNA_min=min(cell2mat(cellfun(@(x)(min(x(:,2),[],1)),RNA_vecs,'uniformoutput',false)));

if ~TIMEPOINT_FLAG
    %if timepoints are not provided determine our own from the histogram of
    %the time data: choose the most commonly reported times so there is
    %less interpolation needed

    RNA_CD_timepoints=cell2mat(cellfun(@(x)(x(:,1)),[RNA_vecs;CD_vecs],'uniformoutput',false));
    [N X]=hist(RNA_CD_timepoints(RNA_CD_timepoints>=0),500);
    [Y I]=sort(N,'descend');
    time_points=sort([0 round(X(I(2:NUM_BINS)))]);
end

%%%%Interpolate the values wherever needed and retrieve the RNA and CD4
%%%%counts
RNA_time_norm=cell2mat(cellfun(@(x)(InterpolateVals(x)),RNA_vecs,'uniformoutput',false));
CD_time_norm=cell2mat(cellfun(@(x)(InterpolateVals(x)),CD_vecs,'uniformoutput',false));


% spot=sum((isnan(RNA_time_norm(:,:))&isnan(CD_time_norm(:,:))),2)==0;
% edges=[{linspace(RNA_min,RNA_max,20)} {linspace(CD_min,CD_max,20)}];
% x_edges=edges{1};
% bin_width=(RNA_max-RNA_min)/20;
% 
% poss_resp=zeros(NUM_BINS,length(x_edges));
% resp_spots=false(length(PAT_STRUCT),NUM_BINS);
% 
% %%%%%Find the most common RNA bin
% for i=1:NUM_BINS
%     N=hist(RNA_time_norm(spot,i),x_edges);
%     [Y I]=sort(N,'descend');
%     poss_resp(i,:)=x_edges(I)';
% end
% 
% %%%%%%Determine which ones are less than the common bin
% for i=1:NUM_BINS
%     resp_spots(:,i)=(((poss_resp(end,1)+bin_width)>RNA_time_norm(:,i))&((poss_resp(end,1)-bin_width)<RNA_time_norm(:,i)))|...
%         (((poss_resp(end,2)+bin_width)>RNA_time_norm(:,i))&((poss_resp(end,2)-bin_width)<RNA_time_norm(:,i)));
% end
% 
% 
% %%%%%%Determine which ones are consitently less than the common bin
% late_responder=(sum(resp_spots(:,2:NUM_BINS),2)>=NUM_TIME_POINTS_NEEDED);
% RESPONDERS=late_responder&spot;





if ~SD_METHOD_FLAG
    RNA_time_diff=RNA_time_norm(:,1:end-1)-RNA_time_norm(:,2:end);
    RESPONDERS=nanmean(RNA_time_diff<0,2)>=0.4&mean(isnan(RNA_time_diff),2)<0.8;
else

    RESPONDERS = cellfun(@SD_Method,RNA_vecs);
    
    
end




%%%%%%%%Setup Output
varargout=SetupOutput(nargout);






%%%%%%%%%%%%%%%SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function FINAL=InterpolateVals(input)
        x_vals=unique(input(:,1));
        if size(x_vals,1)<3
            %not enough input values
            FINAL=ones(size(time_points))*NaN;
        else
            if size(x_vals,1)~=size(input,1)
                %Average repeated vals
                y_vals=zeros(size(x_vals));
                for j=1:length(x_vals)
                    y_vals(j)=mean(input(input(:,1)==x_vals(j),2));
                end
            else
                y_vals=input(:,2);
            end

            FINAL=interp1(x_vals,y_vals,time_points,'spline',NaN);
            FINAL(FINAL<0)=NaN(nnz(FINAL<0),1);
        end
    end

    function output_cell=SetupOutput(input)
        output_cell=cell(input,1);
        start=1;
        if NEED_PAT_STRUCT
            %%%%%Collate data output
            NORMED_DATA=arrayfun(@(x)([time_points; RNA_time_norm(x,:); CD_time_norm(x,:)]'),1:length(RNA_time_norm)','uniformoutput',false);
            
            %%%%%ADD PAT_STRUCT IF DESIRED
            output_cell{start}=PatientStructHelper(PAT_STRUCT,...
                {'NORM_CLINICAL_DATA',(1:length(PAT_STRUCT))',NORMED_DATA'},...
                {'IS_RESPONDER',(1:length(PAT_STRUCT))',RESPONDERS})

            start=2;
        end

        switch input-(start-1)
            case 1
                output_cell{start}=RESPONDERS;
            case 2

                output_cell{start}=RESPONDERS;
                output_cell{start+1}=time_points;
        end
    end

    function status = SD_Method(input)
        x_vals = unique(input(:,1));
        if size(x_vals,1)<3
            %not enough input values
            status=false;
        else
            if size(x_vals,1)~=size(input,1)
                %Average repeated vals
                y_vals=zeros(size(x_vals));
                for j=1:length(x_vals)
                    y_vals(j)=mean(input(input(:,1)==x_vals(j),2));
                end
            else
                y_vals=input(:,2);
            end

            FINAL=interp1(x_vals,y_vals,[0 8],'spline',NaN);
            
            status = diff(FINAL)<=-1.5;
        end
        
        
        
    end
end
