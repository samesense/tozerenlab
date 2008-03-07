function GenerateScatterPlots(PAT_STRUCT,varargin)
%   GenerateScatterPlots
%       Uses the data in a PAT_STRUCT to create a set of time-dependant
%       scatterplots based on the CD4 and Viral-Load.  This function also
%       allows robust control of the set-up and labeling of the
%       scatterplots.
%
%
%
%       GenerateScatterPlots(PAT_STRUCT)
%
%           Creates the following figures:
%       FIG1:
%           A set of histograms of the Viral loads at each normalized
%           time point.
%
%       FIG2:
%           A set of scatter plots of CD4 vs. Viral load at each normalized
%           time point sampled.
%
%       FIG3:
%           A set of scatter plots of the CD4 slope vs. Viral slope between
%           each normalized time point.
%
%       FIG4:
%           A set of histograms of the CD4 slope vs. Viral slope between
%           each normalized time point.
%
%
%
%   Optional Parameters:
%
%       FIG_ONLY        [true]|false        [1:4]
%                       This arguement can be used to control which (if any
%                       figures) to generate as described above.  If the
%                       arguement is TRUE then all figures will be created
%                       while FALSE creates no figures.  You can also
%                       enter a vector indicating which figures to create.
%
%       TIME_POINTS     Controls the normalied time points at which the
%                       scatterplots are displayed.  If this does not match
%                       the data already in the PAT_STRUCT then it will be
%                       re-calculated.
%
%       FORCE_NORM      Forces the algorithm to re-normalize even if data
%                       is present in the PAT_STRUCT.
%
%       NORM_PROPS      Arguements that will be passed to
%                       DetermineResponders if normalization is needed.
%                       See DetermineResponders for more information on
%                       valid properties.
%
%       FORCE_GROUPS    {COLOR_ARRAY, GROUP_MEMBERSHIP, GROUP_LABLES}
%                       Assuming there are N patients in M groups:
%                       COLOR_ARRAY =   Mx3 matrix of each color for each
%                                       group.
%                       GROUP_MEMBERSHIP = Nx1 vector of the membership to
%                                          each group.  Use 0 if there is 
%                                          no group membership.
%                       GROUP_LABELS =  Mx1 cell array of label names for
%                                       each group.
%                   
%




FIG_TOTAL=4;

FIG_PROVIDED=false;
NUM_PLOTS=6;
SUBPLOT_ROWS=2;
SUBPLOT_COLS=3;
PLOT_NEED_CORR=false;
PLOT_BEEN_CORR=false;
FIGURE_FLAG=true;
TIMEPOINT_FLAG=false;
GROUPS_PROVIDED=false;

NEED_NORMALIZATION=false;

NORM_PROPS={};

NUM_HIST_BINS=15;


FIGURES_TO_GENERATE=1:FIG_TOTAL;



if ~isempty(varargin)
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case {'fig_only','fig only','figure_only','figure only'}
                if islogical(varargin{i+1})
                    FIGURES_TO_GENERATE=1:FIG_TOTAL*varargin{i+1};
                elseif isvector(varargin{i+1})&&all(varargin{i+1}>0&&varargin{i+1}<=FIG_TOTAL)
                    FIGURES_TO_GENERATE=varargin{i+1};
                else
                    error('GenerateScatterPlots:BAD_FIG_ONLY_ARG','Arguement to FIG_ONLY must either be a logical scalar or a numerical vector [1 %d].',FIG_TOTAL);
                end

            case {'timepoint','time_points'}
                WANTED_TIME_POINTS=varargin{i+1}(:);
                TIMEPOINT_FLAG=true;
                PLOT_NEED_CORR=true;

                NUM_PLOTS=length(time_points);
                NORM_PROPS=[NORM_PROPS(:);'time_points',WANTED_TIME_POINTS];

            case {'normalization_props','norm_props'}
                if iscell(varargin{i+1})
                    NORM_PROPS=[NORM_PROPS(:); varargin{i+1}(:)];
                    NEED_NORMALIZATION=true;
                else
                    error('GenerateScatterPlots:BAD_NORM_PROPS_ARG','Arguement to NORM_PROPS must be a 1xN cell array.');
                end
            case {'force_norm'}
                if islogical(varargin{i+1})
                    NEED_NORMALIZATION=varargin{i+1};
                else
                    error('GenerateScatterPlots:BAD_FORCE_NORM_ARG','Arguement to FORCE_NORM must be a logical scalar.');
                end

            case {'force groups','force_groups'}
                if iscell(varargin{i+1})
                    GROUPS_PROVIDED=true;
                    temp_arg=varargin{i+1};
                    if iscell(temp_arg)&&size(temp_arg{1},2)==3
                        DESIRED_COLORS=temp_arg{1};
                        NUM_GROUPS=size(DESIRED_COLORS,1);
                    else
                        error('GenerateScatterPlots:BAD_FORCE_GROUPS_ARG','Arguement to FORCE_GROUPS must be cell-array with an mx3 color array and a Nx1 array showing group membership.');
                    end

                    GROUP_MEMBERSHIP=temp_arg{2}(:);
                    if ~any(size(GROUP_MEMBERSHIP)==[length(PAT_STRUCT) 1])||max(GROUP_MEMBERSHIP)>NUM_GROUPS
                        error('GenerateScatterPlots:BAD_FORCE_GROUPS_ARG','Arguement to FORCE_GROUPS must be cell-array with an mx3 color array and a Nx1 array showing group membership.');
                    end

                    if length(temp_arg)==3&&iscell(temp_arg{3})&&length(temp_arg{3})==NUM_GROUPS
                        GROUP_LABELS=temp_arg{3}(:);
                    else
                        GROUP_LABELS=arrayfun(@(x)int2str(x),1:NUM_GROUPS,'uniformoutput',false);
                    end

                    COLOR_ARRAY=zeros(length(GROUP_MEMBERSHIP),3);
                    for k=1:NUM_GROUPS
                        COLOR_ARRAY(GROUP_MEMBERSHIP==k,:)=ones(nnz(GROUP_MEMBERSHIP==k),1)*DESIRED_COLORS(k,:);
                    end

                else
                    error('GenerateScatterPlots:BAD_FORCE_GROUPS_ARG','Arguement to FORCE_GROUPS must be cell-array with an mx3 color array and a Nx1 array showing group membership.');
                end

            otherwise
                error('GenerateScatterPlots:BADPARAM','An unknown parameter was provided: %s',varargin{i})
        end
    end
end


[combined_norm_data combined_norm_index]=PatientStructHelper(PAT_STRUCT,...
    {'NORM_CLINICAL_DATA','explodeCells'});

if ~TIMEPOINT_FLAG
    WANTED_TIME_POINTS=combined_norm_data{1}(:,1);
end

%If the timepoints were provided then we have to check to make sure
%that they correspond with the PAT_STRUCT

all_correct=cellfun(...
    @(x)(size(x,1)==size(WANTED_TIME_POINTS,1)&&all(WANTED_TIME_POINTS==x(:,1))),...
    combined_norm_data);

if any(~all_correct)    %if they don't match then re-normalize
    NEED_NORMALIZATION=true;
else
    %if they match then we need to get the RESPONDER information out
    [RESPONDERS RESP_index]=PatientStructHelper(PAT_STRUCT,...
        {'IS_RESPONDER','leaveNumeric'});
end



if NEED_NORMALIZATION
    PAT_STRUCT=DetermineResponders(PAT_STRUCT,NORM_PROPS{:},'pat_struct_out',true);

    [combined_norm_data combined_norm_index ...
        RESPONDERS RESP_index]=PatientStructHelper(PAT_STRUCT,...
        {'NORM_CLINICAL_DATA','explodeCells'},...
        {'IS_RESPONDER','leaveNumeric'});

end


RNA_time_norm=cell2mat(cellfun(@(x)(x(:,2)'),combined_norm_data,'uniformoutput',false));
CD_time_norm=cell2mat(cellfun(@(x)(x(:,3)'),combined_norm_data,'uniformoutput',false));

RNA_good_rows=sum(isnan(RNA_time_norm),2)<=4;
CD_good_rows=sum(isnan(CD_time_norm),2)<=4;

valid_indexed=combined_norm_index(RNA_good_rows&CD_good_rows);

%%%%%%Collate all data so the GROUP_MEMBERSHIP, normalized data, and color
%%%%%%data are on the same rows.

if ~GROUPS_PROVIDED     %no groups were provided so create them from the responder data
    [junk good_NORM good_RESP]=intersect(valid_indexed,RESP_index);

    GROUP_MEMBERSHIP=RESPONDERS(good_RESP)+1;

    VALID_GROUPS=1:2;

    DESIRED_COLORS=[1 0 0;0 0 1];

    COLOR_ARRAY=zeros(length(GROUP_MEMBERSHIP),3);
    COLOR_ARRAY(GROUP_MEMBERSHIP==1,:)=ones(nnz(GROUP_MEMBERSHIP==1),1)*[1 0 0];
    COLOR_ARRAY(GROUP_MEMBERSHIP==2,:)=ones(nnz(GROUP_MEMBERSHIP==2),1)*[0 0 1];

    GROUP_LABELS={'Non-Responders','Responders'};

else
    %if groups were provided then we need to make sure that the indexes
    %match up with those that have clinical data

    [junk good_NORM good_PROVIDED]=intersect(valid_indexed,(1:length(PAT_STRUCT)).*(GROUP_MEMBERSHIP~=0)');

    GROUP_MEMBERSHIP=GROUP_MEMBERSHIP(good_PROVIDED);
    COLOR_ARRAY=COLOR_ARRAY(good_PROVIDED,:);
    
    %deal with a special case where an entire group is removed
    VALID_GROUPS=unique(GROUP_MEMBERSHIP)';

end
NUM_GROUPS=length(VALID_GROUPS);

RNA_NORM=RNA_time_norm(valid_indexed(good_NORM),:);
CD_NORM=CD_time_norm(valid_indexed(good_NORM),:);
NORM_TIME_POINTS=WANTED_TIME_POINTS;




RNA_NUM_MEMBERS=zeros(NUM_GROUPS,size(RNA_NORM,2));
CD_NUM_MEMBERS=zeros(NUM_GROUPS,size(RNA_NORM,2));
BOTH_NUM_MEMBERS=zeros(NUM_GROUPS,size(RNA_NORM,2));

for i=1:NUM_GROUPS
    RNA_NUM_MEMBERS(i,:)=sum(~isnan(RNA_NORM(GROUP_MEMBERSHIP==VALID_GROUPS(i),:)));
    CD_NUM_MEMBERS(i,:)=sum(~isnan(CD_NORM(GROUP_MEMBERSHIP==VALID_GROUPS(i),:)));
    BOTH_NUM_MEMBERS(i,:)=sum(~isnan(CD_NORM(GROUP_MEMBERSHIP==VALID_GROUPS(i),:))&~isnan(RNA_NORM(GROUP_MEMBERSHIP==VALID_GROUPS(i),:)));
end


% RNA_min=min(RNA_NORM(:));
% RNA_max=max(RNA_NORM(:));
% 
% CD_min=min(CD_NORM(:));
% CD_max=max(CD_NORM(:));


RNA_min=0;
RNA_max=14.5;

CD_min=0;
CD_max=2000;

edges=[{linspace(RNA_min,RNA_max,NUM_HIST_BINS)} {linspace(CD_min,CD_max,NUM_HIST_BINS)}];
x_edges=edges{1};
bin_width=(RNA_max-RNA_min)/20;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FIG1:
%           A set of histograms of the Viral loads at each normalized
%           time point.
if any(FIGURES_TO_GENERATE==1)
    figure
    N=zeros(NUM_GROUPS,length(x_edges));
    for i=1:NUM_PLOTS
        h=subplot(SUBPLOT_ROWS,SUBPLOT_COLS,i);
        hold on

        N_tot=hist(RNA_NORM(:,i),x_edges);
        N_tot_norm=N_tot.*(1/size(RNA_NORM,1));

        %calulate the historgam of each group membership
        N=cell2mat(arrayfun(@(x)(hist(RNA_NORM(GROUP_MEMBERSHIP==x,i),x_edges)),(VALID_GROUPS)','uniformoutput',false));
        %each row is a group and each column is an edge

        %Normalize each row by the number of samples in that group
        N_norm=((1./RNA_NUM_MEMBERS(:,i))*ones(1,length(x_edges))).*N;


        bar_series=bar(h,edges{1},N_norm',1);
        BarSeriesChanger(bar_series);

        %         line_series=plot(h,edges{1},N_norm');
        %         LineSeriesChanger(line_series);

        axis([RNA_min RNA_max 0 1])
        xlabel('Viral RNA count')
        title([num2str(NORM_TIME_POINTS(i)) ' Weeks'])
    end
    legend(GROUP_LABELS{VALID_GROUPS})
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FIG2:
%           A set of scatter plots of CD4 vs. Viral load at each normalized
%           time point sampled.
if any(FIGURES_TO_GENERATE==2)
    figure
    for i=1:NUM_PLOTS
        h=subplot(SUBPLOT_ROWS,SUBPLOT_COLS,i);
        hold on;

        for k=VALID_GROUPS
            scatter(h,RNA_NORM(GROUP_MEMBERSHIP==k,i),CD_NORM(GROUP_MEMBERSHIP==k,i),10*ones(nnz(GROUP_MEMBERSHIP==k),1),COLOR_ARRAY(GROUP_MEMBERSHIP==k,:),'filled')
        end

        title([num2str(NORM_TIME_POINTS(i)) ' Weeks'])
        axis([RNA_min RNA_max CD_min CD_max])
        xlabel('Viral RNA count')
        ylabel('CD4 cell count')
        legend(GROUP_LABELS{VALID_GROUPS})
    end
    
end

%%%%%%%%Calculate the slopes of the patient's normalized data
RNA_vels=[zeros(size(RNA_NORM,1),1) (RNA_NORM(:,2:end)-RNA_NORM(:,1:end-1))./(ones(size(RNA_NORM,1),1)*diff(NORM_TIME_POINTS'))];
CD_vels=[zeros(size(CD_NORM,1),1) (CD_NORM(:,2:end)-CD_NORM(:,1:end-1))./(ones(size(CD_NORM,1),1)*diff(NORM_TIME_POINTS'))];

RNA_vels_min=min(RNA_vels(:));
RNA_vels_max=max(RNA_vels(:));
CD_vels_min=min(CD_vels(:));
CD_vels_max=max(CD_vels(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FIG3:
%           A set of scatter plots of the CD4 slope vs. Viral slope between
%           each normalized time point.
if any(FIGURES_TO_GENERATE==3)
    figure
    for i=2:NUM_PLOTS
        h=subplot(SUBPLOT_ROWS,SUBPLOT_COLS,i);
        hold on;

        for k=VALID_GROUPS
            scatter(h,RNA_vels(GROUP_MEMBERSHIP==k,i),CD_vels(GROUP_MEMBERSHIP==k,i),10*ones(nnz(GROUP_MEMBERSHIP==k),1),COLOR_ARRAY(GROUP_MEMBERSHIP==k,:),'filled')
        end

        title([num2str(NORM_TIME_POINTS(i)) ' Weeks'])
        axis([RNA_vels_min RNA_vels_max CD_vels_min CD_vels_max])
        xlabel('Viral RNA slope')
        ylabel('CD4 cell slope')
    end
    legend(GROUP_LABELS{VALID_GROUPS})
end

RNA_vels_edges=linspace(RNA_vels_min,RNA_vels_max,20);

RNA_VELS_NUM_MEMBERS=zeros(NUM_GROUPS,size(RNA_NORM,2));
for i=VALID_GROUPS
    RNA_VELS_NUM_MEMBERS(i,:)=sum(~isnan(RNA_vels(GROUP_MEMBERSHIP==i,:)));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FIG4:
%           A set of histograms of the CD4 slope vs. Viral slope between
%           each normalized time point.
if any(FIGURES_TO_GENERATE==4)
    figure
    for i=2:NUM_PLOTS
        h=subplot(SUBPLOT_ROWS,SUBPLOT_COLS,i);


        %calulate the historgam of each group membership
        N=cell2mat(arrayfun(@(x)(hist(RNA_vels(GROUP_MEMBERSHIP==x,i),RNA_vels_edges)),(VALID_GROUPS)','uniformoutput',false));
        %each row is a group and each column is an edge

        %Normalize each row by the number of samples in that group
        N_norm=((1./RNA_VELS_NUM_MEMBERS(:,i))*ones(1,length(RNA_vels_edges))).*N;        
        
        bar_series=bar(h,RNA_vels_edges,N_norm',1);
        BarSeriesChanger(bar_series);
        
        
        axis([RNA_vels_min RNA_vels_max 0 1])
        xlabel('Viral RNA Velocity')
        title([num2str(NORM_TIME_POINTS(i)) ' Weeks'])
    end
    legend(GROUP_LABELS{VALID_GROUPS})
end


%%%%%%%%%%%%%%%%%HELPER FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function BarSeriesChanger(bar_handles)
        for p=VALID_GROUPS
            set(bar_handles(p),'FaceColor',DESIRED_COLORS(p,:));
        end
    end

    function LineSeriesChanger(line_series)
        for p=VALID_GROUPS
            set(line_series(p),'Color',DESIRED_COLORS(p,:),'linewidth',2.0);
        end
    end
end














