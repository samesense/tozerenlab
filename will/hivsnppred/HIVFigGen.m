function HIVFigGen(HIV_STRUCT,LOCS,HEIGHTS,NAMES,varargin)
%   HIVFigGen
%       Generates an annotated HIV figure based on an HIV_STRUCT.  It will
%       create a bar of a specified height at a specified location on the
%       x-axis.  Below the x-axis will be an annotated HIV genome.
%
%
%
%   HIVFigGen(HIV_STRUCT,LOCS,HEIGHTS,NAMES)
%
%       HIV_STRUCT      An HIV structure created by HIVrefLoader
%
%       LOCS            The genomic positions to annotate
%
%       HEIGHTS         The height to make the marker.  If left empty then
%                       the all markers are 1.0 units heigh.
%
%       NAMES           A name to annotated the marker.  If left empty then
%                       the marker is annotated with the numeric order it
%                       was provided in.
%
%
%
%   Optional Properties
%
%       USE_AA          true|[false]    If set to TRUE then the LOCS will
%                       be interpreted as a position in the concatenated AA
%                       reference.
%
%       PLOT_LIMITS     This will force the figure to only plot within the
%                       provided region.
%
%       ANNOT_CELL      An annotation cell-array in which the first column
%                       is the BeginPos and EndPos of each gene and the
%                       second column is the GeneName.
%
%

IS_AA_FLAG = false;
MIN_VAL=[];
MAX_VAL=[];
ANNOT_CELL = [];

if ~isempty(varargin)
    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'use_aa'
                if islogical(varagin{i+1})
                    IS_AA_FLAG = varagin{i+1};
                else
                    error('HIVFigGen:BAD_USEAA','Arguement to USE_AA must be a logical');
                end
                
            case 'plot_limits'
                if isdouble(varargin{i+1})&&numel(varargin{i+1})==2
                    MIN_VAL = varargin{i+1}(1);
                    MAX_VAL = varargin{i+1}(2);
                else
                    error('HIVFigGen:BAD_PLOTLIMITS','Arguement to PLOT_LIMITS must be a 2 element double array');
                end
                
            case 'annot_cell'
                if iscell(varargin{i+1})&&size(varargin{i+1},2)==2
                    try
                        temp_vals=cell2mat(varargin{i+1}(:,1));
                        if ~isdouble(temp_vals)
                            error('HIVFigGen:BAD_ANNOTCELL','Bad arguement to ANNOT_CELL');
                        end
                    catch
                        error('HIVFigGen:BAD_ANNOTCELL','Bad arguement to ANNOT_CELL');
                    end
                    ANNOT_CELL = varargin{i+1};
                else
                    error('HIVFigGen:BAD_ANNOTCELL','Bad arguement to ANNOT_CELL');
                end
            
            otherwise
                error('HIVFigGen:BAD_ARG','An unknown input arguements was provided: %s',varargin{i});
        end
    end
end

if isempty(ANNOT_CELL)
    
    if ~IS_AA_FLAG
        ANNOT_CELL = cell(size(HIV_STRUCT.GenePos,1),2);
        ANNOT_CELL(:,1) = num2cell(HIV_STRUCT.GenePos,2);

    else
        ANNOT_CELL = cell(size(HIV_STRUCT.GenePos,1),2);
        temp = cellfun('length',HIV_STRUCT.AAseqs);
        spots = cumsum([1; temp]);

        ANNOT_CELL(:,1)= num2cell([spots(1:end-1) spots(2:end)],2);
    end
    ANNOT_CELL(:,2) = HIV_STRUCT.GeneNames;
end

if isempty(MIN_VAL)
    MIN_VAL = 1;
    MAX_VAL = ANNOT_CELL{end,1}(end);
end



THIS_HANDLE = axes;
hold on
axis([MIN_VAL MAX_VAL 0 1.5*max(HEIGHTS)])

bar(THIS_HANDLE,LOCS,HEIGHTS,0.1)
AnnotFig(MIN_VAL:MAX_VAL,ANNOT_CELL)

[vals, inds]=sort(LOCS);

for i = 1:length(LOCS)
    [epos1 epos2]=dsxy2figxy(gca,LOCS(inds(i)),HEIGHTS(inds(i)));

    thisAnnot = annotation('textarrow');
    
    set(thisAnnot,'string',NAMES{inds(i)},'X',[epos1 epos1],'Y',[epos2+.05 epos2],'textrotation',270,'HorizontalAlignment','center','VerticalAlignment','top')
end

hold off




function varargout = dsxy2figxy(varargin)
% dsxy2figxy -- Transform point or position from axis to figure coords
% Transforms [axx axy] or [xypos] from axes hAx (data) coords into coords
% wrt GCF for placing annotation objects that use figure coords into data
% space. The annotation objects this can be used for are
%    arrow, doublearrow, textarrow
%    ellipses (coordinates must be transformed to [x, y, width, height])
% Note that line, text, and rectangle anno objects already are placed
% on a plot using axes coordinates and must be located within an axes.
% Usage: Compute a position and apply to an annotation, e.g.,
%   [axx axy] = ginput(2);
%   [figx figy] = getaxannopos(gca, axx, axy);
%   har = annotation('textarrow',figx,figy);
%   set(har,'String',['(' num2str(axx(2)) ',' num2str(axy(2)) ')'])

%% Obtain arguments (only limited argument checking is performed).
% Determine if axes handle is specified
if length(varargin{1})== 1 && ishandle(varargin{1}) && ...
  strcmp(get(varargin{1},'type'),'axes') 
 hAx = varargin{1};
 varargin = varargin(2:end);
else
 hAx = gca;
end;
% Parse either a position vector or two 2-D point tuples
if length(varargin)==1 % Must be a 4-element POS vector
 pos = varargin{1};
else
 [x,y] = deal(varargin{:});  % Two tuples (start & end points)
end
%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Need normaized units to do the xform
axpos = get(hAx,'Position');
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));
%% Transform data from figure space to data space
if exist('x','var')     % Transform a and return pair of points
 varargout{1} = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
 varargout{2} = (y-axlim(3))*axpos(4)/axheight + axpos(2);
else                    % Transform and return a position rectangle
 pos(1) = (pos(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
 pos(2) = (pos(2)-axlim(3))/axheight*axpos(4) + axpos(2);
 pos(3) = pos(3)*axpos(3)/axwidth;
 pos(4) = pos(4)*axpos(4)/axheight;
 varargout{1} = pos;
end
%% Restore axes units
set(hAx,'Units',axun)
