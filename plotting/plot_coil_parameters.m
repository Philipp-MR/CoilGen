function  plot_coil_parameters(coil_layouts,plot_title)


% Plot result parameters in one combined plot
valid_layouts=find(arrayfun(@(x) ~isempty(coil_layouts(x).out),1:numel(coil_layouts))); 
pot_levels=arrayfun(@(x) coil_layouts(x).out.num_levels,valid_layouts);
if coil_layouts(1).out.input_data.skip_postprocessing
coil_lengths=arrayfun(@(x) sum([coil_layouts(x).out.coil_parts(:).combined_loop_length]),valid_layouts);
coil_inductances=arrayfun(@(x) 0,valid_layouts).*10^3;
else
coil_lengths=arrayfun(@(x) sum([coil_layouts(x).out.coil_parts(:).coil_length]),valid_layouts);
coil_inductances=arrayfun(@(x) coil_layouts(x).out.coil_parts(1).coil_inductance,valid_layouts).*10^3;
end
coil_currents=arrayfun(@(x) coil_layouts(x).out.potential_step,valid_layouts);
linecolors={'r' [0 .5 0] 'b'};
line_styles={'-' '--' ':'};
marker_style={'o' 'd' '+'};
h3i=plotNy(pot_levels,coil_inductances,1,...
    pot_levels,coil_currents,2,...
    pot_levels,coil_lengths,3,...
    'Linewidth',1.5,...
    'YAxisLabels',{'[mH]' '[A]' '[m]'},...
    'XAxisLabel','n',...
    'Ylim', [min(coil_inductances).*0.9 max(coil_inductances).*1.1; min(coil_currents).*0.9 max(coil_currents).*1.1; min(coil_lengths).*0.9 max(coil_lengths).*1.1],...
    'TitleStr',plot_title,...
    'LineColor',linecolors,...
    'LineStyle',line_styles,...
    'MarkerStyle',marker_style,...
    'FontSize',15,...
    'Fontname','TimesNewRoman',...
    'Grid','on',...
    'LegendString',{'Inductance [mH]' 'Current [A]' 'Coil length [m]'});
%Change the axis colors to match the requested line colors
for ii=1:length(h3i.ax)
    set(h3i.ax(ii),'ycolor',linecolors{ii});
	set(h3i.ax(ii),'Fontweight','bold');
end







function varargout=plotNy(varargin)
% [handles]=plotNy(x1,y1,ax1,x2,y2,ax2,...,[options]) -or-
% [handles]=plotNy({x1,x2,...},{y1,y2,...},[ax1,ax2,...],[options]) -or-
% [handles]=plotNy(x_common,{y1,y2,...},[ax1,ax2,...],[options])
%
% This function plots the set on plots on the N axes specified. The
% following options can be specified *after* all of the plot information.
% Either leave blank [] or do not specify command to leave default
% If an output variable is given, the plot handles are returned. Note that
% the order in which these properties pairs are given does not matter, and
% that properties and responses are not case sensitive.
%
% Options:
%   String arguments:
%   --------------------------------------------------------------------------------
%    'YAxisLabels'	- A cell array of strings for the labels of each y-axis
%          syntax: {Ax1_label,Ax2_label,...} (default 'Y-Ax 1...N')
%    'XAxisLabel'   - Label string for the x-axis
%          syntax: 'Xlabelstr' (defaults to blank '')
%    'TitleStr'     - Label string for the figure's title
%          syntax: 'TitleStr' (defaults to blank '')
%    'LegendString' - A cell array of strings for the labels of each plot in legend
%          syntax: {P1name,P2name,...}(default - Plot N [Ax M])
%   --------------------------------------------------------------------------------
%
%   Graphic arguments:
%   --------------------------------------------------------------------------------
%    'LineWidth'    - Width of each plot line. If one argument is given it is applied
%                     to all plots.
%          syntax: 0.5 (default), or specify line widths [1,.1,3,...]. See 'help plot'
%    'LineColor'    - Color of each line in plot. Can be color strings or RGB triplets.
%                     If one argument is given it is applied to all plots.
%          syntax: {'col1',[r2 g2 b2],'col3',...}. See 'help plot'
%    'LineStyle'    - Line style of each plot line. If one argument is given it
%                     is applied to all plots.
%          syntax: '-' (default), or {'none','-','--',...}. See 'help plot'
%    'MarkerStyle'  - Style of marker to put at each data point in plot. If one
%                     argument is given it is applied to all plots.
%          syntax: 'none' (default), or specify {'.','+',..}. See 'help plot'
%    'MarkerSize'   - Size of marker to put at each data point in plot
%          syntax: 6 (default), or specify marker size. See 'help plot'
%    'Grid'         - Decision on whether or not to include grid. See 'help grid'
%          syntax: 'off' (default), 'minor', 'on'
%    'LegendLoc'    - Decision on where to place the legend. See 'help legend'
%          syntax: 'best' (default), 'none', 'Legloc' (follows 'legend' convention')
%    'XLim'         - Bounds for common X-axis of plot (must be [1x2])
%          syntax: [min max]
%    'Ylim'         - Bounds for each Y-axis of plot (must be [N_axis x 2])
%          syntax: [miny1 maxy1; miny2 maxy2; ...]
%   --------------------------------------------------------------------------------
%
%   Formatting arguments:
%   --------------------------------------------------------------------------------
%    'UpperMarginW' - Margin between axis and upper edge of figure in pixels
%          syntax: default 'auto'/0 fits ax_labels, or specify upper margin width in pixels
%    'LowerMarginW' - Margin between axis and lower edge of figure in pixels
%          syntax: defaults to fit x label, or specify lower margin width in pixels
%    'RightMarginW' - Margin between axis and right edge of figure in pixels
%          syntax: defaults to fit axis, or specify right margin width in pixels
%    'LeftMarginW'  - Margin between axis and left edge of figure in pixels
%          syntax: default 'auto'/0 fits tic labels, or specify left margin width in pixels
%    'AxisWidth'    - Width of each axis in pixels
%          syntax: default to fit tic labels, or specify axis width in pixels
%   --------------------------------------------------------------------------------
%
%   Other arguments:
%   --------------------------------------------------------------------------------
%    'Parent'       - Specify which parent figure/uipanel will host this plot
%          syntax: 0 (new plot,default) , figure or uipanel handle (will delete all other children)
%    'ColorOrd'     - The color order for plots. If n_plot>n_colors, loops back to top with next line style
%          syntax: default is MATLAB axis default, or specify a Nx3 matrix of color triplets. This is
%                  overridden if 'LineColor' vector is provided
%    'LineOrd'      - After using all colors, loops around to color(1,:) with new line style
%          syntax: default is {'-' '--' '-.' ':'}, or specify a cell array of line styles. This is
%                  overridden if 'LineStyle' argument is provided
%   --------------------------------------------------------------------------------
%Check if asking for too many output arguments
if nargout>1
    error('%sToo many output arguments given.\nSee ''help plotNy'' for more information.','')
end
%--------First Check Input arguments, convert to required format--------%
if nargin<3
    error('Not Enough Input arguments given.')
else
    i=1;
    if iscell(varargin{2}) %First two inputs are cell arrays of x and y
        %If we are given vector for x and cell array for y it means all share same x points
        if isnumeric(varargin{1})
            for i=1:numel(varargin{2})
                if ~all(size(varargin{1})==size(varargin{2}{i}))
                    error('%sWith single input vector x and cell array y, all y must be same size as x.\nSee ''help plotNy'' for more information.','')
                end
            end
            varargin{1}=repmat(varargin(1),size(varargin{2}));
        elseif ~iscell(varargin{1})
            error('%sUnrecognized input type for first input argument.\nSee ''help plotNy'' for more information.','')
        end
        
        if ~all(size(varargin{1})==size(varargin{2}))
            error('%sInput cell arrays must be of same size.\nSee ''help plotNy'' for more information.','')
        elseif ~all(size(varargin{2})==size(varargin{3})) || ~isnumeric(varargin{3}) || any(round(varargin{3})~=varargin{3})
            error('%sAxis assignment input array must an integer array of same size as X,Y cell arrays.\nSee ''help plotNy'' for more information.','')
        else
            X=varargin{1};
            Y=varargin{2};
            ax_in=varargin{3};
            nplots=numel(ax_in);
            start_inputs=4;
        end
    elseif iscell(varargin{1}) && ~iscell(varargin{2}) %Can't have only 1 of first two args be cell
        error('%s\nInput arguments must be in triplets or cell array form.\nSee ''help plotNy'' for more information.','')
    else %assume that we have been given triplets, need to loop through triplets
        %Guess at preallocation for cell arrays
        X{1,floor(nargin/3)}='';
        Y=X;
        ax_in=zeros(size(X));
        
        still_trips=true;
        nplots=0;
        
        while still_trips
            if i>nargin || ischar(varargin{i})
                if nplots==0
                    error('%s\nThe first arguments must be plot information.\nSee ''help plotNy'' for more information.','')
                else
                    break;
                end
            elseif ~all(size(varargin{i})==size(varargin{i+1}))
                error('%s\nInput vector x_i and y_i must be of same size.\nSee ''help plotNy'' for more information.','')
            elseif max(size(varargin{i+2}))>1 || ~isnumeric(varargin{i+2}) || round(varargin{i+2})~=varargin{i+2}
                error('%s\nAxis assignment input ax_i must an integer.\nSee ''help plotNy'' for more information.','')
            else
                nplots=nplots+1;
                X{nplots}=varargin{i};
                Y{nplots}=varargin{i+1};
                ax_in(nplots)=varargin{i+2};
                i=i+3;
            end
        end
        X(nplots+1:end)=[];
        Y(nplots+1:end)=[];
        ax_in(nplots+1:end)=[];
        start_inputs=i;
        
    end
end
[axinds,~,retind]=unique(ax_in(:));
nax=length(axinds);
axinds=1:nax;
ax_in=reshape(axinds(retind),size(ax_in));
% %Get monitor resolution for reference in resize
% set(0,'units','pixels')
% monitor_res = get(0,'screensize');
%Set default plotting options
options={...
    'YAxisLabel'	, strcat(repmat({'Axis '},[1 nax]),strsplit(num2str(1:nax),''));...
    'XAxisLabel'	, '';...
    'TitleStr'  	, '';...
    'Parent'        , 0;...
    'LegendString'  , 0;...
    'YRanges'       , '';...
    'XRanges'       , '';...
    'ColorOrd'      , '';...
    'LineOrd'       , 0;...
    };
%Impose Defaults Settings
%         l  r u  d
%marg_w= [40 28 46 46]; %Number of pixels for margins in each direction above
marg_w=  [0 28 0 46]; %Number of pixels for margins in each direction above
ax_w=45;    %Number of pixels for width of each axis
gridstr='grid off';
legloc='best';
linew=repmat(0.5,1,nplots); %line width of plots
msize=repmat(6,1,nplots);
mstyle=repmat({'none'},1,nplots);
fontname=nan;
fontsize=nan;
still_inputs=i<nargin;
i=start_inputs;
while still_inputs
    if ischar(varargin{i})
        if strcmpi(varargin{i},'XAxisLabel')
            if ischar(varargin{i+1})
                options{2,2}=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''XAxisLabels'' must be a string for label of x-axis.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'YAxisLabels')
            if nax==1 && ischar(varargin{i+1})
                options{1,2}=varargin(i+1);
            elseif iscellstr(varargin{i+1}) && numel(varargin{i+1})==nax
                options{1,2}=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''YAxisLabels'' must be a cell array of strings for labels of each y-axis.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'Markerstyle')
            if iscellstr(varargin{i+1}) && numel(varargin{i+1})==nplots
                mstyle=varargin{i+1};
            elseif ischar(varargin{i+1})
                mstyle=repmat(varargin(i+1),1,nplots);
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''Markerstyle'' must be a string for marker style, a cell array of strings.\nSee ''help plotNy'' for more information.','')
			end
		elseif strcmpi(varargin{i},'FontName')
            if ischar(varargin{i+1})
                fontname=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''FontName'' must be a string for figure''s font.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'Linestyle')
            if iscellstr(varargin{i+1}) && numel(varargin{i+1})==nplots
                lstyle=varargin{i+1};
            elseif ischar(varargin{i+1})
                lstyle=repmat(varargin(i+1),1,nplots);
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''Linestyle'' must be a string for line style, a cell array of strings.\nSee ''help plotNy'' for more information.','')
			end
		elseif strcmpi(varargin{i},'FontSize')
            if isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1 && varargin{i+1}>0
                fontsize=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''FontSize'' must be a single positive value for figure''s font size.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'LegendString')
            if iscellstr(varargin{i+1}) && numel(varargin{i+1})==nplots
                options{5,2}=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''LegendString'' must be a cell array of strings for labels of each plot.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'Titlestr')
            if ischar(varargin{i+1})
                options{3,2}=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''TitleStr'' must be a string for plot title.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'Grid')
            if any(strcmpi({'off' 'minor' 'on'},varargin{i+1}))
                %options{3,2}=varargin{i+1};
                gridstr=['grid ' varargin{i+1}];
            elseif ~isempty(varargin{i+1})
                error('%sAllowable grid arguments are ''off'', ''minor'', or ''on''.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'LegendLoc')
            
            possible={'north', 'south', 'east', 'west', 'northeast', 'northwest', 'southeast', 'southwest',...
                'northoutside', 'southoutside', 'eastoutside', 'westoutside', 'northeastoutside', 'northwestoutside', 'southeastoutside',...
                'southwestoutside',	'best',	'bestoutside',	'none'};
            
            if any(strcmpi(possible,varargin{i+1}))
                %options{3,2}=varargin{i+1};
                legloc=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sAllowable lengend locations arguments are available in ''legend'' documentation.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'Parent')
            if ~ishandle(varargin{i+1})
                error('%s\nArgument following ''Parent'' must be a either a figure or uipanel handle.\nSee ''help plotNy'' for more information.','')
            elseif any(strcmpi({'uipanel' 'figure'},get(varargin{i+1},'type')))
                options{4,2}=varargin{i+1};
            elseif ~isempty(varargin{i+1}) && varargin{i+1}~=0
                error('%sHandle type %s not supported as parent for this function, or cannot host this plot.\nSee ''help plotNy'' for more information.','',get(varargin{i+1},'type'))
            end
        elseif strcmpi(varargin{i},'YLim')
            if isnumeric(varargin{i+1}) && all(size(varargin{i+1})==[nax 2])
            	options{6,2}=varargin{i+1}';
            elseif ~isempty(varargin{i+1}) && ~strcmpi(varargin{i+1},'auto')
                error('%sArgument following ''YLim'' must be a 2xN_axes matrix, where first row is lower bounds, next row is upper bounds.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'XLim')
            if isnumeric(varargin{i+1}) && all(size(varargin{i+1}(:))==[2 1])
                options{7,2}=varargin{i+1};
            elseif ~isempty(varargin{i+1}) && ~strcmpi(varargin{i+1},'auto')
                error('%sArgument following ''XLim'' must be a vector where first argument is lower bound, second is upper bound.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'UpperMarginW')
            if isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1 && varargin{i+1}>=0
                marg_w(3)=varargin{i+1};
            elseif ~isempty(varargin{i+1}) && ~strcmpi(varargin{i+1},'auto')
                error('%sArgument following ''UpperMarginW'' must be a single positive value or ''auto''/0 for upper margin.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'LowerMarginW')
            if isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1 && varargin{i+1}>0
                marg_w(4)=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''LowerMarginW'' must be a single positive value for lower margin.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'RightMarginW')
            if isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1 && varargin{i+1}>0
                marg_w(2)=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''RightMarginW'' must be a single positive value for right margin.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'LeftMarginW')
            if isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1 && varargin{i+1}>=0
                marg_w(1)=varargin{i+1};
            elseif ~isempty(varargin{i+1}) && ~strcmpi(varargin{i+1},'auto')
                error('%sArgument following ''LeftMarginW'' must be a single positive value or ''auto''/0 for left margin.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'AxisWidth')
            if isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1 && varargin{i+1}>0
                ax_w=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''AxisWidth'' must be a single positive value for axis margins.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'LineWidth')
            if isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1 && varargin{i+1}>0
                linew=repmat(varargin{i+1},1,nplots);
            elseif isnumeric(varargin{i+1}) && numel(varargin{i+1})==nplots && min(varargin{i+1})>0
                linew=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''LineWidth'' must be a single positive value for line width.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'MarkerSize')
            if isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1 && varargin{i+1}>0
                msize=repmat(varargin{i+1},1,nplots);
            elseif isnumeric(varargin{i+1}) && numel(varargin{i+1})==nplots && min(varargin{i+1})>0
                msize=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''MarkerSize'' must be a single positive value for marker size.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'ColorOrd')
            if isnumeric(varargin{i+1}) && size(varargin{i+1},2)==3 && min(varargin{i+1}(:))>=0 && max(varargin{i+1}(:))<=1
                options{8,2}=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''ColorOrd'' must be a [n_color x 3] matrix of RGB triplets for plotting color (0.0 <= value <= 1.0).\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'LineOrd')
            if iscell(varargin{i+1})
                options{9,2}=varargin{i+1};
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''LineOrd'' must be a cell array of line styles for plotting.\nSee ''help plotNy'' for more information.','')
            end
        elseif strcmpi(varargin{i},'LineColor')
            if iscell(varargin{i+1}) && numel(varargin{i+1})==nplots
                lcolor=varargin{i+1};
            elseif ischar(varargin{i+1}) || (isnumeric(varargin{i+1}) && max(size(varargin{i+1}))==1)
                lcolor=repmat(varargin(i+1),1,nplots);
            elseif ~isempty(varargin{i+1})
                error('%sArgument following ''LineColor'' must be a cell array of line colors, ore one to repeat, for plotting.\nSee ''help plotNy'' for more information.','')
            end
        else
            error('%s is not a valid plot option.\nSee ''help plotNy'' for more information.',varargin{i})
        end
        i=i+2;
        if i>nargin
            break;
        end
    else
        error('%sOptions should be specified as strings (%ith input argument expected to be string).\nSee ''help plotNy'' for more information.','',i)
    end
end
%Do we need to make a figure, or is one specified?
if isnumeric(options{4,2})
    fh=figure('Units','pixels','Visible','off');
    fighost=true;
    
    %Get monitor resolution for reference in resize
    set(0,'units','pixels')
    monitor_res = get(0,'screensize');
else
    fh=options{4,2};
    set(fh,'units','pixels')
    delete(get(fh,'children'))
    if strcmpi('figure',get(fh,'type'))
        fighost=true;
        %Get monitor resolution for reference in resize
        set(0,'units','pixels')
        monitor_res = get(0,'screensize');
    elseif strcmpi('uipanel',get(fh,'type'))
        fighost=false;
        parhand=get(fh,'parent');
        maxnest=100;
        for i=1:maxnest
            if strcmpi('figure',get(parhand,'type'))
                break
            else
                if i==maxnest
                    error('Parent provided is too nested, cannot find figure handle to set data tip callback.')
                else
                    parhand=get(parhand,'parent');
                end
            end
        end
        
        %Get monitor resolution for reference in resize
        set(parhand,'units','pixels')
        monitor_res = get(parhand,'position');
    else
        error('plotNy''s Parent must be either a Figure or UIPanel. Parent ''%s'' not supported.',get(fh,'type'))
    end
end
%------------------------------------------------------------------------%
%Preallocate axis handle structure
h_ax=zeros(1,nax);
h_y=h_ax;
fpos=get(fh,'position');
def_pos=[marg_w(1)+(nax-1)*ax_w marg_w(4) fpos(3)-marg_w(2)-(marg_w(1)+(nax-1)*ax_w) fpos(4)-marg_w(3)-marg_w(4)];
def_pos(3:4)=max([1 1],def_pos(3:4));
for i=1:nax
    
    %Create an axis and offset from main one by i axis widths
    h_ax(i)=axes('Parent',fh,...
        'Unit','Pixels',...
        'Position',def_pos-[ax_w 0 0 0].*(i-1),...
        'ActivePositionProperty','position',...
        'YTickLabelRotation',45);
    
    %Set its units to normalized after doing pixel width math
    set(gca,'Units','normalized')
    if ~isnan(fontname)
        set(h_ax(i),'fontname',fontname);
    end
    if ~isnan(fontsize)
        set(h_ax(i),'fontsize',fontsize);
    end
    
    
    if i>1 %This is a dummy axis object, only need y-axis
        set(get(h_ax(i),'XAxis'),'Visible','off')
        set(h_ax(i),'color','none',...
            'Box','off');
        
        %Ensure that this axis cannot be set to active axis with click
        set(h_ax(i),'HitTest','off')
        
        %Set temporary y-axis lable in place
        h_y(i)=ylabel(h_ax(i),['Y' num2str(i)],...
            'units','normalized',...
            'rotation',45,...
            'HorizontalAlignment','left',...
            'position',[0 titpos(2) 0]);
		
		if ~isnan(fontname)
			set(h_y(i),'fontname',fontname);
		end
		if ~isnan(fontsize)
			set(h_y(i),'fontsize',fontsize);
		end
        
    else %This is the master axis where we will plot
        %Make title, take position as reference for y-axis labels
        hold(h_ax(1),'on')
        h_t=title(h_ax(1),'Title','units','normalized','interpreter', 'none');
        titpos=get(h_t,'Position');
        
        %Offset higher from title in order to miss possible exponents on y-axis
        titpos(2)=titpos(2)*1.01;
        
        %Now can create x axis label object
        h_x=xlabel(h_ax(1),'X-Axis','units','normalized');
		if ~isnan(fontname)
			set(h_x,'fontname',fontname);
		end
		if ~isnan(fontsize)
			set(h_x,'fontsize',fontsize);
		end
        
        %Set temporary y-axis lable in place
        h_y(i)=ylabel(h_ax(1),'Y1',...
            'units','normalized',...
            'rotation',45,...
            'HorizontalAlignment','left',...
            'position',[0 titpos(2) 0]);
		
		if ~isnan(fontname)
			set(h_y(i),'fontname',fontname);
		end
		if ~isnan(fontsize)
			set(h_y(i),'fontsize',fontsize);
		end
        
        %Set the grid setting specified by user
        eval(gridstr);
        
        %Capture zoom and pan event objects, set callbacks later once
        %complete handles list is created
        hzoom_ax1=zoom;
        hpan_ax1=pan;
    end
    
    %We can also override so the rotation of any axis will simply be
    %rejected and view is reset
    htemp=rotate3d;
    htemp.ActionPostCallback=@reject_rotate;
end
%set the current axis
if fighost
    set(fh,'CurrentAxes',h_ax(1))
else
    set(parhand,'CurrentAxes',h_ax(1))
end
%Now impose the correct title, xlabel, and ylabel strings
set(h_x,'string',options{2,2})
set(h_t,'string',options{3,2})
max_labh=-inf;
for i=1:nax
    set(h_y(i),'string',options{1,2}{i})
    
    %%******This is probably the extent in axis FOR not figure FOR!
    set(h_y(i),'Units','pixels')
    tempext=get(h_y(i),'Extent');
    set(h_y(i),'Units','normalized')
    max_labh=max(max_labh,tempext(4));
end
if marg_w(3)==0
    for i=1:nax
        %h_ax(i).Position(3)=h_ax(i).Position(3)-max_labh;
        set(h_ax(i),'Units','pixels')
        ptemp=get(h_ax(i),'Position');
        ptemp(4)=ptemp(4)-max_labh-(titpos(2)-1)*def_pos(4);
        set(h_ax(i),'Position',ptemp);
        set(h_ax(i),'Units','normalized')
    end
    marg_w(3)=max_labh;%-(titpos(2)-1)*def_pos(4);
else
    %Set to 1 so that resize event does not change margin for title offset
    titpos(2)=1;
end
fhtemp=figure('visible','off');
for i=nax:-1:1
    sphand(i)=subplot(1,nax,i);hold on
end
for i=nplots:-1:1
    subplot(1,nax,ax_in(i));
    
    %Plot each pair on their appropriate axis, check to make sure it is successsful
    try
        plot(X{i},Y{i},...
            'LineWidth',linew(i),...
            'markersize',msize(i),...
            'marker',mstyle{i},...
            'linestyle','-');
    catch err
        delete(fhtemp);delete(fh);
        if strcmp(err.identifier,'MATLAB:samelen')
            error('Each X-Y pair must be the same length (Plot %i).\nSee ''help plotNy'' for more information.',i)
        else
            rethrow(err)
        end
    end
end
for i=nax:-1:1
    %Get the default Ylims that matlab would impose
    %benifit from their robust logic
    ranges(:,i)=get(sphand(i),'Ylim');
    hold(h_ax(i),'on');
end
delete(fhtemp);
if isnumeric(options{6,2})
    for i=1:nax
%         if ~any(isnan(options{6,2}(:,i)))
%             ranges(:,i)=options{6,2}(:,i);
%         end
		original=ranges(:,i);
        if ~isnan(options{6,2}(1,i))
            ranges(1,i)=options{6,2}(1,i);
        end
        if ~isnan(options{6,2}(2,i))
            ranges(2,i)=options{6,2}(2,i);
		end
		
		if ranges(2,i)==ranges(1,i)
			ranges(:,i)=original;
		end
    end
end
%Now check how wide the left margin needs to be if left margin is set to auto
if marg_w(1)==0
    if max(abs(ranges(:,end)))>=1e4
        newlw=25;
    elseif max(abs(ranges(:,end)))>=1e3
        newlw=40;
    elseif max(abs(ranges(:,end)))>=1e2
        newlw=35;
    elseif max(abs(ranges(:,end)))>=10
        newlw=30;
    else
        newlw=25;
    end
    
    %Add slightly more room to all fir auto
    newlw=newlw+5;
    
    %offset for the negative sign
    if min(ranges(:,end))<0
        newlw=newlw+5;
    end
    
    %Update the axis left margin depending on how many digits shown in last axis
    for i=1:nax
        %h_ax(i).Position(3)=h_ax(i).Position(3)-max_labh;
        set(h_ax(i),'Units','pixels')
        ptemp=get(h_ax(i),'Position');
        ptemp(3)=ptemp(3)-newlw;
        ptemp(1)=ptemp(1)+newlw;
        set(h_ax(i),'Position',ptemp);
        set(h_ax(i),'Units','normalized')
    end
    marg_w(1)=newlw;
end
%Get Color order/specify line styles for plots
if isnumeric(options{8,2})
    %User has specified color order
    color_ord=options{8,2};
else
    %Use MATLAB default
    color_ord=get(h_ax(1),'ColorOrder');
end
if iscell(options{9,2})
    %User has specified line order
    line_ord=options{9,2};
else
    %Use function default
    line_ord={'-' '--' '-.' ':'};
end
%Do we know the line styles and colors of the plots?
makelstyle=false;
makelcolor=false;
if ~exist('lstyle','var')
    makelstyle=true;
end
if ~exist('lcolor','var')
    makelcolor=true;
end
if and(makelstyle,makelcolor)
    %If neither line style nor color is specified, loop through unique
    %combinations of color and style
    ccnt=1;lcnt=1;
    lstyle{1,nplots}='';
    lcolor{1,nplots}='';
    for i=1:nplots
        lcolor{i}=color_ord(ccnt,:);
        lstyle{i}=line_ord{lcnt};
        
        %Keep running update of unique color,style combinations
        ccnt=ccnt+1;
        if ccnt>size(color_ord,1)
            ccnt=1;lcnt=lcnt+1;
            if lcnt>length(line_ord)
                lcnt=1;
            end
        end
    end
elseif makelstyle
    %If line color specified but not styles, assume they're all solid
    lstyle=repmat({'-'},1,nplots);
elseif makelcolor
    % If line styles specified but not colors, still cycle through
    ccnt=1;
    lcolor{1,nplots}='';
    for i=1:nplots
        lcolor{i}=color_ord(ccnt,:);
        ccnt=ccnt+1;
        if ccnt>size(color_ord,1)
            ccnt=1;
        end
    end
end
h_pl=zeros(1,nplots);
legstr{1,nplots}='';
for i=1:nplots
    if ax_in(i)==1
        ytemp=Y{i};
    else
        ytemp=(Y{i}-ranges(1,ax_in(i)))/(ranges(2,ax_in(i))-ranges(1,ax_in(i)))*(ranges(2,1)-ranges(1,1))+ranges(1,1);
    end
    
    h_pl(i)=plot(X{i},ytemp,...
        'Color',lcolor{i},...
        'Linestyle',lstyle{i},...
        'LineWidth',linew(i),...
        'markersize',msize(i),...
        'marker',mstyle{i},...
        'Markerfacecolor',lcolor{i});
        
    if isnumeric(options{5,2})
        legstr{i}=sprintf('Plot %i [Ax %i]',i,ax_in(i));
    else
        legstr{i}=options{5,2}{i};
    end
    
    
end
for i=nax:-1:1
    %Set  each axis to the default Ylims that matlab would impose
    %benifit from their robust logic
    set(h_ax(i),'Ylim',ranges(:,i)');
    axis(h_ax(i),'manual');
end
if isnumeric(options{7,2})
    %set(h_ax(1),'Xlim',options{7,2});
    xrange=get(h_ax(1),'Xlim');
	orginal=xrange;
    if ~isnan(options{7,2}(1))
        xrange(1)=options{7,2}(1);
    end
    if ~isnan(options{7,2}(2))
        xrange(2)=options{7,2}(2);
	end
	if xrange(2)==xrange(1)
		xrange=orginal;
	end
	
    set(h_ax(1),'Xlim',xrange);
end
if ~strcmpi(legloc,'none')
    %Place legend where requested
    legend(h_pl,legstr,'location',legloc)
end
%Create handles object, pass to output if requested
handles.ax=h_ax;
handles.parent=fh;
handles.xlab=h_x;
handles.ylab=h_y;
handles.title=h_t;
handles.plots=h_pl;
if nargout>0
    varargout{1}=handles;
end
%Final figure is ready, set visable on and pass full information to resize
%callback
set(fh,'units','normalized',...
    'visible','on',...
    'ResizeFcn',{@resize_callback,handles,marg_w,ax_w,nax,monitor_res,titpos})

set(fh,'color','w');

%Assign Zoom/pan callback now that we have complete handles struct
hzoom_ax1.ActionPostCallback={@zoom_event,handles,ranges,nax};
hpan_ax1.ActionPostCallback={@zoom_event,handles,ranges,nax};
if fighost
    %Since we just called function to edit fh, it should be current figure
    set(0, 'CurrentFigure', fh)
    set(fh,'CurrentAxes',h_ax(1))
    
    dcurse=datacursormode(fh);
    dcurse.UpdateFcn={@data_cursor_event,ranges,legstr,ax_in};
else
    
    set(0, 'CurrentFigure', parhand)
    set(parhand,'CurrentAxes',h_ax(1))
    dcurse=datacursormode(parhand);
    dcurse.UpdateFcn={@data_cursor_event,ranges,legstr,ax_in};
end
end
function [] = resize_callback(varargin)
%This function enforces the width of the axes and margins on resize event,
%gives extra space to plot space
req_pos=varargin{1}.Position;
handles=varargin{3};
marg_w=varargin{4};
ax_w=varargin{5};
nax=varargin{6};
monitor_res=varargin{7};
titpos=varargin{8};
%Calculate how much extra room we have after fitting margins and axes
fpos=req_pos.*[monitor_res(3) monitor_res(4) monitor_res(3) monitor_res(4)];
def_pos=[marg_w(1)+(nax-1)*ax_w marg_w(4) fpos(3)-marg_w(2)-(marg_w(1)+(nax-1)*ax_w) fpos(4)-marg_w(3)-marg_w(4)];
def_pos(4)=def_pos(4) - (titpos(2)-1)*def_pos(4);
%Shift each dummy axis to left of master axis for tapered axes look
for i=1:nax
    set(handles.ax(i),'Unit','Pixels')
    set(handles.ax(i),'Position',max([1 1 1 1],def_pos-[ax_w 0 0 0].*(i-1)))
    set(handles.ax(i),'Unit','normalized')
end
end
function [] = reject_rotate(varargin)
%This function resets view to normal after undesired rotate event
view(varargin{2}.Axes,[0 90])
end
function []= zoom_event(varargin)
%This function sets all the dummy axes when main axis range has been
%changed (from zoom or pan)
handles=varargin{3};
ranges=varargin{4};
nax=varargin{5};
%Pull the master axis range to scale with
new_rng=get(handles.ax(1),'Ylim');
%For each dummy axes, scale from original ranges to match current Ylim
for i=2:nax
    temprng(1,2)=(new_rng(2)-ranges(1,1))/(ranges(2,1)-ranges(1,1))*(ranges(2,i)-ranges(1,i))+ranges(1,i);
    temprng(1,1)=(new_rng(1)-ranges(1,1))/(ranges(2,1)-ranges(1,1))*(ranges(2,i)-ranges(1,i))+ranges(1,i);
    
    set(handles.ax(i),'Ylim',temprng);
end
end
function txt = data_cursor_event(varargin)
%This function translate the data plotted on axis 1's range to the range it belongs on
pointData=varargin{2};
pos=pointData.Position;
ranges=varargin{3};
legstr=varargin{4};
ax_in=varargin{5};
%Find line that data cursor is on
linename=get(pointData.Target,'DisplayName');
ind=find(strcmp(legstr,linename));
ind=ind(1); % in case there are duplicates, take only the first
%Which axis does this line below to?
axind=ax_in(ind);
%Scale the y position on ax1 to ax_ind
realy=(pos(2)-ranges(1,1))/(ranges(2,1)-ranges(1,1))*(ranges(2,axind)-ranges(1,axind))+ranges(1,axind);
txt = {['Axis: ',num2str(axind)],...
    ['X: ',num2str(pos(1))],...
    ['Y: ',num2str(realy)]};
end





end

