function [ hcb ] = cbar_custom(varargin)
%[ hcb ] = cbar_custom( varargin)
%   function to draw a custom colormap WITHIN the plot bounds at a location specified with respect
%   to the current axess (n.b. to plot outside axes, can just add another
%   set of axes and specify location w.r.t. those axes
% 
% USAGE:
% 
%  cbar_custom     
%       makes a default colourbar in the southwest
% 
%  cbar_custom(...,'location',[Xleft Xright Ybottom Ytop])
%        default is [0.05 0.5 0.05 0.15]*axlim
%  cbar_custom(...,'lims',[Vmin Vmax])                  
%        default is [0 1]
%  cbar_custom(...,'tickvals',[tick_position_vector])      
%        default 11 ticks evenly spaced
%  cbar_custom(...,'ticklabs',[tick_position_vector])      
%        has to be the same number of elements as in the tickvals
%        default 11 ticks evenly spaced
%  cbar_custom(...,'cmap',colourmap)
%        default is jet
%  cbar_custom(...,'ncols',n_of_colour_bins)            
%        default is 50
%  cbar_custom(...,'tickside','left/right/top/bottom')  
%        adjusts the position of the labels and ticks accordingly
%        default is 'bottom'
%  cbar_custom(...,'title',string)
%        default is none
%  cbar_custom(...,'LineWidth',lw)
%        default is 2
%  cbar_custom(...,'FontSize',fsz)
%        default is 10
%  cbar_custom(...,'FontWeight',fwt)
%        default is 'normal'
%  cbar_custom(...,'FontCol',fcl)
%        default is 'k'
%  cbar_custom(...,'FontAngle',fang)
%        default is 0
%  cbar_custom(...,'interpreter',string)
%        default is 'tex'
%  cbar_custom(ax,...) 
%       makes colorbar in the axes "ax", where this is a handle of the axes
% 
%  hcb = cbar_custom(...)    
%       to get the handle of the colourbar (actually a vector of handles)
%       N.B. not working yet
% 
%  Written by Z. Eilon      05/2015

%% defaults
% location in SW, horizontal, ticks below
axlim = axis;
cscXle = axlim(1) + 0.05*diff(axlim(1:2)); % 5% of way from left 
cscXri = axlim(1) + 0.5*diff(axlim(1:2)); % 20% of way from left 
cscYbo = axlim(3) + 0.05*diff(axlim(3:4)); % 5% of way from bottom 
cscYto = axlim(3) + 0.15*diff(axlim(3:4)); % 8% of way from bottom 
tickside = -1;
orie = 'horizontal';
interpstr = 'tex';
titlestr = '';
% lims
Vmi = 0;
Vma = 1;
ncols = 50;
cscntick = 11;
cmap = jet;
fsz = 12;
fwt = 'normal';
fang = 0;
lw = 2;
fcl = 'k';
ftix = 0.15; %fraction of width of cbar to have labels centred
vtick = linspace(Vmi,Vma,cscntick); % vector of values of ticks
ticklab = [];

%% read inputs
iarg = 1;
while iarg <= length(varargin)
    if iarg==1
        if ishandle(varargin{iarg})
            axes(varargin{iarg});
            iarg = iarg+1;
            continue
        end
    end
            
    switch lower(varargin{iarg})
        case {'location','loc'}
            cscXle = varargin{iarg+1}(1);
            cscXri = varargin{iarg+1}(2);
            cscYbo = varargin{iarg+1}(3);
            cscYto = varargin{iarg+1}(4);
        case {'lims'}
            Vmi = min(varargin{iarg+1});
            Vma = max(varargin{iarg+1});
            if Vmi==Vma, error('Give separate minimum and maximum limits'), end
        case {'tickvals'}
            vtick = varargin{iarg+1};
        case {'ticklabs'}
            ticklab = varargin{iarg+1};
        case {'cmap'}
            cmap = varargin{iarg+1};
        case {'ncols'}
            ncols = varargin{iarg+1};
        case {'tickside'}
            if strcmp(varargin{iarg+1},'left'),  tickside = -1; orie = 'vertical'; end
            if strcmp(varargin{iarg+1},'right'), tickside = +1; orie = 'vertical'; end
            if strcmp(varargin{iarg+1},'bottom'),tickside = -1; orie = 'horizontal'; end
            if strcmp(varargin{iarg+1},'top'),   tickside = +1; orie = 'horizontal'; end
        case {'title'}
            titlestr = varargin{iarg+1};
        case {'fontsize'}
            fsz = varargin{iarg+1};
            if ~isnumeric(fsz), error('FontSize must be numeric'), end
        case {'fontweight'}
            fwt = varargin{iarg+1};
            if isnumeric(fwt), error('FontWeight must be string'), end
        case {'fontcol'}
            fcl = varargin{iarg+1};
            if isnumeric(fcl), error('FontCol must be string'), end
        case {'fontang','fontangle'}
            fang = varargin{iarg+1};
            if ~isnumeric(fang), error('FontAngle must be numeric'), end
        case {'linewidth'}
            lw = varargin{iarg+1};
            if ~isnumeric(fsz), error('LineWdith must be numeric'), end
        case {'interpreter'}
            interpstr = varargin{iarg+1}; 
    end
    iarg = iarg+2;
end

if isempty(ticklab), ticklab = num2str(vtick(:)); end

vtick = vtick(:)'; % must be row

%% calcs - calculate corner positions for each colour bin, depending on orientation
cscv = linspace(Vmi,Vma,ncols+1); % vector of values at edges of colour bins
dX = cscXri-cscXle;
dY = cscYto-cscYbo;

if strcmp(orie,'horizontal')
    cscXv = linspace(cscXle,cscXri,ncols+1); % vector of X positions at edges of colour bins
    cscYv = repmat([cscYbo cscYto],1,ncols); cscYv = cscYv(1:length(cscXv)); % vector of Y positions at edges of colour bins
    % ticks
    csctiXv = interp1(cscv,cscXv,vtick);
    csctiYv = (mean([cscYbo cscYto]) + tickside*dY*(0.65+ftix)) * ones(size(csctiXv));
    % title
    titX = mean([cscXle cscXri]);
    titY = mean([cscYbo cscYto]) - tickside*dY*(0.5 + 1.5*ftix);
    titrot = 0;
    % all strings
    if tickside == 1
        vertal = 'bottom'; titvertal = 'top'; 
    elseif tickside == -1
        vertal = 'top';    titvertal = 'bottom'; 
    end
    horizal = 'center'; tithorizal = 'center';
elseif strcmp(orie,'vertical')
    cscYv = linspace(cscYbo,cscYto,ncols+1); % vector of X positions at edges of colour bins
    cscXv = repmat([cscXle cscXri],1,ncols); cscXv = cscXv(1:length(cscYv)); % vector of Y positions at edges of colour bins
    % ticks
    csctiYv = interp1(cscv,cscYv,vtick);
    csctiXv = (mean([cscXle cscXri]) + tickside*dX*(0.5+2*ftix))  * ones(size(csctiYv));
    % title
    titX = mean([cscXle cscXri]) - tickside*dX*(0.5 + 6*ftix);
    titY = mean([cscYbo cscYto]) ;
    titrot = 90;
    % all strings
    if tickside == 1
        horizal = 'left';  tithorizal = 'center'; 
    elseif tickside == -1
        horizal = 'right'; tithorizal = 'center'; 
    end
    vertal = 'middle'; titvertal = 'middle';
end


%% draw
hcb = hggroup;
hold on
%% colour bins
for ii = 1:ncols
    xx = [cscXv(ii),cscXv(ii+1),cscXv(ii+1),cscXv(ii),cscXv(ii)];
    yy = [cscYv(ii),cscYv(ii),cscYv(ii+1),cscYv(ii+1),cscYv(ii)];
    hp(ii) = patch(xx,yy,colour_get(mean(cscv(ii:ii+1)),Vma,Vmi,cmap)','Parent',hcb);
    set(hp,'LineStyle','none')
end
%% ticks and labels
% labels
text(csctiXv',csctiYv',ticklab,...
        'HorizontalAlignment',horizal,'VerticalAlignment',vertal,...
        'FontSize',fsz,'FontWeight',fwt,'interpreter',interpstr,'color',fcl,'Parent',hcb,'rotation',fang)
% ticks
if strcmp(orie,'horizontal')
plot([1;1]*csctiXv,cscYbo + [0 0.2]*dY,'k','LineWidth',lw/1.5,'Parent',hcb)
plot([1;1]*csctiXv,cscYto - [0 0.2]*dY,'k','LineWidth',lw/1.5,'Parent',hcb)
elseif strcmp(orie,'vertical')
plot(cscXle + [0 0.2]*dX,[1;1]*csctiYv,'k','LineWidth',lw/1.5,'Parent',hcb)
plot(cscXri - [0 0.2]*dX,[1;1]*csctiYv,'k','LineWidth',lw/1.5,'Parent',hcb)
end
%% title
text(titX,titY,titlestr,...
        'HorizontalAlignment',tithorizal,'VerticalAlignment',titvertal,'rotation',titrot,...
        'FontSize',1.3*fsz,'FontWeight',fwt,'interpreter',interpstr,'color',fcl,'Parent',hcb)

plot([cscXle,cscXri,cscXri,cscXle,cscXle],[cscYbo,cscYbo,cscYto,cscYto,cscYbo],'-k','LineWidth',lw,'Parent',hcb)

hold off
end