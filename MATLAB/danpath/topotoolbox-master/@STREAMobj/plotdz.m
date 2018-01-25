function h = plotdz(S,DEM,varargin)

% plot upstream distance version elevation of a stream network
%
% Syntax
%
%     plotdz(S,DEM)
%     plotdz(S,DEM,pn,pv,...)
%     plotdz(S,nal,pn,pv,...)
%     h = ...
%
% Description
%
%     Plot stream distance from the outlet versus elevation. 
%
% Input arguments
%
%     S      instance of STREAMobj
%     DEM    digital elevation model (GRIDobj)
%     nal    node attribute list (as returned by various STREAMobj
%            methods, e.g. STREAMobj/streamorder, STREAMobj/gradient)
%
%     Parameter name/value pairs {default}
%
%     'annotation':  {[]} ix      
%     vector with linear indices of locations into the DEM. The cells
%     referenced by ix must be part of the stream network. Use snap2stream
%     to adjust locations. Annotation is achieved with vertical arrows.
%
%     'annotationtext': cell array of strings
%     if annotated, a cell array of strings can be added to the vertical 
%     arrows
%
%     'distance': {S.distance}
%     node attribute list with custom distances (see STREAMobj/distance) or
%     STREAMobj (see function distance(S,S2))
%
%     'dunit': {'m'} 'km'
%     distance unit. plotdz assumes that distance is given in meter. 
%
%     'doffset': {0}
%     add an offset (scalar) to the distance from outlet to shift the 
%     x-axis of the plot.
%
%     'linewidth', 1
%     scalar that specifies the width of the line. 
%
%     'color': {'b'}
%     line colors as provided to plot. Alternatively, you can supply a node
%     attribute list (nal). The line will then have varying colors based on 
%     nal and h will be a surface object.
%
%     if 'color' is a node attribute list, then following parameter
%     name/values apply
%
%     'colormethod'    {'line'} or 'surface'
%     lines with variable colors can be obtained by different methods.
%     Before introduction of the new graphics system in Matlab, using the
%     edges of surfaces was the standard hack. Since 2014b and the new
%     graphics system, lines have an undocumented edge property. Modifying
%     this edge property results in much smoother and more beautiful lines,
%     but may be bugged in newer versions.
%
%     'colormap'  {'parula'}
%     string that identifies a known colormap (e.g. 'jet','landcolor')
%
%     'colorbar'  {true} or false
%     true adds a colorbar. Note that colors created with the undocumented
%     method can not be changed afterwards by 'caxis'.
%
%     'cbarlabel' {''}
%     string to label colorbar
%     
%
% Output arguments
%
%     h     handle to the line handle. h will be a surface handle if color 
%           is set to a nal and colormethod is surface.
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'mex',true,'preprocess','carve');
%     S  = STREAMobj(FD,flowacc(FD)>1000);
%     S  = klargestconncomps(S);
%     plotdz(S,DEM)
%
% Example 2 (colored plot)
%     
%     z = imposemin(S,DEM);
%     g = gradient(S,z);
%     plotdz(S,DEM,'color',g)
%
% See also: STREAMobj, STREAMobj/plot
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2017

nrnodes = numel(S.x);
ax      = gca;

colororderindex = mod(ax.ColorOrderIndex, size(ax.ColorOrder,1));
if colororderindex==0; colororderindex=size(ax.ColorOrder,1); end

% check input
p = inputParser;         
p.FunctionName = 'plotdz';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'DEM', @(x) isa(x,'GRIDobj') || numel(x) == nrnodes);
addParameter(p,'annotation',[])
addParameter(p,'color',ax.ColorOrder(colororderindex,:));
addParameter(p,'annotationtext',{});
addParameter(p,'distance',[],@(x) isnal(S,x) || isa(x,'STREAMobj'));
addParameter(p,'dunit','m',@(x) ischar(validatestring(x,{'m' 'km'})));
addParameter(p,'doffset',0,@(x) isscalar(x));
addParameter(p,'colormap','parula');
addParameter(p,'linewidth',1);
addParameter(p,'colormethod','line');
addParameter(p,'colorbar',true);
addParameter(p,'cbarlabel','');


parse(p,S,DEM,varargin{:});
S   = p.Results.S;
DEM = p.Results.DEM;

if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    zz = getnal(S,DEM);
elseif isnal(S,DEM);
    zz = DEM;
else
    error('Imcompatible format of second input argument')
end


if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
end

% get dynamic properties of S
order    = S.orderednanlist;

if isempty(p.Results.distance);
    dist = S.distance;
else
    if isa(p.Results.distance,'STREAMobj');
        dist = distance(S,p.Results.distance);
    else
        dist = p.Results.distance;
    end
end


switch lower(p.Results.dunit);
    case 'km'
        dist = dist/1000;
end
        

% apply distance offset
dist = dist + p.Results.doffset;

I     = ~isnan(order);
d     = nan(size(order));
d(I)  = dist(order(I));
z     = nan(size(order));
z(I)  = zz(order(I));

% plot
if ~isnal(S,p.Results.color)
    ht = plot(ax,d,z,'-','Color',p.Results.color,'LineWidth',p.Results.linewidth);
else
    meth = validatestring(p.Results.colormethod,{'line','surface'});
    switch meth
        case 'line'
            %% Plotting colored lines using undocumented Edges property
            % see here:
            % http://undocumentedmatlab.com/blog/plot-line-transparency-and-color-gradient
            ht = plot(ax,d,z,'-');
            c     = zeros(size(order,1),3);
            minc  = min(p.Results.color);
            maxc  = max(p.Results.color);
            
            cmap    = colormap(p.Results.colormap)*255;
            cmapix  = linspace(minc,maxc,size(cmap,1));
            
            col   = interp1(cmapix,cmap,p.Results.color);
            c(I,:) = col(order(I),:);
            c     = c';
            c     = [c;zeros(1,size(c,2))+200];
            c     = uint8(c);
            % nans must be removed
            c     = c(:,I);
            
            ht.LineWidth = p.Results.linewidth;
            % seems that the line must be first drawn ...
            drawnow
            % ... to be colored.
            set(ht.Edge, 'ColorBinding','interpolated', 'ColorData',c);
            
            if p.Results.colorbar
                cc = colorbar;
                caxis([minc maxc]);
            end
            
        otherwise
            %% Plotting colored lines using surface
            
            dummy = z*0;
            c     = nan(size(order));
            c(I)  = p.Results.color(order(I));
            colormap(p.Results.colormap)
            ht = surface([d d],[z z],[dummy dummy],[c c],...
                'facecolor','none',...
                'edgecolor','flat',...
                'linewidth',p.Results.linewidth,...
                'parent',ax);
            if p.Results.colorbar
                cc = colorbar;
            end
            
            
    end
    if p.Results.colorbar && ~isempty(p.Results.cbarlabel);
        cc.Label.String = p.Results.cbarlabel;
    end
end

xlabel(['Distance upstream [' lower(p.Results.dunit) ']'])
ylabel('Elevation [m]')

%% Annotation
if ~isempty(p.Results.annotation);
    ix = p.Results.annotation;
    hold on
    [Lia,Locb] = ismember(ix,S.IXgrid);


    if any(~Lia)
        error('TopoToolbox:STREAMobj',...
            'Some of the annotations are not located on the stream network')
    end
    
    annd = dist(Locb);
    
    if isa(DEM,'GRIDobj')
        annz = DEM.Z(S.IXgrid(Locb));
    else
        annz = zz(Locb);
    end
    
    if ~isempty(p.Results.annotationtext);
        c = p.Results.annotationtext;
        addtext = true;
    else
        addtext = false;
    end
    
    
    for r = 1:numel(ix);
        if addtext
            annotext = [c{r} '\newline\downarrow'];
        else
            annotext = '\downarrow ';
        end
        
        text(ax,'Position',[annd(r), annz(r)],...
             'String', annotext,...
             'VerticalAlignment','bottom',...
             'FontWeight','bold');
    end
    hold off
end


if nargout == 1;
    h = ht;
end
end
