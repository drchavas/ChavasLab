%basin_identify.m -- identify basins in which set of points are located
%Purpose: Given lat/lon points and basin polygons/names, return basin names
%
% Syntax: [basin_pts] = basin_identify(lon_pts,lat_pts,lon_polys,lat_polys,basin_strs,make_plot)  
%
% Inputs:
%   lon_pts [deg E] - longitude of input points [0,360)
%   lat_pts [deg N] - latitude of input points [-90,90]
%   lon_polys [deg E] - cell array of longitude values of polygon vertices
%       (must be ordered either clockwise or counterclockwise)    
%   lat_polys [deg N] - cell array of latitude values of polygon vertices
%       (must be ordered either clockwise or counterclockwise)    
%   basin_strs [] - cell array of basin labels (one per cell in lon_polys)
%   make_plot [0/1] - 1 = make geo plot of basins and points with labels
%
% Outputs:
%   basin_pts [] - cell array of basin labels at input pts; 'none' if
%       outside all input polygons
%
% Example: [basin_pts] = basin_identify(lon_pts,lat_pts,lon_polys,lat_polys,basin_strs,make_plot)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Notes: Lines on output plot do not align perfectly with meridians, so it
% may look like there are points that are in the wrong place. They are not.

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 5 Nov 2014; Last revision:

%------------- BEGIN CODE --------------

function [basin_pts] = basin_identify(lon_pts,lat_pts,lon_polys,lat_polys,basin_strs,make_plot)

if(nargin < 6)  
    
    %%Default: no plot
    make_plot = 0;
    
    if(nargin < 5)
        
        %%Default: number the basins in order
        temp = 1:length(lon_polys);
        basin_strs = strread(num2str(temp),'%s')';
        
    end
    
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTING: boundaries plot %%%%%%%%%%%%%%%%%%%%%%%%%%
if(make_plot)
    fig2 = figure(2);
    clf(fig2)
    set(fig2,'Units','centimeters');
    hpos = [0 0 35 25];
    set(fig2,'Position',hpos);
    set(fig2,'PaperUnits','centimeters');
    set(fig2,'PaperPosition',hpos);
    set(fig2,'PaperSize',hpos(3:4));

    %% Plot land
    ax = worldmap('World')
    getm(ax,'MapProjection')
    setm(ax, 'Origin', [0 180 0])
    % land = shaperead('landareas', 'UseGeoCoords', true);
    % geoshow(ax, land, 'FaceColor', [.5 .5 .5])
    load coast
    plotm(lat, long)
    hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


basin_pts = cell(size(lon_pts));
for ii=1:length(basin_strs)
   
    basin_str = basin_strs{ii};
    lon_polys_temp = lon_polys{ii};
    lat_polys_temp = lat_polys{ii};
    lon_polys_temp = [lon_polys_temp lon_polys_temp(1)];   %to wrap around
    lat_polys_temp = [lat_polys_temp lat_polys_temp(1)];
    
    for jj=1:length(lon_polys_temp)-1
        
        lons_temp = lon_polys_temp(jj:jj+1);
        lats_temp = lat_polys_temp(jj:jj+1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TESTING: boundaries plot %%%%%%%%%%%%%%%%%%%%%%%%%%
        if(make_plot)
            h_line = geoshow(lats_temp,lons_temp,'DisplayType','line');
            set(h_line,'LineWidth',3,'Color','r')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TESTING: boundaries plot %%%%%%%%%%%%%%%%%%%%%%%%%%
    if(make_plot)
        %%Plot label
        textm(mean([max(lat_polys_temp) min(lat_polys_temp)]),...
            mean([max(lon_polys_temp) min(lon_polys_temp)]),basin_str,...
            'FontSize',20,'FontWeight','bold')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    i_basin = inpolygon(lon_pts,lat_pts,lon_polys_temp,lat_polys_temp);
    if(sum(i_basin)>0)
        basin_pts(i_basin) = {basin_str};
    end
    clear i_basin
    
end

basin_pts(cellfun(@isempty,basin_pts)) = {'None'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTING: boundaries plot %%%%%%%%%%%%%%%%%%%%%%%%%%
if(make_plot)
    h_pt = geoshow(lat_pts,lon_pts,'DisplayType','Point');
    set(h_pt,'Marker','x','Color','g','MarkerSize',15)
    textm(lat_pts,lon_pts,basin_pts)
    title('Red X = test points')
    
    sprintf('Note: Lines on output plot do not align perfectly with meridians')
    sprintf('so it may look like there are points that are in the wrong place. They are not.')
    
    %% Save plot %%%%%%%%%%%%%%%%%%
    plot_filename = sprintf('basin_identify_example.pdf');
    saveas(gcf,sprintf('%s',plot_filename),'pdf')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- END CODE --------------