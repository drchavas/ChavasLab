%Vp_monthlyclimo.m - return monthly climatology of Vp
%Purpose: return monthly climatology of Vp
%
% Syntax:  [] = Vp_monthlyclimo(lat,lon,coastal_res)
%
% Inputs:
%   lat [deg N [-90,90]] - vector of latitude values
%   lon [deg E (-180,180]] - vector of corresponding longitude values
%   coastal_res [pts/deg] - resolution of coastline (gridpts per deg);
%       NOTE: higher resolution = more computing time (1 = decent coarse
%           res; 10 = decent high res)
%
% Outputs:
%   is_Ocean - 1 = ocean, 0 = land
%
% Example: 
%   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: pi_grid.mat
%
% Author: Dan Chavas (adapted from Brett Shoelson 9 Feb 2011 -- see above)
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: http://www.princeton.edu/~dchavas/
% 19 Aug 2014; Last revision:
% Revision history:

%------------- BEGIN CODE --------------

clear
clc
close('all')

addpath(genpath('~/Dropbox/Research/MATLAB/'));

%% USER INPUT %%%%%%%%%%%%%%%
Vp_file = sprintf('~/Dropbox/Research/TC_DATA/mpi/Vp_mat_allmonths.mat');
months_in = 8;
use_land_mask = 0;  %0=no mask; 1=land set to NaN
lat_min = -90;   %[-90,90] deg N
lat_max = 90;   %[-90,90] deg N
lon_min = 0;    %[0,360] deg E
lon_max = 360;  %[0,360] deg E

%%Testing
fig_vis = 'on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Saffir-Simpson category info for plots
SS_thresh = [17.5 33 42 49 58 70];   %[m/s]; wind speed thresholds
SS_colors = {[.9 .9 .9] [.5 .75 0] [1 1 0] [1 .5 0] [1 0 0] [.5 0 0]};
SS_str = {'TS/TD (17.5 m/s)' 'Cat1 (33 m/s)' 'Cat2 (42 m/s)' 'Cat3 (49 m/s)' 'Cat4 (58 m/s)' 'Cat5 (70 m/s)'};

for ii=1:length(months_in)

    month_in = months_in(ii);
    
    %% Load Vp data
    [Vp_month_case_mat,lat_Vp_mat,lon_Vp_mat] = Vp_monthlyclimo(Vp_file,month_in,use_land_mask);
    
    %% Subset data to desired region
    ii_keep = lat_Vp_mat <= lat_max & lat_Vp_mat >= lat_min & ...
        lon_Vp_mat <= lon_max & lon_Vp_mat >= lon_min;

    Vp_month_case_mat(~ii_keep) = NaN;
    lat_Vp_mat(~ii_keep) = NaN;
    lon_Vp_mat(~ii_keep) = NaN;

    %%2D interpolation to desired lat/lon pts
    %Vp_TC = griddata(lat_Vp_mat,lon_Vp_mat,Vp_month_case_mat,Lat_BT_TC,Lon_BT_TC);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TESTING: plot Vp %%%%%
    figure(1);set(gcf,'Visible',fig_vis)
    clf(1)
    set(gcf,'Units','centimeters');
    hpos = [0 0 30 30];
    set(gcf,'Position',hpos);
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',hpos);
    set(gcf,'PaperSize',hpos(3:4));

    set(gca,'position',[0.05    0    0.91    1]);


    lats_map = [lat_min lat_max];
    lons_map = [lon_min lon_max];

    if(lons_map(end)<30 && lons_map(end)>0)
        lons_map(end) = -4;  %has issues when going over greenwich meridion
    end

    ax = worldmap(lats_map,lons_map);
    temp = get(ax,'Position');

    %%Add land mask
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [.5 .5 .5])


    [C_Vp,h_Vp] = contourm(lat_Vp_mat,lon_Vp_mat,Vp_month_case_mat,SS_thresh,'LineWidth',3);

    temp = get(h_Vp,'Children');    %first contour lines, then fills
    numthresh = length(SS_thresh);
    for pp=1:numthresh

        %%Update contour color (ordered from highest to lowest value!)
        %%NOTE: all contour thresholds will have children whether or not they actually appear in the plot
        set(temp(numthresh+1-pp),'Color',SS_colors{pp});

        drawnow     %need this or else fill won't work (from MATLAB forum! http://www.mathworks.com/matlabcentral/answers/68373-problem-changing-facecolor-in-contourm)
    end
    drawnow

    tightmap

    %%Add a Legend for Saffir-Simpson scale plus a scale for size
    SS_leg = clegendm(C_Vp,h_Vp,SS_str);

    %%Adjust legend position
    temp = get(SS_leg,'position');
    set(SS_leg,'position',temp+[0 0.2 0 0]);

    title(sprintf('V_p month %i',month_in))

    %%SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_filename = sprintf('Vp_monthlyclimo_test_%s.pdf',sprintf('%2.2d',month_in))
    saveas(gcf,plot_filename,'pdf')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%------------- END OF CODE --------------