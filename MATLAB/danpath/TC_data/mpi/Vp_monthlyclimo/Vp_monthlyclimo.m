%Vp_monthlyclimo.m - return monthly climatology of Vp
%Purpose: return monthly climatology of Vp
%
% Syntax:  [Vp_month_case_mat,lat_Vp_mat,lon_Vp_mat] = Vp_monthlyclimo(Vp_file,month_in,use_land_mask)
%
% Inputs:
%   Vp_file - path to Vp file pi_grid.mat
%   month_in [1-12] - number of desired month
%   use_land_mask - 0=no masking; 1=set land to NaN
%
% Outputs:
%   Vp_month_case_mat - potential intensity matrix
%   lat_Vp_mat - latitude matrix
%   lon_Vp_mat - longitude matrix
%
% Example: [Vp_month_case_mat,lat_Vp_mat,lon_Vp_mat] = Vp_monthlyclimo(Vp_file,month_in,use_land_mask);
%   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: Vp_mat_allmonths.mat
%
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: http://www.princeton.edu/~dchavas/
% 19 Aug 2014; Last revision:
% Revision history:

%------------- BEGIN CODE --------------

function [Vp_month_case_mat,lat_Vp_mat,lon_Vp_mat] = Vp_monthlyclimo(Vp_file,month_in,use_land_mask)

if(nargin==2)
    use_land_mask = 0; %default: no land mask
end

%% Load Vp data
listOfVariables = {'Vp_mat_allmonths','lat_Vp_mat','lon_Vp_mat','time_day_Vp_mat'};
load(Vp_file,listOfVariables{:})
sprintf('Vp monthly data for entire year loaded from %s',Vp_file)

Vp_month_case_mat = Vp_mat_allmonths(:,:,month_in);

%% If desired, set values over land to NaN
if(use_land_mask)
    coastal_res = 1;    %[pts/deg]
    make_plot = 0;
    [isOcean] = land_or_ocean(lat_Vp_mat,lon_Vp_mat,coastal_res,make_plot);

    Vp_month_case_mat(~isOcean) = NaN;
end

%%2D interpolation to desired lat/lon pts
%Vp_TC = griddata(lat_Vp_mat,lon_Vp_mat,Vp_month_case_mat,Lat_BT_TC,Lon_BT_TC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTING: plot Vp %%%%%
%{
%% Saffir-Simpson category info for plots
SS_thresh = [17.5 33 42 49 58 70];   %[m/s]; wind speed thresholds
SS_colors = {[.6 .6 .6] [.9 .9 .9] [.5 .75 0] [1 1 0] [1 .5 0] [1 0 0] [.5 0 0]};
SS_str = {'TS/TD (17.5 m/s)' 'Cat1 (33 m/s)' 'Cat2 (42 m/s)' 'Cat3 (49 m/s)' 'Cat4 (58 m/s)' 'Cat5 (70 m/s)'};


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
    set(temp(numthresh+1-pp),'Color',SS_colors{pp+1});

    drawnow     %need this or else fill won't work (from MATLAB forum! http://www.mathworks.com/matlabcentral/answers/68373-problem-changing-facecolor-in-contourm)
end
drawnow

tightmap

%%Add a Legend for Saffir-Simpson scale plus a scale for size
SS_leg = clegendm(C_Vp,h_Vp,SS_str);

%%Adjust legend position
temp = get(SS_leg,'position');
set(SS_leg,'position',temp+[.2 0.15 0 0]);

%%SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_filename = sprintf('Vp_monthlyclimo_test_%s.pdf',sprintf('%2.2d',month_in))
print('-dpdf','-r300',plot_filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%------------- END OF CODE --------------