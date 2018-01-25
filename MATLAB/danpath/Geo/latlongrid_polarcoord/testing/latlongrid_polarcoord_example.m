%latlongrid_polarcoord_example.m
%Purpose: Return lat/lons for a polar grid
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 2015-10-15; Last revision:

%------------- BEGIN CODE --------------

clc
clear
close all

addpath('../');

dr_inner = 50*1000;   %[m]
r1_inner = dr_inner; %hypothetical inner-core storm inner edge (m)
r2_inner = 100*1000; %hypothetical inner-core storm outer edge (m)
dr_outer = 100*1000;   %[m]
r1_outer = dr_outer; %hypothetical outer annulus environmental inner edge (m)
r2_outer = 300*1000; %hypothetical  outer annulus environmental outer edge (m)
dth = pi/4; %[rad]
lat_cent = 0;


[dlon_inner,dlat_inner,dist_inner] = latlongrid_polarcoord(r1_inner,r2_inner,dr_inner,dth,lat_cent);
[dlon_outer,dlat_outer,dist_outer] = latlongrid_polarcoord(r1_outer,r2_outer,dr_outer,dth,lat_cent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTING
hh=figure(1)
clf(1)
%%Position/size
set(hh,'Units','centimeters');
hpos = [0 0 15 15];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.08    0.07    0.88    0.9]);

scatter(dlon_outer,dlat_outer,70,dist_outer/1000,'x')
hold on
scatter(dlon_inner,dlat_inner,70,dist_inner/1000,'o')
xlabel('dlon [deg]')
ylabel('dlat [deg]')
colorbar
title('color = distances [km]')
axis([-5 5 -5 5])

%% SAVE %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('latlongrid_polarcoord_example_lat%i.pdf',round(lat_cent))
saveas(gcf,plot_filename,'pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

