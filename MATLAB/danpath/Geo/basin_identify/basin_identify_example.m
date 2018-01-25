%basin_identify_example.m -- Example using basin_identify.m with plot

%------------- BEGIN CODE --------------

clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Polygon file
file_polygons = 'example_polygons.mat';

%%Point of interest
N_pts = 100;
lon_pts = 360*rand(N_pts,1);   %[deg E] [0,360)
lat_pts = 90.*sign(rand(N_pts,1)-.5).*rand(N_pts,1);    %[deg N] [-90,90]


%%Make a plot?
make_plot = 1;  %1: makes a plot; ow: no plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listOfVariables = {'basin_strs','lon_poly','lat_poly'};
load(file_polygons,listOfVariables{:});

[basin_pts] = basin_identify(lon_pts,lat_pts,lon_poly,lat_poly,basin_strs,make_plot);


%------------- END OF CODE --------------