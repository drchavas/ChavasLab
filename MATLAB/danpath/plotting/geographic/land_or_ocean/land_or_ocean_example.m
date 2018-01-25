%land_or_ocean_example.m
%example file testing land_or_ocean.m

clear
clc
close('all')


%% Ex1
lat = -90:10:90;
lon = -180:20:180;
coastal_res = 5;
make_plot = 1;  %0 = no plot, 1 = plot
[isOcean1] = land_or_ocean(lat,lon,coastal_res,make_plot)
figure(1)
saveas(gcf,'ex1.pdf','pdf')


%% Ex2
lat = [20 30];
lon = [-70 -100];
coastal_res = 10;    %note: if use 1, both on land; if use >4, first pt is water
make_plot = 1;  %0 = no plot, 1 = plot
[isOcean2] = land_or_ocean(lat,lon,coastal_res,make_plot)
figure(2)
saveas(gcf,'ex2.pdf','pdf')
