%netcdf_plot.m -- template to read in and plot netcdf data

clear
clc
close all

set(0,'defaultaxesfontsize',24,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',3,'DefaultAxesFontName','Helvetica')
            
addpath(genpath('~/Dropbox/Research/MATLAB/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%General
ncdir = sprintf('./');
ncfile = 'air.sfc.mon.ltm.nc';

%%Parameters
month_in = 10;

%%Output
figdir_out = sprintf('./');
topo_min = 100; %[m]; minimum topography contour
dtopo = 100; %[m]; topography contour interval
topo_max = 1000; %[m]; maximum topography contour
clabs = 250:5:310; %[K]; contour interval for temperature
lon_min = -150; %[deg E]
lon_max = -50;  %[deg E]
lat_min = 10;   %[deg N]
lat_max = 60;   %[deg N]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncpath = sprintf('%s/%s',ncdir,ncfile);

%% Read in NetCDF data
time = ncread(ncfile,'time');

Tsfc = ncread(ncfile,'air',[1 1 month_in],[Inf Inf 1]);
Tsfc(Tsfc<0) = NaN; %remove missing values

lat_in = ncread(ncfile,'lat');
lon_in = ncread(ncfile,'lon');

assert(isequal(size(lon_in),size(lat_in)),'lat and lon matrices are not the same size')
assert(isequal(size(Tsfc),size(lat_in)),'lat and Tsfc matrices are not the same size')

%% Block out data east of Greenwich Meridien (fixes plotting issues)
Tsfc(lon_in>0)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
fig1 = figure(1);%set(gcf,'Visible','off')
clf(fig1)
set(fig1,'Units','centimeters');
hpos = [0 0 30 30]; %[cm]; dimensions of the entire plot
%hpos = [0 0 60 30]; 
set(fig1,'Position',hpos);
set(fig1,'PaperUnits','centimeters');
set(fig1,'PaperPosition',hpos);
set(fig1,'PaperSize',hpos(3:4));

%%If single-panel plot only: set fractional dimensions of plot box
set(gca,'position',[0.09    0.09    0.86    0.87]);

%%Multi-panel plot: set spacing (buffers height, buffers width, buffers between subplots) 
%subplot = @(m,n,p) subtightplot(m, n, p, [0.1 0.08], [0.1 0.05], [0.06 0.03]);
%%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)


%% Plot global coastline
load coast
plot(long,lat,'Color',[.9 .9 .9])
hold on

%% Plot US states
states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
lat_states = states.Latitude;
lon_states = states.Longitude;
plot(lon_states,lat_states,'Color',[.7 .7 .7],'LineWidth',1)

%% Plot data
clear hpl hnum input_legend
hnum = 0;

hold on
hnum = hnum + 1;
[c,h] = contour(lon_in,lat_in,Tsfc,clabs,'LineWidth',2);
colorbar
clabel(c,h,'LabelSpacing',200)
caxis([min(clabs) max(clabs)])

%% Plot aesthetics
lon_pl_min = lon_min-5;
lon_pl_max = lon_max+5;
lat_pl_min = lat_min-5;
lat_pl_max = lat_max+5;
axis([lon_pl_min lon_pl_max lat_pl_min lat_pl_max]) %only show this region of the Earth in your plot
xlabel('longitude [deg E]')
ylabel('latitude [deg N]')
input_title = sprintf('Monthly-mean T_{sfc} [K], month %2.2d (range = [%3.1f, %3.1f] K)',month_in,min(Tsfc(:)),max(Tsfc(:)));
title(input_title,'FontSize',12)

    
%% Save plot %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('netcdf_plot.jpg');
saveas(gcf,sprintf('%s/%s',figdir_out,plot_filename),'jpeg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

