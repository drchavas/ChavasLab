%map_mtemplate.m -- template for plotting a basic map without mapping tools
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%
% Outputs:
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Dan Chavas
% EAPS, Purdue University
% email: drchavas@gmail.com
% Website: -
% 2016-02-19; Last revision:

%------------- BEGIN CODE --------------


clear
clc
close all

addpath(genpath('~/Dropbox/Research/MATLAB/'));
set(0,'defaultaxesfontsize',24,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',3,'DefaultAxesFontName','Helvetica')

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%General
datadir_in = sprintf('.');
file1_in = 'map_mtemplate_data.mat';

%%Parameters

%%Plotting
figdir_out = sprintf('.');
topo_min = 100; %[m]; minimum topography contour
dtopo = 100; %[m]; topography contour interval
topo_max = 1000; %[m]; maximum topography contour
pl_clr = 'r';   %color of markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
path_in = sprintf('%s/%s',datadir_in,file1_in);
listOfVariables = {
    'lat_in','lon_in'
    };
load(path_in, listOfVariables{:});
sprintf('Loading random data from %s',path_in)


%% Calculate a few basic things
lon_min = min(lon_in);
lon_max = max(lon_in);
lat_min = min(lat_in);
lat_max = max(lat_in);
N_pts = length(lon_in);

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
set(gca,'position',[0.09    0.09    0.89    0.88]);

%%Multi-panel plot: set spacing (buffers height, buffers width, buffers between subplots) 
%subplot = @(m,n,p) subtightplot(m, n, p, [0.1 0.08], [0.1 0.05], [0.06 0.03]);
%%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)


%% Plot global coastline
load coast
plot(long,lat,'k')
hold on

%% Plot US states
states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
lat_states = states.Latitude;
lon_states = states.Longitude;
plot(lon_states,lat_states,'b')

%% Plot global topography
load topo
lon_topo = -179.5:179.5;                                % longitude
lat_topo = -89.5:89.5;                               % latitude
topo_temp = topo(:,1:end/2);
topo_temp2 = topo(:,end/2+1:end);
topo = [topo_temp2 topo_temp];

topolevs = [topo_min:dtopo:topo_max];
[C_pl,~] = contour(lon_topo,lat_topo,topo,topolevs);
clrmap = colormap('gray');
colormap(flip(clrmap,1));
%h_clab = clabel(C_pl,'FontSize',8,'Rotation',0,'Color',[.3 .3 .3]); % This will produce only one label per contour.

%% Plot data
clear hpl hnum input_legend
hnum = 0;

hnum = hnum + 1;
hpl(hnum) = plot(lon_in,lat_in,'.','MarkerSize',50,'Color',pl_clr);
input_legend{hnum} = sprintf('random points (N=%i)',N_pts);

%% Plot aesthetics
lon_pl_min = lon_min-5;
lon_pl_max = lon_max+5;
lat_pl_min = lat_min-5;
lat_pl_max = lat_max+5;
axis([lon_pl_min lon_pl_max lat_pl_min lat_pl_max]) %only show this region of the Earth in your plot
xlabel('longitude [deg E]')
ylabel('latitude [deg N]')
file1_in_str = strrep(file1_in,'_','\_');
input_title = sprintf('Totally random data (file: %s; lon range = [%2.2d,%2.2d]; lat range = [%2.2d,%2.2d]); N=%i datapoints',file1_in_str,lon_min,lon_max,lat_min,lat_max,N_pts);
title(input_title,'FontSize',12)
legend(hpl,input_legend);% legend boxoff

    
%% Save plot %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('map_mtemplate.jpg');
saveas(gcf,sprintf('%s/%s',figdir_out,plot_filename),'jpeg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- END OF CODE --------------