%sin_fit_temp.m

clear
clc
close all

addpath(genpath('~/Dropbox/Research/MATLAB/'));
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',1,'DefaultAxesFontName','Helvetica')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdata = 1:12;   %months
ydata = 10^3*[1.0174 1.3080 2.1874 3.0870 3.6900 3.8400 3.8800 3.5200 ...
   2.8000 2.1200 1.4600 1.1248];   %aribtrary CAPE values from NARR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Sinusoidal fit
[xdata_fit,ydata_fit,cycleinfo,gof,f,fit1] = sinfit(xdata,ydata);


%% PLOT DATA + MODEL FIT
figure
plot(xdata,ydata,'r.-')
hold on
plot(xdata_fit,ydata_fit,'b.-')
legend('data','model')
plot(cycleinfo.x_maxval*ones(1,2),[0 5000],'g--')
plot([0 13],cycleinfo.amp_offset*ones(1,2),'k--')
plot(cycleinfo.x_amp0*ones(1,2),[0 5000],'g:')
axis([0 13 0 5000])
xlabel('month')
ylabel('CAPE [J/kg]')
title('99ile CAPE from random year in NARR, northern louisiana')