%mfile_template.m -- 
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

%%Parameters

%%Plotting
figdir_out = sprintf('.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
%{
datapath_in = sprintf('%s/%s',datadir_in,datafile_in);
listOfVariables = {
    };
load(datapath_in, listOfVariables{:});
sprintf('Loading data from %s',datapath_in)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure(1);%set(gcf,'Visible','off')
clf(fig1)
set(fig1,'Units','centimeters');
%hpos = [0 0 30 30];
hpos = [0 0 20 20];
%hpos = [0 0 60 30];
set(fig1,'Position',hpos);
set(fig1,'PaperUnits','centimeters');
set(fig1,'PaperPosition',hpos);
set(fig1,'PaperSize',hpos(3:4));

set(gca,'position',[0.09    0.09    0.86    0.86]);
subplot = @(m,n,p) subtightplot(m, n, p, [0.15 0.15], [0.12 0.03], [0.1 0.03]);
    %%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)

subplot(2,2,1)
plot(0:4,0:4)
xlabel('x-axis')
ylabel('y-axis')
%title('TEST')

subplot(2,2,2)
plot(0:4,0:4)
xlabel('x-axis')
ylabel('y-axis')
%title('TEST')

subplot(2,2,3)
plot(0:4,0:4)
xlabel('x-axis')
ylabel('y-axis')
%title('TEST')

subplot(2,2,4)
plot(0:4,0:4)
xlabel('x-axis')
ylabel('y-axis')
%title('TEST')

%% Save plot %%%%%%%%%%%%%%%%%%
%figfile_out = sprintf('TEST.pdf');
%saveas(gcf,sprintf('%s/%s',figdir_out,figfile_out),'pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%------------- END OF CODE --------------