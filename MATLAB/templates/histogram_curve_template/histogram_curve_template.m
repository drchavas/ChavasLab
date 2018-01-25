%histogram_curve_template.m -- 
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
set(0,'defaultaxesfontsize',18,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',3,'DefaultAxesFontName','Helvetica')

            
%% Data
xx = (randn(1000,1)).^2;
dx = 0.2;   %bin width
bin_min = 0;
bin_max = ceil(max(xx));

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
set(gca,'position',[0.14    0.09    0.84    0.86]);

% subplot = @(m,n,p) subtightplot(m, n, p, [0.15 0.15], [0.12 0.03], [0.1 0.03]);
    %%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)

edges = bin_min:dx:bin_max;
[nn] = histcounts(xx,edges,'Normalization','probability');
N_total = length(~isnan(xx));
centers = .2*(edges(1:end-1)+edges(2:end));

plot(centers,nn,'k')
hold on
scatter(centers,nn,30,'r')
xlabel('data x [units]')
ylabel('relative frequency [-]')
input_title = sprintf('relative frequency of x [units] (N=%i), dx=%3.2f',N_total,dx);
title(input_title,'FontSize',12)
    
    
%% Save plot %%%%%%%%%%%%%%%%%%
figfile_out = sprintf('histogram_curve_template.pdf');
saveas(gcf,sprintf('%s',figfile_out),'pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%------------- END OF CODE --------------