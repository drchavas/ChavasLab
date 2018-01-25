%bindata_template.m -- template for plotting a 2D binned dataset
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
% 2017-01-20; Last revision:

%------------- BEGIN CODE --------------

clear
clc
close all

addpath(genpath('~/Dropbox/Research/MATLAB/'));
set(0,'defaultaxesfontsize',14,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',1,'DefaultAxesFontName','Helvetica')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 2*randn(5000,1);
x2 = 3*randn(5000,1);
%x1 = [0 1 2 3] + 0.05;
%x2 = [0 0 0 0] + 0.05;
zz = x1 + x2;
%zz = ones(size(xx.^2));    %pure histogram
x1rg = -5:1:5;
x2rg = -7:1:7;
statistic = 'median';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Bin the data
[zm,x1b,x2b,x1m,x2m,nn] = bindata(zz,x1,x2,x1rg,x2rg,statistic);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot pure histogram
hh=figure(1);set(gcf,'Visible','on')
clf(hh)
set(hh,'Units','centimeters');
hpos = [0 0 30 30];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.08    0.05    0.83   0.9]);

hnum = 0;

%subplot(1,2,1);plot3(x(:,1),x(:,2),y,'.');
%subplot(1,2,2);h = imagesc(xrg,xrg,ym);
xdata = x1b(1,:);  %label for x-axis
ydata = x2b(:,1);   %label for y-axis
nn(nn==0) = NaN;
zdata = nn';    %pure histogram; must rotate data for imagesc
%zdata = zm';    %must rotate data for imagesc
h_pl = imagesc(xdata,ydata,zdata);
    %imagesc turns the matrix directly into a color map, rows to rows and
    %columns to columns.
set(gca,'YDir','normal')
set(h_pl,'alphadata',~isnan(zdata)) %only colors non-NaN values (found on google)

%caxis(clims);  %must set this before redefining colormap!
hclbr = colorbar;
clrbar_pos=get(hclbr,'Position');
clrbar_pos(3)=0.01;
clrbar_pos(1)=clrbar_pos(1)+.05;
clrbar_pos(2)=clrbar_pos(2)-.001;
set(hclbr,'Position',clrbar_pos,'FontSize',18)

hold on

%% Raw data on top -- to check that things are working here
hnum = hnum + 1;
hpl(hnum) = scatter(x1,x2,[],'.','MarkerEdgeColor','k') %just show black dots
input_legend{hnum} = 'raw';
hnum = hnum + 1;
hpl(hnum) = scatter(x1m(:),x2m(:),[],'.','MarkerEdgeColor','m') %(x1,x2) mean within each bin
input_legend{hnum} = '(x,y) bin mean';
plot([-100 100],[0 0],'k') %TEST to make sure this pops up where I want it to
plot([0 0],[-100 100],'k') %TEST to make sure this pops up where I want it to
hnum = hnum + 1;
hpl(hnum) = scatter(4,1,40,'r*') %TEST to make sure this pops up where I want it to
input_legend{hnum} = '(4,1) sanity check';

axis([-10 10 -10 10])
input_xlabel = 'xdat [-]';
xlabel(input_xlabel)
input_ylabel = 'ydat [-]';
ylabel(input_ylabel)
input_title = sprintf('bin statistic (color): count; x bins: [%3.2f,%3.2f), y bins: [%3.2f,%3.2f)',x1rg(1),x1rg(end),x2rg(1),x2rg(end));
title(input_title)
legend(hpl,input_legend)

%% Save plot
plot_filename = sprintf('bindata_template_hist.pdf');
saveas(gcf,plot_filename,'pdf')


%% Plot pure histogram
hh=figure(2);set(gcf,'Visible','on')
clf(hh)
set(hh,'Units','centimeters');
hpos = [0 0 30 30];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.08    0.06    0.83   0.89]);

hnum = 0;

%subplot(1,2,1);plot3(x(:,1),x(:,2),y,'.');
%subplot(1,2,2);h = imagesc(xrg,xrg,ym);
xdata = x1b(1,:);  %label for x-axis
ydata = x2b(:,1);   %label for y-axis
nn(nn==0) = NaN;
%zdata = nn';    %pure histogram; must rotate data for imagesc
zdata = zm';    %must rotate data for imagesc
h_pl = imagesc(xdata,ydata,zdata);
    %imagesc turns the matrix directly into a color map, rows to rows and
    %columns to columns.
set(gca,'YDir','normal')
set(h_pl,'alphadata',~isnan(zdata)) %only colors non-NaN values (found on google)

%caxis(clims);  %must set this before redefining colormap!
hclbr = colorbar;
clrbar_pos=get(hclbr,'Position');
clrbar_pos(3)=0.01;
clrbar_pos(1)=clrbar_pos(1)+.05;
clrbar_pos(2)=clrbar_pos(2)-.001;
set(hclbr,'Position',clrbar_pos,'FontSize',18)

hold on

%% Raw data on top -- to check that things are working here
hnum = hnum + 1;
hpl(hnum) = scatter(x1,x2,[],'.','MarkerEdgeColor','k') %just show black dots
input_legend{hnum} = 'raw';
hnum = hnum + 1;
hpl(hnum) = scatter(x1m(:),x2m(:),[],'.','MarkerEdgeColor','m') %(x1,x2) mean within each bin
input_legend{hnum} = '(x,y) bin mean';
plot([-100 100],[0 0],'k') %TEST to make sure this pops up where I want it to
plot([0 0],[-100 100],'k') %TEST to make sure this pops up where I want it to
hnum = hnum + 1;
hpl(hnum) = scatter(4,1,40,'r*') %TEST to make sure this pops up where I want it to
input_legend{hnum} = '(4,1) sanity check';

axis([-10 10 -10 10])
input_xlabel = 'xdat [-]';
xlabel(input_xlabel)
input_ylabel = 'ydat [-]';
ylabel(input_ylabel)
input_title = sprintf('bin statistic (color): %s; x bins: [%3.2f,%3.2f), y bins: [%3.2f,%3.2f)',statistic,x1rg(1),x1rg(end),x2rg(1),x2rg(end));
title(input_title)
legend(hpl,input_legend)

%% Save plot %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('bindata_template.pdf');
saveas(gcf,plot_filename,'pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%