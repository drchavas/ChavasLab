%bindata2_example.m

clear
clc
close all

x = randn(500,2);
y = sum(x.^2,2) + randn(500,1);
xrg = linspace(-3,3,10)';
[ym,yb] = bindata2(y,x(:,1),x(:,2),xrg,xrg);

hh=figure(1);set(gcf,'Visible','on')
clf(hh)
set(hh,'Units','centimeters');
hpos = [0 0 30 30];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.08    0.05    0.87   0.9]);

subplot(1,2,1);plot3(x(:,1),x(:,2),y,'.');
subplot(1,2,2);h = imagesc(xrg,xrg,ym);
set(h,'AlphaData',~isnan(ym)); box off;

%% Save plot
plot_filename = sprintf('bindata2_example.pdf');
saveas(gcf,plot_filename,'pdf')