%bindata_example.m

%x = randn(100,1);
x = (0:.1:10)+.05;
%y = x.^2 + randn(100,1);
y = -x;
y = y + randn(size(y))/2;
%xrg = linspace(-.5,.5,20)';
xrg = (-2:1:12)';   %bin edges
[ym,yb,xm] = bindata(y,x,xrg);



hh=figure(1);set(gcf,'Visible','on')
clf(hh)
set(hh,'Units','centimeters');
hpos = [0 0 30 30];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.08    0.05    0.87   0.9]);

X = [xrg(1:end-1),xrg(2:end)]';
Y = [ym,ym]';
plot(x,y,'.',X(:),Y(:),'r-');
hold on
plot(xm,ym,'m');
legend('raw data','binned center(x)-mean(y)','binned mean(x)-mean(y)')

%% Save plot
plot_filename = sprintf('bindata_example.pdf');
saveas(gcf,plot_filename,'pdf')