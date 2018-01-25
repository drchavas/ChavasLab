%%plot_basic_template.m -- a template for plotting basic things


%------------- BEGIN CODE --------------

clear
clc
close('all') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx = (1:10)';
yy = xx.^2-2*randn(10,1);
zz = xx.^.5-randn(10,1);
xx2 = 10*rand(100,1);
zz2 = 10*randn(100,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial setup
hh=figure(1);
clf(hh)

%%Position/size
set(hh,'Units','centimeters');
hpos = [0 0 15 15];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.12    0.1    0.84    0.81]);

%%Default options -- as desired
set(0,'defaultaxesfontsize',22,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',3,'DefaultAxesFontName','Helvetica')

%% Plot the data
hold off

%%Begin legend iteration
hnum = 1;
clear input_legend

%%Plot yy
color_in = [.5 0 0];
hh(hnum) = plot(xx,yy);
set(hh(hnum),'Color',color_in,'LineWidth',3);
hold on
input_legend{hnum} = 'yy';
hnum = hnum + 1;

%%Plot zz
color_in = .5*color_in;
hh(hnum) = plot(xx,zz);
set(hh(hnum),'Color',color_in,'LineWidth',3);
hold on
input_legend{hnum} = 'zz';
hnum = hnum + 1;

%%Plot zz2
color_in = [.5 .5 .5];
hh(hnum) = scatter(xx2,zz2);
set(hh(hnum),'Marker','o','SizeData',40,'MarkerEdgeColor',color_in,'MarkerFaceColor','g');
hold on
input_legend{hnum} = 'zz2';
hnum = hnum + 1;


%% Plot aesthetics
%%Define plot boundaries
global_max = max([max(yy) max(zz) max(zz2)]);
global_min = min([min(yy) min(zz) min(zz2)]);
xx_plot_min = min([min(xx) min(xx2) 0]);
xx_plot_max = 1.0*max([max(xx) max(xx2)]);
yy_plot_min = min([global_min 0]);
yy_plot_max = 1.1*global_max;
axis([xx_plot_min xx_plot_max yy_plot_min yy_plot_max])

%%Add a zero line
xx_zero = [xx_plot_min xx_plot_max];
plot(xx_zero,[0 0],'k')

%%Annotate the plot with the global max value
text_location_xx = .5*xx_plot_max;
text_location_yy = .9*yy_plot_max;
text0=text(text_location_xx,text_location_yy,...
    sprintf('$Max = $ %3.2f $ms^{-1}$',global_max),...
    'fontweight','bold','Color','k');
set(text0,'HorizontalAlignment','center','VerticalAlignment',...
    'top','Interpreter','Latex','BackgroundColor','white');

%%Axis labels and title
xlabel('Distance [$m$]','Interpreter','LaTeX')
ylabel('Speed [$ms^{-1}$]','Interpreter','LaTeX')
input_title1 = sprintf('Example plot');
input_title2 = sprintf('Max value = %3.2f',global_max);
temp  =title({input_title1 input_title2},'Interpreter','Latex')

%%Legend
hleg = legend(hh,input_legend,'Interpreter','Latex','FontSize',12);
% set(hleg,'FontSize',14);

%%Other
grid on     %add grid lines
box on  %add box around plot (some output file formats dont do this)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_filename = sprintf('plot_basic_template.pdf')
saveas(gcf,plot_filename,'pdf')

%------------- END OF CODE --------------
