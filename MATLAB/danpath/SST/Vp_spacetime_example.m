%Vp_spacetime_example.m - test Vp_spacetime
%
% Other m-files required: Vp_spacetime, land_or_ocean
% Subfunctions: none
%
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: http://www.princeton.edu/~dchavas/
% 5 Sep 2014; Last revision:
% Revision history:

%------------- BEGIN CODE --------------

clear
clc
close all

addpath(genpath('~/Dropbox/Research/MATLAB/danpath/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Vp data
Vp_file = sprintf('/Volumes/Users/Dan Chavas/DATA/ERAI_Vp_19992009_momean.mat');

%%Point of interest
N_pts = 9;
lats_in = round(linspace(-40,40,N_pts));    %[deg N]
lons_in = 3*lats_in+100;   %[deg E]
% lats_in = 20;    %[deg N]
% lons_in = -50;   %[deg E]
years_in = round(linspace(1999,2009,N_pts));
daysofyear_in = round(linspace(40,320,N_pts));
% years_in = 1971*ones(size(lats_in));
% daysofyear_in = 227*ones(size(lats_in));
% months_in = 8*ones(size(lats_in));
% days_in = 1*ones(size(lats_in));

%%Make a plot?
make_plot = 1;  %1: makes a plot; ow: no plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Load Vp monthly data
listOfVariables={
        'days_since_197001010000_all',...
        'lat_mat','lon_mat','Vp_all'
        };
load(Vp_file,listOfVariables{:})
sprintf('Loading Vp data from %s',Vp_file)

length(lats_in)
[Vp_interp_all,isOcean_all] = data_spacetime(days_since_197001010000_all,lat_mat,lon_mat,Vp_all,lats_in,...
            lons_in,years_in,daysofyear_in);
clear lat_mat lon_mat days_since_197001010000_all Vp_all

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING: Plot points, values on top of closest monthly %%
if(make_plot)
    
    %% Load Vp monthly data
    listOfVariables={
            'year_all','month_all','monthnum','days_since_197001010000_all',...
            'days_since_197001010000_all',...
            'lat_mat','lon_mat','Vp_all'
            };
    load(Vp_file,listOfVariables{:})
    sprintf('Loading Vp data from %s',Vp_file)

    coastal_res = 1;    %[pts/deg]
    make_plot = 0;
    [isOcean_mat] = land_or_ocean(lat_mat,lon_mat,coastal_res,make_plot);

    for jj=1:length(Vp_interp_all)
    
        Vp_interp = Vp_interp_all(jj);
        isOcean = isOcean_all(jj);
        
        lat_in = lats_in(jj);
        lon_in = lons_in(jj);
        year_in = years_in(jj);
        day_in_year = daysofyear_in(jj);

        %% Calculate total days since 1 Jan 1970 0000UTC (Unix standard)
        year0_Unix = 1970;
        month0_Unix = 1;
        day0_Unix = 1;    %can use fractions of days
        days_since_197001010000_in = (day_in_year-1) + days_since(year_in,...
            ones(size(year_in)),ones(size(year_in)),year0_Unix,month0_Unix,day0_Unix);

        %%Extract Vp map data for desired time
        dtime = abs(days_since_197001010000_all-mode(days_since_197001010000_in));
        i_time_pl = find(dtime==min(dtime));
        i_time_pl = i_time_pl(1);
        Vp_temp = Vp_all(:,:,i_time_pl);
        days_since_197001010000_map = days_since_197001010000_all(i_time_pl);

        Vp_temp(~isOcean_mat) = NaN;

        %%Adjust lon to [0,360) deg E
        lon_in_temp = lon_in;
        lon_in_temp(lon_in_temp<0) = lon_in_temp(lon_in_temp<0)+360;

        %%INITIAL SETUP %%%%%%%%
        hh=figure(1001);
        clf(hh)
        set(hh,'units','centimeters');
        hpos = [0 0 60 30];
        set(hh,'Position',hpos);
        set(hh,'PaperUnits','centimeters');
        set(hh,'PaperPosition',hpos);
        set(hh,'PaperSize',hpos(3:4));

        set(gca,'position',[0.08    0.07    0.86    0.91]);

    %     subplot(2,1,1)
    %     ax = worldmap('World');
    %     setm(ax, 'Origin', [0 180 0])
    %     land = shaperead('landareas', 'UseGeoCoords', true);
    %     geoshow(ax, land, 'FaceColor', [0.8 0.9 0.8])
    %     hold on
    %     contourfm(lat_mat,lon_mat,sst_temp,0:1:max(sst_temp(:)));
    %     %caxis([0 25])
    %     colorbar

        %subplot(2,1,2)
        caxis([0 80]);  %must set this before redefining colormap!
        cvals = colormap(bluewhitered(256));   %custom color map: pos = red, neg = blue, 0 = white
        xvals = linspace(0,80,length(cvals)); %values corresponding to color map scale
        ax = worldmap([-90 90],[-180 180]);
        setm(ax, 'Origin', [0 180 0])
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(ax, land, 'FaceColor', [0.8 0.9 0.8])
        hold on
        cint = 10;   %[K]
        [cpl,hpl] = contourfm(lat_mat,lon_mat,Vp_temp,0:cint:80);
        ht = clabelm(cpl,hpl,'LabelSpacing',10^4);
        set(ht,'Color','r','BackgroundColor','white','FontWeight','bold')
        contourm(lat_mat,lon_mat,Vp_temp,[0 0],'Color','g','LineStyle','--','LineWidth',2);  %%Zero contour
        %caxis([0 25])
        colorbar
        Vp_colors = interp1(xvals,cvals,Vp_interp);
        if(isOcean==1)
            lat_pl_temp = lat_in(isOcean);
            lon_pl_temp = lon_in(isOcean);
            Vp_color_pl_temp = Vp_colors(isOcean,:);
            for ii=1:length(lat_pl_temp)
                if(sum(isnan(Vp_color_pl_temp(ii,:)))==0)
                    plotm(lat_pl_temp(ii),lon_pl_temp(ii),'Marker','o','MarkerEdgeColor','g','MarkerFaceColor',Vp_color_pl_temp(ii,:),'MarkerSize',20,'LineWidth',5)
                else
                    plotm(lat_pl_temp(ii),lon_pl_temp(ii),'Marker','o','MarkerEdgeColor','g','MarkerFaceColor','k','MarkerSize',20,'LineWidth',5)
                end
            end
        end
        if(isOcean==0)
            plotm(lat_in(~isOcean),lon_in_temp(~isOcean),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',20,'LineWidth',5)
        end
        title(sprintf('pt: Vp = %3.1f, yr = %i, days since 1970 = %i; map: days since 1970 = %i',Vp_interp,year_in,days_since_197001010000_in,days_since_197001010000_map))

        %% Save plot
        plot_filename = sprintf('Vp_spacetime_example_%i.jpg',jj)
        saveas(gcf,plot_filename,'jpg')
        
    end
    
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
%------------- END OF CODE --------------