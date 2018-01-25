%adiabat_plot.m
%Plot a moist adiabat, both reversible and pseudo-adiabatic

clc
clear
close('all')

addpath(genpath('~/Dropbox/Research/MATLAB/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
T_surf = 300.0;   %[K]
rh_surf = 80/100;   %[-] fraction
p_surf = 1015.00*100;   %[Pa]
dp = 50*100;  %[Pa]  
p = (1000:-dp/100:50)*100;    %[Pa]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%If needed, load constants from CM1
constants_CM1_createdatfile();  %creates constants_CM1_list
load constants_CM1_list


[r_surf] = r_rhinput(p_surf,T_surf,rh_surf); %[kg/kg]
[th_surf] = pottemp(p_surf,T_surf,r_surf);  %[K]
[Tr,r,rl,ri,Tp] = adiabat(T_surf,r_surf,p_surf,p);
% Calculates a reversible and psuedo-adiabat from an input 
% Temperature (K)
% mixing ratio (r_surf)
% pressure     (p_surf)
% 
% The adiabat Temperature and mixing ratios are given at pressures input (p)
%
% The reversible moist adiabat includes ice consistent with its treatment
% in CM1. The pseudo-adiabat is water only - it should give the same result
% as Kerry's subroutine.
%

[th] = pottemp(p,Tr,r);  %[K]
%[Ttest] = temp_frompottemp(p,th,r);  %[K]
[T_temp,rho,rh,the,s,s_sat,gam_m,thv] = thermo(p,th,r,rl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial setup
hh=figure(1);
clf(hh)

%%Position/size
set(hh,'Units','centimeters');
hpos = [0 0 25 25];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.12    0.1    0.84    0.81]);

%%Default options -- as desired
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',3,'DefaultAxesFontName','Helvetica')

subplot(1,2,1)
h_pl = plot(th,p/100,'b',Tr,p/100,'b--');
set(gca, 'Ydir', 'reverse')
xlabel('Potential temperature / Temperature [K]')
ylabel('Pressure [hPa]')
legend('Pot temp','Temp')
hold on
subplot(1,2,2)
h_pl = plot(r*1000,p/100,'g',rl*1000,p/100,'g--',ri*1000,p/100,'g:');
set(gca, 'Ydir', 'reverse')
xlabel('mixing ratios [g/kg]')
ylabel('Pressure [hPa]')
legend('Vapor','Liquid','Ice')

%% Save plot %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('adiabat.pdf')
saveas(gcf,plot_filename,'pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
