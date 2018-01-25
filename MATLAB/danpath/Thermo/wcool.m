%wcool.m
%Calculate reversible moist-adiabatic radiative-subsidence rate
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    ncdir_dr - the directory of the file you'd like
%    ncfile_in - the file you'd like
%    variable_in - the variable you'd like to extract
%    missing_value_flag - value that will be set to NaN
%
% Outputs:
%    data_out - matrix of the desired data
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
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 15 May 2014; Last revision:

%------------- BEGIN CODE --------------

clc
clear
close('all')

addpath(genpath('~/Dropbox/Research/MATLAB/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Inputs for moist adiabat
T_surf = 302.0;   %[K]; air temperature above surface
rh_surf = 80/100;   %[-] fraction
p_surf = 1015.00*100;   %[Pa]
dp = 50*100;  %[Pa]  
p = (1000:-dp/100:50)*100;    %[Pa]

%%Inputs for radiative subsidence rate calculation
Qcool_day = 1;  %[K/day] radiative cooling rate
pmax_wcool = 900*100;   %[Pa]
pmin_wcool = 700*100;   %[Pa]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%If needed, load constants from CM1
constants_CM1_createdatfile();  %creates constants_CM1_list
load constants_CM1_list


[r_surf] = r_rhinput(p_surf,T_surf,rh_surf); %[kg/kg]
[th_surf] = pottemp(p_surf,T_surf,r_surf);  %[K]
[Tr,r,rl,ri,Tp] = calculate_adiabat_DAN(T_surf,r_surf,p_surf,p);
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

[th] = pottemp(p,Tp,r);  %[K]
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
%set(gca,'position',[0.11    0.08    0.85    0.9]);

%%Default options -- as desired
set(0,'defaultaxesfontsize',18,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',3,'DefaultAxesFontName','Helvetica')

subplot(1,2,1)
h_pl = plot(th,p/100,'b',Tp,p/100,'b--');
set(gca, 'Ydir', 'reverse')
xlabel('Potential temperature / Temperature [K]')
ylabel('Pressure [hPa]')
legend('Pot temp','Temp')
box off
hold on
subplot(1,2,2)
h_pl = plot(r*1000,p/100,'g',rl*1000,p/100,'g--',ri*1000,p/100,'g:');
set(gca, 'Ydir', 'reverse')
xlabel('mixing ratios [g/kg]')
ylabel('Pressure [hPa]')
legend('Vapor','Liquid','Ice')
box off

%% Save plot %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('sounding_Tsfc%3.1fK_RHsfc%2.0f_Qcool%2.1fKday.pdf',T_surf,rh_surf,Qcool_day)
saveas(gcf,plot_filename,'pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Estimate dz from hydrostatic balance
dz = dp./(((rho(2:end)+rho(1:end-1))/2).*g);
dthdz = (th(2:end)-th(1:end-1))./dz

%%Estimate radiative subsidence rate
Qcool = Qcool_day/86400;    %[K/s]
wcool_p = Qcool./dthdz;    %[m/s]
dthdz_lowertrop = mean(dthdz(p<=pmax_wcool & p>=pmin_wcool))   %[K/m]
wcool_lowertrop = Qcool/dthdz_lowertrop

hh=figure(2);
clf(hh)

%%Position/size
set(hh,'Units','centimeters');
hpos = [0 0 15 15];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.17    0.12    0.77    0.86]);

%%Default options -- as desired
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',3,'DefaultAxesFontName','Helvetica')
p_mid = (p(2:end)+p(1:end-1))/2;
h_pl = plot(wcool_p,p_mid/100,'b');
set(gca, 'Ydir', 'reverse')
xlabel('radiative subsidence rate [m/s]')
ylabel('Pressure [hPa]')
%legend('Pot temp','Temp')
input_text = sprintf('w_{cool} = %3.2f cm/s',100*wcool_lowertrop)
text_location_xx = .0075;
text_location_yy = .5*max(p_mid)/100;
text(text_location_xx,text_location_yy,input_text,'FontSize',18)
axis([0 .02 0 1000])
box off


%% Save plot %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('wcool_Tsfc%3.1fK_RHsfc%2.0f_Qcool%2.1fKday.pdf',T_surf,rh_surf,Qcool_day)
saveas(gcf,plot_filename,'pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- END OF CODE --------------