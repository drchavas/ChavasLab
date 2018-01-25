%GWB_V_to_P.m -- calculate pressure field for given wind and density profile
%
% Syntax:
%
% Inputs:
%
% Outputs:
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 24 Apr 2015; Last revision:

% Revisions:
%------------- BEGIN CODE --------------

%function [Pmin,PP] = GWB_V_to_P(rr,VV,rho,renv,Penv,Lat_BT)
function [dPmax,dP] = GWB_V_to_P(rr,VV,rho,renv,Lat_BT)

assert(max(rr(~isnan(VV)))>=renv,'renv lies beyond the input radial profile')

%% Constants
omeg = 7.292e-5;    %[s-1]; Earth's rotation rate

%% Calculate some things
fcor = 2*omeg*sind(abs(Lat_BT));

%% Define a high-res radial vector to perform the actual calculation
rr_temp = (0:1000:renv);    %[m]
VV_temp = interp1(rr,VV,rr_temp,'pchip',NaN);
rho_temp = interp1(rr,rho,rr_temp,'pchip',NaN);

%% Calculate Pmin: gradient wind balance
term1 = fcor*VV_temp;
term2 = (VV_temp.^2)./rr_temp;
dPdr = rho_temp.*(term1 + term2);
dr = mean(rr_temp(2:end)-rr_temp(1:end-1));
dP_temp = flip(cumsum(flip(dPdr*dr,2)),2);   %[Pa]; pressure deficit
dP_temp(1) = dP_temp(2)+(dP_temp(2)-dP_temp(3));    %linear extrap to r=0;
%PP_temp = Penv-dP_temp;
%Pmin = min(PP_temp);

%% TESTING
%{
figure(999)
clf(999)
subplot(2,1,1)
plot(rr_temp/1000,VV_temp,'k')
hold on
subplot(2,1,2)
plot(rr_temp/1000,PP_temp/100,'k')
hold on
%}

%% Interpolate to user grid
%PP = interp1(rr_temp,PP_temp,rr,'pchip',NaN);
dP = interp1(rr_temp,dP_temp,rr,'pchip',NaN);
dPmax = max(dP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTING:  V and P
%{
hh=figure(1);
clf(hh)
set(hh,'Units','centimeters');
hpos = [0 0 30 30];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.08    0.09    0.90    0.89]);

clear hpl input_legend
hnum = 0;

subplot(2,1,1)
plot(rr/1000,pp1_grid,'b')
hold on
plot(rr_p2_grid/1000,pp2_grid,'r')
xlabel('radius [km]')
ylabel('pressure [hPa]')

subplot(2,1,2)
plot(rr_v1/1000,VV1,'b')
hold on
%plot(rr_v1/1000,VV1_cycl,'b:')
plot(rr_v2/1000,VV2,'r')
%plot(rr_v2/1000,VV2_cycl,'r:')
xlabel('radius [km]')
ylabel('Wind speed [m/s]')
axis([0 250 0 70])
%}

end

%------------- END OF CODE --------------
