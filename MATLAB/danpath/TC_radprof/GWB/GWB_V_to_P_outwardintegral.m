%GWB_V_to_P_outward integral.m -- calculate dP from gradient wind balance
%integrated radially outwards using centered difference scheme
%
% Syntax:
%
% Inputs:
%
% Outputs:
%
% deltaP_r [Pa]: radial profile of integrated dP moving radially outwards
% up to r=5000 km
%   (deltaP_r(r=0) = 0)
% deltaP_radii [Pa]: deltaP interpolated to input radii vector
% deltaPcent [Pa]: deltaPcent = deltaP_radii(end) (i.e. at largest input radius)
% deltaPmax [Pa}: max(deltaP_r) for entire vector
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
% 22 Mar 2017; Last revision:

% Revisions:
%------------- BEGIN CODE --------------

function [deltaP_r,deltaP_radii,deltaPcent,deltaPmax] = GWB_V_to_P_outwardintegral(rr,VV,rho,r_radii_in,fcor_BT)

assert(max(rr(~isnan(VV)))>=max(r_radii_in),'outermost input radius lies beyond the input radial profile')

%% Calculate some things
fcor = abs(fcor_BT);

%% Define a high-res radial vector to perform the actual calculation
dr = 2*1000;
rr_temp = (0:dr:(max(rr)+dr));    %[m]
rr_temp_midpt = (rr_temp(2:end)+rr_temp(1:end-1))/2;    %[m]
ii_good = ~isnan(VV);
VV_temp_midpt = interp1(rr(ii_good),VV(ii_good),rr_temp_midpt,'pchip',NaN);
ii_good = ~isnan(rho);
rho_temp_midpt = interp1(rr(ii_good),rho(ii_good),rr_temp_midpt,'pchip',NaN);

%% Calculate Pmin: gradient wind balance -- centered finite difference
term1 = fcor*VV_temp_midpt;
term2 = (VV_temp_midpt.^2)./rr_temp_midpt;
dPdr_midpt = rho_temp_midpt.*(term1 + term2);
dr = rr_temp(2:end)-rr_temp(1:end-1);
dP_r_temp = dPdr_midpt.*dr;

%Integrate pressure changes radially outwards from center
deltaP_r_temp = cumsum([0 dP_r_temp]);

%%% TESTING %%%%%%%%%%%%%%%
%{
figure(999)
clf(999)
plot(rr_temp/1000,deltaP_r/100)
'hi'
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dP_temp = flip(cumsum(flip(dPdr.*dr,2)),2);   %[Pa]; pressure deficit
% dP_temp(1) = dP_temp(2)+(dP_temp(2)-dP_temp(3));    %linear extrap to r=0;
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
ii_good = ~isnan(deltaP_r_temp);
deltaP_r = interp1(rr_temp(ii_good),deltaP_r_temp(ii_good),rr,'pchip',NaN);
ii_good = ~isnan(deltaP_r);

%% If necessary, extrapolate out just beyond outermost wind radius (or else you'll get a NaN in deltaP_radii)
rr_temp = rr(ii_good);
deltaP_r_temp = deltaP_r(ii_good);
dr_extra = 1.1*(max(r_radii_in)-max(rr_temp));
if(dr_extra>0)  %outermost radii is beyond interpolated vector (happens rarely)
    rr_temp(end+1) = rr_temp(end) + dr_extra;
    deltaP_r_temp(end+1) = deltaP_r_temp(end) + dr_extra*(deltaP_r_temp(end)-deltaP_r_temp(end-1))/(rr_temp(end)-rr_temp(end-1));
end

%% Interpolate to radii
deltaP_radii = interp1(rr_temp,deltaP_r_temp,r_radii_in,'pchip',NaN);
deltaPcent = deltaP_radii(r_radii_in==max(r_radii_in)); %outermost valid radius
deltaPmax = max(deltaP_r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTING:  V and P
%{
hh=figure(999);
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
plot(rr/1000,deltaP_r/100,'b')
hold on
plot(r_radii_in/1000,deltaP_radii/100,'go')
plot(max(r_radii_in)/1000,deltaPcent/100,'rx')
xlabel('radius [km]')
ylabel('pressure [hPa]')
title(sprintf('max(deltaP) = %2.1f hPa; deltaPcent = %2.1f hPa',deltaPmax/100,deltaPcent/100))
axis([0 5000 0 1.1*deltaPmax/100])

subplot(2,1,2)
plot(rr/1000,VV,'b')
hold on
plot(r_radii_in/1000,0,'go')
plot(max(r_radii_in)/1000,0,'rx')
plot([-10^9 10^9],[0 0],'k-')
xlabel('radius [km]')
ylabel('Wind speed [m/s]')
axis([0 5000 -10 1.1*max(VV)])
%}

end

%------------- END OF CODE --------------
