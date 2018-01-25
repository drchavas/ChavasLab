%thermo.m

%Purpose: calculate some relevant thermodynamic quantities

%Created: 20 May 2011, Dan Chavas


function [T,rho,rh,the,s,s_sat,gam_m,thv] = thermo(p,th,r,rl)


%% Constants (values taken from CM1 model)
% c_CM1 = constants_CM1(); %c_CM1: [g rd cp rv p000 xxlv cpv]
% 
% g=c_CM1(1); %[m/s2]
% rd=c_CM1(2);  %[J/kg/K]
% cp=c_CM1(3); %[J/kg/K]; spec heat of dry air
% rv=c_CM1(4);   %[J/K/kg]
% p00 = c_CM1(5); %[Pa]
% xlv=c_CM1(6);   %[J/kg]
% cpv=c_CM1(7); %[J/kg/K]; spec heat of water
% 
% epsil=rd/rv;

%%If needed, load constants from CM1
constants_CM1_createdatfile();  %creates constants_CM1_list
load constants_CM1_list


%{
%constants
rd = 287;   %[J/K/kg]
cp = 1004; %[J/K/kg]
cpv = 1952; %[J/K/kg]
rv = 461;   %[J/K/kg]
epsil = rd/rv;
p00 = 100000;    %[Pa]
xlv = 2.5e6;    %[J/kg]
g = 9.8067; %[m/s^2]
%}

%%calculate actual temperature from potential temperature
T = temp_frompottemp(p,th,r);

%%calculate virtual temperature from actual temperature and wv mix ratio
Tv = T.*(1+.608*r);

%%calculate density from pressure and virtual temperature
rho = p./(rd*Tv);

%%calculate saturation vapor pressures
T_C = T-273.15;
es = 611.2*exp((17.67*T_C)./(T_C+243.5));  %[Pa]

%%calculate saturation mixing ratios
rs = epsil*(es./(p-es));  %[kg/kg]

%%calculate relative humidity
rh = r./rs;

%%calculate equivalent potential temperature
Cp = cp + r*cpv;
the=th.*exp((xlv.*r)./(Cp.*T));
%the2(i,j)=thet(i,j).*exp(LO.*q(i,j)./1000/(CP*(A*log(q(i,j)./1000)+B))); (from Brian Tang's ASPECH model)

%%calculate virtual potential temperature (BR09b, p.3)
thv=th.*(1+(rv/rd)*r)./(1+r+rl);

%%calculate entropy (eqn (1) Emanuel11)
s = cp*log(T) - rd*log(p) + xlv*r./T;

%%calculate saturation entropy (eqn (1) Emanuel11)
s_sat = cp*log(T) - rd*log(p) + xlv*rs./T;

%%calculate moist adiabatic lapse rate (see eqn 4.7.3 Emanuel Convection
%%book): ignoring liquid water content contribution
gam_d = (g/cp)*(1+r)./(1+r*(cpv/cp));
gam_m = gam_d.*(1+(xlv*r)./(rd.*T))./(1+xlv^2*r.*(1+r/epsil)./(rv*T.^2.*(cp+r*cpv)));

end