% calculate surface (10-m) symetrical wind speed and pressure profile 
% with Holland_80 model
% by Ning Lin Mar. 2011

function [vs, ps] = windpressure_H80(lat,pcs,vms,rmw,r)
% vs=symetrical wind (m/s)
% ps=pressure (mb)
% lat = storm latitude
% long = storm longitude   
% pcs = storm centrol pressure (mb)
% vms = maximum storm symetrical wind speed at 10 m (m/s)
% rms = radius of maximum wind (km)
% ro = outer radius (km)
% r = radius of the point of interest (or a number of points; km);

rmw=rmw*1000;
r=r.*1000;

omega=7.292e-5; %[rad s-1]; rotation rate of the Earth
f=abs(2*omega*sind(lat)); %[s-1]; Coriolis parameter

pcs=pcs*100;  % convert mb to pa

% standard pressure (pa) (may use data)
pns=101325;

% Air density
Rho=1.15;

%%
% Holland80 approximate B parameter used in ADCIRC
% B=vms^2*exp(1)*Rho/(pns-pcs);

% % Holland80 B exact
 vv=(vms+f*rmw/2)^2-f^2*rmw^2/4;
 B=vv*exp(1)*Rho/(pns-pcs);

% calculate pressure 
ps=pcs+(pns-pcs).*exp(-(rmw./r).^B);
ps=ps/100;  %change to mb

% % the cyclostrophic wind profile (need to use the approximate B)
% vs=sqrt((rmw./r).^B.*B.*(pns-pcs)./Rho.*exp(-(rmw./r).^B));

% or gradient wind profile (need to use the exact B)
vs=sqrt((rmw./r).^B.*B.*(pns-pcs)./Rho.*exp(-(rmw./r).^B)+r.^2*f^2/4)-f.*r/2;

vs(vs<0)=0;
