function [vs] = windprofile_E04(lat,vms,rmw,ro,r)
% v = Symetrical wind (m/s)
% lat = storm latitude
% long = storm longitude   
% pcs = storm centrol pressure (mb)
% vms = maximum storm symetrical wind speed at 10 m (m/s)
% rms = radius of maximum wind (km)
% ro = outer radius (km)
% r = radius of the point of interest (or a number of points; km);

% Vortex radial profile  (Emanule et al. 2006)

% Shape parameters  
b=0.25;
nb=0.9;
mb=1.6;
%
fac1=(1-b).*(mb+nb);   %1.875
fac2=b.*(1+2.*mb);     %1.05
fac3=2.*(mb+nb);       %5
fac4=2.*mb+1;          %4.2
%

rat=r./max(rmw,1);
vs=vms.*((ro-r)./(ro-rmw)).*(rat).^mb.*sqrt(fac1./(nb+mb.*rat.^fac3)+fac2./(1+2*mb.*rat.^fac4));
vs=vs.*lat./(abs(lat)+1e-8);

% vs<0 when r>ro
vs(vs<0)=0;


