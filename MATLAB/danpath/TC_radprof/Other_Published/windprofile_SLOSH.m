function [vs] = windprofile_SLOSH(vms,rmw,r)
% vs = Symetrical wind (m/s)
% long = storm longitude   
% pcs = storm centrol pressure (mb)
% vms = maximum storm symetrical wind speed at 10 m (m/s)
% rms = radius of maximum wind (km)
% ro = outer radius (km)
% r = radius of the point of interest (or a number of points; km);

% lat=18.6;
% vms=166;
% rmw=16.76;
% r=1:1:600;

rmw=rmw*1000;
r=r.*1000;

% SLOSH wind profile
vs=vms*2*rmw.*r./(rmw^2+r.^2);

vs(vs<0)=0;

