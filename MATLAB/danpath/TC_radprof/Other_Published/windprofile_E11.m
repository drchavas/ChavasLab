function [vs] = windprofile_E11(lat,vms,rmw,r,C)
% vs = Symetrical wind (m/s)
% lat = storm latitude
% long = storm longitude   
% pcs = storm centrol pressure (mb)
% vms = maximum storm symetrical wind speed at 10 m (m/s)
% rms = radius of maximum wind (km)
% ro = outer radius (km)
% r = radius of the point of interest (or a number of points; km);
% C=Ck/Cd; C~1 in general, but not for intense storms, in which C~0.5

% lat=18.6;
% vms=166;
% rmw=16.76;
% r=1:1:600;

rmw=rmw*1000;
r=r.*1000;

omega=7.292e-5; %[rad s-1]; rotation rate of the Earth
f=abs(2*omega*sind(lat)); %[s-1]; Coriolis parameter

% Vortex radial profile in terms of angular momentum  (Emanuel and Rotunno 2011)
%
Mm=rmw*vms+0.5*f*rmw^2;
rat=r./rmw;
%
% Parameters  Ck/Cd=1;
% M=Mm.*2.*(rat.^2)./(1.+rat.^2);
% vs=(M-0.5.*f.*r.^2)./r;

M=Mm.*(2.*(rat.^2)./(2-C+C.*rat.^2)).^(1/(2-C));
vs=(M-0.5.*f.*r.^2)./r;

% vs<0 when r>ro
vs(vs<0)=0;


