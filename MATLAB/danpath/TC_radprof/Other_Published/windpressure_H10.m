% calculate surface (10-m) symetrical wind speed and pressure at a point 
% with Holland_10 model
% by Ning Lin Nov. 19, 2010

function [vs, ps] = windpressure_H10(lat,pcs,vms,rmw,ro,r,r12)
% vs=symetrical wind (m/s)
% lat = storm latitude
% long = storm longitude   
% pcs = storm centrol pressure (mb)
% vms = maximum storm symetrical wind speed at 10 m (m/s)
% rms = radius of maximum wind (km)
% ro = outer radius (km)
% r = radius of the point of interest (or a number of points; km);
% r12 = radius of v=12m/s;

% lat=18.6;
% pcs=886.8;
% vms=166*0.5144; 
% rmw=16.76;
% ro=400;
% 
% lat=40;
% long=285.93;
% pcs=960.38;
% vms=93**0.5144;
% rmw=33;
% ro=400;

% plat=40.7;
% plong=285;
% % % Distant between storm center and the point
% % % pifac=acos(-1)/180;
% % % dfac=60.*1.852;
% % % dx=dfac.*cos(pifac.*plat).*(plong-long);
% % % dy=dfac.*(plat-lat);
% % % r=sqrt(dx.^2+dy.^2);
% or used a vector of r
% r=1:1:600;

N=length(r);
vs=zeros(N,1);
vso=zeros(N,1);
x=zeros(N,1);
xo=zeros(N,1);

pcs=pcs*100;  % convert mb to pa

% standard pressure (pa) (may use data)
pns=101325;

% pressure at the rmw
pms=pcs+(pns-pcs)*exp(-1);

% SST at rmw (may use data)
Ts=28-3*(lat-10)/20;   
% vapor pressure
qm=0.9*3.802/pms*exp(17.67*Ts/(243.5+Ts));
% virtual surface air temperature
Tvs=(Ts+273.15)*(1+0.61*qm);

% Holland10 bs parameter
Rho=pms/(286.9*Tvs);  % ideal gas law
bs=vms^2*exp(1)*Rho/(pns-pcs);

% calculate pressure 
ps=pcs+(pns-pcs).*exp(-(rmw./r).^bs);

% calculate the exponent of outer wind profile, using ro
xro=log(0.01)/log((rmw/ro)^bs*exp(1-(rmw/ro)^bs));
vro=vms*((rmw/ro)^bs*exp(1-(rmw/ro)^bs))^xro;   %check vro should be 0.01*vms

% % calculate the exponent of outer wind profile, using r12, approximate eq.
% % Constants
% w_rad=.016; %[ms-1]; clear-sky free tropospheric subsidence velocity in tropics
% omega=7.292e-5; %[rad s-1]; rotation rate of the Earth
% f=abs(2*omega*sind(lat)); %[s-1]; Coriolis parameter
% C_D=1e-3;   %drag coefficient
% v12=12; %m/s
% %from Emanuel model
% a=1;b=2*C_D*v12^2/f/w_rad;c=-(ro*1000)^2;
% r12=(-b+sqrt(b^2-4*a*c))/2/a;
% v12c=sqrt(f*w_rad/(2*C_D)*((ro*1000)^2-r12^2)/r12);    %check v12c=v12
% r12=r12/1000; %km

% calculate r12 from Dan's numerical method
v12=12; %m/s
%[r12,residual,niter] = Lilly_LIN(ro,12,lat,0);
 
x12=log(v12/vms)/log((rmw/r12)^bs*exp(1-(rmw/r12)^bs));
vr12=vms*((rmw/r12)^bs*exp(1-(rmw/r12)^bs))^x12;   %check vr12=v12
vrro=vms*((rmw/ro)^bs*exp(1-(rmw/ro)^bs))^x12;   %check vro 

% if x12<0.5           % avoid possible blow up at large radius
%     x12=0.5;
% end

% calculate surface wind
for n=1:N
    if r(n)<=rmw
        x(n)=0.5;
    else
        x(n)=0.5+(r(n)-rmw)*(x12-0.5)/(r12-rmw);
    end
    vs(n)=vms*((rmw/r(n))^bs*exp(1-(rmw/r(n))^bs))^x(n);
end

%if use xo determined from ro
for n=1:N
    if r(n)<=rmw
        xo(n)=0.5;
    else
        xo(n)=0.5+(r(n)-rmw)*(xro-0.5)/(ro-rmw);
    end
    vso(n)=vms*((rmw/r(n))^bs*exp(1-(rmw/r(n))^bs))^xo(n);
end

ps=ps/100;  % convert pa back to mb

%vs=vs/0.5144; % convert from m/s back to kts
%vso=vso/0.5144;

%check if use x=0.5
%vs05=vms.*((rmw./r).^bs.*exp(1-(rmw./r).^bs)).^0.5; 
%vs05=vs05/0.5144;

