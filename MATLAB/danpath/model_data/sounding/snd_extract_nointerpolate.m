%snd_extract.m

%Purpose: extract input_sounding data and linearly interpolate it to the
%desired model levels as is done in CM1. For the pressure field, this code
%outputs only the hydrostatically balanced pressure field associated with
%the distribution of theta, qv, and p_sfc.  (this is the background state CM1 calculates)

%Updated: 09 Feb 2012 DRC

%clear
%clc

function [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract_nointerpolate(dir_in,snd_file)



%% Constants (values taken from CM1 model)
c_CM1 = constants_CM1(); %c_CM1: [g rd cp rv p00 xlv cpv]

g=c_CM1(1); %[m/s2]
Rd=c_CM1(2);  %[J/kg/K]
Cpd=c_CM1(3); %[J/kg/K]; spec heat of dry air
Rv=c_CM1(4);   %[J/K/kg]
p0 = c_CM1(5); %[Pa]
Lv=c_CM1(6);   %[J/kg]
Cpv=c_CM1(7); %[J/kg/K]; spec heat of water

eps=Rd/Rv;

%% EXTRACT DATA FROM input_sounding
sounding_fil = sprintf('%s/%s',dir_in,snd_file);
clear fid
fid = fopen(sounding_fil);

%%header: %1) sfc pres (mb); 2) sfc theta (K); 3) sfc qv (g/kg)
temp=fgetl(fid);
    p_sfc = str2num(temp(1:9));
    th_sfc = str2num(temp(13:20));
    qv_sfc = str2num(temp(22:end));
    clear temp

    
%%initial sounding data: 1) zheight (m); 2) theta (K); 3) qv (g/kg); 4) u (m/s); 5) v (m/s)
[C pos] = textscan(fid,'%f %f %f %f %f');
zz00 = C{1}(1:end);
th00 = C{2}(1:end);
qv00 = C{3}(1:end);
u00  = C{4}(1:end);
v00  = C{5}(1:end);
clear C pos

%%put into proper units
qv00=qv00/1000; %[kg/kg]
p_sfc=p_sfc*100;    %[Pa]
qv_sfc=qv_sfc/1000; %[kg/kg]

dz = zz00(2)-zz00(1); %[km] assumes constant vertical resolution
nz = length(zz00);  %number of levels in sounding
zz00=zz00';
th00=th00';
qv00=qv00';
u00=u00';
v00=v00';

%%calculate virtual potential temperature from theta and wv mix ratio
thv_sfc = th_sfc*(1+qv_sfc/eps)/(1+qv_sfc);
thv00=th00.*(1+qv00/eps)./(1+qv00);

%%calculate pressure from pi (via hydrostatic equation) -- see base.F line 587
pi_sfc = (p_sfc/p0)^(Rd/Cpd);
pi00(1)=pi_sfc-g*zz00(1)/(Cpd*0.5*(thv_sfc+thv00(1)));
for k=2:nz
  pi00(k)=pi00(k-1)-g*(zz00(k)-zz00(k-1))/(Cpd*0.5*(thv00(k)+thv00(k-1)));
end
pp00=p0*(pi00.^(Cpd/Rd));

%load RE87_pres
%%var: p_RE87
%pp00 = p_RE87';

%%calculate actual temperature from potential temperature
T00 = th00.*(pp00/p0).^(Rd/Cpd);

%%calculate virtual temperature from actual temperature and wv mix ratio
Tv00 = T00.*(1+.608*qv00); %virtual potential temperature

%%calculate density from pressure and virtual temperature
rho00 = pp00./(Rd*Tv00);

%%calculate saturation vapor pressures at each level
T00_C = T00-273.15;
es00 = 6.112*exp((17.67*T00_C)./(T00_C+243.5))*100;  %[Pa]

%%calculate saturation mixing ratios at each level
qvs00 = eps*(es00./(pp00-es00));  %[kg/kg]

%%calculate relative humidity at each level
rh00 = qv00./qvs00;

end