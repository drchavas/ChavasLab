%calcVI.m: calculate the ventilation index given input data
%Based on original file 'ventilation_sample.m' by (c) BRIAN TANG (4/15/14)
%Code calculates PI, shear, and chi at all input points
%
% Syntax:  
%
% Inputs:
%
% Outputs: none, file saved in current directory
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
% Notes:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 6 Aug 2014; Last revision:

%------------- BEGIN CODE --------------

%CONSTANTS
Rd = 287.04;
Rv = 461.50;
eps = Rd/Rv;
cpd = 1005.7;
Lv = 2.555E6; %use inflated value of Lv (Bryan 2008) for pseudoadiabatic entropy calculation

%USER DEFINED CONSTANTS
%Shear
plow = 85000; %shear layer lower level (Pa)
pupp = 20000; %shear layer upper level (Pa)

%Entropy deficit
pmid = 60000; %entropy deficit level (Pa)
r1 = 100E3; %hypothetical inner-core storm radius (m)
r2 = 300E3; %hypothetical outer radius for environmental air intrusion (m)

%Setup grid for interpolation to polar coordinates for entropy deficit calculation
th = 0:(pi/4):2*pi-pi/4;
[TH1,R1] = meshgrid(th,(50E3:50E3:r1)./110E3);
TH1 = reshape(TH1,[],1);
R1 = reshape(R1,[],1);
[TH2,R2] = meshgrid(th,(r1:100E3:r2)./110E3);
TH2 = reshape(TH2,[],1);
R2 = reshape(R2,[],1);

[xxi,yyi] = pol2cart(TH1,R1);
[xxo,yyo] = pol2cart(TH2,R2);

%READ IN LATITUDE, LONGITUDE, AND PRESSURE LEVELS HERE
%lat = ;
%lon = ;
%[LON,LAT] = meshgrid(lon,lat);
%p = ;

%READ IN STATE VARIABLES HERE (assuming a single time, loop over time if needed)
%temperature
%t = ; %(pressure, lat, lon)

%water vapor mixing ratio
%r = ; %(pressure, lat, lon)

%sea surface temperature
%sst = ; %(lat, lon)

%mean sea level pressure
%slp = ; %(lat, lon)
       
%zonal wind
%u = ; %(pressure, lat, lon)
       
%meridional wind
%v = ; %(pressure, lat, lon)

%WIND SHEAR
%NOTE: Normally you will have to filter out any TC circulation (if one exists) before calculating the shear
%For seasonal climatologies in GCMs, reanalyses, this is less important
p2 = find(p==pupp);
p1 = find(p==plow);
shear = squeeze(sqrt((u(p2,:,:)-u(p1,:,:)).^2+(v(p2,:,:)-v(p1,:,:)).^2));

%POTENTIAL INTENSITY
%NOTE: Make sure SST data conforms to atmospheric grid
pitemp = zeros(size(sst));
asdeqtemp = zeros(size(sst));
sstC = sst-273.15;
slpmb = slp./100;
pmb = p./100;
tC = t-273.15;
rgkg = r.*1000;
for phi = 1:length(lat)
    for lam = 1:length(lon)
    	if (isnan(sst(phi,lam)))
        	continue;
        end
        %Call modified potential intensity algorithm (pcmin.m) to calculate the potential intensity and denominator of nondimensional entropy deficit
        %See pcmin.m for details
        %[pmin,vmax,airsea,ifl] = pcmin(sstC(phi,lam),slpmb(phi,lam),pmb,tC(:,phi,lam),rgkg(:,phi,lam),0.7,1);

        [pmin,vmax,airsea,~,~]= pcmin_emanuel(Tsst_in,Pmsl_in,prs_in,T_in,qv_in);    %Units: Tsst/T [C]; Pmsl/P [hPa]; qv [g/kg]; lowest model level to highest model level
            
        
        pi(phi,lam) = vmax;
        asdeq(phi,lam) = airsea;
    end
end

%NONDIMENSIONAL ENTROPY DEFICIT
%s*_m
pm = find(p==pmid);
tm = squeeze(t(pm,:,:));
estarm = squeeze(611.2.*exp(17.67.*(tm-273.15)./(243.5+(tm-273.15)))); %Pa
rstarm = eps*estarm./(pmid-estarm);
pd = pmid./(1+rstarm./eps);
sstarm = cpd.*log(tm)-Rd.*log(pd)+Lv.*rstarm./tm;

%Average over (hypothetical) inner-core disc
sstarminnerbar = nan(size(sstarm));
F = TriScatteredInterp(reshape(LON,[],1),reshape(LAT,[],1),reshape(double(sstarm),[],1));
for nn = 1:length(lat)
	for mm = 1:length(lon)
        if (isnan(pi(nn,mm)))
            continue;
        end
        disclon = LON(nn,mm)+xxi;
        disclat = LAT(nn,mm)+yyi;
        sstarminner = F(disclon,disclat);
        sstarminnerbar(nn,mm) = nansum(sstarminner.*R1)/sum(~isnan(sstarminner).*R1);
    end
end

%s_m
rm = squeeze(r(pm,:,:));
pd = pmid./(1+rm./eps);
vappres = rm.*pd./eps;
rh = min(vappres./estarm,1);
sm = cpd.*log(tm)-Rd.*log(pd)+Lv.*rm./tm-Rv.*rm.*log(max(rh,0.01));

%Average over (hypothetical) environmental annulus
smouterbar = nan(size(sm));
F = TriScatteredInterp(reshape(LON,[],1),reshape(LAT,[],1),reshape(double(sm),[],1));
for nn = 1:length(lat)
	for mm = 1:length(lon)
        if (isnan(pi(nn,mm)))
            continue;
        end
        disclon = LON(nn,mm)+xxo;
        disclat = LAT(nn,mm)+yyo;
        smouter = F(disclon,disclat);
        smouterbar(nn,mm) = nansum(smouter.*R2)/sum(~isnan(smouter).*R2);
    end
end

%do not allow supersaturation (negative numerator); also limit smallness of denominator
chim = max(0,sstarminnerbar-smouterbar)./max(1,asdeq);
%eliminate points where net enthalpy flux into ocean
chim(asdeq<0) = NaN;

%VENTILATION INDEX
vi = shear.*chim./max(pi,1);    %PI set to minimum of 1 so it doesn't blow up DRC
vi(isnan(pi)) = NaN;
vi(pi==0) = NaN;

%------------- END OF CODE --------------