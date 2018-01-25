function [shr,chim,PI,PI_p,vi,asdeq,iflag] = insane_calcVent3(basin,yrblock)
% Script for calculating ventilation index: Merging calcVent with (c) BRIAN TANG (4/15/14)
% input: vector of lon,lat,date for each observation
% CONVERT LONS


%CONSTANTS
Rd = 287.04;%J/kg-K, gas const for dry air
Rv = 461.50;%J/kg-K, gas const for water vap
eps = Rd/Rv;
cpd = 1005.7;%J/kg-K, specific heat at constant pressure for dry air
Lv = 2.555E6; %use inflated value of Lv (Bryan 2008) for pseudoadiabatic entropy calculation
                %J/kg, latent heat of vaporization, compensating for neglect of ent. of water vap
% data path:  /Users/renzhijing/Documents/Work/Data/ECMWF  
% data path:  /Users/renzhijing/Documents/NHRA/TE12/ECMWF_Data/  
if strcmp(basin, 'NA')
    switch yrblock
        case 4
%     -2002-2010
load /Users/renzhijing/Documents/Work/Data/ECMWF/T_na_0214_new.mat temp
load /Users/renzhijing/Documents/Work/Data/ECMWF/qa_na_0214_new.mat qa
load /Users/renzhijing/Documents/Work/Data/ECMWF/uv_na_0214_new.mat
load /Users/renzhijing/Documents/Work/Data/ECMWF/sst_na_0214_new.mat sst %hC [global]
load /Users/renzhijing/Documents/Work/Data/ECMWF/psl_na_0214_new.mat psl %Pa [global] 
        case 3
%     -1990-2001
load /Users/renzhijing/Documents/Work/Data/ECMWF/T_na_9001_new.mat temp
load /Users/renzhijing/Documents/Work/Data/ECMWF/qa_na_9001_new.mat qa
load /Users/renzhijing/Documents/Work/Data/ECMWF/uv_na_9001_new.mat
load /Users/renzhijing/Documents/Work/Data/ECMWF/sst_na_9001_new.mat sst %hC [global]
load /Users/renzhijing/Documents/Work/Data/ECMWF/psl_na_9001_new.mat psl %Pa [global] 


        case 2
%     - 1981-1989
load /Users/renzhijing/Documents/Work/Data/ECMWF/T_na_8189_new.mat temp
load /Users/renzhijing/Documents/Work/Data/ECMWF/qa_na_8189_new.mat qa
load /Users/renzhijing/Documents/Work/Data/ECMWF/uv_na_8189_new.mat
load /Users/renzhijing/Documents/Work/Data/ECMWF/sst_na_8189_new.mat sst %hC [global]
load /Users/renzhijing/Documents/Work/Data/ECMWF/psl_na_8189_new.mat psl %Pa [global] 
        case 1
%     - 1970-1980
load /Users/renzhijing/Documents/Work/Data/ECMWF/T_na_7980_new.mat temp
load /Users/renzhijing/Documents/Work/Data/ECMWF/qa_na_7980_new.mat qa 
load /Users/renzhijing/Documents/Work/Data/ECMWF/uv_na_7980_new.mat
load /Users/renzhijing/Documents/Work/Data/ECMWF/sst_na_7980_new.mat sst %hC [global]
load /Users/renzhijing/Documents/Work/Data/ECMWF/psl_na_7980_new.mat psl %Pa [global] 
    end



elseif strcmp(basin,'WNP')
    switch yrblock 
        case 4
%     -2002-2010
load /Users/renzhijing/Documents/Work/Data/ECMWF/T_wnp_0214_new.mat temp
load /Users/renzhijing/Documents/Work/Data/ECMWF/qa_wnp_0214_new.mat qa
load /Users/renzhijing/Documents/Work/Data/ECMWF/uv_wnp_0214_new.mat 
load /Users/renzhijing/Documents/Work/Data/ECMWF/sst_wnp_0214_new.mat sst %hC [global]
load /Users/renzhijing/Documents/Work/Data/ECMWF/psl_wnp_0214_new.mat psl %Pa [global]
        case 3
%     -1990-2001
load /Users/renzhijing/Documents/Work/Data/ECMWF/T_wnp_9001_new.mat temp
load /Users/renzhijing/Documents/Work/Data/ECMWF/qa_wnp_9001_new.mat qa
load /Users/renzhijing/Documents/Work/Data/ECMWF/uv_wnp_9001_new.mat 
load /Users/renzhijing/Documents/Work/Data/ECMWF/sst_wnp_9001_new.mat sst %hC [global]
load /Users/renzhijing/Documents/Work/Data/ECMWF/psl_wnp_9001_new.mat psl %Pa [global]

        case 2
%     -1981-1989
load /Users/renzhijing/Documents/Work/Data/ECMWF/T_wnp_8189_new.mat temp
load /Users/renzhijing/Documents/Work/Data/ECMWF/qa_wnp_8189_new.mat qa
load /Users/renzhijing/Documents/Work/Data/ECMWF/uv_wnp_8189_new.mat 
load /Users/renzhijing/Documents/Work/Data/ECMWF/sst_wnp_8189_new.mat sst %hC [global]
load /Users/renzhijing/Documents/Work/Data/ECMWF/psl_wnp_8189_new.mat psl %Pa [global]

        case 1
%     -1970-1980
load /Users/renzhijing/Documents/Work/Data/ECMWF/T_wnp_7980_new.mat temp
load /Users/renzhijing/Documents/Work/Data/ECMWF/qa_wnp_7980_new.mat qa 
load /Users/renzhijing/Documents/Work/Data/ECMWF/uv_wnp_7980_new.mat 
load /Users/renzhijing/Documents/Work/Data/ECMWF/sst_wnp_7980_new.mat sst %hC [global]
load /Users/renzhijing/Documents/Work/Data/ECMWF/psl_wnp_7980_new.mat psl %Pa [global]  

    end

end

rv = qa.qa./(1 - qa.qa);   % Revised by Renzhi
plev = logical(qa.Z==600);
pm= logical(qa.Z==600);
rv600 = squeeze(rv(:,:,plev,:));    % Renzhi: rv600: kg/kg
rv = rv.*1000.0;        % Renzhi: rv: g/kg
sst.sst(sst.sst==-999)=NaN;
sstC = sst.sst - 273.15;
pmb = temp.Z;
tC = temp.temp-273.15;
slpmb = psl.psl./100.0;

if length(pmb)==12
    pmb(10) = []; tC(:,:,10,:) = [];
elseif length(pmb)== 17
    pmb(10) = []; tC(:,:,10,:) = [];
    pmb(12:16) = []; tC(:,:,12:16,:) = [];
end
sprintf('LOADING DATA FINISHED for yr = %d',yrblock)

lon = double(U.X);
lat = double(U.Y);
time = U.T;
nlon = length(lon);
nlat = length(lat);
ntime = length(time);
[LON,LAT] = meshgrid(lon,lat);

PI_p = nan(size(sstC));
PI   = nan(size(sstC));
iflag = nan(size(sstC));
asdeq = nan(size(sstC));

for j = 1:nlat 
   for i = 1:nlon
        for t = 4:ntime
              rgkg = [reshape(rv(j,i,:,t-3),8,1);0;0;0];
              [pmin,vmax,airsea,~,ifl] = pcmin_new(sstC(j,i,t-3),slpmb(j,i,t-3),double(pmb),tC(j,i,:,t-3),rgkg); %output: pmin,vmax,airsea are numbers,
              PI_p(j,i,t) = pmin;
              PI(j,i,t) = vmax;
              asdeq(j,i,t) = airsea;
              iflag(j,i,t) = ifl;
        end
   end
   sprintf('Calculating MPI for yr = %d,lat = %d',yrblock,j)
end

sprintf('MPI Calculation Finished')

%NONDIMENSIONAL ENTROPY DEFICIT

udiff = U.U(:,:,1,:) - U.U(:,:,2,:);   %850pa-250pa
vdiff = V.V(:,:,1,:) - V.V(:,:,2,:);
shr = squeeze(sqrt(udiff.^2 + vdiff.^2));   %shr(lat,lon,1,time)


r1 = 100E3; %hypothetical inner-core storm radius (m)
r2 = 300E3; %hypothetical outer radius for environmental air intrusion (m)
th = 0:(pi/4):2*pi-pi/4;
[TH1,R1] = meshgrid(th,(50E3:50E3:r1)./110E3);
TH1 = reshape(TH1,[],1);
R1 = reshape(R1,[],1);
[TH2,R2] = meshgrid(th,(r1:100E3:r2)./110E3);
TH2 = reshape(TH2,[],1);
R2 = reshape(R2,[],1);

[xxi,yyi] = pol2cart(TH1,R1);
[xxo,yyo] = pol2cart(TH2,R2);

pmid = 60000; %entropy deficit level (Pa)

sstarminnerbar = nan(size(PI));
smouterbar = nan(size(PI));

for t = 4:ntime
    %NONDIMENSIONAL ENTROPY DEFICIT
    %s*_m
    tm = squeeze(temp.temp(:,:,pm,t));
    estarm = squeeze(611.2.*exp(17.67.*(tm-273.15)./(243.5+(tm-273.15)))); %Pa
    rstarm = eps*estarm./(pmid-estarm);
    pd = pmid./(1+rstarm./eps);
    sstarm = cpd.*log(tm)-Rd.*log(pd)+Lv.*rstarm./tm;

    %Average over (hypothetical) inner-core disc
    F = TriScatteredInterp(reshape(LON,[],1),reshape(LAT,[],1),reshape(double(sstarm),[],1));
    for j = 1:length(lat)
        for i = 1:length(lon)
            if (isnan(PI(j,i,t)))
                continue;
            end
            disclon = LON(j,i)+xxi;
            disclat = LAT(j,i)+yyi;
            sstarminner = F(disclon,disclat);
            sstarminnerbar(j,i,t) = nansum(sstarminner.*R1)/sum(~isnan(sstarminner).*R1);
        end
    end
    
    %s_m
    rm = squeeze(rv600(:,:,t));
    pd = pmid./(1+rm./eps);
    vappres = rm.*pd./eps;
    rh = min(vappres./estarm,1);
    sm = cpd.*log(tm)-Rd.*log(pd)+Lv.*rm./tm-Rv.*rm.*log(max(rh,0.01));

    %Average over (hypothetical) environmental annulus
    F = TriScatteredInterp(reshape(LON,[],1),reshape(LAT,[],1),reshape(double(sm),[],1));
    for j = 1:length(lat)
        for i = 1:length(lon)
            if (isnan(PI(j,i,t)))
                continue;
            end
            disclon = LON(j,i)+xxo;
            disclat = LAT(j,i)+yyo;
            smouter = F(disclon,disclat);
            smouterbar(j,i,t) = nansum(smouter.*R2)/sum(~isnan(smouter).*R2);
        end
    end
    sprintf('vi calculation for yr = %d, time = %d',yrblock,t)

end


%do not allow supersaturation (negative numerator); also limit smallness of denominator
chim = max(0,sstarminnerbar-smouterbar)./max(1,asdeq);
%eliminate points where net enthalpy flux into ocean
chim(asdeq<0) = NaN;

%VENTILATION INDEX
vi = shr.*chim./max(PI,1);
vi(isnan(PI)) = NaN;
vi(PI==0) = NaN;
VREDUC=0.75;
PI = PI.*VREDUC;