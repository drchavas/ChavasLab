%pcmin_direct.m -- input data directly rather than reading sounding file

%Purpose: For a given sounding, calculate the potential intensity (both
%reversible and pseudoadiabatic)

function [Vp_reversible,Vp_pseudoadiabatic]= pcmin_direct(Tsst,p_sfc,th00,qv00,pp00_full,CKCD_in,IDISS_in,b_in,NK_in,VREDUC_in,ptop)

%Calculate both Vp_reversible Vp_pseudoadiabatic
SIG_in = [0 1]; %[0-1], 0=reversible, 1=pseudoadiabatic parcel ascent (can do in between)

%% TESTING %%%%%
%{
clc
clear
close('all')

%% USER INPUT %%%%%%%%
dir_in = '/Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/input_soundings/';
snd_file = 'input_sounding_axCTRLv0qrhSATqdz5000_nx3072_SST310.00K_drag';
Tsst = 310.0;   %[K]

%%Inputs into V_p code (NOTE: code cuts off sounding at ptop = 10 hPa)
CKCD_in = 1.0;    %Ck/Cd; ratio of exchange coefficients
SIG_in = 1; %[0-1], 0=reversible, 1=pseudoadiabatic parcel ascent (can do in between)
IDISS_in = 1;   %0=no dissipative heating, 1=yes (simply multiplies Vmax by Tsst/T0)
b_in = 2; %Exponent, b, in assumed profile of azimuthal velocity in eye, V=V_m(r/r_m)^b. Used only in calculation of central pressure   
NK_in = 1;  %sounding level from which parcels lifted
VREDUC_in = 1.0;    %Factor to reduce gradient wind to 10 m wind
ptop = 10;  %[hPa] checks up to this pressure level; if it blows up, code will increase this until it doesn't (up to 100 hPa, after which gives error)
%%%%%%%%%%%%%%%%%%%%
%}

%% Recalculate T00 using new pressures
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

T00_full = th00(1:length(pp00_full)).*(pp00_full/p0).^(Rd/Cpd);

SST_in = Tsst - 273.15; %[C]; sea surface temperature
PSL_in = p_sfc/100; %[hPa]; surface pressure
PP_in = pp00_full/100;   %[hPa]; pressure vector
TT_in = T00_full-273.15; %[C];   temperature vector
RR_in = qv00*1000;  %[g/kg]; mixing ratio vector

%% Reversible V_p
VMAX = 0;
iter = 0;
ptop_save = ptop;
while(VMAX<=0)  %on occasion the code blows up for some reason when ptop is too low (in magnitude)
    iter = iter + 1;
    ptop = ptop+10*(iter-1);
    [PMIN,VMAX,TOM,IFL]= pcmin_drc20130423(SST_in,PSL_in,PP_in,TT_in,RR_in,CKCD_in,SIG_in(1),IDISS_in,b_in,NK_in,VREDUC_in,ptop);
    if(iter>10)
        assert('ERROR IN CALCULATING MPI')
    end
end
Vp_reversible = VMAX;
ptop = ptop_save;

%% Pseudoadiabatic V_p
VMAX = 0;
iter = 0;
while(VMAX<=0)  %on occasion the code blows up for some reason when ptop is too low (in magnitude)
    iter = iter + 1;
    ptop = ptop+10*(iter-1);
    [PMIN,VMAX,TOM,IFL]= pcmin_drc20130423(SST_in,PSL_in,PP_in,TT_in,RR_in,CKCD_in,SIG_in(2),IDISS_in,b_in,NK_in,VREDUC_in,ptop);
    if(iter>10)
        assert('ERROR IN CALCULATING MPI')
    end
end
Vp_pseudoadiabatic = VMAX;

end