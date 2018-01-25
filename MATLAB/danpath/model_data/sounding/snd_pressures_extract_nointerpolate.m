%snd_pressures_extract_nointerpolate.m

%Purpose: extract input_sounding data and linearly interpolate it to the
%desired model levels as is done in CM1. For the pressure field, this code
%outputs only the hydrostatically balanced pressure field associated with
%the distribution of theta, qv, and p_sfc.  (this is the background state CM1 calculates)

%Updated: 09 Feb 2012 DRC

%clear
%clc

function [pp00_full] = snd_pressures_extract_nointerpolate(dir_in,snd_file)

%dir_in = '/Volumes/CHAVAS_CM1/CM1_output/axisym/TRANSFER/CTRLsnd2540';
%snd_file = 'input_sounding';

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
sounding_fil = sprintf('%s/%s_pressures',dir_in,snd_file);
clear fid
fid = fopen(sounding_fil);

%%no header
    
%%initial sounding data: 1) pressure [Pa]
[C pos] = textscan(fid,'%f');
pp00_full = C{1}(1:end);
clear C pos

end