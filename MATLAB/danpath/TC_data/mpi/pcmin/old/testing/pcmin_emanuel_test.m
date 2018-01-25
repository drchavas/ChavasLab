%pcmin_emanuel_test.m
%Purpose: Test your PI code against a standard!
% The standard is the accompanying file 'pcmin_emanuel_ver2014-10-30.m', whose
% outputs have been modified by DRC 2015-08-24, but the core code is the
% same as the version from Kerry Emanuel's website dated 2014-10-30.
% The input data are totally arbitrary, but give a reasonably high Vp
% value.

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 2015-08-24; Last revision: 2015-09-28
% 2015-09-28 - updated to use mixing ratio (r) instead of specific humidity (q)

%------------- BEGIN CODE --------------

clc
clear
close all

%% USER INPUT %%%%%%%%%%%%%%%
pcmin_emanuel_testfile = 'pcmin_emanuel_testdata.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load input data, saved values for comparison
listOfVariables  = {'Tsst_in','Pmsl_in','prs_in','T_in','r_in',...
    'PMIN_drc','VMAX_drc','AIRSEA_drc','TO_drc','IFL_drc',...
    'CKCD','SIG','IDISS','b','NK','VREDUC',...
    'source','version_date'
    };
load(pcmin_emanuel_testfile,listOfVariables{:})
sprintf('Loading input data and calculated Vp result for comparison from %s',pcmin_emanuel_testfile)


%% Calculate PI
% qv_in = qv_in/1000; %[kg/kg], specific humidity
% r_in = 1000*(qv_in./(1-qv_in));    %[g/kg]; water vapor mixing ratio -- Emanuel (1994) p. 108
[PMIN,VMAX,AIRSEA,TO,IFL]= pcmin_emanuel_ver20141030(Tsst_in,Pmsl_in,prs_in,T_in,r_in); %takes MIXING RATIO as input, not specific humidity!
               
sprintf('Vmax = %3.1f m/s',VMAX)
sprintf('Pmin = %3.1f hPa',PMIN)

%% Compare with saved value
assert(VMAX == VMAX_drc,'Vmax value does not match saved value')
sprintf('Vmax value matches saved value')

assert(PMIN == PMIN_drc,'Pmin value does not match saved value')
sprintf('Pmin value matches saved value')