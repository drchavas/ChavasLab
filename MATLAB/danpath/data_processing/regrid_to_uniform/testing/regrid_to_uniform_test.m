%regrid_to_uniform_test.m
%Purpose: Take input set of data with asymmetric grid-spacing and regrid to
%   symmetric grid-spacing
%
% Syntax:  
%
% Inputs:
%   xx [dist or lon] - matrix of x coordinates
%
% Example: 
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 29 May 2014; Last revision:

% Revision history:

%------------- BEGIN CODE --------------

clear
clc
close all
set(0,'DefaultFigureVisible', 'on');

addpath(genpath('~/Dropbox/Research/MATLAB/'));

%%Load some test data
load('testing/azimdattemp.mat')

rr_mean = nanmean(rr_data)

dr_nativeres = 100*1000;    %[m]; approximate native model resolution
C = 2*pi*rr_mean;
N_bins = floor(C/(sqrt(2)*dr_nativeres))

dthbinedges = (1/N_bins)*360;
%{
tic
thbinmin = 0;
thbinmax = thbinmin+dthbinedges;
for jj=1:N_bins
    
    Vtot_mean_azim(jj) = nanmean(V_data(th_data>=thbinmin & th_data<thbinmax));
    
    thbinmin = thbinmin+dthbinedges;
    thbinmax = thbinmax+dthbinedges;
end
toc
%}

tic

thbinedges = 0:dthbinedges:360;

[th_data i_sort] = sort(th_data);
V_data = V_data(i_sort);
[N_inbin Bin] = histc(th_data,thbinedges);
binnum = 1:N_bins;

mat_temp = Bin*(binnum./(binnum.^2));
mat_temp(mat_temp~=1) = 0;

assert(sum(sum(mat_temp==1))==length(th_data),'Something wrong with azim-mean binning algorithm')

Vtot_mean_azim = (V_data'*mat_temp)./N_inbin(1:end-1)';
%}
%{
%% Slow looping method
icount = 0;
for jj=1:N_bins
    
    N_inbin_temp = N_inbin(jj); %number of storms in this bin
    
    Vtot_mean_azim(jj) = nanmean(V_data(icount+1:icount+N_inbin_temp));
   
    icount = icount + N_inbin_temp;
end

%}
toc
%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vtot_mean_azim
Vtot_mean = nanmean(Vtot_mean_azim)
Vtot_mean_azim_gooddata = ~isnan(Vtot_mean_azim)

figure(1)
scatter(th_data,V_data,'gx')
hold on
plot(thbinedges,0*thbinedges,'m*')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%------------- END OF CODE --------------