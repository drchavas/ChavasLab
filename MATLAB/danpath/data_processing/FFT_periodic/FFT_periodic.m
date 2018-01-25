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

% [th_data i_sort] = sort(th_data);
% V_data = V_data(i_sort);
% [N_inbin Bin] = histc(th_data,thbinedges);
% binnum = 1:N_bins;


h = edft(V_data)

%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
scatter(th_data,V_data,'gx')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%------------- END OF CODE --------------