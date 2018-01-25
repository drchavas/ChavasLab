%weights_periodic_test.m
%Purpose: Assign weights based on spacing between datapoints
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
% i_test = [5 10];
% V_data = V_data(i_test);
% th_data = th_data(i_test);

[th_data,i_sort] = sort(th_data);
V_data = V_data(i_sort);
unsorted = 1:length(th_data);
i_unsort(i_sort) = unsorted;
clear unsorted

th_data_temp = [th_data(end)-360;th_data;th_data(1)+360];

mean_spacing = (th_data_temp(3:end)-th_data_temp(1:end-2))/2;
assert(sum(mean_spacing)==360,'Problem with mean spacing algorithm')
normfac = 360;
data_weights_sorted = mean_spacing/normfac;

data_weights = data_weights_sorted(i_unsort)

assert(sum(data_weights_sorted~=data_weights(i_sort))==0,'problem with unsorting algorithm')

%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
clf(1)
scatter(th_data,V_data,50,data_weights_sorted)
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------- END OF CODE --------------