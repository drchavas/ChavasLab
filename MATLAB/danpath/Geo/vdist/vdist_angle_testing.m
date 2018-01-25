%vdist_angle_testing.m
%Purpose: test output for angles between two points on sphere from function
%   VDIST with direct calculation using spherical law of cosines
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    ncdir_dr - the directory of the file you'd like
%    ncfile_in - the file you'd like
%    variable_in - the variable you'd like to extract
%    missing_value_flag - value that will be set to NaN
%
% Outputs:
%    data_out - matrix of the desired data
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
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 11 Nov 2013; Last revision: 11 Nov 2013

%------------- BEGIN CODE --------------

clear
clc
close all
set(0,'DefaultFigureVisible', 'on');

addpath(genpath('~/Dropbox/Research/MATLAB/'));


yy_center_matrix = 0;  %[deg N]
xx_center_matrix = 0;   %[deg E]
yy_mat = 89.99;  %[deg N]
xx_mat = 20;   %[deg E]

[rr, th_pt1pt2, th_pt2pt1] = vdist(yy_center_matrix,xx_center_matrix,yy_mat,xx_mat);


dlam = xx_mat - xx_center_matrix;

%% Predict th_pt2pt1 using spherical law of cosines
r_unit = 1;  %[m]
r_cent2NP = r_unit*acos(sind(yy_center_matrix));    %arc length from center point to north pole
th_pt2pt1_predict = acosd(-cosd(th_pt1pt2)*cosd(dlam)+sind(th_pt1pt2)*sind(dlam)*cos(r_cent2NP));    %spherical law of cosines

figure(104)
plot(xx_mat,yy_mat,'rx','MarkerSize',20)
hold on
plot(xx_center_matrix,yy_center_matrix,'g.','MarkerSize',20)
xlabel('longitude')
ylabel('latitude')
title(sprintf('r = %5.0f km, th12 = %3.0f deg, th21 = %3.0f deg, th21pred = %3.0f deg',rr/1000,th_pt1pt2,360-th_pt2pt1,th_pt2pt1_predict))


%------------- END OF CODE --------------