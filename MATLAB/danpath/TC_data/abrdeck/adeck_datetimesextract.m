%adeck_datetimesextract.m - extract date/times from adeck forecast file
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
% 12 Sep 2014; Last revision:

%------------- BEGIN CODE --------------
% 
% clear all
% clc
% close all
% addpath(genpath('~/Dropbox/Research/MATLAB/'));

function [TC_datetimes] = adeck_datetimesextract(dir_in,file_in)

%% Read all text in file into single big cell array (columns are variables)
path_in = sprintf('%s/%s',dir_in,file_in);
fid = fopen(path_in,'r');  % Open text file
%%ATCF format: BASIN,CY,YYYYMMDDHH,TECHNUM/MIN,TECH,TAU,LatN/S,LonE/W,VMAX,MSLP,TY,RAD,WINDCODE,RAD1,RAD2,RAD3,RAD4,RADP,RRP,MRD,GUSTS,EYE,SUBREGION,MAXSEAS,INITIALS,DIR,SPEED,STORMNAME,DEPTH,SEAS,SEASCODE,SEAS1,SEAS2,SEAS3,SEAS4,USERDEFINED,userdata
%%37 items
%{
line_format = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
btk_all = textscan(fid, line_format, 'Delimiter', ',');
%}
%%Only need the third line
data_temp = textscan(fid, '%*s %*s %s %*[^\n]', 'Delimiter', ',');
fclose(fid);

%% Extract datetime data
%%YYYYMMDDHH (analysis date-time)
TC_datetimes = data_temp{1};
clear data_temp

%------------- END CODE --------------