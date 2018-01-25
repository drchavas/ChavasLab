%Donelan2004_fit.m - piecewise-linear fit to Donelan et al 2004 Cd(V) data
%
%Data from Donelan2004.m -- email from Mark Donelan 2014-06-17
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
% 17 Jun 2014; Last revision:

%------------- BEGIN CODE --------------

%function [C_d] = Donelan2004_fit(V_in)

clear
clc
close all
set(0,'DefaultFigureVisible', 'on');

addpath(genpath('~/Dropbox/Research/MATLAB/'));

figure(4)
clf(4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%From Donelan2004.m -- email from Mark Donelan 2014-06-17
% Donelan et al 2004 GRL: CD vs U10; 
% Data lifted off Figure 2 by Alex Soloviev.
%Squares
VD1 =  ([263 251 225 223 198 194 174 166 152 139 116 130 111 103   97  92  85  74  69  64   63  60  58]-46)*60/(560-46);
CD1 = -([289 296 309 310 324 326 338 342 350 358 370 362 371 376  379 380 382 382  379 372 367 358 341]-431)*5e-3/(431-26);
%Asterisks
VD2 = ([180 226 276 334 412 400 376]-46)*60/(560-46);
CD2 = -([351 324 278 235 246 246 237]-431)*5e-3/(431-26);
%Circles
VD3 = ([209 502 476 384 254 307 429 356]-46)*60/(560-46);
CD3 = -1.0*([340 244 231 255 305 265 241 241]-431)*5e-3/(431-26);
%Diamonds
VD4 = ([396 433 315 267 183 148 354 220 84]-46)*60/(560-46);
CD4 = -([237 237 248 274 324 341 243 307 381]-431)*5e-3/(431-26);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Piecewise linear fit: constant V<4 m/s, V>33 m/s, linear in between
V_all = [VD1 VD2 VD3 VD4];
Cd_all = [CD1 CD2 CD3 CD4];

V_breakpoints = [6 35.4 60];    %ignore higher Cd data for small V!
Cd_breakpoints = lsq_lut_piecewise( V_all', Cd_all',  V_breakpoints');

linear_slope = (Cd_breakpoints(2)-Cd_breakpoints(1))/(V_breakpoints(2)-V_breakpoints(1))

yintercept = Cd_breakpoints(2) - linear_slope*V_breakpoints(2)

%% Run Cd_Donelan04 to make sure it matches up
V_in = 1:60;
[Cd_D04test] = Cd_Donelan04(V_in);

%% PLOT THINGS
plot(V_all,Cd_all,'.',[0 V_breakpoints],[Cd_breakpoints(1) Cd_breakpoints'],'+-')
hold on
plot(V_in,linear_slope*V_in+yintercept,'m--')
plot(V_in,Cd_D04test,'r')
title(sprintf('Piecewise linear fit (cnst/lin/cnst): breakpoints (%3.1f,%5.5f), (%3.1f,%5.5f); m=%6.6f',...
    V_breakpoints(1),Cd_breakpoints(1),V_breakpoints(2),Cd_breakpoints(2),linear_slope))

%% Save plot %%%%%%%%%%%%%%%%%%
plot_filename = sprintf('Donelan2004_fit.pdf')
saveas(gcf,plot_filename,'pdf')




%------------- END OF CODE --------------