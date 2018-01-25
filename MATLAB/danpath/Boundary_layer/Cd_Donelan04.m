%Cd_Donelan04.m -- wind-dependent drag coefficient
%Purpose: return vectorized drag coefficient based on Donelan (2004)
%
% Syntax:  [C_d] = Cd_Donelan04(V_in)
%
% Inputs:
%   V_in [ms-1] - vector of wind speeds
%
% Outputs:
%   C_d [-] - vector of corresponding drag coefficients
%
% Example: 
%   [C_d] = Cd_Donelan04(V_in)
%
% References: 
%   - Donelan et al. (2004)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 28 Mar 2014; Last revision: 17 Jun 2014

% Revision history:
% 17 Jun 2014 - Updated piecewise linear relationship to best fit data
%------------- BEGIN CODE --------------

function [C_d] = Cd_Donelan04(V_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Piecewise linear fit parameters estimated from Donelan2004_fit.m
C_d_lowV = 6.2e-4;
V_thresh1 = 6;  %m/s; transition from constant to linear increasing
V_thresh2 = 35.4;  %m/s; transition from linear increasing to constant
C_d_highV = 2.35e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linear_slope = (C_d_highV-C_d_lowV)/(V_thresh2-V_thresh1);

%%Vectorize for input data
C_d_lowV = C_d_lowV*ones(size(V_in));
C_d_midV = C_d_lowV + linear_slope*(V_in-V_thresh1);
C_d_highV = C_d_highV*ones(size(V_in));

%%Find true C_d value
C_d = max([min([C_d_midV' C_d_highV'],[],2) C_d_lowV'],[],2);

%%%%%%%%%%%%%%%%%%
%%TESTING
% figure
% plot(V_in,C_d)
%%%%%%%%%%%%%%%%%%%%%
%------------- END OF CODE --------------
