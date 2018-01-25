%Cd_Smith2008Black07.m -- wind-dependent drag coefficient
%Purpose: return drag coefficient based on Smith et al. (2014) QJRMS
%
% Syntax:  [C_d] = Cd_Smith2008Black07(V_in)
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
% 28 Mar 2014; Last revision: 

% Revision history:
%------------- BEGIN CODE --------------

function [C_d] = Cd_Smith2008Black07(V_in)

C_d0 = .7e-3;
C_d1 = 6.5e-5;
C_dmax = 2e-3*ones(size(V_in));
C_d = C_d0 + C_d1*V_in;
C_d = min([C_d;C_dmax],[],1);

%%%%%%%%%%%%%%%%%%
%%TESTING
% figure
% plot(V_in,C_d)
%%%%%%%%%%%%%%%%%%%%%
%------------- END OF CODE --------------
