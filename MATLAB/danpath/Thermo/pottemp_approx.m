%pottemp_approx.m

%Purpose: calculate the potential temperature given temp and pressure

% Inputs:
%
% Outputs:
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
% 8 Apr 2014; Last revision:

%------------- BEGIN CODE --------------

function [theta] = pottemp_approx(p,T)

%%Load constants from CM1
addpath(genpath('~/Dropbox/Research/MATLAB/'));
const = CM1_constants();  %creates constants_CM1_list

%theta = T.*(const.p00./p).^((const.Rd/const.cp).*((1+r/const.eps)./(1+r*(const.cpv/const.cp))));  %Emanuel (1994) p. 111 Eq. 4.2.11 no approximation
%theta = T.*(const.p00./p).^((const.Rd/const.cp).*(1-.24*r));  %Emanuel (1994) p. 111 Eq. 4.2.11 approximation
theta = T.*(const.p00./p).^((const.Rd/const.cp));  %Dry only


end


%------------- END OF CODE --------------