%temp_frompottemp_dry.m

%Purpose: calculate the basic dry potential temperature given temp and pressure

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

function [T] = temp_frompottemp_dry(p,theta)

%%Load constants from CM1
const = CM1_constants();  %creates constants_CM1_list

T = theta.*(p./const.p00).^(const.Rd/const.cp);  %Emanuel (1994) p. 111 Eq. 4.3.2 for virtual temp

end


%------------- END OF CODE --------------