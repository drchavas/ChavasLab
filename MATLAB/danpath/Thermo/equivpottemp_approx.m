%equivpottemp_approx.m

%Purpose: calculate the equivalent potential temperature given temp and pressure
%   ignores contribution to density from water

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

function [thetae] = equivpottemp_approx(theta,T,r)

%%Load constants from CM1
const = CM1_constants();  %creates constants_CM1_list

thetae = theta.*exp((const.Lv0.*r)./(const.cp.*T));   %approximate form

end


%------------- END OF CODE --------------