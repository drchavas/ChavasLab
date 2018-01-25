%temp_frompottemp.m
%[T] = temp_frompottemp(p,th,r)

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

function [T] = temp_frompottemp(p,th,r)

%%If needed, load constants from CM1
constants_CM1_createdatfile();  %creates constants_CM1_list
load constants_CM1_list

T = th.*(p./p00).^((rd/cp).*((1+r/epsil)./(1+r*(cpv/cp))));  %Emanuel (1994) p. 111 Eq. 4.2.11 no approximation
%T = th.*(p./p00).^((rd/cp).*(1-.24*r));  %Emanuel (1994) p. 111 Eq. 4.2.11 approximation

end


%------------- END OF CODE --------------