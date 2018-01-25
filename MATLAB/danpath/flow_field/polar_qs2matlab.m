%polar_qs2matlab.m - vectorized conversion of polar coordinate angle
%convert theta from QuikSCAT format to MATLAB format
%MATLAB: counterclockwise from east in signed radians (-pi,pi]
%QuikSCAT: clockwise from North in positive degrees [0,360)
%
% Syntax:  [th_mat] = polar_qs2matlab(th_qs)
%
% Inputs:
%    th_qs - vector angle in QuikSCAT format
%
% Outputs:
%    th_mat - vector angle in MATLAB format
%
% Example: 
%    [th_mat] = polar_qs2matlab(th_qs);
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
% 13 Nov 2013; Last revision: 13 Nov 2013

%------------- BEGIN CODE --------------

function [th_mat] = polar_qs2matlab(th_qs)

%%deg CW from N --> rad CCW from E signed
th_qs = (360 - th_qs) + 90; %reverse direction and shift by 90 deg
th_qs(th_qs>=360) = th_qs(th_qs>=360) - 360;    %coterminal angle for >= 360
th_qs(th_qs>180) = th_qs(th_qs>180) - 360;   %shift to (-180,180]
th_mat = pi*th_qs/180;  %convert to radians
clear th_qs

assert(min(min(th_mat))>-pi & max(max(th_mat))<=pi,'th_mat out of bounds')

end

%------------- END OF CODE --------------


