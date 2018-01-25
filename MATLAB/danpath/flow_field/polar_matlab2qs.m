%polar_matlab2qs.m - vectorized conversion of polar coordinate angle
%convert theta from MATLAB format to QuikSCAT format
%MATLAB: counterclockwise from east in signed radians (-pi,pi]
%QuikSCAT: clockwise from North in positive degrees [0,360)
%
% Syntax:  [th_qs] = polar_matlab2qs(th_mat)
%
% Inputs:
%    th_mat - vector angle in MATLAB format
%
% Outputs:
%    th_qs - vector angle in QuikSCAT format
%
% Example: 
%    [th_qs] = polar_matlab2qs(th_mat);
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

function [th_qs] = polar_matlab2qs(th_mat)

%%rad CCW from E signed --> deg CCW from E signed
th_mat = 180*th_mat/pi;  %convert to degrees (-180,180]
th_qs = 360 - (th_mat - 90);    %shift by 90 degrees and reverse direction
clear th_mat
th_qs(th_qs>=360) = th_qs(th_qs>=360)-360; %coterminal angle for >= 360

assert(min(min(th_qs))>=0 & max(max(th_qs))<360,'th_qs out of bounds')

end

%------------- END OF CODE --------------


