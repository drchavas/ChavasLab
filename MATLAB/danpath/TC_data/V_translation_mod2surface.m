%V_translation_mod2surface.m - surface-modified TC translation vector
%reduce magnitude and rotate input vector
%
% Syntax:  [th_translation_surface V_translation_surface] = ...
%           V_translation_mod2surface(V_translation,th_translation,...
%           V_translation_reduction_factor,V_translation_rotation_angle,fcor_sign)
%
% Inputs:
%    V_translation - vector magnitude
%    th_translation - vector angle [deg CW from N (0,360]]
%    V_translation_reduction_factor - magnitude reduction factor
%    V_translation_rotation_angle - vector rotation angle (positive = CCW)
%    fcor_sign - sign of fcor (if neg (i.e. SH), will switch sign of V_translation_rotation_angle)
%
% Outputs:
%    th_translation_surface - surface-modified vector angle [deg CW from N (0,360]]
%    V_translation_surface - surface-modified vector magnitude
%
% Example: 
%   [th_translation_surface V_translation_surface] = ...
%       V_translation_mod2surface(V_translation,th_translation,...
%       V_translation_reduction_factor,V_translation_rotation_angle,fcor_sign);
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
% 12 Nov 2013; Last revision: 12 Nov 2013

%------------- BEGIN CODE --------------

function [th_translation_surface,V_translation_surface] = ...
    V_translation_mod2surface(V_translation,th_translation,...
    V_translation_reduction_factor,V_translation_rotation_angle,fcor_sign)

%%Account for cyclonic in southern hemisphere
if(fcor_sign<0)
    V_translation_rotation_angle = -1*V_translation_rotation_angle;
end

V_translation_surface = V_translation_reduction_factor*V_translation;
th_translation_surface = th_translation - V_translation_rotation_angle;
if(th_translation_surface>360)
    th_translation_surface = th_translation_surface - 360;
elseif(th_translation_surface<0)
    th_translation_surface = th_translation_surface + 360;
end

assert(min(min(th_translation_surface))>=0 & max(max(th_translation_surface))<360,...
    'th_translation_surface out of bounds')
assert(min(min(V_translation_surface))>=0,'V_translation_surface out of bounds')

end

%------------- END OF CODE --------------


