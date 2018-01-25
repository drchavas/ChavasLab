%distances_from_point_sphere.m -- distances from a point on a sphere
%Purpose: calculate distances of matrix of lat/lon points from single lat/lon
%point on sphere
%
% Syntax:  [xx, yy, rr, th_pt1pt2] = distances_from_point_sphere(lat_mat,...
%               lon_mat,lat_center,lon_center)
%
% Inputs:
%   lat_mat [deg N [-90,90]] - matrix of latitude values
%   lon_mat [deg E [0,360)] - matrix of longitude values
%   lat_center [deg N [-90,90]] - latitude of point from which distances are calculated
%   lon_center [deg E [0,360)] - longitude of point from which distances are calculated
%
% Outputs:
%   xx - matrix of zonal distances from center point [m]
%   yy - matrix of meridional distances from center point [m]
%   rr - matrix of total distances from center point [m]
%   th_pt1pt2 [deg CCW from E (-pi,pi]]- matrix of azimuthal angles from center point
%
% Example: 
%   [xx, yy, rr, th_pt1pt2] = distances_from_point_sphere(lat_mat,lon_mat,lat_center,lon_center)
%
% Other m-files required: none
% Subfunctions: vdist
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 9 Dec 2013; Last revision: 9 Dec 2013

%------------- BEGIN CODE --------------


function [xx, yy, rr, th_pt1pt2] = distances_from_point_sphere(lat_mat,lon_mat,lat_center,lon_center)

%% Calculate great circle distances and angle from center in deg CCW from E (-180,180]
%% for all lat/lon points to TC center point
lat_center_matrix = lat_center*ones(size(lat_mat));
lon_center_matrix = lon_center*ones(size(lat_mat));

[rr, th_pt1pt2_temp] = vdist(lat_center_matrix,lon_center_matrix,lat_mat,lon_mat);
    %distances [m]; angle deg CW from N [0,360)
rr = real(rr);    %occasionally calculation can give tiny imaginary part
assert(isreal(th_pt1pt2_temp),'th_pt1pt2 has imaginary part')

%%convert th_pt1pt2_temp to deg CCW from E (-pi,pi]
[th_pt1pt2] = polar_qs2matlab(th_pt1pt2_temp); 
clear th_pt1pt2_temp

clear lat_center_matrix lon_center_matrix

assert(min(min(th_pt1pt2))>-pi & max(max(th_pt1pt2))<=pi,'th_pt1pt2 out of bounds')

%%convert to cartesian grid
[xx, yy] = pol2cart(th_pt1pt2,rr);

%------------- END OF CODE --------------