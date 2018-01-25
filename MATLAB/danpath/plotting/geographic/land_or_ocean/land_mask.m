%land_mask.m -- return land mask for desired region
%Purpose: land mask: 1 = land; 0 = ocean
%
% Syntax:  [isOcean] = land_mask(lat_minmax,lon_minmax)
%
% Inputs:
%   lat_minmax [deg N [-90,90]] - bounding latitude values
%   lon_minmax [deg E (-180,180] or [0,360)] - bounding longitude values
%
% Outputs:
%   isLand - 1 = land, 0 = ocean
%   lon_grid - corresponding grid of longitudes
%   lat_grid - corresponding grid of latitudes
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: coast.mat (included in MATLAB)

% Original code: http://www.mathworks.com/matlabcentral/answers/1065-determining
%   -whether-a-point-on-earth-given-latitude-and-longitude-is-on-land-or-ocean

% Author: Dan Chavas (adapted from Brett Shoelson 9 Feb 2011 -- see above)
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: http://www.princeton.edu/~dchavas/
% 5 Mar 2014; Last revision: 5 Mar 2014

%------------- BEGIN CODE --------------


function [isLand,lon_grid,lat_grid] = land_mask(lat_minmax,lon_minmax)

load('landmaskdat.mat') %global land mask data, 10 pts/deg
    %vars: isLand,lon_grid,lat_grid
%%Zoom into region of interest
j_zoom_0 = find(lon_grid(1,:) >= lon_minmax(1),1)-1;
j_zoom_f = find(lon_grid(1,:) >= lon_minmax(2),1);
i_zoom_0 = find(lat_grid(:,1) >= lat_minmax(1),1)-1;
i_zoom_f = find(lat_grid(:,1) >= lat_minmax(2),1);
lon_grid = lon_grid(i_zoom_0:i_zoom_f,j_zoom_0:j_zoom_f);
lat_grid = lat_grid(i_zoom_0:i_zoom_f,j_zoom_0:j_zoom_f);
isLand = isLand(i_zoom_0:i_zoom_f,j_zoom_0:j_zoom_f);
% toc

end

%------------- END OF CODE --------------