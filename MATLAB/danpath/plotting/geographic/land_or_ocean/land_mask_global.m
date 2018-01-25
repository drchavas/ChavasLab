%land_mask_global.m -- return land/ocean for input region
%Purpose: land mask: 1 = land; 0 = ocean
%
% Syntax:  [isOcean] = land_mask_global(coastal_res,make_plot)
%
% Inputs:
%   coastal_res [pts/deg] - resolution of coastline (gridpts per deg);
%       NOTE: higher resolution = more computing time (1 = decent coarse
%           res; 10 = decent high res)
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


function [isLand,lon_grid,lat_grid] = land_mask_global(coastal_res,make_plot)

% tic

%%Global coordinates
lat_minmax = [-90 90];
lon_minmax = [0 360];

switch nargin
    case 2
        coastal_res = 1;
        make_plot = 0;
    case 3
        make_plot = 0;
end

%% Load coastal data
coast = load('coast.mat');

%% Define land inside of coast
[Z, R] = vec2mtx(coast.lat, coast.long, ...
    coastal_res, lat_minmax, lon_minmax, 'filled');
    %output: 1= land, 2 = ocean; rows = latitude, columns = longitude;

%% Define corresponding latitude/longitude
numpts_lon = size(Z,2);    %columns
numpts_lat = size(Z,1);    %rows
dlon = (lon_minmax(2)-lon_minmax(1))/(numpts_lon-1);
dlat = (lat_minmax(2)-lat_minmax(1))/(numpts_lat-1);
lon_vec = lon_minmax(1):dlon:lon_minmax(2);
lat_vec = lat_minmax(1):dlat:lat_minmax(2);
lon_grid = repmat(lon_vec,numpts_lat,1);
lat_grid = repmat(lat_vec',1,numpts_lon);

%% Determine land/ocean at same set of gridpoints
val = ltln2val(Z, R, lat_grid, lon_grid);
isOcean = val == 2;
isLand = ~isOcean;
clear isOcean

assert(sum(size(isLand) == size(lon_grid))==2,'Problem with longitude grid creation')
assert(sum(size(isLand) == size(lat_grid))==2,'Problem with latitude grid creation')


%% Plot the points on geographic map
if(make_plot)
    figure
    contourf(lon_grid,lat_grid,isLand)
    colormap([0 1 0;0 0 0])
end

% toc

end

%------------- END OF CODE --------------