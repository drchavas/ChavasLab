%land_or_ocean.m -- return land/ocean for input points
%Purpose: determine if input points are over land or ocean
%
% Syntax:  [isOcean] = land_or_ocean(lat,lon,coastal_res)
%
% Inputs:
%   lat [deg N [-90,90]] - vector of latitude values
%   lon [deg E (-180,180]] - vector of corresponding longitude values
%   coastal_res [pts/deg] - resolution of coastline (gridpts per deg);
%       NOTE: higher resolution = more computing time (1 = decent coarse
%           res; 10 = decent high res)
%
% Outputs:
%   is_Ocean - 1 = ocean, 0 = land
%
% Example: (see land_or_ocean_example.m)
%   lat = -90:10:90;
%   lon = -180:20:180;
%   coastal_res = 5;
%   [isOcean] = land_or_ocean(lat,lon,coastal_res)
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
% 27 Jan 2014; Last revision: 19 Aug 2014

% Revision history:
% 14 May 2014 - updated min/max to allow for matrix input
% 19 Aug 2014 - fixed if statement for 180-->360 adjustment to work for
%   matrices
% 1 Jul 2015 - commented out print to screen about longitude adjustment
% 2017-04-19 - increased range of lat/lon values for search box to 3 from 2

%------------- BEGIN CODE --------------


function [isOcean] = land_or_ocean(lat,lon,coastal_res,make_plot)

% tic

switch nargin
    case 2
        coastal_res = 1;
        make_plot = 0;
    case 3
        make_plot = 0;
end

if(sum(lon(:)>180)>0)
    lon(lon>180) = lon(lon>180) - 360;  %adjust if using [0,360) lon values
%    sprintf('Adjusting lon values from [0,360) to (-180,180]')
end

%% Load coastal data
coast = load('coast.mat');

%% Define search region (want as small as possible to minimize computation)
lat_search_min = max([min(min(lat))-3 -90]); %deg N (-90,90]
lat_search_max = min([max(max(lat))+3 90]);  %deg N (-90,90]
lon_search_min = max([min(min(lon))-3 -180]);    %deg W (-180,180]
lon_search_max = min([max(max(lon))+3 180]); %deg W (-180,180]

%% Define land inside of coast
[Z, R] = vec2mtx(coast.lat, coast.long, ...
    coastal_res, [lat_search_min lat_search_max], [lon_search_min lon_search_max], 'filled');

%% Return land/ocean for each input point
val = ltln2val(Z, R, lat, lon);
isOcean = val == 2;
%isLand = ~isOcean;

%% Plot the points on geographic map
if(make_plot)
    figure; worldmap(Z, R)
    geoshow(Z, R, 'DisplayType', 'texturemap')
    colormap([0 1 0;0 0 0;0 1 0;0 0 1])
    plotm(lat(isOcean),lon(isOcean),'ro')
    plotm(lat(~isOcean),lon(~isOcean),'mx')
end

% toc

end

%------------- END OF CODE --------------