%%latlon2distfrompt.m
%Purpose: Input lat/lon data + center point and return gridpoints whose
%   great-circle distance from center point is below a threshold
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%   lon_in [deg E [0,360)] - input matrix of longitudes
%   lat_in [deg N [-90,90)] - input matrix of latitudes
%   lon_cent [deg E [0,360)] - input center point longitude
%   lat_cent [deg N [-90,90)] - input center point latitude
%   rdist_max [m] - only keep points with distances less than this
%
% Outputs:
%   indices_out [] - linear indices for all saved gridpoints
%   rdist_out [m] - matrix of distances (NaN if too far)
%   th_grid_out [deg CW from N [0,360)]- matrix of azimuthal angles FROM
%      center point  (NaN if too far)
%   th_tocenter_out [deg CW from N [0,360)]- matrix of azimuthal angles TO
%      center point  (NaN if too far)
%   lon_out [deg E [0,360)] - output matrix of longitudes (NaN if too far)
%   lat_out [deg N [-90,90)] - output matrix of latitudes (NaN if too far)
%   xx_out [m] - matrix of zonal distances from center point
%   yy_out [m] - matrix of meridional distances from center point
%
% Example:
%   [ii_out,jj_out,rdist_out,th_grid,lon_out,lat_out,xx_out,yy_out] = ...
%       latlon2distfrompt(lon_mat,lat_mat,lon_cent_SLP_in,...
%       lat_cent_SLP_in,rdist_max);
%
% Other m-files required: distances_from_point
% Subfunctions: none
% MAT-files required: none
%
% See also:
%
% Note: Angles on spherical triangle do not add up to 180! Crazy!

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 27 May 2014; Last revision:

%------------- BEGIN CODE --------------


function [indices_out,rdist_out,th_grid_out,th_tocenter_out,lon_out,lat_out,xx_out,yy_out] = latlon2distfrompt(lon_in,lat_in,lon_cent,lat_cent,rdist_max,Earth_sphere);

%% Calculate great circle distances from (lon_cent,lat_cent) to all points 
%%(on sphere: WGS84 great circle)
grid_type = 'lonlat';
[xx_in, yy_in, rdist_in, th_grid_in, th_tocenter_in] = distances_from_point(grid_type,lon_in,lat_in,...
    lon_cent,lat_cent,Earth_sphere);
    %for 'lonlat' grid_type: xx, yy, rr will be in distances in [m]
        
%% Keep those points with distances r<=rdist_max
[indices_out] = find(rdist_in<=rdist_max);
rdist_out = NaN(size(rdist_in));
rdist_out(indices_out) = rdist_in(indices_out);
th_grid_out = NaN(size(th_grid_in));
th_grid_out(indices_out) = th_grid_in(indices_out);
th_tocenter_out = NaN(size(th_tocenter_in));
th_tocenter_out(indices_out) = th_tocenter_in(indices_out);
lon_out = NaN(size(lon_in));
lon_out(indices_out) = lon_in(indices_out);
lat_out = NaN(size(lat_in));
lat_out(indices_out) = lat_in(indices_out);
xx_out = NaN(size(xx_in));
xx_out(indices_out) = xx_in(indices_out);
yy_out = NaN(size(yy_in));
yy_out(indices_out) = yy_in(indices_out);

%% TESTING: th_grid %%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1234);set(gcf,'Visible','on')
hp=pcolor(lon_out,lat_out,th_grid_out);set(hp, 'EdgeColor', 'none');
%pcolor(xx_out,yy_out,th_grid_out)
hold on
plot(lon_cent,lat_cent,'g.','MarkerSize',20)
%plot(0,0,'g.','MarkerSize',20)
colorbar
figure(1235);set(gcf,'Visible','on')
hp=pcolor(lon_out,lat_out,th_tocenter_out);set(hp, 'EdgeColor', 'none');
%pcolor(xx_out,yy_out,th_grid_out)
hold on
plot(lon_cent,lat_cent,'g.','MarkerSize',20)
%plot(0,0,'g.','MarkerSize',20)
colorbar
figure(1236);set(gcf,'Visible','on')
hp=pcolor(lon_out,lat_out,th_grid_out+th_tocenter_out);set(hp, 'EdgeColor', 'none');
%pcolor(xx_out,yy_out,th_grid_out)
hold on
plot(lon_cent,lat_cent,'g.','MarkerSize',20)
%plot(0,0,'g.','MarkerSize',20)
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TESTING: rdist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1236);set(gcf,'Visible','on')
contour(lon_out,lat_out,rdist_out/1000)
colorbar
hold on
plot(lon_cent,lat_cent,'m*','MarkerSize',20)
lon_rdistmax = lon_out(rdist_out==max(max(rdist_out)));
lat_rdistmax = lat_out(rdist_out==max(max(rdist_out)));
plot(lon_rdistmax,lat_rdistmax,'gx','MarkerSize',10)
xlabel('longitude')
ylabel('latitude')
title(sprintf('great circle distances [km]; max = %5.0f',max(max(rdist_out))/1000))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

%------------- END OF CODE --------------