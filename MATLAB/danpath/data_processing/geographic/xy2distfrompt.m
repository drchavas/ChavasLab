%%xy2distfrompt.m
%Purpose: Input x/y data + center point and return gridpoints whose
%   distance from center point is below a threshold
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%   xx_in [-] - input matrix of longitudes
%   yy_in [-] - input matrix of latitudes
%   x_cent [-] - input center point longitude
%   y_cent [-] - input center point latitude
%   rdist_max [-] - only keep points with distances less than this
%
% Outputs:
%   indices_out [] - linear indices for all saved gridpoints
%   rdist_out [-] - matrix of distances (NaN if too far)
%   th_grid_out [deg CW from N [0,360)]- matrix of azimuthal angles FROM
%      center point  (NaN if too far)
%   th_tocenter_out [deg CW from N [0,360)]- matrix of azimuthal angles TO
%      center point  (NaN if too far)
%   xx_out [-] - matrix of zonal distances from center point
%   yy_out [-] - matrix of meridional distances from center point
%
% Example:
%   [ii_out,jj_out,rdist_out,th_grid,xx_out,yy_out,xx_out,yy_out] = ...
%       xy2distfrompt(xx_mat,yy_mat,x_cent_SLP_in,...
%       y_cent_SLP_in,rdist_max);
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
% 17 Jun 2017; Last revision:

%------------- BEGIN CODE --------------


function [indices_out,rdist_out,th_grid_out,th_tocenter_out,xx_out,yy_out] = xy2distfrompt(xx_in,yy_in,x_cent,y_cent,rdist_max);

%% Calculate distances from (x_cent,y_cent) to all points 
grid_type = 'xy';
Earth_sphere = NaN;
[~, ~, rdist_in, th_grid_in, th_tocenter_in] = distances_from_point(grid_type,xx_in,yy_in,...
    x_cent,y_cent,Earth_sphere);
    %for 'lonlat' grid_type: xx, yy, rr will be in distances in [-]
        
%% Keep those points with distances r<=rdist_max
[indices_out] = find(rdist_in<=rdist_max);
% rdist_out = NaN(size(rdist_in));
% rdist_out(indices_out) = rdist_in(indices_out);
% th_grid_out = NaN(size(th_grid_in));
% th_grid_out(indices_out) = th_grid_in(indices_out);
% th_tocenter_out = NaN(size(th_tocenter_in));
% th_tocenter_out(indices_out) = th_tocenter_in(indices_out);
% xx_out = NaN(size(xx_in));
% xx_out(indices_out) = xx_in(indices_out);
% yy_out = NaN(size(yy_in));
% yy_out(indices_out) = yy_in(indices_out);
% xx_out = NaN(size(xx_in));

rdist_out = rdist_in(indices_out);
th_grid_out = th_grid_in(indices_out);
th_tocenter_out = th_tocenter_in(indices_out);
xx_out = xx_in(indices_out) - x_cent;
yy_out = yy_in(indices_out) - y_cent;

%% TESTING: th_grid %%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1234);set(gcf,'Visible','on')
scatter(xx_out,yy_out,2,th_grid_out)
%pcolor(xx_out,yy_out,th_grid_out)
hold on
plot(x_cent,y_cent,'g.','MarkerSize',20)
%plot(0,0,'g.','MarkerSize',20)
colorbar
figure(1235);set(gcf,'Visible','on')
scatter(xx_out,yy_out,2,th_tocenter_out)
%pcolor(xx_out,yy_out,th_grid_out)
hold on
plot(x_cent,y_cent,'g.','MarkerSize',20)
%plot(0,0,'g.','MarkerSize',20)
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TESTING: rdist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1236);set(gcf,'Visible','on')
scatter(xx_out,yy_out,2,rdist_out)
colorbar
hold on
plot(0,0,'m*','MarkerSize',20)
xx_rdistmax = xx_out(rdist_out==max(max(rdist_out)));
yy_rdistmax = yy_out(rdist_out==max(max(rdist_out)));
plot(xx_rdistmax,yy_rdistmax,'gx','MarkerSize',10)
xlabel('xdist')
ylabel('ydist')
title(sprintf('distances [-]; max = %5.0f',max(max(rdist_out))/1000))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

%------------- END OF CODE --------------