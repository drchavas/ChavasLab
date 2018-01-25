%TC_center_find_vortmax.m -- find center of circulation in 2d flow
%Purpose: find center of circulation in a given 2D (xy) rotating flow based
%on smoothed vorticity maximum
%
% Syntax:  [xx_center, yy_center] = TC_center_find_vortmax(xx,yy,...
%                uu,vv)
%
% Inputs:
%   xx [m] - matrix of x-distances
%   yy [m] - matrix of y-distances
%   uu [ms-1] - x-direction flow speeds
%   vv [ms-1] - y-direction flow speeds
%
% Outputs:
%   xx_center [m] - x-coordinate of rotating flow center
%   yy_center [m] - y-coordinate of rotating flow center
%
% Example: 
%   [xx_center, yy_center] = TC_center_find_vortmax(xx,yy,uu,vv)
%
% Other m-files required: none
% Subfunctions: filter2, gridddata
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 9 Dec 2013; Last revision: 10 Dec 2013

%------------- BEGIN CODE --------------


function [xx_center, yy_center] = TC_center_find_vortmax(xx,yy,uu,vv)

%% Parameters used below
num_smooth = 1; %number of passes by 9-pt (3x3) smoother
L_zoombox = 100*1000;  %[m] side length of box to zoom in on
highres_dxy = 1*1000;   %[m]; resolution of high-res zoomed in box

%% Calculate vertical vorticity
%%curl() function only works if grid is regular
%%vorticity = dv/dx - du/dy (positive = CCW)
zvort = NaN*uu;
dvdx = (vv(3:end,2:end-1) - vv(1:end-2,2:end-1))./(xx(3:end,2:end-1) - xx(1:end-2,2:end-1));
dudy = (uu(2:end-1,3:end) - uu(2:end-1,1:end-2))./(yy(2:end-1,3:end) - yy(2:end-1,1:end-2));
zvort_temp = dvdx - dudy;
zvort(2:end-1,2:end-1) = zvort_temp;
clear zvort_temp

%% Smooth data
h = 1/9*ones(3);
zvort_smooth = zvort;
for ii=1:num_smooth
    zvort_smooth = filter2(h,zvort_smooth);
end

%% Interpolate central region of smoothed field to high res

%%Begin with location of maximum value of vertical vorticity
[ii_zvortmax, jj_zvortmax] = find(zvort_smooth == max(max(zvort_smooth)));
%TC_xx_zvortmax = xx(ii_zvortmax,jj_zvortmax);  %FOR TESTING ONLY (BELOW)
%TC_yy_zvortmax = yy(ii_zvortmax,jj_zvortmax);  %FOR TESTING ONLY (BELOW)

%%Zoom into box centered around point of max value
dr_x_mean = nanmean(nanmean(xx(2:end,2:end-1) - xx(1:end-1,2:end-1)));
dr_y_mean = nanmean(nanmean(xx(2:end-1,2:end) - xx(2:end-1,1:end-1)));
dr_mean = sqrt(dr_x_mean.^2+dr_y_mean.^2);
L_gridpts_zoombox = ceil(L_zoombox/dr_mean);
dii = floor(L_gridpts_zoombox/2);
ii_interp = ii_zvortmax-dii:ii_zvortmax+dii;    %odd number of points always
jj_interp = jj_zvortmax-dii:jj_zvortmax+dii;

clear ii_zvortmax jj_zvortmax

xx_zoom = xx(ii_interp,jj_interp);
yy_zoom = yy(ii_interp,jj_interp);
zvort_smooth_zoom = zvort_smooth(ii_interp,jj_interp);

clear ii_interp jj_interp

%%Interpolate vertical vorticity values in zoomed box to high-resolution grid
xx_highres_temp = min(min(xx_zoom)):highres_dxy:max(max(xx_zoom));
yy_highres_temp = min(min(yy_zoom)):highres_dxy:max(max(yy_zoom));

xx_highres = repmat(xx_highres_temp',1,length(yy_highres_temp));
yy_highres = repmat(yy_highres_temp,length(xx_highres_temp),1);

clear xx_highres_temp yy_highres_temp

zvort_smooth_interp = griddata(xx_zoom(:),yy_zoom(:),zvort_smooth_zoom(:),xx_highres,yy_highres,'cubic');

%%Identify (x,y) value for point of maximum vertical vorticity
%%This is the center of circulation
[ii_temp, jj_temp] = find(zvort_smooth_interp == max(max(zvort_smooth_interp)));
xx_center = xx_highres(ii_temp,jj_temp);
yy_center = yy_highres(ii_temp,jj_temp);

clear ii_temp jj_temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING: plot comparison of raw and smoothed vertical vorticity fields
%%with locations of peak values %%%%%
%{
figure(1003)
hold off
contourf(xx_zoom,yy_zoom,zvort_smooth_zoom)
colorbar
hold on
plot(TC_xx_zvortmax,TC_yy_zvortmax,'mx','MarkerSize',20)

figure(1004)
hold off
contourf(xx_highres,yy_highres,zvort_smooth_interp)
colorbar
hold on
plot(xx_center,yy_center,'mx','MarkerSize',20)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear xx_zoom yy_zoom zvort_smooth_zoom 
clear xx_highres yy_highres zvort_smooth_interp

%------------- END OF CODE --------------