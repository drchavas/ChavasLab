%distances_from_point.m -- distances from a point on grid or sphere
%Purpose: calculate distances of matrix of points from single point either
%in xy grid (arbitrary units) or WGS-84 great distances ([m]) between
%lon/lat on a sphere.
%Note: Angles on spherical triangle do not add up to 180! Crazy!
%
% Syntax:  [xx, yy, rr, th_pt1pt2, th_pt2pt1] = distances_from_point(grid_type,...
%               xx_mat,yy_mat,xx_center,yy_center)
%
% Inputs:
%   grid_type - 'lonlat' for sphere or 'xy' for regular grid
%   xx_mat [deg E [0,360)] or [distance units] - matrix of x-coordinates
%   yy_mat [deg N [-90,90]] or [distance units] - matrix of y-coordinates
%   xx_center [deg E [0,360)] or [distance units] - x-coordinate of center point
%   yy_center [deg N [-90,90]] or [distance units] - y-coordinate of center point
%
% Outputs:
%   xx [m] or [distance units] - matrix of x-distances from center point
%   yy [m] or [distance units] - matrix of meridional distances from center point
%   rr [m] or [distance units] - matrix of total distances from center point
%   th_pt1pt2 [deg CW from N [0,360)]- matrix of azimuthal angles FROM center point
%   th_pt2pt1 [deg CW from N [0,360)]- matrix of azimuthal angles TO center point
%
% Example: 
%   [xx, yy, rr, th_grid, th_tocenter] = distances_from_point(grid_type,xx_in,yy_in,...
%       TC_lon_spline,TC_lat_spline)
%
% Other m-files required: polar_qs2matlab, polar_matlab2qs
% Subfunctions: vdist
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 9 Dec 2013; Last revision: 18 Dec 2013

%------------- BEGIN CODE --------------


function [xx, yy, rr, th_pt1pt2, th_pt2pt1] = distances_from_point(grid_type,xx_mat,yy_mat,xx_center,yy_center,Earth_sphere)

switch grid_type
    case 'lonlat'   %all inputs assumed to be lat/lon
        assert(max(max(xx_mat))<360 & min(min(xx_mat)) >= 0,...
            'Input longitude grid value out of bounds')
        assert(max(max(yy_mat))<=180 & min(min(yy_mat)) > -180,...
            'Input latitude grid value out of bounds')
        assert(xx_center<360 & xx_center >= 0,...
            'Input longitude center value out of bounds')
        assert(yy_center<=180 & yy_center > -180,...
            'Input latitude center value out of bounds')
        
        %%Calculate great circle distances and angle from center in deg CCW from E (-180,180]
        %%for all lon/lat points to TC center point
        xx_center_matrix = xx_center*ones(size(yy_mat));
        yy_center_matrix = yy_center*ones(size(yy_mat));

        if(Earth_sphere==1)     %use sphere with fixed Earth's radius
            [rr, th_pt1pt2, th_pt2pt1] = vdist_sphere(yy_center_matrix,xx_center_matrix,yy_mat,xx_mat);
        else    %use WGS84 ellipsoid
            [rr, th_pt1pt2, th_pt2pt1] = vdist(yy_center_matrix,xx_center_matrix,yy_mat,xx_mat);
        end
        %distances [m]; angle deg CW from N [0,360)
        %Note: Angles on spherical triangle do not add up to 180! Crazy!
        rr = real(rr);    %occasionally calculation can give tiny imaginary part
        
        assert(isreal(th_pt1pt2),'th_pt1pt2 has imaginary part')
        clear yy_center_matrix xx_center_matrix
        
        %%For pol2cart, need angle rad CCW from E (-pi,pi]
        [th_pt1pt2_radCWfromE] = polar_qs2matlab(th_pt1pt2);
        
        assert(min(min(th_pt1pt2))>=0 & max(max(th_pt1pt2))<360,'th_pt1pt2 out of bounds')
        assert(min(min(th_pt2pt1))>=0 & max(max(th_pt2pt1))<360,'th_pt2pt1 out of bounds')
        assert(min(min(rr))>=0,'rr out of bounds')

        %%convert to cartesian grid
        [xx, yy] = pol2cart(th_pt1pt2_radCWfromE,rr);
        clear th_pt1pt2_radCWfromE

            
    case 'xy'   %simple cartesian grid
        
        xdist=xx_mat-xx_center; %[-] distance between point and TC center
        ydist=yy_mat-yy_center; %[-] distance between point and TC center
        [th_grid_radCCWfromE,rr] = cart2pol(xdist,ydist);  %[rad CW of E, (-pi,pi]]
        %convert to [deg CW of N, [0,360)]
        
        %%angle of vector from center towards grid point
        th_pt1pt2 = polar_matlab2qs(th_grid_radCCWfromE);   %th_grid_degCWfromN
        
        %%angle of vector from grid point towards center (for cartesian
        %%this is simply reverse of th_pt1pt2)
        th_pt2pt1 = th_pt1pt2 - 180;
        th_pt2pt1(th_pt2pt1<0) = th_pt2pt1(th_pt2pt1<0) + 360;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        figure(1234)
        contourf(xdist,ydist,th_grid_degCWfromN)
        colorbar
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        xx = xdist;
        yy = ydist;
        
    otherwise
        assert(1==2,'Must give grid_type "lonlat" or "xy"')
end


%------------- END OF CODE --------------