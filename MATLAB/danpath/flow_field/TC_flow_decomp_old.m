%TC_flow_decomp_old.m -- complete decomposition of 2d flow field
%Purpose: Isolation of vortex flow, helmholtz decomposition, wind radii
%
% Syntax:  [flow_center, dist_from_input_center, ...
%    grid_from_center, V_helmholtz_2d, rmaxVmaxxy_2d, ...
%    rr_mean, V_helmholtz_radprof, n_r, asym_r, rmaxVmax_radprof, radii_radprof, ...
%    rr_mean_quadrant, VVazim_TC_r_quadrant, n_r_quadrant, uv_translation_mod] = ...
%    TC_flow_decomp(use_flow_center,grid_type,...
%    xx_mat,yy_mat,uu,vv,xx_center,yy_center,u_translation,v_translation,...
%    V_translation_reduction_factor,V_translation_rotation_angle,fcor_sign,...
%    ring_width,V_radii_in,n_r_minthresh)
%
% Inputs:
%   use_flow_center [] - 0 = just use input center point; 1 = determine center from flow field
%   grid_type [] -  'lonlat' = data is geographic; 'xy' = data is simple cartesian
%   xx_mat [lon or dist] - matrix of x coordinates
%   yy_mat [lat or dist] - matrix of y coordinates
%   uu [m/s] - matrix of x-direction wind speeds
%   vv [m/s] - matrix of y-direction wind speeds
%   xx_center [lon or dist] - x coordinate of circulation center
%   yy_center [lat or dist] - y coordinate of circulation center
%   u_translation [m/s] - x-direction of translation vector
%   v_translation [m/s] - y-direction of translation vector
%   V_translation_reduction_factor [-] - factor reduction of V_translation
%   V_translation_rotation_angle [deg CCW] - angle rotation of V_translation
%   fcor_sign [-] - sign of Coriolis parameter (+: cyclonic = CCW)
%   ring_width [m or dist] - bin width for radial profile calculation
%   V_radii_in [m/s] - vector of wind speeds for radii calculation
%   n_r_minthresh [] - minimum number of valid points for inclusion in
%       radial profile of FULL wind field (not by quadrant)
%
% Outputs:
%   %%Center point
%   flow_center [lat/lon or dist] - flow-based center of circulation
%       [xx_flowcenter yy_flowcenter]
%   dist_from_input_center [lat/lon or dist] - difference in center estimates
%       [xx_diff yy_diff]    
%
%   %%2d data
%   grid_from_center [m or dist] - grid relative to center point
%       {xx,yy,rr,th_grid}
%   V_helmholtz_2d [m/s] - 2d flow field decomposition
%       {VV_TC,VVazim_TC,UUrad_TC,vvazimx_TC,vvazimy_TC,uuradx_TC,uurady_TC}
%   rmaxVmaxxy_2d [r: m or dist xy: lat/lon or dist] - 2d (rmax,Vmax,x,y) data
%       {[r_TCmax_xy V_TCmax_xy xx_mat_r_TCmax_xy yy_mat_r_TCmax_xy],...
%           [razim_TCmax_xy Vazim_TCmax_xy xx_mat_razim_TCmax_xy yy_mat_razim_TCmax_xy]}
%
%   %%Radprof data
%   rr_mean [m or dist] - vector of radii for mean radial profile
%   V_helmholtz_radprof [m/s] - radial profiles of 2d flow field decomposition
%       {VV_TC_r_mean,VVazim_TC_r_mean,UUrad_TC_r_mean}
%   n_r [-] - number of datapoints in each radial bin
%   asym_r [-; deg CCW from N [0,360)] - normalized asymmetry vectors (mag
%       [0,1] and angle)
%   rmaxVmax_radprof [m or dist; m/s] - radprof (rmax,Vmax) data
%       {[r_TCmax_r V_TCmax_r],[razim_TCmax_r Vazim_TCmax_r]}
%   rmaxVmax_azimmean [m or dist; m/s] - azimuthal-mean (rmax,Vmax) point data
%       {[razim_TCmax_azimmean Vazim_TCmax_azimmean]}
%   radii_radprof [m or dist] - radii of V_radii_in from radial profiles
%       {r_TCradii_r,razim_TCradii_r}
%
%   %%Vazim_TC by quadrant
%   rr_mean_quadrant [m or dist] - vectors of radii for VVazim_TC_r_quadrant
%       {rr_mean_RF,rr_mean_RR,rr_mean_LR,rr_mean_LF} 
%   VVazim_TC_r_quadrant [m/s] - radial profiles of Vazim_TC by quadrant
%       {VVazim_TC_r_mean_RF,VVazim_TC_r_mean_RR,VVazim_TC_r_mean_LR,VVazim_TC_r_mean_LF}
%   n_r_quadrant [-] - number of datapoints in each radial bin for VVazim_TC_r_quadrant
%       {n_r_RF,n_r_RR,n_r_LR,n_r_LF}
%
%   %%Modified translation vector
%   uv_translation_mod [m/s] - modified translation vector components
%
% Example: 
%   [qs_flow_center, qs_dist_from_BT, ...
%     qs_grid_from_center, qs_V_helmholtz_2d, qs_rmaxVmaxxy_2d, ...
%     qs_rr_mean, qs_V_helmholtz_radprof, qs_n_r, qs_asym_r, qs_rmaxVmax_radprof, qs_radii_radprof, ...
%     qs_rr_mean_quadrant, qs_VVazim_TC_r_quadrant, qs_n_r_quadrant, qs_uv_translation_mod] = ...
%     TC_flow_decomp(use_flow_center,grid_type,...
%     qs_lon,qs_lat,qs_u,qs_v,TC_lon_spline,TC_lat_spline,u_translation,v_translation,...
%     V_translation_reduction_factor,V_translation_rotation_angle,sign(fcor),...
%     ring_width,V_radii_in);
%
% Other m-files required: polar_qs2matlab, polar_matlab2qs, ...
%   V_translation_mod2surface, distances_from_point, TC_center_find_vortmax...
%   helmholtz_decomp, radprof, TC_radprof_quadrants, TC_radprof_radii
% Subfunctions: vdist
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 17 Dec 2013; Last revision: 27 May 2014

% Revision history:
%  29 Jan 2014: added input of n_r_minthresh passed on to radprof()
%  3 Feb 2014: added output variable uv_translation_mod 
%  27 May 2014: optimized for-loop code to skip NaNs

%------------- BEGIN CODE --------------

function [flow_center, dist_from_input_center, ...
    grid_from_center, V_helmholtz_2d, rmaxVmaxxy_2d, ...
    rr_mean, V_helmholtz_radprof, n_r, asym_r, rmaxVmax_radprof, rmaxVmax_azimmean, radii_radprof, ...
    rr_mean_quadrant, VVazim_TC_r_quadrant, n_r_quadrant, uv_translation_mod] = ...
    TC_flow_decomp_old(use_flow_center,grid_type,...
    xx_mat,yy_mat,uu,vv,xx_center,yy_center,u_translation,v_translation,...
    V_translation_reduction_factor,V_translation_rotation_angle,fcor_sign,...
    ring_width,V_radii_in,n_r_minthresh)
        
%% Removed modified translation vector from flow
%%Convert vector to polar coordinates: [0,360) deg CW from N
if(~isnan(u_translation))
    [th_temp, V_translation] = cart2pol(u_translation,v_translation); %theta deg CCW from E (-pi,pi]
    [th_translation] = polar_matlab2qs(th_temp);   %theta deg CW from N [0,360)
    clear th_temp

    %%Calculate modified translation vector
    [th_translation_mod,V_translation_mod] = ...
        V_translation_mod2surface(V_translation,th_translation,...
        V_translation_reduction_factor,V_translation_rotation_angle,yy_center);

    %%Convert to cartesian
    [th_temp] = polar_qs2matlab(th_translation_mod);   %theta deg CW from N [0,360)
    [u_translation_mod, v_translation_mod] = ...
        pol2cart(th_temp,V_translation_mod); %theta deg CCW from E (-pi,pi]

    assert(abs(V_translation_reduction_factor*V_translation-V_translation_mod)<...
        abs(V_translation_mod)/10^6,'Translation reduction calculation not working')

else
    u_translation_mod = 0;
    v_translation_mod = 0;
    th_translation = 0;
end

%%Subtract surface translation vector from wind field
uu_TC = uu - u_translation_mod;
vv_TC = vv - v_translation_mod;

%%Calculate flow magnitudes
VV_TC = sqrt(uu_TC.^2 + vv_TC.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING: above steps -- these should make sense given the translation speed %%%%%
%{
VV = sqrt(uu.^2 + vv.^2);
figure(1011)
contourf(xx_mat',yy_mat',VV');
colorbar
title(sprintf('max = %3.1f',max(max(VV))))
figure(1012)
contourf(xx_mat',yy_mat',VV_TC');
title(sprintf('max = %3.1f',max(max(VV_TC))))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING: plot vector flow field %%%%%
%{
hh=figure(1002);
hold off
clf(hh)
set(hh,'Units','centimeters');
hpos = [0 0 30 30];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
quiver(xx_mat',yy_mat',uu_TC',vv_TC');
hold on
plot(xx_center,yy_center,'mx','MarkerSize',20)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate grid of distances between points relative to
%%  input center point (on sphere: WGS84 great circle)
[xx, yy, rr, th_grid] = distances_from_point(grid_type,xx_mat,yy_mat,...
    xx_center,yy_center);
    %for 'lonlat' grid_type: xx, yy, rr will be in distances in [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING great circle distance, angle calculation, cartesian conversion %%
%{
figure(1000)
subplot(2,1,1)
pcolor(xx_mat,yy_mat,rr)
%pcolor(xx_mat,yy_mat,xx/1000)
%pcolor(xx_mat,yy_mat,yy/1000)
%pcolor(xx_mat,yy_mat,(sqrt(xx.^2+yy.^2)/1000))
    %coastline is visible because QuikSCAT gridpoint spacing is irregular
colorbar
subplot(2,1,2)
pcolor(xx_mat,yy_mat,th_grid)
    %coastline is visible because QuikSCAT gridpoint spacing is irregular
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine desired center of circulation and update distances %%%%%%%%
switch use_flow_center
    case 0
        %use input flow center point
        xx_flowcenter = NaN;
        yy_flowcenter = NaN;
        xx_diff = NaN;
        yy_diff = NaN;
    case 1
        %%Find center of circulation from flow field directly %%%%
        %%%STEP 2: Find location of maximum in smoothed vertical vorticity field
        [xx_flowcenter, yy_flowcenter] = TC_center_find_vortmax(xx,yy,uu,vv);

        %%Map center location back to grid (only useful for lon/lat)
        ii_good_temp = ~isnan(yy_mat);
        [xx_flowcenter_temp] = griddata(xx(ii_good_temp),yy(ii_good_temp),...
            xx_mat(ii_good_temp),xx_flowcenter,yy_flowcenter,'cubic');
        [yy_flowcenter_temp] = griddata(xx(ii_good_temp),yy(ii_good_temp),...
            yy_mat(ii_good_temp),xx_flowcenter,yy_flowcenter,'cubic');                
        xx_flowcenter = xx_flowcenter_temp;
        yy_flowcenter = yy_flowcenter_temp;
        clear xx_flowcenter_temp yy_flowcenter_temp

        %%Recalculate grid of distances between points relative to
        %%new flow-center point (on sphere: WGS84 great circle)
        [xx, yy, rr, th_grid] = distances_from_point('lonlat',xx_mat,yy_mat,...
            xx_flowcenter,yy_flowcenter);
        
        %%Record distance between two center estimates (in input units)
        xx_diff = xx_flowcenter-xx_center;
        yy_diff = yy_flowcenter-yy_center;
        %         assert(abs(TC_yy_diff)<=TC_cent_deg_max_diff,...
        %             sprintf('QuikSCAT-estimated TC lat is %3.2f deg off from Best Track',...
        %             TC_yy_diff))
        %         assert(abs(yy_flowcenter-TC_lat)<=TC_cent_deg_max_diff,...
        %             sprintf('QuikSCAT-estimated TC lon is %3.2f deg off from Best Track',...
        %             TC_xx_diff))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING: plot vector wind field with estimated center and compare to
%%original Best Track-derived center%%%%%
%{
figure(1005)
hold off
quiver(xx_mat,yy_mat,uu_TC,vv_TC)
hold on
plot(xx_flowcenter,yy_flowcenter,'mx','MarkerSize',20)
plot(xx_center,yy_center,'m*','MarkerSize',20)
title('Input wind vectors, input center (*), flow center (X)','FontSize',12)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Helmholtz decomposition into azimuthal/radial flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vvazimx_TC,vvazimy_TC,VVazim_TC,uuradx_TC,...
    uurady_TC,UUrad_TC] = helmholtz_decomp(th_grid,uu_TC,vv_TC,fcor_sign);
    %VVazim_TC: signed azimuthal flow (+ = cyclonic)
    %UUrad_TC: signed radial flow (+ = outward)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RADIAL PROFILES AND RADII %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate radial profile of interesting things
indices_good = ~isnan(VV_TC);
rr_good = rr(indices_good);
th_grid_good = th_grid(indices_good);
VVazim_TC_good = VVazim_TC(indices_good);

[rr_mean, VV_TC_r_mean, n_r, asym_mag_r, asym_th_r] = radprof(rr_good,VV_TC(indices_good),ring_width,n_r_minthresh,th_grid_good);
[~, VVazim_TC_r_mean, ~] = radprof(rr_good,VVazim_TC_good,ring_width,n_r_minthresh);
[~, UUrad_TC_r_mean, ~] = radprof(rr_good,UUrad_TC(indices_good),ring_width,n_r_minthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Breakdown radial profile by quadrant
[~, ~, ~, rr_mean_RF, VVazim_TC_r_mean_RF, n_r_RF, rr_mean_LF, ...
    VVazim_TC_r_mean_LF, n_r_LF, rr_mean_RR, VVazim_TC_r_mean_RR, n_r_RR, rr_mean_LR, ...
    VVazim_TC_r_mean_LR, n_r_LR] = TC_radprof_quadrants(rr_good,...
    th_grid_good,VVazim_TC_good,th_translation,ring_width);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract values of interest from 2D field
%%Vmax -- point-estimate from 2d data
V_TCmax_xy = max(max(VV_TC));
Vazim_TCmax_xy = max(max(VVazim_TC));

%%r_TCmax_xy -- point-estimate from 2d data
[i_r_TCmax_xy, j_r_TCmax_xy] = find(VV_TC==V_TCmax_xy);
[i_razim_TCmax_xy, j_razim_TCmax_xy] = find(VVazim_TC==Vazim_TCmax_xy);
r_TCmax_xy = rr(i_r_TCmax_xy,j_r_TCmax_xy);
razim_TCmax_xy = rr(i_razim_TCmax_xy,j_razim_TCmax_xy);
xx_mat_r_TCmax_xy = xx_mat(i_r_TCmax_xy,j_r_TCmax_xy);
yy_mat_r_TCmax_xy = yy_mat(i_r_TCmax_xy,j_r_TCmax_xy);
xx_mat_razim_TCmax_xy = xx_mat(i_razim_TCmax_xy,j_razim_TCmax_xy);
yy_mat_razim_TCmax_xy = yy_mat(i_razim_TCmax_xy,j_razim_TCmax_xy);
clear i_r_TCmax_xy j_r_TCmax_xy i_razim_TCmax_xy j_razim_TCmax_xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract radii of interest from radial profiles
%%Parameters for the radial wind profile of V_TC
[V_TCmax_r, r_TCmax_r, r_TCradii_r] = TC_radprof_radii(rr_mean,VV_TC_r_mean,V_radii_in);

%%Parameters for the radial wind profile of Vazim_TC
[Vazim_TCmax_r, razim_TCmax_r, razim_TCradii_r] = TC_radprof_radii(rr_mean,VVazim_TC_r_mean,V_radii_in);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate azimuthal mean (rmax,Vmax) (rather than (rmax,Vmax) from the azimuthal mean V(r))
num_azim_bins = 8;
[razim_TCmax_azimmean,Vazim_TCmax_azimmean] = ...
    rmaxVmax_points_azimuthalmean(rr_good,VVazim_TC_good,th_grid_good,num_azim_bins);

%% Cleaning up the output %%%%%%%%%%%%%%%%%%
%%Center point
flow_center = [xx_flowcenter yy_flowcenter];
dist_from_input_center = [xx_diff yy_diff];

%%2d grid data relative to center point
grid_from_center{1} = xx;
grid_from_center{2} = yy;
grid_from_center{3} = rr;
grid_from_center{4} = th_grid;

%%2d decomposed flow fields
V_helmholtz_2d{1} = VV_TC;
V_helmholtz_2d{2} = VVazim_TC;
V_helmholtz_2d{3} = UUrad_TC;
V_helmholtz_2d{4} = vvazimx_TC;
V_helmholtz_2d{5} = vvazimy_TC;
V_helmholtz_2d{6} = uuradx_TC;
V_helmholtz_2d{7} = uurady_TC;

%%2d (rmax,Vmax, xmax, ymax)
rmaxVmaxxy_2d{1} = [r_TCmax_xy V_TCmax_xy xx_mat_r_TCmax_xy yy_mat_r_TCmax_xy];
rmaxVmaxxy_2d{2} = [razim_TCmax_xy Vazim_TCmax_xy xx_mat_razim_TCmax_xy yy_mat_razim_TCmax_xy];

%%Radial profiles of decomposed flow fields
V_helmholtz_radprof{1} = VV_TC_r_mean;
V_helmholtz_radprof{2} = VVazim_TC_r_mean;
V_helmholtz_radprof{3} = UUrad_TC_r_mean;

%%Radprof (rmax,Vmax)
rmaxVmax_radprof{1} = [r_TCmax_r V_TCmax_r];
rmaxVmax_radprof{2} = [razim_TCmax_r Vazim_TCmax_r];

%%Azimuthal-mean (rmax,Vmax)
rmaxVmax_azimmean = [razim_TCmax_azimmean Vazim_TCmax_azimmean];

%%Radii
radii_radprof{1} = r_TCradii_r;
radii_radprof{2} = razim_TCradii_r;

%%Breakdown of VVazim_TC_r by quadrant
%%Radii by quadrant
rr_mean_quadrant{1} = rr_mean_RF;
rr_mean_quadrant{2} = rr_mean_RR;
rr_mean_quadrant{3} = rr_mean_LR;
rr_mean_quadrant{4} = rr_mean_LF;

%%Vazim_TC by quadrant
VVazim_TC_r_quadrant{1} = VVazim_TC_r_mean_RF;
VVazim_TC_r_quadrant{2} = VVazim_TC_r_mean_RR;
VVazim_TC_r_quadrant{3} = VVazim_TC_r_mean_LR;
VVazim_TC_r_quadrant{4} = VVazim_TC_r_mean_LF;

%%n_r by quadrant
n_r_quadrant{1} = n_r_RF;
n_r_quadrant{2} = n_r_RR;
n_r_quadrant{3} = n_r_LR;
n_r_quadrant{4} = n_r_LF;

%%asymmetry vector
asym_r{1} = asym_mag_r;
asym_r{2} = asym_th_r;

%%modified translation vector
uv_translation_mod = [u_translation_mod v_translation_mod];

%------------- END OF CODE --------------