%TC_flow_decomp_fast.m -- faster decomposition of 2d flow field
%Purpose: Same as TC_flow_decomp, but assumes xy data as input, skips
%   quadrant calculations
%
% Syntax:  [V_helmholtz_2d, rmaxVmaxxy_2d, rr_mean, V_helmholtz_radprof,...
%    n_r, asym_r, rmaxVmax_radprof, rmaxVmax_azimmean, radii_radprof, ...
%    uv_translation_mod] = ...
%    TC_flow_decomp_fast(xx,yy,rr,th_tocenter,uu,vv,u_translation,v_translation,...
%    V_translation_reduction_factor,V_translation_rotation_angle,fcor_sign,...
%    ring_width,V_radii_in,n_r_minthresh)
%
% Inputs:
%   xx [dist or lon] - matrix of x coordinates
%   yy [dist or lat] - matrix of y coordinates
%   rr [dist] - matrix of great circle distances from center
%   th_tocenter [deg CW from N [0,360)] - angle from gridpoints TO TC center
%   th_grid [deg CW from N [0,360)] - angle from TC center to gridpoints
%   uu [m/s] - matrix of x-direction wind speeds
%   vv [m/s] - matrix of y-direction wind speeds
%   u_translation [m/s] - x-direction of translation vector
%   v_translation [m/s] - y-direction of translation vector
%   V_translation_reduction_factor [-] - factor reduction of V_translation
%   V_translation_rotation_angle [deg CCW] - angle rotation of V_translation
%   fcor_sign [-] - sign of Coriolis parameter (+: cyclonic = CCW)
%   ring_width [m or dist] - bin width for radial profile calculation
%   V_radii_in [m/s] - vector of wind speeds for radii calculation
%   n_r_minthresh [] - minimum number of valid points for inclusion in
%       radial profile of FULL wind field (not by quadrant)
%   account_for_azimasym [1/0] - indicator of whether to use pure azimuthal
%       mean or rebin first

%
% Outputs:
%   %%2d data
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
%   %%Modified translation vector
%   uv_translation_mod [m/s] - modified translation vector components
%
% Example: 
%    [V_helmholtz_2d, rmaxVmaxxy_2d, ...
%    rr_mean, V_helmholtz_radprof, n_r, asym_r, ...
%    rmaxVmax_radprof, rmaxVmax_azimmean, radii_radprof, ~] = ...
%    TC_flow_decomp_fast(xx_out,yy_out,rdist_out,th_tocenter_out,ubot_out,...
%    vbot_out,u_translation,v_translation,...
%    V_translation_reduction_factor,V_translation_rotation_angle,sign(fcor),...
%    ring_width,V_radii_in,n_r_minthresh);
%
% Other m-files required: polar_qs2matlab, polar_matlab2qs, ...
%   V_translation_mod2surface, TC_center_find_vortmax...
%   helmholtz_decomp, radprof, TC_radprof_quadrants, TC_radprof_radii
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 27 May 2014; Last revision: 1 Aug 2017

% Revision history:
% 2 Jun 2014 - added restriction on data for 2D (rmax,Vmax) extraction to
%   only look at data for r < 3*rmax_azimmean
% 4 Sep 2015 - cleaned up code, removed unused input parameter
% 20 Jul 2016 - added code to retain only smallest r_TCmax_r value in case
%   there are multiple
% 1 Aug 2017 - added code to retain arbitrary value of r_TCmax_xy and
%   razim_TCmax_xy

%------------- BEGIN CODE --------------

function [V_helmholtz_2d, rmaxVmaxxy_2d, rr_mean, V_helmholtz_radprof,...
    n_r, asym_r, rmaxVmax_radprof, rmaxVmax_azimmean, radii_radprof, ...
    uv_translation_mod] = ...
    TC_flow_decomp_fast(xx,yy,rr,th_tocenter,th_grid,uu,vv,u_translation,v_translation,...
    V_translation_reduction_factor,V_translation_rotation_angle,fcor_sign,...
    ring_width,V_radii_in,n_r_minthresh,account_for_azimasym)
        
%% Removed modified translation vector from flow
%%Convert vector to polar coordinates: [0,360) deg CW from N
if(~isnan(u_translation) && V_translation_reduction_factor > 0)
    [th_temp, V_translation] = cart2pol(u_translation,v_translation); %theta deg CCW from E (-pi,pi]
    [th_translation] = polar_matlab2qs(th_temp);   %theta deg CW from N [0,360)
    clear th_temp

    %%Calculate modified translation vector
    [th_translation_mod,V_translation_mod] = ...
        V_translation_mod2surface(V_translation,th_translation,...
        V_translation_reduction_factor,V_translation_rotation_angle,fcor_sign);

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
figure(1011);set(gcf,'Visible','on')
contourf(xx',yy',VV');
colorbar
title(sprintf('max = %3.1f',max(max(VV))))
figure(1012)
contourf(xx',yy',VV_TC');
title(sprintf('max = %3.1f',max(max(VV_TC))))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING: plot vector flow field %%%%%
%{
hh=figure(1002);set(gcf,'Visible','on')
hold off
clf(hh)
set(hh,'Units','centimeters');
hpos = [0 0 30 30];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
quiver(xx',yy',uu_TC',vv_TC');
hold on
plot(0,0,'g.','MarkerSize',20)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Helmholtz decomposition into azimuthal/radial flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Define angle of unit vector pointing radially away from TC center (along which u_r>0)
th_awayfromcenter = th_tocenter-180;
th_awayfromcenter(th_awayfromcenter<0) = ...
    th_awayfromcenter(th_awayfromcenter<0)+360;    %[deg CW from N [0,360)]

%% TESTING: th_awayfromcenter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(5678);set(gcf,'Visible','on')
hp=pcolor(xx,yy,th_tocenter);set(hp, 'EdgeColor', 'none');
hold on
plot(0,0,'g.','MarkerSize',20)
colorbar
figure(5679);set(gcf,'Visible','on')
hp=pcolor(xx,yy,th_awayfromcenter);set(hp, 'EdgeColor', 'none');
hold on
plot(0,0,'g.','MarkerSize',20)
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Decompose flow onto radial and azimuthal directions relative to TC center
[vvazimx_TC,vvazimy_TC,VVazim_TC,uuradx_TC,...
    uurady_TC,UUrad_TC] = helmholtz_decomp(th_awayfromcenter,uu_TC,vv_TC,fcor_sign);
    %VVazim_TC: signed azimuthal flow (+ = cyclonic)
    %UUrad_TC: signed radial flow (+ = outward)

%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(2003)
clf(2003)
scatter(xx/1000,yy/1000,2,sqrt(uu_TC.^2+vv_TC.^2))
hold on
plot(0,0,'mx')
title('storm-centered wind speed')
saveas(gcf,sprintf('4wspd_origcode.jpg'),'jpeg')

figure(2004)
clf(2004)
quiver(xx,yy,uu_TC,vv_TC);
hold on
plot(0,0,'mx')
title('storm-centered wind vector')
saveas(gcf,sprintf('5wndvec_origcode.jpg'),'jpeg')

figure(2005)
clf(2005)
scatter(xx/1000,yy/1000,10,rr/1000)
colorbar
hold on
plot(0,0,'mx')
title('distance from storm center')
saveas(gcf,sprintf('6distance_origcode.jpg'),'jpeg')

figure(2006)
clf(2006)
scatter(xx/1000,yy/1000,2,UUrad_TC)
colorbar
hold on
plot(0,0,'mx')
title('signed radial wind speed (pos=inward)')
saveas(gcf,sprintf('7urad_origcode.jpg'),'jpeg')


figure(2007)
clf(2007)
scatter(xx/1000,yy/1000,2,VVazim_TC)
colorbar
hold on
plot(0,0,'mx')
title('signed azimuthal wind speed (pos=CCW)')
saveas(gcf,sprintf('8vazim_origcode.jpg'),'jpeg')
'hi'
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    
%%FOR TESTING (NEED TO PASS IN GRID DATA (either xy or latlon))
% lat_cent_SLP_in = 78.7500
% lon_cent_SLP_in = 299.5312
% [vvazimx_TC,vvazimy_TC,VVazim_TC,uuradx_TC,...
%     uurady_TC,UUrad_TC] = helmholtz_decomp(th_awayfromcenter,uu_TC,vv_TC,fcor_sign,xx,yy,lon_cent_SLP_in,lat_cent_SLP_in);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RADIAL PROFILES AND RADII %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use good data
indices_good = ~isnan(VV_TC);
rr_good = rr(indices_good);
th_awayfromcenter_good = th_awayfromcenter(indices_good);
th_grid_good = th_grid(indices_good);
VVazim_TC_good = VVazim_TC(indices_good);
VV_TC_good = VV_TC(indices_good);

%% Calculate radial profile of interesting things
[rr_mean, VV_TC_r_mean, n_r, asym_mag_r, asym_th_r] = radprof(rr_good,VV_TC_good,ring_width,n_r_minthresh,th_grid_good,account_for_azimasym);
[~, VVazim_TC_r_mean, ~] = radprof(rr_good,VVazim_TC_good,ring_width,n_r_minthresh,th_grid_good,account_for_azimasym);
[~, UUrad_TC_r_mean, ~] = radprof(rr_good,UUrad_TC(indices_good),ring_width,n_r_minthresh,th_grid_good,account_for_azimasym);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Breakdown radial profile by quadrant
% [~, ~, ~, rr_mean_RF, VVazim_TC_r_mean_RF, n_r_RF, rr_mean_LF, ...
%     VVazim_TC_r_mean_LF, n_r_LF, rr_mean_RR, VVazim_TC_r_mean_RR, n_r_RR, rr_mean_LR, ...
%     VVazim_TC_r_mean_LR, n_r_LR] = TC_radprof_quadrants(rr_good,...
%     th_awayfromcenter_good,VVazim_TC_good,th_translation,ring_width);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Extract radii of interest from radial profiles
%%Parameters for the radial wind profile of V_TC
[V_TCmax_r, r_TCmax_r, r_TCradii_r] = TC_radprof_radii(rr_mean,VV_TC_r_mean,V_radii_in);
if(length(r_TCmax_r)>1)
    r_TCmax_r = min(r_TCmax_r); %in case there are multiple, keep both
    sprintf('WARNING MORE THAN ONE r_TCmax_r VALUE; KEEPING SMALLEST VALUE')
end
ii_inner_core = rr_good < 3*r_TCmax_r;  %only those points at sufficiently small radii, to avoid 2nd storm in vicinity
rr_good_inner_core = rr_good(ii_inner_core);
VV_TC_good_inner_core = VV_TC_good(ii_inner_core);
VVazim_TC_good_inner_core = VVazim_TC_good(ii_inner_core);
th_awayfromcenter_good_inner_core = th_awayfromcenter_good(ii_inner_core);

%%Parameters for the radial wind profile of Vazim_TC
[Vazim_TCmax_r, razim_TCmax_r, razim_TCradii_r] = TC_radprof_radii(rr_mean,VVazim_TC_r_mean,V_radii_in);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate azimuthal mean (rmax,Vmax) (rather than (rmax,Vmax) from the azimuthal mean V(r))
num_azim_bins = 8;
[razim_TCmax_azimmean,Vazim_TCmax_azimmean] = ...
    rmaxVmax_points_azimuthalmean(rr_good_inner_core,VVazim_TC_good_inner_core,th_awayfromcenter_good_inner_core,num_azim_bins);

%% Extract values of interest from 2D field
%%Vmax -- point-estimate from 2d data
V_TCmax_xy = max(max(VV_TC_good_inner_core));
Vazim_TCmax_xy = max(max(VVazim_TC_good_inner_core));

%%r_TCmax_xy -- point-estimate from 2d data
[i_r_TCmax_xy, j_r_TCmax_xy] = find(VV_TC==V_TCmax_xy);
if(length(i_r_TCmax_xy)>1)
    i_r_TCmax_xy = i_r_TCmax_xy(1); %in case there are multiple, keep both
    j_r_TCmax_xy = j_r_TCmax_xy(1); %in case there are multiple, keep both
    sprintf('WARNING MORE THAN ONE r_TCmax_xy VALUE; KEEPING ARBITRARY ONE')
end
[i_razim_TCmax_xy, j_razim_TCmax_xy] = find(VVazim_TC==Vazim_TCmax_xy);
if(length(i_razim_TCmax_xy)>1)
    i_razim_TCmax_xy = i_razim_TCmax_xy(1); %in case there are multiple, keep both
    j_razim_TCmax_xy = j_razim_TCmax_xy(1); %in case there are multiple, keep both
    sprintf('WARNING MORE THAN ONE razim_TCmax_xy VALUE; KEEPING ARBITRARY ONE')
end
r_TCmax_xy = rr(i_r_TCmax_xy,j_r_TCmax_xy);
razim_TCmax_xy = rr(i_razim_TCmax_xy,j_razim_TCmax_xy);
xx_r_TCmax_xy = xx(i_r_TCmax_xy,j_r_TCmax_xy);
yy_r_TCmax_xy = yy(i_r_TCmax_xy,j_r_TCmax_xy);
xx_razim_TCmax_xy = xx(i_razim_TCmax_xy,j_razim_TCmax_xy);
yy_razim_TCmax_xy = yy(i_razim_TCmax_xy,j_razim_TCmax_xy);
clear i_r_TCmax_xy j_r_TCmax_xy i_razim_TCmax_xy j_razim_TCmax_xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleaning up the output %%%%%%%%%%%%%%%%%%
%%2d decomposed flow fields
V_helmholtz_2d{1} = VV_TC;
V_helmholtz_2d{2} = VVazim_TC;
V_helmholtz_2d{3} = UUrad_TC;
V_helmholtz_2d{4} = vvazimx_TC;
V_helmholtz_2d{5} = vvazimy_TC;
V_helmholtz_2d{6} = uuradx_TC;
V_helmholtz_2d{7} = uurady_TC;

%%2d (rmax,Vmax, xmax, ymax)
rmaxVmaxxy_2d{1} = [r_TCmax_xy V_TCmax_xy xx_r_TCmax_xy yy_r_TCmax_xy];
rmaxVmaxxy_2d{2} = [razim_TCmax_xy Vazim_TCmax_xy xx_razim_TCmax_xy yy_razim_TCmax_xy];

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
%%asymmetry vector
asym_r{1} = asym_mag_r;
asym_r{2} = asym_th_r;

%%modified translation vector
uv_translation_mod = [u_translation_mod v_translation_mod];

%------------- END OF CODE --------------