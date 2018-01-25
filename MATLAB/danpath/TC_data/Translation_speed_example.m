%Translation_speed_example.m 
%
% Syntax:  [u_translation,v_translation] = ...
%                Translation_speed(lat_in,lon_in,t_tot_tc_in,t_interp)
%
% Inputs:
%   lat_in - vector of track latitudes (any continuous format)
%   lon_in - vector of track longitudes (any continuous format)
%   t_tot_tc_in [day] - vector of track times (any single-number format)
%   t_interp [day] - desired interpolation time (same format as t_tot_tc_in)
%
% Outputs:
%   u_translation - x-direction translation speed [ms-1]
%   v_translation - y-direction translation speed [ms-1]
%


lat = 10:15;
lon = -80:-75;
tt = 10:15;   %whatever units you like
tt_interp = 13.3;   %same units

[u_translation,v_translation] = Translation_speed(lat,lon,tt,tt_interp)
sprintf('this is the translation velocity vector at the interpolation time.')
sprintf('u is in the zonal (west-east) direction (positive = towards the east).')
sprintf('v is in the meridional (south-north) direction (positive = towards the north).')