%Translation_speed.m -- calculate TC translation velocity from track data
%Purpose: calculate TC translation velocity at input time/location from
%input track data
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
% Example: 
%	[u_translation v_translation] = Translation_speed(Lat_BT_TC,Lon_BT_TC,...
%       tt_daysince1858111700UTC_TC,qs_tt_daysince1858111700UTC);
%
% Other m-files required: none
% Subfunctions: spline
% MAT-files required: none
%
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 4 Dec 2013; Last revision: 23 Oct 2014

% Revision history:
% 13 Jan 2014: changed spline() to pchip()
% 23 Oct 2014: vectorized (a single period added to u_translation)
% 06 Apr 2017: switched pchip to interp1(...,'pchip',NaN) to avoid
%   extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [u_translation,v_translation] = Translation_speed(lat_in,lon_in,t_tot_tc_in,t_interp)

%% Interpolate translation speeds to desired time

%set any constants
km_per_deg=111.325; %[km/deg]

%%Lat (y) NOTE: MUST DO LAT FIRST!
t = t_tot_tc_in; y = lat_in;
lat_TC = interp1(t,y,t_interp,'pchip',NaN);   %latitude at interpolated point is needed to calculate u later!
pp = pchip(t,y);   %piecewise cubic polynomial is created; (one between each pair of consecutive points)

if(~isnan(lat_TC))

    %redefine pp to represent the DERIVATIVE of the actual spline
    A=zeros(size(pp.coefs));    %pp.coefs gives matrix of coefficients
    if(size(A,2)==4)  %>=4 input best track data points (cubic)
        A(:,2)=3*pp.coefs(:,1);
        A(:,3)=2*pp.coefs(:,2);
        A(:,4)=pp.coefs(:,3);
    elseif(size(A,2)==3)  %only 3 input best track data points (can only do quadratic)
        A(:,2)=2*pp.coefs(:,1);
        A(:,3)=pp.coefs(:,2);
    elseif(size(A,2)==2)  %only 2 input best track data points (linear)
        A(:,2)=pp.coefs(:,1);
    end     %only 1 input best track data point (constant)
    pp.coefs=A; %replace coefs --> now pp represents the spline of the DERIVATIVE

    %evaluate the derivative at the desired time (d(lat)/dt)
    dlat_dt=ppval(pp,t_interp);

    %calculate meridional speed [m/s]: v=(d(lat)/dt) x (111.325 km/deg)
    v_translation=dlat_dt*km_per_deg;   %[km/hr]
    v_translation=v_translation*1000/(24*60*60);   %[m/s]

    clear pp A dlat_dt

    %%Lon (x)
    x = lon_in;
    pp = spline(t,x);   %piecewise cubic polynomial is created; (one between each pair of consecutive points)

    %redefine pp to represent the DERIVATIVE of the actual spline
    A=zeros(size(pp.coefs));    %pp.coefs gives matrix of coefficients
    if(size(A,2)==4)  %cubic
        A(:,2)=3*pp.coefs(:,1);
        A(:,3)=2*pp.coefs(:,2);
        A(:,4)=pp.coefs(:,3);
    elseif(size(A,2)==3)  %quadratic
        A(:,2)=2*pp.coefs(:,1);
        A(:,3)=pp.coefs(:,2);
    elseif(size(A,2)==2)  %linear
        A(:,2)=pp.coefs(:,1);
    end     %constant --> A=0 (deriv = 0)
    pp.coefs=A; %replace coefs --> now pp represents the spline of the DERIVATIVE

    %evaluate the derivative at the desired time (d(lon)/dt)
    dlon_dt=ppval(pp,t_interp);

    %calculate zonal speed [m/s]: u=(d(lon)/dt) x (111.325 km/deg)*(cos(lat))
    u_translation=dlon_dt.*km_per_deg.*cosd(lat_TC);
    u_translation=u_translation*1000/(24*60*60);   %[m/s]
else
    u_translation=NaN;
    v_translation=NaN;
end

end