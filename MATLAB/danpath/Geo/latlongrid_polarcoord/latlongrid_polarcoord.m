%latlongrid_polarcoord.m
%Purpose: Return lat/lons for a polar grid
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 2015-10-15; Last revision:

%------------- BEGIN CODE --------------

function [dlon_out,dlat_out,dist_out] = latlongrid_polarcoord(r_in,r_out,dr,dth,lat_cent)

%% Calculate non-dimensional entropy deficit

%%Make polar grid of points for inner disk and outer annulus
% dr = 50*1000;   %[m]
% dth = pi/4; %[rad]
% r_in = dr; %hypothetical inner-core storm radius (m)
% r_out = 100*1000; %hypothetical inner-core storm radius (m)

if(~exist('lat_cent'))
    lat_cent = 0;   %assume equator
end

th = 0:dth:2*pi-dth;
r_E_eq = 6378137;   %[m]; equatorial radius, WGS84
eccsquared = 0.0066944; %square of the Earth's eccentricity, WGS84

%%ref: https://en.wikipedia.org/wiki/Latitude#Length_of_a_degree_of_latitude
m_per_deglat = pi*r_E_eq*(1-eccsquared)./(180*(1-eccsquared*sind(lat_cent)^2)^(3/2)); %[m/deg lat] at a given latitude, WGS84
m_per_deglon = pi*r_E_eq*cosd(lat_cent)./(180*(1-eccsquared*sind(lat_cent)^2)^(1/2)); %[m/deg lon] at a given latitude, WGS84

[TH1,R1] = meshgrid(th,r_in:dr:r_out);  %regular grid in distances
TH1 = reshape(TH1,[],1);
R1 = reshape(R1,[],1); 

[dx_out,dy_out] = pol2cart(TH1,R1);
dlat_out = dy_out./m_per_deglat;
dlon_out = dx_out./m_per_deglon;
dist_out = R1;



