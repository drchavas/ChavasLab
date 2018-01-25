%TC_radprof_distribution.m
%Purpose: Given TC wind field, plot all data and radial profiles
%
% Syntax:  [rr_mean, V_r_mean, n_r, rr_mean_RF, V_r_mean_RF, n_r_RF, rr_mean_LF, ...
%           V_r_mean_LF, n_r_LF, rr_mean_RR, V_r_mean_RR, n_r_RR, rr_mean_LR, ...
%           V_r_mean_LR, n_r_LR] = TC_radprof_quadrants(rmat,thmat,V_xy,...
%           th_translation,ring_width)
%
% Inputs:
%    rmat - matrix of radii from center
%    thmat - matrix of angles of gridpoints from center [0,360) CW from N (points in this direction)
%    V_xy - matrix of quantity at gridpoints (e.g. wind speed)
%    th_translation - translation vector angle [0,360) CW from N (points in this direction)
%    ring_width - width of radial bin for averaging (same units as rmat)
%
% Outputs:
%    rr_mean - vector of mean radii within each bin
%    V_r_mean - vector of mean values of V_xy within each radial bin
%    n_r - vector of number of datapoints within each bin
%    rr_mean_[RF/LF/RR/LR] - vector of mean radii within each bin
%       for quadrant 
%    V_r_mean_[RF/LF/RR/LR] - vector of mean values of V_xy within each
%       radial bin for quadrant
%    n_r_[RF/LF/RR/LR] - vector of number of datapoints within each bin    
%       for quadrant
%
% Example: 
%   [rr_mean, V_r_mean, n_r, rr_mean_RF, V_r_mean_RF, n_r_RF, rr_mean_LF, ...
%       V_r_mean_LF, n_r_LF, rr_mean_RR, V_r_mean_RR, n_r_RR, rr_mean_LR, ...
%       V_r_mean_LR, n_r_LR] = TC_radprof_quadrants(rr,...
%       th_pt1pt2_degCWfromN0360,qs_Vazim_TC,th_translation,ring_width);
%
% Other m-files required: radprof
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 5 Dec 2013; Last revision: 13 Jan 2014

% Revision history:
% 13 Jan 2014: added if statements to deal with quadrants with no data

%------------- BEGIN CODE --------------

function [rr_mean, V_r_mean, n_r, rr_mean_RF, V_r_mean_RF, n_r_RF, rr_mean_LF, V_r_mean_LF, n_r_LF, rr_mean_RR, V_r_mean_RR, n_r_RR, rr_mean_LR, V_r_mean_LR, n_r_LR] = TC_radprof_quadrants(rmat,thmat,V_xy,th_translation,ring_width,n_r_minthresh,account_for_azimasym,nativeres)


%% Calculate radial wind profile
assert(th_translation>=0 & th_translation<360,sprintf('Translation vector theta out of bounds: %3.2i',th_translation))
assert(min(min(thmat))>=0 & max(max(thmat))<360,sprintf('Grid theta value out of bounds'))
assert(min(min(rmat))>=0,sprintf('Radius value negative makes no sense'))
    
%%Identify indices of quadrants relative to translation vector (must account for cross-over 0 --> 360)
if(th_translation<90) %NE quadrant
    i_RF = thmat>=th_translation & thmat<th_translation+90;         %Right-front
    i_RR = thmat>=th_translation+90 & thmat<th_translation+180;     %Right-rear
    i_LR = thmat>=th_translation+180 & thmat<th_translation+270;    %Left-rear
    i_LF = thmat>=270+th_translation | thmat<th_translation;    %Left-front
elseif(th_translation<180) %SE quadrant
    i_RF = thmat>=th_translation & thmat<th_translation+90;         %Right-front
    i_RR = thmat>=th_translation+90 & thmat<th_translation+180;     %Right-rear
    i_LR = thmat>=270+(th_translation-90) | thmat<th_translation-90;    %Left-rear
    i_LF = thmat>=th_translation-90 & thmat<th_translation;    %Left-front
elseif(th_translation<270)    %SW quadrant
    i_RF = thmat>=th_translation & thmat<th_translation+90;         %Right-front
    i_RR = thmat>=270+(th_translation-180) | thmat<th_translation-180;     %Right-rear
    i_LR = thmat>=th_translation-180 & thmat<th_translation-90;    %Left-rear
    i_LF = thmat>=th_translation-90 & thmat<th_translation;    %Left-front
else        %NW quadrant
    i_RF = thmat>=270+(th_translation-270) | thmat<th_translation-270;    %Right-front
    i_RR = thmat>=th_translation-270 & thmat<th_translation-180;    %Right-rear
    i_LR = thmat>=th_translation-180 & thmat<th_translation-90;    %Left-rear
    i_LF = thmat>=th_translation-90 & thmat<th_translation;    %Left-front
end
%To plot individual quadrants, simply set values to NaN

%% Calculate mean radial wind profile for each quadrant

%%RIGHT FRONT
rmat_RF = rmat;
V_xy_RF = V_xy;
thmat_RF = thmat;
rmat_RF(i_RF==0) = NaN;
V_xy_RF(i_RF==0) = NaN;
thmat_RF(i_RF==0) = NaN;
% sprintf('Calculating radial profile for RF quadrant')
if(~isnan(max(max(V_xy_RF))))
    [rr_mean_RF,V_r_mean_RF,n_r_RF] = radprof(rmat_RF,V_xy_RF,ring_width,n_r_minthresh,thmat_RF,account_for_azimasym,nativeres);
else
    rr_mean_RF = NaN;
    V_r_mean_RF = NaN;
    n_r_RF = NaN;
end

%%RIGHT REAR
rmat_RR = rmat;
V_xy_RR = V_xy;
thmat_RR = thmat;
rmat_RR(i_RR==0) = NaN;
V_xy_RR(i_RR==0) = NaN;
thmat_RR(i_RR==0) = NaN;
% sprintf('Calculating radial profile for RR quadrant')
if(~isnan(max(max(V_xy_RR))))
    [rr_mean_RR,V_r_mean_RR,n_r_RR] = radprof(rmat_RR,V_xy_RR,ring_width,n_r_minthresh,thmat_RR,account_for_azimasym,nativeres);    %[m] [ms-1]
else
    rr_mean_RR = NaN;
    V_r_mean_RR = NaN;
    n_r_RR = NaN;
end

%%LEFT REAR
rmat_LR = rmat;
V_xy_LR = V_xy;
thmat_LR = thmat;
rmat_LR(i_LR==0) = NaN;
V_xy_LR(i_LR==0) = NaN;
thmat_LR(i_LR==0) = NaN;
% sprintf('Calculating radial profile for LR quadrant')
if(~isnan(max(max(V_xy_LR))))
    [rr_mean_LR,V_r_mean_LR,n_r_LR] = radprof(rmat_LR,V_xy_LR,ring_width,n_r_minthresh,thmat_LR,account_for_azimasym,nativeres);    %[m] [ms-1]
else
    rr_mean_LR = NaN;
    V_r_mean_LR = NaN;
    n_r_LR = NaN;
end

%%LEFT FRONT
rmat_LF = rmat;
V_xy_LF = V_xy;
thmat_LF = thmat;
rmat_LF(i_LF==0) = NaN;
V_xy_LF(i_LF==0) = NaN;
thmat_LF(i_LF==0) = NaN;
% sprintf('Calculating radial profile for LF quadrant')
if(~isnan(max(max(V_xy_LF))))
    [rr_mean_LF,V_r_mean_LF,n_r_LF] = radprof(rmat_LF,V_xy_LF,ring_width,n_r_minthresh,thmat_LF,account_for_azimasym,nativeres);    %[m] [ms-1]
else
    rr_mean_LF = NaN;
    V_r_mean_LF = NaN;
    n_r_LF = NaN;
end

%%Calculate mean radial wind profile for all data
% sprintf('Calculating radial profile for all data')
if(~isnan(max(max(V_xy))))
    [rr_mean,V_r_mean,n_r] = radprof(rmat,V_xy,ring_width,n_r_minthresh,thmat,account_for_azimasym,nativeres);    %[m] [ms-1]
else
    rr_mean = NaN;
    V_r_mean = NaN;
    n_r = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING: Plot all the radial wind profiles %%%%
%{
figure(1008)
plot(rr_mean,V_r_mean,'k')
hold on
plot(rr_mean_RF/1000,V_r_mean_RF,'r')
plot(rr_mean_RR/1000,V_r_mean_RR,'b')
plot(rr_mean_LF/1000,V_r_mean_LF,'m')
plot(rr_mean_LR/1000,V_r_mean_LR,'g')
xlabel('radius [km]')
ylabel('V [m/s]')
%}
%%%%%%%%%%%%%%%%%%

%%Put all data into a single giant vector (NOT REALLY NEEDED)
%r_vec = rmat(:);    %one giant vector
%V_xy_vec = V_xy(:);
%[r_sort i_sort] = sort(r_vec);
%V_xy_sort = V_xy_vec(i_sort);
    

%------------- END OF CODE --------------