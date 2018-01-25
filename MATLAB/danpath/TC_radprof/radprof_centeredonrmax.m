%radprof_centeredonrmax.m -- calculate azimuthal-average radial profile
%Purpose: same as radprof, but radial grid is first centered on rmax
%
% Syntax:
%   [rvals, v_r, n_r] = radprof(rmat,Vtot,ring_width)
%   or
%   [rvals, v_r, n_r, asym_mag_r, asym_th_r] = radprof(rmat,Vtot,ring_width,thmat)
%   or
%   [rvals, v_r, n_r, asym_mag_r, asym_th_r] = radprof(rmat,Vtot,ring_width,thmat,n_r_minthresh)
%
% Inputs:
%    rmat [m] - matrix of radii from center
%    Vtot [ms-1] - quantity at each radii (e.g. wind speed)
%    ring_width [m] - width of radial bin for averaging (same units as rmat)
%    n_r_minthresh [] - minimum number of valid points for inclusion in radial profile
%    thmat [deg CCW from N [0,360)] - matrix of angles for gridpoints
%
% Outputs:
%    rvals [m] - vector of mean radii within each bin
%    v_r [ms-1] - vector of mean values of Vtot within each radial bin
%    n_r [-] - vector of number of datapoints within each bin
%    asym_mag_r [-] - normalized ([0,1]) magnitude of asymmetry vector
%    asym_th_r [deg CCW from N [0,360)] - angle of asymmetry vector
%
% Example: 
%   [rr_mean, VV_TC_r_mean, n_r, asym_mag, asym_th] = ...
%       radprof(rr,VV_TC,ring_width,th_grid,n_r_minthresh);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 5 Dec 2013; Last revision: 29 Jan 2014

% Revision history:
%   21 Jan 2014: added option to output asymmetry vector; fixed error in
%       calculating n_r
%   29 Jan 2014: added option for minimum n_r threshold for inclusion in
%       radial profile

%------------- BEGIN CODE --------------

function [rvals, v_r, n_r, asym_mag_r, asym_th_r] = radprof_centeredonrmax(rmat,Vtot,ring_width,n_r_minthresh,thmat)

%% Find 2d-point rmax
rmax_2d = rmat(Vtot==max(max(Vtot)));

if(nargin == 5) %thmat is input
    %%Convert gridpoints to cartesian (if data available)
    [th_temp] = polar_qs2matlab(thmat);   %input: theta deg CW from N [0,360)
    [xmat, ymat] = pol2cart(th_temp,rmat); %input: theta deg CCW from E (-pi,pi]
    clear th_temp
end

%% Create radial profiles

%%Create radial grid with bin centered around rmax_2d
clear num_rings
num_rings=floor(max(max(rmat))/ring_width);  %number of rings

rtemp = rmax_2d - ring_width/2;
r1 = mod(rtemp,ring_width);
clear rtemp
rr_grid = r1+ ring_width*(0:num_rings);   %radial grid with bin centered on rmax_2d

assert(length(rr_grid)==num_rings+1,'Problem with radial grid for radprof')

%% TESTING %%%%%%
%{
figure(1043)
clf(1043)
rdat = rr_grid(rr_grid<1.3*rmax_2d)/1000;
plot(rdat,max(max(Vtot))*ones(length(rdat),1),'rx')
hold on
plot(rmax_2d/1000,max(max(Vtot)),'go')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Iterate ring by ring to calculate radial profile
clear indices n_r v_r rvals x_r y_r n asym_mag_r asym_th_r
for n=1:num_rings   %iterate through rings
    %%Calculate mean wind in ring
    indices=find(rmat >= rr_grid(n) & rmat < rr_grid(n+1) & ~isnan(Vtot));    %1=in ring, 0=not in ring
    n_r(n)=length(indices);  %number of gridpoints within ring that have good data   

    v_r(n) = nanmean(Vtot(indices));
    rvals(n) = nanmean(rmat(indices));
    
    %%Asymmetry vector components ( = (0,0) if perfectly symmetric)
    if(nargin == 5 && nargout == 5)   %return asymmetry vector
        asym_x = nanmean(xmat(indices))/rvals(n);
        if(~isnan(asym_x))
            asym_y = nanmean(ymat(indices))/rvals(n);

            [th_temp, asym_mag_r(n)] = cart2pol(asym_x,asym_y); %input: theta deg CCW from E (-pi,pi]
            [asym_th_r(n)] = polar_matlab2qs(th_temp);   %input: theta deg CW from N [0,360)
        else
            asym_mag_r(n) = NaN;
            asym_th_r(n) = NaN;
        end
        clear th_temp asym_x asym_y
    end

end
clear xmat ymat
    
if(nargin >= 4)     %apply n_r minimum threshold
    indices_bad = find(n_r < n_r_minthresh);
    v_r(indices_bad) = NaN;
    clear indices_bad
end
%{
    %% Plot radial profile
    clf(1)
    figure(1)
    plot(rvals,v_r)
    hold on
    plot(rmax,v_r(i_rmax),'rx')
    plot(rvals,N,'r')
%}    
    
end

%------------- END OF CODE --------------
