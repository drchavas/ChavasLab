%radprof.m -- calculate azimuthal-average radial profile using
%Purpose: azimuthal-average radial profile using moving average
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
%    account_for_azimasym [1/0] - indicator of whether to use pure azimuthal
%       mean or rebin first
%    nativeres [m] - native resolution of data (determines azimuthal bin
%       size for azimuthal average, useful when azimuthal data coverage is
%       asymmetric)
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
% 24 Feb 2014; Last revision: 24 Sep 2014

% Revision history:
% 24 Feb 2014: new version of radprof.m, now using moving average; old
%   one moved to radprof_binned.m
% 19 Mar 2014: fixed bug with radial profile and NaN in if statement of
%   line 93
% 29 May 2014: added option for azimuthal re-binning to account for asymmetric data coverage
% 24 Sep 2014: removed 'nativeres' from input argument list since ignored

%------------- BEGIN CODE --------------

function [rvals, v_r, n_r, asym_mag_r, asym_th_r] = radprof(rmat,Vtot,ring_width,n_r_minthresh,thmat,account_for_azimasym,nativeres)
    
%% Create radial profiles

if(nargin >= 5) %thmat is input
    %%Convert gridpoints to cartesian (if data available)
    [th_temp] = polar_qs2matlab(thmat);   %input: theta deg CW from N [0,360)
    [xmat, ymat] = pol2cart(th_temp,rmat); %input: theta deg CCW from E (-pi,pi]
    clear th_temp
else
    thmat = NaN;
    account_for_azimasym = 0;
%    nativeres = NaN;
end

%%Initialize some things
clear num_rings
rmat_max=max(max(rmat));  %maximum radius in data
r_ringcenter = 0; %[m]
dr_ringcenter = ring_width/4;   %[m] -- radial increment for outward movement of ring center

%%Preallocate memory
rvals = NaN(1,ceil(rmat_max/dr_ringcenter));
n_r = NaN(1,ceil(rmat_max/dr_ringcenter));
v_r = NaN(1,ceil(rmat_max/dr_ringcenter));
asym_mag_r = NaN(1,ceil(rmat_max/dr_ringcenter));
asym_th_r = NaN(1,ceil(rmat_max/dr_ringcenter));

%%Iterate ring by ring to calculate radial profile
iter = 0;
rvals_temp = 0; %[m]
while(r_ringcenter<=rmat_max)
    
    %%Update radius of ring center
    r_ringcenter = r_ringcenter + dr_ringcenter; %[m]
    
    %%Update ring boundaries
    r_ring_inner = r_ringcenter-ring_width/2;   %[m]
    r_ring_inner = max([r_ring_inner 0]);   %[m]; ensure its at least 0
    r_ring_outer = r_ringcenter+ring_width/2;   %[m]
    
    %%Calculate mean wind in ring [r_ringcenter-ring_width/2,r_ringcenter+ring_width/2)
    indices=find(rmat >= r_ring_inner & rmat < r_ring_outer & ~isnan(Vtot));    %1=in ring, 0=not in ring

    %%mean radius of data in ring
    rvals_old = rvals_temp; %save old value
    rvals_temp = nanmean(rmat(indices));    %determine new value
        
    %%Only record data if mean radius is different than previous value
    if(rvals_temp ~= rvals_old && ~isnan(rvals_temp))
        iter = iter+1;
        rvals(iter) = rvals_temp;
        n_r(iter)=length(indices);  %number of gridpoints within ring that have good data   
        
        %%Calculate azimuthal average
        if(account_for_azimasym==1)

            %%Weight azimuthal data by mean distance between adjacent points
            thmax = 360;    %[deg]
            [data_weights] = weights_periodic(thmat(indices),thmax);
            assert(abs(sum(data_weights)-1)<10^-5,'Problem with data weights')
            v_r(iter) = sum(Vtot(indices).*data_weights);
            
            %%Smart azimuthal-average that accounts for data coverage asymmetry in azimuth
%             [Vtot_mean_azim_temp] = regrid_to_uniform(rvals_temp,Vtot(indices),thmat(indices),nativeres);
%             v_r(iter) = nanmean(Vtot_mean_azim_temp);
            
        else        
            %%Simple mean of all gridpoint values
            v_r(iter) = nanmean(Vtot(indices));
        end

        %%Asymmetry vector components ( = (0,0) if perfectly symmetric)
        if(nargin >= 5 && nargout == 5)   %return asymmetry vector
            asym_x = nanmean(xmat(indices))/rvals(iter);
            if(~isnan(asym_x))
                asym_y = nanmean(ymat(indices))/rvals(iter);

                [th_temp, asym_mag_r(iter)] = cart2pol(asym_x,asym_y); %input: theta deg CCW from E (-pi,pi]
                [asym_th_r(iter)] = polar_matlab2qs(th_temp);   %input: theta deg CW from N [0,360)
            else
                asym_mag_r(iter) = NaN;
                asym_th_r(iter) = NaN;
            end
            clear th_temp asym_x asym_y
        end
    end
    
end
    
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
