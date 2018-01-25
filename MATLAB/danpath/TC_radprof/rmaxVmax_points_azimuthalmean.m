%rmaxVmax_points_azimuthalmean.m -- calculate true azimuthal-mean (rmax,Vmax)
%Purpose: calculate true azimuthal-mean (rmax,Vmax) from 2d flow data
%
% Syntax:
%   [rmax_azimuthalmean,Vmax_azimuthalmean] = ...
%       rmaxVmax_points_azimuthalmean(rmat,Vtot,thmat,num_azim_bins)
%
% Inputs:
%    rmat [m] - matrix of radii from center
%    Vtot [ms-1] - quantity at each radii (e.g. wind speed)
%    thmat [deg CCW from N [0,360)] - matrix of angles for gridpoints
%    num_azim_bins [] - number of azimuthal bins (equal sized [0,360))
%
% Outputs:
%    rmax_azimuthalmean [m] - azimuthal mean point rmax across all bins
%    Vmax_azimuthalmean [ms-1] - azimuthal mean point Vmax across all bins
%
% Example: 
%   [rmax_azimuthalmean,Vmax_azimuthalmean] = ...
%       radprof(rr,VV_TC,th_grid,num_azim_bins);
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
% 19 Mar 2014; Last revision: 19 May 2014

% Revision history:
% 13 May 2014 - added if statement checking if Vtot is empty
% 19 May 2014 - added check for data coverage sufficiency in inner core

%------------- BEGIN CODE --------------

function [rmax_azimuthalmean,Vmax_azimuthalmean] = rmaxVmax_points_azimuthalmean(rmat,Vtot,thmat,num_azim_bins)
    
%% Loop over azimuthal bins and calculate (rmax,Vmax) within each bin

%%Geometry: 8 nearest neighbors (incl. diagonal) = 8 azimuthal bins to
%%ensure at least one nearest neighbor in every bin
% num_azim_bins = 8;    %[]

th_min_min = 0;
th_max_max = 360;
dth_all = th_max_max - th_min_min;
dth_bin = dth_all/num_azim_bins;
rmat_all = rmat;
thmat_all = thmat;
Vtot_all = Vtot;

%%Preallocate memory
Vmax = NaN(1,num_azim_bins);
rmax = NaN(1,num_azim_bins);
for ii=1:num_azim_bins

    %%Define azimuthal bin
    th_min_bin = th_min_min + (ii-1)*dth_bin;
    th_max_bin = th_min_min + ii*dth_bin;
    
    %%Extract data within azimuthal bin
    r_thresh = 200*1000;    %[m] distance within which to check data coverage
    max_pts_possible=sum(sum(thmat_all >= th_min_bin & thmat_all < th_max_bin & rmat_all < r_thresh));
    num_pts_good=sum(sum(thmat_all >= th_min_bin & thmat_all < th_max_bin & rmat_all < r_thresh & ~isnan(Vtot_all)));
    indices=thmat_all >= th_min_bin & thmat_all < th_max_bin & ~isnan(Vtot_all);
    thmat = thmat_all(indices);
    rmat = rmat_all(indices);
    Vtot = Vtot_all(indices);
    
    %%Take (rmax,Vmax) point from data within azimuthal bin
    data_coverage = num_pts_good/max_pts_possible;
    data_coverage_thresh = .5;  %require at least 50% of data coverage 
    if(~isempty(Vtot) & data_coverage > data_coverage_thresh)
        Vmax(ii) = max(max(Vtot));
        rmax(ii) = rmat(Vtot==Vmax(ii));
    else
        Vmax(ii) = NaN;
        rmax(ii) = NaN;
    end
    
    
    %%TESTING %%%%%%%%%%%%%%%%%%%%%
    %{
    color_plot = [0 0 ii/num_azim_bins];
    
    figure(100)
%     hold off
    plot(rmat(:)/1000,Vtot(:),'.','Color',color_plot)
    hold on
    plot(rmax(ii)/1000,Vmax(ii),'*','MarkerEdgeColor',color_plot([3 1 2]),'MarkerSize',20)
    
    figure(800)
%     hold off
    hold on
    h_pl = polar(thmat/360*2*pi,rmat/1000,'.');
    set(h_pl,'MarkerEdgeColor',color_plot,'MarkerSize',10)
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
end

rmax_azimuthalmean = nanmean(rmax);
Vmax_azimuthalmean = nanmean(Vmax);

end

%------------- END OF CODE --------------
