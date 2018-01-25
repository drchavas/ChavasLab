%Find_TC_from_IBTracs_latlontime.m -- Extract single track based on point
%Purpose: Input single point in time/lat/lon and return TC lifecycle
%
% Syntax:  [tt_daysince1858111700UTC_TC, Lat_BT_TC, Lon_BT_TC, Vmax_BT_TC,...
%    tt_daysince1858111700UTC_in, name_BT_TC, storm_found_in_BT, ...
%    dist_BT_median_min] = ...
%    Find_TC_from_IBTracs_latlontime(year,month,day,Lat_in,Lon_in,...
%    sources_IBTracs)
%
% Inputs:
%   year - year of case
%   month - month of case
%   day - day (with fraction) of case
%   Lat_in [deg N (-90,90]] - latitude of case center
%   Lon_in [deg E [0,360)]- longitude of case center
%   sources_IBTracs - IBTracs sources to check for Best Track match
%
% Outputs:
%   WMOID_BT_TC [] - WMO ID serial number for the storm
%   tt_daysince1858111700UTC_TC [day] - vector of times for BT lifecycle since 1858 11 17 00UTC
%   Lat_BT_TC [deg N (-90,90]] - center lat spline-interpolated from BT 
%   Lon_BT_TC [deg E [0,360)] - center lon spline-interpolated from BT 
%   Vmax_BT_TC [ms-1] - Vmax spline-interpolated from BT 
%   Pmin_BT_TC [hPa] - minimum central pressure spline-interpolated from BT
%   tt_daysince1858111700UTC_in [day] - center time spline-interpolated from BT 
%   name_BT_TC [] - Official name in IBTracs data source (often "UNNAMED")
%   storm_found_in_BT - indicator of whether a storm was found or not
%   dist_BT_median_min [deg] - median space-time distance of the 5 SMALLEST
%       distances between input data and ANY point on found TC track
%
% Example: 
%    [tt_daysince1858111700UTC_TC, Lat_BT_TC, Lon_BT_TC, Vmax_BT_TC,...
%        qs_tt_daysince1858111700UTC, name_BT_TC, storm_found_in_BT, ...
%        dist_BT_median_min] = ...
%        Find_TC_from_IBTracs_latlontime(...
%        year,month,day_hourmin,TC_lat,TC_lon,sources_IBTracs)
%
% Other m-files required: days_since
% Subfunctions: none
% MAT-files required: none
%
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 3 Dec 2013; Last revision: 9 May 2014

% Revision history:
%  22 Jan 2014 - fixed bug to now find TC based jointly on lat/lon/time
%  23 Jan 2014 - updated storm_found_in_BT to give '2' if find multiple TCs
%  23 Jan 2014 - remove data across greenwich meriden (where lon 360-->0)
%  3 Feb 2014 - updated search algorithm to minimize total space-time
%   distance between vector of input lats/lons/times and entire track
%   of every TC. This should be very robust for identifying TC.
%  4 Feb 2014 - fixed error, now will search every input IBTracs database
%   and find case with minimum distance
%  5 Feb 2014 - added output of dist_BT_ave_min, name_BT_TC
%  6 Feb 2014 - changed dist_BT_ave to dist_BT_median for cases that cross
%   prime meridian and thus include a few points with very large dLon
%  7 Feb 2014 - changed dist_BT_median to dist_BT_min to find closest
%   single point to avoid rare cases with lots of data far from BT
%  10 Feb 2014 - merged the two -- median of 5 minimum distances; fixed bug
%   in median/sort calculation for case with single input point -- now it
%   explicitly states the dimension along which to apply median/sort
%  9 May 2014 - added output for Pmin_BT_TC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- BEGIN CODE --------------

function [WMOID_BT_TC, tt_daysince1858111700UTC_TC, Lat_BT_TC, Lon_BT_TC, Vmax_BT_TC, Pmin_BT_TC, tt_daysince1858111700UTC_in, name_BT_TC, storm_found_in_BT, dist_BT_median_min] = Find_TC_from_IBTracs_latlontime(year,month,day,Lat_in,Lon_in,sources_IBTracs)

ii = 0;
Lat_BT_TC = [];
dist_BT_median_min = 10^6; %[deg] just something large
while(ii<length(sources_IBTracs))

    ii = ii + 1;

%     assert(ii<=length(sources_IBTracs),'No storm found in any of input IBTracs sources')

    source_IBTracs = sources_IBTracs{ii};

    %% Calculate total time in days since 1858-11-17 00UTC to match IBTracs
    year0_IBTracs = 1858;
    month0_IBTracs = 11;
    day0_IBTracs = 17.0;    %can use fractions of days

    tt_daysince1858111700UTC_in = days_since(year,month,day,year0_IBTracs,month0_IBTracs,day0_IBTracs);

    %% Load IBTracs data
    matdir_in = '~/Dropbox/Research/TC_DATA/IBTracs -- 2013-12-03/Data/MATLAB';
    matpath_in = sprintf('%s/IBTracs_%s.mat',matdir_in,source_IBTracs);
    load(matpath_in)
    %{
    tt_1858_11_17 = datenum(1858,11,17);
    %%A serial date number of 1 corresponds to Jan-1-0000.
    datenum(1858,11,18)-tt_1858_11_17
    %}

    %% Sum space-time distance from all input points
    num_cases = length(Lat_in);
    num_storms = size(Lat_BT,2);
    clear dist_BT_all
    
    %%Preallocate memory
    dist_BT_all = NaN(num_cases,num_storms);
    
    for mm=1:num_cases
    
        %%Distance in space
        dlat_BT = abs(Lat_BT - Lat_in(mm)); %[deg]
        dlon_BT = abs(Lon_BT - Lon_in(mm)); %[deg]

        %%Distance in time -- convert time difference to distance using characteristic translation speed
        dt_BT = abs(tt_daysince1858111700UTC_BT - tt_daysince1858111700UTC_in(mm)); %[day]
        km_per_deg=111.325; %[km/deg] at equator
        v_trans_characteristic_ms = 5;  %[ms-1]
        v_trans_characteristic_degday = v_trans_characteristic_ms*(60*60*24)/(1000*km_per_deg);  %[deg day-1]
        dt_BT_deg = v_trans_characteristic_degday*dt_BT;   %[deg]

        %%Distances to ANY point in TC Best Track (in space-time) -- minimize this
        dist_BT_all(mm,:) = nanmin(sqrt(dt_BT_deg.^2 + dlat_BT.^2 + dlon_BT.^2));    %minimum distances to ANY point along a given TC's Best Track
    end
    %%Keep the 5 smallest distances -- this focuses on those points most
        %%likely to be immediately on the Best Track, while ignoring points
        %%before/after Best Track
    dim_cases = 1;
    dist_BT_mins = sort(dist_BT_all,dim_cases,'ascend'); 
    if(num_cases>=5)    %at least 5 input points; otherwise use all data
        dist_BT_mins = dist_BT_mins(1:5,:); %five closest points only
    end
    
    %%Take the median (avoids too strongly weighting weird outliers, e.g. crossing prime meridian)
    dist_BT_median = nanmedian(dist_BT_mins,dim_cases);

    %% Find minimum space-time distance TC in Best Track database
    dist_BT_median_min_temp = min(min(dist_BT_median));
    
    %%Only update guess if find a better match than old one
    if(dist_BT_median_min_temp<dist_BT_median_min)    %use minimum as metric
        
        %%Update value
        dist_BT_median_min = dist_BT_median_min_temp;
        
        [jj_TC] = find(dist_BT_median==dist_BT_median_min);    %index of relevant TCs (columns)
        jj_TC = unique(jj_TC);

        %%TESTING: finding a lat/lon/time point near multiple TCs %%%%
        %{
        % jj_TC = [1604 1606];  %manually choose TCs from Best Track to plot
        temp = (dlat_BT <= dlat_max) & (dlon_BT <= dlon_max) & (dt_BT <= days_max);   %[deg]
        [kk_goodlatlontime] = find(temp==1);    %index of relevant TCs (columns)
        lat = Lat_BT(kk_goodlatlontime)
        lon = Lon_BT(kk_goodlatlontime)
        time = tt_daysince1858111700UTC_BT(kk_goodlatlontime)
        figure(1)
        i_TC = 1
        plot3(Lat_BT(:,jj_TC(i_TC)),Lon_BT(:,jj_TC(i_TC)),tt_daysince1858111700UTC_BT(:,jj_TC(i_TC)))
        grid on
        i_TC = 2
        hold on
        plot3(Lat_BT(:,jj_TC(i_TC)),Lon_BT(:,jj_TC(i_TC)),tt_daysince1858111700UTC_BT(:,jj_TC(i_TC)),'r')
        lat = Lat_BT(kk_goodlatlontime)
        lon = Lon_BT(kk_goodlatlontime)
        time = tt_daysince1858111700UTC_BT(kk_goodlatlontime)
        scatter3(lat,lon,time,'m')
        scatter3(Lat_in,Lon_in,tt_daysince1858111700UTC_in,'g')
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Check that it found exactly one TC
        assert(length(jj_TC)<=1,sprintf('Multiple TCs found in %s IBTracs database!',source_IBTracs))

        %% Extract relevant data for the TC
        Lon_BT_TC = Lon_BT(:,jj_TC);

        %%Check for crossings of greenwich meridien (359 --> 0 --> 1 deg E)
        dLon_BT_TC = abs(Lon_BT_TC(2:end)-Lon_BT_TC(1:end-1));
        temp = find(dLon_BT_TC > 180);  %index of last point before jump
        if(isempty(temp))   %no crossing
            ii_TC = ~isnan(Lon_BT_TC);  %use all non-nan data
        else    %crossing
            ii_TC = 1:temp;  %use all non-nan data until just before crossing
        end

        tt_daysince1858111700UTC_TC = tt_daysince1858111700UTC_BT(ii_TC,jj_TC);
        Lon_BT_TC = Lon_BT(ii_TC,jj_TC);
        Lat_BT_TC = Lat_BT(ii_TC,jj_TC);
        Vmax_BT_TC = Vmax_BT(ii_TC,jj_TC);
        Pmin_BT_TC = Pmin_BT(ii_TC,jj_TC);
        %basin_BT_TC = basin_BT(ii_TC,jj_TC);
        %classification_BT_TC = classification_BT(ii_TC,jj_TC);
        name_BT_TC = name_BT{jj_TC};
        %season_BT_TC = season_BT(jj_TC);
        WMOID_BT_TC = serialnum_BT(jj_TC);

        %basin_gen_BT_TC = basin_gen_BT(jj_TC);
    end

end

if(isempty(Lat_BT_TC))    
    storm_found_in_BT = 0;
else
    storm_found_in_BT = 1;
end

%% TEST: Plot the recovered track %%%%%
%{
figure
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
hold on
geoshow(Lat_BT_TC,Lon_BT_TC)
%}

%------------- END OF CODE --------------



    