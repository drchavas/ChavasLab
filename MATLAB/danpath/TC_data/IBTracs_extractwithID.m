%IBTracs_extractwithID.m -- Extract IBTracs track given input ID
%
% Syntax:  [tt_daysince1858111700UTC_TC, Lat_BT_TC, Lon_BT_TC, Vmax_BT_TC, ...
%           Pmin_BT_TC, name_BT_TC, classification_BT_TC] = ...
%           IBTracs_extractwithID(IBTracs_ID,sources_IBTracs)
%
% Inputs:
%   IBTracs_ID - official IBTracs ID (e.g. '2005236N23285' for Katrina)
%   sources_IBTracs - IBTracs sources to check IN ORDER for Best Track match
%
% Outputs:
%   tt_daysince1858111700UTC_TC [day] - vector of times for BT lifecycle since 1858 11 17 00UTC
%   Lat_BT_TC [deg N (-90,90]] - center lat spline-interpolated from BT 
%   Lon_BT_TC [deg E [0,360)] - center lon spline-interpolated from BT 
%   Vmax_BT_TC [ms-1] - Vmax spline-interpolated from BT 
%   Pmin_BT_TC [hPa] - minimum central pressure spline-interpolated from BT
%   name_BT_TC [] - Official name in IBTracs data source (often "UNNAMED")
%   classification_BT_TC [] - storm type classification (0 = TS - Tropical;
%       1 = SS - Subtropical; 2 = ET - Extratropical; 3 = DS - Disturbance)
%
% Example: 
%    [tt_daysince1858111700UTC_TC, Lat_BT_TC, Lon_BT_TC, Vmax_BT_TC, ...
%           Pmin_BT_TC, name_BT_TC, classification_BT_TC] = ...
%           IBTracs_extractwithID(IBTracs_ID,sources_IBTracs)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: IBTracs_[basin].mat for each basin
%
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 4 Nov 2014; Last revision: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- BEGIN CODE --------------

function [tt_daysince1858111700UTC_TC, Lat_BT_TC, Lon_BT_TC, Vmax_BT_TC, Pmin_BT_TC, name_BT_TC, classification_BT_TC] = IBTracs_extractwithID(IBTracs_ID,sources_IBTracs)

ii = 0;
Lat_BT_TC = [];
dist_BT_median_min = 10^6; %[deg] just something large
while(ii<length(sources_IBTracs))

    ii = ii + 1;

%     assert(ii<=length(sources_IBTracs),'No storm found in any of input IBTracs sources')

    source_IBTracs = sources_IBTracs{ii};

    %% Load IBTracs data
    matdir_in = '~/Dropbox/Research/TC_DATA/IBTracs -- 2013-12-03/Data/MATLAB';
    matpath_in = sprintf('%s/IBTracs_%s.mat',matdir_in,source_IBTracs);
    load(matpath_in)
    %{
    tt_1858_11_17 = datenum(1858,11,17);
    %%A serial date number of 1 corresponds to Jan-1-0000.
    datenum(1858,11,18)-tt_1858_11_17
    %}


    jj_TC = find(strcmp(IBTracs_ID,serialnum_BT)==1);
    
    %% Check that it found exactly one TC
    assert(length(jj_TC)<=1,sprintf('Multiple TCs found in %s IBTracs database!',source_IBTracs))

    %% Extract relevant data for the TC
    Lon_BT_TC = Lon_BT(:,jj_TC);
    
    if(sum(~isnan(Lon_BT_TC))>0)

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
        %basin_gen_BT_TC = basin_gen_BT(jj_TC);
        classification_BT_TC = classification_BT(ii_TC,jj_TC);
        name_BT_TC = name_BT{jj_TC};

        %% Found it!
        break 
        
    end


end

if(sum(~isnan(Lat_BT_TC))==0)
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
title(sprintf('%s',name_BT_TC))
%}

%------------- END OF CODE --------------



    