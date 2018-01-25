%adeck2mat.m - extract data from adeck forecast file for a specified time
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    ncdir_dr - the directory of the file you'd like
%    ncfile_in - the file you'd like
%    variable_in - the variable you'd like to extract
%    missing_value_flag - value that will be set to NaN
%
% Outputs:
%    data_out - matrix of the desired data
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
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
% 12 Sep 2014; Last revision:

%------------- BEGIN CODE --------------
% 
% clear all
% clc
% close all
% addpath(genpath('~/Dropbox/Research/MATLAB/'));

% function [TC_datetime,TC_technum,TC_tech,TC_Lat,TC_Lon,TC_Vmaxms,TC_PminhPa,TC_rmaxkm,...
%     TC_Vradiuskt,TC_rwindkm,TC_firstquad,TC_rOCIkm,TC_POCIhPa,...
%     TC_Name,TC_num,TC_basin,TC_tau,TC_type] = adeck2mat(dir_in,file_in,datetime_in)

function [TC_datetime,TC_technum_out,TC_tech_out,TC_tau_out,...
    TC_Lat_out,TC_Lon_out,TC_Vmaxms_out,TC_PminhPa_out,TC_rmaxkm_out,...
    TC_r34km_out,TC_r50km_out,TC_r64km_out,...
    TC_firstquad_out,TC_rOCIkm_out,TC_POCIhPa_out,...
    TC_basin_out,TC_type_out,TC_Name_out] = ...
    adeck2mat(dir_in,file_in,datetime_in,models_in)

if(nargin==3)
    models_in = 'all';
end

%% Constants
ms_kt = .5144444;   %1 kt = .514444 m/s
km_nautmi = 1.852;
m_ft = 1/3.282084;

%% Extract data from adeck forecast file for specified datetime
[TC_datetimes] = adeck_datetimesextract(dir_in,file_in);    %vector of all times
TC_datetimes_unique = unique(TC_datetimes);

%%Special time option, if desired
switch datetime_in
    case 'last'
        datetime_in = TC_datetimes_unique(end);
    case 'first'
        datetime_in = TC_datetimes_unique(1);
end

%%Check if time exists
if(sum(strcmp(TC_datetimes_unique,datetime_in))==0)
    assert(1==2,'datetime_in not found for this storm!')
end

%%Identify indices for specified datetime
indices_datetime = strcmp(datetime_in,TC_datetimes);
N_lines_datetime = sum(indices_datetime);
istart_datetime = find(indices_datetime==1,1,'first');
iend_datetime = find(indices_datetime==1,1,'last');
assert(iend_datetime - istart_datetime + 1 == N_lines_datetime,'datetime_in found in non-consecutive entries!')


%% Read all text in file into single big cell array (columns are variables)
path_in = sprintf('%s/%s',dir_in,file_in);
fid = fopen(path_in,'r');  % Open text file
%%ATCF format: BASIN,CY,YYYYMMDDHH,TECHNUM/MIN,TECH,TAU,LatN/S,LonE/W,VMAX,MSLP,TY,RAD,WINDCODE,RAD1,RAD2,RAD3,RAD4,RADP,RRP,MRD,GUSTS,EYE,SUBREGION,MAXSEAS,INITIALS,DIR,SPEED,STORMNAME,DEPTH,SEAS,SEASCODE,SEAS1,SEAS2,SEAS3,SEAS4,USERDEFINED,userdata
%%37 items

%%Skip over lines before specified time
line_format_junk = '%*[^\n]';
junk = textscan(fid, line_format_junk, istart_datetime-1, 'Delimiter', ',');
clear junk

%%Extract data for specified time
line_format = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
data_all_datetime = textscan(fid, line_format, N_lines_datetime, 'Delimiter', ',');

fclose(fid);

%% Extract data and convert to desired units as needed
%% BASIN
TC_basin = data_all_datetime{1};

%%CY (annual cyclone number 1-99)
TC_num = data_all_datetime{2};
TC_num = unique(TC_num);
assert(length(TC_num)==1,'multiple TC numbers found')

%%YYYYMMDDHH (warning date-time)
TC_datetime = data_all_datetime{3};
TC_datetime = unique(TC_datetime);
assert(length(TC_datetime)==1,'More than one date found!')

%%TECHNUM/MIN (objective technique sorting number, minutes for best track: 00 - 99)
TC_technum = data_all_datetime{4};

%%TECH (acronym for each objective technique or CARQ or WRNG, BEST for best track)
TC_tech = data_all_datetime{5};

%%TAU (forecast period: -24 through 240 hours, 0 for best-track, negative taus used for CARQ and WRNG records)
TC_tau = cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{6}, 'UniformOutput',false));  %0 = analysis time

%%LatN/S (Latitude (tenths of degrees) for the DTG: 0 through 900, N/S is the hemispheric index)
TC_LatNS = data_all_datetime{7};
TC_Lat = .1*cell2mat(cellfun(@(x) str2num(x(1:end-1)), TC_LatNS, 'UniformOutput',false));
TC_HemNS = cellfun(@(x) x(end), TC_LatNS, 'UniformOutput',false);
i_HemN = strcmp('N',TC_HemNS);
TC_Lat(i_HemN) = abs(TC_Lat(i_HemN));
TC_Lat(~i_HemN) = -1*abs(TC_Lat(~i_HemN)); %deg N [-90,90]
clear TC_LatNS i_HemN TC_HemNS

%%LonE/W (Longitude (tenths of degrees) for the DTG: 0 through 1800, E/W is the hemispheric index)
TC_LonEW = data_all_datetime{8};
TC_Lon = .1*cell2mat(cellfun(@(x) str2num(x(1:end-1)), TC_LonEW, 'UniformOutput',false));
TC_HemEW = cellfun(@(x) x(end), TC_LonEW, 'UniformOutput',false);
i_HemW = strcmp('W',TC_HemEW);
TC_Lon(i_HemW) = -1*TC_Lon(i_HemW); %deg E [-180,180)
clear TC_LonEW i_HemW TC_HemEW

%%VMAX (Maximum sustained wind speed in knots: 0 through 300)
TC_Vmaxkt = cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{9}, 'UniformOutput',false));
TC_Vmaxkt(TC_Vmaxkt==0) = NaN;  %[kt]
TC_Vmaxms = ms_kt*TC_Vmaxkt;    %[m/s]

%%MSLP (Minimum sea level pressure, 1 through 1100 MB)
TC_PminhPa = cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{10}, 'UniformOutput',false));
TC_PminhPa(TC_PminhPa==0) = NaN;    %[hPa]

%%TY (Level of tc development; DB/TD/TS/HU)
TC_type = data_all_datetime{11};

%%RAD (Wind intensity (kts) for the radii defined in this record: 34, 50, 64)
TC_Vradiuskt = cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{12}, 'UniformOutput',false));
TC_Vradiuskt(TC_Vradiuskt==0) = NaN;

%%WINDCODE (Radius code: AAA - full circle; QQQ - first quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)
TC_firstquad = data_all_datetime{13};

%%RAD1 (wind radius for region specified by WINDCODE. 0 - 1200 nm)
TC_rwindkm1 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{14}, 'UniformOutput',false));
TC_rwindkm1(TC_rwindkm1==0) = NaN;

%%RRAD2 (Full circle: not used; Quadrant: wind radius for 2nd quadrant (moving clockwise from first quadrant). 0 - 1200 nm.)
TC_rwindkm2 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{15}, 'UniformOutput',false));
TC_rwindkm2(TC_rwindkm2==0) = NaN;

%%RRAD3 (Full circle: not used; Quadrant: wind radius for 3rd quadrant (moving clockwise from first quadrant). 0 - 1200 nm.)
TC_rwindkm3 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{16}, 'UniformOutput',false));
TC_rwindkm3(TC_rwindkm3==0) = NaN;

%%RRAD4 (Full circle: not used; Quadrant: wind radius for 4th quadrant (moving clockwise from first quadrant). 0 - 1200 nm.)
TC_rwindkm4 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{17}, 'UniformOutput',false));
TC_rwindkm4(TC_rwindkm4==0) = NaN;

%%RADP (pressure in millibars of the last closed isobar, 900 - 1050 mb)
ix=cellfun(@isempty,data_all_datetime{18});
data_all_datetime{18}(ix)={'-9999'}; 
TC_POCIhPa = cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{18}, 'UniformOutput',false));
TC_POCIhPa(TC_POCIhPa<=0) = NaN;

%%RRP (radius of the last closed isobar in nm, 0 - 9999 nm)
ix=cellfun(@isempty,data_all_datetime{19});
data_all_datetime{19}(ix)={'-9999'}; 
TC_rOCIkm = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{19}, 'UniformOutput',false));
TC_rOCIkm(TC_rOCIkm<=0) = NaN;

%%MRD (radius of max winds, 0 - 999 nm)
ix=cellfun(@isempty,data_all_datetime{20});
data_all_datetime{20}(ix)={'-9999'}; 
TC_rmaxkm = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{20}, 'UniformOutput',false));
TC_rmaxkm(TC_rmaxkm<=0) = NaN;

%%GUSTS (gusts, 0 through 995 kts)
ix=cellfun(@isempty,data_all_datetime{21});
data_all_datetime{21}(ix)={'-9999'}; 
TC_Vgustkt = cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{21}, 'UniformOutput',false));
TC_Vgustkt(TC_Vgustkt<=0) = NaN;
TC_Vgustms = ms_kt*TC_Vgustkt;

%%EYE (eye diameter, 0 through 999 nm)
ix=cellfun(@isempty,data_all_datetime{22});
data_all_datetime{22}(ix)={'-9999'}; 
TC_reyekm = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{22}, 'UniformOutput',false));
TC_reyekm(TC_reyekm<=0) = NaN;

%%SUBREGION (subregion code: W, A, B, S, P, C, E, L, Q)
TC_subregion = data_all_datetime{23};

%%MAXSEAS (max seas: 0 through 999 ft)
% ix=cellfun(@isempty,data_all_datetime{24});
% data_all_datetime{24}(ix)={'-9999'}; 
% TC_maxseam = m_ft*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{24}, 'UniformOutput',false));
% TC_maxseam(TC_maxseam<0) = NaN;

%%INITIALS (Forecaster's initials, used for tau 0 WRNG, up to 3 chars)
%TC_forecaster = data_all_datetime{25};

%%DIR (storm direction in compass coordinates, 0 - 359 degrees)
%TC_thtrans = cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{26}, 'UniformOutput',false));

%%SPEED (storm speed, 0 - 999 kts)
%TC_Vtrans = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{27}, 'UniformOutput',false));

%%STORMNAME (literal storm name, NONAME or INVEST)
TC_Name = data_all_datetime{28};

%%DEPTH (system depth, D-deep, M-medium, S-shallow, X-unknown)
%TC_depthtype = data_all_datetime{29}; %wtf does this even mean?

%%SEAS (Wave height for radii defined in SEAS1-SEAS4, 0-99 ft)
% ix=cellfun(@isempty,data_all_datetime{30});
% data_all_datetime{30}(ix)={'-9999'}; 
% TC_wavehtm = m_ft*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{30}, 'UniformOutput',false));
% TC_wavehtm(TC_wavehtm<0) = NaN;

%%SEASCODE (Radius code: AAA - full circle; QQQ - first quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)
% TC_firstquadsea = data_all_datetime{31};

%%SEAS1 (first quadrant seas radius as defined by SEASCODE, 0 through 999 nm)
% ix=cellfun(@isempty,data_all_datetime{32});
% data_all_datetime{32}(ix)={'-9999'}; 
% TC_radiusseakm1 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{32}, 'UniformOutput',false));
% TC_radiusseakm1(TC_radiusseakm1<0) = NaN;

%%SEAS2 (second quadrant seas radius as defined by SEASCODE, 0 through 999 nm)
% ix=cellfun(@isempty,data_all_datetime{33});
% data_all_datetime{33}(ix)={'-9999'}; 
% TC_radiusseakm2 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{33}, 'UniformOutput',false));
% TC_radiusseakm2(TC_radiusseakm2<0) = NaN;

%%SEAS3 (third quadrant seas radius as defined by SEASCODE, 0 through 999 nm)
% ix=cellfun(@isempty,data_all_datetime{34});
% data_all_datetime{34}(ix)={'-9999'}; 
% TC_radiusseakm3 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{34}, 'UniformOutput',false));
% TC_radiusseakm3(TC_radiusseakm3<0) = NaN;

%%SEAS4 (fourth quadrant seas radius as defined by SEASCODE, 0 through 999 nm)
% ix=cellfun(@isempty,data_all_datetime{35});
% data_all_datetime{35}(ix)={'-9999'}; 
% TC_radiusseakm4 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,data_all_datetime{35}, 'UniformOutput',false));
% TC_radiusseakm4(TC_radiusseakm4<0) = NaN;

%%USERDEFINED (20 character description of format to follow) userdata
TC_userinfo1 = data_all_datetime{36};

%%userdata
TC_userinfo2 = data_all_datetime{37};


%% Store wind radii
TC_rwindkm = [TC_rwindkm1 TC_rwindkm2 TC_rwindkm3 TC_rwindkm4];

%% Extract data for all desired models
[TC_technum_out,TC_tech_out,TC_tau_out,...
    TC_Lat_out,TC_Lon_out,TC_Vmaxms_out,TC_PminhPa_out,TC_rmaxkm_out,...
    TC_r34km_out,TC_r50km_out,TC_r64km_out,...
    TC_firstquad_out,TC_rOCIkm_out,TC_POCIhPa_out,...
    TC_basin_out,TC_type_out,TC_Name_out] = ...
    adeck_modelextract(models_in,TC_technum,TC_tech,TC_tau,...
    TC_Lat,TC_Lon,TC_Vmaxms,TC_PminhPa,TC_rmaxkm,...
    TC_Vradiuskt,TC_rwindkm,...
    TC_firstquad,TC_rOCIkm,TC_POCIhPa,...
    TC_basin,TC_type,TC_Name);

%------------- END CODE --------------