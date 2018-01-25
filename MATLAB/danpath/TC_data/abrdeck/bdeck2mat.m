%bdeck2mat.m - extract data from bdeck best track file
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

function [TC_datetime,TC_technum,TC_tech,TC_Lat,TC_Lon,TC_Vmaxms,TC_PminhPa,TC_rmaxkm,...
    TC_r34km,TC_r50km,TC_r64km,TC_firstquad,TC_rOCIkm,TC_POCIhPa,...
    TC_Name,TC_num,TC_basin,TC_tau,TC_type] = bdeck2mat(dir_in,file_in)

%% Constants
ms_kt = .5144444;   %1 kt = .514444 m/s
km_nautmi = 1.852;
m_ft = 1/3.282084;

%% Read all text in file into single big cell array (columns are variables)
path_in = sprintf('%s/%s',dir_in,file_in);
fid = fopen(path_in,'r');  % Open text file
%%ATCF format: BASIN,CY,YYYYMMDDHH,TECHNUM/MIN,TECH,TAU,LatN/S,LonE/W,VMAX,MSLP,TY,RAD,WINDCODE,RAD1,RAD2,RAD3,RAD4,RADP,RRP,MRD,GUSTS,EYE,SUBREGION,MAXSEAS,INITIALS,DIR,SPEED,STORMNAME,DEPTH,SEAS,SEASCODE,SEAS1,SEAS2,SEAS3,SEAS4,USERDEFINED,userdata
%%37 items
%line_format = '%s%f32%f32%s%s%f32%f32%f32%f32%f32%s%f32%s%f32%f32%f32%f32%f32%f32%f32%f32%f32%s%f32%s%f32%f32%s%s%f32%s%f32%f32%f32%f32%s';
%line_format = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
line_format = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
btk_all = textscan(fid, line_format, 'Delimiter', ',');
fclose(fid);

%% Extract data and convert to desired units as needed
%% BASIN
TC_basin = btk_all{1}{1};

%%CY (annual cyclone number 1-99)
TC_num = btk_all{2}{1};

%%YYYYMMDDHH (warning date-time)
TC_datetime = btk_all{3};

%%TECHNUM/MIN (objective technique sorting number, minutes for best track: 00 - 99)
TC_technum = btk_all{4};

%%TECH (acronym for each objective technique or CARQ or WRNG, BEST for best track)
TC_tech = btk_all{5};

%%TAU (forecast period: -24 through 240 hours, 0 for best-track, negative taus used for CARQ and WRNG records)
TC_tau = cell2mat(cellfun(@(x) str2num(x) ,btk_all{6}, 'UniformOutput',false));  %0 = analysis time

%%LatN/S (Latitude (tenths of degrees) for the DTG: 0 through 900, N/S is the hemispheric index)
TC_LatNS = btk_all{7};
TC_Lat = .1*cell2mat(cellfun(@(x) str2num(x(1:end-1)), TC_LatNS, 'UniformOutput',false));
TC_HemNS = cellfun(@(x) x(end), TC_LatNS, 'UniformOutput',false);
i_HemN = strcmp('N',TC_HemNS);
TC_Lat(i_HemN) = abs(TC_Lat(i_HemN));
TC_Lat(~i_HemN) = -1*abs(TC_Lat(~i_HemN)); %deg N [-90,90]
clear TC_LatNS i_HemN TC_HemNS

%%LonE/W (Longitude (tenths of degrees) for the DTG: 0 through 1800, E/W is the hemispheric index)
TC_LonEW = btk_all{8};
TC_Lon = .1*cell2mat(cellfun(@(x) str2num(x(1:end-1)), TC_LonEW, 'UniformOutput',false));
TC_HemEW = cellfun(@(x) x(end), TC_LonEW, 'UniformOutput',false);
i_HemW = strcmp('W',TC_HemEW);
TC_Lon(i_HemW) = -1*TC_Lon(i_HemW); %deg E [-180,180)
clear TC_LonEW i_HemW TC_HemEW

%%VMAX (Maximum sustained wind speed in knots: 0 through 300)
TC_Vmaxkt = cell2mat(cellfun(@(x) str2num(x) ,btk_all{9}, 'UniformOutput',false));
TC_Vmaxkt(TC_Vmaxkt==0) = NaN;  %[kt]
TC_Vmaxms = ms_kt*TC_Vmaxkt;    %[m/s]

%%MSLP (Minimum sea level pressure, 1 through 1100 MB)
TC_PminhPa = cell2mat(cellfun(@(x) str2num(x) ,btk_all{10}, 'UniformOutput',false));
TC_PminhPa(TC_PminhPa==0) = NaN;    %[hPa]

%%TY (Level of tc development; DB/TD/TS/HU)
TC_type = btk_all{11};

%%RAD (Wind intensity (kts) for the radii defined in this record: 34, 50, 64)
TC_Vradiuskt = cell2mat(cellfun(@(x) str2num(x) ,btk_all{12}, 'UniformOutput',false));
TC_Vradiuskt(TC_Vradiuskt==0) = NaN;

%%WINDCODE (Radius code: AAA - full circle; QQQ - first quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)
TC_firstquad = btk_all{13};
NEQ_test = TC_firstquad(cellfun(@(x) ~isempty(x),TC_firstquad));
assert(sum(strcmp('NEQ',NEQ_test))==length(NEQ_test),...
    'There is a starting wind radius quadrant that is not NEQ -- code cannot account for this right now!')

%%RAD1 (wind radius for region specified by WINDCODE. 0 - 1200 nm)
TC_rwindkm1 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{14}, 'UniformOutput',false));
TC_rwindkm1(TC_rwindkm1==0) = NaN;

%%RRAD2 (Full circle: not used; Quadrant: wind radius for 2nd quadrant (moving clockwise from first quadrant). 0 - 1200 nm.)
TC_rwindkm2 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{15}, 'UniformOutput',false));
TC_rwindkm2(TC_rwindkm2==0) = NaN;

%%RRAD3 (Full circle: not used; Quadrant: wind radius for 3rd quadrant (moving clockwise from first quadrant). 0 - 1200 nm.)
TC_rwindkm3 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{16}, 'UniformOutput',false));
TC_rwindkm3(TC_rwindkm3==0) = NaN;

%%RRAD4 (Full circle: not used; Quadrant: wind radius for 4th quadrant (moving clockwise from first quadrant). 0 - 1200 nm.)
TC_rwindkm4 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{17}, 'UniformOutput',false));
TC_rwindkm4(TC_rwindkm4==0) = NaN;

%%RADP (pressure in millibars of the last closed isobar, 900 - 1050 mb)
ix=cellfun(@isempty,btk_all{18});
btk_all{18}(ix)={'-9999'}; 
TC_POCIhPa = cell2mat(cellfun(@(x) str2num(x) ,btk_all{18}, 'UniformOutput',false));
TC_POCIhPa(TC_POCIhPa<=0) = NaN;

%%RRP (radius of the last closed isobar in nm, 0 - 9999 nm)
ix=cellfun(@isempty,btk_all{19});
btk_all{19}(ix)={'-9999'}; 
TC_rOCIkm = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{19}, 'UniformOutput',false));
TC_rOCIkm(TC_rOCIkm<=0) = NaN;

%%MRD (radius of max winds, 0 - 999 nm)
ix=cellfun(@isempty,btk_all{20});
btk_all{20}(ix)={'-9999'}; 
TC_rmaxkm = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{20}, 'UniformOutput',false));
TC_rmaxkm(TC_rmaxkm<=0) = NaN;

%%GUSTS (gusts, 0 through 995 kts)
ix=cellfun(@isempty,btk_all{21});
btk_all{21}(ix)={'-9999'}; 
TC_Vgustkt = cell2mat(cellfun(@(x) str2num(x) ,btk_all{21}, 'UniformOutput',false));
TC_Vgustkt(TC_Vgustkt<=0) = NaN;
TC_Vgustms = ms_kt*TC_Vgustkt;

%%EYE (eye diameter, 0 through 999 nm)
ix=cellfun(@isempty,btk_all{22});
btk_all{22}(ix)={'-9999'}; 
TC_reyekm = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{22}, 'UniformOutput',false));
TC_reyekm(TC_reyekm<=0) = NaN;

%%SUBREGION (subregion code: W, A, B, S, P, C, E, L, Q)
TC_subregion = btk_all{23};

%%MAXSEAS (max seas: 0 through 999 ft)
% ix=cellfun(@isempty,btk_all{24});
% btk_all{24}(ix)={'-9999'}; 
% TC_maxseam = m_ft*cell2mat(cellfun(@(x) str2num(x) ,btk_all{24}, 'UniformOutput',false));
% TC_maxseam(TC_maxseam<0) = NaN;

%%INITIALS (Forecaster's initials, used for tau 0 WRNG, up to 3 chars)
%TC_forecaster = btk_all{25};

%%DIR (storm direction in compass coordinates, 0 - 359 degrees)
%TC_thtrans = cell2mat(cellfun(@(x) str2num(x) ,btk_all{26}, 'UniformOutput',false));

%%SPEED (storm speed, 0 - 999 kts)
%TC_Vtrans = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{27}, 'UniformOutput',false));

%%STORMNAME (literal storm name, NONAME or INVEST)
TC_Name = btk_all{28};

%%DEPTH (system depth, D-deep, M-medium, S-shallow, X-unknown)
%TC_depthtype = btk_all{29}; %wtf does this even mean?

%%SEAS (Wave height for radii defined in SEAS1-SEAS4, 0-99 ft)
% ix=cellfun(@isempty,btk_all{30});
% btk_all{30}(ix)={'-9999'}; 
% TC_wavehtm = m_ft*cell2mat(cellfun(@(x) str2num(x) ,btk_all{30}, 'UniformOutput',false));
% TC_wavehtm(TC_wavehtm<0) = NaN;

%%SEASCODE (Radius code: AAA - full circle; QQQ - first quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)
% TC_firstquadsea = btk_all{31};

%%SEAS1 (first quadrant seas radius as defined by SEASCODE, 0 through 999 nm)
% ix=cellfun(@isempty,btk_all{32});
% btk_all{32}(ix)={'-9999'}; 
% TC_radiusseakm1 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{32}, 'UniformOutput',false));
% TC_radiusseakm1(TC_radiusseakm1<0) = NaN;

%%SEAS2 (second quadrant seas radius as defined by SEASCODE, 0 through 999 nm)
% ix=cellfun(@isempty,btk_all{33});
% btk_all{33}(ix)={'-9999'}; 
% TC_radiusseakm2 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{33}, 'UniformOutput',false));
% TC_radiusseakm2(TC_radiusseakm2<0) = NaN;

%%SEAS3 (third quadrant seas radius as defined by SEASCODE, 0 through 999 nm)
% ix=cellfun(@isempty,btk_all{34});
% btk_all{34}(ix)={'-9999'}; 
% TC_radiusseakm3 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{34}, 'UniformOutput',false));
% TC_radiusseakm3(TC_radiusseakm3<0) = NaN;

%%SEAS4 (fourth quadrant seas radius as defined by SEASCODE, 0 through 999 nm)
% ix=cellfun(@isempty,btk_all{35});
% btk_all{35}(ix)={'-9999'}; 
% TC_radiusseakm4 = km_nautmi*cell2mat(cellfun(@(x) str2num(x) ,btk_all{35}, 'UniformOutput',false));
% TC_radiusseakm4(TC_radiusseakm4<0) = NaN;

%%USERDEFINED (20 character description of format to follow) userdata
TC_userinfo1 = btk_all{36};

%%userdata
TC_userinfo2 = btk_all{37};


%% Combine wind radii data, which are stored as separate entries
[TC_datetime_unique,i_datetime_unique] = unique(TC_datetime);   %index points to first appearance
N_datetime_unique = length(TC_datetime_unique);   %number of unique times in file

%%Initialize values to NaN
TC_r34km1 = NaN(N_datetime_unique,1);
TC_r34km2 = NaN(N_datetime_unique,1);
TC_r34km3 = NaN(N_datetime_unique,1);
TC_r34km4 = NaN(N_datetime_unique,1);
TC_r50km1 = NaN(N_datetime_unique,1);
TC_r50km2 = NaN(N_datetime_unique,1);
TC_r50km3 = NaN(N_datetime_unique,1);
TC_r50km4 = NaN(N_datetime_unique,1);
TC_r64km1 = NaN(N_datetime_unique,1);
TC_r64km2 = NaN(N_datetime_unique,1);
TC_r64km3 = NaN(N_datetime_unique,1);
TC_r64km4 = NaN(N_datetime_unique,1);

for ii=1:N_datetime_unique
    
    %%Find indices of same datetime
    TC_datetime_temp = TC_datetime_unique{ii};
    i_match = strcmp(TC_datetime_temp,TC_datetime);
    indices_match = find(i_match==1);

    %%Extract corresponding wind speeds
    TC_Vradiikt_temp = TC_Vradiuskt(indices_match);
    TC_rwindskm1_temp = TC_rwindkm1(indices_match);
    TC_rwindskm2_temp = TC_rwindkm2(indices_match);
    TC_rwindskm3_temp = TC_rwindkm3(indices_match);
    TC_rwindskm4_temp = TC_rwindkm4(indices_match);

    for jj=1:length(TC_Vradiikt_temp)

        TC_Vradiuskt_temp = TC_Vradiikt_temp(jj);
        TC_rwindkm1_temp =TC_rwindskm1_temp(jj);
        TC_rwindkm2_temp =TC_rwindskm2_temp(jj);
        TC_rwindkm3_temp =TC_rwindskm3_temp(jj);
        TC_rwindkm4_temp =TC_rwindskm4_temp(jj);

        %%Store wind radii in separate vectors
        switch TC_Vradiuskt_temp
            case 34
                TC_r34km1(ii) = TC_rwindkm1_temp;
                TC_r34km2(ii) = TC_rwindkm2_temp;
                TC_r34km3(ii) = TC_rwindkm3_temp;
                TC_r34km4(ii) = TC_rwindkm4_temp;
            case 50
                TC_r50km1(ii) = TC_rwindkm1_temp;
                TC_r50km2(ii) = TC_rwindkm2_temp;
                TC_r50km3(ii) = TC_rwindkm3_temp;
                TC_r50km4(ii) = TC_rwindkm4_temp;
            case 64
                TC_r64km1(ii) = TC_rwindkm1_temp;
                TC_r64km2(ii) = TC_rwindkm2_temp;
                TC_r64km3(ii) = TC_rwindkm3_temp;
                TC_r64km4(ii) = TC_rwindkm4_temp;
        end

    end
    
end

%% Keep only data at unique datetimes
TC_technum = TC_technum(i_datetime_unique);
TC_tech = TC_tech(i_datetime_unique);
TC_tau = TC_tau(i_datetime_unique);
TC_Lat = TC_Lat(i_datetime_unique);
TC_Lon = TC_Lon(i_datetime_unique);
TC_Vmaxkt = TC_Vmaxkt(i_datetime_unique);
TC_Vmaxms = TC_Vmaxms(i_datetime_unique);
TC_PminhPa = TC_PminhPa(i_datetime_unique);
TC_type = TC_type(i_datetime_unique);
TC_firstquad = TC_firstquad(i_datetime_unique);
TC_POCIhPa = TC_POCIhPa(i_datetime_unique);
TC_rOCIkm = TC_rOCIkm(i_datetime_unique);
TC_rmaxkm = TC_rmaxkm(i_datetime_unique);
TC_Vgustkt = TC_Vgustkt(i_datetime_unique);
TC_Vgustms = TC_Vgustms(i_datetime_unique);
TC_reyekm = TC_reyekm(i_datetime_unique);
TC_subregion = TC_subregion(i_datetime_unique);
% TC_maxseam = TC_maxseam(i_datetime_unique);
%TC_forecaster = TC_forecaster(i_datetime_unique);
%TC_thtrans = TC_thtrans(i_datetime_unique);
%TC_Vtrans = TC_Vtrans(i_datetime_unique);
TC_Name = TC_Name(i_datetime_unique);
%TC_depthtype = TC_depthtype(i_datetime_unique);
% TC_wavehtm = TC_wavehtm(i_datetime_unique);
% TC_firstquadsea = TC_firstquadsea(i_datetime_unique);
% TC_radiusseakm1 = TC_radiusseakm1(i_datetime_unique);
% TC_radiusseakm2 = TC_radiusseakm2(i_datetime_unique);
% TC_radiusseakm3 = TC_radiusseakm3(i_datetime_unique);
% TC_radiusseakm4 = TC_radiusseakm4(i_datetime_unique);
TC_userinfo1 = TC_userinfo1(i_datetime_unique);
TC_userinfo2 = TC_userinfo2(i_datetime_unique);

TC_r34km = [TC_r34km1 TC_r34km2 TC_r34km3 TC_r34km4];
TC_r50km = [TC_r50km1 TC_r50km2 TC_r50km3 TC_r50km4];
TC_r64km = [TC_r64km1 TC_r64km2 TC_r64km3 TC_r64km4];



%------------- END CODE --------------