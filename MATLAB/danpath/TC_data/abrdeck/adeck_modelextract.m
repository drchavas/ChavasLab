%adeck_modelextract.m - extract model-specific data from input adeck data
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

function [TC_technum_out,TC_tech_out,TC_tau_out,...
    TC_Lat_out,TC_Lon_out,TC_Vmaxms_out,TC_PminhPa_out,TC_rmaxkm_out,...
    TC_r34km_out,TC_r50km_out,TC_r64km_out,...
    TC_firstquad_out,TC_rOCIkm_out,TC_POCIhPa_out,...
    TC_basin_out,TC_type_out,TC_Name_out] = ...
    adeck_modelextract(models_in,TC_technum,TC_tech,TC_tau,...
    TC_Lat,TC_Lon,TC_Vmaxms,TC_PminhPa,TC_rmaxkm,...
    TC_Vradiuskt,TC_rwindkm,...
    TC_firstquad,TC_rOCIkm,TC_POCIhPa,...
    TC_basin,TC_type,TC_Name)

%% Extract data for individual models from datetime-specific data
switch models_in{1}
    case 'all'
        Model_names = unique(TC_tech);
    otherwise
        Model_names = models_in;
end

N_models = length(Model_names);

for jj=1:N_models

    Model_model = Model_names{jj};
    indices_model = strcmp(Model_model,TC_tech);

    TC_technum_model = TC_technum(indices_model);
    TC_tech_model = TC_tech(indices_model);
    TC_Lat_model = TC_Lat(indices_model);
    TC_Lon_model = TC_Lon(indices_model);
    TC_Vmaxms_model = TC_Vmaxms(indices_model);
    TC_PminhPa_model = TC_PminhPa(indices_model);
    TC_rmaxkm_model = TC_rmaxkm(indices_model);
    TC_Vradiuskt_model = TC_Vradiuskt(indices_model);
    TC_rwindkm_model = TC_rwindkm(indices_model,:);
    TC_firstquad_model = TC_firstquad(indices_model);
    TC_rOCIkm_model = TC_rOCIkm(indices_model);
    TC_POCIhPa_model = TC_POCIhPa(indices_model);

    TC_basin_model = TC_basin(indices_model);
    TC_tau_model = TC_tau(indices_model);
    TC_type_model = TC_type(indices_model);
    TC_Name_model = TC_Name(indices_model);

    %% IMPLEMENT bdeck2mat METHODOLOGY FOR EXTRACTING r34, r50, r64 and compressing data to 1 per timestep
    %% Combine wind radii data, which are stored as separate entries
    [TC_tau_model_unique,i_tau_model_unique] = unique(TC_tau_model);   %index points to first appearance
    N_tau_model_unique = length(TC_tau_model_unique);   %number of unique times in file

    %%Initialize values to NaN
    TC_r34km_model = NaN(N_tau_model_unique,4);
    TC_r50km_model = NaN(N_tau_model_unique,4);
    TC_r64km_model = NaN(N_tau_model_unique,4);

    for kk=1:N_tau_model_unique

        %%Find indices of same tau
        TC_tau_model_tau = TC_tau_model_unique(kk);
        indices_match = find(TC_tau_model==TC_tau_model_tau);

        %%Extract corresponding wind speeds (34, 50, and/or 64 kt)
        TC_Vradiikt_model_tau = TC_Vradiuskt_model(indices_match);
        TC_rwindskm_model_tau = TC_rwindkm_model(indices_match,:);

        %%Loop over each wind speed entry
        for ll=1:length(TC_Vradiikt_model_tau)

            %%Extract wind speed entry, corresponding wind radii
            TC_Vradiuskt_model_tau = TC_Vradiikt_model_tau(ll);
            TC_rwindkm_model_tau = TC_rwindskm_model_tau(ll,:);

            %%Store wind radii in separate vectors (34, 50, 64 kt)
            switch TC_Vradiuskt_model_tau
                case 34
                    TC_r34km_model(kk,:) = TC_rwindkm_model_tau;
                case 50
                    TC_r50km_model(kk,:) = TC_rwindkm_model_tau;
                case 64
                    TC_r64km_model(kk,:) = TC_rwindkm_model_tau;
            end

        end

    end

    %% Some tests
    assert(length(TC_tau_model(i_tau_model_unique))==size(TC_r34km_model,1),'Matrix dimensions of _out data do not line up!')
    assert(length(TC_Vmaxms_model(i_tau_model_unique))==size(TC_r34km_model,1),'Matrix dimensions of _out data do not line up!')

    %% Keep only data at unique taus
    TC_technum_out{jj} = TC_technum_model(i_tau_model_unique);
    %TC_tech_out{jj} = TC_tech_model(i_tau_model_unique);   %defined below
    TC_tau_out{jj} = TC_tau_model_unique;
    TC_Lat_out{jj} = TC_Lat_model(i_tau_model_unique);
    TC_Lon_out{jj} = TC_Lon_model(i_tau_model_unique);
    TC_Vmaxms_out{jj} = TC_Vmaxms_model(i_tau_model_unique);
    TC_PminhPa_out{jj} = TC_PminhPa_model(i_tau_model_unique);
    TC_rmaxkm_out{jj} = TC_rmaxkm_model(i_tau_model_unique);
    TC_firstquad_out{jj} = TC_firstquad_model(i_tau_model_unique);
    TC_rOCIkm_out{jj} = TC_rOCIkm_model(i_tau_model_unique);
    TC_POCIhPa_out{jj} = TC_POCIhPa_model(i_tau_model_unique);

    TC_basin_out{jj} = TC_basin(i_tau_model_unique);
    TC_type_out{jj} = TC_type(i_tau_model_unique);
    TC_Name_out{jj} = TC_Name(i_tau_model_unique);
    
    TC_r34km_out{jj} = TC_r34km_model;
    TC_r50km_out{jj} = TC_r50km_model;
    TC_r64km_out{jj} = TC_r64km_model;

    
    
end
TC_tech_out = Model_names;

%------------- END CODE --------------