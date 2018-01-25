function [data_out] = nc_extract_var(ncdir_in,ncfile_in,variable_in,missing_value_flag,start_vec,count_vec)
%nc_extract_var.m - Extract data for a single variable from a netcdf file
%Also sets a missing value flag to NaN
%
% Syntax:  [data_out] = nc_extract_var(ncdir_in,ncfile_in,variable_in,missing_value_flag)
%
% Inputs:
%    ncdir_in - the directory of the file you'd like
%    ncfile_in - the file you'd like
%    variable_in - the variable you'd like to extract
%    missing_value_flag - value that will be set to NaN
%    start_vec - vector of starting indices along each dimension
%    count_vec - vector of number of entries to extract along each dimension
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
% 11 Nov 2013; Last revision: 3 Jun 2014

% Revision history:
%  3 Jun 2014 - allowed for subsets of data to be extracted


%------------- BEGIN CODE --------------

%% Define full path to file
ncpath_in = sprintf('%s/%s',ncdir_in,ncfile_in);

%% Access the desired data
if(nargin==4)   %read in all data
    data_out = ncread(ncpath_in,variable_in);
elseif(nargin==6)
    data_out = ncread(ncpath_in,variable_in,start_vec,count_vec);
end

%% Set missing value flag to NaN
data_out(data_out == missing_value_flag) = NaN;

%------------- END OF CODE --------------

end


