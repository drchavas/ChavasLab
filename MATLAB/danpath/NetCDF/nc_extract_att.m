function [attribute_out] = nc_extract_att(ncdir_in,ncfile_in,variable_in,attribute_in)
%nc_extract_att.m - Extract attribute for a variable (or global) from a netcdf file
%
%
% Syntax:  [attribute_out] = nc_extract_att(ncdir_in,ncfile_in,variable_in,attribute_in)
%
% Inputs:
%    ncdir_in - the directory of the file you'd like
%    ncfile_in - the file you'd like
%    variable_in - the variable whose attribute you'd like to extract
%    attribute_in - the attribute you'd like to extract
%
% Outputs:
%    attribute_out - matrix of the desired data
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
% 11 Nov 2013; Last revision: 10 Jan 2014
%
% Revision history:
% 10 Jan 2014: added check for existence of attribute

%------------- BEGIN CODE --------------

%% Define full path to file
ncpath_in = sprintf('%s/%s',ncdir_in,ncfile_in);

%% Access the desired data
if(strcmp(variable_in,'global'))
    variable_in = '/';
end

%%Check if attribute exists
temp = ncinfo(ncpath_in,variable_in);
temp = struct2cell(temp.Attributes(:));
temp = temp(1,:);
temp = strcmp(attribute_in,temp);
if(sum(temp)==1)
    attribute_out = ncreadatt(ncpath_in,variable_in,attribute_in);
else
    attribute_out = 'NOATT';
end

%------------- END OF CODE --------------

end


