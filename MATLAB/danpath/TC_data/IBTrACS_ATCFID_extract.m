%IBTrACS_ATCFID_extract.m
%Purpose: Return IBTrACS info for input IBTrACS ID
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%
% Outputs:
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
% 9 Jun 2014; Last revision:

%------------- BEGIN CODE --------------

function [source_ID_out,source_out] = IBTrACS_ATCFID_extract(dir_in,IBTrACS_ID_in)

%% IBTracs filename
file_in = sprintf('%s.ibtracs.v03r05.nc',IBTrACS_ID_in);

%% Extract source info
source_info = cellstr(nc_extract_var(dir_in,file_in,'track_from_source',-9999)');

%% Check ATCF first, JTWC second
i_good = regexp(source_info,'b\w\w\d\d\d\d\d\d.dat');
i_good = ~cellfun(@isempty,i_good);

if(sum(i_good)>=1)   %ATCF entry found  (if multiple, use first one)
    source_ID_out = source_info{find(i_good==1,1)};
    source_out = 'ATCF';
else    %search for JTWC entry
    i_good = regexp(source_info,'b\w\w\d\d\d\d\d\d.txt');
    i_good = ~cellfun(@isempty,i_good);
    if(sum(i_good)>=1)   %JTWC entry found (if multiple, use first one)
        source_ID_out = source_info{find(i_good==1,1)};
        source_out = 'JTWC';
    else    %no ID found
        source_ID_out = 'NOTFOUND';
    end
end


end

%------------- END OF CODE --------------

