%fixedwidthtxt_to_cell.m
%
% Syntax:  get_ebtr(basin)
%
% Inputs:
%   basin [] -- string of desired basin ('atl, 'epac')
%
% Outputs: none, file saved in current directory
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: get_ftp_files
% Subfunctions: none
% MAT-files required: none
%
% Source: http://www.mathworks.com/matlabcentral/answers/67819-how-do-a-read-a-text-file-as-fixed-width-columns-as-in-excel

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 10 Apr 2014; Last revision:

%------------- BEGIN CODE --------------


function [data_cell] = fixedwidthtxt_to_cell(filename,colw)

buffer = fileread(filename) ;
buffer(buffer<' ') = [] ;                           % Remove \n,\r, etc. [EDITED]
buffer = reshape(buffer, sum(colw), []).' ;

%%Check that things work
assert(sum(colw)==size(buffer,2),sprintf('Column widths (%i) do sum to line length (%i)!',sum(colw),size(buffer,2)))

data_cell = strtrim(mat2cell(buffer, ones(size(buffer,1),1), colw));   %strtrim removes white space at start and end of strings

end

%------------- END OF CODE --------------