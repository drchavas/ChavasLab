%get_file_from_web.m -- download file from http
%Purpose:
%
% Syntax:  get_ftp_files(ftp_address,ftp_dir,ftp_file,target_dir,output_file)
%
% Inputs:
%   ftp_address - pure ftp address (e.g. mwsci.jpl.nasa.gov)
%   ftp_dir - ftp directory (e.g. subdir1/subdir2/subdir3)
%   ftp_file - desired filename (e.g. file.m)
%   target_dir - download location (e.g. ~/Documents/)
%   output_file - download filename (e.g. file_out.m)
%
% Outputs:
%
% Example: 
%    get_ftp_files(ftp_address,ftp_dir,ftp_file,target_dir,output_file)
%
% Other m-files required: none
% Subfunctions: wget
% MAT-files required: none
%
% See also:
%
% Notes: 
%
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 10 Apr 2014; Last revision:

%------------- BEGIN CODE --------------

function [] = get_file_from_web(url,filename)

urlwrite(url,filename)

end

%------------- END OF CODE --------------