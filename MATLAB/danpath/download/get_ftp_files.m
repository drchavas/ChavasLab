%get_ftp_files.m -- download file from ftp server
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
% Notes: -wc is used, so it will not clobber (i.e. copy over) existing files
%
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 9 Jan 2014; Last revision: 9 Jan 2014

%------------- BEGIN CODE --------------

function [] = get_ftp_files(ftp_address,ftp_dir,ftp_file,target_dir,output_file)

if(nargin==3)
    target_dir = '.';
    output_file = ftp_file;
end

%{
%%MATLAB's built in mget -- a little slower (~10%) than system() method
ftp_address = 'mwsci.jpl.nasa.gov';
ftp_dir = 'outgoing/north_atlantic/2005/katrina/quikscat/data';
ftp_file = 'SEAWINDS_QUIKSCAT_L2_WIND_20050822_2305.nc';
output_dir = '~/Dropbox/Research/2013/QuikSCAT_new/RAW_DATA/';
mw = ftp(ftp_address);
cd(mw,ftp_dir)
mget(mw,ftp_file);
close(mw);
%}

%%System() method -- do via unix
ftp_path = sprintf('ftp://%s/%s/%s',ftp_address,ftp_dir,ftp_file);

home_dir = pwd;
cd(target_dir)

ftp_download_command =sprintf('/opt/local/bin/wget ''%s'' -O %s',ftp_path,output_file) %will overwrite
system(ftp_download_command);

cd(home_dir)

%------------- END OF CODE --------------