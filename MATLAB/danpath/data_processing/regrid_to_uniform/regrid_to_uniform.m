%regrid_to_uniform.m
%Purpose: Take input set of data with asymmetric azimuthal grid-spacing and regrid to
%   symmetric grid-spacing
%
% Syntax:  
%
% Inputs:
%   xx [dist or lon] - matrix of x coordinates
%
% Example: 
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 29 May 2014; Last revision:

% Revision history:

%------------- BEGIN CODE --------------

function [VV_mean_azim] = regrid_to_uniform(rr_mean,VV,th,nativeres)

%nativeres = 100*1000;    %[m]; approximate native model resolution

C = 2*pi*rr_mean;   %circumference of circle of data
N_bins = floor(C/(sqrt(2)*nativeres)); %number of bins given circumference, native data resolution
dthbinedges = (1/N_bins)*360;

%% Vectorized mean calculation in azimuthal bins
thbinedges = 0:dthbinedges:360; %bin boundaries
[th,i_sort] = sort(th); %sort by theta
VV = VV(i_sort);
[N_inbin,Bin] = histc(th,thbinedges);
N_inbin = N_inbin(1:end-1); %last entry is always zero
binnum = 1:N_bins;

mat_temp = Bin*(binnum./(binnum.^2));
mat_temp = abs(mat_temp-1)<10^-5;   %matrix of zeroes and ones

%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(sum(sum(mat_temp==1))~=length(th))
%     'hi'    
%     mat_temp = abs(mat_temp-1)<10^-5;
% else
%     mat_temp = abs(mat_temp-1)<10^-5;
% %     mat_temp(mat_temp~=1) = 0;    %this approach has precision errors!
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(sum(sum(mat_temp==1))==length(th),'Something wrong with azim-mean binning algorithm')

VV_mean_azim = (VV'*mat_temp)./N_inbin';

%{
%% Loop around azimuth in fixed bin widths
icount = 0;
for jj=1:N_bins
    
    N_inbin_temp = N_inbin(jj); %number of storms in this bin
    
    VV_mean_azim(jj) = nanmean(VV(icount+1:icount+N_inbin_temp));
   
    icount = icount + N_inbin_temp;
    
end
%}
%}




%------------- END OF CODE --------------