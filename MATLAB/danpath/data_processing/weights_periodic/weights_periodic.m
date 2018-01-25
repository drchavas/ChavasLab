%weights_periodic.m
%Purpose: Assign weights based on spacing between datapoints
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

function [data_weights] = weights_periodic(th_data,thmax)

if(size(th_data,1)==1)
    th_data = th_data';
end

[th_data,i_sort] = sort(th_data);
unsorted = 1:length(th_data);
i_unsort(i_sort) = unsorted;
clear unsorted

th_data_temp = [th_data(end)-thmax;th_data;th_data(1)+thmax];

mean_spacing = (th_data_temp(3:end)-th_data_temp(1:end-2))/2;
assert(abs(sum(mean_spacing)-thmax)<10^5,'Problem with mean spacing algorithm')
normfac = thmax;
data_weights_sorted = mean_spacing/normfac;

data_weights = data_weights_sorted(i_unsort);

assert(sum(data_weights_sorted~=data_weights(i_sort))==0,'problem with unsorting algorithm')

%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1)
clf(1)
scatter(th_data,ones(size(th_data)),50,data_weights_sorted)
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- END OF CODE --------------