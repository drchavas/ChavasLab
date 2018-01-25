function [i_cents_SLP,j_cents_SLP] = LocalMinFinder(pp_pert_temp,num_minima,neighborhood_size,num_smooth);
%LocalMinFinder.m
%Purpose: find local minima in input 2D field
%
% Syntax:  [i_cents_SLP,j_cents_SLP] = ...
%           LocalMinFinder(pp_pert_temp,num_minima,neighborhood_size,num_smooth)
%
% Inputs:
%   pp_pert_temp -- input 2D field
%   num_minima -- maximum number of minima to return
%   neighborhood_size [# gridpts] -- side length of box for search
%   num_smooth [#] -- number of applications of 9-pt (3x3) smoother before
%       minima are extracted
%
% Outputs:
%   i_cents_SLP -- row numbers (i.e. index of first dim in pp_pert_temp) of center
%   j_cents_SLP -- column numbers (i.e. index of second dim in pp_pert_temp) of center
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
% 26 May 2014; Last revision:

%------------- BEGIN CODE --------------

pp_pert_neg = -pp_pert_temp;

%% Need to copy pp_pert field in all directions to account for cases with TC near edge of domain
pp_pert_neg_repeat = repmat(pp_pert_neg,3,3);

%% Smooth the data [num_smooth]-times with a 9-point 2D filter
h = 1/9*ones(3);
for ii=1:num_smooth
    pp_pert_neg_repeat = filter2(h,pp_pert_neg_repeat);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure(1000+ii)
    contourf(pp_pert_neg_repeat)
    max(max(pp_pert_neg_repeat))
    colorbar
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%{
%% TESTING: Did copied field work? %%%%%
figure(95)
set(gcf,'Visible','on')
contourf(-1*pp_pert_neg_repeat)
title(sprintf('TEST: p-pert [Pa] < -3-sig; copied then smoothed'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

hLocalMax = vision.LocalMaximaFinder;
hLocalMax.MaximumNumLocalMaxima = num_minima;    %find the [num_minima] largest local maxima
hLocalMax.NeighborhoodSize = [neighborhood_size neighborhood_size];   %need to look at a broad area to find it
hLocalMax.Threshold = 100/1000;   %[hPa]; just need something a bit larger than 0 (all values below threshold were zeroed out)

%% Find locations of all local maxes in large copied domain
location = step(hLocalMax, pp_pert_neg_repeat);
i_cents_SLP = location(:,2);
j_cents_SLP = location(:,1);

%% If desired: algorithm that checks if multiple points very close together. If yes, then take mean
%%But this should really be implicit already in [neighborhood_size]

%% Extract from here only the locations of local maxes within our original domain
ni_orig = size(pp_pert_neg,1);
nj_orig = size(pp_pert_neg,2);
i_orig_0 = ni_orig+1;   %first x-point in original domain
i_orig_f = 2*ni_orig;   %last x-point in original domain
% i_orig_0 = 1;   %first y-point in original domain
% i_orig_f = ni_orig;   %last y-point in original domain
j_orig_0 = nj_orig+1;   %first y-point in original domain
j_orig_f = 2*nj_orig;   %last y-point in original domain
indices_orig = find(i_cents_SLP>=i_orig_0 & i_cents_SLP<=i_orig_f & j_cents_SLP>=j_orig_0 & j_cents_SLP<=j_orig_f);
i_cents_SLP = i_cents_SLP(indices_orig)-ni_orig;
%i_cents_SLP = i_cents_SLP(indices_orig);
j_cents_SLP = j_cents_SLP(indices_orig)-nj_orig;

%}

%{
%% TESTING: Did it extract the correct centers? %%%%%
figure(96)
clf(96)
set(gcf,'Visible','on')
contourf(-1*pp_pert_neg)
%contourf(-1*pp_pert_neg_repeat)
colorbar
hold on
plot(j_cents_SLP,i_cents_SLP,'w*','MarkerSize',20)
title(sprintf('TEST: p-pert [Pa] < -3-sig; TC centers (%ix 9-pt smooth)',num_smooth))
%title(sprintf('loc1 = (%i,%i); loc2 = (%i,%i)',location(1,1),location(1,2),location(2,1),location(2,2)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%end
