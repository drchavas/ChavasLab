function [i_cents_SLP,j_cents_SLP,pp_pert_neg_cents_SLP] = LocalMinFinder_new(pp_pert_temp,N_gridpoints_connected_minthresh,num_smooth)
%LocalMinFinder.m
%Purpose: find local minima in input 2D field
%
% Syntax:  [i_cents_SLP,j_cents_SLP] = ...
%            LocalMinFinder_new(pp_pert_temp,N_gridpoints_connected_minthresh,num_smooth)
%
% Inputs:
%   pp_pert_temp -- input 2D field
%   N_gridpoints_connected_minthresh -- # gridpoints required to count as
%       connected region, also # gridpoints used to define data-weighted center
%   num_smooth [#] -- number of applications of 9-pt (3x3) smoother before
%       minima are extracted (default = 0)
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

%%Default num_smooth = 0
if(nargin < 3)
    num_smooth = 0;
end

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

% hLocalMax = vision.LocalMaximaFinder;
% %hLocalMax.MaximumNumLocalMaxima = num_minima;    %find the [num_minima] largest local maxima
% hLocalMax.NeighborhoodSize = [neighborhood_size neighborhood_size];   %need to look at a broad area to find it
% hLocalMax.Threshold = 100/1000;   %[hPa]; just need something a bit larger than 0 (all values below threshold were zeroed out)

%% Find locations of all local maxes in large copied domain
% location = step(hLocalMax, pp_pert_neg_repeat);
% i_cents_SLP = location(:,2);
% j_cents_SLP = location(:,1);

%regmax = imregionalmax(pp_pert_neg_repeat); %uses 8-connected neighborhood as default
                % (i.e. must be greater than adjacent including corners)
%[i_cents_SLP,j_cents_SLP] = find(regmax==1);


pp_pert_binary = pp_pert_neg_repeat>0;
regmax_binary = imregionalmax(pp_pert_binary); %uses 8-connected neighborhood as default
                % (i.e. must be greater than adjacent including corners)
                
CC = bwconncomp(regmax_binary);
N_connectedregions = length(CC.PixelIdxList);
%i_cents_SLP = NaN(N_connectedregions,1);
%j_cents_SLP = NaN(N_connectedregions,1);
i_cents_SLP = [];
j_cents_SLP = [];
pp_pert_neg_cents_SLP = [];
iter = 0;
for jj=1:N_connectedregions
    [ii_connected,jj_connected] = ind2sub(size(regmax_binary),CC.PixelIdxList{jj});

    %% Check that region has sufficient # of gridpoints to keep
    if(length(ii_connected)>=N_gridpoints_connected_minthresh)
 
        iter = iter + 1;
        %Find centroid index of connected region
        %%THIS ASSUMES CONSTANT HORIZONTAL GRID SPACING
%        i_cents_SLP(jj) = mean(ii_connected);
%        j_cents_SLP(jj) = mean(jj_connected);
%        i_cents_SLP(iter) = mean(ii_connected);
%        j_cents_SLP(iter) = mean(jj_connected);
        
        %Find PPPERT-WEIGHTED centroid index of connected region
        %%THIS ASSUMES CONSTANT HORIZONTAL GRID SPACING
%        i_cents_SLP(jj) = mean(ii_connected);
%        j_cents_SLP(jj) = mean(jj_connected);

        pp_pert_neg_repeat_ijconnected = pp_pert_neg_repeat(sub2ind(size(pp_pert_neg_repeat),ii_connected,jj_connected));
        assert(isequal(size(pp_pert_neg_repeat_ijconnected),size(jj_connected)),'pp_pert_neg_repeat_ijconnected different size from jj_connected')
        
        %% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        figure(999)
        hist(pp_pert_neg_repeat_ijconnected)    %this should be have tail skewed towards large values!
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%Define centroid ONLY from the N_gridpts_centroid
        [~,kk_nearcenter] = sort(pp_pert_neg_repeat_ijconnected,'descend');
        kk_nearcenter = kk_nearcenter(1:N_gridpoints_connected_minthresh);
        ii_connected_nearcenter = ii_connected(kk_nearcenter);
        jj_connected_nearcenter = jj_connected(kk_nearcenter);
        pp_pert_neg_repeat_ijconnected_nearcenter = pp_pert_neg_repeat_ijconnected(kk_nearcenter);
        
        %% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        figure(998)
        hist(pp_pert_neg_repeat_ijconnected_nearcenter)    %this should be ONLY large values!
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        i_cent_SLP_temp = sum(pp_pert_neg_repeat_ijconnected_nearcenter.*ii_connected_nearcenter)/sum(pp_pert_neg_repeat_ijconnected_nearcenter);
        j_cent_SLP_temp = sum(pp_pert_neg_repeat_ijconnected_nearcenter.*jj_connected_nearcenter)/sum(pp_pert_neg_repeat_ijconnected_nearcenter);
        pp_pert_neg_cent_SLP_temp = interp2(1:size(pp_pert_neg_repeat,1),1:size(pp_pert_neg_repeat,2),pp_pert_neg_repeat,j_cent_SLP_temp,i_cent_SLP_temp,'linear');

        %% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        figure(905)
        scatter(ii_connected,jj_connected)
        hold on
        scatter(ii_connected_nearcenter,jj_connected_nearcenter,'rx')
        plot(i_cent_SLP_temp,j_cent_SLP_temp,'gx')
        input_title = sprintf('pppertcent = -%3.3f hPa; min pppert grid value = -%3.3f hPa (should be close)',pp_pert_neg_cent_SLP_temp,max(pp_pert_neg_repeat_ijconnected_nearcenter));
        title(input_title)
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% If valid data, then save (interp2 can give a NaN for pp_pert_neg_cent_SLP_temp at very edge of domain in some versions of MATLAB)
        if(~isnan(pp_pert_neg_cent_SLP_temp))
            dP_test = abs(pp_pert_neg_cent_SLP_temp-max(pp_pert_neg_repeat_ijconnected_nearcenter));
            %assert(dP_test<3,sprintf('central pressure estimate differs from min nearby grid value by more than 3 hPa (dP = %3.3f hPa)',dP_test))
            if(dP_test<3)   %only save the data if it passes this test; otherwise it is probably from two storms very close in each other's orbit
        
                %% Save the data
                i_cents_SLP(iter) = i_cent_SLP_temp;
                j_cents_SLP(iter) = j_cent_SLP_temp;
                pp_pert_neg_cents_SLP(iter) = pp_pert_neg_cent_SLP_temp;
                
            else
                
                sprintf('STORM CENTER IGNORED: central pressure estimate differs from min nearby grid value by more than 3 hPa (dP = %3.3f hPa)',dP_test)
                
            end
                
        end
        
        
%    else
%        i_cents_SLP(jj) = NaN;
%        j_cents_SLP(jj) = NaN;
    end
    
    
    
end


%% TESTING: PLOT pp_pert < threshold + minima %%%%
%{
figure(923)
contourf(pp_pert_neg_repeat)
colorbar
hold on
plot(j_cents_SLP,i_cents_SLP,'mx')
title(sprintf('TEST: TILED pressure perturbation with centers'))
%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% Extract from here only the locations of local maxes within our original domain
ni_orig = size(pp_pert_neg,1);
nj_orig = size(pp_pert_neg,2);
%% For domain centered on (0,0)
%i_orig_0 = ni_orig+0.5;   %halfway between last x-point in left tile and first x-point in original domain
%i_orig_f = 2*ni_orig+0.5;   %halfway between last x-point in original domain and first x-point in right tile
%j_orig_0 = nj_orig+0.5;   %halfway between last y-point in bottom tile and first y-point in original domain
%j_orig_f = 2*nj_orig+0.5;   %halfway between last y-point in original domain and first y-point in top tile
%% For domain with distances that are positive in each direction starting at (0,0)
i_orig_0 = ni_orig+1;   %first x-point in original domain
i_orig_f = 2*ni_orig+1;   %one point beyond last x-point in original domain
j_orig_0 = nj_orig+1;   %first y-point in original domain
j_orig_f = 2*nj_orig+1;   %one point beyond last y-point in original domain

indices_orig = find(i_cents_SLP>=i_orig_0 & i_cents_SLP<i_orig_f & j_cents_SLP>=j_orig_0 & j_cents_SLP<j_orig_f);
i_cents_SLP = i_cents_SLP(indices_orig)-ni_orig;
j_cents_SLP = j_cents_SLP(indices_orig)-nj_orig;
pp_pert_neg_cents_SLP = pp_pert_neg_cents_SLP(indices_orig);

assert(isequal(size(pp_pert_neg_cents_SLP),size(i_cents_SLP)),'mismatch in size pp_pert_neg_cents_SLP and i_cents_SLP')

%}


%% TESTING: Did it extract the correct centers? %%%%%
%{
figure(96)
clf(96)
set(gcf,'Visible','on')
contourf(-1*pp_pert_neg)
colorbar
hold on
plot(j_cents_SLP,i_cents_SLP,'w*','MarkerSize',20)
title(sprintf('TEST: p-pert [Pa] < -3-sig; TC centers (%ix 9-pt smooth)',num_smooth))
%title(sprintf('loc1 = (%i,%i); loc2 = (%i,%i)',location(1,1),location(1,2),location(2,1),location(2,2)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%'hi'
%end
