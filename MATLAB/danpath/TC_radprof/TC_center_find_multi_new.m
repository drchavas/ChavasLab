%%TC_center_find_multi_new.m

%% Created 2017-09-26, Dan Chavas

%% Purpose: Find centers by
%% 1) First defining connected regions with data values < [ppert_abs_minthresh] and N_gridpoints > [N_gridpoints_connected_minthresh]
%% 2) Second locating data-weighted center of mass from [N_gridpoints_connected_minthresh] smallest values

%% ivals and jvals can have any units

function [ivals_cent,jvals_cent,pppertnegvals_cent] = TC_center_find_multi_new(ivals,jvals,pp_pert,ppert_abs_minthresh,N_gridpoints_connected_minthresh)

%{
%% TESTING %%%%%%
clear
clc
close('all')
load temp200K.mat
%load temp250K.mat
%%vars: ivals,jvals,pp_pert,R_search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% Step 1: Find points of minimum pressure

%{
%% TESTING: PLOT pp_pert FIELD AND pp_pert < threshold
figure(91)
contourf(pp_pert-mean(pp_pert(:)))
colorbar
title('TEST: pressure perturbation [hPa]')
%%%%%%%%%%%%%%
%}

%%Define threshold as multiple of the standard deviation of the field
%threshold = mean(pp_pert(:))+3*std(pp_pert(:));    %[hPa]; MAGNITUDE of negative pressure perturbation above which to search for
pp_pert_temp = pp_pert;
pp_pert_temp(pp_pert>=-1*ppert_abs_minthresh)=0;


%% TESTING: PLOT pp_pert < threshold %%%%
%{
figure(922)
contourf(pp_pert_temp)
colorbar
title(sprintf('TEST: pressure perturbation [hPa] < -%3.2f Pa, numgridpts >= %i',ppert_abs_minthresh,N_gridpoints_connected_minthresh))
%%%%%%%%%%%%%%%%%%%%%%%%
%}

num_smooth = 0;    %number of times to apply smoother to pp_pert_temp for defining TC center
[i_cents_SLP,j_cents_SLP,pppertnegvals_cent] = LocalMinFinder_new(pp_pert_temp,N_gridpoints_connected_minthresh,num_smooth);




%% TESTING: PLOT pp_pert < threshold + minima %%%%
%{
figure(92)
contourf(pp_pert_temp)
colorbar
hold on
plot(j_cents_SLP,i_cents_SLP,'mx')
title(sprintf('TEST: pressure perturbation [hPa] < -%3.2f Pa, numgridpts >= %i',ppert_abs_minthresh,N_gridpoints_connected_minthresh))
%%%%%%%%%%%%%%%%%%%%%%%%
%}


%% Extract the actual x- and y-distances of these points
ivals_cent = [];
jvals_cent = [];
if(~isempty(i_cents_SLP))

%% Need to account for center point being at very edge of domain "between" first and last points
d_ivals = ivals(2)-ivals(1);
d_jvals = jvals(2)-jvals(1);
ivals_extended = [ivals(1)-d_ivals;ivals;ivals(end)+d_ivals];
jvals_extended = [jvals(1)-d_jvals;jvals;jvals(end)+d_jvals];
    
for ii=1:length(i_cents_SLP)
    i_cent_SLP=i_cents_SLP(ii);
    j_cent_SLP=j_cents_SLP(ii);
%    ivals_cent(ii) = ivals(j_cent_SLP);
%    jvals_cent(ii) = jvals(i_cent_SLP);
    ivals_cent(ii) = interp1(0:1:length(ivals)+1,ivals_extended,j_cent_SLP);
    jvals_cent(ii) = interp1(0:1:length(jvals)+1,jvals_extended,i_cent_SLP);
end
   
end

%% TESTING: PLOT pp_pert < threshold + minima %%%%
%{
figure(93)
contourf(pp_pert_temp)
colorbar
hold on
plot(j_cents_SLP,i_cents_SLP,'mx')
title(sprintf('TEST: pressure perturbation [hPa] < -%3.2f Pa, numgridpts >= %i',ppert_abs_minthresh,N_gridpoints_connected_minthresh))
%}
%%%%%%%%%%%%%%%%%%%%%%%%


%% For cases less than X distance from one another, choose the one with the lower pressure
%% DONT NEED TO DO THIS, USING MINIMUM # GRIDPOINTS REQUIREMENT INSTEAD
%{
%[D] = pdist2([j_cents_SLP i_cents_SLP],[j_cents_SLP i_cents_SLP],'euclidean')
if(length(jvals_cent)>1)   %only check if more than one center found
    %[D_closest,I_closest] = pdist2([j_cents_SLP i_cents_SLP],[j_cents_SLP i_cents_SLP],'euclidean','Smallest',2)
    %First closest is distance to itself (=0); thus second is what we want
    %D_closest = D_closest(2,:);
    %I_closest = I_closest(2,:);
    
    %%Need to keep ALL points and then save those that are less than
    %%R_search (in case there are more than 2 points within R_search of
    %%each other)
    D_closest = pdist2([jvals_cent ivals_cent],[jvals_cent ivals_cent],'euclidean');
    find(D_closest>0 & D_closest < R_search)
    
end
%}


end