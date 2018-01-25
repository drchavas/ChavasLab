%%TC_center_find_multi.m

%% Created 18 Dec 2012, Dan Chavas
%% Modified 6 May 2016, Dan Chavas
%% Modified 16 Jun 2017, Dan Chavas -- distances all in the same units now

%% Purpose: To find multiple TC centers by
%% 1) First determining locations of minimum pressure below a threshold
%% 2) Second locating the center of mass near the pressure minimum
%% The code also quantifies the difference between the two, giving an error message if large

function [ivals_cent,jvals_cent] = TC_center_find_multi(ivals,jvals,pp_pert,R_search);
    %R_search [m]

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

di = ivals(2)-ivals(1);
dj = jvals(2)-jvals(1);

%% Step 1: Find points of minimum pressure
num_minima = 100;   %maximum number of minima it will search for
%km_per_deg = 111.325; %[km/deg] at equator
TC_search_radius = R_search;   %[deg]; use Vp/f -- need to scale with size, but if too small it will find multiple centers within the same storm sometimes
neighborhood_size = floor(2*TC_search_radius/di); %side-length, in gridpoints, of
if(mod(neighborhood_size,2)==0)
    neighborhood_size = neighborhood_size + 1;  %must be an odd number
end

%{
%% TESTING: PLOT pp_pert FIELD AND pp_pert < threshold
figure(91)
contourf(pp_pert-mean(pp_pert(:)))
colorbar
title('TEST: pressure perturbation [Pa]')
%%%%%%%%%%%%%%
%}

%%Define threshold as multiple of the standard deviation of the field
%threshold = mean(pp_pert(:))+3*std(pp_pert(:));    %[Pa]; MAGNITUDE of negative pressure perturbation above which to search for
threshold = 0;    %[Pa]; MAGNITUDE of negative pressure perturbation above which to search for
pp_pert_temp = pp_pert;
pp_pert_temp(pp_pert>=-threshold)=0;

%{
%% TESTING: PLOT pp_pert < threshold %%%%
figure(92)
contourf(pp_pert_temp)
colorbar
title(sprintf('TEST: pressure perturbation [Pa] < -%3.2f Pa',threshold))
%%%%%%%%%%%%%%%%%%%%%%%%
%}

num_smooth = 5;    %number of times to apply smoother to pp_pert_temp for defining TC center
%neighborhood_size = 49; 
num_minima = 100;
[i_cents_SLP,j_cents_SLP] = LocalMinFinder(pp_pert_temp,num_minima,neighborhood_size,num_smooth);

%{
%% TESTING: PLOT pp_pert < threshold + minima %%%%
figure(92)
contourf(pp_pert_temp)
colorbar
hold on
plot(j_cents_SLP,i_cents_SLP,'mx')
title(sprintf('TEST: pressure perturbation [Pa] < -%3.2f Pa',threshold))
%%%%%%%%%%%%%%%%%%%%%%%%
%}

ivals_cent = [];
jvals_cent = [];
if(~isempty(i_cents_SLP))
    
for ii=1:length(i_cents_SLP)
    i_cent_SLP=i_cents_SLP(ii);
    j_cent_SLP=j_cents_SLP(ii);
    ivals_cent(ii) = ivals(j_cent_SLP);
    jvals_cent(ii) = jvals(i_cent_SLP);
end
    
end

end