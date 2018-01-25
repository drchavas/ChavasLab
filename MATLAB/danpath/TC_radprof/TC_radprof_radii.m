%%TC_radprof_radii.m -- calculate desired radii from input radial profile
%Purpose: Given input magnitude values, find radii beyond rmax
%
% Syntax:  [Vmax_r, rmax, rradii] = TC_radprof_radii(rr_mean,...
%               V_r_mean,V_radii_in)
%
% Inputs:
%    rr_mean [m] - vector of mean radii within each bin
%    V_r_mean [ms-1] - vector of mean values of V_xy within each radial bin
%    V_radii_in [ms-1] - vector of magnitudes (e.g. wind speed) for desired radii
%
% Outputs:
%    Vmax_r [ms-1] - maximum value of in radial profile
%    rmax [m] - radius of Vmax_r
%    rradii [m] - radii of V_radii_in
%
% Example: 
%   [Vmax_r, rmax, rradii] = TC_radprof_radii(rr_mean,V_r_mean,V_radii)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 5 Dec 2013; Last revision: 5 Dec 2013

%------------- BEGIN CODE --------------

%------------- BEGIN CODE --------------
    
function [Vmax_r, rmax, rradii] = TC_radprof_radii(rr_mean,V_r_mean,V_radii_in)

%% Vmax -- from radial profile
Vmax_r = max(V_r_mean);    %[ms-1]

%% rmax -- from radial profile
i_rmax = find(V_r_mean==Vmax_r);
rmax = rr_mean(i_rmax);  %[km]

%% Other radii -- from radial profile
V_r_mean_beyondrmax = V_r_mean(i_rmax:end);   %ignore data for r<rmax
rr_mean_beyondrmax = rr_mean(i_rmax:end);   %ignore data for r<rmax
    
for ii=1:length(V_radii_in)

    V_radii_in_temp = V_radii_in(ii);

    ii_temp_out = find(V_r_mean_beyondrmax<V_radii_in(ii),1);
    if(ii_temp_out>1)
        ii_temp_in = find(V_r_mean_beyondrmax<V_radii_in(ii),1)-1;
        v_temp_out = V_r_mean_beyondrmax(ii_temp_out);
        v_temp_in = V_r_mean_beyondrmax(ii_temp_in);
        dr_temp = rr_mean_beyondrmax(ii_temp_out)-rr_mean_beyondrmax(ii_temp_in); %mean dr in radial profile, should be very close to ring_width for reasonably large sample size

        %interpolate to radius of V=V_radii_in_temp
        rradii(ii) = rr_mean_beyondrmax(ii_temp_out) - dr_temp*((v_temp_out-V_radii_in_temp)/(v_temp_out-v_temp_in));
        clear ii_temp_in v_temp_out v_temp_in dr_temp
    else    %this wind speed does not exist
        rradii(ii) = NaN;
    end

    clear ii_temp_out

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING: Plot radial profile + radii %%%%%%%
%{
%PLOT MEAN RADIAL PROFILE
figure(1010)
hold on
plot(rr_mean/1000,V_r_mean,'g','LineWidth',2)

%%plot various points of interest in v-r space
plot(rmax/1000,Vmax_r,'rx')
symbols_radprof = {'bo' 'bs' 'bd' 'b+'};
for ii=1:length(V_radii_in)
    plot(rradii(ii)/1000,V_radii_in(ii),symbols_radprof{ii})
end

input_xlabel=sprintf('Radius');
xlabel(input_xlabel);
input_ylabel=sprintf('Wind speed');
ylabel(input_ylabel);
%}
%%END TESTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- END OF CODE --------------