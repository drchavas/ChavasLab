%height_hydrostatic.m
%Calculate heights given pressures, temperatures, mixing ratios
%values start from the surface and move upwards
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%
% Outputs:
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
% 6 Jul 2016; Last revision:

%------------- BEGIN CODE --------------

function [zz] = height_hydrostatic(T,r,rl,ri,p,T_surf,r_surf,p_surf)

assert(p_surf >= max(p),'p_surf cannot be smaller than largest value in p')

%% Get all the constants we need
const = CM1_constants;

%% Calculate densities at each pressure level
[Trho_surf] = densitytemp(T_surf,r_surf,0,0);
[Trho] = densitytemp(T,r,rl,ri);

rho_surf = p_surf./(const.Rd*Trho_surf);
rho = p./(const.Rd*Trho);

rho_all = [rho_surf rho];
p_all = [p_surf p];

dp_all = p_all(2:end)-p_all(1:end-1);
rhobar_all = ((p_all(2:end).*rho_all(2:end)+p_all(1:end-1).*rho_all(1:end-1)))./(p_all(2:end) + p_all(1:end-1));

dz_all = -dp_all./(rhobar_all*const.g);
zz = cumsum(dz_all);

end

function const = CM1_constants


    % Make plots quickly by only using a few data points?
    const.doquick = 1;

    
    const.cp        = 1005.7;
    const.Rd        = 287.04;
    const.g         = 9.81;
    
    const.cpv       = 1870.0;
    const.Rv        = 461.5;
    
    const.cpl       = 4190.0;
    const.cpi       = 2106.0;

    const.cv        = const.cp-const.Rd;
    const.cvv       = const.cpv-const.Rv; 
    
    const.eps       = const.Rd/const.Rv;
    
    const.T0        = 273.15;

    % Lv = Lv1 - Lv2*T
    const.Lv0       = 2501000.0;
    const.Lv1       = const.Lv0+(const.cpl-const.cpv)*const.T0;
    const.Lv2       = const.cpl-const.cpv;
    
    % Ls = Ls1 - Ls2*T
    const.Ls0       = 2834000.0;
    const.Ls1       = const.Ls0+(const.cpi-const.cpv)*const.T0;
    const.Ls2       = const.cpi-const.cpv;
    
    
    % Goddard constants
    const.T00k      = 233.15;
    const.T0k       = 273.15;
    
    
    % reference pressure
    const.p00 = 100000;
    const.e0  = 611.2;
    
end

%------------- END OF CODE --------------