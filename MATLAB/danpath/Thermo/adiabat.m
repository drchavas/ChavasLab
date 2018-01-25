%adiabat.m - calculate moist adiabat reversible andpseudoadiabatically
%
%From Marty Singh "calculate_adiabat.m" 2012
%
% Syntax:  [Tr,r,rl,ri,Tp] = adiabat(T_surf,r_surf,p_surf,p)
%
% Inputs:
%   T_surf [K] - surface temperature
%   r_surf [kg/kg] - surface water vapor mixing ratio
%   p_surf [Pa] - surface pressure
%   p [Pa] - vector of output pressure levels
%
% Outputs:
%   Tr [K] - temperature along reversible moist adiabat
%   r [kg/kg] - water vapor mass mixing ratio
%   rl [kg/kg] - liquid water mass mixing ratio
%   ri [kg/kg] - ice water mass mixing ratio
%   Tp [K] - temperature along pseudo-adiabatic moist adiabat
%
% Example:
%    [Tr,r,rl,ri,Tp] = adiabat(T_surf,r_surf,p_surf,p);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 8 Apr 2014; Last revision:

%------------- BEGIN CODE --------------

function [Tr,r,rl,ri,Tp] = adiabat(T_surf,r_surf,p_surf,p)

% Calculates a reversible and psuedo-adiabat from an input 
% Temperature (K)
% mixing ratio (r_surf)
% pressure     (p_surf)
% 
% The adiabat Temperature and mixing ratios are given at pressures input (p)
%
% The reversible moist adiabat includes ice consistent with its treatment
% in CM1. The pseudo-adiabat is water only - it should give the same result
% as Kerry's subroutine.
%


    %% Get all the constants we need
    c = CM1_constants;


    %% Calculate the potential temperatures
    entropy = calc_entropy(p_surf,r_surf,T_surf);
    thetaep = calc_thetaep(p_surf,r_surf,T_surf);
    disp(thetaep)
    disp(entropy)
    
    %% Calculate adiabats
    [Tr,r,rl,ri] = calc_reversible(entropy,p,r_surf);
    Tp = calc_pseudo(thetaep,p,r_surf);

    
    if 0
        pi = linspace(15000,p_surf,30);
        Tpi = interp1(p,Tp,pi);
%        bfig('a5l')
        subplot(121)
        plot(Tr,p./100,'-k')
        hold on
        plot(Tpi,pi./100,'.k')
        set(gca,'ylim',[150 1000])
        set(gca,'ydir','reverse')
        ylabel('pressure (hPa)')
        xlabel('Temperature (K)')
        title('(a)')
        d = annotation('textbox','string','Reversible and psuedo-adiabtaic ascents','position',[0.3 0.9 0.4 0.1]);
        set(d,'linestyle','none','horizontalalignment','center')
        
        
        
        subplot(122)
        plot(Tr-Tp,p./100,'k','linewidth',2); set(gca,'ydir','reverse')
        xlabel('\Delta T (K)')
        title('(b)')
        set(gca,'xlim',[-1 3],'ylim',[150 1000])
    end
        
        
    
    
    
end

function [T,r,rl,ri] = calc_reversible(entropy,p,r_t)

    T = zeros(size(p));
    Tguess = 270;
    for i = 1:length(p)
        T(i) = fzero(@(x) calc_entropy(p(i),r_t,x)-entropy,Tguess);
        if T(i)>400 | T(i)<100; T(i) = nan; end
        if isfinite(T(i)); Tguess = T(i); end
    end
    
    [r,rl,ri] = goddard_microphysics(p,T,r_t);
    
    
end

function T = calc_pseudo(thetaep,p,r_t)


    T = zeros(size(p));
    Tguess = 270;
    for i = 1:length(p)
        T(i) = fzero(@(x) calc_thetaep(p(i),r_t,x)-thetaep,Tguess);
        if T(i)>400 | T(i)<100; T(i) = nan; end
        if isfinite(T(i)); Tguess = T(i); end
    end

end



function thetaep = calc_thetaep(p,r_t,T)

    %% Define some constants
    c = CM1_constants;


    [es,esl] = e_sat(T);
        
    rs = c.eps.*(es./(p-es));
    
    r = min(r_t,rs);
    e = p.*r./(c.eps+r);
    H = e./es;


    % Boltons pseudo-equivelant potential temperature (ignores ice)
    t1 = T.*(c.p00./p).^(0.2854.*(1-0.28.*r));
    Tstar = 2840./(3.5.*log(T) - log(e./100) - 4.805) + 55;
    t2 = exp(r.*(1+0.81.*r).*(3376/Tstar-2.54));

    thetaep = t1*t2;
    
end


function s = calc_entropy(p,r_t,T)

    %% Define some constants
    c = CM1_constants;

    [Lcond,Lfrz,Ldep] = Latent_heats(T);
    [es,esl,esi] = e_sat(T);
    [r,rl,ri] = goddard_microphysics(p,T,r_t);
    
    e = p.*r./(c.eps+r);
    pd = p-e;   

    % There are two ways of calculating the entropy
    % They give slightly different results, mostly because the CM1
    % microphysics are not exaclty thermodynamically consistent.
    
    % calculate specific entropy:
    % s = s_d + r_v*s_v + r_l*s_l + r_i*s_i
    %s_d = c.cp   .* log(T./c.T0) - c.Rd.* log(pd./c.p00);
    %s_v = c.cpv  .* log(T./c.T0) - c.Rv.* log(e./c.e0);
    %s_l = c.cpl  .* log(T./c.T0) - c.Lv0./c.T0;
    %s_i = c.cpi  .* log(T./c.T0) - c.Ls0./c.T0;
    %s = s_d + r.*s_v + rl.*s_l + ri.*s_i;

    
    % Derived from first law of thermodynamics 
    if(r_t>0)
        s = (c.cp + r_t.*c.cpv).*log(T) - c.Rd.*log(pd) ...
             - Ldep.*ri./T - Lcond.*rl./T ...
             - c.Rv.*( rl.*log(esl) + ri.*log(esi) + r.*log(e) );
    else
        s = c.cp.*log(T) - c.Rd.*log(pd);
    end
    
    
    
end

function [r,rl,ri] = goddard_microphysics(p,T,r_t)
c = CM1_constants;
    es = e_sat(T);
    
    rs = c.eps.*(es./(p-es));

    r = min(r_t,rs);
    
    r(r<0) = 0;
    rl = 0;
    ri = 0;
    
    % Goddard ice scheme
    fliq = ( T-c.T00k )./(c.T0k - c.T00k);
    fliq(fliq<0) = 0;
    fliq(fliq>1) = 1;
    fice = 1-fliq;

    
    rl = fliq.*(r_t-r);
    ri = fice.*(r_t-r);

    rl(abs(r_t-r)<1e-8) = 0;
    ri(abs(r_t-r)<1e-8) = 0;
    
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

function [Lcond,Lfrz,Ldep] = Latent_heats(T)

c = CM1_constants;


Lcond = c.Lv1  -  c.Lv2.*T;
Ldep =  c.Ls1  -  c.Ls2.*T;
Lfrz =    Ldep -    Lcond;

end

function [es,varargout] = e_sat(T,varargin)
  
  % Get constants
  c = CM1_constants;

  % Method we will use
  method = 'goddard';
  if nargin>1; method = varargin{1}; end


  switch method
    
    
    case 'goddard'
      
      fliq = (T-c.T00k) ./( c.T0k - c.T00k);
      fliq(T>c.T0k) = 1;
      fliq(T<c.T00k) = 0;
     
      fice = 1 - fliq;  
     
      esl = 611.2.*exp( 17.67      .* ( T  - 273.15 ) ./ ( T  - 29.65 ) );  
      esi = 611.2.*exp( 21.8745584 .* ( T  - 273.15 ) ./ ( T  - 7.66  ) );

      es = fliq.*esl + fice.*esi;

        
        
    otherwise
        error('Unknown saturation calculation method')
  end

  
  if nargout>1; varargout{1} = esl; end
  if nargout>2; varargout{2} = esi; end
  
end

%------------- END OF CODE --------------
