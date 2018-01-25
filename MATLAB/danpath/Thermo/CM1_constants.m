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