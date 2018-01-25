function [esat] = satvappres(T_K)
% [esat] = satvappres(T_K)
% Calculates the saturation vapor pressure [Pa] for water, Bolton formulation
% (Emanuel Eq. 4.4.14)
% Temperature (K)

    %% Get all the constants we need
    T_C = T_K - 273.15;
    esat = 611.2*exp(17.67*T_C./(T_C+243.5)); %[Pa]
    
end
