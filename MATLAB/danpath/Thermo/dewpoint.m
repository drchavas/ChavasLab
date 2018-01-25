function [T_d] = dewpoint(e)
% [T_d] = dewpoint(e)
% Calculates the dew point [K] when given the vapor pressure [Pa], Bolton
% formulation inverted
% (Emanuel Eq. 4.6.2)
% Temperature (K)

    T_d = 273.15 + 243.5./((17.67./log(e/611.2))-1);
    
end
