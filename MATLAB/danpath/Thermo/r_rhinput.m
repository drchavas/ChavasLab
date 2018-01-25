%r_rhinput.m
%[r] = r_rhinput(p,T_K,rh); p [Pa], T_K [K], rh [0-1]
%Purpose: calculate r for given relative humidity, p, T

%Created: 8 April 2014, Dan Chavas


function [r] = r_rhinput(p,T_K,RH)

const = CM1_constants;


%%calculate saturation vapor pressures
[es] = satvappres(T_K);  %[Pa]

%%calculate actual vapor pressure
e = RH.*es;

%%calculate actual mixing ratio
r = const.eps*(e/(p-e));

end