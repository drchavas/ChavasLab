

clear
clc
close all
addpath(genpath('~/Dropbox/Research/MATLAB/'));

year0 = 1979;
month0 = 1;
day0 = 1;
month = 3;
day = 1;
years = (1979:2001)';

for ii = 1:length(years)
    year = years(ii);
    d(ii) = days_since(year,month,day,year0,month0,day0);
    d_noleap(ii) = days_since_noleap(year,month,day,year0,month0,day0);
end

sprintf('year / mo / da / true days_since / days_since no leap year / difference; relative to start of %i',year0)
[years month*ones(size(years)) day*ones(size(years)) d' d_noleap' d'-d_noleap']