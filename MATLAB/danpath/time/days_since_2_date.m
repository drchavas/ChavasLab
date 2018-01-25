%%days_since_2_date.m -- convert time offset to actual time
%Syntax: d = days_since_2_date(tt_days_since,year0,month0,day0)
%Output date/time format: 2012-06-15 23:28:28

function d = days_since_2_date(tt_days_since,year0,month0,day0)

d = datestr(tt_days_since+datenum(year0,month0,day0),31);
    %format 31:  'yyyy-mm-dd HH:MM:SS'  (ex: 2000-03-01 15:45:17)