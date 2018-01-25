function d = days_since(year,month,day,year0,month0,day0)
%d = days_since(year,month,day,year0,month0,day0)

d = datenum(year,month,day)-datenum(year0,month0,day0);