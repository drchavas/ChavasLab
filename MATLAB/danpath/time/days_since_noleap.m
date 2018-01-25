function d = days_since_noleap(year,month,day,year0,month0,day0)

d = datenum(year,month,day)-datenum(year0,month0,day0); %with leap days included

years_test = year0:max(year(:));    %account for all leap days

assert(~isempty(years_test),'year0 must be less than or equal to year')

iter = 0;
d_leapday = [];
for ii=1:length(years_test)
    
    year_temp = years_test(ii);
    
    if(leapyear_drc(year_temp)) %note leapyear_drc() is identical to leapyear()
        
        iter = iter + 1;
    
        d_leapday(iter) = datenum(year_temp,2,29)-datenum(year0,month0,day0); %with leap days included
        
    end
    
    
end

bins = [d_leapday max(d(:))+1];

for ii=1:length(bins)-1

    ii_good = d>=bins(ii) & d<bins(ii+1);
    d(ii_good) = d(ii_good)-ii;
    
end





