function success=tskew(pz,tz,rhz)
%
% pz = pressure (hPa), tz = temperature (C), rhz = relative humidity (0-1)
%
ez=6.112.*exp(17.67.*tz./(243.5+tz));
qz=rhz.*0.622.*ez./(pz-ez);
chi=log(pz.*qz./(6.112.*(0.622+qz)));
tdz=243.5.*chi./(17.67-chi);
%
p=[1050:-25:100];
pplot=transpose(p);
t0=[-48:2:50];
[ps1,ps2]=size(p);
ps=max(ps1,ps2);
[ts1,ts2]=size(t0);
ts=max(ts1,ts2);
for i=1:ts,
   for j=1:ps,
      tem(i,j)=t0(i)+30.*log(0.001.*p(j));
      thet(i,j)=(273.15+tem(i,j)).*(1000./p(j)).^.287;
      es=6.112.*exp(17.67.*tem(i,j)./(243.5+tem(i,j)));
      q(i,j)=622.*es./(p(j)-es);
      thetaea(i,j)=thet(i,j).*exp(2.5.*q(i,j)./(tem(i,j)+273.15));
   end
end
p=transpose(p);
t0=transpose(t0);
temp=transpose(tem);
theta=transpose(thet);
thetae=transpose(thetaea);
qs=transpose(sqrt(q));
h=contour(t0,pplot,temp,16,'k');
hold on
set(gca,'ytick',[1000:100:100])
set(gca,'yscale','log','ydir','reverse')
set(gca,'fontweight','bold')
set(gca,'ytick',[100:100:1000])
set(gca,'ygrid','on')
hold on
h=contour(t0,pplot,theta,24,'b'); %dry adiabate
h=contour(t0,pplot,qs,24,'g'); %constant mixing ratio
h=contour(t0,pplot,thetae,24,'r'); %moist adiabate
%tsound=30.+43.5.*log(0.001.*p);
%tsoundm=tsound-30.*log(0.001.*p);
tzm=tz-30.*log(0.001.*pz);
tdzm=tdz-30.*log(0.001.*pz);
g=9.80665; M=0.0289644; R=8.31447; L=0.0065; ex_10=(R*L)/(g*M);
T0=(tzm(find(min(pz)))+273.15);
altz=(1-(pz/1013.25).^ex_10)*T0/L;
if 0
h=plot(tzm,pz,'k',tdzm,pz,'k--');
else
[ax, h] = plot2axes(tzm,pz,'k',...
    tdzm,pz,'k--','yscale', @(pz) (1-(pz/1013.25).^ex_10)*T0/L/1000);
    set(ax,'fontweight','bold','ydir','reverse')
    ylabel(ax(2), 'Altitude (km)')
end
set(h,'linewidth',2)

hold off
xlabel('Temperature (C)','fontweight','bold')
ylabel('Pressure (mb)','fontweight','bold')
success='yes';