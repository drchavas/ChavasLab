%C15_asymptotic_ke.m
%Purpose: Kerry's code for analytical approximation

%Reference: "Approximate Solutions for Outer, Inner, and Merged Wind
%Fields" Kerry Emanuel email 05/2013

%------------- BEGIN CODE --------------

function [v,r0,ra,va,r] = C15_asymptotic_ke(vm,rm,f,chi)

vp=80;
chi = chi*vp;
rm = rm/1000;

gam=1000.*f.*rm/vm;
r0=rm.*sqrt(chi.*vm.*(1+0.5.*gam).^2/(vp*(1+gam))).*((2*(2+gam).*(1+gam)).^(2/3).*gam.^(-1/3)-1-gam).^0.75;
r=0:1.2*floor(r0);
ra=(rm*r0.^2*vp.*(1+gam)/(chi*vm*(1+0.5*gam).^2))^(1/3);
va=2*vm*(ra/rm)*(1+0.5*gam).^2./(1+gam+(ra/rm)^2);
q=max(size(r));
m=zeros(1,q);
for i=1:q,
    if r(i) <= ra
        m(i)=2*rm*vm*(1+0.5*gam)^2*(r(i)/rm)^2/(1+gam+(r(i)/rm)^2);
        m(i)=m(i)-0.5.*f.*1000.*r(i).^2;
        m(i)=max(m(i),0);
    else
        %m(i)=ra*va*(sqrt(1+2*chi*ra*va*(r(i)-ra)*(r0^2-r(i)^2)/(vp*r0^4))-((r(i)-ra)/(r0-ra))^2);
        m(i)=ra*va*(sqrt(1+2*chi*ra*va*(r(i)-ra)*(r0^2-r(i)^2)/(vp*r0^4)));
        m(i)=m(i)-0.5.*f.*1000.*r(i).^2./(1+(r(i)./ra-1).^1.5);
        %m(i)=m(i)-offset.*r(i);
    end
end
%m=m-0.1.*f.*1000.*r.^2;
%m=max(m,0);
v=zeros(1,q);
v(1,2:end)=m(2:end)./r(2:end);

ra = ra*1000;   %convert to [m]
r = r*1000;   %convert to [m]
r0 = r0*1000;   %convert to [m]
%
%
%{
h=plot(rg,vg,'b',r,v,'--b');
set(h,'linewidth',2)
set(gca,'fontweight','bold')
if max(size(pf)) >= 3 && strcmp(pf(1:3),'kat')
    legend('Observations','Analytic Model')
else
    legend('Numerical Model','Analytic Model')
end    
hold on
ia=q-1;
for i=1:q,
    if r(i) < ra
        ia=i;
    end
end
w=(ra-r(ia))/(r(2)-r(1));
rra=w*r(ia+1)+(1-w)*r(ia);
vra=w*v(ia+1)+(1-w)*v(ia);
j=errorbar(rra,vra,5);
set(j,'linewidth',2,'color','r')
hold off
set(gca,'xlim',[0 r0+50])
xlabel('Radius (km)')
ylabel('Azimuthal Wind Speed (m/s)')
title(pf)
%}

%------------- END OF CODE --------------