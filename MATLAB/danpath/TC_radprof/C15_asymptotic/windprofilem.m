function [ V ] = windprofilem( vm, rm, vm2, rm2, r, wp)
save temp2 vm rm vm2 rm2 r wp
%
% This function return radial profiles of azimuthal wind, V, given matrices
% containing the maximum circular wind speeds, vm and vm2 (knots), the radii of maximum
% wind, rm  and rm2 (km), associated with primary and secondary wind maxima;
% distances r (km) of each event from the points of interest, 
% and the wind profile wp (see below; must be equal to 1, 2, 3, or 4). 
%
% Copyright WindRiskTech, L.L.C., 2013
% Last revised November, 2014
%-----------------------------------------------------------------------
%
wprofile=wp; % Use holland (1) or emanuel (2) or er2011 (3) or ec2013 (4) wind profile
%
vm=vm.*1852./3600;  % Convert maximum wind speed to m/s
vm2=vm2.*1852./3600;  % Convert maximum wind speed to m/s
%
se=sum(nonzeros(rm2));  % Test if there are any secondary eyewalls
%
if wprofile == 1
    %
    % Holland 2010 wind profile model
    %
    bs=1.8;
    rn=300/20;     % Holland parameters
    xn=0.8;
    %
    rh=r./rm;
    x=0.5+(rh-1).*(xn-0.5)./(rn-1);
    x=max(x,0.5);
    V=vm.*(rh.^-bs.*exp(1-rh.^-bs)).^x;   
    %
    if se ~= 0
        rh=r./rm2;
        x=0.5+(rh-1).*(xn-0.5)./(rn-1);
        x=max(x,0.5);
        V2=vm2.*(rh.^-bs.*exp(1-rh.^-bs)).^x;   
    end
    %
elseif wprofile == 2
    %
    % Outer radius (km)
    %
    r0=1000;
    %
    %  Re-scale radius of maximum winds by random number drawn from log-normal
    %  distribution
    %
    %r0=r0.*rfac;    
    %
    % Shape parameters
    %
    b=0.25;
    nb=0.9;
    mb=1.6;  
    %
    mb2=2.*mb;
    fac1=(1-b).*(mb+nb);
    fac2=b.*(1+2.*mb);
    fac3=2.*(mb+nb);
    fac4=2.*mb+1; 
    %
    rat=r./max(rm,1);
    V=vm.*(max((r0-r),0)./(r0-rm)).*sqrt(rat.^mb2.* ...
        (fac1./(nb+mb.*rat.^fac3)+fac2./(1+mb2.*rat.^fac4)));
    %
    if se ~= 0
        rat=r./max(rm2,1);
        V2=vm2.*(max((r0-r),0)./(r0-rm2)).*sqrt(rat.^mb2.* ...
            (fac1./(nb+mb.*rat.^fac3)+fac2./(1+mb2.*rat.^fac4)));
     end
    %
elseif wprofile == 3 
    crat=1;
    f=5.0e-5;
    f=f.*1000;  % effectively converts radii from kms to meters
    %
    Mm=rm.*vm+0.5.*f.*rm.^2;
    rn=r./rm;
    if crat == 1
        M=Mm.*(2.*rn.^2./(1+rn.^2));
    else    
        M=Mm.*(2.*rn.^2./(2-crat+crat.*rn.^2)).^(1./(2.-crat));
    end
    %V=(M-0.5.*f.*r.^2)./(r+1e-8);
    a=0.6*f./Mm-1/600^2;
    V=(M-0.5*f*r.^2./(1+a.*r.^2))./(r+1e-8);
    %V=M./(r+1e-8);  % (Add long tail to V to avoid discontinuity in vorticity at outer radius) 3/2013
    V=max(V,0);
    %
    if se ~= 0
        Mm=rm2.*vm2+0.5.*f.*rm2.^2;
        rn=r./rm2;
        if crat == 1
            M=Mm.*(2.*rn.^2./(1+rn.^2));
        else    
            M=Mm.*(2.*rn.^2./(2-crat+crat.*rn.^2)).^(1./(2.-crat));
        end
        %V2=(M-0.5.*f.*r.^2)./(r+1e-8);
        a=0.6*f./Mm-1/600^2;
        V2=(M-0.5*f*r.^2./(1+a.*r.^2))./(r+1e-8);
        %V2=M./(r+1e-8);  % (Add long tail to V to avoid discontinuity in vorticity at outer radius) 3/2013
        V2=max(V2,0);
    end
elseif wprofile == 4  % (Experimental; not recommended)
    chi=80;  %  Chi parameter for outer profile
    r0=600;  %  Outer radius (but winds not forced to zero there in this implementation)
    vp=80;   %  Nominal potential intensity. But only chi/vp matters. 
    %
    %r0=rm.*sqrt(chi.*vm./vp).*((16.*vm./(1000.*f.*rm)).^(1/3)-1).^0.75;
    ra=(rm.*r0.^2.*vp./(chi*vm)).^(1/3);
    va=2.*vm.*(ra./rm)./(1+(ra./rm).^2);
    u=max(sign(ra-r),0);
    m1=2.*rm.*vm.*(r./rm).^2./(1+(r./rm).^2); % Inner (ER2011) profile
    m2=ra.*va.*(sqrt(1+2.*chi.*ra.*va.*(r-ra).*max((r0.^2-r.^2),0)./(vp.*r0.^4))); % Outer (EC2013) profile
    V=(u.*m1+(1-u).*m2)./max(r,0.1); % Join inner and outer profiles
    %
    if se ~= 0
        ra=(rm2.*r0.^2.*vp./(chi*vm2)).^(1/3);
        va=2.*vm2.*(ra./rm2)./(1+(ra./rm2).^2);
        u=max(sign(ra-r),0);
        m1=2.*rm2.*vm2.*(r./rm2).^2./(1+(r./rm2).^2); % Inner (ER2011) profile
        m2=ra.*va.*(sqrt(1+2.*chi.*ra.*va.*(r-ra).*max((r0.^2-r.^2),0)./(vp.*r0.^4))); % Outer (EC2013) profile
        V2=(u.*m1+(1-u).*m2)./max(r,0.1); % Join inner and outer profiles    
    end
end
%
% Merge primary and secondary wind profiles
%
if se ~= 0
    u=max(r,1)./rm;
    hu=max(sign(u-1),0);
    V=V.*(1-hu)+hu.*vm./u;
    u=r./max(rm2,1);
    hu=max(sign(u-1),0);
    V2=V2.*hu+(1-hu).*vm2.*u;
    V=max(V,V2);
end    
%
V=V.*3600./1852;    % Convert wind speed to knots
%
end