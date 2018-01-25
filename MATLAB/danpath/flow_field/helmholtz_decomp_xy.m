%helmholtz_decomp_xy.m

%Purpose: Take as input a (u,v) cartesian wind field and convert it to
%polar coordinates. Assumes center is at (0,0)

function [vazimx,vazimy,vazim,uradx,...
     urady,urad] = helmholtz_decomp_xy(xmat_in,ymat_in,u_in,v_in,fcor_sign)

    %%2) Create matrix of angles (- 90 deg) for all points from center
    x_TC = 0;   %[km]
    y_TC = 0;   %[km]
    xdist=xmat_in-x_TC; %[km] distance between point and TC center
    ydist=ymat_in-y_TC; %[km] distance between point and TC center

    %%angle of vector from center to each gridpoint
%    [th,~]=cart2pol(xdist,ydist);  %E=0; N=pi/2; S=-pi/2; W=+pi

    %% Define a radial unit vector pointing outwards from the center towards each gridpoint
    rdist = sqrt(xdist.^2 + ydist.^2);
    u_unitvec = xdist./rdist;
    v_unitvec = ydist./rdist;
    assert(max(abs(sqrt(u_unitvec.^2+v_unitvec.^2)-ones(size(u_unitvec))))<10^-4,'problem with radial unit vectors')
    
%     figure(1001)
%     quiver(xdist/1000,ydist/1000,u_unitvec,v_unitvec);
    
    %% Project each wind vector onto this radial unit vector
    %%proj = (A.B)/(B.B) * Bvec
    AA = [u_in v_in];
    BB = [u_unitvec v_unitvec];
    uradvec = repmat((dot(AA,BB,2)./dot(BB,BB,2)),1,2).*BB; %positive outwards from center
    uradx = uradvec(:,1);
    urady = uradvec(:,2);
    
    %%vazim is the residual vector
    vazimvec = AA - uradvec;
    vazimx = vazimvec(:,1);
    vazimy = vazimvec(:,2);

    %% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure(1002)
    clf(1002)
    %quiver(xdist,ydist,u_in,v_in);
    quiver(xdist/1000,ydist/1000,u_in,v_in,'k');
    hold on
    quiver(xdist/1000,ydist/1000,uradx,urady,'r');
    quiver(xdist/1000,ydist/1000,vazimx,vazimy,'g');
    title('wind vectors and radial/azimuthal decomposition')
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define unsigned magnitudes
    urad = sqrt(uradx.^2+urady.^2);
    vazim = sqrt(vazimx.^2 + vazimy.^2);
    
    %% Check that the magnitudes align
    assert(max(abs(sqrt(urad(:).^2 + vazim(:).^2) - sqrt(u_in(:).^2 + v_in(:).^2)))<10^-4,'decomposed wind magnitude does not match input wind magnitude')
            
    %% Assign signs to urad (+=outward) and vazim (+ = counterclockwise)    
    %% urad
    %%urad (+=outward)
    %%north of storm (ydist > 0)
    ii_makeneg = ydist > 0 & urady<0;
    urad(ii_makeneg) = -1*urad(ii_makeneg);
    
    %%south of storm (ydist < 0)
    ii_makeneg = ydist < 0 & urady>0;
    urad(ii_makeneg) = -1*urad(ii_makeneg);
    
    %%DUE east of storm (xdist > 0 & ydist == 0)
    ii_makeneg = xdist > 0 & ydist == 0 & uradx<0;
    urad(ii_makeneg) = -1*urad(ii_makeneg);
    
    %%DUE west of storm (xdist < 0 & ydist == 0)
    ii_makeneg = xdist < 0 & ydist == 0 & uradx>0;
    urad(ii_makeneg) = -1*urad(ii_makeneg);

    %% vvazim
    %%vazim (+ = counterclockwise for f>0) -- depends on relative location
    %%north of storm (ydist > 0)
    ii_makeneg = ydist > 0 & vazimx>0;
    vazim(ii_makeneg) = fcor_sign*-1*vazim(ii_makeneg);
    
    %%south of storm (ydist < 0)
    ii_makeneg = ydist < 0 & vazimx<0;
    vazim(ii_makeneg) = fcor_sign*-1*vazim(ii_makeneg);
    
    %%DUE east of storm (xdist > 0 & ydist == 0)
    ii_makeneg = xdist > 0 & ydist == 0 & vazimy<0;
    vazim(ii_makeneg) = fcor_sign*-1*vazim(ii_makeneg);
    
    %%DUE west of storm (xdist < 0 & ydist == 0)
    ii_makeneg = xdist < 0 & ydist == 0 & vazimy>0;
    vazim(ii_makeneg) = fcor_sign*-1*vazim(ii_makeneg);
    
    %% set center to have NaN urad and vazim
    ii_makeNaN = xdist == 0 & ydist == 0;
    urad(ii_makeNaN) = NaN;
    uradx(ii_makeNaN) = NaN;
    urady(ii_makeNaN) = NaN;
    vazim(ii_makeNaN) = NaN;
    vazimx(ii_makeNaN) = NaN;
    vazimy(ii_makeNaN) = NaN;
    
    
    
    %% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure(1004)
    clf(1004)
    scatter(xdist/1000,ydist/1000,2,urad)
    colorbar
    title('signed radial wind speed (pos=outward)')
    saveas(gcf,sprintf('urad_test.jpg'),'jpeg')
    
    figure(1005)
    clf(1005)
    scatter(xdist/1000,ydist/1000,2,vazim)
    colorbar
    title('signed azimuthal wind speed (pos=CCW)')
    saveas(gcf,sprintf('vazim_test.jpg'),'jpeg')
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


end