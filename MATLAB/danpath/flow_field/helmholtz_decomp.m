%helmholtz_decomp.m
%Purpose: decompose 2D flow into radial and azimuthal components
%
% Syntax:  [Vazimx,Vazimy,Vazim,Uradx,...
%            Urady,Urad] = helmholtz_decomp(th_awayfromcenter,uu,vv,fcor_sign)
%
% Inputs:
%   th_awayfromcenter [deg CW from N [0,360)]- angle from gridpoint radially away from center
%   uu [ms-1] - x-component of flow vector
%   vv [ms-1] - y-component of flow vector
%   fcor_sign [s-1] - sign of Coriolis parameter
%
% Outputs:
%   Vazimx [ms-1] - x-component of azimuthal flow vector
%   Vazimy [ms-1] - y-component of azimuthal flow vector
%   Vazim [ms-1] - SIGNED azimuthal flow magnitude (pos = cyclonic)
%   Uradx [ms-1] - x-component of radial flow vector
%   Urady [ms-1] - y-component of radial flow vector
%   Urad [ms-1] - SIGNED radial flow magnitude (pos = outward)
%
% Example: 
%   [qs_vazimx_TC,qs_vazimy_TC,qs_Vazim_TC,qs_uradx_TC,...
%       qs_urady_TC,qs_Urad_TC] = helmholtz_decomp(th_awayfromcenter,qs_u_TC,qs_v_TC,fcor_sign);
%
% Other m-files required: none
% Subfunctions: vdist
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: --
% 10 Dec 2013; Last revision: 22 Oct 2014
%
% Revision history:
% 14 May 2014 - updated final assert to look at sum of errors rather than
%   error between sums (why the former didn't match is still a mystery!)
% 22 Oct 2014 - fixed bug that occurs when wind speed is exactly zero at a
%   point, which tripped the assert statement below since theta of flow
%   vector cant be defined if no 

%------------- BEGIN CODE --------------

%% FOR TESTING PLOTS COMMENTED OUT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Vazimx,Vazimy,Vazim,Uradx,...
%     Urady,Urad] = helmholtz_decomp(th_awayfromcenter,uu,vv,fcor_sign,xx,yy,loncent,latcent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vazimx,Vazimy,Vazim,Uradx,...
    Urady,Urad] = helmholtz_decomp(th_awayfromcenter,uu,vv,fcor_sign)

%%Convert angle to rad CCW from E (-pi,pi]
[th_awayfromcenter] = polar_qs2matlab(th_awayfromcenter);

%% TESTING: th_awayfromcenter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1040);set(gcf,'Visible','on')
hp=pcolor(xx,yy,th_awayfromcenter);set(hp, 'EdgeColor', 'none');
hold on
plot(0,0,'g.','MarkerSize',20)
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Create matrix of unit vectors parallel to azimuth CCW positive
%%Angle will be 90 deg (pi/2 rad) to the LEFT of the azimuth itself
th_azim = th_awayfromcenter + pi/2;
th_azim(th_azim>pi) = th_azim(th_azim>pi)-2*pi;

assert(min(min(th_azim))>-pi & max(max(th_azim))<=pi,'th_azim out of bounds')

V_azim = 1; %unit vector

%%Convert to cartesian
[u_azim, v_azim] = pol2cart(th_azim,V_azim); %theta deg CCW from E (-pi,pi]

%%%%%%%%%%%%%%%
%%TESTING azimuthal unit vector field%%%
%{
figure(1041);set(gcf,'Visible','on')
quiver(xx',yy',u_azim',v_azim');
hold on 
%quiver(uu',vv','m');
%plot(loncent,latcent,'g.','MarkerSize',30)
%}
%%%%%%%%%%%%%%%

%%Project wind field onto local azimuthal unit vector field

%%%%%%%%%%%%%%%
%%TESTING input vector field%%%
%{
figure(10412);set(gcf,'Visible','on')
quiver(xx',yy',uu',vv');
hold on 
%quiver(uu',vv','m');
%plot(loncent,latcent,'g.','MarkerSize',30)
%}
%%%%%%%%%%%%%%%

%%Change only points with valid data
indices_good = find(~isnan(uu));
uu_good = uu(indices_good);
vv_good = vv(indices_good);
uu_azimunit_good = u_azim(indices_good);
vv_azimunit_good = v_azim(indices_good);

%%Arbitrary initialization with same structure and NaN field as uu
Vazimx = uu;
Vazimy = vv;
qs_Vazimsign_TC = uu;
qs_Uradsign_TC = uu;

%%Note: unit vectors so don't need to divide by magnitude
if(size(uu_good,2)>1 && size(uu_good,1)==1) %make a column vector
    uu_good = uu_good';
end
if(size(vv_good,2)>1 && size(vv_good,1)==1) %make a column vector
    vv_good = vv_good';
end
if(size(uu_azimunit_good,2)>1 && size(uu_azimunit_good,1)==1) %make a column vector
    uu_azimunit_good = uu_azimunit_good';
end
if(size(vv_azimunit_good,2)>1 && size(vv_azimunit_good,1)==1) %make a column vector
    vv_azimunit_good = vv_azimunit_good';
end
    
uv_dot_uvazim = sum([uu_good.*uu_azimunit_good vv_good.*vv_azimunit_good],2);

%%Azimuthal flow vector components
Vazimx(indices_good) = uv_dot_uvazim.*uu_azimunit_good;
Vazimy(indices_good) = uv_dot_uvazim.*vv_azimunit_good;

%%Radial flow vector components
Uradx = uu-Vazimx;
Urady = vv-Vazimy;

%%Sign of azimuthal flow (cyclonic = positive)
qs_Vazimsign_TC(indices_good) = fcor_sign*sign(uv_dot_uvazim);

%%%%%%%%%%%%%%%
%%TESTING Vazim field%%%
%{
figure(10413);set(gcf,'Visible','on')
quiver(xx',yy',Vazimx',Vazimy');
hold on 
%quiver(uu',vv','m');
%plot(loncent,latcent,'g.','MarkerSize',30)
%}
%%%%%%%%%%%%%%%

%% TESTING: Vazimsign plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1042)
h=pcolor(xx,yy,qs_Vazimsign_TC);
set(h,'EdgeColor','none');
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Sign of radial flow (away from TC center = positive)
[Urad_th, ~] = cart2pol(Uradx,Urady); %theta deg CCW from E (-pi,pi]

%%Account for unlikely case that wind speed is exactly zero (i.e. Urad_th=0, but really its undefined)
indices_zerowind = uu==0 & vv==0;
Urad_th(indices_zerowind) = th_awayfromcenter(indices_zerowind);    %just set it to be radially away; it doesn't matter

dth_Urad_Uradazimunit = abs(Urad_th - th_awayfromcenter);   %difference between flow Urad and unit vector Urad -- should be either 0 or pi!

%% TESTING: dth_Urad_Uradazimunit plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1043)
h=pcolor(dth_Urad_Uradazimunit);
set(h,'EdgeColor','none');
colorbar

figure(1044)
hist(dth_Urad_Uradazimunit(:),50)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%check that difference between flow Urad and unit vector Urad is indeed either 0 or pi!
assert(sum(sum(dth_Urad_Uradazimunit > pi/16 & ...
    dth_Urad_Uradazimunit < 15*pi/16))==0,'Urad vector not aligned with radial axis')

%%Assign radial flow sign (+/- 1)
qs_Uradsign_TC(indices_good) = 1;   %set all non-NaN values to 1
qs_Uradsign_TC(dth_Urad_Uradazimunit>pi/16) = -1;

%% TESTING: Uradsign plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(12345);set(gcf,'Visible','on')
h=pcolor(qs_Uradsign_TC);
set(h,'EdgeColor','none');
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Calculate SIGNED azimuthal and radial flow magnitudes (pos = cyclonic, outward)
Vazim = qs_Vazimsign_TC.*sqrt(Vazimx.^2 + Vazimy.^2);
Urad = qs_Uradsign_TC.*sqrt(Uradx.^2 + Urady.^2);
clear qs_Vazimsign_TC qs_Uradsign_TC

%%Check for error within 1:10^6
qs_V_TC_check = sqrt(uu.^2 + vv.^2);
assert(nansum(nansum(abs(sqrt(Vazim.^2 + Urad.^2) - ...
    qs_V_TC_check))) < nansum(nansum(qs_V_TC_check))/10^6,...
    'Urad_TC and VazimTC dont match V_TC')
clear qs_V_TC_check


%% TESTING: Vazim and Urad SIGNED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1);set(gcf,'Visible','on')
clf(1)
subplot(2,1,1)
h=pcolor(Vazim);
set(h,'EdgeColor','none');
colorbar
hold on
subplot(2,1,2)
h=pcolor(Urad);
set(h,'EdgeColor','none');
colorbar
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%a_again = reshape(b, size(uu))

%{
for ii=1:length(indices_good)
    
    index_good = indices_good(ii);
    
        uu_temp = uu(index_good);
        vv_temp = vv(index_good);
        
%        if(~isnan(uu_temp))
            
            u_azim_temp = u_azim(index_good);
            v_azim_temp = v_azim(index_good);
            th_awayfromcenter_temp = th_awayfromcenter(index_good);

            %%Azimuthal flow vector -- projection
            qs_vazimvec_TC_temp = proj([uu_temp vv_temp],[u_azim_temp v_azim_temp]);

            %%Sign of azimuthal flow vector (pos = cyclonic (i.e. CCW in NH))
            qs_Vazimsign_TC_temp = fcor_sign.*sign(dot([uu_temp vv_temp],[u_azim_temp v_azim_temp]));

            %%Radial flow vector -- vector difference
            qs_uradvec_TC_temp = [uu_temp vv_temp] - qs_vazimvec_TC_temp;

            %%Sign of radial flow vector (pos = outward)
            [qs_uradth_TC_temp, ~] = cart2pol(qs_uradvec_TC_temp(1),qs_uradvec_TC_temp(2)); %theta deg CCW from E (-pi,pi]
            if(abs(qs_uradth_TC_temp - th_awayfromcenter_temp)<pi/16) 
                %angle is same as outward pointing gridpoint location vector
                qs_Uradsign_TC_temp = 1;
            else
                qs_Uradsign_TC_temp = -1;
            end

            Vazimx(index_good) = qs_vazimvec_TC_temp(1);
            Vazimy(index_good) = qs_vazimvec_TC_temp(2);
            qs_Vazimsign_TC(index_good) = qs_Vazimsign_TC_temp;
            Uradx(index_good) = qs_uradvec_TC_temp(1);
            Urady(index_good) = qs_uradvec_TC_temp(2);
            qs_Uradsign_TC(index_good) = qs_Uradsign_TC_temp;
           
end
clear uu_temp vv_temp u_azim_temp v_azim_temp th_awayfromcenter_temp ...
    qs_vazimvec_TC_temp qs_Vazimsign_TC_temp qs_uradvec_TC_temp ...
    qs_Uradsign_TC_temp
%}

%%TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Check sign -- should be mostly cyclonic and radially-inwards
%{
pcolor(qs_Vazimsign_TC)
pcolor(qs_Uradsign_TC)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- END OF CODE --------------