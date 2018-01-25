function [cape, cin]=getcape( nk , p_in , t_in , td_in);
%-----------------------------------------------------------------------

% original code in fortran by G. Bryan - see below
% 2011 fortran version was used to translated to matlab
% matlab input variables: 
% Input variables: nk - number of levels in the sounding (integer)
%        p_in - one-dimensional array of pressure (mb) (real)
%        t_in - one-dimensional array of temperature (C) (real)
%        td_in - one-dimensional array of dewpoint temperature (C) (real)
%
%  Output variables:  cape - Convective Available Potential Energy (J/kg) (real)
%           cin - Convective Inhibition (J/kg) (real)

% Some general notes on CAPE 
%
% 1-1500 J/kg: positive
% 1500-2500 J/kg: large
% 2500+ J/kg: extreme
% example webresources:
% http://www.theweatherprediction.com/habyhints/305/
% http://www.atmos.albany.edu/deas/atmclasses/atm301/CAPE.htm
% other parameters: wmax = sqrt(2*cape);

%  original code and notes from here on:
%  getcape - a fortran90 subroutine to calculate Convective Available
%            Potential Energy (CAPE) from a sounding.
%
%  Version 1.02                           Last modified:  10 October 2008
%
%  Author:  George H. Bryan
%           Mesoscale and Microscale Meteorology Division
%           National Center for Atmospheric Research
%           Boulder, Colorado, USA
%           gbryan@ucar.edu
%
%  Disclaimer:  This code is made available WITHOUT WARRANTY.
%
%  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
%               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
%
%-----------------------------------------------------------------------
%
%  Input:     nk - number of levels in the sounding (integer)
%
%           p_in - one-dimensional array of pressure (mb) (real)
%
%           t_in - one-dimensional array of temperature (C) (real)
%
%          td_in - one-dimensional array of dewpoint temperature (C) (real)
%
%  Output:  cape - Convective Available Potential Energy (J/kg) (real)
%
%            cin - Convective Inhibition (J/kg) (real)
%
%-----------------------------------------------------------------------
%  User options:
% Pressure increment (Pa)
clear global; clear functions;

persistent adiabat avgqv avgth b1 b2 cloud converge cp cpdg cpdrd cpi cpl cpm cpv debug_level doit dp dz eps fice fliq frac g i ice k kmax lfc lhf lhs lhv ls1 ls2 lv1 lv2 maxthe ml_depth n narea nloop not_converged orec p p00 p1 p2 parea pb pc pi pi1 pi2 pinc pn ps pt ptv q qi1 qi2 qibar ql1 ql2 qlbar qt qv1 qv2 qvbar rd rddcp reps rm rp00 rv source t t0 t1 t2 tbar td th th1 th2 the thlast thv thv1 thv2 xls xlv z ; 

if isempty(pinc), pinc = 100.0; end;
% (smaller number yields more accurate
%  results,larger number makes code
%  go faster)
% Source parcel:
if isempty(source), source = 1; end;
% 1 = surface
% 2 = most unstable (max theta-e)
% 3 = mixed-layer (specify ml_depth)
% depth (m) of mixed layer
if isempty(ml_depth), ml_depth =  200.0; end;
% for source=3
% Formulation of moist adiabat:
if isempty(adiabat), adiabat = 1; end;
% 1 = pseudoadiabatic, liquid only
% 2 = reversible, liquid only
% 3 = pseudoadiabatic, with ice
% 4 = reversible, with ice
%-----------------------------------------------------------------------
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%-----------------------------------------------------------------------
%            No need to modify anything below here:
%-----------------------------------------------------------------------
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%-----------------------------------------------------------------------
if isempty(doit), doit=false; end;
if isempty(ice), ice=false; end;
if isempty(cloud), cloud=false; end;
if isempty(not_converged), not_converged=false; end;
if isempty(k), k=0; end;
if isempty(kmax), kmax=0; end;
if isempty(n), n=0; end;
if isempty(nloop), nloop=0; end;
if isempty(i), i=0; end;
if isempty(orec), orec=0; end;
if isempty(p), p=zeros(1,nk); end;
if isempty(t), t=zeros(1,nk); end;
if isempty(td), td=zeros(1,nk); end;
if isempty(pi), pi=zeros(1,nk); end;
if isempty(q), q=zeros(1,nk); end;
if isempty(th), th=zeros(1,nk); end;
if isempty(thv), thv=zeros(1,nk); end;
if isempty(z), z=zeros(1,nk); end;
if isempty(pt), pt=zeros(1,nk); end;
if isempty(pb), pb=zeros(1,nk); end;
if isempty(pc), pc=zeros(1,nk); end;
if isempty(pn), pn=zeros(1,nk); end;
if isempty(ptv), ptv=zeros(1,nk); end;
if isempty(the), the=0; end;
if isempty(maxthe), maxthe=0; end;
if isempty(parea), parea=0; end;
if isempty(narea), narea=0; end;
if isempty(lfc), lfc=0; end;
if isempty(th1), th1=0; end;
if isempty(p1), p1=0; end;
if isempty(t1), t1=0; end;
if isempty(qv1), qv1=0; end;
if isempty(ql1), ql1=0; end;
if isempty(qi1), qi1=0; end;
if isempty(b1), b1=0; end;
if isempty(pi1), pi1=0; end;
if isempty(thv1), thv1=0; end;
if isempty(qt), qt=0; end;
if isempty(dp), dp=0; end;
if isempty(dz), dz=0; end;
if isempty(ps), ps=0; end;
if isempty(frac), frac=0; end;
if isempty(th2), th2=0; end;
if isempty(p2), p2=0; end;
if isempty(t2), t2=0; end;
if isempty(qv2), qv2=0; end;
if isempty(ql2), ql2=0; end;
if isempty(qi2), qi2=0; end;
if isempty(b2), b2=0; end;
if isempty(pi2), pi2=0; end;
if isempty(thv2), thv2=0; end;
if isempty(thlast), thlast=0; end;
if isempty(fliq), fliq=0; end;
if isempty(fice), fice=0; end;
if isempty(tbar), tbar=0; end;
if isempty(qvbar), qvbar=0; end;
if isempty(qlbar), qlbar=0; end;
if isempty(qibar), qibar=0; end;
if isempty(lhv), lhv=0; end;
if isempty(lhs), lhs=0; end;
if isempty(lhf), lhf=0; end;
if isempty(rm), rm=0; end;
if isempty(cpm), cpm=0; end;
if isempty(avgth), avgth=0; end;
if isempty(avgqv), avgqv=0; end;
%-----------------------------------------------------------------------
if isempty(g), g     = 9.81; end;
if isempty(p00), p00   = 100000.0; end;
if isempty(cp), cp    = 1005.7; end;
if isempty(rd), rd    = 287.04; end;
if isempty(rv), rv    = 461.5; end;
if isempty(xlv), xlv   = 2501000.0; end;
if isempty(xls), xls   = 2836017.0; end;
if isempty(t0), t0    = 273.15; end;
if isempty(cpv), cpv   = 1875.0; end;
if isempty(cpl), cpl   = 4190.0; end;
if isempty(cpi), cpi   = 2118.636; end;
if isempty(lv1), lv1   = xlv+(cpl-cpv).*t0; end;
if isempty(lv2), lv2   = cpl-cpv; end;
if isempty(ls1), ls1   = xls+(cpi-cpv).*t0; end;
if isempty(ls2), ls2   = cpi-cpv; end;
if isempty(rp00), rp00  = 1.0./p00; end;
if isempty(eps), eps   = rd./rv; end;
if isempty(reps), reps  = rv./rd; end;
if isempty(rddcp), rddcp = rd./cp; end;
if isempty(cpdrd), cpdrd = cp./rd; end;
if isempty(cpdg), cpdg  = cp./g; end;
if isempty(converge), converge = 0.2; end;
if isempty(debug_level), debug_level =   0; end;
%-----------------------------------------------------------------------
%---- convert p,t,td to mks units; get pi,q,th,thv ----!
for k=1:nk;
p(k) = 100.0.*p_in(k);
t(k) = 273.15+t_in(k);
td(k) = 273.15+td_in(k);
pi(k) =(p(k).*rp00).^rddcp;
[q(k) ,p(k),td(k)]=getqvs(p(k),td(k));
th(k) = t(k)./pi(k);
thv(k) = th(k).*(1.0+reps.*q(k))./(1.0+q(k));
end; k=fix(nk+1);
%---- get height using the hydrostatic equation ----!
z(1) = 0.0;
for k=2:nk;
dz = -cpdg.*0.5.*(thv(k)+thv(k-1)).*(pi(k)-pi(k-1));
z(k) = z(k-1) + dz;
end; k=fix(nk+1);
%---- find source parcel ----!
if(source==1)
% use surface parcel
kmax = 1;
elseif(source==2);
% use most unstable parcel (max theta-e)
if(p(1)<50000.0)
% first report is above 500 mb ... just use the first level reported
kmax = 1;
[maxthe ,p(1),t(1),td(1),q(1)]=getthe(p(1),t(1),td(1),q(1));
else;
% find max thetae below 500 mb
maxthe = 0.0;
for k=1:nk;
if(p(k)>=50000.0)
[the ,p(k),t(k),td(k),q(k)]=getthe(p(k),t(k),td(k),q(k));
if( the>maxthe )
maxthe = the;
kmax = fix(k);
end;
end;
end; k=fix(nk+1);
end;
if(debug_level>=100)
print .*,'  kmax,maxthe = ',kmax,maxthe;
end;
elseif(source==3);
% use mixed layer
if((z(2)-z(1))>ml_depth )
% the second level is above the mixed-layer depth:  just use the
% lowest level
avgth = th(1);
avgqv = q(1);
kmax = 1;
elseif( z(nk)<ml_depth );
% the top-most level is within the mixed layer:  just use the
% upper-most level
avgth = th(nk);
avgqv = q(nk);
kmax = fix(nk);
else;
% calculate the mixed-layer properties:
avgth = 0.0;
avgqv = 0.0;
k = 2;
if(debug_level>=100)
print .*,'  ml_depth = ',ml_depth;
end;
if(debug_level>=100)
print .*,'  k,z,th,q:';
end;
if(debug_level>=100)
print .*,1,z(1),th(1),q(1);
end;
while((z(k)<=ml_depth) &(k<=nk) );
if(debug_level>=100)
print .*,k,z(k),th(k),q(k);
end;
avgth = avgth + 0.5.*(z(k)-z(k-1)).*(th(k)+th(k-1));
avgqv = avgqv + 0.5.*(z(k)-z(k-1)).*(q(k)+q(k-1));
k = fix(k + 1);
end;
th2 = th(k-1)+(th(k)-th(k-1)).*(ml_depth-z(k-1))./(z(k)-z(k-1));
qv2 =  q(k-1)+( q(k)- q(k-1)).*(ml_depth-z(k-1))./(z(k)-z(k-1));
if(debug_level>=100)
print .*,999,ml_depth,th2,qv2;
end;
avgth = avgth + 0.5.*(ml_depth-z(k-1)).*(th2+th(k-1));
avgqv = avgqv + 0.5.*(ml_depth-z(k-1)).*(qv2+q(k-1));
if(debug_level>=100)
print .*,k,z(k),th(k),q(k);
end;
avgth = avgth./ml_depth;
avgqv = avgqv./ml_depth;
kmax = 1;
end;
if(debug_level>=100)
print .*,avgth,avgqv;
end;
else;
writef(1,['%0.15g \n']);
writef(1,['%s \n'], '  Unknown value for source');
writef(1,['%0.15g \n']);
writef(1,['%s %0.15g \n'], '  source = ',source);
writef(1,['%0.15g \n']);
error(['stop encountered in original fortran code  ',char(10),';']);
end;
%---- define parcel properties at initial location ----!
narea = 0.0;
if((source==1)||(source==2) )
k    = fix(kmax);
th2  = th(kmax);
pi2  = pi(kmax);
p2   = p(kmax);
t2   = t(kmax);
thv2 = thv(kmax);
qv2  = q(kmax);
b2   = 0.0;
elseif( source==3 );
k    = fix(kmax);
th2  = avgth;
qv2  = avgqv;
thv2 = th2.*(1.0+reps.*qv2)./(1.0+qv2);
pi2  = pi(kmax);
p2   = p(kmax);
t2   = th2.*pi2;
b2   = g.*( thv2-thv(kmax) )./thv(kmax);
end;
ql2 = 0.0;
qi2 = 0.0;
qt  = qv2;
cape = 0.0;
cin  = 0.0;
lfc  = 0.0;
doit = true;
cloud = false;
if(adiabat==1||adiabat==2)
ice = false;
else;
ice = true;
end;
t2_orig=t2;    [the ,p2,t2,dumvar4,qv2]=getthe(p2,t2,t2,qv2);    t2(dumvar4~=t2_orig)=dumvar4(dumvar4~=t2_orig);
if(debug_level>=100)
print .*,'  the = ',the;
end;
%---- begin ascent of parcel ----!
if(debug_level>=100)
writef(1,['%s \n'], '  Start loop:');
writef(1,['%s %0.15g %0.15g %0.15g \n'], '  p2,th2,qv2 = ',p2,th2,qv2);
end;
while( doit &(k<nk) );
k = fix(k+1);
b1 =  b2;
dp = p(k-1)-p(k);
if( dp<pinc )
nloop = 1;
else;
nloop = fix(1 + fix( dp./pinc ));
dp = dp./(nloop);
end;
for n=1:nloop;
p1 =  p2;
t1 =  t2;
pi1 = pi2;
th1 = th2;
qv1 = qv2;
ql1 = ql2;
qi1 = qi2;
thv1 = thv2;
p2 = p2 - dp;
pi2 =(p2.*rp00).^rddcp;
thlast = th1;
i = 0;
not_converged = true;
while( not_converged );
i = fix(i + 1);
t2 = thlast.*pi2;
if(ice)
fliq = max(min((t2-233.15)./(273.15-233.15),1.0),0.0);
fice = 1.0-fliq;
else;
fliq = 1.0;
fice = 0.0;
end;
qv2 = min( qt , fliq.*getqvs(p2,t2) + fice.*getqvi(p2,t2) );
qi2 = max( fice.*(qt-qv2) , 0.0 );
ql2 = max( qt-qv2-qi2 , 0.0 );
tbar  = 0.5.*(t1+t2);
qvbar = 0.5.*(qv1+qv2);
qlbar = 0.5.*(ql1+ql2);
qibar = 0.5.*(qi1+qi2);
lhv = lv1-lv2.*tbar;
lhs = ls1-ls2.*tbar;
lhf = lhs-lhv;
rm=rd+rv.*qvbar;
cpm=cp+cpv.*qvbar+cpl.*qlbar+cpi.*qibar;
th2=th1.*exp(  lhv.*(ql2-ql1)./(cpm.*tbar)+lhs.*(qi2-qi1)./(cpm.*tbar)+(rm./cpm-rd./cp).*log(p2./p1) );
if(i>90)
print .*,i,th2,thlast,th2-thlast;
end;
if(i>100)
writef(1,['%0.15g \n']);
writef(1,['%s \n'], '  Error:  lack of convergence');
writef(1,['%0.15g \n']);
writef(1,['%s \n'], '  ... stopping iteration ');
writef(1,['%0.15g \n']);
error(['stop encountered in original fortran code  ',char(10),' 1001;']);
end;
if( abs(th2-thlast)>converge )
thlast=thlast+0.3.*(th2-thlast);
else;
not_converged = false;
end;
end;
% Latest pressure increment is complete.  Calculate some
% important stuff:
if( ql2>=1.0e-10 )
cloud = true;
end;
if(adiabat==1||adiabat==3)
% pseudoadiabat
qt  = qv2;
ql2 = 0.0;
qi2 = 0.0;
elseif(adiabat<=0||adiabat>=5);
writef(1,['%0.15g \n']);
writef(1,['%s \n'], '  Undefined adiabat');
writef(1,['%0.15g \n']);
error(['stop encountered in original fortran code  ',char(10),' 10000;']);
end;
end; n=fix(nloop+1);
thv2 = th2.*(1.0+reps.*qv2)./(1.0+qv2+ql2+qi2);
b2 = g.*( thv2-thv(k) )./thv(k);
dz = -cpdg.*0.5.*(thv(k)+thv(k-1)).*(pi(k)-pi(k-1));
t2_orig=t2;    [the ,p2,t2,dumvar4,qv2]=getthe(p2,t2,t2,qv2);    t2(dumvar4~=t2_orig)=dumvar4(dumvar4~=t2_orig);
% Get contributions to CAPE and CIN:
if((b2>=0.0) &&(b1<0.0) )
% first trip into positive area
ps = p(k-1)+(p(k)-p(k-1)).*(0.0-b1)./(b2-b1);
frac = b2./(b2-b1);
parea =  0.5.*b2.*dz.*frac;
narea = narea-0.5.*b1.*dz.*(1.0-frac);
if(debug_level>=200)
writef(1,['%s %0.15g %0.15g \n'], '      b1,b2 = ',b1,b2);
writef(1,['%s %0.15g %0.15g %0.15g \n'], '      p1,ps,p2 = ',p(k-1),ps,p(k));
writef(1,['%s %0.15g \n'], '      frac = ',frac);
writef(1,['%s %0.15g \n'], '      parea = ',parea);
writef(1,['%s %0.15g \n'], '      narea = ',narea);
end;
cin  = cin  + narea;
narea = 0.0;
elseif((b2<0.0) &&(b1>0.0) );
% first trip into neg area
ps = p(k-1)+(p(k)-p(k-1)).*(0.0-b1)./(b2-b1);
frac = b1./(b1-b2);
parea =  0.5.*b1.*dz.*frac;
narea = -0.5.*b2.*dz.*(1.0-frac);
if(debug_level>=200)
writef(1,['%s %0.15g %0.15g \n'], '      b1,b2 = ',b1,b2);
writef(1,['%s %0.15g %0.15g %0.15g \n'], '      p1,ps,p2 = ',p(k-1),ps,p(k));
writef(1,['%s %0.15g \n'], '      frac = ',frac);
writef(1,['%s %0.15g \n'], '      parea = ',parea);
writef(1,['%s %0.15g \n'], '      narea = ',narea);
end;
elseif( b2<0.0 );
% still collecting negative buoyancy
parea =  0.0;
narea = narea-0.5.*dz.*(b1+b2);
else;
% still collecting positive buoyancy
parea =  0.5.*dz.*(b1+b2);
narea =  0.0;
end;
cape = cape + max(0.0,parea);
if(debug_level>=200)
writef(1,[repmat(['%13.4f'] ,1,5),repmat(' ',1,2),'%1f' ' \n'], p2,b1,b2,cape,cin,cloud);
%format(5(f13.4),2x,l1);
end;
if((p(k)<=10000.0)&&(b2<0.0) )
% stop if b < 0 and p < 100 mb
doit = false;
end;
end;
%---- All done ----!
return;
end %subroutine getcape
%-----------------------------------------------------------------------
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%-----------------------------------------------------------------------
function [getqvsresult,p,t]=getqvs(p,t);
getqvsresult=[];
persistent eps es ; 

;
if isempty(es), es=0; end;
if isempty(eps), eps = 287.04./461.5; end;
es = 611.2.*exp(17.67.*(t-273.15)./(t-29.65));
getqvsresult = eps.*es./(p-es);
csnil=dbstack(1); csnil=csnil(1).name(1)~='@';
if csnil&&~isempty(inputname(2)), assignin('caller','FUntemp',t); evalin('caller',[inputname(2),'=FUntemp;']); end
if csnil&&~isempty(inputname(1)), assignin('caller','FUntemp',p); evalin('caller',[inputname(1),'=FUntemp;']); end
return;
end %function getqvs
%-----------------------------------------------------------------------
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%-----------------------------------------------------------------------
function [getqviresult,p,t]=getqvi(p,t);
getqviresult=[];
persistent eps es ; 

;
if isempty(es), es=0; end;
if isempty(eps), eps = 287.04./461.5; end;
es = 611.2.*exp(21.8745584.*(t-273.15)./(t-7.66));
getqviresult = eps.*es./(p-es);
csnil=dbstack(1); csnil=csnil(1).name(1)~='@';
if csnil&&~isempty(inputname(2)), assignin('caller','FUntemp',t); evalin('caller',[inputname(2),'=FUntemp;']); end
if csnil&&~isempty(inputname(1)), assignin('caller','FUntemp',p); evalin('caller',[inputname(1),'=FUntemp;']); end
return;
end %function getqvi
%-----------------------------------------------------------------------
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%-----------------------------------------------------------------------
function [gettheresult,p,t,td,q]=getthe(p,t,td,q);
gettheresult=[];
persistent tlcl ; 

;
if isempty(tlcl), tlcl=0; end;
if((td-t)>=-0.1 )
tlcl = t;
else;
tlcl = 56.0 +((td-56.0).^(-1) + 0.00125.*log(t./td) ).^(-1);
end;
gettheresult=t.*((100000.0./p).^(0.2854.*(1.0-0.28.*q)) ).*exp(((3376.0./tlcl)-2.54).*q.*(1.0+0.81.*q) );
csnil=dbstack(1); csnil=csnil(1).name(1)~='@';
if csnil&&~isempty(inputname(3)), assignin('caller','FUntemp',td); evalin('caller',[inputname(3),'=FUntemp;']); end
if csnil&&~isempty(inputname(2)), assignin('caller','FUntemp',t); evalin('caller',[inputname(2),'=FUntemp;']); end
if csnil&&~isempty(inputname(4)), assignin('caller','FUntemp',q); evalin('caller',[inputname(4),'=FUntemp;']); end
if csnil&&~isempty(inputname(1)), assignin('caller','FUntemp',p); evalin('caller',[inputname(1),'=FUntemp;']); end
return;
end %function getthe
%-----------------------------------------------------------------------
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%-----------------------------------------------------------------------




function out=writef(fid,varargin)
% function out=writef(fid,varargin)
%  Catches fortran stdout (6) and reroutes in to Matlab's stdout (1)
%  Catches fortran stderr (0) and reroutes in to Matlab's stderr (2)
if isnumeric(fid)
 if fid==6,      out=fprintf(1,varargin{:});
 elseif fid==0,  out=fprintf(2,varargin{:});
 elseif isempty(fid) %% treat empty array like a string array [sethg 2008-03-03]
  out=sprintf(varargin{:});
  if nargin>2 %set the calling var to out
   if ~isempty(inputname(1)), assignin('caller',inputname(1),out); end
  end
 else,           out=fprintf(fid,varargin{:});
 end
elseif ischar(fid)
 out=sprintf(varargin{:});
 if nargin>2 %set the calling var to out
  if ~isempty(inputname(1)), assignin('caller',inputname(1),out); end
 end
else,            out=fprintf(fid,varargin{:});
end
end