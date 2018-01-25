%cape_emanuel.m -- cape function from pcmin_emanuel_ver20141030.m
    %Input vectors indexed from high pressure to low pressure


     function [CAPED,TOB,IFLAG]= cape_emanuel(TP,RP,PP,T,R,P,SIG)
%
%     This function calculates the CAPE of a parcel with pressure PP (mb), 
%       temperature TP (K) and mixing ratio RP (gm/gm) given a sounding
%       of temperature (T in K) and mixing ratio (R in gm/gm) as a function
%       of pressure (P in mb). CAPED is
%       the calculated value of CAPE and TOB is the temperature at the
%       level of neutral buoyancy.  IFLAG is a flag
%       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
%       not run owing to improper sounding (e.g.no water vapor at parcel level).
%       IFLAG=2 indicates that routine did not converge.                 
%-------------------------------------------------------------------------
      ptop=70;   %  Pressure below which sounding is ignored
%------------------------------------------------------------------------      
      Nold=max(size(P));
      N=1;
      for i=Nold:-1:1,
          if P(i) > ptop
              N=max(N,i);
              break
          end    
      end   
      if N < Nold
          P(N+1:Nold)=[];
          T(N+1:Nold)=[];
          R(N+1:Nold)=[];
      end    
      TVRDIF=zeros(1,N);
%
%   ***   Default values   ***
%      
      CAPED=0.0;
      TOB=T(1);
      IFLAG=1;
%
%   ***   Check that sounding is suitable    ***
%
      if RP < 1e-6 || TP < 200
       IFLAG=0;
       return
      end            
%
%   ***   Assign values of thermodynamic constants     ***
%
      CPD=1005.7;
      CPV=1870.0;
%     CL=4190.0;
      CL=2500.0;
      CPVMCL=CPV-CL;
      RV=461.5;
      RD=287.04;
      EPS=RD./RV;
      ALV0=2.501e6;
%
%   ***  Define various parcel quantities, including reversible   ***
%   ***                       entropy, S.                         ***
%                           
      TPC=TP-273.15;
      ESP=6.112*exp(17.67.*TPC/(243.5+TPC));
      EVP=RP*PP/(EPS+RP);
      RH=EVP/ESP;
      RH=min(RH,1.0);
      ALV=ALV0+CPVMCL*TPC;
      S=(CPD+RP*CL)*log(TP)-RD*log(PP-EVP)+...
         ALV*RP./TP-RP*RV*log(RH);            
%
%   ***  Find lifted condensation pressure, PLCL   ***
%     
	CHI=TP/(1669.0-122.0*RH-TP);
	PLCL=PP*(RH^CHI);
%
%   ***  Begin updraft loop   ***
%
	NCMAX=0;
%
	JMIN=1e6;
%    
    for J=1:N,
%
	JMIN=min(JMIN,J);
%
%    ***  Parcel quantities below lifted condensation level   ***
%	 
     if P(J) >= PLCL
	  TG=TP*(P(J)./PP)^(RD/CPD);
	  RG=RP;
%
%   ***   Calculate buoyancy   ***
%  
	  TLVR=TG*(1.+RG/EPS)./(1.+RG);
	  TVRDIF(J)=TLVR-T(J).*(1.+R(J)/EPS)/(1+R(J));
     else
%
%   ***  Parcel quantities above lifted condensation level  ***
%	 
	  TGNEW=T(J);          
	  TJC=T(J)-273.15; 
	  ES=6.112*exp(17.67*TJC/(243.5+TJC));
	  RG=EPS*ES/(P(J)-ES);
%
%   ***  Iteratively calculate lifted parcel temperature and mixing   ***
%   ***                ratio for reversible ascent                    ***
%
	  NC=0;
      TG=0.0;
%    
      while (abs(TGNEW-TG)) > 0.001
%
	   TG=TGNEW;
	   TC=TG-273.15;
	   ENEW=6.112*exp(17.67*TC./(243.5+TC));
       RG=EPS*ENEW/(P(J)-ENEW);   
%          
	  NC=NC+1;
%
%   ***  Calculate estimates of the rates of change of the entropy    ***
%   ***           with temperature at constant pressure               ***
%  
	  ALV=ALV0+CPVMCL*(TG-273.15);
	  SL=(CPD+RP*CL+ALV*ALV*RG./(RV*TG*TG))/TG;
	  EM=RG*P(J)/(EPS+RG);
	  SG=(CPD+RP*CL)*log(TG)-RD*log(P(J)-EM)+ ...
          ALV*RG/TG;
          if NC < 3
	   AP=0.3;
          else
	   AP=1.0;
          end
	  TGNEW=TG+AP*(S-SG)/SL;  
%
%   ***   Bail out if things get out of hand   ***
%
       if NC > 500 || ENEW > (P(J)-1)
            IFLAG=2;
            return
       end
%       
       end
%       
       NCMAX=max(NC,NCMAX);
%
%   *** Calculate buoyancy   ***
%
      RMEAN=SIG*RG+(1-SIG)*RP;
	  TLVR=TG*(1.+RG/EPS)/(1.+RMEAN);
	  TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J));
     end
    end
%
%  ***  Begin loop to find NA, PA, and CAPE from reversible ascent ***
%
	NA=0.0;
	PA=0.0;
%
%   ***  Find maximum level of positive buoyancy, INB    ***
%
	INB=1;
    for J=N:-1:JMIN;
        if TVRDIF(J) > 0
            INB=max(INB,J);
        end
    end    
    if INB == 1
        return
    end    
%
%   ***  Find positive and negative areas and CAPE  ***
%
    if INB > 1
        for J=(JMIN+1):INB
            PFAC=RD*(TVRDIF(J)+TVRDIF(J-1))*(P(J-1)-P(J))/(P(J)+P(J-1));
            PA=PA+max(PFAC,0.0);
            NA=NA-min(PFAC,0.0);
        end    
  
%   ***   Find area between parcel pressure and first level above it ***
%
        PMA=(PP+P(JMIN)) ;
        PFAC=RD*(PP-P(JMIN))/PMA;
        PA=PA+PFAC*max(TVRDIF(JMIN),0.0);
        NA=NA-PFAC*min(TVRDIF(JMIN),0.0);
%
%   ***   Find residual positive area above INB and TO  ***
%
       PAT=0.0;
       TOB=T(INB);
       if INB < N
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/ ...
         (TVRDIF(INB)-TVRDIF(INB+1));
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB);
	    TOB=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/ ...
         (P(INB)-P(INB+1));
       end
%
%   ***   Find CAPE  ***
%            
	 CAPED=PA+PAT-NA;
	 CAPED=max(CAPED,0.0);
    end