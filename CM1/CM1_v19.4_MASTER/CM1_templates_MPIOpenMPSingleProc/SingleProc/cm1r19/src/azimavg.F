  MODULE azimavg_module

  implicit none

  private
  public :: azimavg

  integer, parameter :: varmax = 1000

  integer, parameter :: fnums = 67
  integer, parameter :: fnumw = 68

  CONTAINS

    subroutine azimavg(mtime,nwritea,arecs,arecw,qname,dt,dosfcflx,            &
                   icrs,icenter,jcenter,xcenter,ycenter,                       &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,            &
                   yh,vh,rvh,yf,vf,rvf,                                        &
                   xfref,yfref,rds,sigma,rdsf,sigmaf,                          &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,                   &
                   rho0s,pi0s,prs0s,rth0s,                                     &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,u0,v0,                  &
                   qi0,rr0,rf0,rrf0,                                           &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,        &
                   tsk,znt,ust,thflux,qvflux,cd,ch,cq,u1,v1,s1,xland,psfc,tlh, &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                    &
                   divx,rho,rf,prs,t11,t12,t13,t22,t23,t33,                    &
                   rru,u3d,uten,us   ,                                         &
                   rrv,v3d,vten,vs   ,                                         &
                   rrw,w3d,wten,ws   ,                                         &
                   pp3d,ppten,sten,th3d,thv ,thten,thten1,                     &
                   q3d,qten,kmh,kmv,khh,khv,tkea,tke3d,tketen,                 &
                   nm,defv,defh,dissten,                                       &
                   rain,u10,v10,s10,                                           &
                   tdiag,qdiag,udiag,vdiag,wdiag,pdiag,out3d,getdbz,getvt,     &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,flag)
    use input
    use constants
    use cm1libs , only : rslf,rsif
#ifdef MPI
    use mpi
#endif
    implicit none

    double precision, intent(in) :: mtime
    integer, intent(inout) :: nwritea,arecs,arecw
    character(len=3), intent(in), dimension(maxq) :: qname
    real, intent(inout) :: dt
    logical, intent(in) :: dosfcflx
    integer, intent(in) :: icrs
    integer, intent(inout) :: icenter,jcenter
    real,    intent(inout) :: xcenter,ycenter
    real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
    real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
    real, intent(in), dimension(jb:je) :: yh,vh,rvh
    real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
    real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
    real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
    real, intent(in), dimension(kb:ke) :: rds,sigma
    real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh,c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
    real, intent(in), dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s,rth0s
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: qi0,rr0,rf0,rrf0
    real, intent(in), dimension(ib:ie,jb:je) :: zs
    real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
    real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
    real, intent(in), dimension(ib:ie,jb:je) :: tsk,znt,ust,thflux,qvflux,cd,ch,cq,  &
                                                u1,v1,s1,xland,psfc,tlh
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: divx,rho,rf,prs
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
    real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d,uten,us
    real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,v3d,vten,vs
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,w3d,wten,ws
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: pp3d,ppten,sten
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: th3d,thv,thten,thten1
    real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d,qten
    real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
    real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh,dissten
    real, intent(in), dimension(ib:ie,jb:je,nrain) :: rain
    real, intent(in), dimension(ibl:iel,jbl:jel) :: u10,v10,s10
    real, intent(in) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
    real, intent(in) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
    real, intent(in) , dimension(ibdp:iedp,jbdp:jedp,kbdp:kedp,npdiag) :: pdiag
    real, intent(in) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
    logical, intent(in) :: getdbz,getvt
    real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
    real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
    logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag

    !---------------------------------------------------------------------------

    integer :: i,j,k,n,nn,nq,nvar,nvar2d,nfile,ntmp
    integer :: ip,im
    real :: rddr
    character(len=80) :: newstring
    character(len=80), dimension(varmax) :: varname,vardesc,varunit
    integer, dimension(varmax) :: varlvls
    character(len=1),  dimension(varmax) :: vargrid
    character(len=80) :: a1,a2
    character(len=16) :: a16

    integer, dimension(icrs) :: navg
    real, dimension(icrs) :: ravg
    double precision, dimension(icrs) :: savg2d
    double precision, dimension(icrs,nk+1) :: savg
    real, dimension(icrs,nk) :: uavg,vavg,wavg,wsavg,pavg,rhoavg,thavg,qvavg,mavg,thvavg,ppiavg
    integer, dimension(ni,nj) :: ir
    real, dimension(ni,nj) :: rr,angle
    real :: pisfc,prssfc,thsfc,qvsfc,ssfc,pd,p1,p2,p3,tx,qx,px,tlcl,ee,ql,qi,qvs
    real :: ff,qvsl,qvsi,qt
    logical :: doit

    if( myid.eq.0 ) print *,' ..... begin azimavg code ..... '

    rddr = 1.0/ddr
    nvar = 0

    ! by default, assume 3d and scalar levels:
    varlvls = nk
    vargrid = 's'

    savg2d = 0.0
    savg = 0.0

    !---------------------------------------------------------------------------
    !  Get center point:

    call       getcenter(rddr,xh,yh,zh,xfref,yfref,u3d,us  ,v3d,vs  ,  &
                         nwritea,icenter,jcenter,xcenter,ycenter)

    !---------------------------------------------------------------------------
    ! open file:

    IF( myid.eq.0 )THEN

        print *,'  nwritea = ',nwritea

        ! grads-format
        ! (by default, for now)

        newstring = string

        newstring(totlen+1:totlen+14) = '_azimavg_s.dat'
        print *,newstring
        open(unit=fnums,file=newstring,form='unformatted',access='direct',recl=4*icrs)

        newstring(totlen+1:totlen+14) = '_azimavg_w.dat'
        print *,newstring
        open(unit=fnumw,file=newstring,form='unformatted',access='direct',recl=4*icrs)

    ENDIF

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      navg = 0.0

      do j=1,nj
      do i=1,ni
        rr(i,j) = sqrt( (xh(i)-xcenter)**2 + (yh(j)-ycenter)**2 )
        ir(i,j) = 1 + int( rr(i,j)*rddr )
        if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
          navg(ir(i,j)) = navg(ir(i,j)) + 1
        endif
        if( (xh(i)-xcenter).ge.0.0 )then
          angle(i,j) =      asin( (yh(j)-ycenter)/max(smeps,rr(i,j)) )
        else
          angle(i,j) = pi - asin( (yh(j)-ycenter)/max(smeps,rr(i,j)) )
        endif
      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE,navg,icrs,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

      ravg = 0.0

      do i=1,icrs
        ravg(i) = 0.5*ddr + (i-1)*ddr
      enddo


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !  2d variables:

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'navg'
      vardesc(nvar) = 'number of data points in each bin'
      varunit(nvar) = 'count'
      varlvls(nvar) = 0

      do i=1,icrs
        savg2d(i) = navg(i)
      enddo

      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ravg'
      vardesc(nvar) = 'avg radius'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = rr(i,j)
      enddo
      enddo

      call     getavg2d(dum1(ib,jb,1),icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ravg2'
      vardesc(nvar) = 'actual radius'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      do i=1,icrs
        savg2d(i) = ravg(i)
      enddo

      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'cd'
      vardesc(nvar) = 'surface drag coefficient'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(cd,icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ch'
      vardesc(nvar) = 'surface sensible heat exchange coeff'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(ch,icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'cq'
      vardesc(nvar) = 'surface latent heat exchange coeff'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(cq,icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'znt'
      vardesc(nvar) = 'surface roughness length'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      call     getavg2d(znt,icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ust'
      vardesc(nvar) = 'surface friction velocity'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(ust,icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'thflux'
      vardesc(nvar) = 'surface potential temperature flux'
      varunit(nvar) = 'K m/s'
      varlvls(nvar) = 0

      call     getavg2d(thflux,icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'qvflux'
      vardesc(nvar) = 'surface water vapor flux'
      varunit(nvar) = 'g/g m/s'
      varlvls(nvar) = 0

      call     getavg2d(qvflux,icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      do j=1,nj
      do i=1,ni
        p1 = pi0(1,1,1)+pp3d(i,j,1)
        p2 = pi0(1,1,2)+pp3d(i,j,2)
        p3 = pi0(1,1,3)+pp3d(i,j,3)
        pisfc=cgs1*p1+cgs2*p2+cgs3*p3
        prssfc = p00*(pisfc**cpdrd)
        thsfc = tsk(i,j)/pisfc
        qvsfc = rslf(prssfc,tsk(i,j))
        pd = prssfc/(1.0+qvsfc*reps)
        ssfc = cp*alog(tsk(i,j))   &
             + 2555000.0*qvsfc/tsk(i,j)   &
             - rd*alog(pd)
        dum1(i,j,1) = pisfc
        dum1(i,j,2) = prssfc
        dum1(i,j,3) = thsfc
        dum1(i,j,4) = qvsfc
        dum1(i,j,5) = ssfc
      enddo
      enddo

      nvar = nvar+1
      varname(nvar) = 'pisfc'
      vardesc(nvar) = 'surface nondimensional pressure'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,1),icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'prssfc'
      vardesc(nvar) = 'surface pressure'
      varunit(nvar) = 'Pa'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,2),icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'tsk'
      vardesc(nvar) = 'surface (ocean/land) temperature'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(tsk,icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'thsfc'
      vardesc(nvar) = 'surface (ocean/land) potential temperature'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,3),icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'qvsfc'
      vardesc(nvar) = 'surface (ocean/land) water vapor mixing ratio'
      varunit(nvar) = 'g/g'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,4),icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ssfc'
      vardesc(nvar) = 'surface (ocean/land) moist entropy'
      varunit(nvar) = 'J/kg/K'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,5),icrs,savg2d,ir,rr,navg)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar2d = nvar

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! 3d vars:  averages

              ! dum1 = u at s pts
              ! dum2 = v at s pts
              ! dum3 = ur at s pts
              ! dum4 = vt at s pts
              ! dum5 = w at s pts

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = 0.5*(u3d(i,j,k)+u3d(i+1,j,k))
          dum2(i,j,k) = 0.5*(v3d(i,j,k)+v3d(i,j+1,k))
          if( abs(rr(i,j)).lt.1.0e-4 )then
            dum3(i,j,k)=0.0
            dum4(i,j,k)=0.0
          else
            dum3(i,j,k)=dum2(i,j,k)*sin(angle(i,j))+dum1(i,j,k)*cos(angle(i,j))
            dum4(i,j,k)=dum2(i,j,k)*cos(angle(i,j))-dum1(i,j,k)*sin(angle(i,j))
          endif
          dum5(i,j,k) = 0.5*(w3d(i,j,k)+w3d(i,j,k+1))
          dum6(i,j,k) = sqrt( dum1(i,j,k)**2 + dum2(i,j,k)**2 )
        enddo
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'u'
        vardesc(nvar) = 'radial velocity'
        varunit(nvar) = 'm/s'

        call getavg3d(dum3,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save uavg:
        do k=1,nk
        do i=1,icrs
          uavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'v'
        vardesc(nvar) = 'tangential velocity'
        varunit(nvar) = 'm/s'

        call getavg3d(dum4,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save vavg:
        do k=1,nk
        do i=1,icrs
          vavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'm'
        vardesc(nvar) = 'angular momentum'
        varunit(nvar) = 'm^2/s'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = rr(i,j)*dum4(i,j,k) + 0.5*fcor*rr(i,j)*rr(i,j)
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save mavg:
        do k=1,nk
        do i=1,icrs
          mavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'wsp'
        vardesc(nvar) = 'horizontal wind speed'
        varunit(nvar) = 'm/s'

        call getavg3d(dum6,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        ! w interpolated to scalar levels:
        nvar = nvar+1
        varname(nvar) = 'w'
        vardesc(nvar) = 'vertical velocity'
        varunit(nvar) = 'm/s'

        call getavg3d(dum5,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save wavg:
        do k=1,nk
        do i=1,icrs
          wavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        ! w at staggered level:
        nvar = nvar+1
        varname(nvar) = 'w'
        vardesc(nvar) = 'vertical velocity'
        varunit(nvar) = 'm/s'
        vargrid(nvar) = 'w'   ! w levels

        call getavg3d(w3d(ib,jb,kb),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save wsavg:
        do k=1,nk
        do i=1,icrs
          wsavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'rho'
        vardesc(nvar) = 'dry air density'
        varunit(nvar) = 'kg/m3'

        call getavg3d(rho,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save rhoavg:
        do k=1,nk
        do i=1,icrs
          rhoavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'prs'
        vardesc(nvar) = 'pressure'
        varunit(nvar) = 'Pa'

        call getavg3d(prs,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save pavg:
        do k=1,nk
        do i=1,icrs
          pavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'ppi'
        vardesc(nvar) = 'nondimensional pressure perturbation'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pp3d,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save ppiavg:
        do k=1,nk
        do i=1,icrs
          ppiavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'th'
        vardesc(nvar) = 'potential temperature'
        varunit(nvar) = 'K'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = th0(i,j,k)+th3d(i,j,k)
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save thavg:
        do k=1,nk
        do i=1,icrs
          thavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 't'
        vardesc(nvar) = 'temperature'
        varunit(nvar) = 'K'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'rh'
        vardesc(nvar) = 'relative humidity wrt liquid'
        varunit(nvar) = 'nondimensional'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = q3d(i,j,k,nqv)/rslf(prs(i,j,k),dum1(i,j,k))
        enddo
        enddo
        enddo

        call getavg3d(dum2,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. iice.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'rhi'
        vardesc(nvar) = 'relative humidity wrt ice'
        varunit(nvar) = 'nondimensional'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = q3d(i,j,k,nqv)/rsif(prs(i,j,k),dum1(i,j,k))
        enddo
        enddo
        enddo

        call getavg3d(dum2,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qv'
        vardesc(nvar) = 'water vapor mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqv),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save qvavg:
        do k=1,nk
        do i=1,icrs
          qvavg(i,k) = savg(i,k)
        enddo
        enddo

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqc.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qc'
        vardesc(nvar) = 'cloud water mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqc),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqr.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qr'
        vardesc(nvar) = 'rain water mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqr),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqi.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qi'
        vardesc(nvar) = 'cloud ice mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqi),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqs.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qs'
        vardesc(nvar) = 'snow mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqs),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqg.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qg'
        vardesc(nvar) = 'graupel/hail mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqg),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. output_dbz.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'dbz'
        vardesc(nvar) = 'reflectivity'
        varunit(nvar) = 'dBZ'

        call getavg3d(qdiag(ibdq,jbdq,kbdq,qd_dbz),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'thv'
        vardesc(nvar) = 'virtual potential temperature'
        varunit(nvar) = 'K'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          ql = 0.0
          do nq=nql1,nql2
            ql=ql+q3d(i,j,k,nq)
          enddo
          qi = 0.0
          if(iice.eq.1)then
          do nq=nqs1,nqs2
            qi=qi+q3d(i,j,k,nq)
          enddo
          endif
          dum2(i,j,k) = (th0(i,j,k)+th3d(i,j,k))*(1.0+reps*q3d(i,j,k,nqv)) &
                                                /(1.0+q3d(i,j,k,nqv)+ql+qi)
        enddo
        enddo
        enddo

        call getavg3d(dum2,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save thvavg:
        do k=1,nk
        do i=1,icrs
          thvavg(i,k) = savg(i,k)
        enddo
        enddo

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'the'
        vardesc(nvar) = 'equivalent potential temperature'
        varunit(nvar) = 'K'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          tx = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
          qx = q3d(i,j,k,nqv)
          px = p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
          if(q3d(i,j,k,nqc).ge.clwsat)then
            tlcl=tx
          else
            ee=0.01*px*qx/(eps+qx)
            tlcl=55.0+2840.0/(3.5*alog(tx)-log(1.0e-20+ee)-4.805)
          endif
          dum1(i,j,k) = tx*((p00/px)**(0.2854*(1.0-0.28*qx)))   &
                          *exp(((3376.0/tlcl)-2.54)*qx*(1.0+0.81*qx))
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'satfrac'
        vardesc(nvar) = 'fraction of points that are saturated'
        varunit(nvar) = 'nondimensional'

        dum1 = 0.0

        IF( iice.eq.0 )THEN
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tx = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
            qvs = rslf(prs(i,j,k),tx)
            ql = max( q3d(i,j,k,nqv)+q3d(i,j,k,nqc) - qvs , 0.0 )
            if( ql.gt.1.0e-10 ) dum1(i,j,k) = 1.0
          enddo
          enddo
          enddo
        ELSE
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tx = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
            qt = q3d(i,j,k,nqv)+q3d(i,j,k,nqc)+q3d(i,j,k,nqi)
            if( tx.ge.273.15 )then
              qvs = rslf(prs(i,j,k),tx)
              ql = max( qt-qvs , 0.0 )
              if( ql.gt.1.0e-10 ) dum1(i,j,k) = 1.0
            elseif( tx.le.233.15 )then
              qvs = rsif(prs(i,j,k),tx)
              qi = max( qt-qvs , 0.0 )
              if( qi.gt.1.0e-10 ) dum1(i,j,k) = 1.0
            else
              ff = (tx-233.15)/(273.15-233.15)
              qvsl = rslf(prs(i,j,k),tx)
              qvsi = rsif(prs(i,j,k),tx)
              qvs = ff*qvsl + (1.0-ff)*qvsi
              qi = max( 0.0 , (1.0-ff)*(qt-qvs) )
              ql = max( 0.0 , ff*(qt-qvs) )
              if( (ql+qi).gt.1.0e-10 ) dum1(i,j,k) = 1.0
            endif
          enddo
          enddo
          enddo
        ENDIF

        call getavg3d(dum1,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'vgrad'
        vardesc(nvar) = 'gradient wind speed (from azim-avg fields)'
        varunit(nvar) = 'm/s'

        do k=1,nk
        do i=1,icrs
          ip = min( i+1 , icrs )
          im = max( i-1 , 1 )
          
          savg(i,k) = -0.5*fcor*ravg(i) + sqrt( max(0.0,               &
                            0.25*fcor*fcor*ravg(i)*ravg(i)               &
           +ravg(i)*cp*thvavg(i,k)*(ppiavg(ip,k)-ppiavg(im,k))/(ravg(ip)-ravg(im))   &
                                           ) )
        enddo
        enddo

        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      IF( cm1setup.eq.1 .and. sgsmodel.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'tke'
        vardesc(nvar) = 'subgrid turbulence kinetic energy'
        varunit(nvar) = 'm^2/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(tke3d(ibt,jbt,kbt),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( iusekm.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'kmh'
        vardesc(nvar) = 'eddy viscosity, horiz direction'
        varunit(nvar) = 'm^2/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(kmh(ibc,jbc,kbc),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( iusekm.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'kmv'
        vardesc(nvar) = 'eddy viscosity, vert direction'
        varunit(nvar) = 'm^2/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(kmv(ibc,jbc,kbc),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( iusekh.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'khh'
        vardesc(nvar) = 'eddy diffusivity, horiz direction'
        varunit(nvar) = 'm^2/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(khh(ibc,jbc,kbc),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( iusekh.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'khv'
        vardesc(nvar) = 'eddy diffusivity, vert direction'
        varunit(nvar) = 'm^2/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(khv(ibc,jbc,kbc),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( pdcomp )THEN

        nvar = nvar+1
        varname(nvar) = 'pipb'
        vardesc(nvar) = 'diagnosed pi-prime: buoyancy component'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pdiag(ib,jb,kb,1),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

        nvar = nvar+1
        varname(nvar) = 'pipdl'
        vardesc(nvar) = 'diagnosed pi-prime: linear dynamic component'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pdiag(ib,jb,kb,2),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

        nvar = nvar+1
        varname(nvar) = 'pipdn'
        vardesc(nvar) = 'diagnosed pi-prime: nonlinear dynamic component'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pdiag(ib,jb,kb,3),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

      if( icor.eq.1 )then
        nvar = nvar+1
        varname(nvar) = 'pipc'
        vardesc(nvar) = 'diagnosed pi-prime: Coriolis component'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pdiag(ib,jb,kb,4),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
      endif

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      if( wd_buoy.ge.1 )then

        nvar = nvar+1
        varname(nvar) = 'wb_buoy'
        vardesc(nvar) = 'w budget: buoyancy'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

        call getavg3d(wdiag(ib,jb,kb,wd_buoy),icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      endif

      !c-c-c-c-c-c-c-c-c-c

      IF( pdcomp )THEN

        nvar = nvar+1
        varname(nvar) = 'pgradb'
        vardesc(nvar) = 'vert pres grad: buoyancy component'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do j=1,nj
          do i=1,ni
            wten(i,j,1) = 0.0
            wten(i,j,nk+1) = 0.0
          enddo
          enddo

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,1)-pdiag(i,j,k-1,1))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

        nvar = nvar+1
        varname(nvar) = 'pgraddl'
        vardesc(nvar) = 'vert pres grad: linear dynamic component'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,2)-pdiag(i,j,k-1,2))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

        nvar = nvar+1
        varname(nvar) = 'pgraddn'
        vardesc(nvar) = 'vert pres grad: nonlinear dynamic component'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,3)-pdiag(i,j,k-1,3))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

      if( icor.eq.1 )then
        nvar = nvar+1
        varname(nvar) = 'pgradc'
        vardesc(nvar) = 'vert pres grad: Coriolis component'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,4)-pdiag(i,j,k-1,4))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
      endif

        !-----

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( pdcomp .and. wd_buoy.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'btot'
        vardesc(nvar) = 'total buoyant acceleration'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,1)-pdiag(i,j,k-1,1))*rdz*mf(i,j,k)  &
                          +wdiag(i,j,k,wd_buoy)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      endif

      !c-c-c-c-c-c-c-c-c-c

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! close output file:

    id0:  IF( myid.eq.0 )THEN

        ! grads-format
        ! (by default, for now)

        ! close binary output file:
        close(unit=fnums)
        close(unit=fnumw)

        ! write descriptor file:
        do nfile=1,2
          newstring = string
          if( nfile.eq.1 )then
            newstring(totlen+1:totlen+14) = '_azimavg_s.ctl'
          elseif( nfile.eq.2 )then
            newstring(totlen+1:totlen+14) = '_azimavg_w.ctl'
          endif
          open(unit=66,file=newstring)
          if( nfile.eq.1 )then
            sstring(baselen+1:baselen+1+14) = '_azimavg_s.dat'
          elseif( nfile.eq.2 )then
            sstring(baselen+1:baselen+1+14) = '_azimavg_w.dat'
          endif
          write(66,101) sstring
          if( nfile.eq.1 )then
            write(66,103)
          elseif( nfile.eq.2 )then
            write(66,203)
          endif
          write(66,104)
          write(66,105) icrs,0.001*0.5*ddr,0.001*ddr
          write(66,106)
          if( nfile.eq.1 )then
            write(66,107) nk
            do k=1,nk
              write(66,217) 0.001*zh(1,1,k)
            enddo
          elseif( nfile.eq.2 )then
            write(66,107) nk+1
            do k=1,nk+1
              write(66,217) 0.001*zf(1,1,k)
            enddo
          endif
          write(66,108) nwritea
          ntmp = 0
          if( nfile.eq.1 )then
            do n=1,nvar
              if( vargrid(n).eq.'s' ) ntmp = ntmp+1
            enddo
          else
            do n=1,nvar
              if( vargrid(n).eq.'w' ) ntmp = ntmp+1
            enddo
          endif
          write(66,109) ntmp
          do n=1,nvar
            doit = .false.
            if( nfile.eq.1 .and. vargrid(n).eq.'s' ) doit = .true.
            if( nfile.eq.2 .and. vargrid(n).eq.'w' ) doit = .true.
            if( doit )then
              a1 = varname(n)
              a2 = vardesc(n)
              !---
              a16 = '                '
              nn = len(trim(varunit(n)))
              write(a16(2:15),314) varunit(n)
              write(a16(1:1),301 )       '('
              write(a16(nn+2:nn+2),301 ) ')'
              !---
              if(     nfile.eq.1 )then
                write(66,110) a1(1:12),varlvls(n),a2(1:40),a16
              elseif( nfile.eq.2 )then
                write(66,110) a1(1:12),(varlvls(n)+1),a2(1:40),a16
              endif
            endif
          enddo
          write(66,111)
          close(unit=66)
        enddo

 301    format(a1)
 314    format(a14)

 101    format('dset ^',a70)
 102    format('options template')
 103    format('title CM1 azimuthally averaged output, scalar levels')
 203    format('title CM1 azimuthally averaged output, w levels')
 104    format('undef -99999999.')
 105    format('xdef ',i6,' linear ',f13.6,1x,f13.6)
 106    format('ydef 1 linear 0 1')
 107    format('zdef ',i6,' levels')
 217    format(2x,f13.6)
 108    format('tdef ',i10,' linear 00:00Z03JUL0001 1YR')
 109    format('vars ',i6)
 110    format(a12,2x,i6,' 99 ',a40,1x,a16)
 111    format('endvars')


    ENDIF  id0

    if( myid.eq.0 ) print *,' ..... end azimavg code ..... '

    end subroutine azimavg


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getavg2d(s2d,icrs,savg2d,ir,rr,navg)
    use input
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je) :: s2d
    integer, intent(in) :: icrs
    double precision, intent(inout), dimension(icrs) :: savg2d
    integer, intent(in), dimension(ni,nj) :: ir
    real, intent(in), dimension(ni,nj) :: rr
    integer, intent(in), dimension(icrs) :: navg

    integer :: i,j

    savg2d = 0.0

    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg2d(ir(i,j)) = savg2d(ir(i,j)) + s2d(i,j)
      endif
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg2d,icrs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    do i=1,icrs
      savg2d(i) = savg2d(i)/dble(max(1,navg(i)))
    enddo
   
    end subroutine getavg2d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getavg3d(s,icrs,savg,ir,rr,navg)
    use input
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: s
    integer, intent(in) :: icrs
    double precision, dimension(icrs,nk+1) :: savg
    integer, intent(in), dimension(ni,nj) :: ir
    real, intent(in), dimension(ni,nj) :: rr
    integer, intent(in), dimension(icrs) :: navg

    integer :: i,j,k

    savg = 0.0

    do k=1,nk+1
    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg(ir(i,j),k) = savg(ir(i,j),k) + s(i,j,k)
      endif
    enddo
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg,icrs*(nk+1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    do k=1,nk+1
    do i=1,icrs
      savg(i,k) = savg(i,k)/dble(max(1,navg(i)))
    enddo
    enddo

    end subroutine getavg3d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    use input
    implicit none

    integer, intent(in) :: icrs
    double precision, intent(in), dimension(icrs) :: savg2d
    integer, intent(inout) :: arecs,arecw,nvar
    character(len=80), intent(in), dimension(varmax) :: varname,vardesc,varunit
    character(len=1), intent(in), dimension(varmax) :: vargrid

    integer :: i
      
    if( myid.eq.0 )then

      print *,nvar,trim(varname(nvar))

        ! grads-format file:

        ! note:  all 2d vars go into scalar output file
        write(fnums,rec=arecs) (sngl(savg2d(i)),i=1,icrs)
        arecs = arecs+1

    endif

    end subroutine write2d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    use input
    implicit none

    integer, intent(in) :: icrs
    double precision, intent(in), dimension(icrs,nk+1) :: savg
    integer, intent(inout) :: arecs,arecw,nvar
    character(len=80), intent(in), dimension(varmax) :: varname,vardesc,varunit
    character(len=1), intent(in), dimension(varmax) :: vargrid

    integer :: i,k
      
    if( myid.eq.0 )then

      print *,nvar,vargrid(nvar),' ',trim(varname(nvar))

        ! grads-format file:

      if( vargrid(nvar).eq.'s' )then
        do k=1,nk
          write(fnums,rec=arecs) (sngl(savg(i,k)),i=1,icrs)
          arecs = arecs+1
        enddo
      elseif( vargrid(nvar).eq.'w' )then
        do k=1,nk+1
          write(fnumw,rec=arecw) (sngl(savg(i,k)),i=1,icrs)
          arecw = arecw+1
        enddo
      else
        print *,' 67833 '
        call stopcm1
      endif

    endif

    end subroutine write3d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getcenter(rddr,xh,yh,zh,xfref,yfref,u3d,utmp,v3d,vtmp,  &
                         nwritea,icenter,jcenter,xcenter,ycenter)
    use input
    use constants
#ifdef MPI
      use mpi
#endif
    implicit none

      !----------------------------------------------------------------
      ! for center-finding code:
      ! (the code searches for the point that yields maximum possible
      !  azimuthally averaged azimuthal velocity, and defines this
      !  point as "storm center")
      !----------------------------------------------------------------
      !  other schemes can be coded in the future
      !----------------------------------------------------------------
      ! (note: for axisymmetric simulations, this code isn't used)
      !----------------------------------------------------------------

      real, intent(in) :: rddr
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, intent(in   ), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: utmp
      real, intent(in   ), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vtmp
      integer, intent(inout) :: nwritea,icenter,jcenter
      real, intent(inout) :: xcenter,ycenter

    integer :: i,j,k,ii,jj,ir,proc,ifoo,jfoo,myifoo,myjfoo,arecs,arecw,irec1,nr,kmax
    integer :: ictest,jctest,imin,imax,jmin,jmax,ninc
    real :: rr,vmax,mostmax,mostx,mosty,mostr,mostz,rmax,zmax,u,v,angle
    real :: xctest,yctest

    logical, parameter :: doit = .true.

#ifdef MPI
    integer status(MPI_STATUS_SIZE)
#endif

             ! assume rmw <= rcrit
             ! assume zmw <= zcrit
      real, parameter :: rcrit = 100000.0       ! max possible rmw (m)
      real, parameter :: zcrit =   2000.0       ! max possible zmw (m)
                                                ! (Note:  code uses last known location 
                                                !  of storm center as a starting point)

      real, parameter :: search_radius  =  50000.0   ! to reduce cost, search only within 
                                                     ! this radius (m) from previous known 
                                                     ! position

      integer, dimension(:), allocatable :: navg
      double precision, dimension(:), allocatable :: ravg
      double precision, dimension(:,:), allocatable :: davg
      character(len=80) :: newstring
      character(len=70) :: rstring

      integer, parameter :: nvars = 7

!-----------------------------------------------------------------------

  docenter:  IF(doit)THEN

    if(myid.eq.0) print *,'  Entering getcenter .... '

      nr = int( rcrit*rddr + 1 )

      kmax = 0
      do while( zh(1,1,kmax).le.zcrit )
        kmax = kmax+1
      enddo

    ! print xc,yc from last time:
    if(myid.eq.0) print *,'  Previous values: '
    if(myid.eq.0) print *,'    xcenter,ycenter = ',xcenter,ycenter
    if(myid.eq.0) print *,'    icenter,jcenter = ',icenter,jcenter
    if(myid.eq.0) print *,'  nr,kmax           = ',nr,kmax

    allocate( navg(nr) )
    navg = 0
    allocate( ravg(nr) )
    ravg = 0
    allocate( davg(nr,kmax) )
    davg = 0.0

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
    do k=1,kmax
    do j=1,nj
    do i=1,ni
      utmp(i,j,k) = 0.5*( u3d(i,j,k)+u3d(i+1,j,k) )
      vtmp(i,j,k) = 0.5*( v3d(i,j,k)+v3d(i,j+1,k) )
    enddo
    enddo
    enddo

      mostmax = 0.0
      mostx = 0.0
      mosty = 0.0
      mostr = 0.0
      mostz = 0.0
      ictest = nx/2 + 1
      jctest = ny/2 + 1

      ifoo = int(search_radius*rddr) + 1

      imin = icenter - ifoo
      imax = icenter + ifoo

      jmin = jcenter - ifoo
      jmax = jcenter + ifoo

      ! to reduce cost, skip some grid points:
      !   (set ninc=1 to search every gridpoint)
      !   (set ninc=2 to search every other gridpoint)
      !   (set ninc=3 to search every third gridpoint)
      ninc  =  1

    if(myid.eq.0)then
      print *,'  imin,imax,ninc = ',imin,imax,ninc
      print *,'  jmin,jmax,ninc = ',jmin,jmax,ninc
    endif

    if( imin.lt.1 .or. imax.gt.nx .or.  &
        jmin.lt.1 .or. jmax.gt.ny )then
      print *,'  56091 '
      call stopcm1
    endif

      jloop:  do jj=jmin,jmax,ninc
      iloop:  do ii=imin,imax,ninc
        navg = 0
        ravg = 0
        davg = 0.0
        xctest = xfref(ii)
        yctest = yfref(jj)
        checkit: IF( sqrt( (xfref(ii)-xcenter)**2 + (yfref(jj)-ycenter)**2 ) .le. search_radius )THEN
          ! Calculate azim-avg vt:
          do j=1,nj
          do i=1,ni
            rr = sqrt( (xh(i)-xctest)**2 + (yh(j)-yctest)**2 )
            if( rr.lt.rcrit )then
              ir = 1 + int( rr*rddr )
              navg(ir)=navg(ir)+1
              if( abs(rr).lt.1.0e-4 )then
                u=0.0
                v=0.0
              else
                if( (xh(i)-xctest).ge.0.0 )then
                  angle =      asin( (yh(j)-yctest)/rr )
                else
                  angle = pi - asin( (yh(j)-yctest)/rr )
                endif
                do k=1,kmax
                  v=vtmp(i,j,k)*cos(angle)-utmp(i,j,k)*sin(angle)
                  davg(ir,k)=davg(ir,k)+v
                enddo
              endif
            endif
          enddo
          enddo
!-----------------------------------------------------------------------
!  communicate:
#ifdef MPI
        if( myid.eq.0 )then
          call MPI_REDUCE(mpi_in_place,davg,nr*kmax,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
          call MPI_REDUCE(mpi_in_place,navg,nr     ,mpi_integer         ,mpi_sum,0,mpi_comm_world,ierr)
        else
          call MPI_REDUCE(  davg      ,davg,nr*kmax,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
          call MPI_REDUCE(  navg      ,navg,nr     ,mpi_integer         ,mpi_sum,0,mpi_comm_world,ierr)
        endif
#endif
!-----------------------------------------------------------------------
        IF(myid.eq.0)THEN
          ! Search for maximum value for this center point:
          vmax = 0.0
          do i=1,nr
            ravg(i) = 1.0/dble(max(1,navg(i)))
          enddo
          do k=1,kmax
          do i=1,nr
            davg(i,k) = davg(i,k)*ravg(i)
            if( davg(i,k).gt.vmax )then
              vmax = davg(i,k)
              rmax = (i-1)*ddr + 0.5*ddr
              zmax = zh(1,1,k)
            endif
          enddo
          enddo
          if( vmax.gt.mostmax )then
            mostmax = vmax
            mostx = xctest
            mosty = yctest
            mostr = rmax
            mostz = zmax
            ictest = ii
            jctest = jj
          endif
        ENDIF
        ENDIF  checkit
      enddo  iloop
      enddo  jloop

      IF(myid.eq.0)THEN
        xcenter = mostx
        ycenter = mosty
        icenter = ictest
        jcenter = jctest
      ENDIF

    if(myid.eq.0) print *,'  New values: '
    if(myid.eq.0) print *,'    xcenter,ycenter = ',xcenter,ycenter
    if(myid.eq.0) print *,'    icenter,jcenter = ',icenter,jcenter
    if(myid.eq.0) print *,'    vmax,rmax,zmax  = ',mostmax,mostr,mostz


  ELSE  docenter

    ! just use center of domain:
    icenter = nx/2 + 1
    jcenter = ny/2 + 1
    xcenter = minx + 0.5*(maxx-minx)
    ycenter = miny + 0.5*(maxy-miny)

  ENDIF  docenter

!-----------------------------------------------------------------------

    IF(myid.eq.0)THEN
      print *,'  writing azimavg stats file'

      irec1 = 1 + nvars*(nwritea-1)

        newstring = string

        newstring(totlen+1:totlen+18) = '_azimavg_stats.dat'
        print *,newstring
        open(unit=67,file=newstring,form='unformatted',access='direct',recl=4,status='unknown')

        write(67,rec=irec1) float(icenter)
        irec1=irec1+1
        write(67,rec=irec1) float(jcenter)
        irec1=irec1+1
        write(67,rec=irec1) xcenter
        irec1=irec1+1
        write(67,rec=irec1) ycenter
        irec1=irec1+1
        write(67,rec=irec1) mostmax
        irec1=irec1+1
        write(67,rec=irec1) mostr
        irec1=irec1+1
        write(67,rec=irec1) mostz
        irec1=irec1+1

        close(unit=67)

      print *,'  done'
    ENDIF

!-----------------------------------------------------------------------
!  more comms:

#ifdef MPI
    call MPI_BCAST(icenter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(jcenter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(xcenter,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ycenter,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

!-----------------------------------------------------------------------
!  GrADS descriptor file:

    IF( myid.eq.0 )THEN
      newstring = string
      newstring(totlen+1:totlen+18) = '_azimavg_stats.ctl'
      open(unit=50,file=newstring,status='unknown')
      rstring = '                                                                      '
      rstring(1:baselen) = output_basename(1:baselen)
      rstring(baselen+1:baselen+1+18) = '_azimavg_stats.dat'
      write(50,231) rstring(1:70)
231   format('dset ^',a70)
      write(50,202)
      write(50,203)
      write(50,204)
      write(50,205)
      write(50,206)
      write(50,207) nwritea
      write(50,208) nvars
      write(50,209) 'icenter ','i value of center point                           '
      write(50,209) 'jcenter ','j value of center point                           '
      write(50,209) 'xcenter ','x (m) value of center point                       '
      write(50,209) 'ycenter ','y (m) value of center point                       '
      write(50,209) 'vmax    ','max azimuthally averaged tangential velocity (m/s)'
      write(50,209) 'rmax    ','radius of Vmax (m)                                '
      write(50,209) 'zmax    ','height (m ASL) of Vmax                            '
      write(50,210)
      close(unit=50)
202   format('title CM1 azimuthally averaged output, stats and info')
203   format('undef -99999999.')
204   format('xdef 1 linear 0.0 1.0')
205   format('ydef 1 linear 0.0 1.0')
206   format('zdef 1 linear 0.0 1.0')
207   format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
208   format('vars ',i4)
209   format(a8,' 0 99 ',a50)
210   format('endvars')
    ENDIF

!-----------------------------------------------------------------------

    deallocate( navg )
    deallocate( ravg )
    deallocate( davg )

    if(myid.eq.0) print *,'  .... leaving getcenter '

    return
    end subroutine getcenter


  END MODULE azimavg_module
