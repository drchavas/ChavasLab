  MODULE comm_module

  implicit none

  CONTAINS

!-----------------------------------------------------------------------
!  message passing routines
!-----------------------------------------------------------------------


      integer function nabor(i,j,nx,ny)
      implicit none
      integer i,j,nx,ny
      integer newi,newj

      newi=i
      newj=j

      if ( newi .lt.  1 ) newi = nx
      if ( newi .gt.  nx) newi = 1

      if ( newj .lt.  1 ) newj = ny
      if ( newj .gt.  ny) newj = 1

      nabor = (newi-1) + (newj-1)*nx

      end function nabor

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prepcorners(s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,comm)
      use input
      use bc_module
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      integer, intent(inout), dimension(rmp) :: reqs_p
      integer, intent(in) :: comm

      integer :: i,j

!--------------------------------------------!
!  This subroutine is ONLY for parcel_interp !
!--------------------------------------------!

      IF( comm.eq.1 )THEN
        call bcs(s)
      ENDIF

      IF( bbc.eq.1 .or. bbc.eq.2 .or. bbc.eq.3 )THEN
        ! extrapolate:
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          s(i,j,0) = cgs1*s(i,j,1)+cgs2*s(i,j,2)+cgs3*s(i,j,3)
        enddo
        enddo
      ENDIF

      IF( tbc.eq.1 .or. tbc.eq.2 )THEN
        ! extrapolate:
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          s(i,j,nk+1) = cgt1*s(i,j,nk)+cgt2*s(i,j,nk-1)+cgt3*s(i,j,nk-2)
        enddo
        enddo
      ENDIF

      end subroutine prepcorners


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prepcornert(t,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,comm)
      use input
      use bc_module
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: t
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      integer, intent(inout), dimension(rmp) :: reqs_p
      integer, intent(in) :: comm

      integer :: i,j
      real :: c1,c2

!--------------------------------------------!
!  This subroutine is ONLY for parcel_interp !
!--------------------------------------------!

      IF( comm.eq.1 )THEN
        call bcw(t,0)
      ENDIF

      end subroutine prepcornert



  END MODULE comm_module
