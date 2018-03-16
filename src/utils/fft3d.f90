!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: fft3d
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: omputils
!+
!--------------------------------------------------------------------------
module fft3d
 use omputils, only:limits_omp

contains
!**********************************************************************
subroutine fft3d_k2 (k2,dx,dy,dz,mx,my,mz)
 implicit none
 integer mx, my, mz
 real :: k2(mx,my,mz)
 real dx, dy, dz, kx, ky, kz
 real, parameter :: pi=3.14519265
 integer i, j, k, izs, ize
!----------------------------------------------------------------------
 call limits_omp(1,mz,izs,ize)
 do k=izs,ize
    kz=2.*pi/(dz*mz)*(k/2)
    do j=1,my
       ky=2.*pi/(dy*my)*(j/2)
       do i=1,mx
          kx=2.*pi/(dx*mx)*(i/2)
          k2(i,j,k)=kx**2+ky**2+kz**2
       enddo
    enddo
 enddo
 !$omp barrier
end subroutine

!**********************************************************************
subroutine fft3df (r,rr,mx,my,mz)
 implicit none
 integer mx,my,mz
 real :: r(mx,my,mz), rr(mx,my,mz)
 integer, parameter :: lw=4096, lwx=2*lw+15, lwy=2*lw+15, lwz=2*lw+15
 real c, wx(lwx), wy(lwy), wz(lwz), fx(mx), fy(my), fz(mz,mx)
 integer i, j, k, iys, iye, izs, ize
!----------------------------------------------------------------------
 if (mx > lw .or. my > lw .or. mz > lw) then
    print *,'fft3db: dim too small',mx,my,mz,lw
    stop
 endif
 call srffti (mx,wx)
 call srffti (my,wy)
 call srffti (mz,wz)

 !$omp barrier
 call limits_omp(1,mz,izs,ize)
 do k=izs,ize
    do j=1,my
       fx(:)=r(:,j,k)
       call srfftf (mx,fx,wx)
       c=1./mx
       rr(:,j,k)=c*fx(:)
    enddo
    do i=1,mx
       fy=rr(i,:,k)
       call srfftf (my,fy,wy)
       c=1./my
       rr(i,:,k)=c*fy
    enddo
 enddo
 !$omp barrier

 call limits_omp(1,my,iys,iye)
 do j=iys,ize
    do k=1,mz
       fz(k,:)=rr(:,j,k)
    enddo
    do i=1,mx
       call srfftf (mz,fz(:,i),wz)
    enddo
    c=1./mz
    do k=1,mz
       rr(:,j,k)=c*fz(k,:)
    enddo
 enddo
 !$omp barrier

end subroutine

!**********************************************************************
subroutine fft3db (rr,r,mx,my,mz)
 implicit none
 integer mx,my,mz
 real :: r(mx,my,mz), rr(mx,my,mz)
 integer, parameter :: lw=4096, lwx=2*lw+15, lwy=2*lw+15, lwz=2*lw+15
 real wx(lwx), wy(lwy), wz(lwz), fx(mx), fy(my), fz(mz,mx)
 integer i, j, k, iys, iye, izs, ize
!----------------------------------------------------------------------

 if (mx > lw .or. my > lw .or. mz > lw) then
    print *,'fft3db: dim too small',mx,my,mz,lw
    stop
 endif
 call srffti (mx,wx)
 call srffti (my,wy)
 call srffti (mz,wz)

 !$omp barrier
 call limits_omp(1,mz,izs,ize)
 do k=izs,ize
    do j=1,my
       fx=rr(:,j,k)
       call srfftb (mx,fx,wx)
       r(:,j,k)=fx
    enddo
    do i=1,mx
       fy=r(i,:,k)
       call srfftb (my,fy,wy)
       r(i,:,k)=fy
    enddo
 enddo
 !$omp barrier

 call limits_omp(1,my,iys,iye)
 do j=iys,ize
    do k=1,mz
       fz(k,:)=r(:,j,k)
    enddo
    do i=1,mx
       call srfftb (mz,fz(:,i),wz)
    enddo
    do k=1,mz
       r(:,j,k)=fz(k,:)
    enddo
 enddo
 !$omp barrier

end subroutine
end module
