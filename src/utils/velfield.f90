!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: velfield
!
!  DESCRIPTION:
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module velfield
 public :: set_velfield

 private

contains

subroutine set_velfield(xyzh,vxyzu,npart)
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: vxyzu(:,:)
 real(kind=8), parameter :: pi=3.1415927d0
 integer, parameter :: ngrid = 64

 integer :: i,j,k,iseed,iseed2,ii,jj,kk
 real(kind=8) :: nindex,kmod,kdotq,lf
 real(kind=8) :: amp_phi(6,2*ngrid,2*ngrid,2*ngrid)
 real(kind=8) :: pow(2*ngrid,2*ngrid,2*ngrid),vgridxyz(3,2*ngrid,2*ngrid,2*ngrid)
 real(kind=8) :: xk,yk,zk,xi,yi,zi
 real(kind=8) :: rtot,aa,ampxi,ampyi,ampzi,sinphix,sinphiy,sinphiz
 real(kind=8) :: velx,vely,velz,contrib,powsum,sigma2,phix,phiy,phiz
 real :: tprev,t2

 rtot=1.0d0
 lf=2.0d0*rtot
 aa=0.5d-1/32.52d0
 nindex=-6.0d0
!
!--standard case:
!
 iseed=17
 iseed2=11
 !print*,'enter 2 seeds (integers, e.g. 17,11)'
 !read (*,*) iseed
 !read (*,*) iseed2
 !print*,'enter index (e.g. n=-17/3 for kolmogorov)'
 !read (*,*) nindex
!
!--gives p(k)=k^{nindex) power spectrum
!     note: rayldev returns the square root of the -log of a random number
!
 powsum=0.d0
 do k=1,ngrid
    do j=1,2*ngrid
       do i=1,2*ngrid
          amp_phi(4,i,j,k)=(-pi + 2.d0*pi*ran3(iseed))
          amp_phi(5,i,j,k)=(-pi + 2.d0*pi*ran3(iseed))
          amp_phi(6,i,j,k)=(-pi + 2.d0*pi*ran3(iseed))
          kmod=sqrt(real((i-ngrid)**2+(j-ngrid)**2+k*k))
          pow(i,j,k)=aa*(kmod)**(nindex)
          powsum=powsum+pow(i,j,k)
          amp_phi(1,i,j,k)=rayldev(iseed2)*sqrt(pow(i,j,k))
          amp_phi(2,i,j,k)=rayldev(iseed2)*sqrt(pow(i,j,k))
          amp_phi(3,i,j,k)=rayldev(iseed2)*sqrt(pow(i,j,k))
       enddo
    enddo
 enddo

 sigma2=0.d0
 tprev = 0.

!$omp parallel do default(none) &
!$omp shared(lf,amp_phi,npart,ngrid) &
!$omp shared(xyzh,vxyzu) &
!$omp private(ipart,kk,jj,ii,i,j,k,kx,ky,kz,kdotq,contrib,velx,vely,velz) &
 overz: do kk = 1, 2*ngrid
    !print*,'kk ',kk
    call cpu_time(t2)
    if (kk > 1) print*,kk,' time = ',t2-tprev,' to finish = ',(2*ngrid-kk)*(t2-tprev)
    tprev = t2

    zi = (kk-ngrid-0.5)/real(ngrid)
    overy: do jj = 1, 2*ngrid
       !print*,'jj ',jj
       yi = (jj-ngrid-0.5)/real(ngrid)
       overx: do ii = 1, 2*ngrid
          xi = (ii-ngrid-0.5)/real(ngrid)

          !do ipart=1,npart
          !   xi = xyzh(1,ipart)
          !   yi = xyzh(2,ipart)
          !   zi = xyzh(3,ipart)
          !if (mod(ipart,1000)==0) then
          !   call cpu_time(t2)
          !   print*,ipart,' time = ',t2-tprev,' per part = ',(t2-tprev)/1000.,' to finish = ',(npart-ipart)/1000.*(t2-tprev)
          !   tprev = t2
          !endif

          velx = 0.
          vely = 0.
          velz = 0.

          do k=1,ngrid
             do j=1,2*ngrid
                do i=1,2*ngrid
                   xk = real(i-ngrid)
                   yk = real(j-ngrid)
                   zk = real(k)
!!                        kmod=(2.d0*pi/lf)*
!!     &                       sqrt(real(kx**2+ky**2+kz**2))
                   kdotq=(2.d0*pi/lf)*(xk*xi + yk*yi + zk*zi)
                   ampxi = amp_phi(1,i,j,k)
                   ampyi = amp_phi(2,i,j,k)
                   ampzi = amp_phi(3,i,j,k)
                   phix = amp_phi(4,i,j,k)
                   phiy = amp_phi(5,i,j,k)
                   phiz = amp_phi(6,i,j,k)
                   sinphix = sin(kdotq+phix)
                   sinphiy = sin(kdotq+phiy)
                   sinphiz = sin(kdotq+phiz)
                   contrib = ampzi*(2.d0*pi/lf)*yk*2.d0*sinphiz &
                       - ampyi*(2.d0*pi/lf)*zk*2.d0*sinphiy
                   velx = velx + contrib
                   contrib = ampxi*(2.d0*pi/lf)*zk*2.d0*sinphix &
                   &                      - ampzi*(2.d0*pi/lf)*xk*2.d0*sinphiz
                   vely = vely + contrib
                   contrib = ampyi*(2.d0*pi/lf)*xk*2.d0*sinphiy &
                       - ampxi*(2.d0*pi/lf)*yk*2.d0*sinphix
                   velz = velz + contrib
                enddo
             enddo
          enddo
          vgridxyz(1,ii,jj,kk) = velx
          vgridxyz(2,ii,jj,kk) = vely
          vgridxyz(3,ii,jj,kk) = velz

!      vxyzu(1,ipart) = velx
!      vxyzu(2,ipart) = vely
!      vxyzu(3,ipart) = velz
       enddo overx
    enddo overy
 enddo overz
!$omp end parallel do

end subroutine set_velfield

!*** the hopefully better random generator ****
!** note that this function is double precision!!

!
!
!     the function subroutine ran3.
!
!      this suroutine generates a uniform random deviate on (0,1).
!
double precision function ran3(idum)
!
 implicit double precision (m)
 double precision fac
 parameter (mbig=4000000., mseed=1618033., mz=0., fac=1./mbig)
 dimension ma(55)
 save
 data iff/0/
 if (idum < 0.or.iff==0) then
    iff=1
    mj=mseed-iabs(idum)
    mj=mod(mj,mbig)
    ma(55)=mj
    mk=1
    do 100 i=1,54
       ii=mod(21*i,55)
       ma(ii)=mk
       mk=mj-mk
       if (mk < mz)mk=mk+mbig
       mj=ma(ii)
100 continue
    do 200 k=1,4
       do 300 i=1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if (ma(i) < mz)ma(i)=ma(i)+mbig
300    continue
200 continue
    inext=0
    inextp=31
    idum=1
 endif
 inext=inext+1
 if (inext==56)inext=1
 inextp=inextp+1
 if (inextp==56)inextp=1
 mj=ma(inext)-ma(inextp)
 if (mj < mz)mj=mj+mbig
 ma(inext)=mj
 ran3=mj*fac
 return
end function ran3

real*8 function rayldev(idum)
 integer, intent(in) :: idum
 double precision random

 random = ran4(idum)
 do while (random==0.0)
    random = ran4(idum)
    write (*,*) 'done ',random
 enddo
 rayldev = dsqrt(-log(random))

 return
end function rayldev

!*** the hopefully better random generator ****
!** note that this function is double precision!!

!
!
!     the function subroutine ran4.
!
!      this suroutine generates a uniform random deviate on (0,1).
!
double precision function ran4(idum)
 implicit double precision (m)
 double precision fac
 parameter (mbig=4000000., mseed=1618033., mz=0., fac=1./mbig)
 dimension ma(55)
 save
 data iff/0/
 if (idum < 0.or.iff==0) then
    iff=1
    mj=mseed-iabs(idum)
    mj=mod(mj,mbig)
    ma(55)=mj
    mk=1
    do i=1,54
       ii=mod(21*i,55)
       ma(ii)=mk
       mk=mj-mk
       if (mk < mz)mk=mk+mbig
       mj=ma(ii)
    enddo
    do k=1,4
       do i=1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if (ma(i) < mz)ma(i)=ma(i)+mbig
       enddo
    enddo
    inext=0
    inextp=31
    idum=1
 endif
 inext=inext+1
 if (inext==56)inext=1
 inextp=inextp+1
 if (inextp==56)inextp=1
 mj=ma(inext)-ma(inextp)
 if (mj < mz)mj=mj+mbig
 ma(inext)=mj
 ran4=mj*fac
 return
end function ran4

end module velfield
