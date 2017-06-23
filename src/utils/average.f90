! $Id$

!***********************************************************************
REAL FUNCTION average(mx,my,mz,f)
!
! Averaging routine, safe from round off errors
!
!-----------------------------------------------------------------------
 implicit none
 integer mx,my,mz,ix,iy,iz
 real :: f(mx,my,mz)
 real(kind=8) sumx,sumy,sumz

 sumz=0.                                                               ! sum over volume
 !$omp parallel do private(ix,iy,iz,sumx,sumy), reduction(+:sumz)
 do iz=1,mz
    sumy=0.                                                             ! sum over planes
    do iy=1,my
       sumx=0.                                                           ! sum along lines
       do ix=1,mx
          sumx=sumx+f(ix,iy,iz)
       enddo
       sumy=sumy+sumx
    enddo
    sumz=sumz+sumy
 enddo
 average=sumz/(real(mx)*real(my)*real(mz))                             ! safe above 2G points
END function average

!***********************************************************************
REAL FUNCTION rms(mx,my,mz,f,aver)
!
! Averaging routine, safe from round off errors
!
!-----------------------------------------------------------------------
 implicit none
 integer mx,my,mz,ix,iy,iz
 real :: f(mx,my,mz)
 real(kind=8) sumx,sumy,sumz
 real aver,average

 aver=average(mx,my,mz,f)
 sumz=0.                                                               ! sum over volume
 !$omp parallel do private(ix,iy,iz,sumx,sumy), reduction(+:sumz)
 do iz=1,mz
    sumy=0.                                                             ! sum over planes
    do iy=1,my
       sumx=0.                                                           ! sum along lines
       do ix=1,mx
          sumx=sumx+(f(ix,iy,iz)-aver)**2
       enddo
       sumy=sumy+sumx
    enddo
    sumz=sumz+sumy
 enddo
 rms=sqrt(sumz/(real(mx)*real(my)*real(mz)))                           ! safe above 2G points
END function rms
