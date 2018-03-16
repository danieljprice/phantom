!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: ekin
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: ekin [no arguments]
!
!  DEPENDENCIES: dustywaves
!+
!--------------------------------------------------------------------------
program ekin
 use dustywaves
 implicit none
 integer, parameter :: nx = 200
 integer, parameter :: npts = 200
 integer, parameter :: nk = 7
 real :: xplot(nx),vgas(nx),vdust(nx),rhogas(nx),rhodust(nx)
 real :: dt,tmax,tmin,time,ampl,cs,Kdragin,lambda,x0,ymean_gas,ymean_dust
 real :: ekintot,ekinm,vxi,rhoi,rhodi,vdxi,ekindust,ekingas,ekingasm,ekindustm
 real :: xmin,xmax,dx,gam,eint,etot,totmass,dk,rhomax
 integer :: ierr,i,j,k
 character(len=120) :: filename,string

 Kdragin     = 1.
 ampl        = 1.e-4
 ymean_gas   = 1.
 ymean_dust  = 1.
 cs          = 1.
 x0          = 0.
 lambda      = 1.
 tmin        = 0.
 tmax        = 10.
 gam         = 5./3.
 totmass     = 1.

 print*,'Enter dust-to-gas ratio'
 read*,ymean_dust
 write(string,"(f10.3)") ymean_dust
 ymean_dust = ymean_dust*ymean_gas

 xmin = 0.
 xmax = 1.
 dx = (xmax - xmin)/real(nx)
 do i=1,nx
    xplot(i) = xmin + (i-1)*dx
 enddo

 dt = (tmax - tmin)/real(npts-1)

 dk = 6./real(nk-1)
 overk: do k=1,nk
    Kdragin= 10**(-3. + (k-1)*dk)

    write(filename,"(f10.3)") Kdragin
    filename = 'dustywave-dgr'//trim(adjustl(string))//'-K'//trim(adjustl(filename))//'.out'

    open(unit=1,file=trim(filename),status='replace',form='formatted')
    do i=1,npts
       time = tmin + (i-1)*dt
       call exact_dustywave(time,ampl,cs,Kdragin,lambda,x0,ymean_gas,ymean_dust,xplot,vgas,vdust,rhogas,rhodust,ierr)
       ekingas = 0.
       ekingasm = 0.
       ekindustm = 0.
       ekindust = 0.
       eint = 0.
       rhomax = -1.
       do j=1,nx
          vxi = vgas(j)
          vdxi = vdust(j)
          rhoi = rhogas(j)
          rhodi = rhodust(j)
          ekingas = ekingas + rhoi*vxi*vxi*dx
          ekindust = ekindust + rhodi*vdxi*vdxi*dx
          !pri = cs**2*rhoi/gam
          !ui  = pri/((gam-1.)*rhoi)
          !eint = eint + rhoi*ui
          rhomax = max(rhoi,rhomax)
       enddo
       ekintot = 0.5*(ekingas + ekindust)
!    eint = eint/real(npts)
       eint = 1.5
       etot = ekintot + eint

       if (abs(ekintot) < 1.e-3 .and. ekintot /= 0.) then
          write(1,*) time,ekintot,eint,0.,0.,etot,0.,0.,rhomax
       endif
    enddo
    close(1)

 enddo overk

end program ekin
