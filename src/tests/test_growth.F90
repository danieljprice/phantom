!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testgrowth
!
!  DESCRIPTION:
!   Unit tests of the growth module
!
!  REFERENCES:
!
!  OWNER: Arnaud Vericel
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, deriv, dim, eos, growth, io, kernel, mpiutils,
!    options, part, physcon, step_lf_global, testutils, timestep, unifdis,
!    units
!+
!--------------------------------------------------------------------------
module testgrowth
 use testutils, only:checkval
 use io,        only:id,master
 implicit none
 public :: test_growth

 private

contains

subroutine test_growth(ntests,npass)
#ifdef DUST
#ifdef DUSTGROWTH
 use growth,  only:init_growth,get_growth_rate,ifrag,isnow
 use physcon,     only:solarm,au
 use units,       only:set_units
 use mpiutils,    only:barrier_mpi
#endif
#endif
 integer, intent(inout) :: ntests,npass

#ifdef DUST
#ifdef DUSTGROWTH
 integer :: nfailed(5),ierr,iregime !don't forget the dimension of nfailed

 if (id==master) write(*,"(/,a)") '--> TESTING DUSTGROWTH MODULE'

 call set_units(mass=solarm,dist=au,G=1.d0)

 if (id==master) write(*,"(/,a)") '--> testing growth initialisation'

 nfailed = 0
 ntests = ntests + 1
 do ifrag=0,2
    do isnow=0,2
       call init_growth(ierr)
       call checkval(ierr,0,0,nfailed(ifrag+isnow+1),'growth initialisation')
    enddo
 enddo
 if (all(nfailed==0)) npass = npass + 1

 !
 ! GROWINGBOX test
 !
 call test_growingbox(ntests,npass)
 call barrier_mpi()

 if (id==master) write(*,"(/,a)") '<-- DUSTGROWTH TEST COMPLETE'
#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING DUSTGROWTH TEST (REQUIRES -DDUST -DDUSTGROWTH)'
#endif
#endif

end subroutine test_growth

#ifdef DUST
#ifdef DUSTGROWTH
!----------------------------------------------------
!+
!  Growingbox test
!  This tests the dust growth algorithm
!+
!----------------------------------------------------
subroutine test_growingbox(ntests,npass)
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:idust,npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,rhoh,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,dustevol,ddustfrac,temperature,iphase,iamdust,maxtypes
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use unifdis,        only:set_unifdis
 use eos,            only:ieos,polyk,gamma
 use options,        only:alpha,alphamax
 use physcon,        only:au,solarm,Ro
 use dim,            only:periodic,mhd
 use timestep,       only:dtmax
 use io,             only:iverbose
 use units,          only:set_units,udist,unit_density
 use mpiutils,       only:reduceall_mpi
 use growth,         only:ifrag,get_vrelonvfrag,vfrag,isnow,vfragin,vfragout,rsnow
 integer, intent(inout) :: ntests,npass
 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: nx, itype, npart_previous, i, j, nsteps, ncheck(6), nerr(6)
 real :: deltax, dz, hfact, totmass, rhozero, errmax(6), dtext_dum
 real :: t, dt, dtext, dtnew
 real :: csj, Stj, s, Vt,vrelonvfrag = 10.
 real :: slast = 0., si = 0., so = 0., r = 0.
 real :: sinit = 1.e-2, dens = 1.
 integer :: switch
 real, parameter :: tols = 2.e-5

 if (periodic) then
    if (id==master) write(*,"(/,a)") '--> testing GROWINGBOX'
 else
    if (id==master) write(*,"(/,a)") '--> skipping GROWINGBOX (need -DPERIODIC)'
    return
 endif
 !
 ! setup for growingbox problem
 !
 nx = 32
 deltax = 1./nx
 dz = 2.*sqrt(6.)/nx
 call set_boundary(-0.5,0.5,-0.25,0.25,-dz,dz)
 hfact = hfact_default
 rhozero = 1.
 totmass = rhozero*dxbound*dybound*dzbound
 npart = 0

 npart_previous = npart
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,verbose=.false.)
 call set_units(mass=solarm,dist=au,G=1.d0)
 do i=npart_previous+1,npart
    call set_particle_type(i,idust)
    vxyzu(:,i) = 0.
    fxyzu(:,i) = 0.
    fext(:,i) = 0.
    dustprop(1,i) = sinit
    dustprop(2,i) = dens
    dustprop(3,i) = 0.
    dustprop(4,i) = 0.
    dustprop(5,i) = 0.
 enddo
 npartoftype(idust) = npart - npart_previous
 npartoftypetot(idust) = reduceall_mpi('+',npartoftype(idust))
 massoftype(idust) = totmass/npartoftypetot(idust)
 !
 ! runtime parameters
 !
 ifrag = 0
 isnow = 0
 alpha = 1.
 iverbose = 0
 ieos = 1
 polyk = 1.
 gamma = 1.
 !
 ! call deriv the first time around
 !
 dt = 1.e-3
 nsteps = 100
 t = 0
 dtmax = nsteps*dt
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)
 !
 ! run growingbox problem
 !
 csj = 1.
 Stj = 1.
 ncheck(:) = 0
 nerr(:) = 0
 errmax(:) = 0.
 Vt = sqrt(2**(0.5)*alpha*Ro)*csj

 !
 ! ifrag = 0 : pure growth
 !
 write(*,"(/,a)")'------------------ pure growth (ifrag = 0) ------------------'

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    do j=1,npart
       call get_vrelonvfrag(xyzh(:,j),dustprop(:,j),csj,Stj,0.) !--get vrel
    enddo
    t = t + dt
    dtext = dt

    call step(npart,npart,t,dt,dtext,dtnew)

    do j=1,npart
       s = sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(t-dt) !- scor
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(1),ncheck(1),errmax(1))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(1),nerr(1),errmax(1),tols)

 ntests = ntests + 1

 write(*,"(/,a)")'------------------ pure fragmentation (ifrag = 1) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = 100*sinit
 dustprop(2,:) = dens
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 dustprop(5,:) = 0.
 vfrag = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! 0.01*vrel to assure fragmentation
 ifrag = 1

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    do j=1,npart
       call get_vrelonvfrag(xyzh(:,j),dustprop(:,j),csj,Stj,0.) !--get vrel, vrel/vfrag
    enddo
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew) !--integrate dust size

    do j=1,npart
       s = 100*sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(t-dt)
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(2),ncheck(2),errmax(2))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(2),nerr(2),errmax(2),tols)

 write(*,"(/,a)")'------------------ pure fragmentation (ifrag = 2) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = 100*sinit
 dustprop(2,:) = dens
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 dustprop(5,:) = 0.
 ifrag = 2

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    do j=1,npart
       call get_vrelonvfrag(xyzh(:,j),dustprop(:,j),csj,Stj,0.) !--get vrel, vrel/vfrag
    enddo
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew) !--integrate dust size

    do j=1,npart
       s = 100*sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*vrelonvfrag**2/(vrelonvfrag**2+1)*(t-dt)
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(3),ncheck(3),errmax(3))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(3),nerr(3),errmax(3),tols)

 write(*,"(/,a)")'------------------ growth-fragmentation switch ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = 100*sinit
 dustprop(2,:) = dens
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 dustprop(5,:) = 0.
 vfrag = vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel < vfrag : growth
 ifrag = 1

 t = 0
 switch = abs(nsteps/2)
 slast = 100*sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(switch*dt)

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    do j=1,npart
       call get_vrelonvfrag(xyzh(:,j),dustprop(:,j),csj,Stj,0.) !--get vrel, vrel/vfrag
    enddo
    t = t + dt
    dtext = dt
    if (i==switch-1) vfrag = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)  ! at nsteps/2-1, vrel > vfrag : fragmentation
    if (i<=switch) then
       s = 100*sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(t-dt)
    else
       s = slast - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(t-(switch-1)*dt)
    endif
    call step(npart,npart,t,dt,dtext,dtnew)
    do j=1,npart
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(4),ncheck(4),errmax(4))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(4),nerr(4),errmax(4),tols)

write(*,"(/,a)")'------------------ box with a position based snow line ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = 100*sinit
 dustprop(2,:) = dens
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 dustprop(5,:) = 0.
 vfragin = vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel < vfrag : growth
 vfragout = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel > vfrag : fragmentation
 ifrag = 1
 isnow = 1
 rsnow = 0.5

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    do j=1,npart
       call get_vrelonvfrag(xyzh(:,j),dustprop(:,j),csj,Stj,0.) !--get vrel, vrel/vfrag
    enddo
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    si = 100*sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(t-dt)
    so = 100*sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(t-dt)
    do j=1,npart
       r = sqrt(xyzh(1,j)**2+xyzh(2,j)**2)
       if (r < rsnow) call checkvalbuf(dustprop(1,j),si,tols,'size',nerr(5),ncheck(5),errmax(5))
       if (r > rsnow) call checkvalbuf(dustprop(1,j),so,tols,'size',nerr(6),ncheck(6),errmax(6))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution (in)',ncheck(5),nerr(5),errmax(5),tols)
 call checkvalbuf_end('size match exact solution (out)',ncheck(6),nerr(6),errmax(6),tols)

 if (all(nerr(1:6)==0)) npass = npass + 1

end subroutine test_growingbox

#endif
#endif

end module testgrowth
