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
!  DEPENDENCIES: boundary, deriv, dim, dust, energies, eos, growth, io,
!    kernel, mpiutils, options, part, physcon, step_lf_global, testutils,
!    timestep, unifdis, units
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
 !
 ! check stokes number interpolation
 !
 call check_stokes_number(ntests,npass)
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
                          dustfrac,dustevol,ddustfrac,temperature,iphase,iamdust,maxtypes,St
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use unifdis,        only:set_unifdis
 use eos,            only:ieos,polyk,gamma
 use options,        only:alpha,alphamax
 use physcon,        only:au,solarm,Ro
 use dim,            only:periodic
 use timestep,       only:dtmax
 use io,             only:iverbose
 use units,          only:set_units
 use mpiutils,       only:reduceall_mpi
 use growth,         only:ifrag,get_vrelonvfrag,vfrag,isnow,vfragin,vfragout,rsnow,iinterpol
 integer, intent(inout) :: ntests,npass
 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: nx, itype, npart_previous, i, j, nsteps, ncheck(7), nerr(7)
 real :: deltax, dz, hfact, totmass, rhozero, errmax(7), dtext_dum
 real :: t, dt, dtext, dtnew
 real :: csj = 1., Stj = 1., s, Vt,vrelonvfrag = 10.
 real :: slast = 0., si = 0., so = 0., r = 0.
 real :: T08,tau
 real :: sinit = 1., dens = 1.
 integer :: switch
 real, parameter :: tols = 2.e-5
 logical :: do_output = .false.

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
    St(i)         = Stj
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
 iinterpol = .false.
 !
 !
 !
 dt = 1.e-3
 nsteps = 100
 t = 0
 dtmax = nsteps*dt
 !
 ! run growingbox problem
 !
 ncheck(:) = 0
 nerr(:) = 0
 errmax(:) = 0.
 !
 ! usefull variables used along the tests
 !
 Vt = sqrt(2**(0.5)*alpha*Ro)*csj
 vfrag = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vfrag < vrel : fragmentation
 tau = 1/(sqrt(2**1.5*Ro*alpha))
 switch = abs(nsteps/2)
 slast = sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(switch*dt)
 !
 ! testing pure growth with St=cst & St=f(size)
 !
 write(*,"(/,a)")'------------------ pure growth (ifrag = 0, St = const) ------------------'
 !
 ! call deriv the first time around
 !
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
            Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)
 do i=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    s = sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    do j=1,npart
    !print*,dustprop(1,142),s
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(1),ncheck(1),errmax(1))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(1),nerr(1),errmax(1),tols)
 write(*,"(/,a)")'------------------ ifrag = 0, St=f(size) inspired from Laibe et al. (2008) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 St(:) = Stj

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t = t + dt
    dtext = dt
    T08 = t/tau+2*sqrt(sinit)*(1+sinit/3)
    call step(npart,npart,t,dt,dtext,dtnew) !--integrate dust size
    s = (8+9*T08**2+3*T08*sqrt(16+9*T08**2))**(1./3.)/2 + 2/(8+9*T08**2+3*T08*sqrt(16+9*T08**2))**(1./3.) - 2
    St(:) = Stj*s
    do j=1,npart
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(2),ncheck(2),errmax(2))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(2),nerr(2),errmax(2),tols)
 !
 ! testing pure fragmentation (ifrag = 1 & ifrag = 2)
 !
 write(*,"(/,a)")'------------------ pure fragmentation (ifrag = 1, St = const) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 St(:)         = Stj
 ifrag         = 1

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew) !--integrate dust size
    s = sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    do j=1,npart
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(3),ncheck(3),errmax(3))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(3),nerr(3),errmax(3),tols)

 write(*,"(/,a)")'------------------ pure fragmentation (ifrag = 2, St = const) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 ifrag         = 2

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew) !--integrate dust size
    s = sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*vrelonvfrag**2/(vrelonvfrag**2+1)*t
    do j=1,npart
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(4),ncheck(4),errmax(4))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(4),nerr(4),errmax(4),tols)
 !
 ! testing pure growth then at half steps switch to pure fragmentation
 !
 write(*,"(/,a)")'------------------ growth-fragmentation switch ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 vfrag = vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel < vfrag : growth
 ifrag = 1

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t = t + dt
    dtext = dt
    if (i==switch) vfrag = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)  ! vrel > vfrag : fragmentation
    if (i<=switch) then
       s = sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    else
       s = slast - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*(t-switch*dt)
    endif
    call step(npart,npart,t,dt,dtext,dtnew)
    if (do_output) call write_file(i,dt,xyzh,dustprop/sinit,npart,'switch_')
    do j=1,npart
       call checkvalbuf(dustprop(1,j),s,tols,'size',nerr(5),ncheck(5),errmax(5))
    enddo
 enddo

 call checkvalbuf_end('size match exact solution',ncheck(5),nerr(5),errmax(5),tols)
 !
 ! testing growth inside the snow line and fragmentation outside of it
 !
 write(*,"(/,a)")'------------------ position based snow line (ifrag = (in:0, out:1)) ------------------'
 !
 ! initialise again
 !
 dustprop(1,:) = sinit
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 vfragin = vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel < vfrag : growth
 vfragout = 1/vrelonvfrag*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1) ! vrel > vfrag : fragmentation
 isnow = 1
 rsnow = 0.2

 t = 0

 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)

 call init_step(npart,t,dtmax)

 do i=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    si = sinit + rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    so = sinit - rhozero/dens*sqrt(2.)*Vt*sqrt(Stj)/(Stj+1)*t
    if (do_output) call write_file(i,dt,xyzh,dustprop/sinit,npart,'snowline_')
    do j=1,npart
       r = sqrt(xyzh(1,j)**2+xyzh(2,j)**2)
       if (r < rsnow) call checkvalbuf(dustprop(1,j),si,tols,'size',nerr(6),ncheck(6),errmax(6))
       if (r > rsnow) call checkvalbuf(dustprop(1,j),so,tols,'size',nerr(7),ncheck(7),errmax(7))
    enddo
 enddo
 call checkvalbuf_end('size match exact solution (in)',ncheck(6),nerr(6),errmax(6),tols)
 call checkvalbuf_end('size match exact solution (out)',ncheck(7),nerr(7),errmax(7 ),tols)

 if (all(nerr(1:7)==0)) npass = npass + 1
 ntests = ntests + 1

end subroutine test_growingbox

subroutine check_stokes_number(ntests,npass) 
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:igas,idust,npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,dustevol,ddustfrac,temperature,iphase,iamdust,maxtypes,St,xyzmh_ptmass
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use energies,       only:compute_energies,ekin
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use eos,            only:ieos,polyk,gamma
 use dust,           only:K_code,idrag
 use growth,         only:ifrag,iinterpol
 use options,        only:alpha,alphamax,use_dustfrac
 use unifdis,        only:set_unifdis
 use dim,            only:periodic,mhd,use_dust,ndusttypes
 use timestep,       only:dtmax
 use io,             only:iverbose
 use mpiutils,       only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: nx, itype, npart_previous, i, j, nsteps, ncheck(1), nerr(1)
 real :: deltax, dz, hfact, totmass, rhozero, errmax(1), dtext_dum
 real :: Stcomp, r, sinit = 1., dens = 1.,s
 real :: t, dt, dtext, dtnew
 real, parameter :: tolst = 2.e-4
 
 write(*,"(/,a)")'--> testing STOKES NUMBER INTERPOLATION'
 write(*,"(/,a)")'------------------ ts = const ------------------'

 !
 ! initialise
 !
 dustprop(1,:) = sinit
 dustprop(2,:) = dens
 dustprop(3,:) = 0.
 dustprop(4,:) = 0.
 xyzmh_ptmass(4,1) = 1.
 !
 ! setup for dustybox problem
 !
 nx = 32
 deltax = 1./nx
 dz = 2.*sqrt(6.)/nx
 call set_boundary(-0.5,0.5,-0.25,0.25,-dz,dz)
 hfact = hfact_default
 rhozero = 1.
 totmass = rhozero*dxbound*dybound*dzbound
 npart = 0

 do itype=1,2
    npart_previous = npart
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                     deltax,hfact,npart,xyzh,verbose=.false.)
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
       vxyzu(:,i) = 0.
       fext(:,i) = 0.
       if (mhd) Bevol(:,i) = 0.
       if (use_dust) then
          dustevol(:,i) = 0.
          dustfrac(:,i) = 0.
       endif
    enddo
    npartoftype(itype) = npart - npart_previous
    npartoftypetot(itype) = reduceall_mpi('+',npartoftype(itype))
    massoftype(itype) = totmass/npartoftypetot(itype)
 enddo
 !
 ! runtime parameters
 !
 K_code = 1.
 ieos = 1
 idrag = 2
 ifrag = 0
 polyk = 1.
 gamma = 1.
 alpha = 0.
 alphamax = 0.
 iverbose = 0
 iinterpol = .true.
 !
 ! call derivs the first time around
 !
 dt = 1.e-3
 nsteps = 100
 t = 0
 dtmax = nsteps*dt
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,t,0.,dtext_dum)
 !
 ! run dustybox problem
 !
 ncheck(:) = 0
 nerr(:) = 0
 errmax(:) = 0.
 call init_step(npart,t,dtmax)
 do i=1,nsteps
    t = t + dt
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)

    do j=1,npart
       if (iamdust(iphase(j))) then
          r      = sqrt(xyzh(1,j)**2+xyzh(2,j)**2)
          Stcomp = 1/(2*K_code*r**(1.5))
          call checkvalbuf(St(j),Stcomp,tolst,'St',nerr(1),ncheck(1),errmax(1))
       endif
    enddo
 enddo

 call checkvalbuf_end('Stokes number interpolation match exact solution',ncheck(1),nerr(1),errmax(1),tolst)

 ntests = ntests + 1
 if (all(nerr(1:1)==0)) npass = npass + 1

 end subroutine check_stokes_number
!---------------------------------------------------
!+
!  write an output file with x, y, z ,
!  dustprop(1) (size) and dustprop(4) (vrel/vfrag)
!+
!---------------------------------------------------
subroutine write_file(step,dt,xyzh,dustprop,npart,prefix)
 real, intent(in)              :: dt
 real, intent(in)              :: xyzh(:,:),dustprop(:,:)
 character(len=*), intent(in)  :: prefix
 integer, intent(in)           :: npart,step
 character(len=30)             :: filename,str
 integer                       :: i,lu
 real                          :: r2

 write(str,"(i000.4)") step
 filename = prefix//trim(adjustl(str))//'.txt'
 open(newunit=lu,file=filename,status='replace')
 write(lu,*) step*dt
 do i=1,npart
    write(lu,*) xyzh(1,i),xyzh(2,i),dustprop(1,i),dustprop(4,i)
 enddo
 close(lu)

end subroutine write_file
#endif
#endif

end module testgrowth
