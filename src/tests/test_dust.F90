!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testdust
!
!  DESCRIPTION:
!   Unit tests of the dust module
!
!  REFERENCES:
!   Laibe & Price (2011),  MNRAS 418, 1491
!   Laibe & Price (2012a), MNRAS 420, 2345
!   Laibe & Price (2012b), MNRAS 420, 2365
!   Price & Laibe (2015),  MNRAS 451, 5332
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, deriv, dim, dust, energies, eos, growth, io,
!    kernel, mpiutils, options, part, physcon, random, set_dust,
!    step_lf_global, table_utils, testutils, timestep, unifdis, units,
!    vectorutils
!+
!--------------------------------------------------------------------------
module testdust
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_dust

#ifdef DUST
#ifdef DUSTGROWTH
 public :: test_dustybox
#endif
#endif

 private

contains

subroutine test_dust(ntests,npass)
#ifdef DUST
 use dust,        only:idrag,init_drag,get_ts
 use set_dust,    only:set_dustbinfrac
 use physcon,     only:solarm,au
 use units,       only:set_units,unit_density,udist
 use eos,         only:gamma
 use dim,         only:use_dust
 use mpiutils,    only:barrier_mpi
 use options,     only:use_dustfrac
 use table_utils, only:logspace
#ifdef DUSTGROWTH
 use growth,      only:init_growth
#endif
#endif
 integer, intent(inout) :: ntests,npass
#ifdef DUST
 integer :: nfailed(3),ierr,iregime
 real    :: rhoi,rhogasi,rhodusti,spsoundi,tsi,grainsizei,graindensi

 if (id==master) write(*,"(/,a)") '--> TESTING DUST MODULE'

 call set_units(mass=solarm,dist=au,G=1.d0)

 if (id==master) write(*,"(/,a)") '--> testing drag initialisation'

 nfailed = 0
 gamma = 5./3.
 do idrag=1,2
    call init_drag(ierr)
    call checkval(ierr,0,0,nfailed(idrag),'drag initialisation')
 enddo
#ifdef DUSTGROWTH
 call init_growth(ierr)
 call checkval(ierr,0,0,nfailed(3),'growth initialisation')
#endif
 call update_test_scores(ntests,nfailed,npass)

 idrag = 1
 rhoi = 1.e-13/unit_density
 spsoundi = 1.
 grainsizei = 1./udist
 graindensi = 1./unit_density
 rhogasi  = 0.5*rhoi
 rhodusti = 0.5*rhoi
 call get_ts(idrag,grainsizei,graindensi,rhogasi,rhodusti,spsoundi,0.,tsi,iregime)
 call checkval(iregime,1,0,nfailed(1),'deltav=0 gives Epstein drag')
 call update_test_scores(ntests,nfailed(1:1),npass)

 !
 ! Test transition between Epstein/Stokes drag
 !
 call test_epsteinstokes(ntests,npass)
 call barrier_mpi()

 !
 ! Test that drag conserves momentum and energy
 !
 use_dustfrac = .false.
 call test_drag(ntests,npass)
 call barrier_mpi()

 !
 ! DUSTYBOX test
 !
 use_dustfrac = .false.
 call test_dustybox(ntests,npass)
 call barrier_mpi()

 !
 ! DUSTYDIFFUSE test
 !
 use_dustfrac = .true.
 call test_dustydiffuse(ntests,npass)
 call barrier_mpi()

 if (id==master) write(*,"(/,a)") '<-- DUST TEST COMPLETE'
#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING DUST TEST (REQUIRES -DDUST)'
#endif

end subroutine test_dust

#ifdef DUST
!----------------------------------------------------
!+
!  Dustybox test from Laibe & Price (2011, 2012a,b)
!  This tests the *two fluid* dust algorithm
!+
!----------------------------------------------------
subroutine test_dustybox(ntests,npass)
 use dim,            only:maxp,maxalpha
 use boundary,       only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use kernel,         only:hfact_default
 use part,           only:igas,idust,npart,xyzh,vxyzu,npartoftype,massoftype,set_particle_type,&
                          fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                          dustfrac,dustevol,ddustevol,temperature,iphase,iamdust,maxtypes,&
                          ndusttypes,alphaind
 use step_lf_global, only:step,init_step
 use deriv,          only:derivs
 use energies,       only:compute_energies,ekin
 use testutils,      only:checkvalbuf,checkvalbuf_end
 use eos,            only:ieos,polyk,gamma
 use dust,           only:K_code,idrag
 use options,        only:alpha,alphamax
 use unifdis,        only:set_unifdis
 use dim,            only:periodic,mhd,use_dust
 use timestep,       only:dtmax
 use io,             only:iverbose
 use mpiutils,       only:reduceall_mpi
 use kernel,         only:kernelname
#ifdef DUSTGROWTH
 use part,           only:dustgasprop,dustprop
 use growth,         only:ifrag
#endif
 integer, intent(inout) :: ntests,npass
 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: nx, itype, npart_previous, i, j, nsteps
 real :: deltax, dz, hfact, totmass, rhozero, dtext_dum
#ifdef DUSTGROWTH
 integer         :: ncheck(6), nerr(6)
 real            :: errmax(6)
 real, parameter :: toldv = 2.e-4
#else
 integer :: ncheck(5), nerr(5)
 real    :: errmax(5)
#endif
 real :: t, dt, dtext, dtnew
 real :: vg, vd, deltav, ekin_exact, fd
 real :: tol,tolvg,tolfg,tolfd

 if (index(kernelname,'quintic') /= 0) then
    tol = 1.e-5; tolvg = 2.5e-5; tolfg = 3.3e-4; tolfd = 3.3e-4
 else
    tol = 1.e-4; tolvg = 1.e-4; tolfg = 3.e-3; tolfd = 3.e-3
 endif

 if (periodic .and. use_dust) then
    if (id==master) write(*,"(/,a)") '--> testing DUSTYBOX'
 else
    if (id==master) write(*,"(/,a)") '--> skipping DUSTYBOX (need -DPERIODIC and -DDUST)'
    return
 endif

#ifdef DUSTGROWTH
 if (id==master) write(*,"(/,a)") '--> Adding dv interpolation test'
#endif

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
 fxyzu = 0.
 dustprop = 0.
 ddustprop = 0.
 ddustevol = 0.
 dBevol = 0.
 if (maxalpha==maxp) alphaind(:,:) = 0.

 itype = igas
 npart_previous = npart
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,verbose=.false.)
 do i=npart_previous+1,npart
    call set_particle_type(i,itype)
    vxyzu(:,i) = 0.
    if (iamdust(iphase(i))) then
       vxyzu(1,i) = 1.
    endif
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


 ndusttypes = 1         ! only works with one dust type currently
 do j=1,ndusttypes
    itype = idust + j - 1
    npart_previous = npart
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                     deltax,hfact,npart,xyzh,verbose=.false.)
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
       vxyzu(:,i) = 0.
       if (iamdust(iphase(i))) then
          vxyzu(1,i) = 1.
       endif
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
 K_code = 0.35
 ieos = 1
 idrag = 2
#ifdef DUSTGROWTH
 ifrag = -1
#endif
 polyk = 1.
 gamma = 1.
 alpha = 0.
 alphamax = 0.
 iverbose = 0
 !
 ! call derivs the first time around
 !
 dt = 1.e-3
 nsteps = 100
 t = 0
 dtmax = nsteps*dt
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum)
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
    !print*,'t = ',t
    call step(npart,npart,t,dt,dtext,dtnew)
    call compute_energies(t)

    deltav = exp(-2.*K_code*t)
    vg = 0.5*(1. - deltav)
    vd = 0.5*(1. + deltav)
    fd = K_code*(vg - vd)
    do j=1,npart
       if (iamdust(iphase(j))) then
          call checkvalbuf(vxyzu(1,j),vd,tol,'vd',nerr(1),ncheck(1),errmax(1))
          call checkvalbuf(fxyzu(1,j),fd,tolfd,'fd',nerr(2),ncheck(2),errmax(2))
#ifdef DUSTGROWTH
          call checkvalbuf(dustgasprop(4,j),deltav,toldv,'dv',nerr(6),ncheck(6),errmax(6))
#endif
       else
          call checkvalbuf(vxyzu(1,j),vg,tolvg,'vg',nerr(3),ncheck(3),errmax(3))
          call checkvalbuf(fxyzu(1,j),-fd,tolfg,'fg',nerr(4),ncheck(4),errmax(4))
       endif
    enddo
    !call checkval(npart/2-1,vxyzu(1,1:npart),vg,tolvg,nerr(2),'vg')
    ekin_exact = 0.5*totmass*(vd**2 + vg**2)
    !print*,' step ',i,'t = ',t,' ekin should be ',ekin_exact, ' got ',ekin,(ekin-ekin_exact)/ekin_exact
    call checkvalbuf(ekin,ekin_exact,tol,'ekin',nerr(5),ncheck(5),errmax(5))
 enddo
 call checkvalbuf_end('dust velocities match exact solution',ncheck(1),nerr(1),errmax(1),tol)
 call checkvalbuf_end('dust accel matches exact solution',   ncheck(2),nerr(2),errmax(2),tolfd)
 call checkvalbuf_end('gas velocities match exact solution',ncheck(3),nerr(3),errmax(3),tolvg)
 call checkvalbuf_end('gas accel matches exact solution',   ncheck(4),nerr(4),errmax(4),tolfg)
 call checkvalbuf_end('kinetic energy decay matches exact',ncheck(5),nerr(5),errmax(5),tol)
#ifdef DUSTGROWTH
 call checkvalbuf_end('interpolated dv matches exact solution',ncheck(6),nerr(6),errmax(6),toldv)
#endif

#ifdef DUSTGROWTH
 call update_test_scores(ntests,nerr(1:6),npass)
#else
 call update_test_scores(ntests,nerr(1:5),npass)
#endif

end subroutine test_dustybox

!----------------------------------------------------
!+
!  3D dust diffusion test from Price & Laibe (2015)
!+
!----------------------------------------------------
subroutine test_dustydiffuse(ntests,npass)
 use dim,       only:maxp,periodic,maxtypes,mhd,use_dust,maxdustsmall,maxalpha
 use part,      only:hfact,npart,npartoftype,massoftype,igas,dustfrac,ddustevol,dustevol,&
                     xyzh,vxyzu,Bevol,dBevol,divcurlv,divcurlB,fext,fxyzu,&
                     set_particle_type,rhoh,temperature,dustprop,ddustprop,&
                     ndusttypes,ndustsmall,alphaind
 use kernel,    only:hfact_default
 use eos,       only:gamma,polyk,ieos
 use dust,      only:K_code,idrag
 use boundary,  only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax,set_boundary
 use io,        only:iverbose
 use unifdis,   only:set_unifdis
 use deriv,     only:derivs
 use testutils, only:checkvalbuf,checkvalbuf_end
 use mpiutils,  only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: nx,j,i,n,nsteps
 integer :: nerr(1),ncheck(1)
 integer :: eps_type
 real    :: errmax(1)
 real    :: deltax,rhozero,totmass,dt,dtnew,time,tmax
 real    :: epstot,epsi(maxdustsmall),rc,rc2,r2,A,B,eta
 real    :: erri,exact,errl2,term,tol
 real,allocatable   :: ddustevol_prev(:,:)
 logical, parameter :: do_output = .false.
 real,    parameter :: t_write(5) = (/0.1,0.3,1.0,3.0,10.0/)

 if (use_dust .and. periodic) then
    if (id==master) write(*,"(/,a)") '--> testing DUSTYDIFFUSE'
 else
    if (id==master) write(*,"(/,a)") '--> skipping DUSTYDIFFUSE (need -DDUST and -DPERIODIC)'
    return
 endif
 !
 ! setup uniform box
 !
 nx = 32
 deltax = 1./nx
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.5,0.5)
 hfact = hfact_default
 rhozero = 3.
 totmass = rhozero*dxbound*dybound*dzbound
 time  = 0.
 npart = 0
 npartoftype(:) = 0
 ndustsmall = maxdustsmall
 ndusttypes = ndustsmall
 iverbose = 2
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,verbose=.false.)
 npartoftype(igas) = npart
 npartoftypetot(igas) = reduceall_mpi('+',npartoftype(igas))
 massoftype(igas)  = totmass/npartoftypetot(igas)
 allocate(ddustevol_prev(ndustsmall,npart))
 vxyzu = 0.
 if (maxalpha==maxp) alphaind(:,:) = 0.
 if (mhd) Bevol = 0.
 !
 ! runtime options
 !
 K_code = 0.1
 ieos = 1
 idrag = 3
 polyk = 1.
 gamma = 1.
 iverbose = 0

 !
 ! setup dust fraction in the box
 !
 epstot = 0.1

 eps_type = 2
 select case(eps_type)
 case(1)
    !--Equal dust fractions
    epsi(:) = epstot/real(ndustsmall)
 case(2)
    !--Unequal dust fractions
    epsi = 0.
    do i=1,ndustsmall
       epsi(i) = 1./real(i)
    enddo
    if (ndustsmall > 0) epsi = epstot/sum(epsi)*epsi
 case default
    print*,'ERROR: eps_type not valid!'
    return
 end select

 !--check that individual dust fractions add up to the total dust fraction
 nerr = 0
 call checkval(sum(epsi),epstot,1.e-14,nerr(1),'sum(epsilon_k) = epsilon')
 ntests = ntests + 1
 if (nerr(1)==0) npass = npass + 1

 rc   = 0.25
 rc2  = rc**2
 dustfrac = 0.
 do i=1,npart
    r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
    if (r2 < rc2) then
       dustfrac(1:ndustsmall,i) = epsi(:)*(1. - r2/rc2)
    endif
    call set_particle_type(i,igas)
 enddo

 ! factors in exact solution (Eq. 51 in PL15)
 B = rc2/epstot
 A = epstot*B**0.6
 eta = 0.1

 !
 ! evolve only the dust fraction with a simple predictor-corrector scheme
 !
 dt = 0.05
 tmax = 10.
 nsteps = nint(tmax/dt)
 dt = tmax/nsteps
 fxyzu = 0.
 fext = 0.
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
             dustfrac,ddustevol,temperature,time,dt,dtnew)

 if (do_output) call write_file(time,xyzh,dustfrac,npart)
 do i=1,npart
!------------------------------------------------
!--sqrt(rho*epsilon) method
!    dustevol(:,i) = sqrt(dustfrac(1:ndustsmall,i)*rhoh(xyzh(4,i),massoftype(igas)))
!------------------------------------------------
!--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
    dustevol(:,i) = sqrt(dustfrac(1:ndustsmall,i)/(1.-dustfrac(1:ndustsmall,i)))
!------------------------------------------------
!--asin(sqrt(epsilon)) method
!    dustevol(:,i) = asin(sqrt(dustfrac(1:ndustsmall,i)))
!------------------------------------------------
 enddo

 nerr = 0
 ncheck = 0
 errmax = 0.
 do j=1,nsteps
    time = j*dt
    !$omp parallel do private(i)
    do i=1,npart
       ddustevol_prev(:,i) = ddustevol(:,i)
       dustevol(:,i) = dustevol(:,i) + dt*ddustevol(:,i)
!------------------------------------------------
!--sqrt(rho*epsilon) method
!       dustfrac(1:ndustsmall,i) = dustevol(:,i)**2/rhoh(xyzh(4,i),massoftype(igas))
!------------------------------------------------
!--sqrt(epsilon/1-epsilon) method (Ballabio et al. 2018)
       dustfrac(1:ndustsmall,i) = dustevol(:,i)**2/(1.+dustevol(:,i)**2)
!------------------------------------------------
!--asin(sqrt(epsilon)) method
!       dustfrac(1:ndustsmall,i) = sin(dustevol(:,i))**2
!------------------------------------------------
    enddo
    !$omp end parallel do
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
                dustprop,ddustprop,dustfrac,ddustevol,temperature,time,dt,dtnew)
    !$omp parallel do private(i)
    do i=1,npart
       dustevol(:,i) = dustevol(:,i) + 0.5*dt*(ddustevol(:,i) - ddustevol_prev(:,i))
    enddo
    !$omp end parallel do

    !
    ! check solution matches the analytic solution at each timestep
    !
    term = 10.*eta*time + B
    n = 0
    errl2 = 0.
    !$omp parallel do private(i,r2,exact,erri) reduction(+:errl2,n)
    do i=1,npart
       r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
       exact = A*abs(term)**(-0.6) - r2/term
       if (exact > 0.) then
          erri  = sum(dustfrac(:,i)) - exact
          errl2 = errl2 + erri*erri
          n = n + 1
       endif
    enddo
    !$omp end parallel do
    errl2 = sqrt(errl2/n)
    tol = 2.6e-3 !1.5e-3/(1. + time)  ! take tolerance down with time
    call checkvalbuf(errl2,0.,tol,'L2 err',nerr(1),ncheck(1),errmax(1))
    !
    ! write solution to file if necessary
    !
    if (do_output .and. any(abs(t_write-time) < 0.01*dt)) call write_file(time,xyzh,dustfrac,npart)
 enddo
 call checkvalbuf_end('dust diffusion matches exact solution',ncheck(1),nerr(1),errmax(1),tol)
 call update_test_scores(ntests,nerr(1:1),npass)

 !
 ! clean up dog poo
 !
 dustevol  = 0.
 dustfrac  = 0.
 ddustevol = 0.

end subroutine test_dustydiffuse

!---------------------------------------------------------------------------------
!+
!  check that drag implementation conserves momentum, angular momentum and energy
!+
!---------------------------------------------------------------------------------
subroutine test_drag(ntests,npass)
 use dim,         only:maxp,periodic,maxtypes,mhd,maxvxyzu,maxdustlarge,maxalpha,use_dustgrowth
 use part,        only:hfact,npart,npartoftype,massoftype,igas,dustfrac,ddustevol,&
                       xyzh,vxyzu,Bevol,dBevol,divcurlv,divcurlB,fext,fxyzu,&
                       set_particle_type,rhoh,temperature,dustprop,ddustprop,&
                       idust,iphase,iamtype,ndusttypes,grainsize,graindens,alphaind
 use options,     only:use_dustfrac
 use eos,         only:polyk,ieos
 use kernel,      only:hfact_default
 use dust,        only:K_code,idrag
 use boundary,    only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax,set_boundary
 use io,          only:iverbose
 use unifdis,     only:set_unifdis
 use deriv,       only:derivs
 use mpiutils,    only:reduceall_mpi
 use random,      only:ran2
 use vectorutils, only:cross_product3D
 use units,       only:udist,unit_density
 integer, intent(inout) :: ntests,npass
 integer(kind=8) :: npartoftypetot(maxtypes)
 integer :: nx,i,j,nfailed(7),itype,iseed,npart_previous,iu
 real    :: da(3),dl(3),temp(3)
 real    :: psep,time,rhozero,totmass,dtnew,dekin,deint

 if (id==master) write(*,"(/,a)") '--> testing DUST DRAG'
!
! set up particles in random distribution
!
 nx = 50
 psep = 1./nx
 iseed= -14255
 call set_boundary(xmin,xmax,ymin,ymax,zmin)
 hfact = hfact_default
 rhozero = 3.
 totmass = rhozero*dxbound*dybound*dzbound
 time  = 0.
 npart = 0
 npartoftype(:) = 0
 if (maxvxyzu < 4) then
    ieos = 1
    polyk = 1.
 else
    ieos = 2
 endif
 fxyzu(:,:) = 0.
 if (maxalpha==maxp) alphaind(:,:) = 0.

 iverbose = 2
 use_dustfrac = .false.
 iu = 4

 call set_unifdis('random',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                      psep,hfact,npart,xyzh,verbose=.false.)
 npartoftype(igas) = npart
 npartoftypetot(igas) = reduceall_mpi('+',npartoftype(igas))
 massoftype(igas) = totmass/npartoftypetot(igas)

 do i=1,npart
    call set_particle_type(i,igas)
    vxyzu(1:3,i) = (/ran2(iseed),ran2(iseed),ran2(iseed)/)
    if (maxvxyzu >= 4) vxyzu(iu,i) = ran2(iseed)
 enddo

 ndusttypes = maxdustlarge
 do j=1,ndusttypes
    grainsize(j) = j*1./udist
    graindens(j) = 1./unit_density
    npart_previous = npart
    call set_unifdis('random',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                         3.*psep,hfact,npart,xyzh,verbose=.false.)

    itype = idust + j - 1
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
       vxyzu(1:3,i) = (/ran2(iseed),ran2(iseed),ran2(iseed)/)
       if (maxvxyzu >= 4) vxyzu(iu,i) = 0.
    enddo
    npartoftype(itype) = npart - npart_previous
    npartoftypetot(itype) = reduceall_mpi('+',npartoftype(itype))
    massoftype(itype) = totmass/npartoftypetot(itype)
 enddo

 if (mhd) Bevol = 0.
 if (use_dustgrowth) then
    dustprop(:,:) = 0.
    dustprop(1,:) = grainsize(1)
    dustprop(2,:) = graindens(1)
 endif
!
! call derivatives
!
 idrag=1
 if (idrag==2) K_code = 100.

 fext = 0.
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
             dustfrac,ddustevol,temperature,time,0.,dtnew)

!
! check that momentum and energy are conserved
!
 da(:) = 0.
 dl(:) = 0.
 dekin = 0.
 deint = 0.
 do i=1,npart
    itype = iamtype(iphase(i))
    da(:) = da(:) + massoftype(itype)*fxyzu(1:3,i)
    if (.not.periodic) then !the angular momentum is not conserved for particle systems with periodic boundary conditions
       call cross_product3D(xyzh(1:3,i),fxyzu(1:3,i),temp)
       dl(:) = dl(:) + massoftype(itype)*temp(:)
    endif
    if (maxvxyzu >= 4) then
       dekin  = dekin  + massoftype(itype)*dot_product(vxyzu(1:3,i),fxyzu(1:3,i))
       deint  = deint  + massoftype(itype)*fxyzu(iu,i)
    endif
 enddo

 da = reduceall_mpi('+', da)
 dl = reduceall_mpi('+', dl)
 dekin = reduceall_mpi('+', dekin)
 deint = reduceall_mpi('+', deint)

 nfailed=0
 call checkval(da(1),0.,7.e-7,nfailed(1),'acceleration from drag conserves momentum(x)')
 call checkval(da(2),0.,7.e-7,nfailed(2),'acceleration from drag conserves momentum(y)')
 call checkval(da(3),0.,7.e-7,nfailed(3),'acceleration from drag conserves momentum(z)')
 if (.not.periodic) then
    call checkval(dl(1),0.,1.e-9,nfailed(4),'acceleration from drag conserves angular momentum(x)')
    call checkval(dl(2),0.,1.e-9,nfailed(5),'acceleration from drag conserves angular momentum(y)')
    call checkval(dl(3),0.,1.e-9,nfailed(6),'acceleration from drag conserves angular momentum(z)')
 endif
 if (maxvxyzu >= 4) call checkval(dekin+deint,0.,7.e-7,nfailed(7),'acceleration from drag conserves energy')

 call update_test_scores(ntests,nfailed,npass)

end subroutine test_drag

!---------------------------------------------------------
!+
!  check that the Epstein/Stokes transition is continuous
!+
!---------------------------------------------------------
subroutine test_epsteinstokes(ntests,npass)
 use dust,      only:idrag,get_ts
 use units,     only:unit_density,unit_velocity,utime,udist
 use physcon,   only:years,kb_on_mh,pi
 use testutils, only:checkval,checkvalbuf,checkvalbuf_end
 integer, intent(inout) :: ntests,npass
 integer :: iregime,i,j,nfailed(1),ncheck
 integer, parameter :: npts=1001, nrhopts = 11
 real :: rhogas,spsoundi,tsi,ts1,deltav,tol,grainsizei,graindensi
 real :: smin,smax,ds,rhomin,rhomax,drho,psi,exact,err,errmax
 logical :: write_output = .false.
 character(len=60) :: filename
 integer, parameter :: lu = 36

 if (id==master) write(*,"(/,a)") '--> testing Epstein/Stokes drag transition'

 idrag = 1
 spsoundi = 6.e4/unit_velocity
 deltav = 1.e-2*spsoundi
 graindensi = 1./unit_density
 !print*,' T = ',(6.e4)**2/kb_on_mh*2.,' delta v/cs = ',deltav/spsoundi
 smin = 1.e-6 ! cgs units
 smax = 1.e8
 rhomin = 1.e-19  ! cgs units
 rhomax = 1.e-09
 ds = (log10(smax) - log10(smin))/real(npts-1)
 drho = (log10(rhomax) - log10(rhomin))/real(nrhopts-1)

 do j=1,nrhopts
    rhogas = rhomin*10**((j-1)*drho)/unit_density
    if (write_output) then
       write(filename,"(a,1pe8.2,a)") 'ts-rho',rhogas*unit_density,'.out'
       open(unit=lu,file=filename,status='replace')
       print "(a)",' writing '//trim(filename)
    endif
    write(filename,"(a,1pe8.2)") 'rho=',rhogas*unit_density
    ncheck = 0
    nfailed = 0
    errmax = 0.
    tol = 6.3e-2
    do i=1,npts
       grainsizei = smin*10**((i-1)*ds)/udist
       !--no need to test drag transition 'ndusttypes' times...once is enough
       call get_ts(idrag,grainsizei,graindensi,rhogas,0.,spsoundi,deltav**2,tsi,iregime)
       !print*,'s = ',grainsizei,' ts = ',tsi*utime/years,',yr ',iregime

       if (i > 1) call checkvalbuf((tsi-ts1)/abs(tsi),0.,tol,'ts is continuous into Stokes regime',nfailed(1),ncheck,errmax)
       ts1 = tsi
       if (write_output) write(lu,*) grainsizei,tsi*utime/years,iregime
    enddo
    if (write_output) close(lu)
    call checkvalbuf_end('ts is continuous into Stokes regime: '//trim(filename),ncheck,nfailed(1),errmax,tol)
    call update_test_scores(ntests,nfailed(1:1),npass)
 enddo

 !
 ! graph of variation with delta v / cs
 !
 grainsizei = 0.1/udist
 rhogas = 1.e-13/unit_density
 smin = 0.
 smax = 10.
 ds = (smax-smin)/real(npts-1)
 err = 0.
 if (write_output) open(unit=lu,file='ts-deltav.out',status='replace')
 do i=1,npts
    deltav = (smin + (i-1)*ds)*spsoundi
    call get_ts(idrag,grainsizei,graindensi,rhogas,0.,spsoundi,deltav**2,tsi,iregime)
    psi = sqrt(0.5)*deltav/spsoundi
    if (i==1) then
       ts1 = tsi
    else
       ! exact non-linear Epstein drag formula
       exact = 20./3./sqrt(pi)/(deltav/spsoundi* &
                ((1./psi + 1./(2.*psi**3))*exp(-psi**2) &
               + (1. + 1./psi**2 - 1./(4.*psi**4))*sqrt(pi)*erf(psi)))
       err = err + (tsi/ts1 - exact)**2
    endif
    if (write_output) write(lu,*) deltav/spsoundi,tsi/ts1,iregime
 enddo
 err = sqrt(err/npts)
 call checkval(err,0.,3.9e-3,nfailed(1),'Epstein drag formula matches non-linear solution')
 call update_test_scores(ntests,nfailed(1:1),npass)

 if (write_output) close(lu)

end subroutine test_epsteinstokes

!---------------------------------------------------
!+
!  write an output file with r, dustfrac in table
!  this is what is shown in Phantom paper
!+
!---------------------------------------------------
subroutine write_file(time,xyzh,dustfrac,npart)
 use dim,  only:maxp
 use part, only:ndusttypes
 real, intent(in)     :: time
 real, intent(in)    :: xyzh(:,:),dustfrac(:,:)
 integer, intent(in) :: npart
 character(len=30)   :: filename,str1,str2,fmt1
 integer :: i,lu
 real    :: r2,dustfracsum(maxp)

 write(str1,"(f5.1)") time
 if (ndusttypes>1) then
    write(str2,"(I5)") ndusttypes+2
    dustfracsum = sum(dustfrac,1)
 else
    write(str2,"(I5)") ndusttypes+1
 endif
 filename = 'dustfrac_t'//trim(adjustl(str1))//'.txt'
 fmt1 = '('//trim(adjustl(str2))//'F15.8)'
 open(newunit=lu,file=filename,status='replace')
 print*,' writing '//filename
 write(lu,*) time
 do i=1,npart
    r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
    if (ndusttypes>1) then
       write(lu,fmt1) sqrt(r2),dustfracsum(i),dustfrac(:,i)
    else
       write(lu,fmt1) sqrt(r2),dustfrac(:,i)
    endif
 enddo
 close(lu)

end subroutine write_file

#endif

end module testdust
