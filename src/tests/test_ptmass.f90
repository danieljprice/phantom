!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testptmass
!
! Unit tests of the ptmass/sink particles module
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: HIIRegion, boundary, centreofmass, checksetup, cons2prim,
!   deriv, dim, energies, eos, eos_HIIR, extern_binary, extern_gr,
!   externalforces, gravwaveutils, io, kdtree, kernel, metric,
!   metric_tools, mpiutils, options, part, physcon, ptmass, random,
!   setbinary, setdisc, setup_params, spherical, step_lf_global,
!   stretchmap, subgroup, substepping, testutils, timestep, timing, units
!
 use testutils, only:checkval,update_test_scores
 implicit none
 public :: test_ptmass

 private

contains

subroutine test_ptmass(ntests,npass,string)
 use io,      only:id,master,iskfile
 use eos,     only:polyk,gamma
 use part,    only:nptmass,gr
 use options, only:iexternalforce,alpha
 use ptmass,  only:use_fourthorder,set_integration_precision
 character(len=*), intent(in) :: string
 character(len=20) :: filename
 character(len=40) :: stringf
 integer, intent(inout) :: ntests,npass
 integer :: itmp,ierr,itest,istart
 logical :: do_test_binary,do_test_accretion,do_test_createsink,do_test_softening
 logical :: do_test_chinese_coin,do_test_merger,do_test_potential,do_test_HII,do_test_SDAR
 logical :: do_test_binary_gr
 logical :: testall

 if (id==master) write(*,"(/,a,/)") '--> TESTING PTMASS MODULE'

 do_test_binary = .false.
 do_test_accretion = .false.
 do_test_createsink = .false.
 do_test_softening = .false.
 do_test_merger = .false.
 do_test_potential = .false.
 do_test_chinese_coin = .false.
 do_test_HII = .false.
 do_test_SDAR = .false.
 do_test_binary_gr = .false.
 testall = .false.
 istart = 1
 select case(trim(string))
 case('ptmassbinary')
    do_test_binary = .true.
 case('ptmassgenrel')
    do_test_binary_gr = .true.
 case('ptmassaccrete')
    do_test_accretion = .true.
 case('ptmasscreatesink')
    do_test_createsink = .true.
 case('ptmasssoftening')
    do_test_softening = .true.
 case('ptmassmerger')
    do_test_merger = .true.
 case('ptmasspotential')
    do_test_potential = .true.
 case('ptmasschinchen','ptmasscoin','chinchen','coin','chinesecoin')
    do_test_chinese_coin = .true.
 case('ptmassfsi','fsi')
    istart = 2
    do_test_binary = .true.
    do_test_softening = .true.
    do_test_merger = .true.
 case('ptmassHII')
    do_test_HII = .true.
 case('ptmassSDAR')
    do_test_SDAR = .true.
 case default
    testall = .true.
 end select
 !
 !--general settings
 !
 polyk = 0.
 gamma = 1.
 iexternalforce = 0
 alpha = 0.01
 !
 !  Test for sink particles in GR
 !
 if ((do_test_binary_gr .or. testall) .and. gr) then
    call test_sink_binary_gr(ntests,npass,string)
    return
 endif

 do itest=istart,2
    !
    !  select order of integration
    !
    if (itest == 2) then
       use_fourthorder = .true.
       stringf = ' with Forward Symplectic Integrator'
    else
       use_fourthorder = .false.
       stringf = ' with Leapfrog integrator'
    endif
    call set_integration_precision
    !
    !  Tests of a sink particle binary
    !
    if (do_test_binary .or. testall) call test_binary(ntests,npass,stringf)
    !
    !  Test of softening between sinks
    !
    if (do_test_softening .or. testall) call test_softening(ntests,npass)
    !
    !  Test of Chinese Coin problem
    !
    if (do_test_chinese_coin .or. testall) call test_chinese_coin(ntests,npass,stringf)
    !
    !  Test sink particle mergers
    !
    if (do_test_merger .or. testall) call test_merger(ntests,npass)

 enddo
 !
 !  Test of sink particle potentials
 !
 if (do_test_potential .or. testall) call test_sink_potential(ntests,npass)
 !
 !  Tests of accrete_particle routine
 !
 if (do_test_accretion .or. testall) then
    do itest=1,2
       call test_accretion(ntests,npass,itest)
    enddo
 endif
 !
 !  Test sink particle creation
 !
 if (do_test_createsink .or. testall) call test_createsink(ntests,npass)

 if (do_test_SDAR .or. testall) call test_SDAR(ntests,npass)

 if (do_test_HII) call test_HIIregion(ntests,npass)


 !reset stuff and clean up temporary files
 itmp    = 201
 nptmass = 0
 close(iskfile,iostat=ierr)
 write(filename,"(i3)") iskfile
 filename = 'fort.'//trim(adjustl(filename))
 open(unit=itmp,file=filename,status='old',iostat=ierr)
 close(itmp,status='delete',iostat=ierr)
 open(unit=itmp,file='Sink00.sink',status='old',iostat=ierr)
 close(itmp,status='delete',iostat=ierr)

 if (id==master) write(*,"(/,a)") '<-- PTMASS TEST COMPLETE'

end subroutine test_ptmass

!-----------------------------------------------------------------------
!+
!  Unit tests of a sink particle binary orbit
!+
!-----------------------------------------------------------------------
subroutine test_binary(ntests,npass,string)
 use dim,        only:periodic,gravity,ind_timesteps,use_sinktree
 use io,         only:id,master,iverbose
 use physcon,    only:pi,deg_to_rad
 use ptmass,     only:get_accel_sink_sink,h_soft_sinksink, &
                      get_accel_sink_gas,f_acc,use_fourthorder
 use part,       only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,fext,&
                      npart,npartoftype,massoftype,xyzh,vxyzu,fxyzu,&
                      hfact,igas,epot_sinksink,init_part,iJ2,ispinx,ispiny,ispinz,iReff,istar,&
                      shortsinktree,fxyz_ptmass_tree,ihsoft
 use energies,   only:angtot,etot,totmom,compute_energies,hp,hx
 use timestep,   only:dtmax,C_force,tolv
 use kdtree,     only:tree_accuracy
 use eos,        only:gamma,ieos,polyk
 use setbinary,  only:set_binary
 use setdisc,    only:set_disc
 use units,      only:set_units
 use mpiutils,   only:bcast_mpi,reduce_in_place_mpi
 use step_lf_global, only:init_step,step
 use gravwaveutils,  only:get_strain_from_circular_binary,get_G_on_dc4,calc_gravitwaves
 use testutils,      only:checkvalf,checkvalbuf,checkvalbuf_end
 use checksetup,     only:check_setup
 use deriv,          only:get_derivs_global
 use timing,         only:getused,printused
 use options,        only:ipdv_heating,ishock_heating,iexternalforce
 use externalforces, only:iext_corotate,omega_corotate,externalforce_vdependent
 integer,          intent(inout) :: ntests,npass
 character(len=*), intent(in)    :: string
 integer :: i,ierr,itest,nfailed(3),nsteps,nerr,nwarn,norbits
 integer :: merge_ij(2),merge_n,nparttot,nfailgw(2),ncheckgw(2)
 integer, parameter :: nbinary_tests = 6
 real :: m1,m2,a,ecc,hacc1,hacc2,dt,dtext,t,dtnew,tolen,tolmom,tolang,hp_exact,hx_exact
 real :: angmomin,etotin,totmomin,dum,dum2,omega,errmax,dtsinksink,fac,errgw(2)
 real :: angle,rin,rout
 real :: fxyz_sinksink(4,2),dsdt_sinksink(3,2) ! we only use 2 sink particles in the tests here
 real(kind=4) :: t1
 character(len=20) :: dumpfile
 real, parameter :: tolgw = 1.2e-2
 !
 !--no gas particles
 !
 call init_part()
 iverbose = 0
 tree_accuracy = 0.
 h_soft_sinksink = 0.
 calc_gravitwaves = .true.
 ipdv_heating = 0
 ishock_heating = 0

 tolv = 1e-2

 binary_tests: do itest = 1,nbinary_tests
    select case(itest)
    case(4)
       if (use_fourthorder) then
          if (id==master) write(*,"(/,a)") '--> skipping integration of binary orbit with oblateness with FSI'
          cycle binary_tests
       else
          if (id==master) write(*,"(/,a)") '--> testing integration of binary orbit with oblateness'//trim(string)
       endif
    case(2,3,5)
       if (periodic) then
          if (id==master) write(*,"(/,a)") '--> skipping circumbinary disc test (-DPERIODIC is set)'
          cycle binary_tests
       elseif (use_fourthorder .and. itest==5) then
          if (id==master) write(*,"(/,a)") '--> skipping circumbinary disc around oblate star test with FSI'
          cycle binary_tests
       elseif (use_sinktree .and. itest==5) then
          if (id==master) write(*,"(/,a)") '--> skipping circumbinary disc around oblate star test with Sinktree'
          cycle binary_tests
       else
          if (itest==5) then
             if (id==master) write(*,"(/,a)") '--> testing integration of disc around oblate star'//trim(string)
          elseif (itest==3) then
             if (id==master) write(*,"(/,a)") '--> testing integration of disc around eccentric binary'//trim(string)
          else
             if (id==master) write(*,"(/,a)") '--> testing integration of circumbinary disc'//trim(string)
          endif
       endif
    case(6)
       if (id==master) write(*,"(/,a)") '--> testing integration of binary orbit in a corotating frame'//trim(string)
    case default
       if (id==master) write(*,"(/,a)") '--> testing integration of binary orbit'//trim(string)
    end select
    !
    !--setup sink-sink binary (no gas particles)
    !
!      time = 0.
    npart = 0
    npartoftype = 0
    nptmass = 0
    m1    = 1.
    m2    = 1.
    a     = 1.
    rin   = 1.5*a
    rout  = 15.*a
    if (itest==5) then
       m2 = 0.0
       rin = 1.
       rout = 5.
    endif
    if (itest==3 .or. itest==4) then
       ecc = 0.5
    else
       ecc = 0.
    endif
    hacc1  = 0.35
    hacc2  = 0.35
    C_force = 0.25
    t = 0.
    if (itest==3) C_force = 0.25
    omega = sqrt((m1+m2)/a**3)
    call set_units(mass=1.d0,dist=1.d0,G=1.d0)
    if (itest==6) then
       if (use_fourthorder) cycle binary_tests ! corotating frame currently does not work with 4th order scheme
       iexternalforce = iext_corotate
       call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                       verbose=.false.)
    else
       iexternalforce = 0
       call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,verbose=.false.)
    endif

    if (ierr /= 0) nerr = nerr + 1

    if (itest==2 .or. itest==3 .or. itest==5) then
       !  add a circumbinary gas disc around it
       nparttot = 1000
       call set_disc(id,master,nparttot=nparttot,npart=npart,rmin=rin,rmax=rout,p_index=1.0,q_index=0.75,&
                     HoverR=0.1,disc_mass=0.01*m1,star_mass=m1+m2,gamma=gamma,&
                     particle_mass=massoftype(igas),hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,&
                     polyk=polyk,verbose=.false.)
       npartoftype(igas) = npart
    endif

    if (use_sinktree) then
       shortsinktree = 1
       fxyz_ptmass_tree = 0.
       if (itest /= 4) then
          xyzmh_ptmass(ihsoft,1) = hacc1
          xyzmh_ptmass(ihsoft,2) = hacc1
       endif
    endif

    if (itest==4 .or. itest==5) then
       if (itest==5) nptmass = 1
       ! set oblateness
       xyzmh_ptmass(iJ2,1) = 0.01629
       angle = 45.*deg_to_rad
       xyzmh_ptmass(ispinx,1) = 1e2*sin(angle)
       xyzmh_ptmass(ispiny,1) = 0.
       xyzmh_ptmass(ispinz,1) = 1e2*cos(angle)
       xyzmh_ptmass(iReff,1) = hacc1
    endif
    !
    ! check that no errors occurred when setting up initial conditions
    !
    nfailed = 0
    call check_setup(nerr,nwarn)
    call checkval(nerr,0,0,nfailed(1),'no errors during disc setup')
    call update_test_scores(ntests,nfailed,npass)

    tolv = 1.e-2
    iverbose = 0
    ieos = 3
    fac = 1./get_G_on_dc4()
    !
    !--compute SPH forces
    !
    if (npart > 0) then
       fxyzu(:,:) = 0.
       call get_derivs_global()
    endif
    !
    ! initialise forces
    !
    if (id==master) then
       call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_sinksink,epot_sinksink,&
                                dtsinksink,iexternalforce,0.,merge_ij,merge_n,dsdt_sinksink)
    endif
    fxyz_ptmass(:,1:nptmass) = 0.
    dsdt_ptmass(:,1:nptmass) = 0.
    call bcast_mpi(epot_sinksink)
    call bcast_mpi(dtsinksink)

    fext(:,:) = 0.
    if (.not. use_sinktree) then
       do i=1,npart
          call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                fext(1,i),fext(2,i),fext(3,i),dum,massoftype(igas),fxyz_ptmass,dsdt_ptmass,dum,dum2)
       enddo
    endif
    if (id==master) then
       fxyz_ptmass(:,1:nptmass) = fxyz_ptmass(:,1:nptmass) + fxyz_sinksink(:,1:nptmass)
       dsdt_ptmass(:,1:nptmass) = dsdt_ptmass(:,1:nptmass) + dsdt_sinksink(:,1:nptmass)
    endif
    call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))
    call reduce_in_place_mpi('+',dsdt_ptmass(:,1:nptmass))

    !
    !--take the sink-sink timestep specified by the get_forces routine
    !
    dt = C_force*dtsinksink
    if (m2 <= 0.) then
       dt = min(C_force*dtsinksink,4.e-3*sqrt(2.*pi/omega))
    elseif (itest==6) then
       dt = 1.25e-2    !time step of the system on a not corotating frame
    endif

    dtmax = dt  ! required prior to derivs call, as used to set ibin
    !
    !--evolve this for a number of orbits
    !
    call compute_energies(t)
    etotin   = etot
    totmomin = totmom
    angmomin = angtot

    !
    !--check that initial potential on the two sinks is correct
    !
    nfailed(:) = 0
    if (itest==1) then
       call checkval(epot_sinksink,-m1*m2/a,epsilon(0.),nfailed(1),'potential energy')
       call update_test_scores(ntests,nfailed,npass)
       !
       !--check initial angular momentum on the two sinks is correct
       !
       call checkval(angtot,m1*m2*sqrt(a/(m1 + m2)),1e6*epsilon(0.),nfailed(1),'angular momentum')
       call update_test_scores(ntests,nfailed,npass)
    endif
    !
    !--determine number of steps per orbit for information
    !
    nsteps  = int(2.*pi/omega/dt) + 1
    if (itest==2 .or. itest==3 .or. itest==5) then
       norbits = 10
    else
       norbits = 100
    endif
    if (id==master) print*,'steps/orbit = ',nsteps,' norbits = ',norbits,' dt = ',dt
    nsteps = nsteps*norbits
    errmax = 0.; errgw = 0.
    nfailgw = 0; ncheckgw = 0
    dumpfile='test_00000'
    f_acc = 1.
    if (id==master) call getused(t1)
    call init_step(npart,t,dtmax)
    do i=1,nsteps
       dtext = dt
       if (id==master .and. iverbose > 2) write(*,*) ' t = ',t,' dt = ',dt
       call step(npart,npart,t,dt,dtext,dtnew)
       call compute_energies(t)
       errmax = max(errmax,abs(etot - etotin))
       !if (itest==3) print*,t,abs(angtot-angmomin)/angmomin
       !
       ! Check the gravitational wave strain if the binary is circular.
       ! There is a phase error that grows with time, so only check the first 10 orbits
       !
       if (calc_gravitwaves .and. abs(ecc) < epsilon(ecc) .and. itest==1 .and. t < 20.*pi/omega) then
          call get_strain_from_circular_binary(t+dt,m1,m2,a,0.,hx_exact,hp_exact)
          call checkvalbuf(10.+hx(1)*fac,10.+hx_exact*fac,tolgw,&
                           'gw strain (x)',nfailgw(1),ncheckgw(1),errgw(1))
          call checkvalbuf(10.+hp(1)*fac,10.+hp_exact*fac,tolgw,&
                           'gw strain (+)',nfailgw(2),ncheckgw(2),errgw(2))
       endif
       t = t + dt
    enddo
    call compute_energies(t)
    if (id==master) call printused(t1)
    nfailed(:) = 0
    tolmom = 2.e-14
    tolang = 2.e-14
    select case(itest)
    case(5)
       tolen = 9.e-1
    case(4)
       tolmom = 1.e-14
       tolen = 1.6e-2
    case(3)
       if (ind_timesteps) then
          tolang = 2.1e-6
       else
          tolang = 6.e-10
       endif
       tolen = 1.2e-2
    case(2)
       tolen = 1.2e-3
       if (gravity) tolen = 3.1e-3

       if (use_fourthorder) then
          tolang = 2.e-11
       endif
    case default
       if (calc_gravitwaves .and. itest==1) then
          call checkvalbuf_end('grav. wave strain (x)',ncheckgw(1),nfailgw(1),errgw(1),tolgw)
          call checkvalbuf_end('grav. wave strain (+)',ncheckgw(2),nfailgw(2),errgw(2),tolgw)
          call update_test_scores(ntests,nfailgw(1:2),npass)
       endif
       if (use_fourthorder) then
          tolen = 1.e-13
       else
          tolen = 3.e-8
       endif
    end select
    !
    !--check energy conservation
    !
    call checkval(angtot,angmomin,tolang,nfailed(1),'angular momentum')
    call checkval(totmom,totmomin,tolmom,nfailed(2),'linear momentum')
    call checkval(etotin+errmax,etotin,tolen,nfailed(3),'total energy')
    do i=1,3
       call update_test_scores(ntests,nfailed(i:i),npass)
    enddo

    ! reset iexternalforce
    iexternalforce = 0
 enddo binary_tests

end subroutine test_binary

!-----------------------------------------------------------------------
!+
!  Test that binary setup in GR using sink particles is OK.
!+
!-----------------------------------------------------------------------
subroutine test_sink_binary_gr(ntests,npass,string)
 use io,             only:id,master,iverbose
 use part,           only:init_part,npart,npartoftype,nptmass,xyzmh_ptmass,vxyz_ptmass,&
                          epot_sinksink,metrics_ptmass,metricderivs_ptmass,pxyzu_ptmass,&
                          fxyz_ptmass,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,&
                          fext
 use timestep,       only:C_force,dtextforce
 use physcon,        only:solarm,pi
 use units,          only:set_units
 use setbinary,      only:set_binary
 use metric,         only:mass1
 use checksetup,     only:check_setup
 use testutils,      only:checkval,checkvalf,update_test_scores
 use ptmass,         only:get_accel_sink_sink
 use metric_tools,   only:init_metric
 use cons2prim,      only:prim2consall
 use extern_gr,      only:get_grforce_all
 use substepping,    only:combine_forces_gr
 use energies,       only:angtot,etot,totmom,compute_energies,epot
 use substepping,    only:substep_gr
 integer, intent(inout)          :: ntests,npass
 character(len=*), intent(in)    :: string
 real :: fxyz_sinksink(4,2),dsdt_sinksink(3,2) ! we only use 2 sink particles in the tests here
 real    :: m1,m2,a,ecc,hacc1,hacc2,t,dt,tol_en
 real    :: dtsinksink,tol,omega,errmax,dis
 real    :: angmomin,etotin,totmomin,dtsph,dtorb,vphi
 integer :: ierr,nerr,nfailed(6),nwarn,nsteps,i,ntypes
 integer :: merge_ij(2),merge_n,norbits
 character(len=20) :: dumpfile
 !
 !--no gas particles
 !
 call init_part()
 !
 !--set quantities
 !
 npartoftype = 0
 npart   = 0
 nptmass = 0
 m1      = 1.e-6
 m2      = 1.e-6
 a       = 2.35 ! udist in GR is 1.48E+11. 5 Rsun in code units
 ecc     = 0.   ! eccentricity of binary orbit
 hacc1   = 0.75 ! 0.35 rsun in code units
 hacc2   = 0.75
 mass1   = 0.   ! set BH mass as 0. So the metric becomes Minkowski
 t       = 0.
 iverbose = 0
 ! chose a very small value because a value of 0.35 was resulting in distance - distance_init of 1.e-3
 ! but using a small timestep resulted in values smaller than equal to 1.e-4
 C_force = 0.25
 tol     = epsilon(0.)
 omega   = sqrt((m1+m2)/a**3)
 vphi    = a*omega
 ! set sinks around each other
 call set_units(mass=1.e6*solarm,c=1.d0,G=1.d0)
 call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,verbose=.false.)
 dis = norm2(xyzmh_ptmass(1:3,1) - xyzmh_ptmass(1:3,2))

 if (ierr /= 0) nerr = nerr + 1

 ! check the setup is ok
 nfailed = 0
 call check_setup(nerr,nwarn)
 call checkval(nerr,0,0,nfailed(1),'no errors during setting sink binary orbit')
 call update_test_scores(ntests,nfailed,npass)
 !
 !--initialise forces and test that the curvature contribution is 0. when mass1 is 0.
 !
 if (id==master) then

    call init_metric(nptmass,xyzmh_ptmass,metrics_ptmass,metricderivs_ptmass)
    call prim2consall(nptmass,xyzmh_ptmass,metrics_ptmass,&
                     vxyz_ptmass,pxyzu_ptmass,use_dens=.false.,use_sink=.true.)
    ! sinks in GR, provide external force due to metric to determine the sink total force
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_sinksink,epot_sinksink,&
                           dtsinksink,0,0.,merge_ij,merge_n,dsdt_sinksink)
    call get_grforce_all(nptmass,xyzmh_ptmass,metrics_ptmass,metricderivs_ptmass,&
                     vxyz_ptmass,fxyz_ptmass,dtextforce,use_sink=.true.)
    call combine_forces_gr(nptmass,fxyz_sinksink,fxyz_ptmass)

    ! Test the force calculated is same as sink-sink because there is no curvature.

    call checkval(fxyz_sinksink(1,1), fxyz_ptmass(1,1),tol,nfailed(1),'x force term for sink 1')
    call checkval(fxyz_sinksink(2,1), fxyz_ptmass(2,1),tol,nfailed(2),'y force term for sink 1')
    call checkval(fxyz_sinksink(3,1), fxyz_ptmass(3,1),tol,nfailed(3),'z force term for sink 1')
    call checkval(fxyz_sinksink(1,2), fxyz_ptmass(1,2),tol,nfailed(4),'x force term for sink 2')
    call checkval(fxyz_sinksink(2,2), fxyz_ptmass(2,2),tol,nfailed(5),'y force term for sink 2')
    call checkval(fxyz_sinksink(3,2), fxyz_ptmass(3,2),tol,nfailed(6),'z force term for sink 2')

    call update_test_scores(ntests,nfailed(1:3),npass)
    call update_test_scores(ntests,nfailed(3:6),npass)

 endif
 !
 !--check energy and angular momentum of the system
 !
 dtextforce =  min(C_force*dtsinksink,dtextforce)
 dt    = dtextforce
 call compute_energies(t)
 etotin   = etot
 totmomin = totmom
 angmomin = angtot

 call checkval(epot,-m1*m2/a,epsilon(0.),nfailed(1),'potential energy')
 call update_test_scores(ntests,nfailed,npass)
 !
 !--check initial angular momentum on the two sinks is correct
 !
 call checkval(angtot,m1*m2*sqrt(a/(m1 + m2)),1e6*epsilon(0.),nfailed(2),'angular momentum')
 call update_test_scores(ntests,nfailed,npass)
 !
 !--check initial total energy of the two sinks is correct
 !--using Virial Theorem for the test
 !
 call checkval(etot,epot*0.5,epsilon(0.),nfailed(3),'total energy')
 call update_test_scores(ntests,nfailed,npass)
 !
 !--determine number of steps per orbit for information
 !
 dtorb = 2.*pi/omega
 dt = dtorb
 norbits = 100
 nsteps = norbits*nint(dtorb/dt)
 errmax = 0.
 dumpfile='test_00000'
 ntypes = 2

 do i=1,nsteps
    dtsph = dt
    call substep_gr(npart,nptmass,ntypes,dtsph,dtextforce,xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext,t,&
                       xyzmh_ptmass,vxyz_ptmass,pxyzu_ptmass,metrics_ptmass,metricderivs_ptmass,fxyz_ptmass)
    call compute_energies(t)
    errmax = max(errmax,abs(etot - etotin))
    t = t + dt
    dis = norm2(xyzmh_ptmass(1:3,1) - xyzmh_ptmass(1:3,2))
 enddo
 !
 !--check the radius of the orbit does not change
 !
 call checkval(dis,a,7.e-4,nfailed(1),"radius of orbit")
 call update_test_scores(ntests,nfailed,npass)
 !
 !--check energy, linear and angular momentum conservation
 !
 tol_en = 1.e-13
 call compute_energies(t)
 call checkval(angtot,angmomin,tol_en,nfailed(1),'angular momentum')
 call checkval(totmom,totmomin,tol_en,nfailed(2),'linear momentum')
 call checkval(etotin+errmax,etotin,tol_en,nfailed(3),'total energy')
 do i=1,3
    call update_test_scores(ntests,nfailed(i:i),npass)
 enddo

end subroutine test_sink_binary_gr
!-----------------------------------------------------------------------
!+
!  Test softening between sink particles. Run a binary orbit
!  and check conservation, also check the potential energy is correct
!+
!-----------------------------------------------------------------------
subroutine test_softening(ntests,npass)
 use io,         only:id,master,iverbose
 use physcon,    only:pi
 use testutils,  only:checkval,checkvalf,update_test_scores
 use ptmass,     only:get_accel_sink_sink,h_soft_sinksink, &
                      get_accel_sink_gas
 use part,       only:npart,npartoftype,epot_sinksink,&
                      nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,&
                      shortsinktree,fxyz_ptmass_tree,poten
 use energies,   only:angtot,etot,totmom,compute_energies,epot
 use timestep,   only:dtmax,C_force
 use setbinary,  only:set_binary
 use units,      only:set_units
 use mpiutils,   only:bcast_mpi,reduce_in_place_mpi
 use step_lf_global, only:init_step,step
 use kernel,         only:kernel_softening
 use dim,            only:use_sinktree,maxpsph,maxp
 integer, intent(inout) :: ntests,npass
 integer :: i,ierr,nfailed(4),nsteps,norbits,merge_ij(2),merge_n
 real :: m1,m2,a,ecc,hacc1,hacc2,t,dt,dtext,dtnew,dtsinksink
 real :: q,phisoft,fsoft,mu,v_c1,v_c2,r1,r2,omega1,omega2,omega
 real :: etotin,totmomin,angmomin,errmax

 if (id==master) write(*,"(/,a)") '--> testing softening in sink particle binary'
 nptmass = 0
 npart = 0
 npartoftype = 0
 m1    = 1.
 m2    = 1.
 a     = 1.
 ecc   = 0.
 hacc1 = 0.
 hacc2 = 0.
 t     = 0.
 if (use_sinktree) then
    shortsinktree = 1
    fxyz_ptmass_tree = 0.
    poten(maxpsph+1:maxp) = 0.
 endif

 h_soft_sinksink = 0.8*a

 call set_units(mass=1.d0,dist=1.d0,G=1.d0)
 call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,verbose=.false.)

 q   = a/h_soft_sinksink
 call kernel_softening(q*q,q,phisoft,fsoft)

 ! Test energy and momentum conservation
 mu = (m1*m2)/(m1+m2)
 r1 = sqrt(xyzmh_ptmass(1,1)**2+xyzmh_ptmass(2,1)**2+xyzmh_ptmass(3,1)**2)
 r2 = sqrt(xyzmh_ptmass(1,2)**2+xyzmh_ptmass(2,2)**2+xyzmh_ptmass(3,2)**2)
 omega1 = sqrt(m2*fsoft/(r1*h_soft_sinksink**2))
 omega2 = sqrt(m1*fsoft/(r2*h_soft_sinksink**2))
 v_c1 = omega1*r1
 v_c2 = omega2*r2
 vxyz_ptmass(1,1) = 0.
 vxyz_ptmass(2,1) = v_c1
 vxyz_ptmass(3,1) = 0.
 vxyz_ptmass(1,2) = 0.
 vxyz_ptmass(2,2) = -v_c2
 vxyz_ptmass(3,2) = 0.
 call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                          dtsinksink,0,0.,merge_ij,merge_n,dsdt_ptmass)
 call compute_energies(t)
 etotin   = etot
 totmomin = totmom
 angmomin = angtot

 call checkval(epot,m1*m2*(phisoft)/h_soft_sinksink,2.*epsilon(0.),nfailed(1),'potential energy')
 call update_test_scores(ntests,nfailed(1:1),npass)

 C_force = 0.25
 dt      = 0.3*C_force*dtsinksink
 !if (id==master) print*,' dt for sinks = ',dt
 dtmax   = dt
 omega   = omega1
 nsteps  = int(2.*pi/omega/dt) + 1
 norbits = 10
 if (id==master) print*,' nsteps per orbit = ',nsteps,' norbits = ',norbits
 nsteps = nsteps*norbits
 errmax = 0.
 iverbose = 0
 do i=1,nsteps
    t = t + dt
    dtext = dt
    if (id==master .and. iverbose > 2) write(*,*) ' t = ',t,' dt = ',dt
    call step(npart,npart,t,dt,dtext,dtnew)
!    write(1,'(10(es10.3,1x))') xyzmh_ptmass(:,1)
!    write(1,'(10(es10.3,1x))') xyzmh_ptmass(:,2)
    call compute_energies(t)
    errmax = max(errmax,abs(etot - etotin))
 enddo
 call compute_energies(t)
 nfailed(:) = 0
 call checkval(angtot,angmomin,2.e-14,nfailed(1),'angular momentum')
 call checkval(totmom,totmomin,tiny(0.),nfailed(2),'linear momentum')
 call checkval(etotin+errmax,etotin,2.e-9,nfailed(3),'total energy')
!  call checkval(      ,r_max,1.e-10,nfailed(4),'radius')

 call update_test_scores(ntests,nfailed(1:4),npass)

 ! Reset sink softening
 h_soft_sinksink = 0.0

end subroutine test_softening

!-----------------------------------------------------------------------
!+
!  Test Chinese Coin problem from Chin & Chen (2005)
!+
!-----------------------------------------------------------------------
subroutine test_chinese_coin(ntests,npass,string)
 use io,             only:id,master,iverbose
 use part,           only:xyzmh_ptmass,vxyz_ptmass,ihacc,nptmass,npart,npartoftype,fxyz_ptmass,dsdt_ptmass,&
                          shortsinktree,fxyz_ptmass_tree
 use dim,            only:use_sinktree
 use extern_binary,  only:mass1,mass2
 use options,        only:iexternalforce
 use externalforces, only:iext_binary,update_externalforce
 use physcon,        only:pi
 use step_lf_global, only:step
 use ptmass,         only:use_fourthorder,get_accel_sink_sink
 integer,          intent(inout) :: ntests,npass
 character(len=*), intent(in)    :: string
 character(len=10) :: tag
 integer :: nfailed(3),merge_ij(1),merge_n,norbit
 real :: t,dtorb,dtnew,dtext,tmax,epot_sinksink,y0,v0
 real :: tol_per_orbit_y,tol_per_orbit_v

 if (id==master) write(*,"(/,a)") '--> testing Chinese coin problem'//trim(string)//' (coin)'

 ! no gas
 npart = 0
 npartoftype = 0
 if (use_sinktree) then
    shortsinktree = 1
    fxyz_ptmass_tree = 0.
 endif

 ! add  a single sink particle
 y0 = 0.0580752367; v0 = 0.489765446
 nptmass = 1
 xyzmh_ptmass = 0.
 xyzmh_ptmass(2,1) = y0
 xyzmh_ptmass(4,1) = 1.0
 xyzmh_ptmass(ihacc,1) = 0.1
 vxyz_ptmass = 0.
 vxyz_ptmass(1,1) = v0

 ! external binary
 iexternalforce = iext_binary
 mass1 = 0.5
 mass2 = mass1
 dtorb = 9.*pi
 tmax = 3.*dtorb

 t = 0.
 dtext = 1.e-15
 iverbose = 1
 call update_externalforce(iexternalforce,t,0.)
 call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,&
                          dtext,iexternalforce,t,merge_ij,merge_n,dsdt_ptmass)

 dtext = 1.e-15 ! take small first step
 norbit = 0
 nfailed(:) = 0
 tol_per_orbit_y = 2.5e-2
 tol_per_orbit_v = 1.15e-2
 if (use_fourthorder) then
    tol_per_orbit_y = 1.1e-3
    tol_per_orbit_v = 3.35e-4
 endif
 do while (t < tmax)
    ! do a whole orbit but with the substepping handling how many steps per orbit
    call step(npart,npart,t,dtorb,dtext,dtnew)
    t = t + dtorb
    norbit = norbit + 1

    write(tag,"(a,i1,a)") '(orbit ',norbit,')'
    call checkval(xyzmh_ptmass(2,1),y0,norbit*tol_per_orbit_y,nfailed(1),'y pos of sink '//trim(tag))
    call checkval(vxyz_ptmass(1,1),v0,norbit*tol_per_orbit_v,nfailed(2),'x vel of sink '//trim(tag))
 enddo

 call update_test_scores(ntests,nfailed(1:2),npass)
 iverbose = 0  ! reset verbosity
 iexternalforce = 0

end subroutine test_chinese_coin

!-----------------------------------------------------------------------
!+
!  Test accretion of gas particles onto sink particles
!+
!-----------------------------------------------------------------------
subroutine test_accretion(ntests,npass,itest)
 use io,        only:id,master
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,massoftype, &
                     npart,npartoftype,xyzh,vxyzu,fxyzu,igas,ihacc,&
                     isdead_or_accreted,set_particle_type,ndptmass,hfact
 use ptmass,    only:ptmass_accrete,update_ptmass
 use energies,  only:compute_energies,angtot,etot,totmom
 use mpiutils,  only:bcast_mpi,reduce_in_place_mpi
 use testutils, only:checkval,update_test_scores
 use kernel,    only:hfact_default
 use eos,       only:polyk
 use setdisc,   only:set_disc
 integer, intent(inout) :: ntests,npass
 integer, intent(in)    :: itest
 integer :: i,nfailed(11),np_disc
 integer(kind=1) :: ibin_wakei
 character(len=20) :: string
 logical :: accreted
 real :: t
 real :: dptmass(ndptmass,1)
 real :: dptmass_thread(ndptmass,1)
 real :: angmomin,etotin,totmomin

 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:)  = 0.

 string = 'of two particles'
 if (itest==2) string = 'of a whole disc'
 if (id==master) write(*,"(/,a)") '--> testing accretion '//trim(string)//' onto sink particles'
 nptmass = 1
 !--setup 1 point mass at (-5,-5,-5)
 xyzmh_ptmass(1:3,1)   = 1.
 xyzmh_ptmass(4,1)     = 10. ! mass of sink
 xyzmh_ptmass(ihacc,1) = 20. ! accretion radius
 vxyz_ptmass(1:3,1)    = -40.
 fxyz_ptmass(1:3,1)    = 40.
 hfact = hfact_default

 if (itest==1) then
    !--setup 2 SPH particles at (5,5,5)
    if (id==master) then
       call set_particle_type(1,igas)
       call set_particle_type(2,igas)
       npartoftype(igas) = 2
       npart        = 2
       xyzh(1:3,1:2)  = 5.
       xyzh(4,1:2)    = 0.01
       vxyzu(1:3,1) = [40.,40.,-10.]
       vxyzu(1:3,2) = [120.,120.,-30.]
       fxyzu(1:3,1:2) = 20.
       massoftype(1)  = 5.
    else
       npartoftype(igas) = 0
       npart        = 0
    endif
 else
    ! eat a large portion of a disc
    np_disc = 1000
    call set_disc(id,master,nparttot=np_disc,npart=npart,rmin=1.,rmax=2.*xyzmh_ptmass(ihacc,1),p_index=1.0,q_index=0.75,&
                  HoverR=0.1,disc_mass=0.5*xyzmh_ptmass(4,1),star_mass=xyzmh_ptmass(4,1),gamma=1.,&
                  particle_mass=massoftype(igas),hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,&
                  polyk=polyk,verbose=.false.)
    npartoftype(igas) = npart
 endif

 !--perform a test of the accretion of the SPH particle by the point mass
 nfailed(:)  = 0
 !--check energies before accretion event
 t=0.
 call compute_energies(t)
 etotin   = etot
 totmomin = totmom
 angmomin = angtot
 ibin_wakei = 0

 dptmass(:,1:nptmass) = 0.
 !$omp parallel default(shared) private(i) firstprivate(dptmass_thread)
 dptmass_thread(:,1:nptmass) = 0.
 !$omp do
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call ptmass_accrete(1,nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                           vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),fxyzu(1,i),fxyzu(2,i),fxyzu(3,i), &
                           igas,massoftype(igas),xyzmh_ptmass,vxyz_ptmass, &
                           accreted,dptmass_thread,t,1.0,ibin_wakei,ibin_wakei)
    endif
 enddo
 !$omp enddo
 !$omp critical(dptmassadd)
 dptmass(:,1:nptmass) = dptmass(:,1:nptmass) + dptmass_thread(:,1:nptmass)
 !$omp end critical(dptmassadd)
 !$omp end parallel

 call reduce_in_place_mpi('+',dptmass(:,1:nptmass))

 if (id==master) call update_ptmass(dptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass)

 call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
 call bcast_mpi(vxyz_ptmass(:,1:nptmass))
 call bcast_mpi(fxyz_ptmass(:,1:nptmass))

 if (itest==1) then
    call bcast_mpi(accreted)
    call bcast_mpi(xyzh(4,1:2))
    call checkval(accreted,.true.,nfailed(1),'accretion flag')
    !--check that h has been changed to indicate particle has been accreted
    call checkval(isdead_or_accreted(xyzh(4,1)),.true.,nfailed(2),'isdead_or_accreted flag(1)')
    call checkval(isdead_or_accreted(xyzh(4,2)),.true.,nfailed(2),'isdead_or_accreted flag(2)')
    call checkval(xyzmh_ptmass(1,1),3.,tiny(0.),nfailed(3),'x(ptmass) after accretion')
    call checkval(xyzmh_ptmass(2,1),3.,tiny(0.),nfailed(4),'y(ptmass) after accretion')
    call checkval(xyzmh_ptmass(3,1),3.,tiny(0.),nfailed(5),'z(ptmass) after accretion')
    call checkval(vxyz_ptmass(1,1),20.,tiny(0.),nfailed(6),'vx(ptmass) after accretion')
    call checkval(vxyz_ptmass(2,1),20.,tiny(0.),nfailed(7),'vy(ptmass) after accretion')
    call checkval(vxyz_ptmass(3,1),-30.,tiny(0.),nfailed(8),'vz(ptmass) after accretion')
    call checkval(fxyz_ptmass(1,1),30.,tiny(0.),nfailed(9), 'fx(ptmass) after accretion')
    call checkval(fxyz_ptmass(2,1),30.,tiny(0.),nfailed(10),'fy(ptmass) after accretion')
    call checkval(fxyz_ptmass(3,1),30.,tiny(0.),nfailed(11),'fz(ptmass) after accretion')

    call update_test_scores(ntests,nfailed(1:2),npass)
    call update_test_scores(ntests,nfailed(3:5),npass)
    call update_test_scores(ntests,nfailed(6:8),npass)
    call update_test_scores(ntests,nfailed(9:11),npass)
 endif

 !--compute conserved quantities after accretion event
 nfailed(:) = 0
 call compute_energies(t)
 call checkval(angtot,angmomin,1.e-14,nfailed(3),'angular momentum')
 call checkval(totmom,totmomin,2.*epsilon(0.),nfailed(2),'linear momentum')
 !call checkval(etot,etotin,1.e-6,'total energy',nfailed(1))
 call update_test_scores(ntests,nfailed(3:3),npass)
 call update_test_scores(ntests,nfailed(2:2),npass)

end subroutine test_accretion

!-----------------------------------------------------------------------
!+
!  Test sink particle creation
!+
!-----------------------------------------------------------------------
subroutine test_createsink(ntests,npass)
 use dim,        only:gravity,maxp,maxphase
 use boundary,   only:set_boundary
 use deriv,      only:get_derivs_global
 use eos,        only:ieos,polyk
 use kdtree,     only:tree_accuracy
 use units,      only:set_units
 use io,         only:id,master,iverbose
 use part,       only:init_part,npart,npartoftype,igas,xyzh,massoftype,hfact,rhoh,&
                      iphase,isetphase,fext,divcurlv,vxyzu,fxyzu,poten, &
                      nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,ndptmass, &
                      dptmass,fxyz_ptmass_sinksink,sf_ptmass
 use ptmass,     only:ptmass_accrete,update_ptmass,icreate_sinks,&
                      ptmass_create,finish_ptmass,ipart_rhomax,h_acc,rho_crit,rho_crit_cgs, &
                      ptmass_create_stars,tmax_acc,tseeds,ipart_createseeds,ipart_createstars,&
                      ptmass_create_seeds,get_accel_sink_sink,n_max
 use energies,   only:compute_energies,angtot,etot,totmom
 use mpiutils,   only:bcast_mpi,reduce_in_place_mpi,reduceloc_mpi,reduceall_mpi
 use spherical,  only:set_sphere
 use stretchmap, only:rho_func
 use units,      only:set_units
 use physcon,    only:solarm,pc
 use setup_params, only:npart_total
 integer, intent(inout) :: ntests,npass
 integer :: i,j,itest,itestp,nfailed(6),imin(1)
 integer :: id_rhomax,ipart_rhomax_global
 real :: psep,totmass,r2min,r2,t,coremass,starsmass
 real :: etotin,angmomin,totmomin,rhomax,rhomax_test
 real :: ke,pe,pei,d2,d1,rmax,ri(3)
 logical :: rtest,stest
 procedure(rho_func), pointer :: density_func

 call set_units(mass=1.d0,dist=1.d0,G=1.d0)
 density_func => gaussianr
 t        = 0.
 iverbose = 1
 rho_crit = rho_crit_cgs
 ieos     = 1
 polyk    = 0.

 do itest=1,3
    select case(itest)
    case(3)
       if (id==master) write(*,"(/,a)") '--> testing sink particle creation (cores and stars prescription)'
    case(2)
       if (id==master) write(*,"(/,a)") '--> testing sink particle creation (sin)'
    case default
       if (id==master) write(*,"(/,a)") '--> testing sink particle creation (uniform density)'
    end select
    !
    ! initialise arrays to zero
    !
    call set_units(mass=solarm,dist=pc,G=1.d0)
    call init_part()
    vxyzu(:,:) = 0.
    fxyzu(:,:) = 0.
    fext(:,:)  = 0.

    !
    ! set a boundary that is larger than the sphere size, so test still works with periodic boundaries
    !
    call set_boundary(-1.,1.,-1.,1.,-1.,1.)
    !
    ! set up gas particles in a uniform sphere with radius R=0.2
    !
    psep = 0.05  ! required as a variable since this may change under conditions not requested here
    npart_total = 0
    if (id == master) then
       if (itest==2) then
          ! use random so particle with maximum density is unique
          call set_sphere('cubic',id,master,0.,0.2,psep,hfact,npartoftype(igas),xyzh,nptot=npart_total,rhofunc=density_func)
       else
          call set_sphere('cubic',id,master,0.,0.2,psep,hfact,npartoftype(igas),xyzh,nptot=npart_total)
       endif
    else
       npartoftype(igas) = 0
    endif
    totmass = 1.0
    massoftype(igas) = totmass/real(reduceall_mpi('+',npart_total))  ! reduceall because only setup particles on master thread
    npart = npartoftype(igas)

    if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
    !
    ! set up tree for neighbour finding
    ! and make sure that gravitational potential energy has been computed
    !
    tree_accuracy = 0.
    if (itest==3) then
       icreate_sinks = 2
       sf_ptmass = 0.
       tmax_acc = 0.
       tseeds = 0.
       ipart_createseeds = 1
       ipart_createstars = 1
    else
       icreate_sinks = 1
    endif

    call get_derivs_global()

    !
    ! calculate itest after calling derivs because particles will
    ! rebalance across tasks
    !
    r2min = huge(r2min)
    itestp = npart
    do i=1,npart
       r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
       if (r2 < r2min) then
          itestp = i
          r2min = r2
       endif
    enddo

    !
    ! check that particle being tested is at the maximum density
    !
    if (itest==2 .and. gravity) then
       imin = minloc(xyzh(4,1:npart))
       itestp = imin(1)
       rhomax_test = rhoh(xyzh(4,itestp),massoftype(igas))
       !
       ! only check on the thread that has rhomax
       !
       ipart_rhomax_global = ipart_rhomax
       call reduceloc_mpi('max',ipart_rhomax_global,id_rhomax)
       if (id == id_rhomax) then
          rhomax = rhoh(xyzh(4,ipart_rhomax),massoftype(igas))
          call checkval(rhomax,rhomax_test,epsilon(0.),nfailed(1),'rhomax',thread_id=id)
       else
          itestp = -1 ! set itest = -1 on other threads
          call checkval(ipart_rhomax,-1,0,nfailed(1),'ipart_rhomax',thread_id=id)
       endif
       call update_test_scores(ntests,nfailed(1:1),npass)
    endif
    !
    ! check energies before insertion of sink
    !
    call compute_energies(t)
    etotin   = etot
    totmomin = totmom
    angmomin = angtot
    !
    ! now create point mass by accreting these particles
    !
    h_acc = 0.15

    !
    ! if gravity is not enabled, then need to choose a particle to create ptmass from
    !
    if (.not. gravity) then
       ipart_rhomax_global = itestp
       call reduceloc_mpi('max',ipart_rhomax_global,id_rhomax)
    endif
    call ptmass_create(nptmass,npart,itestp,xyzh,vxyzu,fxyzu,fext,divcurlv,poten,&
                       massoftype,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,fxyz_ptmass_sinksink,sf_ptmass,dptmass,0.)
    if (itest==3) then
       coremass = 0.
       starsmass = 0.
       ke = 0.
       pe = 0.
       rmax = epsilon(rmax)
       coremass = xyzmh_ptmass(4,1)
       ri(3)    = xyzmh_ptmass(3,1)
       ri(2)    = xyzmh_ptmass(2,1)
       ri(1)    = xyzmh_ptmass(1,1)
       call ptmass_create_seeds(nptmass,ipart_createseeds,sf_ptmass,0.)
       call ptmass_create_stars(nptmass,ipart_createstars,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass, &
                                fxyz_ptmass_sinksink,sf_ptmass,0.)
       do i=1,nptmass
          pei = 0.
          do j=1,nptmass
             if (j/=i) then
                d2 = (xyzmh_ptmass(1,i)-xyzmh_ptmass(1,j))**2+&
                     (xyzmh_ptmass(2,i)-xyzmh_ptmass(2,j))**2+&
                     (xyzmh_ptmass(3,i)-xyzmh_ptmass(3,j))**2
                d1 = 1./sqrt(d2)
                pei = pei + xyzmh_ptmass(4,j)*d1
             endif
          enddo
          pe = pe + 0.5*pei*xyzmh_ptmass(4,i)
          starsmass = starsmass + xyzmh_ptmass(4,i)
          ke  = ke + 0.5*xyzmh_ptmass(4,i)*(vxyz_ptmass(1,i)**2 + vxyz_ptmass(2,i)**2 + vxyz_ptmass(3,i)**2)
          rmax = max(sqrt((xyzmh_ptmass(1,i)-ri(1))**2+(xyzmh_ptmass(2,i)-ri(2))**2+(xyzmh_ptmass(3,i)-ri(3))**2),rmax)
       enddo



    endif
    !
    ! check that creation succeeded
    !
    nfailed(:) = 0
    if (itest == 3) then
       rtest = rmax < h_acc
       stest = nptmass < n_max
       call checkval(stest,.true.,nfailed(1),'nptmass< nseeds max')
       call checkval(starsmass-coremass,0.,6e-17,nfailed(4),'Mass conservation')
       call checkval(ke/pe,0.5,5e-16,nfailed(5),'Virialised system')
       call checkval(rtest,.true.,nfailed(6),'rmax < h_acc')
    else
       call checkval(nptmass,1,0,nfailed(1),'nptmass=1')
    endif
    call update_test_scores(ntests,nfailed,npass)
    !
    ! check that linear and angular momentum and energy is conserved
    !
    nfailed(:) = 0
    call compute_energies(t)
    if (itest /= 3) call checkval(angtot,angmomin,1.e-10,nfailed(3),'angular momentum')
    call checkval(totmom,totmomin,epsilon(0.),nfailed(2),'linear momentum')
    !call checkval(etot,etotin,1.e-6,nfailed(1),'total energy')
    call update_test_scores(ntests,nfailed,npass)

    call finish_ptmass(nptmass)
 enddo

 ! reset options
 iverbose = 0
 icreate_sinks  = 0

end subroutine test_createsink

!-----------------------------------------------------------------------
!+
!  Test sink particle mergers
!+
!-----------------------------------------------------------------------
subroutine test_merger(ntests,npass)
 use dim,            only:periodic,use_sinktree
 use io,             only:id,master,iverbose
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass, &
                          npart,ihacc,itbirth,epot_sinksink,dsdt_ptmass,&
                          sf_ptmass,shortsinktree,fxyz_ptmass_tree
 use ptmass,         only:h_acc,h_soft_sinksink,get_accel_sink_sink, &
                          r_merge_uncond,r_merge_cond,r_merge_uncond2,&
                          r_merge_cond2,r_merge2,icreate_sinks,n_max
 use random,         only:ran2
 use step_lf_global, only:init_step,step
 use timestep,       only:dtmax
 use mpiutils,       only:bcast_mpi,reduce_in_place_mpi
 use energies,       only:compute_energies,angtot,totmom,mtot

 integer, intent(inout) :: ntests,npass
 integer, parameter :: max_to_test = 100
 logical, parameter :: print_sink_paths = .false. ! print sink paths in the merger test
 integer :: i,j,iseed,itest,nfailed(80),merge_ij(max_to_test),merge_n
 integer :: nsink0,nsinkf,nsteps
 logical :: merged,merged_expected
 real :: t,dt,dtext,dtnew,dtsinksink,r2,v2
 real :: angmom0,mtot0,mv0,dx(3),dv(3)
 real :: fxyz_sinksink(4,max_to_test)

 iseed           = -74205
 nfailed(:)      = 0
 iverbose        = 0
 nptmass         = 2
 npart           = 0
 h_acc           = 0.1
 h_soft_sinksink = h_acc
 r_merge_uncond  = 2.*h_acc    ! sinks will unconditionally merge if they touch
 r_merge_cond    = 4.*h_acc    ! sinks will merge if bound within this radius
 r_merge_uncond2 = r_merge_uncond**2
 r_merge_cond2   = r_merge_cond**2
 r_merge2        = max(r_merge_uncond2,r_merge_cond2)
 if (use_sinktree) then
    shortsinktree = 1
    fxyz_ptmass_tree = 0.
 endif
 do itest=1,10
    t                 = 0.
    xyzmh_ptmass(:,:) = 0.
    xyzmh_ptmass(4,:) = 1.
    xyzmh_ptmass(ihacc,:) = h_acc
    vxyz_ptmass(:,:)  = 0.
    icreate_sinks = 1
    select case(itest)
    case(1)
       if (id==master) write(*,"(/,a)") '--> testing fast flyby: no merger'
       ! fast flyby within r_merge_uncond < r < r_merge_cond
       xyzmh_ptmass(1,1) =  1.
       xyzmh_ptmass(2,1) =  1.5*h_acc
       vxyz_ptmass(1,1)  = -10.
       merged_expected   = .false.
    case(2)
       if (id==master) write(*,"(/,a)") '--> testing fast flyby: impact so merger'
       ! fast flyby within r < r_merge_uncond
       xyzmh_ptmass(1,1) =  1.
       xyzmh_ptmass(2,1) =  0.5*h_acc
       vxyz_ptmass(1,1)  = -10.
       merged_expected   = .true.
    case(3)
       if (id==master) write(*,"(/,a)") '--> testing slow flyby: capture and merger'
       ! slow flyby within r_merge_uncond < r < r_merge_cond
       xyzmh_ptmass(1,1) =  1.
       xyzmh_ptmass(2,1) =  1.5*h_acc
       vxyz_ptmass(1,1)  = -1.
       merged_expected   = .true.
    case(4)
       if (id==master) write(*,"(/,a)") '--> testing slow flyby: impact and merger'
       ! slow flyby within r < r_merge_cond
       xyzmh_ptmass(1,1) =  1.
       xyzmh_ptmass(2,1) =  0.5*h_acc
       vxyz_ptmass(1,1)  = -1.
       merged_expected   = .true.
    case(5)
       if (id==master) write(*,"(/,a)") '--> testing flyby: slingshot & no merger'
       ! flyby within r_merge_uncond < r < r_merge_cond
       xyzmh_ptmass(1,1) =  1.
       xyzmh_ptmass(2,1) =  1.5*h_acc
       vxyz_ptmass(1,1)  = -5.
       merged_expected   = .false.
    case(6)
       if (id==master) write(*,"(/,a)") '--> testing orbit: stable & no merger'
       ! stable orbit within r >  r_merge_cond
       xyzmh_ptmass(1,1) = 2.5*h_acc
       vxyz_ptmass(2,1)  = sqrt(0.25*xyzmh_ptmass(4,1)/xyzmh_ptmass(1,1))
       merged_expected   = .false.
    case(7)
       if (id==master) write(*,"(/,a)") '--> testing orbit: decaying & merger'
       ! decaying orbit within r >  r_merge_cond
       xyzmh_ptmass(1,1) = 2.5*h_acc
       vxyz_ptmass(2,1)  = 0.9*sqrt(0.25*xyzmh_ptmass(4,1)/xyzmh_ptmass(1,1))
       merged_expected   = .true.
    case(8)
       if (id==master) write(*,"(/,a)") '--> testing multiple sink interations'
       nptmass = max_to_test
       do i = 1,nptmass
          xyzmh_ptmass(1:3,i) = ( (/ran2(iseed),ran2(iseed),ran2(iseed)/) - 0.5) * 2.  ! in range (-1,1)
          vxyz_ptmass(1:3,i)  = ( (/ran2(iseed),ran2(iseed),ran2(iseed)/) - 0.5) * 6.  ! in range (-3,3)
       enddo
       merged_expected   = .true. ! this logical does not have meaning here
    case(9)
       if (id==master) write(*,"(/,a)") '--> testing release during merging with icreate_sinks == 2'
       nptmass = 2
       xyzmh_ptmass(1,1) =  1.
       xyzmh_ptmass(2,1) =  0.5*h_acc
       vxyz_ptmass(1,1)  = -10.
       xyzmh_ptmass(1,itbirth) = 0.2
       xyzmh_ptmass(2,itbirth) = 0.4
       n_max = 5
       icreate_sinks     = 2
       sf_ptmass(1,:)    = 1
       sf_ptmass(2,1)    = 4
       sf_ptmass(2,2)    = 3
       merged_expected   = .true.
    case(10)
       if (id==master) write(*,"(/,a)") '--> testing merging with icreate_sinks == 2 (one sink is only gas)'
       nptmass = 2
       xyzmh_ptmass(1,1) =  1.
       xyzmh_ptmass(2,1) =  0.5*h_acc
       vxyz_ptmass(1,1)  = -10.
       xyzmh_ptmass(1,itbirth) = 0.01
       xyzmh_ptmass(2,itbirth) = 0.4
       n_max = 5
       icreate_sinks     = 2
       sf_ptmass(1,:)    = 1
       sf_ptmass(2,1)    = 0
       sf_ptmass(2,2)    = 3
       merged_expected   = .true.


    end select
    if (itest /= 8) then
       xyzmh_ptmass(1:3,2) = -xyzmh_ptmass(1:3,1)
       vxyz_ptmass(1:3,2)  = -vxyz_ptmass(1:3,1)
    endif
    !
    ! get initial values
    !
    call compute_energies(0.)
    nsink0 = nptmass; angmom0 = angtot; mv0 = totmom; mtot0 = mtot
    !
    ! initialise forces
    !
    if (id==master) then
       call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_sinksink,epot_sinksink,&
                                dtsinksink,0,0.,merge_ij,merge_n,dsdt_ptmass)
    endif
    fxyz_ptmass(:,:) = 0.
    call bcast_mpi(epot_sinksink)
    call bcast_mpi(dtsinksink)

    if (id==master) fxyz_ptmass(:,1:nptmass) = fxyz_ptmass(:,1:nptmass) + fxyz_sinksink(:,1:nptmass)
    call reduce_in_place_mpi('+',fxyz_ptmass)
    !
    ! integrate
    !
    nsteps = 1000
    dt     = huge(dt)
    do i = 1,nptmass-1
       do j = i+1,nptmass
          dx = xyzmh_ptmass(1:3,i)-xyzmh_ptmass(1:3,j)
          r2   = dot_product(dx,dx)
          dv = vxyz_ptmass(1:3,i)-vxyz_ptmass(1:3,j)
          v2   = dot_product(dv,dv)+epsilon(v2)
          dt   = min(dt,sqrt(r2/v2))
       enddo
    enddo
    dt     = 2.*dt/nsteps
    dtmax  = dt*nsteps
    t      = 0.
    call init_step(npart,t,dtmax)
    if (print_sink_paths) then
       write(333,*) itest,0,xyzmh_ptmass(1:4,1),xyzmh_ptmass(1:4,2)
       if (itest==8) then
          do j = 1,nptmass
             if (xyzmh_ptmass(4,j) > 0.) write(334,*) t,j,xyzmh_ptmass(1:4,j)
          enddo
       endif
    endif
    do i=1,nsteps
       t = t + dt
       dtext = dt
       if (id==master .and. iverbose > 2) write(*,*) ' t = ',t,' dt = ',dt
       call step(npart,npart,t,dt,dtext,dtnew)
       if (print_sink_paths) then
          write(333,*) itest,i,xyzmh_ptmass(1:4,1),xyzmh_ptmass(1:4,2)
          if (itest==8) then
             do j = 1,nptmass
                if (xyzmh_ptmass(4,j) > 0.) write(334,*) t,j,xyzmh_ptmass(1:4,j)
             enddo
          endif
       endif
    enddo
    !
    ! check results
    !
    call compute_energies(t)
    nsinkF = count(xyzmh_ptmass(4,1:nptmass) > 0.) ! only count non-merged sinks
    !
    if (xyzmh_ptmass(4,2) < 0.) then
       merged = .true.
    else
       merged = .false.
    endif
    if (itest==8) then
       call checkval(nsinkF,41,0,nfailed(itest),'final number of sinks')
    else
       call checkval(merged,merged_expected,nfailed(itest),'merger')
       if (merged_expected .and. itest/=9) then
          call checkval(xyzmh_ptmass(1,1),0.,epsilon(0.),nfailed(2*itest),'final x-position')
          call checkval(xyzmh_ptmass(2,1),0.,epsilon(0.),nfailed(3*itest),'final y-position')
          v2 = dot_product(vxyz_ptmass(1:2,1),vxyz_ptmass(1:2,1))
          call checkval(sqrt(v2),0.,epsilon(0.),nfailed(4*itest),'final velocity')
       endif
    endif
    call checkval(totmom,    mv0,1.e-13,nfailed(5*itest),'conservation of linear momentum')
    if ( itest/=9) then
       call checkval(angtot,angmom0,1.e-13,nfailed(6*itest),'conservation of angular momentum')
    endif
    call checkval(mtot,    mtot0,1.e-13,nfailed(7*itest),'conservation of mass')
    if (itest==9) then
       call checkval(sf_ptmass(2,1)+nsinkF,8,0,nfailed(8*itest),'conservation of star seeds')
    elseif (itest==10) then
       call checkval(sf_ptmass(2,2)+nsinkF,4,0,nfailed(8*itest),'conservation of star seeds')
    endif
 enddo
 call update_test_scores(ntests,nfailed(1:80),npass)

 ! reset options
 r_merge_uncond = 0.
 r_merge_cond   = 0.
 icreate_sinks  = 0

end subroutine test_merger

!-----------------------------------------------------------------------
!+
!  Test HII region expansion around sink particles
!+
!-----------------------------------------------------------------------
subroutine test_HIIregion(ntests,npass)
 use dim,            only:maxp,maxphase,maxvxyzu
 use io,             only:id,master,iverbose,iprint
 use eos_HIIR,       only:polykion,init_eos_HIIR
 use eos,            only:gmw,ieos,polyk,gamma
 use deriv,          only:get_derivs_global
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass, &
                            npart,ihacc,irstrom,xyzh,vxyzu,hfact,igas, &
                            npartoftype,fxyzu,massoftype,isionised,init_part,&
                            iphase,isetphase,irateion,irstrom
 use ptmass,         only:h_acc
 use step_lf_global, only:init_step,step
 use spherical,      only:set_sphere
 use units,          only:set_units,utime,unit_velocity,udist,umass
 use physcon,        only:pc,solarm,years,pi,kboltz,mass_proton_cgs
 use kernel,         only:hfact_default
 use kdtree,         only:tree_accuracy
 use testutils,      only:checkval,update_test_scores
 use HIIRegion,      only:initialize_H2R,update_ionrate,HII_feedback,iH2R,nHIIsources,ar,mH
 use setup_params,   only:npart_total
 integer, intent(inout) :: ntests,npass
 integer        :: np,i,nfailed(1)
 real           :: totmass,psep
 real           :: Rstrom,ci,k,rho0
 real           :: totvol,nx,rmin,rmax,temp
 if (id==master) write(iprint,"(/,a)") '--> testing HII region expansion around massive stars...'

 call set_units(dist=pc,mass=solarm,G=1.d0)
 call init_eos_HIIR()
 iverbose = 0
 !
 ! initialise arrays to zero
 !
 call init_part()
 gmw = 1.0

 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:)  = 0.

 h_acc = 0.002

 xyzmh_ptmass(4,1) = -1.
 xyzmh_ptmass(irateion,1) = 49. ! rate_ion [s^-1]
 nptmass = 1
 nHIIsources = 1

 hfact = 1.2
 gamma = 1.
 rmin  = 0.
 rmax  = 2.91*pc/udist
 ieos  = 21
 tree_accuracy = 0.5
 temp = 1000.
!
!--setup particles
!
 np       = 1000000
 totvol   = 4./3.*pi*rmax**3
 nx       = int(np**(1./3.))
 psep     = totvol**(1./3.)/real(nx)
 npart    = 0
 npart_total = 0
 ! only set up particles on master, otherwise we will end up with n duplicates
 if (id==master) then
    call set_sphere('cubic',id,master,rmin,rmax,psep,hfact,npart,xyzh,nptot=npart_total,np_requested=np)
 endif
 np       = npart


!
!--set particle properties
!
 totmass        = 8.e3*solarm/umass
 npartoftype(:) = 0
 npartoftype(igas) = npart
 massoftype(:)  = 0.0
 massoftype(igas)  = totmass/npartoftype(igas)
 if (maxphase==maxp) then
    do i=1,npart
       iphase(i) = isetphase(igas,iactive=.true.)
    enddo
 endif


 iH2R = 1
 if (id==master) then
    call initialize_H2R
    !call HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyz_ptmass,isionised)
 endif

 rho0 = totmass/totvol

 Rstrom = 10**((1./3)*(log10(((3*mH**2)/(4*pi*ar*rho0**2)))+xyzmh_ptmass(irateion,1)+log10(utime)))
 xyzmh_ptmass(irstrom,1) = -1.
 ci   = sqrt(polykion)
 k = 0.005

 polyk = (kboltz*temp)/(gmw*mass_proton_cgs)*((unit_velocity)**2)
 vxyzu(:,:) = 0.
 fxyzu(:,:) = 0.
 if (maxvxyzu >= 4) then
    vxyzu(4,:) = polyk
    ieos = 22
 endif

 call get_derivs_global()

 call HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyz_ptmass,isionised)

 call checkval(xyzmh_ptmass(irstrom,1),Rstrom,1.e-2,nfailed(1),'Initial strömgren radius')

 call update_test_scores(ntests,nfailed,npass)

end subroutine test_HIIregion

!-----------------------------------------------------------------------
!+
!  Test SDAR integration method on a stable triple system
!+
!-----------------------------------------------------------------------
subroutine test_SDAR(ntests,npass)
 use dim,        only:periodic,gravity,ind_timesteps
 use io,         only:id,master,iverbose
 use physcon,    only:pi,deg_to_rad
 use ptmass,     only:get_accel_sink_sink,h_soft_sinksink, &
                        get_accel_sink_gas,f_acc,use_fourthorder,use_regnbody
 use part,       only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass,fext,&
                        npart,npartoftype,massoftype,xyzh,vxyzu,&
                        igas,epot_sinksink,init_part,iJ2,ispinx,ispiny,ispinz,iReff,istar
 use part,       only:group_info,bin_info,n_group,n_ingroup,n_sing,nmatrix
 use energies,   only:angtot,etot,totmom,compute_energies
 use timestep,   only:dtmax,C_force,tolv
 use kdtree,     only:tree_accuracy
 use eos,        only:ieos
 use setbinary,  only:set_binary
 use units,      only:set_units
 use mpiutils,   only:bcast_mpi,reduce_in_place_mpi
 use step_lf_global, only:init_step,step
 use testutils,      only:checkvalf,checkvalbuf,checkvalbuf_end
 use checksetup,     only:check_setup
 use deriv,          only:get_derivs_global
 use timing,         only:getused,printused
 use options,        only:ipdv_heating,ishock_heating
 use subgroup,       only:group_identify,r_neigh
 use centreofmass,   only:reset_centreofmass
 integer,          intent(inout) :: ntests,npass
 integer :: i,ierr,nfailed(4),nerr,nwarn
 integer :: merge_ij(3),merge_n
 real :: m1,m2,a,ecc,incl,hacc1,hacc2,dt,dtext,t,dtnew,tolen,tolmom,tolang,tolecc
 real :: angmomin,etotin,totmomin,dum,dum2,omega,errmax,dtsinksink,tmax,eccfin,decc
 real :: fxyz_sinksink(4,3),dsdt_sinksink(3,3) ! we only use 3 sink particles in the tests here
 real :: xsec(3),vsec(3)
 real(kind=4) :: t1
 if (id==master) write(*,"(/,a)") '--> testing SDAR module : Kozai-Lidov effect'
 !
 !--no gas particles
 !
 call init_part()
 iverbose = 0
 tree_accuracy = 0.
 h_soft_sinksink = 0.
 ipdv_heating = 0
 ishock_heating = 0
 use_regnbody = .true.
 r_neigh = 10.
 use_fourthorder = .true.

 tolv = 1e-2

 !
 !--setup triple system with Kozai-Lidov resonance
 !
 npart = 0
 npartoftype = 0
 nptmass = 0
 m1    = 2.0
 m2    = 1.0
 a     = 1.0000
 ecc = 0.990000
 incl = 0.10/deg_to_rad
 hacc1  = 1e-4
 hacc2  = 1e-4
 C_force = 0.25
 omega = sqrt((m1+m2)/a**3)
 t = 0.
 call set_units(mass=1.d0,dist=1.d0,G=1.d0)
 call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                 posang_ascnode=0.0,arg_peri=0.0,incl=incl,mean_anomaly=179.999999,verbose=.false.)


 xsec(1:3) = xyzmh_ptmass(1:3,2)
 vsec(1:3) = vxyz_ptmass(1:3,2)
 m1 = 0.90
 m2 = 0.10
 a  = 0.00099431556644
 ecc = 0.90000
 incl = 1.5/deg_to_rad

 nptmass = nptmass - 1
 xyzmh_ptmass(:,2) = 0.
 vxyz_ptmass(:,2)  = 0.

 call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                 posang_ascnode=0.0,arg_peri=0.0,incl=incl,mean_anomaly=179.999999,verbose=.false.)

 xyzmh_ptmass(1:3,2) =  xyzmh_ptmass(1:3,2) + xsec(1:3)
 vxyz_ptmass(1:3,2)  =  vxyz_ptmass(1:3,2) + vsec(1:3)
 xyzmh_ptmass(1:3,3) =  xyzmh_ptmass(1:3,3) + xsec(1:3)
 vxyz_ptmass(1:3,3)  =  vxyz_ptmass(1:3,3) + vsec(1:3)


 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)



 if (ierr /= 0) nerr = nerr + 1

 !
 ! check that no errors occurred when setting up initial conditions
 !
 nfailed = 0
 call check_setup(nerr,nwarn)
 call checkval(nerr,0,0,nfailed(1),'no errors during setup')
 call update_test_scores(ntests,nfailed,npass)

 tolv = 1.e-2
 iverbose = 0
 ieos = 1
 !
 ! initialise forces
 !
 if (id==master) then
    call group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,vxyz_ptmass,&
                        group_info,bin_info,nmatrix)
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_sinksink,epot_sinksink,&
                                  dtsinksink,0,0.,merge_ij,merge_n,dsdt_sinksink,&
                                  group_info=group_info,bin_info=bin_info)
 endif
 fxyz_ptmass(:,1:nptmass) = 0.
 dsdt_ptmass(:,1:nptmass) = 0.
 call bcast_mpi(epot_sinksink)
 call bcast_mpi(dtsinksink)

 fext(:,:) = 0.
 do i=1,npart
    call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                  fext(1,i),fext(2,i),fext(3,i),dum,massoftype(igas),fxyz_ptmass,dsdt_ptmass,dum,dum2)
 enddo
 if (id==master) then
    fxyz_ptmass(:,1:nptmass) = fxyz_ptmass(:,1:nptmass) + fxyz_sinksink(:,1:nptmass)
    dsdt_ptmass(:,1:nptmass) = dsdt_ptmass(:,1:nptmass) + dsdt_sinksink(:,1:nptmass)
 endif
 call reduce_in_place_mpi('+',fxyz_ptmass(:,1:nptmass))
 call reduce_in_place_mpi('+',dsdt_ptmass(:,1:nptmass))


 dt = 0.01

 dtmax = dt  ! required prior to derivs call, as used to set ibin


 !
 !--evolve this for a number of orbits
 !
 call compute_energies(t)
 etotin   = etot
 totmomin = totmom
 angmomin = angtot
 call bcast_mpi(etotin)
 call bcast_mpi(totmomin)
 call bcast_mpi(angmomin)
 ecc      = bin_info(2,2)
 call bcast_mpi(ecc)
 decc     = 0.09618

 tmax = 7.*3.63 ! 7 out binary periods
 t    = 0.
 errmax = 0.
 f_acc = 1.
 !
 !--integration loop
 !
 if (id==master) call getused(t1)
 call init_step(npart,t,dtmax)
 do while (t < tmax)
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    call compute_energies(t)
    errmax = max(errmax,abs(etot - etotin))
    t = t + dt
 enddo

 call compute_energies(t)

 if (id==master) call printused(t1)
 nfailed(:) = 0
 eccfin = 0.99617740539553523
 tolecc = 3e-5
 tolmom = 2.e-11
 tolang = 3.e-11
 tolen  = 8.e-6
 !
 !--check energy conservation
 !
 call checkval(angtot,angmomin,tolang,nfailed(1),'angular momentum')
 call checkval(totmom,totmomin,tolmom,nfailed(2),'linear momentum')
 call checkval(etotin+errmax,etotin,tolen,nfailed(3),'total energy')
 call checkval(eccfin-ecc,decc,tolecc,nfailed(4),'eccentricity')
 do i=1,4
    call update_test_scores(ntests,nfailed(i:i),npass)
 enddo

 use_regnbody = .false.
 call init_part()

end subroutine test_SDAR

!-----------------------------------------------------------------------
!+
!  Test sink particle surface force, simply that the acceleration
!  is the gradient of the potential
!+
!-----------------------------------------------------------------------
subroutine test_sink_potential(ntests,npass)
 use io,         only:id,master
 use testutils,  only:checkval,update_test_scores
 use ptmass,     only:get_accel_sink_gas,isink_potential
 use part,       only:npart,npartoftype,nptmass,xyzmh_ptmass,ihacc,iReff
 use units,      only:set_units
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(1)
 real :: phi1,phi,eps,x0(3)
 real :: dphidx,hi,xi,yi,zi,dumxi,dumyi,dumzi,fxi,fyi,fzi,rp

 if (id==master) write(*,"(/,a)") '--> testing sink particle surface force'
 nptmass = 1
 npart = 0
 npartoftype = 0
 hi = 0.
 x0 = [100.,100.,100.]
 rp = 2.
 isink_potential = 1
 ! place a single point mass at a random location
 xyzmh_ptmass(:,:)     = 0.
 xyzmh_ptmass(1:3,1)   = x0
 xyzmh_ptmass(4,1)     = 3.14159
 xyzmh_ptmass(ihacc,1) = 0.
 xyzmh_ptmass(iReff,1) = rp  ! surface radius = 2

 call set_units(mass=1.d0,dist=1.d0,G=1.d0)

 ! evaluate sink-gas acceleration at some position
 xi = x0(1) + 1.00001*rp
 yi = x0(2) + 1.*rp
 zi = x0(3) + 1.*rp
 fxi = 0.; fyi = 0.; fzi = 0.; phi = 0.
 call get_accel_sink_gas(nptmass,xi,yi,zi,hi,xyzmh_ptmass,fxi,fyi,fzi,phi)
 ! evaluate sink-gas acceleration at some position + epsilon
 eps = 1.e-6
 dumxi = 0.; dumyi = 0.; dumzi = 0.; phi1 = 0.
 call get_accel_sink_gas(nptmass,xi+eps,yi,zi,hi,xyzmh_ptmass,dumxi,dumyi,dumzi,phi1)
 ! get the derivative of phi and check it equals the acceleration
 dphidx = -(phi1 - phi)/eps

 call checkval(dphidx,fxi,3.3e-8,nfailed(1),'dphi/dx = acceleration')
 call update_test_scores(ntests,nfailed(1:1),npass)

 ! reset options
 isink_potential = 0

end subroutine test_sink_potential

!-----------------------------------------------------------------------
!+
!  Helper function used in sink particle creation test
!+
!-----------------------------------------------------------------------
real function gaussianr(r)
 real, intent(in) :: r

 gaussianr = exp(-(r/0.05)**2) !1./(r**2 + 0.0001**2)

end function gaussianr

end module testptmass
