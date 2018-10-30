!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testptmass
!
!  DESCRIPTION:
!   Unit tests of the ptmass/sink particles module
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, checksetup, deriv, dim, energies, eos,
!    fileutils, io, kdtree, kernel, mpiutils, options, part, physcon,
!    ptmass, setbinary, setdisc, spherical, step_lf_global, testutils,
!    timestep, units
!+
!--------------------------------------------------------------------------
module testptmass
 implicit none
 public :: test_ptmass

 private

contains

subroutine test_ptmass(ntests,npass)
 use dim,      only:maxp,mhd,periodic,gravity,maxptmass
 use io,       only:id,master,iverbose,iskfile
 use part,     only:npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu,fext,&
                    xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,nptmass,epot_sinksink,&
                    ihacc,isdead_or_accreted,igas,divcurlv,iphase,isetphase,maxphase,&
                    Bevol,dBevol,dustfrac,ddustevol,temperature,divcurlB,fxyzu,set_particle_type,&
                    ispinx,ispiny,ispinz,dustprop,ddustprop,poten,rhoh
 use part,            only:pxyzu,dens,metrics
 use eos,             only:gamma,polyk
 use timestep,        only:dtmax,C_force,tolv
 use testutils,       only:checkval,checkvalf
 use setbinary,       only:set_binary
 use step_lf_global,  only:step,init_step
 use energies,        only:compute_energies,etot,totmom,epot,angtot !,accretedmass
 use ptmass,          only:get_accel_sink_sink,ptmass_accrete,h_soft_sinksink, &
                           ptmass_create,h_acc,get_accel_sink_gas,f_acc,finish_ptmass, &
                           ipart_rhomax,icreate_sinks, &
                           idxmsi,idymsi,idzmsi,idmsi,idspinxsi,idspinysi,idspinzsi, &
                           idvxmsi,idvymsi,idvzmsi,idfxmsi,idfymsi,idfzmsi, &
                           ndptmass,update_ptmass
 use physcon,         only:pi
 use setdisc,         only:set_disc
 use spherical,       only:set_sphere
 use boundary,        only:set_boundary
 use deriv,           only:derivs
 use kdtree,          only:tree_accuracy
 !use readwrite_dumps, only:write_fulldump,write_smalldump
 use fileutils,       only:getnextfilename
 use checksetup,      only:check_setup
 use options,         only:ieos,iexternalforce
 use units,           only:set_units
 use kernel,          only:kernel_softening
#ifdef IND_TIMESTEPS
 use part,            only:ibin
#endif
 use mpiutils,        only:bcast_mpi,reduce_in_place_mpi,reduceloc_mpi
 integer, intent(inout) :: ntests,npass
 integer                :: i,nsteps,nbinary_tests,itest,nerr,nwarn,itestp
 integer                :: nparttot
 logical                :: test_binary,test_accretion,test_createsink, test_softening
 logical                :: accreted
 real                   :: massr,m1,a,ecc,hacc1,hacc2,dt,dtext,t,dtnew,dr
 real                   :: etotin,totmomin,dtsinksink,omega,mred,errmax,angmomin
 real                   :: r2,r2min,dtext_dum,xcofm(3),totmass,dum,dum2,psep,tolen
 real                   :: xyzm_ptmass_old(4,1), vxyz_ptmass_old(3,1)
 real                   :: q,phisoft,fsoft,m2,mu,v_c1,v_c2,r1,omega1,omega2
 real                   :: dptmass(ndptmass,maxptmass)
 real                   :: dptmass_thread(ndptmass,maxptmass)
 real                   :: fxyz_sinksink(4,maxptmass),rhomax_test,rhomax
 integer                :: norbits,itmp,ierr
 integer                :: nfailed(11),imin(1)
 integer                :: id_rhomax,ipart_rhomax_global
 integer(kind=1)        :: ibin_wakei
 character(len=20)      :: dumpfile,filename

 if (id==master) write(*,"(/,a,/)") '--> TESTING PTMASS MODULE'

 test_binary = .true.
 test_accretion = .true.
 test_createsink = .true.
 test_softening = .true.
 nbinary_tests = 3
 !
 !--general settings
 !
 polyk = 0.
 gamma = 1.
 iexternalforce = 0
!
!  Test 1: orbit of a single binary
!  Test 2: with gas disc around it
!
 testbinary: if (test_binary) then
    !
    !--no gas particles
    !
    npart = 0
    npartoftype(:) = 0
    massoftype = 0.

    xyzh(:,:)  = 0.
    vxyzu(:,:) = 0.
    fxyzu(:,:) = 0.
    fext(:,:)  = 0.
    Bevol(:,:) = 0.

#ifdef IND_TIMESTEPS
    ibin(:) = 0_1
#endif
    iverbose = 0
    tree_accuracy = 0.

    binary_tests: do itest = 1,nbinary_tests
       select case(itest)
       case(2,3)
          if (periodic) then
             if (id==master) write(*,"(/,a)") '--> skipping circumbinary disc test (-DPERIODIC is set)'
             cycle binary_tests
          else
             if (itest==3) then
                if (id==master) write(*,"(/,a)") '--> testing integration of disc around eccentric binary'
             else
                if (id==master) write(*,"(/,a)") '--> testing integration of circumbinary disc'
             endif
          endif
       case default
          if (id==master) write(*,"(/,a)") '--> testing integration of binary orbit'
       end select
       !
       !--setup sink-sink binary (no gas particles)
       !
!       time = 0.
       nptmass = 0
       m1    = 1.
       massr = 1.
       a     = 1.
       if (itest==3) then
          ecc = 0.5
       else
          ecc   = 0.
       endif
       hacc1  = 0.35
       hacc2  = 0.35
       C_force = 0.25
       t = 0.
       call set_units(mass=1.d0,dist=1.d0,G=1.d0)
       call set_binary(m1,massr,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,verbose=.false.)
       if (itest==2 .or. itest==3) then
          !  add a circumbinary gas disc around it
          nparttot = 1000
          call set_disc(id,master,nparttot=nparttot,npart=npart,rmin=1.5*a,rmax=15.*a,p_index=1.5,q_index=0.75,&
                        HoverR=0.1,disc_mass=0.01*m1,star_mass=m1+massr*m1,gamma=gamma,&
                        particle_mass=massoftype(igas),hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,&
                        polyk=polyk,verbose=.false.)
          npartoftype(1) = npart

          !
          ! check that no errors occurred when setting up disc
          !
          ntests = ntests + 1
          nfailed = 0
          call check_setup(nerr,nwarn)
          call checkval(nerr,0,0,nfailed(1),'no errors during disc setup')
          !call checkval(nwarn,0,0,nfailed(2),'no warnings during disc setup')
          if (all(nfailed(1:2)==0)) npass = npass + 1

       endif

       tolv = 1.e3
       iverbose = 0
       ieos = 3
       !print*,'initial v = ',vxyz_ptmass(:,1)
       !
       ! initialise forces
       !
       if (id==master) then
          call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_sinksink,epot_sinksink,dtsinksink,0,0.)
       endif
       fxyz_ptmass(:,:) = 0.
       call bcast_mpi(epot_sinksink)
       call bcast_mpi(dtsinksink)

       fext(:,:) = 0.
       do i=1,npart
          call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                   fext(1,i),fext(2,i),fext(3,i),dum,massoftype(igas),fxyz_ptmass,dum,dum2)
       enddo
       if (id==master) fxyz_ptmass(:,:) = fxyz_ptmass(:,:) + fxyz_sinksink(:,:)
       call reduce_in_place_mpi('+',fxyz_ptmass)
       !
       !--take the sink-sink timestep specified by the get_forces routine
       !
       if (id==master) print*,' dt for sinks = ',C_force*dtsinksink
       dt      = C_force*dtsinksink !2.0/(nsteps)
       dtmax   = dt  ! required prior to derivs call, as used to set ibin
       !
       !--compute SPH forces
       !
       if (itest==2 .or. itest==3) then
          fxyzu(:,:) = 0.
          call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                      Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,t,0.,dtext_dum,pxyzu,dens,metrics)
       endif
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
          ntests = ntests + 1
          call checkval(epot_sinksink,-m1*m1*massr/a,epsilon(0.),nfailed(1),'potential energy')
          if (nfailed(1)==0) npass = npass + 1
          !
          !--check initial angular momentum on the two sinks is correct
          !
          ntests = ntests + 1
          call checkval(angtot,m1*m1*massr*sqrt(a/(m1 + m1*massr)),1.e6*epsilon(0.),nfailed(1),'angular momentum')
          if (nfailed(1)==0) npass = npass + 1
       endif
       !
       !--determine number of steps per orbit for information
       !
       mred    = m1*m1*massr/(m1 + m1*massr)
       omega   = sqrt(mred/a**3)
       nsteps  = int(2.*pi/omega/dt) + 1
       if (itest==2 .or. itest==3) then
          norbits = 10
       else
          norbits = 100
       endif
       if (id==master) print*,' nsteps per orbit = ',nsteps,' norbits = ',norbits
       nsteps = nsteps*norbits
       errmax = 0.
       dumpfile='test_00000'
       f_acc = 1.
       call init_step(npart,t,dtmax)
       do i=1,nsteps
          t = t + dt
          dtext = dt
          if (id==master .and. iverbose > 2) write(*,*) ' t = ',t,' dt = ',dt
          call step(npart,npart,t,dt,dtext,dtnew)
          call compute_energies(t)
          errmax = max(errmax,abs(etot - etotin))
          if (itest==2) then
             !   write(1,*) i,t,angtot,totmom,etot,accretedmass,epot
             !   print*,i,t,angtot,totmom,etot,accretedmass,epot
          endif
       enddo
       call compute_energies(t)
       nfailed(:) = 0
       select case(itest)
       case(3)
#ifdef IND_TIMESTEPS
          call checkval(angtot,angmomin,2.1e-6,nfailed(3),'angular momentum')
          call checkval(totmom,totmomin,5.e-6,nfailed(2),'linear momentum')
#else
          call checkval(angtot,angmomin,1.1e-6,nfailed(3),'angular momentum')
          call checkval(totmom,totmomin,4.e-14,nfailed(2),'linear momentum')
#endif
          call checkval(etotin+errmax,etotin,1.2e-2,nfailed(1),'total energy')
       case(2)
          call checkval(angtot,angmomin,3.e-7,nfailed(3),'angular momentum')
          call checkval(totmom,totmomin,6.e-14,nfailed(2),'linear momentum')
          tolen = 2.e-3
          if (gravity) tolen = 3.1e-3
          call checkval(etotin+errmax,etotin,tolen,nfailed(1),'total energy')
       case default
          call checkval(angtot,angmomin,3.e-14,nfailed(3),'angular momentum')
          call checkval(totmom,totmomin,epsilon(0.),nfailed(2),'linear momentum')
          call checkval(etotin+errmax,etotin,3.e-8,nfailed(1),'total energy')
       end select
       !
       !--check energy conservation
       !
       ntests = ntests + 3
       do i=1,3
          if (nfailed(i)==0) npass = npass + 1
       enddo
    enddo binary_tests

 endif testbinary

!
!  Test of softening between sinks
!
 testsoftening: if (test_softening) then
    if (id==master) write(*,"(/,a)") '--> testing softening in sink particle binary'
    nptmass = 0
    npart = 0
    npartoftype = 0
    m1    = 1.
    massr = 1.
    a     = 1.
    ecc   = 0.
    hacc1 = 0.
    hacc2 = 0.
    t     = 0.

    h_soft_sinksink = 0.8*a

    call set_units(mass=1.d0,dist=1.d0,G=1.d0)
    call set_binary(m1,massr,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,verbose=.false.)

    q   = a/h_soft_sinksink
    call kernel_softening(q*q,q,phisoft,fsoft)

    ! Test energy and momentum conservation
    m2 = m1*massr
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
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epot_sinksink,dtsinksink,0,0.)
    call compute_energies(t)
    etotin   = etot
    totmomin = totmom
    angmomin = angtot
    ntests = ntests + 1

    call checkval(epot,m1*m1*massr*(phisoft)/h_soft_sinksink,2.*epsilon(0.),nfailed(1),'potential energy')
    if (nfailed(1)==0) npass = npass + 1

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
!      write(1,'(10(es10.3,1x))') xyzmh_ptmass(:,1)
!      write(1,'(10(es10.3,1x))') xyzmh_ptmass(:,2)
       call compute_energies(t)
       errmax = max(errmax,abs(etot - etotin))
    enddo
    call compute_energies(t)
    nfailed(:) = 0
    call checkval(angtot,angmomin,2.e-14,nfailed(1),'angular momentum')
    call checkval(totmom,totmomin,tiny(0.),nfailed(2),'linear momentum')
    call checkval(etotin+errmax,etotin,2.e-9,nfailed(3),'total energy')
!    call checkval(      ,r_max,1.e-10,nfailed(4),'radius')

    ntests = ntests + 1
    if (all(nfailed(1:4)==0)) npass = npass + 1

    ! Reset sink softening
    h_soft_sinksink = 0.0

 endif testsoftening

!
!  Tests of accrete_particle routine
!
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:)  = 0.

 testaccretion: if (test_accretion) then
    if (id==master) write(*,"(/,a)") '--> testing accretion onto sink particles'
    nptmass = 1
    !--setup 1 point mass at (-5,-5,-5)
    xyzmh_ptmass(1:3,1)   = 1.
    xyzmh_ptmass(4,1)     = 10. ! mass of sink
    xyzmh_ptmass(ihacc,1) = 20. ! accretion radius
    vxyz_ptmass(1:3,1)    = -40.
    fxyz_ptmass(1:3,1)    = 40.
    massoftype(1)   = 10.
    !--setup 1 SPH particle at (5,5,5)
    if (id==master) then
       call set_particle_type(1,igas)
       npartoftype(igas) = 1
       npart        = 1
       xyzh(1:3,1)  = 5.
       xyzh(4,1)    = 0.01
       vxyzu(1:3,1) = 80.
       fxyzu(1:3,1) = 20.
    else
       npartoftype(igas) = 0
       npart        = 0
    endif
    xyzm_ptmass_old = xyzmh_ptmass(1:4,1:nptmass)
    vxyz_ptmass_old = vxyz_ptmass (1:3,1:nptmass)
    dr = sqrt(dot_product(xyzh(1:3,1) - xyzmh_ptmass(1:3,1),xyzh(1:3,1) - xyzmh_ptmass(1:3,1)))
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
       call ptmass_accrete(1,nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                           vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),fxyzu(1,i),fxyzu(2,i),fxyzu(3,i), &
                           igas,massoftype(igas),xyzmh_ptmass,vxyz_ptmass, &
                           accreted,dptmass_thread,t,1.0,ibin_wakei,ibin_wakei)
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

    if (id==master) then
       call checkval(accreted,.true.,nfailed(1),'accretion flag')
       !--check that h has been changed to indicate particle has been accreted
       call checkval(isdead_or_accreted(xyzh(4,1)),.true.,nfailed(2),'isdead_or_accreted flag')
    endif
    call checkval(xyzmh_ptmass(1,1),3.,tiny(0.),nfailed(3),'x(ptmass) after accretion')
    call checkval(xyzmh_ptmass(2,1),3.,tiny(0.),nfailed(4),'y(ptmass) after accretion')
    call checkval(xyzmh_ptmass(3,1),3.,tiny(0.),nfailed(5),'z(ptmass) after accretion')
    call checkval(vxyz_ptmass(1,1),20.,tiny(0.),nfailed(6),'vx(ptmass) after accretion')
    call checkval(vxyz_ptmass(2,1),20.,tiny(0.),nfailed(7),'vy(ptmass) after accretion')
    call checkval(vxyz_ptmass(3,1),20.,tiny(0.),nfailed(8),'vz(ptmass) after accretion')
    call checkval(fxyz_ptmass(1,1),30.,tiny(0.),nfailed(9), 'fx(ptmass) after accretion')
    call checkval(fxyz_ptmass(2,1),30.,tiny(0.),nfailed(10),'fy(ptmass) after accretion')
    call checkval(fxyz_ptmass(3,1),30.,tiny(0.),nfailed(11),'fz(ptmass) after accretion')

    ntests = ntests + 4
    if (all(nfailed(1:2)==0)) npass = npass + 1
    if (all(nfailed(3:5)==0)) npass = npass + 1
    if (all(nfailed(6:8)==0)) npass = npass + 1
    if (all(nfailed(9:11)==0)) npass = npass + 1

    !--compute energies after accretion event
    nfailed(:) = 0
    call compute_energies(t)
    call checkval(angtot,angmomin,1.e-10,nfailed(3),'angular momentum')
    call checkval(totmom,totmomin,epsilon(0.),nfailed(2),'linear momentum')
    !call checkval(etot,etotin,1.e-6,'total energy',nfailed(1))
    ntests = ntests + 2
    if (nfailed(3)==0) npass = npass + 1
    if (nfailed(2)==0) npass = npass + 1

 endif testaccretion
!
!  Test sink particle creation
!
 testcreatesink: if (test_createsink) then
    do itest=1,2
       select case(itest)
       case(2)
          if (id==master) write(*,"(/,a)") '--> testing sink particle creation (sin)'
       case default
          if (id==master) write(*,"(/,a)") '--> testing sink particle creation (uniform density)'
       end select
       xyzmh_ptmass(:,:) = 0.
       vxyz_ptmass(:,:) = 0.
       nptmass = 0

       npart = 0
       npartoftype(:) = 0
       massoftype = 0.

       xyzh(:,:)  = 0.
       vxyzu(:,:) = 0.
       fxyzu(:,:) = 0.
       fext(:,:)  = 0.
       Bevol(:,:) = 0.
       iverbose   = 1
       call set_boundary(-1.,1.,-1.,1.,-1.,1.)
       !
       ! set up gas particles in a uniform sphere with radius R=0.2
       !
       psep = 0.05  ! required as a variable since this may change under conditions not requested here
       if (itest==2) then
          ! use random so particle with maximum density is unique
          call set_sphere('cubic',id,master,0.,0.2,psep,hfact,npartoftype(igas),xyzh,rhofunc=gaussianr)
       else
          call set_sphere('cubic',id,master,0.,0.2,psep,hfact,npartoftype(igas),xyzh)
       endif
       totmass = 1.0
       massoftype(igas) = totmass/real(npartoftype(igas))
       npart = npartoftype(igas)
       !
       ! give inward radial velocities
       !
       itestp = npart
       r2min = huge(r2min)
       xcofm(:) = 0.
       do i=1,npart
          r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
          if (r2 < r2min) then
             itestp = i
             r2min = r2
          endif
          xcofm = xcofm + xyzh(1:3,i)
       enddo

       xcofm = massoftype(igas)*xcofm/totmass
       if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
       !
       ! set up tree for neighbour finding
       ! and make sure that gravitational potential energy has been computed
       !
       tree_accuracy = 0.
       icreate_sinks = 1
       iverbose = 1
       call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                   Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,0.,0.,dtext_dum,pxyzu,dens,metrics)
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
             call checkval(rhomax,rhomax_test,epsilon(0.),nfailed(1),'rhomax')
          else
             itestp = -1 ! set itest = -1 on other threads
             call checkval(ipart_rhomax,-1,0,nfailed(1),'ipart_rhomax')
          endif
          ntests = ntests + 1
          if (nfailed(1)==0) npass = npass + 1
       endif
       !
       ! check energies before insertion of sink
       !
       call compute_energies(t)
       !print*,' got ETOT = ',etot,totmom,epot
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
                          massoftype,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,0.)
       !
       ! check that creation succeeded
       !
       nfailed(:) = 0
       call checkval(nptmass,1,0,nfailed(1),'nptmass=1')
       ntests = ntests + 1
       if (all(nfailed==0)) npass = npass + 1

       !
       ! check that linear and angular momentum and energy is conserved
       !
       nfailed(:) = 0
       call compute_energies(t)
       !print*,' got ETOT = ',etot,totmom,epot
       call checkval(angtot,angmomin,1.e-10,nfailed(3),'angular momentum')
       call checkval(totmom,totmomin,epsilon(0.),nfailed(2),'linear momentum')
       !call checkval(etot,etotin,1.e-6,nfailed(1),'total energy')
       ntests = ntests + 1
       if (all(nfailed==0)) npass = npass + 1
       iverbose = 0
       call finish_ptmass(nptmass)
    enddo
 endif testcreatesink

 !--reset stuff
 nptmass = 0

 ! clean up temporary files
 itmp = 201
 do i=iskfile+1,iskfile+2    ! we used 2 point masses in tests above
    close(i,iostat=ierr)     ! close file unit numbers if they are already open
    write(filename,"(i3)") i
    filename = 'fort.'//trim(adjustl(filename))
    open(unit=itmp,file=filename,status='old',iostat=ierr)
    close(itmp,status='delete',iostat=ierr)
 enddo
 open(unit=itmp,file='SinkSink0001N00.ev',status='old',iostat=ierr)
 close(itmp,status='delete',iostat=ierr)

 if (id==master) write(*,"(/,a)") '<-- PTMASS TEST COMPLETE'

end subroutine test_ptmass

real function gaussianr(r)
 real, intent(in) :: r

 gaussianr = exp(-(r/0.05)**2) !1./(r**2 + 0.0001**2)

end function gaussianr

end module testptmass
