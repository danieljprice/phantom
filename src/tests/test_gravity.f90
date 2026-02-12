!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testgravity
!
! Unit tests of self-gravity
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: checksetup, deriv, dim, directsum, energies, eos, io,
!   kdtree, kernel, mpibalance, mpidomain, mpiutils, neighkdtree, options,
!   part, physcon, ptmass, random, setplummer, setup_params,
!   sort_particles, sortutils, spherical, table_utils, testapr, testutils,
!   timing, units
!
 use io, only:id,master
 implicit none
 public :: test_gravity
 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests for Newtonian gravity (i.e. Poisson solver)
!+
!-----------------------------------------------------------------------
subroutine test_gravity(ntests,npass,string)
 use dim, only:gravity
 use testapr, only:setup_apr_region_for_test
 integer,          intent(inout) :: ntests,npass
 character(len=*), intent(in)    :: string
 logical :: testdirectsum,test_mom,testtaylorseries,testall,test_plummer
 logical :: plot_plummer

 testdirectsum    = .false.
 testtaylorseries = .false.
 test_mom         = .false.
 testall          = .false.
 test_plummer     = .false.
 plot_plummer     = .false.
 select case(string)
 case('taylorseries')
    testtaylorseries = .true.
 case('directsum')
    testdirectsum = .true.
 case('fmm')
    test_mom = .true.
 case('spheres','plummer','hernquist')
    test_plummer = .true.
 case('plotplummer')
    plot_plummer = .true.
 case default
    testall = .true.
 end select
 if (gravity) then
    if (id==master) write(*,"(/,a,/)") '--> TESTING SELF-GRAVITY'
    !
    !--unit test of Taylor series expansions in the treecode
    !
    if (testtaylorseries .or. testall) call test_taylorseries(ntests,npass)
    !
    !--unit tests of treecode gravity by direct summation
    !
    if (testdirectsum .or. testall) call test_directsum(ntests,npass)
    !
    !--unit tests of FMM momentum conservation
    !
    if (test_mom .or. testall) call test_FMM(ntests,npass)
    !
    !--unit tests of Plummer and Hernquist spheres
    !
    if (test_plummer .or. testall) call test_spheres(ntests,npass)
    !
    !--Plot routine of Plummer and Homogeneous sphere (store data to be plotted)
    !
    if (plot_plummer) call plot_SFMM()

    if (id==master) write(*,"(/,a)") '<-- SELF-GRAVITY TESTS COMPLETE'
 else
    if (id==master) write(*,"(/,a)") '--> SKIPPING SELF-GRAVITY TESTS (need -DGRAVITY)'
 endif

end subroutine test_gravity

!-----------------------------------------------------------------------
!+
!  Unit tests of the Taylor series expansion about local and distant nodes
!+
!-----------------------------------------------------------------------
subroutine test_taylorseries(ntests,npass)
 use kdtree,    only:compute_M2L,expand_fgrav_in_taylor_series
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(18),i,npnode
 real :: xposi(3),xposj(3),x0(3),dx(3),fexact(3),f0(3)
 real :: xposjd(3,3)
 real :: fnode(20),quads(6)
 real :: dr,dr2,phi,phiexact,pmassi,totmass

 if (id==master) write(*,"(/,a)") '--> testing taylor series expansion about current node'
 totmass = 5.
 xposi = (/0.05,-0.04,-0.05/)  ! position to evaluate the force at
 xposj = (/1., 1., 1./)       ! position of distant node
 x0 = 0.                      ! position of nearest node centre

 call get_dx_dr(xposi,xposj,dx,dr)
 fexact = -totmass*dr**3*dx   ! exact force between i and j
 phiexact = -totmass*dr       ! exact potential between i and j

 call get_dx_dr(x0,xposj,dx,dr)
 fnode = 0.
 quads = 0.
 call compute_M2L(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

 dx = xposi - x0   ! perform expansion about x0
 call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),phi)
 !print*,'           exact force = ',fexact,' phi = ',phiexact
 !print*,'       force at origin = ',fnode(1:3), ' phi = ',fnode(20)
 !print*,'force w. taylor series = ',f0, ' phi = ',phi
 nfailed(:) = 0
 call checkval(f0(1),fexact(1),3.e-4,nfailed(1),'fx taylor series about f0')
 call checkval(f0(2),fexact(2),1.1e-4,nfailed(2),'fy taylor series about f0')
 call checkval(f0(3),fexact(3),9.e-5,nfailed(3),'fz taylor series about f0')
 call checkval(phi,phiexact,8.e-4,nfailed(4),'phi taylor series about f0')
 call update_test_scores(ntests,nfailed,npass)

 if (id==master) write(*,"(/,a)") '--> testing taylor series expansion about distant node'
 totmass = 5.
 npnode = 3
 xposjd(:,1) = (/1.03, 0.98, 1.01/)       ! position of distant particle 1
 xposjd(:,2) = (/0.95, 1.01, 1.03/)       ! position of distant particle 2
 xposjd(:,3) = (/0.99, 0.95, 0.95/)       ! position of distant particle 3
 xposj = 0.
 pmassi = totmass/real(npnode)
 do i=1,npnode
    xposj = xposj + pmassi*xposjd(:,i)     ! centre of mass of distant node
 enddo
 xposj = xposj/totmass
 !print*,' centre of mass of distant node = ',xposj
 !--compute quadrupole moments
 quads = 0.
 do i=1,npnode
    dx(:) = xposjd(:,i) - xposj
    dr2   = dot_product(dx,dx)
    quads(1) = quads(1) + pmassi*(dx(1)*dx(1))
    quads(2) = quads(2) + pmassi*(dx(1)*dx(2))
    quads(3) = quads(3) + pmassi*(dx(1)*dx(3))
    quads(4) = quads(4) + pmassi*(dx(2)*dx(2))
    quads(5) = quads(5) + pmassi*(dx(2)*dx(3))
    quads(6) = quads(6) + pmassi*(dx(3)*dx(3))
 enddo

 x0 = 0.      ! position of nearest node centre
 xposi = x0   ! position to evaluate the force at
 fexact = 0.
 phiexact = 0.
 do i=1,npnode
    dx = xposi - xposjd(:,i)
    dr = 1./sqrt(dot_product(dx,dx))
    fexact = fexact - dr**3*dx   ! exact force between i and j
    phiexact = phiexact - dr     ! exact force between i and j
 enddo
 fexact = fexact*pmassi
 phiexact = phiexact*pmassi

 call get_dx_dr(x0,xposj,dx,dr)
 fnode = 0.
 call compute_M2L(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

 dx = xposi - x0   ! perform expansion about x0
 call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),phi)
 !print*,'           exact force = ',fexact,' phi = ',phiexact
 !print*,'       force at origin = ',fnode(1:3), ' phi = ',fnode(20)
 !print*,'force w. taylor series = ',f0, ' phi = ',phi
 nfailed(:) = 0
 call checkval(f0(1),fexact(1),8.7e-5,nfailed(1),'fx taylor series about f0')
 call checkval(f0(2),fexact(2),1.5e-6,nfailed(2),'fy taylor series about f0')
 call checkval(f0(3),fexact(3),1.6e-5,nfailed(3),'fz taylor series about f0')
 call checkval(phi,phiexact,5.9e-6,nfailed(4),'phi taylor series about f0')
 call update_test_scores(ntests,nfailed,npass)

 if (id==master) write(*,"(/,a)") '--> testing taylor series expansion about both current and distant nodes'
 x0 = 0.                      ! position of nearest node centre
 xposi = (/0.05,0.05,-0.05/)  ! position to evaluate the force at
 fexact = 0.
 phiexact = 0.
 do i=1,npnode
    dx = xposi - xposjd(:,i)
    dr = 1./sqrt(dot_product(dx,dx))
    fexact = fexact - dr**3*dx   ! exact force between i and j
    phiexact = phiexact - dr     ! exact force between i and j
 enddo
 fexact = fexact*pmassi
 phiexact = phiexact*pmassi

 dx = x0 - xposj
 dr = 1./sqrt(dot_product(dx,dx))  ! compute approx force between node and j
 fnode = 0.
 call compute_M2L(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

 dx = xposi - x0   ! perform expansion about x0
 call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),phi)
 !print*,'           exact force = ',fexact,' phi = ',phiexact
 !print*,'       force at origin = ',fnode(1:3), ' phi = ',fnode(20)
 !print*,'force w. taylor series = ',f0, ' phi = ',phi
 nfailed(:) = 0
 call checkval(f0(1),fexact(1),1.3e-4,nfailed(1),'fx taylor series about f0')
 call checkval(f0(2),fexact(2),1.4e-4,nfailed(2),'fy taylor series about f0')
 call checkval(f0(3),fexact(3),3.2e-4,nfailed(3),'fz taylor series about f0')
 call checkval(phi,phiexact,9.7e-4,nfailed(4),'phi taylor series about f0')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_taylorseries

!-----------------------------------------------------------------------
!+
!  Unit tests of the tree code gravity, checking it agrees with
!  gravity computed via direct summation
!+
!-----------------------------------------------------------------------
subroutine test_directsum(ntests,npass)
 use io,              only:id,master,nprocs
 use dim,             only:maxp,maxptmass,mpi,use_apr,use_sinktree,maxpsph
 use part,            only:init_part,npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu, &
                           gradh,poten,iphase,isetphase,maxphase,labeltype,&
                           nptmass,xyzmh_ptmass,fxyz_ptmass,dsdt_ptmass,ibelong,&
                           fxyz_ptmass_tree,istar,shortsinktree
 use eos,             only:polyk,gamma
 use options,         only:ieos,alpha,alphau,alphaB,tolh
 use spherical,       only:set_sphere
 use deriv,           only:get_derivs_global
 use physcon,         only:pi
 use timing,          only:getused,printused
 use directsum,       only:directsum_grav
 use energies,        only:compute_energies,epot
 use kdtree,          only:tree_accuracy
 use testutils,       only:checkval,checkvalbuf_end,update_test_scores
 use ptmass,          only:get_accel_sink_sink,get_accel_sink_gas,h_soft_sinksink
 use mpiutils,        only:reduceall_mpi,bcast_mpi
 use neighkdtree,     only:build_tree
 use sort_particles,  only:sort_part_id
 use mpibalance,      only:balancedomains
 use testapr,         only:setup_apr_region_for_test
 use setup_params,    only:npart_total

 integer, intent(inout) :: ntests,npass
 integer :: nfailed(18),boundi,boundf
 integer :: maxvxyzu,nx,np,i,k,merge_n,merge_ij(maxptmass),nfgrav
 real :: psep,totvol,totmass,rhozero,tol,pmassi,fsum(3)
 real :: time,rmin,rmax,phitot,dtsinksink,fonrmax,phii,epot_gas_sink
 real(kind=4) :: t1,t2
 real :: epoti,tree_acc_prev
 real, allocatable :: fgrav(:,:),fxyz_ptmass_gas(:,:)

 maxvxyzu = size(vxyzu(:,1))
 tree_acc_prev = tree_accuracy
 do k = 1,6
    if (labeltype(k)/='bound') then
       if (id==master) write(*,"(/,3a)") '--> testing gravity force in densityforce for ',labeltype(k),' particles'
!
!--general parameters
!
       time  = 0.
       hfact = 1.2
       gamma = 5./3.
       rmin  = 0.
       rmax  = 1.
       ieos  = 2
       tree_accuracy = 0.5
!
!--setup particles
!
       call init_part()
       np       = 1000
       totvol   = 4./3.*pi*rmax**3
       nx       = int(np**(1./3.))
       psep     = totvol**(1./3.)/real(nx)
       psep     = 0.18
       npart    = 0
       npart_total = 0
       ! only set up particles on master, otherwise we will end up with n duplicates
       if (id==master) then
          call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart,xyzh,npart_total,np_requested=np)
       endif
       np       = npart
!
!--set particle properties
!
       totmass        = 1.
       rhozero        = totmass/totvol
       npartoftype(:) = 0
       npartoftype(k) = int(reduceall_mpi('+',npart),kind=kind(npartoftype))
       massoftype(:)  = 0.0
       massoftype(k)  = totmass/npartoftype(k)
       if (maxphase==maxp) then
          do i=1,npart
             iphase(i) = isetphase(k,iactive=.true.)
          enddo
       endif
!
!--call apr setup if using it - this must be called after massoftype is set
!       we're not using this right now, this test fails as is
!       if (use_apr) call setup_apr_region_for_test()

!
!--set thermal terms and velocity to zero, so only force is gravity
!
       polyk      = 0.
       vxyzu(:,:) = 0.
!
!--make sure AV is off
!
       alpha  = 0.
       alphau = 0.
       alphaB = 0.
       tolh = 1.e-5

       fxyzu = 0.0
!
!--call derivs to get everything initialised
!
       call get_derivs_global()

!
!--reset force to zero
!
       fxyzu = 0.0
!
!--move particles to master and sort for direct summation
!
       if (mpi) then
          ibelong(:) = 0
          call balancedomains(npart)
       endif
       call sort_part_id
!
!--allocate array for storing direct sum gravitational force
!
       allocate(fgrav(maxvxyzu,npart))
       fgrav = 0.0
!
!--compute gravitational forces by direct summation
!
       if (id == master) then
          call directsum_grav(xyzh,gradh,fgrav,phitot,npart)
       endif
!
!--send phitot to all tasks
!
       call bcast_mpi(phitot)
!
!--calculate derivatives
!
       call getused(t1)
       call get_derivs_global()
       call getused(t2)
       if (id==master) call printused(t1)
!
!--move particles to master and sort for test comparison
!
       if (mpi) then
          ibelong(:) = 0
          call balancedomains(npart)
       endif
       call sort_part_id
!
!-- sum all the force for conservation checks
!
       fsum=0.
       do i=1,npart
          fsum(1) = fsum(1) + fxyzu(1,i)
          fsum(2) = fsum(2) + fxyzu(2,i)
          fsum(3) = fsum(3) + fxyzu(3,i)
       enddo

       fsum(:) = fsum(:)*massoftype(k)
!
!--compare the results
!
       call checkval(npart,fxyzu(1,:),fgrav(1,:),7.2e-3,nfailed(1),'fgrav(x)')
       call checkval(npart,fxyzu(2,:),fgrav(2,:),6.e-3,nfailed(2),'fgrav(y)')
       call checkval(npart,fxyzu(3,:),fgrav(3,:),9.4e-3,nfailed(3),'fgrav(z)')
       call checkval(fsum(1), 0., 2.5e-17, nfailed(4),'fsum(x)')
       call checkval(fsum(2), 0., 2.5e-17, nfailed(5),'fsum(y)')
       call checkval(fsum(3), 0., 2.9e-17, nfailed(6),'fsum(z)')
       deallocate(fgrav)
       epoti = 0.
       do i=1,npart
          epoti = epoti + poten(i)
       enddo
       epoti = reduceall_mpi('+',epoti)
       call checkval(epoti,phitot,5.2e-4,nfailed(7),'potential')
       call checkval(epoti,-3./5.*totmass**2/rmax,3.6e-2,nfailed(8),'potential=-3/5 GMM/R')
       ! check that potential energy computed via compute_energies is also correct
       call compute_energies(0.)
       call checkval(epot,phitot,5.2e-4,nfailed(9),'epot in compute_energies')
       call update_test_scores(ntests,nfailed(1:9),npass)
    endif
 enddo

!--test that the same results can be obtained from a cloud of sink particles
!  with softening lengths equal to the original SPH particle smoothing lengths
!
 !
 !--general parameters
 !
 time  = 0.
 hfact = 1.2
 gamma = 5./3.
 rmin  = 0.
 rmax  = 1.
 ieos  = 2
 tree_accuracy = 0.5
 !
 !--setup particles
 !
 call init_part()
 np       = 1000
 totvol   = 4./3.*pi*rmax**3
 nx       = int(np**(1./3.))
 psep     = totvol**(1./3.)/real(nx)
 psep     = 0.18
 npart    = 0

 ! only set up particles on master, otherwise we will end up with n duplicates
 if (id==master) then
    call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart,xyzh,npart_total,np_requested=np)
 endif
 np       = npart
 !
 !--set particle properties
 !
 totmass        = 1.
 rhozero        = totmass/totvol
 npartoftype(:) = 0
 npartoftype(istar) = int(reduceall_mpi('+',npart),kind=kind(npartoftype))
 massoftype(:)  = 0.0
 massoftype(istar)  = totmass/npartoftype(istar)
 if (maxphase==maxp) then
    do i=1,npart
       iphase(i) = isetphase(istar,iactive=.true.) ! set all particles to star to avoid comp gas force (only grav here)
    enddo
 endif
 if (maxptmass >= npart) then
    if (use_sinktree) then
       if (id==master) write(*,"(/,3a)") '--> testing gravity in uniform cloud of softened sink particles (SinkInTree)'
    else
       if (id==master) write(*,"(/,3a)") '--> testing gravity in uniform cloud of softened sink particles (direct)'
    endif
!
!--move particles to master for sink creation
!
    if (mpi) then
       ibelong(:) = 0
       call balancedomains(npart)
    endif
!
!--sort particles so that they can be compared at the end
!
    call sort_part_id

    pmassi = totmass/reduceall_mpi('+',npart)
    call copy_gas_particles_to_sinks(npart,nptmass,xyzh,xyzmh_ptmass,pmassi)
    h_soft_sinksink = hfact*psep
!
!--compute direct sum for comparison, but with fixed h and hence gradh terms switched off
!
    do i=1,npart
       xyzh(4,i)  = h_soft_sinksink
       gradh(1,i) = 1.
       gradh(2,i) = 0.
       vxyzu(:,i) = 0.
    enddo
    allocate(fgrav(maxvxyzu,npart))
    fgrav = 0.0
    call directsum_grav(xyzh,gradh,fgrav,phitot,npart)
    call bcast_mpi(phitot)
!
!--compute gravity on the sink particles
!
    shortsinktree(1:nptmass,1:nptmass) = 1
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epoti,&
                             dtsinksink,0,0.,merge_ij,merge_n,dsdt_ptmass)
    call bcast_mpi(epoti)
!
!--compare the results
!
    tol = 1.e-14
    call checkval(npart,fxyz_ptmass(1,:),fgrav(1,:),tol,nfailed(1),'fgrav(x)')
    call checkval(npart,fxyz_ptmass(2,:),fgrav(2,:),tol,nfailed(2),'fgrav(y)')
    call checkval(npart,fxyz_ptmass(3,:),fgrav(3,:),tol,nfailed(3),'fgrav(z)')
    call checkval(epoti,phitot,8e-3,nfailed(4),'potential')
    call checkval(epoti,-3./5.*totmass**2/rmax,4.1e-2,nfailed(5),'potential=-3/5 GMM/R')
    call update_test_scores(ntests,nfailed(1:5),npass)

!
!--now perform the same test, but with HALF the cloud made of sink particles
!  and HALF the cloud made of gas particles. Do not re-evaluate smoothing lengths
!  so that the results should be identical to the previous test
!
    if (use_sinktree) then
       if (id==master) write(*,"(/,3a)") &
       '--> testing softened gravity in uniform sphere with half sinks and half gas (SinkInTree)'
    else
       if (id==master) write(*,"(/,3a)") &
       '--> testing softened gravity in uniform sphere with half sinks and half gas (direct)'
    endif

!--sort the particles by ID so that the first half will have the same order
!  even after half the particles have been converted into sinks. This sort is
!  not really necessary because the order shouldn't have changed since the
!  last test because derivs hasn't been called since.
    call sort_part_id
    call copy_half_gas_particles_to_sinks(npart,nptmass,xyzh,xyzmh_ptmass,pmassi,hfact*psep)
    !nptmass = 0

    if (mpi) then
       if (use_sinktree) then
          ibelong((maxpsph)+1:maxp) = -1
          boundi = (maxpsph)+(nptmass / nprocs)*id
          boundf = (maxpsph)+(nptmass / nprocs)*(id+1)
          if (id == nprocs-1) boundf = boundf + mod(nptmass,nprocs)
          ibelong(boundi+1:boundf) = id
          fxyz_ptmass_tree = 0.
       endif
    endif

    print*,' Using ',npart,' SPH particles and ',nptmass,' point masses'
    call get_derivs_global(icall=0) ! icall = 0 refresh tree cache used for h1j in the force routine

    epoti = 0.0
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,epoti,&
                             dtsinksink,0,0.,merge_ij,merge_n,dsdt_ptmass)
!
!--prevent double counting of sink contribution to potential due to MPI
!
    if (id /= master) epoti = 0.0

    if (use_sinktree) then
       epot_gas_sink = 0.
       do i=1,npart
          epoti = epoti + poten(i)
       enddo
       do i=1,nptmass
          epoti = epoti + poten(i+maxpsph)
       enddo

       fxyz_ptmass(1:3,1:nptmass) = fxyz_ptmass(1:3,1:nptmass) + fxyz_ptmass_tree(1:3,1:nptmass)
       epoti  = reduceall_mpi('+',epoti)
    else
!
!--allocate an array for the gas contribution to sink acceleration
!
       allocate(fxyz_ptmass_gas(size(fxyz_ptmass,dim=1),nptmass))
       fxyz_ptmass_gas = 0.0

       epot_gas_sink = 0.0
       do i=1,npart
          call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),&
                               xyzmh_ptmass,fxyzu(1,i),fxyzu(2,i),fxyzu(3,i),&
                               phii,pmassi,fxyz_ptmass_gas,dsdt_ptmass,fonrmax,dtsinksink)
          epot_gas_sink = epot_gas_sink + pmassi*phii
          epoti = epoti + poten(i)
       enddo
!
!--the gas contribution to sink acceleration has to be added afterwards to
!  prevent double counting the sink contribution when calling reduceall_mpi
!
       fxyz_ptmass_gas = reduceall_mpi('+',fxyz_ptmass_gas)
       fxyz_ptmass(:,1:nptmass) = fxyz_ptmass(:,1:nptmass) + fxyz_ptmass_gas(:,1:nptmass)
       deallocate(fxyz_ptmass_gas)
!
!--sum up potentials across MPI tasks
!
       epoti         = reduceall_mpi('+',epoti)
       epot_gas_sink = reduceall_mpi('+',epot_gas_sink)
    endif
!
!--move particles to master for comparison
!
    if (mpi) then
       ibelong(:) = 0
       call balancedomains(npart)
    endif
    call sort_part_id

    call checkval(npart,fxyzu(1,:),fgrav(1,:),5.e-2,nfailed(1),'fgrav(x)')
    call checkval(npart,fxyzu(2,:),fgrav(2,:),6.e-2,nfailed(2),'fgrav(y)')
    call checkval(npart,fxyzu(3,:),fgrav(3,:),9.4e-2,nfailed(3),'fgrav(z)')

!
!--fgrav doesn't exist on worker tasks, so it needs to be sent from master
!
    call bcast_mpi(npart)
    if (id == master) nfgrav = size(fgrav,dim=2)
    call bcast_mpi(nfgrav)
    if (id /= master) then
       deallocate(fgrav)
       allocate(fgrav(maxvxyzu,nfgrav))
    endif
    call bcast_mpi(fgrav)
    call checkval(nptmass,fxyz_ptmass(1,1:nptmass),fgrav(1,npart+1:2*npart),2.3e-2,nfailed(4),'fgrav(xsink)')
    call checkval(nptmass,fxyz_ptmass(2,1:nptmass),fgrav(2,npart+1:2*npart),2.9e-2,nfailed(5),'fgrav(ysink)')
    call checkval(nptmass,fxyz_ptmass(3,1:nptmass),fgrav(3,npart+1:2*npart),3.7e-2,nfailed(6),'fgrav(zsink)')

    call checkval(epoti+epot_gas_sink,phitot,8e-3,nfailed(7),'potential')
    call checkval(epoti+epot_gas_sink,-3./5.*totmass**2/rmax,4.1e-2,nfailed(8),'potential=-3/5 GMM/R')
    call update_test_scores(ntests,nfailed(1:8),npass)
    deallocate(fgrav)
 endif
!
!--clean up doggie-doos
!
 npartoftype(:) = 0
 massoftype(:) = 0.
 tree_accuracy = tree_acc_prev
 fxyzu = 0.
 vxyzu = 0.

end subroutine test_directsum

!-----------------------------------------------------------------------
!+
!   test that we conserve linear momentum with the symmetrical FMM
!+
!-----------------------------------------------------------------------
subroutine test_FMM(ntests,npass)
 use io,        only:id,master,iverbose
 use part,      only:npart,npartoftype,xyzh,massoftype,hfact,&
                       init_part,fxyzu,istar,set_particle_type,&
                       ibelong
 use mpidomain, only:i_belong
 use mpiutils,  only:reduceall_mpi
 use mpibalance,only:balancedomains
 use options,   only:ieos
 use physcon,   only:solarr,solarm,pi
 use units,     only:set_units
 use eos,       only:gamma
 use kdtree,    only:tree_accuracy
 use checksetup, only:check_setup
 use spherical, only:set_sphere
 use deriv,     only: get_derivs_global
 use testutils,       only:checkval,checkvalbuf_end,update_test_scores
 use sort_particles,  only:sort_part_id
 use dim, only:maxp,maxphase,mpi

 integer, intent(inout) :: ntests,npass
 real :: x0(3),rmin,rmax,nx,psep,totvol,time,fsum(3)
 integer(kind=8) ::npart_total
 integer :: np
 integer :: nfail(3),i

 if (id==master) write(*,"(/,a)") '--> testing linear momentum conservation with symmetric fmm'
 if (mpi) then
    if (id==master) write(*,"(/,a)") '--> skipped... No sym FMM with MPI'
    return
 endif
 npart = 0
 npartoftype = 0
 massoftype = 0.
 iverbose = 0

 x0 = 0.
 fsum = 0.

 !
 !--general parameters
 !
 time  = 0.
 hfact = 1.2
 gamma = 5./3.
 rmin  = 0.
 rmax  = 1.
 ieos  = 2
 tree_accuracy = 0.5
 !
 !--setup particles
 !
 call init_part()
 np       = 10000
 totvol   = 4./3.*pi*rmax**3
 nx       = int(np**(1./3.))
 psep     = totvol**(1./3.)/real(nx)
 psep     = 0.18
 npart    = 0
 npart_total = 0

 ! do this test twice, to check the second star relaxes...
 do i=1,2
    if (i==2) x0 = [20.,0.,0.]
    ! only set up particles on master, otherwise we will end up with n duplicates
    if (id==master) then
       call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart,xyzh,npart_total,np_requested=np,xyz_origin=x0)
    endif
 enddo
 npartoftype(:) = 0
 npartoftype(istar) = int(reduceall_mpi('+',npart),kind=kind(npartoftype))
 massoftype(:)  = 0.
 massoftype(istar)  = 1./npartoftype(istar)

 if (maxphase==maxp) then
    do i=1,npart
       call set_particle_type(i,istar)
    enddo
 endif

 call get_derivs_global()

 !
 !--move particles to master and sort for test comparison
 !
 if (mpi) then
    ibelong(:) = 0
    call balancedomains(npart)
 endif
 call sort_part_id

 do i=1,npart
    fsum(1) = fsum(1) + fxyzu(1,i)
    fsum(2) = fsum(2) + fxyzu(2,i)
    fsum(3) = fsum(3) + fxyzu(3,i)
 enddo
 fsum = fsum*massoftype(istar)
 call checkval(fsum(1),0.,2.e-16,nfail(1),"momentum conservation x")
 call checkval(fsum(2),0.,2.e-16,nfail(2),"momentum conservation y")
 call checkval(fsum(3),0.,2.e-16,nfail(3),"momentum conservation z")
 call update_test_scores(ntests,nfail,npass)

end subroutine test_FMM

!-----------------------------------------------------------------------
!+
!   Unit tests of Plummer and Hernquist spheres
!+
!-----------------------------------------------------------------------
subroutine test_spheres(ntests,npass)
 use testutils,  only:checkval,update_test_scores
 use setplummer, only:iprofile_plummer,iprofile_hernquist
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a)") '--> testing Plummer and Hernquist spheres'

 call test_sphere(ntests,npass,iprofile_plummer)
 !call test_sphere(ntests,npass,iprofile_hernquist)

 if (id==master) write(*,"(/,a)") '<-- Plummer and Hernquist spheres test complete'

end subroutine test_spheres

!-----------------------------------------------------------------------
!+
!  Monte Carlo MASE test for a specified spherical density profile
!+
!-----------------------------------------------------------------------
subroutine test_sphere(ntests,npass,iprofile)
 use dim,         only:maxp
 use deriv,       only:get_derivs_global
 use eos,         only:gamma,polyk
 use mpiutils,    only:reduceall_mpi
 use options,     only:ieos,alpha,alphau,alphaB,tolh
 use part,        only:init_part,npart,xyzh,fxyzu,hfact,&
                       npartoftype,massoftype,istar,maxphase,iphase,isetphase
 use setup_params,only:npart_total
 use testutils,   only:checkval,update_test_scores
 use setplummer,  only:get_accel_profile,profile_label,radius_from_mass,density_profile
 use spherical,   only:set_sphere,iseed_mc
 use kdtree,      only:tree_accuracy
 use kernel,      only:hfact_default
 use table_utils, only:linspace
 use mpidomain,   only:i_belong
 use io,          only:id,master,iverbose
 integer, intent(inout) :: ntests,npass
 integer, intent(in)    :: iprofile
 integer :: nfailed(1)
 integer :: npart_target,nrealisations,i,ireal
 real :: err_sum,ref_sum,err_local,ref_local
 real :: mase,mase_tol,total_samples
 real :: rsoft,mass_total,cut_fraction
 real :: rmin,rmax,psep
 real :: acc_exact(3),diff(3)
 character(len=32) :: label
 integer, parameter :: ntab = 1000
 real :: rgrid(ntab),rhotab(ntab)

 label = profile_label(iprofile)

 npart_target = 10000
 total_samples = 1.0e5
 nrealisations = int(total_samples/real(npart_target))

 if (id==master) then
    write(*,"(1x,a,i8,a,i8)") 'Monte Carlo '//trim(label)//' test: N = ',npart_target, &
                                     ', nrealisations = ',nrealisations
 endif

 call init_part()
 hfact = hfact_default
 tree_accuracy = 0.5
 gamma = 5./3.
 polyk = 0.
 ieos  = 11
 alpha  = 0.; alphau = 0.; alphaB = 0.
 tolh = 1.e-5
 rsoft = 1.0
 mass_total = 1.0

 ! construct tables for radius and density
 cut_fraction = 0.999
 rmin = 0.
 rmax = radius_from_mass(iprofile,cut_fraction,rsoft)
 call linspace(rgrid,0.,rmax)
 do i=1,ntab
    rhotab(i) = density_profile(iprofile,rgrid(i),rsoft,mass_total)
 enddo

 psep = rmax/real(ntab) ! this is not used for random placement anyway
 iverbose = 0
 err_sum = 0.
 ref_sum = 0.

 do ireal=1,nrealisations
    iseed_mc = ireal
    npart = 0
    npart_total = 0
    call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart, &
                    xyzh,npart_total,rhotab=rhotab,rtab=rgrid,exactN=.true.,&
                    np_requested=npart_target,mask=i_belong,verbose=.false.)

    massoftype(istar) = mass_total/real(npart_total)
    npartoftype(istar) = npart
    if (maxphase==maxp) then
       iphase(1:npart) = isetphase(istar,iactive=.true.)
    endif

    call get_derivs_global()

    err_local = 0.
    ref_local = 0.
    do i=1,npart
       call get_accel_profile(iprofile,xyzh(1:3,i),rsoft,mass_total,acc_exact)
       diff = fxyzu(1:3,i) - acc_exact
       err_local = err_local + dot_product(diff,diff)
       ref_local = ref_local + dot_product(acc_exact,acc_exact)
    enddo
    err_local = reduceall_mpi('+',err_local)
    ref_local = reduceall_mpi('+',ref_local)
    if (iverbose > 0 .and. id==master) then
       print*,' realisation ',ireal,' mase_local = ',err_local/npart
    endif
    err_sum = err_sum + err_local
    ref_sum = ref_sum + ref_local
 enddo

 if (ref_sum > tiny(0.)) then
    mase = err_sum/(nrealisations*npart_total)
 else
    mase = 0.
 endif

 mase_tol = 8.5e-4
 nfailed = 0
 call checkval(mase,0.,mase_tol,nfailed(1),'MASE '//trim(label))
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_sphere

!-----------------------------------------------------------------------
!+
! Generate plot data showing the perf and accuracy of self-gravity solver
! in Phantom (see Bernard et al. 2026)
!+
!-----------------------------------------------------------------------
subroutine plot_SFMM()
 use setplummer, only:iprofile_plummer
 integer :: ntarg(7),i

 ntarg = (/1000,3000,10000,30000,100000,300000,1000000/)

 if (id==master) write(*,*) '--> Plot routine : Plummer sphere tests with different Npart'
 do i=1,size(ntarg)
    if (id==master) write(*,*) 'Test with Npart = ',ntarg(i)
    call get_plummer_prec_perf(ntarg(i),iprofile_plummer)
 enddo

end subroutine plot_SFMM

!-----------------------------------------------------------------------
!+
! test function to compare FMM,SFMM to directsum for Plummer and Homo
! sphere. It tests precision and perf for different theta max and Npart
!+
!-----------------------------------------------------------------------
subroutine get_plummer_prec_perf(npart_target,iprofile)
 use dim,         only:maxp
 use deriv,       only:get_derivs_global
 use eos,         only:gamma,polyk
 use mpiutils,    only:reduceall_mpi
 use options,     only:ieos,alpha,alphau,alphaB,tolh
 use part,        only:init_part,npart,xyzh,fxyzu,hfact,&
                       npartoftype,massoftype,istar,maxphase,iphase,isetphase,rhoh
 use setup_params,only:npart_total
 use testutils,   only:checkval,update_test_scores
 use setplummer,  only:get_accel_profile,profile_label,radius_from_mass,density_profile
 use spherical,   only:set_sphere,iseed_mc
 use random,      only:ran2
 use kdtree,      only:tree_accuracy
 use kernel,      only:hfact_default
 use table_utils, only:linspace
 use mpidomain,   only:i_belong
 use io,          only:id,master,iverbose
 use neighkdtree, only:use_dualtree
 use timing,      only:get_timings
 use sortutils,   only:indexx
 integer, intent(in)    :: iprofile,npart_target
 integer :: i,it,itest
 integer, parameter :: niter=10
 real, allocatable :: fxyz_dir(:,:),err_rel(:)
 integer, allocatable :: erridx(:)
 real :: rsoft,mass_total,cut_fraction,rmin,rmax,psep,theta_crit
 character(len=64) :: label,filename_max,type
 integer, parameter :: ntab = 1000
 real :: rgrid(ntab),rhotab(ntab)
 real :: maxerr(3,niter+1),minerr(3,niter+1),meanerr(3,niter+1)
 integer :: iunit
 real(kind=4) :: t1,t2,tcpu1,tcpu2,timings(3,niter+1)

 label = profile_label(iprofile)

 write(filename_max,'("data_plot_plummer_sphere_",i8.8,".ev")') npart_target
 filename_max = adjustl(filename_max)

 open(newunit=iunit,file=trim(filename_max),action='write',status='replace')
 write(iunit,"(a)") '# \theta, emax_SFMM, emax_FMM, emin_SFMM, emin_FMM, &
 &tcpu_SFMM, tcpu_FMM, tcpu_direct'

 call init_part()
 hfact = hfact_default
 gamma = 5./3.
 polyk = 0.
 ieos  = 11
 alpha  = 0.; alphau = 0.; alphaB = 0.
 tolh = 1.e-5
 rsoft = 1.0
 mass_total = 1.0

 ! construct tables for radius and density
 cut_fraction = 0.999
 rmin = 0.
 rmax = radius_from_mass(iprofile,cut_fraction,rsoft)
 call linspace(rgrid,0.,rmax)

 do i=1,ntab
    rhotab(i) = density_profile(iprofile,rgrid(i),rsoft,mass_total)
 enddo

 psep = rmax/real(ntab) ! this is not used for random placement anyway
 iverbose = 1

 iseed_mc = 1
 npart = 0
 npart_total = 0

 call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart, &
                  xyzh,npart_total,rhotab=rhotab,rtab=rgrid,exactN=.true.,&
                  np_requested=npart_target,mask=i_belong,verbose=.false.)
 !call set_sphere('random',id,master,rmin,rmax,psep,hfact,npart,xyzh,npart_total,np_requested=npart_target)

 massoftype(istar) = mass_total/real(npart_total)
 npartoftype(istar) = npart

 if (maxphase==maxp) then
    iphase(1:npart) = isetphase(istar,iactive=.true.)
 endif

 allocate(fxyz_dir(3,npart))
 allocate(err_rel(npart))
 allocate(erridx(npart))

 call get_derivs_global(icall=1)

 tree_acc: do it=0,niter
    theta_crit = 0.1 + it*0.05
    do itest=3,1,-1
       if (itest==1) then
          type = "SFMM"
          use_dualtree = .true.
          tree_accuracy = theta_crit
       elseif (itest==2) then
          type = "FMM"
          use_dualtree = .false.
          tree_accuracy = theta_crit
       else
          type = "direct"
          use_dualtree = .false.
          tree_accuracy = 0.
       endif

       if (itest==3 .and. it>0) cycle
       call get_timings(t1,tcpu1)
       call get_derivs_global(icall=2)
       call get_timings(t2,tcpu2)

       if (itest==3) fxyz_dir = fxyzu(0:3,1:npart)

       err_rel = norm2(fxyzu(1:3,1:npart)-fxyz_dir,1)/norm2(fxyz_dir,1)
       call indexx(npart, err_rel, erridx)

       maxerr(itest,it+1)  = maxval(err_rel)
       minerr(itest,it+1)  = err_rel(erridx(npart/10))
       meanerr(itest,it+1) = sum(err_rel)/npart

       if (itest==3) then
          timings(itest,1:niter+1) = tcpu2-tcpu1
       else
          timings(itest,it+1) = tcpu2-tcpu1
       endif
    enddo
 enddo tree_acc

 do it=0,niter
    write(iunit,*) 0.1 + it*0.05, maxerr(1:2,it+1), minerr(1:2,it+1), meanerr(1:2,it+1), timings(1:3,it+1)
 enddo

 close(iunit)

 use_dualtree = .true.

end subroutine get_plummer_prec_perf

!-----------------------------------------------------------------------
!+
!   Copy gas particles to sinks
!+
!-----------------------------------------------------------------------
subroutine copy_gas_particles_to_sinks(npart,nptmass,xyzh,xyzmh_ptmass,massi)
 integer, intent(in)  :: npart
 integer, intent(out) :: nptmass
 real, intent(in)  :: xyzh(:,:),massi
 real, intent(out) :: xyzmh_ptmass(:,:)
 integer :: i

 nptmass = npart
 do i=1,npart
    ! make a sink particle with the position of each SPH particle
    xyzmh_ptmass(1:3,i) = xyzh(1:3,i)
    xyzmh_ptmass(4,i) =  massi ! same mass as SPH particles
    xyzmh_ptmass(5:,i) = 0.
 enddo

end subroutine copy_gas_particles_to_sinks

!-----------------------------------------------------------------------
!+
!   Copy half of the gas particles to sinks
!+
!-----------------------------------------------------------------------
subroutine copy_half_gas_particles_to_sinks(npart,nptmass,xyzh,xyzmh_ptmass,massi,hi)
 use io,       only: id,master,fatal
 use mpiutils, only:bcast_mpi
 integer, intent(inout) :: npart
 integer, intent(out)   :: nptmass
 real, intent(in)  :: xyzh(:,:),massi,hi
 real, intent(out) :: xyzmh_ptmass(:,:)
 integer :: i, nparthalf

 nptmass = 0
 nparthalf = npart/2

 call bcast_mpi(nparthalf)

 if (id==master) then
    ! Assuming all gas particles are already on master,
    ! create sinks here and send them to other tasks

    ! remove half the particles by changing npart
    npart = nparthalf

    do i=npart+1,2*npart
       nptmass = nptmass + 1
       call bcast_mpi(nptmass)
       ! make a sink particle with the position of each SPH particle
       xyzmh_ptmass(1:3,nptmass) = xyzh(1:3,i)
       xyzmh_ptmass(4,nptmass)  =  massi ! same mass as SPH particles
       xyzmh_ptmass(5:,nptmass) = hi
       xyzmh_ptmass(6,nptmass)  = hi
       call bcast_mpi(xyzmh_ptmass(1:6,nptmass))
    enddo
 else
    ! Assuming there are no gas particles here,
    ! get sinks from master

    if (npart /= 0) call fatal("copy_half_gas_particles_to_sinks","there are particles on a non-master task")

    ! Get nparthalf from master, but don't change npart from zero
    do i=nparthalf+1,2*nparthalf
       call bcast_mpi(nptmass)
       call bcast_mpi(xyzmh_ptmass(1:6,nptmass))
    enddo
 endif

end subroutine copy_half_gas_particles_to_sinks

!-----------------------------------------------------------------------
!+
!   Get the distance between two points
!+
!-----------------------------------------------------------------------
subroutine get_dx_dr(x1,x2,dx,dr)
 real, intent(in) :: x1(3),x2(3)
 real, intent(out) :: dx(3),dr

 dx = x1 - x2
 dr = 1./sqrt(dot_product(dx,dx))

end subroutine get_dx_dr

end module testgravity
