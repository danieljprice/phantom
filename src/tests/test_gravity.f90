!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testgravity
!
! Unit tests of self-gravity
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: deriv, dim, directsum, energies, eos, io, kdtree,
!   linklist, mpibalance, mpiutils, options, part, physcon, ptmass,
!   setup_params, sort_particles, spherical, testapr, testutils, timing
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
 logical :: testdirectsum,testpolytrope,testtaylorseries,testall

 testdirectsum    = .false.
 testtaylorseries = .false.
 testpolytrope    = .false.
 testall          = .false.
 select case(string)
 case('taylorseries')
    testtaylorseries = .true.
 case('directsum')
    testdirectsum = .true.
 case('polytrope')
    testpolytrope = .true.
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
 use kdtree,    only:compute_fnode,expand_fgrav_in_taylor_series
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(18),i,npnode
 real :: xposi(3),xposj(3),x0(3),dx(3),fexact(3),f0(3)
 real :: xposjd(3,3),dfdx_approx(3,3),d2f(3,3),dpot(3)
 real :: fnode(20),quads(6)
 real :: dr,dr2,phi,phiexact,pmassi,tol,totmass

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
 call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

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
    quads(1) = quads(1) + pmassi*(3.*dx(1)*dx(1) - dr2)
    quads(2) = quads(2) + pmassi*(3.*dx(1)*dx(2))
    quads(3) = quads(3) + pmassi*(3.*dx(1)*dx(3))
    quads(4) = quads(4) + pmassi*(3.*dx(2)*dx(2) - dr2)
    quads(5) = quads(5) + pmassi*(3.*dx(2)*dx(3))
    quads(6) = quads(6) + pmassi*(3.*dx(3)*dx(3) - dr2)
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
 call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

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

 if (id==master) write(*,"(/,a)") '--> checking results of compute_fnode routine'
 !
 ! check that components of fnode are derivatives of each other
 !
 tol = 1.e-6
 call get_finite_diff(3,x0,xposj,totmass,quads,fnode,dfdx_approx,dpot,d2f,tol)
 nfailed(:) = 0
 call checkval(fnode(1),dpot(1),tol,nfailed(1),'fx=-dphi/dx')
 call checkval(fnode(2),dpot(2),tol,nfailed(2),'fy=-dphi/dy')
 call checkval(fnode(3),dpot(3),tol,nfailed(3),'fz=-dphi/dz')
 call checkval(fnode(4),dfdx_approx(1,1),tol,nfailed(4),'dfx/dx')
 call checkval(fnode(5),dfdx_approx(1,2),tol,nfailed(5),'dfx/dy')
 call checkval(fnode(6),dfdx_approx(1,3),tol,nfailed(6),'dfx/dz')
 call checkval(fnode(7),dfdx_approx(2,2),tol,nfailed(7),'dfy/dy')
 call checkval(fnode(8),dfdx_approx(2,3),tol,nfailed(8),'dfx/dz')
 call checkval(fnode(9),dfdx_approx(3,3),tol,nfailed(9),'dfz/dz')
 call checkval(fnode(10),d2f(1,1),1.e-3,nfailed(10),'d^2fx/dx^2')
 call checkval(fnode(13),d2f(1,2),1.1e-3,nfailed(11),'d^2fx/dy^2')
 call checkval(fnode(15),d2f(1,3),1.e-3,nfailed(12),'d^2fx/dz^2')
 call checkval(fnode(11),d2f(2,1),1.e-3,nfailed(13),'d^2fy/dx^2')
 call checkval(fnode(16),d2f(2,2),1.e-3,nfailed(14),'d^2fy/dy^2')
 call checkval(fnode(18),d2f(2,3),1.e-3,nfailed(15),'d^2fy/dz^2')
 call checkval(fnode(12),d2f(3,1),1.e-3,nfailed(16),'d^2fz/dx^2')
 call checkval(fnode(17),d2f(3,2),1.2e-3,nfailed(17),'d^2fz/dy^2')
 call checkval(fnode(19),d2f(3,3),1.e-3,nfailed(18),'d^2fz/dz^2')
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
 call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

 dx = xposi - x0   ! perform expansion about x0
 call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),phi)
 !print*,'           exact force = ',fexact,' phi = ',phiexact
 !print*,'       force at origin = ',fnode(1:3), ' phi = ',fnode(20)
 !print*,'force w. taylor series = ',f0, ' phi = ',phi
 nfailed(:) = 0
 call checkval(f0(1),fexact(1),4.3e-5,nfailed(1),'fx taylor series about f0')
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
 use linklist,        only:set_linklist
 use sort_particles,  only:sort_part_id
 use mpibalance,      only:balancedomains
 use testapr,         only:setup_apr_region_for_test
 use setup_params,    only:npart_total

 integer, intent(inout) :: ntests,npass
 integer :: nfailed(18),boundi,boundf
 integer :: maxvxyzu,nx,np,i,k,merge_n,merge_ij(maxptmass),nfgrav
 real :: psep,totvol,totmass,rhozero,tol,pmassi
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
!--compare the results
!
       call checkval(npart,fxyzu(1,:),fgrav(1,:),5.e-3,nfailed(1),'fgrav(x)')
       call checkval(npart,fxyzu(2,:),fgrav(2,:),6.e-3,nfailed(2),'fgrav(y)')
       call checkval(npart,fxyzu(3,:),fgrav(3,:),9.4e-3,nfailed(3),'fgrav(z)')
       deallocate(fgrav)
       epoti = 0.
       do i=1,npart
          epoti = epoti + poten(i)
       enddo
       epoti = reduceall_mpi('+',epoti)
       call checkval(epoti,phitot,5.2e-4,nfailed(4),'potential')
       call checkval(epoti,-3./5.*totmass**2/rmax,3.6e-2,nfailed(5),'potential=-3/5 GMM/R')
       ! check that potential energy computed via compute_energies is also correct
       call compute_energies(0.)
       call checkval(epot,phitot,5.2e-4,nfailed(6),'epot in compute_energies')
       call update_test_scores(ntests,nfailed(1:6),npass)
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

    if(mpi) then
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

    if(use_sinktree) then
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

subroutine copy_half_gas_particles_to_sinks(npart,nptmass,xyzh,xyzmh_ptmass,massi,hi)
 use io,       only: id,master,fatal
 use mpiutils, only: bcast_mpi
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
       xyzmh_ptmass(5:,nptmass) = 0.
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

subroutine get_dx_dr(x1,x2,dx,dr)
 real, intent(in) :: x1(3),x2(3)
 real, intent(out) :: dx(3),dr

 dx = x1 - x2
 dr = 1./sqrt(dot_product(dx,dx))

end subroutine get_dx_dr

subroutine get_finite_diff(ndim,x0,xposj,totmass,quads,fnode,dfdx,dpot,d2f,eps)
 use kdtree,    only:compute_fnode
 integer, intent(in)  :: ndim
 real,    intent(in)  :: x0(ndim),xposj(ndim),totmass,quads(6),fnode(20),eps
 real,    intent(out) :: dfdx(ndim,ndim),dpot(ndim),d2f(ndim,ndim)
 integer :: i,j
 real :: dx(ndim),x0_plus(ndim),x0_minus(ndim)
 real :: dr,fnode_plus(20),fnode_minus(20)

 do j=1,ndim
    x0_plus     = x0
    x0_plus(j)  = x0(j) + eps
    x0_minus    = x0
    x0_minus(j) = x0(j) - eps
    do i=1,ndim
       call get_dx_dr(x0_plus,xposj,dx,dr)
       fnode_plus = 0.
       call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode_plus)

       call get_dx_dr(x0_minus,xposj,dx,dr)
       fnode_minus = 0.
       call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode_minus)

       dfdx(i,j) = (fnode_plus(i) - fnode_minus(i))/(2.*eps)
       d2f(i,j) = (fnode_plus(i) - 2.*fnode(i) + fnode_minus(i))/(eps*eps)
    enddo
    dpot(j) = -(fnode_plus(20) - fnode_minus(20))/(2.*eps)
 enddo

end subroutine get_finite_diff

end module testgravity
