!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
! :Dependencies: deriv, dim, directsum, energies, eos, io, kdtree, options,
!   part, physcon, spherical, testutils, timing
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
 use dim,       only:maxp
 use part,      only:init_part,npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu, &
                     gradh,poten,iphase,isetphase,maxphase,labeltype
 use eos,       only:polyk,gamma
 use options,   only:ieos,alpha,alphau,alphaB,tolh
 use spherical, only:set_sphere
 use deriv,     only:get_derivs_global
 use physcon,   only:pi
 use timing,    only:getused,printused
 use directsum, only:directsum_grav
 use energies,  only:compute_energies,epot
 use kdtree,    only:tree_accuracy
 use testutils, only:checkval,checkvalbuf_end,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(18)
 integer :: maxvxyzu,nx,np,i,k
 real :: psep,totvol,totmass,rhozero
 real :: time,rmin,rmax,phitot
 real(kind=4) :: t1,t2
 real :: epoti,tree_acc_prev

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
       maxvxyzu = size(vxyzu(:,1))
       totvol   = 4./3.*pi*rmax**3
       nx       = int(np**(1./3.))
       psep     = totvol**(1./3.)/real(nx)
       !print*,' got psep = ',nx,psep
       psep     = 0.18
       npart    = 0
       call set_sphere('cubic',id,master,rmin,rmax,psep,hfact,npart,xyzh)
       !print*,' using npart = ',npart
       np       = npart
       !iverbose = 5
!
!--set particle properties
!
       totmass        = 1.
       rhozero        = totmass/totvol
       npartoftype(:) = 0
       npartoftype(k) = npart
       massoftype(:)  = 0.0
       massoftype(k)  = totmass/npart
       if (maxphase==maxp) then
          do i=1,npart
             iphase(i) = isetphase(k,iactive=.true.)
          enddo
       endif
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
!
!--calculate derivatives
!
       call getused(t1)
       call get_derivs_global()
       call getused(t2)
       if (id==master) call printused(t1)
!
!--compute gravitational forces by direct summation
!
       call directsum_grav(xyzh,gradh,vxyzu,phitot,npart)
!
!--compare the results
!
       call checkval(np,fxyzu(1,:),vxyzu(1,:),5.e-3,nfailed(1),'fgrav(x)')
       call checkval(np,fxyzu(2,:),vxyzu(2,:),6.e-3,nfailed(2),'fgrav(y)')
       call checkval(np,fxyzu(3,:),vxyzu(3,:),9.4e-3,nfailed(3),'fgrav(z)')
       epoti = 0.
       do i=1,npart
          epoti = epoti + poten(i)
       enddo
       call checkval(epoti,phitot,5.1e-4,nfailed(4),'potential')
       call checkval(epoti,-3./5.*totmass**2/rmax,3.6e-2,nfailed(5),'potential=-3/5 GMM/R')
       ! check that potential energy computed via compute_energies is also correct
       call compute_energies(0.)
       call checkval(epot,phitot,5.1e-4,nfailed(6),'epot in compute_energies')
       call update_test_scores(ntests,nfailed(1:6),npass)
    endif
 enddo
!
!--clean up doggie-doos
!
 npartoftype(:) = 0
 massoftype(:) = 0.
 tree_accuracy = tree_acc_prev
 fxyzu = 0.
 vxyzu = 0.

end subroutine test_directsum

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
