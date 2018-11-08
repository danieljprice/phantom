!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testgravity
!
!  DESCRIPTION:
!  Unit tests of self-gravity
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: deriv, dim, directsum, eos, io, kdtree, options, part,
!    physcon, spherical, testutils, timing
!+
!--------------------------------------------------------------------------
module testgravity
 implicit none
 public :: test_gravity

 private

contains

subroutine test_gravity(ntests,npass,string)
 use io,        only:id,master
#ifdef GRAVITY
 use dim,       only:maxp
 use part,      only:npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu,fext,Bevol,mhd, &
                     alphaind,maxalpha,dustprop,ddustprop, &
                     divcurlv,divcurlB,dBevol,gradh,poten,&
                     iphase,isetphase,maxphase,dustfrac,ddustevol,temperature,labeltype, &
                     pxyzu,dens,metrics
 use eos,       only:polyk,gamma
 use options,   only:ieos,alpha,alphau,alphaB,tolh
 use testutils, only:checkval,checkvalf,checkvalbuf_start,checkvalbuf,checkvalbuf_end
 use spherical, only:set_sphere
 use deriv,     only:derivs
 use physcon,   only:pi
 use timing,    only:getused,printused
 use directsum, only:directsum_grav
 use kdtree,    only:compute_fnode,expand_fgrav_in_taylor_series,tree_accuracy
#endif
 integer,          intent(inout) :: ntests,npass
 character(len=*), intent(in)    :: string
#ifdef GRAVITY
 integer :: nfailed(18)
 logical                :: testdirectsum,testpolytrope,testtaylorseries,testall
 integer :: maxvxyzu,nx,np,i,npnode,k
 real :: psep,totvol,totmass,rhozero
 real :: time,rmin,rmax,dtext_dum,phitot
 real(kind=4) :: t1,t2
 real :: xposi(3),xposj(3),x0(3),dx(3),fexact(3),f0(3)
 real :: xposjd(3,3),dfdx_approx(3,3),d2f(3,3),dpot(3)
 real :: fnode(20)
 real :: quads(6)
 real :: dr,phiexact,phi,dr2,pmassi,epot,tol,tree_acc_prev

 if (id==master) write(*,"(/,a,/)") '--> TESTING SELF-GRAVITY'

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
!
!--unit test of the taylor series expansion
!  about the node centre
!
 if (testtaylorseries .or. testall) then
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
    ntests = ntests + 1
    if (all(nfailed==0)) npass = npass + 1

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
    ntests = ntests + 1
    if (all(nfailed==0)) npass = npass + 1

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
    ntests = ntests + 1
    if (all(nfailed==0)) npass = npass + 1

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
    ntests = ntests + 1
    if (all(nfailed==0)) npass = npass + 1

 endif

 testsum: if (testdirectsum .or. testall) then

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
          if (mhd) Bevol(:,:) = 0.
!
!--make sure AV is off
!
          alpha  = 0.
          alphau = 0.
          alphaB = 0.
          tolh = 1.e-5
          if (maxalpha==maxp)  alphaind(:,:) = 0.

          ntests = ntests + 1
!
!--calculate derivatives
!
          call getused(t1)
          call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                      Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum,pxyzu,dens,metrics)
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
          epot = 0.
          do i=1,npart
             epot = epot + poten(i)
          enddo
          call checkval(epot,phitot,4.8e-4,nfailed(4),'potential')

          if (all(nfailed(1:4)==0)) npass = npass + 1
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

 endif testsum

 if (id==master) write(*,"(/,a)") '<-- SELF-GRAVITY TESTS COMPLETE'

#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING SELF-GRAVITY TESTS (need -DGRAVITY)'
 return
#endif

end subroutine test_gravity

#ifdef GRAVITY
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
#endif

end module testgravity
