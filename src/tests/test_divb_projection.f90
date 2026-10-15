!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testdivbproj
!
! Unit tests of the adjoint divergence projection: the projection
! divergence matches the one computed in the force loop, the difference
! operator and the symmetric gradient are adjoints, and one projection
! reduces div B without increasing the discrete magnetic energy.
! Repeated with unequal particle masses when using adaptive particle
! refinement.
!
! :References:
!   Tsukamoto (2026), arXiv:2606.28197
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: apr, boundary, deriv, dim, divb_projection, eos, io,
!   mpidomain, options, part, physcon, testutils, unifdis
!
 use physcon, only:pi
 implicit none
 public :: test_divb_projection

 private

contains
!-----------------------------------------------------------------------
!+
!  Set up a random particle distribution in a periodic box with a
!  non-solenoidal field, then run the checks with equal masses and,
!  under adaptive particle refinement, with two mass levels
!+
!-----------------------------------------------------------------------
subroutine test_divb_projection(ntests,npass)
 use dim,       only:mhd,mpi,periodic,use_apr,maxphase,maxp,maxvxyzu
 use io,        only:id,master
 use boundary,  only:xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use part,      only:npart,npartoftype,massoftype,hfact,xyzh,vxyzu,Bevol,&
                     iphase,igas,isetphase,apr_level,init_part
 use unifdis,   only:set_unifdis
 use mpidomain, only:i_belong
 use apr,       only:sync_aprmassoftype
 use eos,       only:gamma,init_eos
 use options,   only:ieos
 integer, intent(inout) :: ntests,npass
 integer :: i,ierr
 real    :: psep

 if (.not. mhd .or. mpi) then
    if (id == master) write(*,"(/,a)") '--> SKIPPING DIV B PROJECTION TEST (requires MHD, no MPI)'
    return
 endif
 if (id == master) write(*,"(/,a)") '--> TESTING DIV B PROJECTION'

 call init_part()
 psep = dxbound/12.
 npart = 0
 call set_unifdis('random',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                  hfact,npart,xyzh,periodic,inputiseed=-24680,mask=i_belong)
 npartoftype(igas) = npart
 massoftype(igas)  = dxbound*dybound*dzbound/real(npart)
 call sync_aprmassoftype
 if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
 vxyzu(:,1:npart) = 0.
 if (maxvxyzu >= 4) vxyzu(4,1:npart) = 1.
 ! adiabatic, which is also compatible with radiation
 ieos  = 2
 gamma = 5./3.
 call init_eos(ieos,ierr)
 Bevol(:,1:npart) = 0.

 call run_checks('equal masses',ntests,npass)
 call check_zero_field(ntests,npass)

 if (use_apr) then
    ! two refinement levels, so neighbouring masses differ by a factor of 2
    do i=1,npart
       if (mod(i,2)==0) apr_level(i) = 2
    enddo
    call run_checks('APR, unequal masses',ntests,npass)
 endif

end subroutine test_divb_projection

!-----------------------------------------------------------------------
!+
!  Checks on the current particle distribution
!+
!-----------------------------------------------------------------------
subroutine run_checks(label,ntests,npass)
 use io,              only:id,master
 use part,            only:npart,xyzh,rho,gradh,Bevol,Bxyz,dBevol,divcurlB
 use deriv,           only:get_derivs_global
 use divb_projection, only:project_divb,get_divb,divb_adjoint_mismatch
 use testutils,       only:checkval,update_test_scores
 character(len=*), intent(in)    :: label
 integer,          intent(inout) :: ntests,npass
 integer :: i,ierr,nfailed(4)
 real    :: mismatch,divb_before,divb_after,emag_ratio,errmax
 real, allocatable :: divb(:),xvec(:,:),piscal(:)

 if (id == master) write(*,"(/,a)") '--> checking projection with '//trim(label)
 allocate(divb(npart),xvec(3,npart),piscal(npart))

 ! smoothing lengths and densities, then B/rho from the target field
 call get_derivs_global()
 do i=1,npart
    Bxyz(1:3,i) = bfunc(xyzh(1:3,i))
    Bevol(1:3,i) = Bxyz(1:3,i)/rho(i)
    Bevol(4,i) = 0.
 enddo
 call get_derivs_global()
 ! B as seen by the force loop, with the densities it used
 do i=1,npart
    Bxyz(1:3,i) = Bevol(1:3,i)*rho(i)
 enddo

 ! projection divergence against the difference divergence from force
 nfailed = 0
 call get_divb(npart,xyzh,rho,gradh,Bxyz,divb,ierr)
 errmax = maxval(abs(divb(1:npart) - divcurlB(1,1:npart)))/maxval(abs(divb(1:npart)))
 call checkval(ierr,0,0,nfailed(1),'get_divb ierr')
 call checkval(errmax,0.,1.e-5,nfailed(2),'D B matches force-loop div B')
 call update_test_scores(ntests,nfailed(1:2),npass)

 ! adjoint identity for smooth test fields
 nfailed = 0
 do i=1,npart
    xvec(1:3,i) = bfunc(xyzh(1:3,i)) + 0.3*xyzh(3,i)
    piscal(i) = sin(2.*pi*(xyzh(1,i) + xyzh(2,i)))
 enddo
 call divb_adjoint_mismatch(npart,xyzh,rho,gradh,xvec,piscal,mismatch,ierr)
 if (id == master) write(*,"(a,es10.3)") ' adjoint mismatch = ',mismatch
 call checkval(ierr,0,0,nfailed(1),'adjoint ierr')
 call checkval(mismatch,0.,1.e-10,nfailed(2),'D and G are adjoints')
 call update_test_scores(ntests,nfailed(1:2),npass)

 ! one projection of the non-solenoidal field
 nfailed = 0
 Bevol(4,1:npart)  = 1.
 dBevol(4,1:npart) = 1.
 call project_divb(npart,xyzh,rho,gradh,Bevol,Bxyz,dBevol,divcurlB,1.e-8,1000,1,ierr,&
                   divb_before,divb_after,emag_ratio)
 if (id == master) write(*,"(a,es10.3,a,es10.3,a,es10.3)") &
    ' |div B|_V ',divb_before,' -> ',divb_after,', E_B ratio ',emag_ratio
 call checkval(ierr,0,0,nfailed(1),'project_divb ierr')
 call checkval(divb_after/divb_before,0.,1.e-3,nfailed(2),'div B reduced')
 call checkval(max(emag_ratio - 1.,0.),0.,1.e-10,nfailed(3),'E_B not increased')
 call checkval(maxval(abs(Bevol(4,1:npart))) + maxval(abs(dBevol(4,1:npart))),0.,0.,nfailed(4),'psi and dpsi/dt reset')
 call update_test_scores(ntests,nfailed,npass)

 ! Bxyz holds the projected field, consistent with B/rho
 nfailed = 0
 errmax = 0.
 do i=1,npart
    errmax = max(errmax,maxval(abs(Bxyz(1:3,i) - Bevol(1:3,i)*rho(i))))
 enddo
 call checkval(errmax,0.,1.e-12,nfailed(1),'Bxyz = rho*Bevol after projection')
 call update_test_scores(ntests,nfailed(1:1),npass)

end subroutine run_checks

!-----------------------------------------------------------------------
!+
!  Field that vanishes in part of the box: the practical stop on
!  h|div B|/|B| must still fire, rather than running to the cap
!+
!-----------------------------------------------------------------------
subroutine check_zero_field(ntests,npass)
 use io,              only:id,master
 use boundary,        only:xmin,dxbound
 use part,            only:npart,xyzh,rho,gradh,Bevol,Bxyz,dBevol,divcurlB
 use deriv,           only:get_derivs_global
 use divb_projection, only:project_divb
 use testutils,       only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer, parameter :: maxit = 1000
 integer :: i,ierr,nits,nfailed(2)
 real    :: divb_before,divb_after,emag_ratio

 if (id == master) write(*,"(/,a)") '--> checking projection with a zero-field region'
 do i=1,npart
    Bxyz(1:3,i) = bfunc(xyzh(1:3,i))
    if (xyzh(1,i) < xmin + 0.25*dxbound) Bxyz(1:3,i) = 0.
    Bevol(1:3,i) = Bxyz(1:3,i)/rho(i)
    Bevol(4,i) = 0.
 enddo
 call get_derivs_global()

 nfailed = 0
 ! zero relative tolerance, so only the practical stop can end the solve
 call project_divb(npart,xyzh,rho,gradh,Bevol,Bxyz,dBevol,divcurlB,0.,maxit,1,ierr,&
                   divb_before,divb_after,emag_ratio,nits)
 call checkval(ierr,0,0,nfailed(1),'project_divb ierr')
 call checkval(nits < maxit,.true.,nfailed(2),'practical stop before the cap')
 call update_test_scores(ntests,nfailed,npass)

end subroutine check_zero_field

!-----------------------------------------------------------------------
!+
!  Smooth periodic test field with non-zero divergence
!+
!-----------------------------------------------------------------------
function bfunc(x) result(b)
 use boundary, only:xmin,ymin,zmin,dxbound,dybound,dzbound
 real, intent(in) :: x(3)
 real :: b(3)

 b(1) = 1. + 0.1*sin(2.*pi*(x(1) - xmin)/dxbound)
 b(2) = 0.5 + 0.2*cos(2.*pi*(x(2) - ymin)/dybound)
 b(3) = 0.1*sin(2.*pi*(x(3) - zmin)/dzbound)

end function bfunc

end module testdivbproj
