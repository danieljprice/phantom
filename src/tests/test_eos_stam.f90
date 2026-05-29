!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module test_eos_stam
!
! Unit tests of the Stamatellos/Lombardi equation of state
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: eos, eos_stamatellos, physcon, testutils, units
!

 use testutils,     only:checkval,update_test_scores,checkvalbuf,checkvalbuf_end,checkvalbuf_start
 use eos_stamatellos, only:read_optab,getintenerg_opdep,getopac_opdep,eos_file
 use eos_stamatellos, only:optable
 implicit none
 public :: run_test_stam

 integer, parameter :: nrho = 100, nT = 100

 private

contains

subroutine run_test_stam(ntests,npass)
 integer, intent(inout) :: ntests,npass
 integer :: ierr,nfail(2)
 logical :: got_phantom_dir
 character(len=20) :: pdir

 call get_environment_variable('PHANTOM_DIR',pdir)
 got_phantom_dir = (len_trim(pdir) > 0)
 if (.not. got_phantom_dir) return
 call read_optab(eos_file,ierr)
 nfail(:) = 0
 if (ierr  /=  0) then
    nfail(1) = 1
 else
    call test_interp_optab(nfail(1),npass)
 endif
 ntests = ntests + 1
 call update_test_scores(ntests,nfail(:),npass)
 call finish_test_stam
end subroutine run_test_stam

subroutine test_interp_optab(nfail,npass)
 use units,  only:unit_density,unit_ergg,unit_velocity
 use eos, only:equationofstate
 use physcon, only:kb_on_mh
 integer, intent(out) :: nfail,npass
 real :: logrhomin,logrhomax,logtmin,logtmax,tol,errmax
 real :: dlogtemp,dlogrho,rhoi,Ti,ui,Tref,spsoundi
 integer  :: irho,itemp,ndiff,ncheck
 real :: ponrhoi,xi,yi,zi,spsoundrefi,gammai
 real :: kappaBar,kappaPart,mui
 real :: rhoi_cgs, ui_cgs
 character(len=30) :: label
 label = 'eos_stamatellos interpolation'
 xi = 0.; yi = 0.; zi = 0. !These aren't used but needed for eos call
 ndiff = 0; ncheck = 0; nfail = 0
 tol = 0.001 !tolerance
 call checkvalbuf_start(label)
 ! create temperature array and rho array
 !1e-24 to 1 g cm^-3
 logtmin = log10(4.); logtmax = 4. !K to 10^4 K
 logrhomax = -0.045; logrhomin = -24.
 dlogtemp = (logtmax-logtmin)/nT
 dlogrho = (logrhomax - logrhomin)/nrho
 do irho=1,nrho
    rhoi = 10.**(logrhomin+(irho*dlogrho))
    do itemp=1,nT
       Ti = 10**(logtmin+(itemp*dlogtemp))
       !print *, irho,itemp, rhoi, Ti
       call getintenerg_opdep(Ti,rhoi,ui)
       Tref = Ti
       rhoi_cgs = rhoi/real(unit_density) ! to prevent compiler errors in equationofstate
       ui_cgs = ui/real(unit_ergg)
       call equationofstate(24,ponrhoi,spsoundi,rhoi_cgs,xi,yi,zi,Ti,ui_cgs) !this calls getopac...
       call getopac_opdep(ui,rhoi,kappaBar,kappaPart,Ti,mui)
       gammai = 1. + (ponrhoi/(ui/unit_ergg))
       spsoundrefi = sqrt(gammai*kb_on_mh*Tref/mui)
       spsoundi = spsoundi * unit_velocity
       call checkvalbuf(spsoundi,spsoundrefi,tol,label,ndiff,ncheck,errmax,use_rel_tol=.true.)
    enddo
 enddo
 call checkvalbuf_end(label,ncheck,ndiff,errmax,tol)
 if (ndiff >0) then
    nfail = nfail + 1
 else
    npass = npass + 1
 endif
end subroutine test_interp_optab

subroutine finish_test_stam
 if (allocated(optable)) deallocate(optable)
end subroutine finish_test_stam

end module test_eos_stam
