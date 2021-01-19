!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module checkconserved
!
! Utility routines to perform runtime checks that
! conservation laws are appropriately satisfied
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, externalforces, io, options, part
!
 use dim, only:maxdusttypes
 implicit none
 real, public :: get_conserv = 1.0 ! to track when we have initial values for conservation laws
 real, public :: etot_in,angtot_in,totmom_in,mdust_in(maxdusttypes)

 public :: init_conservation_checks, check_conservation_error

 private

contains

!----------------------------------------------------------------
!+
!  check if conservation of various properties *should* be
!  possible given the range of physics selected
!+
!----------------------------------------------------------------
subroutine init_conservation_checks(should_conserve_energy,should_conserve_momentum,&
                                    should_conserve_angmom,should_conserve_dustmass)
 use options, only:icooling,ieos,ipdv_heating,ishock_heating,&
                   iresistive_heating,use_dustfrac,iexternalforce
 use dim,     only:mhd,maxvxyzu,periodic,particles_are_injected
 use part,    only:iboundary,npartoftype
 logical, intent(out) :: should_conserve_energy,should_conserve_momentum
 logical, intent(out) :: should_conserve_angmom,should_conserve_dustmass

 !
 ! should conserve energy if using adiabatic equation of state with no cooling
 ! as long as all heating terms are included
 !
 should_conserve_energy = (maxvxyzu==4 .and. ieos==2 .and. &
                          icooling==0 .and. ipdv_heating==1 .and. ishock_heating==1 &
                          .and. (.not.mhd .or. iresistive_heating==1))
 !
 ! code should conserve momentum unless boundary particles are employed
 !
 if (iexternalforce/=0) then
    should_conserve_momentum = .false.
 else
    should_conserve_momentum = (npartoftype(iboundary)==0)
 endif
 !
 ! code should conserve angular momentum as long as no boundaries (fixed or periodic)
 ! and as long as there are no non-radial forces (iexternalforce > 1)
 !
 should_conserve_angmom = (npartoftype(iboundary)==0 .and. .not.periodic &
                          .and. iexternalforce <= 1)
 !
 ! should always conserve dust mass
 !
 should_conserve_dustmass = use_dustfrac
 !
 ! Each injection routine will need to bookeep conserved quantities, but until then...
 !
 if (particles_are_injected) then
    should_conserve_energy   = .false.
    should_conserve_momentum = .false.
    should_conserve_angmom   = .false.
 endif

end subroutine init_conservation_checks

!----------------------------------------------------------------
!+
!  routine to check conservation errors during the calculation
!  and stop if it is too large
!+
!----------------------------------------------------------------
subroutine check_conservation_error(val,ref,tol,label,decrease)
 use io,             only:error,fatal,iverbose
 use options,        only:iexternalforce
 use externalforces, only:iext_corot_binary
 real, intent(in) :: val,ref,tol
 character(len=*), intent(in) :: label
 logical, intent(in), optional :: decrease
 real :: err
 character(len=20) :: string

 if (abs(ref) > 1.e-3) then
    err = (val - ref)/abs(ref)
 else
    err = (val - ref)
 endif
 if (present(decrease)) then
    err = max(err,0.) ! allow decrease but not increase
 else
    err = abs(err)
 endif
 if (err > tol) then
    if ((trim(label) == 'angular momentum' .or. trim(label) == 'energy') &
        .and. iexternalforce == iext_corot_binary) then
       call error('evolve',trim(label)//' is not being conserved due to corotating frame',var='err',val=err)
    else
       call error('evolve','Large error in '//trim(label)//' conservation ',var='err',val=err)
       call get_environment_variable('I_WILL_NOT_PUBLISH_CRAP',string)
       if (.not. (trim(string)=='yes')) then
          print "(2(/,a))",' You can ignore this error and continue by setting the ',&
                           ' environment variable I_WILL_NOT_PUBLISH_CRAP=yes to continue'
          call fatal('evolve',' Conservation errors too large to continue simulation')
       endif
    endif
 else
    if (iverbose >= 2) print "(a,es10.3)",trim(label)//' error is ',err
 endif

end subroutine check_conservation_error

end module checkconserved
