!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module injection
!
! wrapper module as interface to swap-in injection modules,
! contains options common to all injection modules
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - rkill : *deactivate particles outside this radius (<0 is off)*
!
! :Dependencies: dim, infile_utils, inject
!
 use dim, only:inject_parts
 implicit none

 real, public :: rkill = -1.

 public :: write_options_injection,read_options_injection

contains

!-----------------------------------------------------------------------
!+
!  Wrapper routine for writing injection options
!+
!-----------------------------------------------------------------------
subroutine write_options_injection(iunit)
 use infile_utils, only:write_inopt
 use inject,       only:write_options_inject,inject_type,update_injected_par
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for injecting/removing particles'
 if (inject_parts) call write_options_inject(iunit)
 if (inject_parts .and. inject_type=='sim') call update_injected_par()
 call write_inopt(rkill,'rkill','deactivate particles outside this radius (<0 is off)',iunit)

end subroutine write_options_injection

!-----------------------------------------------------------------------
!+
!  Wrapper routine for reading injection options
!+
!-----------------------------------------------------------------------
subroutine read_options_injection(db,nerr)
 use inject,       only:read_options_inject
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(rkill,'rkill',db,errcount=nerr,default=rkill)
 if (inject_parts) call read_options_inject(db,nerr)

end subroutine read_options_injection

end module injection
