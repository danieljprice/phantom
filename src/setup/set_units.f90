!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setunits
!
! this module contains utilities for writing/reading units
! to/from .setup files
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - dist_unit : *distance unit (e.g. au,pc,kpc,0.1pc)*
!   - mass_unit : *mass unit (e.g. solarm,jupiterm,1e6*solarm)*
!
! :Dependencies: infile_utils, io, prompting, units
!
 implicit none
 public :: write_options_units
 public :: read_options_units,read_options_and_set_units
 public :: set_units_interactive

 character(len=20), public :: dist_unit='au'
 character(len=20), public :: mass_unit='solarm'

 private

contains

!-------------------------------------------------------------------
!+
!  interactively prompt for units
!+
!------------------------------------------------------------------
subroutine set_units_interactive(gr)
 use prompting, only:prompt
 use units,     only:select_unit
 logical, intent(in), optional :: gr
 real(kind=8) :: umass,udist
 logical :: nogr
 integer :: ierr

 nogr = .true.
 if (present(gr)) nogr = .not.gr

 ! units
 ierr = 1
 do while (ierr /= 0)
    call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
    call select_unit(mass_unit,umass,ierr)
    if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
 enddo
 if (nogr) then
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo
 endif

end subroutine set_units_interactive

!-------------------------------------------------------------------
!+
!  write units options to the .setup file
!+
!------------------------------------------------------------------
subroutine write_options_units(iunit,gr)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit
 logical, intent(in), optional :: gr
 logical :: nogr

 nogr = .true.
 if (present(gr)) nogr = .not.gr

 ! units
 write(iunit,"(/,a)") '# units'
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm,jupiterm,1e6*solarm)',iunit)
 if (nogr) call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au,pc,kpc,0.1pc)',iunit)

end subroutine write_options_units

!-------------------------------------------------------------------
!+
!  read units options from the .setup file
!+
!------------------------------------------------------------------
subroutine read_options_units(db,umass,udist,nerr,gr)
 use units,        only:select_unit
 use infile_utils, only:read_inopt,inopts
 use io,           only:error
 type(inopts), allocatable, intent(inout) :: db(:)
 real(kind=8),              intent(out)   :: umass,udist
 integer,                   intent(inout) :: nerr
 logical,                   intent(in), optional :: gr
 logical :: nogr
 integer :: ierr

 nogr = .true.
 if (present(gr)) nogr = .not.gr

 ! units
 call read_inopt(mass_unit,'mass_unit',db,errcount=nerr,default=trim(mass_unit))
 if (nogr) call read_inopt(dist_unit,'dist_unit',db,errcount=nerr,default=trim(dist_unit))

 !
 ! parse units
 !
 call select_unit(mass_unit,umass,ierr)
 if (ierr /= 0) then
    call error('units','mass unit not recognised')
    nerr = nerr + 1
 endif
 if (nogr) then
    call select_unit(dist_unit,udist,ierr)
    if (ierr /= 0) then
       call error('units','length unit not recognised')
       nerr = nerr + 1
    endif
 endif

end subroutine read_options_units

!-------------------------------------------------------------------
!+
!  combined read-and-set units routine
!+
!------------------------------------------------------------------
subroutine read_options_and_set_units(db,nerr,gr)
 use units,        only:set_units
 use infile_utils, only:inopts
 type(inopts), allocatable, intent(inout) :: db(:)
 integer, intent(inout) :: nerr
 logical, intent(in), optional :: gr
 real(kind=8) :: umass,udist
 logical :: mygr
 integer :: nerr_was

 mygr = .false.
 if (present(gr)) mygr = gr

 nerr_was = nerr
 call read_options_units(db,umass,udist,nerr,gr=mygr)

 if (nerr == nerr_was) then
    if (mygr) then
       call set_units(mass=umass,c=1.d0,G=1.d0) ! use geometric units for gr
    else
       call set_units(dist=udist,mass=umass,G=1.d0)
    endif
 endif

end subroutine read_options_and_set_units

end module setunits
