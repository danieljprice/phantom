!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Change the units of a dumpfile.
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - change_udist : *adjust the length unit*
!   - change_umass : *adjust the mass unit*
!   - change_utime : *adjust the time unit*
!   - udist_factor : *factor to scale the length unit by (if udist_fixed=T)*
!   - udist_fixed  : *adjust the length unit by a fixed factor*
!   - umass_factor : *factor to scale the mass unit by (if umass_fixed=T)*
!   - umass_fixed  : *adjust the mass unit by a fixed factor*
!   - utime_factor : *factor to scale the time unit by (if utime_fixed=T)*
!   - utime_fixed  : *adjust the time unit by a fixed factor*
!
! :Dependencies: infile_utils, moddump_utils, prompting, units
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 logical :: change_utime = .false., utime_fixed = .false.
 logical :: change_udist = .false., udist_fixed = .false.
 logical :: change_umass = .false., umass_fixed = .false.
 real    :: utime_factor = 1.0, udist_factor = 1.0, umass_factor = 1.0

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use units, only:umass,udist,utime

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: fixed_tot,total_units_to_change
 real    :: umass_tmp,utime_tmp,udist_tmp,grav_const

 grav_const = udist**3/(utime**2*umass)

 print*,'Current time unit is ',utime,'.'
 print*,'Current length unit is ',udist,'.'
 print*,'Current mass unit is ',umass,'.'
 !
 ! prompt only when the driver did not find a parameter file
 !
 if (prompt_for_params) call read_interactive_moddumpfile()

 ! Cannot fix all three units simultaneously (mass is the one dropped)
 fixed_tot = 0
 if (utime_fixed) fixed_tot = fixed_tot + 1
 if (udist_fixed) fixed_tot = fixed_tot + 1
 if (umass_fixed .and. fixed_tot >= 2) then
    print*,'Cannot change the mass unit by a fixed value since two other units are being adjusted by fixed value.'
    umass_fixed = .false.
 endif

 ! Check to make sure that more than one unit is being changed
 total_units_to_change = 0
 if (change_utime) total_units_to_change = total_units_to_change + 1
 if (change_udist) total_units_to_change = total_units_to_change + 1
 if (change_umass) total_units_to_change = total_units_to_change + 1

 if (total_units_to_change<=1) then
    print*,'You must allow more than one unit to be changed'
    stop
 endif

 ! Begin writing tmp code units
 utime_tmp = utime
 udist_tmp = udist
 umass_tmp = umass
 if (change_utime .and. utime_fixed) utime_tmp = utime_factor*utime
 if (change_udist .and. udist_fixed) udist_tmp = udist_factor*udist
 if (change_umass .and. umass_fixed) umass_tmp = umass_factor*umass

 print*,'Temporary time unit is ',utime_tmp,'.'
 print*,'Temporary length unit is ',udist_tmp,'.'
 print*,'Temporary mass unit is ',umass_tmp,'.'

 ! Begin changing code units that are allowed to vary, but are not fixed
 if ((.not.utime_fixed) .AND. change_utime) then
    utime_tmp = (udist_tmp**3/(grav_const*umass_tmp))**(1.0/2.0)
 endif

 if ((.not.udist_fixed) .AND. change_udist) then
    udist_tmp = (grav_const*utime_tmp**2*umass_tmp)**(1.0/3.0)
 endif

 if ((.not.umass_fixed) .AND. change_umass) then
    umass_tmp = udist_tmp**3/(grav_const*utime_tmp**2)
 endif

 print*,'Temporary time unit is ',utime_tmp,'.'
 print*,'Temporary length unit is ',udist_tmp,'.'
 print*,'Temporary mass unit is ',umass_tmp,'.'

 ! Check that everything has been changed properly
 print*,'Gravitational constant is ',grav_const,'.'
 print*,'With new units, it is now ', udist_tmp**3/(utime_tmp**2*umass_tmp),'.'

 ! Write new code units to header
 umass = umass_tmp
 udist = udist_tmp
 utime = utime_tmp

end subroutine modify_dump

!
!---Interactively set the moddump parameters--------------------------------
!
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt
 integer :: fixed_tot

 call prompt('Would you like time unit to be adjusted',change_utime)
 if (change_utime) then
    call prompt('Would you like time unit to be adjusted by a fixed value',utime_fixed)
    if (utime_fixed) call prompt('Enter in value you want to scale time unit by',utime_factor,0.)
 endif

 call prompt('Would you like length unit to be adjusted',change_udist)
 if (change_udist) then
    call prompt('Would you like length unit to be adjusted by a fixed value',udist_fixed)
    if (udist_fixed) call prompt('Enter in value you want to scale length unit by',udist_factor,0.)
 endif

 fixed_tot = 0
 if (utime_fixed) fixed_tot = fixed_tot + 1
 if (udist_fixed) fixed_tot = fixed_tot + 1

 call prompt('Would you like mass unit to be adjusted',change_umass)
 if (change_umass) then
    if (fixed_tot<2) then
       call prompt('Would you like mass unit to be adjusted by a fixed value',umass_fixed)
    else
       print*,'Cannot change the mass unit by a fixed value since two other units are being adjusted by fixed value.'
    endif
    if (umass_fixed) call prompt('Enter in value you want to scale mass unit by',umass_factor,0.)
 endif

end subroutine read_interactive_moddumpfile

!
!---Write the moddump parameter file----------------------------------------
!
subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 call write_inopt(change_utime,'change_utime','adjust the time unit',iunit)
 call write_inopt(utime_fixed,'utime_fixed','adjust the time unit by a fixed factor',iunit)
 call write_inopt(utime_factor,'utime_factor','factor to scale the time unit by (if utime_fixed=T)',iunit)
 call write_inopt(change_udist,'change_udist','adjust the length unit',iunit)
 call write_inopt(udist_fixed,'udist_fixed','adjust the length unit by a fixed factor',iunit)
 call write_inopt(udist_factor,'udist_factor','factor to scale the length unit by (if udist_fixed=T)',iunit)
 call write_inopt(change_umass,'change_umass','adjust the mass unit',iunit)
 call write_inopt(umass_fixed,'umass_fixed','adjust the mass unit by a fixed factor',iunit)
 call write_inopt(umass_factor,'umass_factor','factor to scale the mass unit by (if umass_fixed=T)',iunit)
 close(iunit)

end subroutine write_moddump

!
!---Read the moddump parameter file-----------------------------------------
!
subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(change_utime,'change_utime',db,errcount=nerr)
 call read_inopt(utime_fixed,'utime_fixed',db,errcount=nerr)
 call read_inopt(utime_factor,'utime_factor',db,min=0.,errcount=nerr)
 call read_inopt(change_udist,'change_udist',db,errcount=nerr)
 call read_inopt(udist_fixed,'udist_fixed',db,errcount=nerr)
 call read_inopt(udist_factor,'udist_factor',db,min=0.,errcount=nerr)
 call read_inopt(change_umass,'change_umass',db,errcount=nerr)
 call read_inopt(umass_fixed,'umass_fixed',db,errcount=nerr)
 call read_inopt(umass_factor,'umass_factor',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

end module moddump

