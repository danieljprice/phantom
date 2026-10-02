!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Changes the dustfrac around sink particles, prevents spurious shadows in RT
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - outer_radius : *radius within which to taper the dust fraction [code units]*
!
! :Dependencies: dim, infile_utils, moddump_utils, part, prompting
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter (written to / read from the prefix.moddump file)
 real :: outer_radius = 10.   ! radius within which to taper the dust fraction [code units]

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,           only:use_dust
 use part,          only:igas,ndusttypes,dustfrac

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,np_gas
 real    :: dust_to_gas,r_g

 if (.not. use_dust) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes'
    stop
 endif

 dust_to_gas = 0.01

 if (prompt_for_params) call read_interactive_moddumpfile()

 !- grainsize and graindens already set if convert from one fluid to two fluid with growth
 np_gas = npartoftype(igas)

 do i=1,np_gas
    r_g = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
    if (r_g < outer_radius) then
       r_g = r_g/outer_radius
       dustfrac(1:ndusttypes,i) = dust_frac_func(r_g)*dustfrac(1:ndusttypes,i)
    endif
 enddo

 ! massoftype(igas) = massoftype(igas)*(1. + dust_to_gas)
 ! npart = np_gas

end subroutine modify_dump

real function dust_frac_func(x)
 real, intent(in) :: x
 real, parameter :: pi = 4.*atan(1.)

 dust_frac_func = sin(x*pi/2)

end function dust_frac_func

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter outer radius within which to taper dustfrac (code units)',outer_radius,0.)

end subroutine read_interactive_moddumpfile

subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 23
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(outer_radius,'outer_radius',db,errcount=nerr,min=0.)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# dustsinkfrac parameters'
 call write_inopt(outer_radius,'outer_radius','radius within which to taper the dust fraction [code units]',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
