!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Changes particle mass
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - disc_mass : *desired total disc mass [code units]*
!
! :Dependencies: infile_utils, moddump_utils, part, prompting, units
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter (written to / read from the prefix.moddump file)
 real :: disc_mass = 0.05   ! desired total disc mass [code units]

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:igas,isdead_or_accreted,kill_particle,shuffle_part
 use units,         only:umass

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real     :: current_disc_mass, mass_factor
 integer  :: i

 ! Remove particles that are dead or accreted
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       call kill_particle(i)
    endif
 enddo

 call shuffle_part(npart)
 npartoftype(igas) = npart

 if (prompt_for_params) call read_interactive_moddumpfile()

 current_disc_mass = npartoftype(igas)*massoftype(igas)
 mass_factor = disc_mass/current_disc_mass

 massoftype(igas) = mass_factor*massoftype(igas)
 print*,'Particle mass is now ', massoftype(igas)*umass, ' g'
 print*,'Total disc mass is now ', npartoftype(igas)*massoftype(igas)*umass, ' g'
 print*,'Total disc mass is now ', npartoftype(igas)*massoftype(igas), 'Msun'

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter desired total disc mass in code units',disc_mass,0.)

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
 call read_inopt(disc_mass,'disc_mass',db,errcount=nerr,min=0.)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# changemass parameters'
 call write_inopt(disc_mass,'disc_mass','desired total disc mass [code units]',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
