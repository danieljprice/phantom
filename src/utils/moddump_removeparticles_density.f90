!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! None
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - rho_threshold : *delete particles with density below this [code units]*
!
! :Dependencies: infile_utils, io, moddump_utils, part, prompting
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter (written to / read from the prefix.moddump file)
 real :: rho_threshold = 5e-8   ! delete particles with density below this [code units]

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:rho,kill_particle,shuffle_part
 use io,           only:fatal
 integer, intent(inout) :: npart
 integer, dimension(:), intent(inout) :: npartoftype
 real, dimension(:), intent(inout) :: massoftype
 real, dimension(:,:), intent(inout) :: xyzh,vxyzu
 real   :: rhoi
 integer :: i, compt

 if (prompt_for_params) call read_interactive_moddumpfile()

 compt = 0
 do i=1,npart
    rhoi = rho(i)
    ! write(*,*) rhoi      ! uncomment to have an idea of the density
    if (rhoi < rho_threshold) then
       call kill_particle(i,npartoftype)
       compt = compt+1
       !xyzh(4,i) = -abs(hi)    !call kill_particle(i,npoftype)
    endif
 enddo
 write(*,*) 'Particles deleted :', compt
 call shuffle_part(npart)
 if (npart /= sum(npartoftype)) call fatal('del_dead_part_dens','particles not conserved')

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter density threshold below which to delete particles (code units)',rho_threshold,0.)

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
 call read_inopt(rho_threshold,'rho_threshold',db,errcount=nerr,min=0.)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# remove-low-density-particles parameters'
 call write_inopt(rho_threshold,'rho_threshold','delete particles with density below this [code units]',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
