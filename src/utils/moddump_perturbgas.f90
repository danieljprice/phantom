!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Give velocity perturbation to gas particles
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - perturb_factor      : *fractional velocity perturbation for gas*
!   - perturb_sink        : *also perturb sink particle velocities*
!   - sink_perturb_factor : *fractional velocity perturbation for sinks*
!
! :Dependencies: infile_utils, moddump_utils, part, prompting
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters
 real    :: perturb_factor      = 0.5      ! fractional velocity perturbation for gas
 real    :: sink_perturb_factor = 0.5      ! fractional velocity perturbation for sinks
 logical :: perturb_sink        = .false.  ! also perturb sink particle velocities

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:vxyz_ptmass,nptmass

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i

 if (prompt_for_params) call read_interactive_moddumpfile()

 print "(a,g0)", '>>> perturbing gas with factor ',perturb_factor
 do i=1,npart
    vxyzu(1:3,i) = (1. + perturb_factor) * vxyzu(1:3,i)
 enddo

 if (perturb_sink) then
    print "(a,g0)", '>>> perturbing sink particle with factor ',sink_perturb_factor
    do i=1,nptmass
       vxyz_ptmass(1:3,i) = (1. + sink_perturb_factor) * vxyz_ptmass(1:3,i)
    enddo
 endif

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter fractional velocity perturbation for gas',perturb_factor)
 call prompt('Also perturb sink particle velocities?',perturb_sink)
 if (perturb_sink) call prompt('Enter fractional velocity perturbation for sinks',sink_perturb_factor)

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
 call read_inopt(perturb_factor,'perturb_factor',db,errcount=nerr)
 call read_inopt(perturb_sink,'perturb_sink',db,errcount=nerr)
 if (perturb_sink) call read_inopt(sink_perturb_factor,'sink_perturb_factor',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# gas velocity perturbation parameters'
 call write_inopt(perturb_factor,'perturb_factor','fractional velocity perturbation for gas',iunit)
 call write_inopt(perturb_sink,'perturb_sink','also perturb sink particle velocities',iunit)
 if (perturb_sink) call write_inopt(sink_perturb_factor,'sink_perturb_factor','fractional velocity perturbation for sinks',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
