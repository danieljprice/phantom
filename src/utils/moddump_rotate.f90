!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! set solid body rotation for all gas particles about the z-axis given angular frequency
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - omega : *angular frequency of solid body rotation about the z-axis*
!
! :Dependencies: infile_utils, moddump_utils, prompting
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 real :: omega = 7.92e-3   ! angular frequency of solid body rotation about the z-axis

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: vphi,R
 integer                :: i

 if (prompt_for_params) call read_interactive_moddumpfile()

 ! Assume rotation axis is the z-axis
 print*,"Adding solid body rotation with omega = ",omega

 do i = 1,npart
    R = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    vphi = omega*R
    vxyzu(1,i) = -omega*xyzh(2,i)
    vxyzu(2,i) = omega*xyzh(1,i)
    vxyzu(3,i) = 0.
 enddo

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter angular frequency omega for solid body rotation about z',omega)

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
 call read_inopt(omega,'omega',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# solid body rotation parameters'
 call write_inopt(omega,'omega','angular frequency of solid body rotation about the z-axis',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
