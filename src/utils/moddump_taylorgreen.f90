!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Add Taylor-Green velocity field
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - vzero : *amplitude of the Taylor-Green velocity field*
!
! :Dependencies: infile_utils, moddump_utils, physcon, prompting
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 real :: vzero = 0.1   ! amplitude of the Taylor-Green velocity field

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use physcon,       only:pi

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i

 print*,' Adding velocity field for Taylor-Green problem '
 if (prompt_for_params) call read_interactive_moddumpfile()

 print*,' vzero = ',vzero
 do i=1,npart
    vxyzu(1:3,i) = 0.
    !--velocity field for the Taylor-Green problem
    vxyzu(1,i) =  vzero*sin(2.*pi*xyzh(1,i))*cos(2.*pi*xyzh(2,i))
    vxyzu(2,i) = -vzero*cos(2.*pi*xyzh(1,i))*sin(2.*pi*xyzh(2,i))
 enddo

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter amplitude of the Taylor-Green velocity field',vzero)

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
 call read_inopt(vzero,'vzero',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# Taylor-Green velocity field'
 call write_inopt(vzero,'vzero','amplitude of the Taylor-Green velocity field',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
