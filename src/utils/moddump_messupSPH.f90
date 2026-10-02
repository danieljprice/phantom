!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! default moddump routine: does not make any modifications
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - factor  : *gap randomization factor (ioption=2)*
!   - inside  : *randomize inside (T) or outside (F) the Hill sphere (ioption=2)*
!   - ioption : *operation (1=randomize azimuth, 2=randomize gap, 3=delete Hill sphere)*
!
! :Dependencies: infile_utils, mess_up_SPH, moddump_utils, part, prompting,
!   units
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer :: ioption = 1            ! 1=randomize azimuth, 2=randomize gap, 3=delete Hill sphere
 real    :: factor  = 1.0          ! gap randomization factor (ioption=2)
 logical :: inside  = .false.      ! randomize inside/outside the Hill sphere (ioption=2)

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use mess_up_SPH ! module from MCFOST
 use part,      only:xyzmh_ptmass,nptmass,kill_particle,shuffle_part
 use units,     only:udist

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, allocatable :: mask(:)
 integer :: i

 print*,'udist=',udist
 !
 ! prompt only when the driver did not find a parameter file
 !
 if (prompt_for_params) call read_interactive_moddumpfile()

 allocate(mask(npart))
 mask = 0
 select case(ioption)
 case(3)
    call mask_Hill_sphere(npart,nptmass,xyzh,xyzmh_ptmass,udist,mask)
 case(2)
    print*,' randomizing gap...'
    call randomize_gap(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,udist,real(factor,kind=8),inside)
 case(1)
    print*,' randomizing azimuths...'
    call randomize_azimuth(npart,xyzh,vxyzu)
 case default
    print*,' no modifications performed '
 end select
 print*,' done ',count(mask==1),' particles killed'
 if (any(mask==1)) then
    do i=1,npart
       if (mask(i) == 1) call kill_particle(i,npartoftype)
    enddo
    call shuffle_part(npart)
    print*,' NEW NUMBER OF PARTICLES = ',npart
 endif
 deallocate(mask)

end subroutine modify_dump

!
!---Interactively set the moddump parameters--------------------------------
!
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 print "(3(/,a))",' 1) randomize azimuth ',&
                  ' 2) randomize gap ',&
                  ' 3) delete Hill sphere'
 call prompt('enter option ',ioption,1,3)
 if (ioption == 2) then
    call prompt('enter factor ',factor)
    call prompt('do you want to randomize inside or outside the Hill sphere? ',inside)
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
 call write_inopt(ioption,'ioption','operation (1=randomize azimuth, 2=randomize gap, 3=delete Hill sphere)',iunit)
 call write_inopt(factor,'factor','gap randomization factor (ioption=2)',iunit)
 call write_inopt(inside,'inside','randomize inside (T) or outside (F) the Hill sphere (ioption=2)',iunit)
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
 call read_inopt(ioption,'ioption',db,min=1,max=3,errcount=nerr)
 call read_inopt(factor,'factor',db,errcount=nerr)
 call read_inopt(inside,'inside',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

end module moddump

