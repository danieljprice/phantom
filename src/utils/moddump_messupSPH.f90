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
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - factor  : *gap randomization factor (ioption=2)*
!   - inside  : *randomize inside (T) or outside (F) the Hill sphere (ioption=2)*
!   - ioption : *operation (1=randomize azimuth, 2=randomize gap, 3=delete Hill sphere)*
!
! :Dependencies: mess_up_SPH, moddump_utils, part,
!   prompting, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer :: ioption = 1            ! 1=randomize azimuth, 2=randomize gap, 3=delete Hill sphere
 real    :: factor  = 1.0          ! gap randomization factor (ioption=2)
 logical :: inside  = .false.      ! randomize inside/outside the Hill sphere (ioption=2)

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use mess_up_SPH ! module from MCFOST
 use part,          only:xyzmh_ptmass,nptmass,kill_particle,shuffle_part
 use units,         only:udist
 use moddump_utils, only:get_moddump_options

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, allocatable :: mask(:)
 integer :: i

 print*,'udist=',udist
 ! read parameters from the .moddump file, prompting if it does not exist
 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

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

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 print "(3(/,a))",' 1) randomize azimuth ',&
                  ' 2) randomize gap ',&
                  ' 3) delete Hill sphere'
 call prompt('enter option ',ioption,1,3)
 if (ioption == 2) then
    call prompt('enter factor ',factor)
    call prompt('do you want to randomize inside or outside the Hill sphere? ',inside)
 endif

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(ioption,'ioption','1=randomize azimuth, 2=randomize gap, 3=delete Hill sphere',iunit)
 call write_inopt(factor,'factor','gap randomization factor (ioption=2)',iunit)
 call write_inopt(inside,'inside','randomize inside (T) or outside (F) the Hill sphere (ioption=2)',iunit)

end subroutine write_options_moddump

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(ioption,'ioption',db,min=1,max=3,errcount=nerr)
 call read_inopt(factor,'factor',db,errcount=nerr)
 call read_inopt(inside,'inside',db,errcount=nerr)

end subroutine read_options_moddump

end module moddump

