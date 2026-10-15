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
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - vzero : *amplitude of the Taylor-Green velocity field*
!
! :Dependencies: moddump_utils, physcon, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 real :: vzero = 0.1   ! amplitude of the Taylor-Green velocity field

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use moddump_utils, only:get_moddump_options
 use physcon,       only:pi
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i

 print*,' Adding velocity field for Taylor-Green problem '
 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 print*,' vzero = ',vzero
 do i=1,npart
    vxyzu(1:3,i) = 0.
    !--velocity field for the Taylor-Green problem
    vxyzu(1,i) =  vzero*sin(2.*pi*xyzh(1,i))*cos(2.*pi*xyzh(2,i))
    vxyzu(2,i) = -vzero*cos(2.*pi*xyzh(1,i))*sin(2.*pi*xyzh(2,i))
 enddo

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter amplitude of the Taylor-Green velocity field',vzero)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  read moddump options from the .moddump file
!+
!----------------------------------------------------------------
subroutine read_options_moddump(db,nerr)
 use moddump_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(vzero,'vzero',db,errcount=nerr)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# Taylor-Green velocity field'
 call write_inopt(vzero,'vzero','amplitude of the Taylor-Green velocity field',iunit)

end subroutine write_options_moddump

end module moddump
