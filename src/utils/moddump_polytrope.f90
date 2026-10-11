!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! perturb a star with a radial velocity perturbation to excite
! the fundamental oscillation mode
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - amp : *amplitude of the radial velocity perturbation*
!
! :Dependencies: moddump_utils, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 real :: amp = 1.e-4   ! amplitude of the radial velocity perturbation

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use moddump_utils, only:get_moddump_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)

 do i=1,npart
    vxyzu(1:3,i) = amp*xyzh(1:3,i)
 enddo

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter amplitude of the radial velocity perturbation',amp)

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

 call read_inopt(amp,'amp',db,errcount=nerr)

end subroutine read_options_moddump

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# fundamental-mode perturbation parameters'
 call write_inopt(amp,'amp','amplitude of the radial velocity perturbation',iunit)

end subroutine write_options_moddump

end module moddump
