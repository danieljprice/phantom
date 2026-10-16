!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Set rigid rotation about the mass-weighted centre, retaining bulk velocity.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - omega : *angular frequency of solid body rotation about the z-axis*
!
! :Dependencies: centreofmass, damping, io, moddump_utils, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 real :: omega = 7.92e-3   ! angular frequency of solid body rotation about the z-axis

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use moddump_utils, only:get_moddump_options
 use centreofmass, only:get_centreofmass
 use damping,     only:rigid_rotation_velocity
 use part,        only:nptmass,igas,isdead_or_accreted
 use io,          only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: xcom(3),vcom(3)
 integer                :: i

 call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)
 if (npart <= 0 .or. nptmass > 0 .or. npartoftype(igas) /= npart) &
    call fatal('moddump_rotate','rigid spin requires an isolated gas star without sinks')
 call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu)

 ! Assume rotation axis is the z-axis
 print*,"Adding solid body rotation with omega = ",omega

 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    vxyzu(1:3,i) = vcom + rigid_rotation_velocity(xyzh(1:3,i)-xcom,[0.,0.,omega])
 enddo

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactive prompting for moddump parameters
!+
!----------------------------------------------------------------
subroutine moddump_interactive()
 use prompting, only:prompt

 call prompt('Enter angular frequency omega for solid body rotation about z',omega)

end subroutine moddump_interactive

!----------------------------------------------------------------
!+
!  write moddump options to the .moddump file
!+
!----------------------------------------------------------------
subroutine write_options_moddump(iunit)
 use moddump_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# solid body rotation parameters'
 call write_inopt(omega,'omega','angular frequency of solid body rotation about the z-axis',iunit)

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

 call read_inopt(omega,'omega',db,errcount=nerr)

end subroutine read_options_moddump

end module moddump
