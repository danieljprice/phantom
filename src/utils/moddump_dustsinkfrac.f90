!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Changes the dustfrac around sink particles, prevents spurious shadows in RT
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters: None
!
! :Dependencies: dim, part, prompting
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,          only:use_dust
 use part,         only:igas,idust,set_particle_type,ndusttypes,dustfrac
 use prompting,    only:prompt
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,np_gas
 real    :: dust_to_gas,outer_radius,r_g

 if (.not. use_dust) then
    print*,' DOING NOTHING: COMPILE WITH DUST=yes'
    stop
 endif

 dust_to_gas = 0.01

 outer_radius = 10.
 !- grainsize and graindens already set if convert from one fluid to two fluid with growth
 np_gas = npartoftype(igas)

 do i=1,np_gas
    r_g = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
    if (r_g < outer_radius) then
       r_g = r_g/outer_radius
       dustfrac(1:ndusttypes,i) = dust_frac_func(r_g)*dustfrac(1:ndusttypes,i)
    endif
 enddo

 ! massoftype(igas) = massoftype(igas)*(1. + dust_to_gas)
 ! npart = np_gas

end subroutine modify_dump

real function dust_frac_func(x)
 real, intent(in) :: x
 real, parameter :: pi = 4.*atan(1.)

 dust_frac_func = sin(x*pi/2)

end function dust_frac_func

end module moddump
