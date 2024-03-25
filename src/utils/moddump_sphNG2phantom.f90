!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
! :Runtime parameters: None
!
! :Dependencies: boundary, eos, kernel, part, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use boundary, only:set_boundary
 use eos,      only:polyk,polyk2
 use units,    only:unit_velocity
 use part,     only:hfact
 use kernel,   only:hfact_default
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: xmini,xmaxi,ymini,ymaxi,zmini,zmaxi
 real :: cs_sphere, cs_medium, density_contrast

 print*,' *** modifying sphNG dump file ***'
 print*,' the sound speed in code units is ',sqrt(polyk)
 print*,' the sound speed in cm/s       is ',sqrt(polyk)*unit_velocity

 print*,' enter new boundary settings (xmin,xmax,ymin,ymax,zmin,zmax)'
 read*,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi
 call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)

 print*,' enter density contrast between sphere and medium (needed to set polyk2)'
 read*,density_contrast

 cs_sphere = sqrt(polyk)
 cs_medium   = cs_sphere*sqrt(density_contrast)
 polyk2 = cs_medium**2
 print*,' the sound speed in the medium is ',sqrt(polyk2)
 print*,'                          in cm/s ',sqrt(polyk2)*unit_velocity

 hfact = hfact_default

end subroutine modify_dump

end module moddump
