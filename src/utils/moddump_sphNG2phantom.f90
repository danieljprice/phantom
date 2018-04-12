!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  default moddump routine: does not make any modifications
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, eos, part, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use boundary, only:set_boundary
 use eos,      only:polyk,polyk2
 use units,    only:unit_velocity
 use part,     only:hfact
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

 hfact = 1.2

 return
end subroutine modify_dump

end module moddump

