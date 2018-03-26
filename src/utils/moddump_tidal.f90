!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:                                                                   !
!  Puts star into a parabolic orbit                                            !
!                                                                           !
!  REFERENCES: None                                                           !
!                                                                           !
!  OWNER: Daniel Price                                                           !
!                                                                           !
!  $Id$                           !
!                                                                           !
!  RUNTIME PARAMETERS: None                                                   !
!                                                                           !
!  DEPENDENCIES: None                                                           !
!
!  REFERENCES: None                                                           !
!                                                                           !
!  OWNER: Daniel Price                                                           !
!                                                                           !
!  $Id$                           !
!                                                                           !
!  RUNTIME PARAMETERS: None                                                   !
!                                                                           !
!  DEPENDENCIES: None                                                           !
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, externalforces, options, prompting
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use centreofmass
 use externalforces, only:mass1
 use externalforces, only:accradius1
 use options,        only:iexternalforce
 use prompting,      only:prompt
 integer,  intent(inout) :: npart
 integer,  intent(inout) :: npartoftype(:)
 real,     intent(inout) :: massoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                 :: i
 real                    :: beta, b, rt, rp, r, rs, Ms, Mh

 !--Reset center of mass
 call reset_centreofmass(npart,xyzh,vxyzu)

 !--Defaults
 beta = 1.0                   ! penetration factor
 Mh = 1.e6                    ! BH mass
 Ms = 1.0                     ! stellar mass
 rs = 1.0                     ! stellar radius

 !--User enter values
 call prompt(' Enter a value for the penetration factor (beta): ',beta,0.)
 call prompt(' Enter a value for blackhole mass (in code units): ',Mh,0.)
 call prompt(' Enter a value for the stellar mass (in code units): ',Ms,0.)
 call prompt(' Enter a value for the stellar radius (in code units): ',rs,0.)

 rt = (Mh/Ms)**(1./3.) * rs   ! tidal radius
 rp = rt/beta                 ! pericenter distance
 b = sqrt(2.)*rp              ! impact parameter (when b=x)

 !--Set input file parameters
 mass1 = Mh
 iexternalforce = 1
 accradius1 = (2*Mh*rs)/((6.8565e2)**2) ! R_sch = 2*G*Mh*rs/c**2

 !--Putting star into orbit
 do i = 1, npart
    !--translate star by x=-b, y=b
    !
    xyzh(1,i) = xyzh(1,i) - b
    xyzh(2,i) = xyzh(2,i) + b
    xyzh(3,i) = xyzh(3,i)
    xyzh(4,i) = xyzh(4,i)
    !
    r = sqrt(2.)*b ! when b=x
    !--giving star a velocity
    !--velocity of star in parabolic orbit is sqrt(2GM/r)
    !
    vxyzu(1,i) = sqrt(2.*Mh/r)
    vxyzu(2,i) = 0.
    vxyzu(3,i) = 0.
    !
 enddo

 write(*,'(a)') "======================================================================"
 write(*,'(a,Es12.5,a)') ' Pericenter distance = ',rp,' R_sun'
 write(*,'(a,Es12.5,a)') ' Tidal radius        = ',rt,' R_sun'
 write(*,'(a,Es12.5,a)') ' Impact parameter    = ',b,' R_sun'
 write(*,'(a,Es12.5,a)') ' Radius of star      = ',rs,' R_sun'
 write(*,'(a,Es12.5,a)') ' Stellar mass        = ',Ms,' M_sun'
 write(*,'(a)') "======================================================================"

 return
end subroutine modify_dump

end module moddump

