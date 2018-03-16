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
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ii
 real    :: phi,A,m,x,y,r
 real    :: rcyl2,rcyl,rsph,Rtorus,v2onr,omegai

 ! Implementing the m=1 density perturbation from Price & Bate 2007

 A = 0.05
 m = 3.0
 Rtorus = 1.0 !Need to check this manually

 ! Adding velocities
 do ii = 1,npart
    rcyl2 = dot_product(xyzh(1:2,ii),xyzh(1:2,ii))
    rcyl  = sqrt(rcyl2)
    rsph  = sqrt(rcyl2 + xyzh(3,ii)*xyzh(3,ii))
    v2onr = 1./(Rtorus)*(-Rtorus*rcyl/rsph**3 + Rtorus**2/(rcyl2*rcyl)) + rcyl/rsph**3

    omegai = sqrt(v2onr/rcyl)
    vxyzu(1,ii) = -omegai*xyzh(2,ii)
    vxyzu(2,ii) = omegai*xyzh(1,ii)
    vxyzu(3,ii) = 0.
 enddo
 print*,'Velocities added.'


 do ii=1,npart
    x=xyzh(1,ii)
    y=xyzh(2,ii)
    r = sqrt(x**2 + y**2)
    phi = atan2(y,x)
    phi = phi - 0.5*A*sin(m*phi)
    xyzh(1,ii) = r*cos(phi)
    xyzh(2,ii) = r*sin(phi)
 enddo

 print*,'Density perturbation added.'

 return
end subroutine modify_dump

end module moddump

