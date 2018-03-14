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
!  rndisc moddump routine: adds a warp in the disc/adds magnetic field
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, physcon, setdisc, setup_params
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use setdisc, only:set_incline_or_warp
 use physcon, only:pi
 use part, only:Bxyz,mhd,rhoh,igas
 use setup_params, only:ihavesetupB
 integer, intent(in)    :: npartoftype(:)
 real,    intent(in)    :: massoftype(:)
 integer, intent(inout) :: npart
 real :: R_warp,H_warp
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: npart_start_count,npart_tot,ii,i
 real    :: beta,rhosum,Bzero,pmassii,phi,incl,posangl
 real    :: rhoc,r2,r2cyl,r,omega,cs,HonR,pressure,psimax
 real    :: vphiold2,vphiold,vadd,vphicorr2

! ihavesetupB=.true.
 HonR = 0.02       ! Must check that this is the same as in setup_rndisc

! Set warp parameters
 R_warp = 2.321!80.
 H_warp = 0.0!20.
 incl = 0.5 ! sine of inclination angle, 0->1
 posangl = 0.

! Similar to that in set_disc
 npart_start_count=1
 npart_tot=npart
!
!---------------------------------------------
! Call setwarp to actually calculate the warp
 call set_incline_or_warp(xyzh,vxyzu,npart_tot,npart_start_count,posangl,incl,&
                          R_warp,H_warp,psimax)
!---------------------------------------------
 do i=npart_start_count,npart_tot
    xyzh(1,i)=xyzh(1,i)
    xyzh(2,i)=xyzh(2,i)
    xyzh(3,i)=xyzh(3,i)
    xyzh(4,i)=xyzh(4,i)

    vxyzu(1,i)=vxyzu(1,i)
    vxyzu(2,i)=vxyzu(2,i)
    vxyzu(3,i)=vxyzu(3,i)
 enddo
 print*,' Disc is now warped '

! Add magnetic field
 if (mhd) then
    beta=10.

! Set up a magnetic field just in Bphi
    do ii = 1,npart
       r2 = xyzh(1,ii)**2 + xyzh(2,ii)**2 + xyzh(3,ii)**2
       r = sqrt(r2)
       phi = atan2(xyzh(2,ii),xyzh(1,ii))
       omega = r**(-1.5)
       cs = HonR*r*omega
       pmassii = massoftype(igas)
       pressure = cs**2*rhoh(xyzh(4,ii),pmassii)
       Bzero = sqrt(2.*pressure/beta)
       Bxyz(1,ii) = -Bzero*sin(phi)
       Bxyz(2,ii) = Bzero*cos(phi)

       ! Calculate correction in v_phi due to B
       vphiold = (-xyzh(2,ii)*vxyzu(1,ii) + xyzh(1,ii)*vxyzu(2,ii))/r
       vphiold2 = vphiold**2
       vphicorr2 = -2.*cs**2
!    if (vphicorr2 > vphi
       vadd = sqrt(vphiold2 + vphicorr2)
       vxyzu(1,ii) = vxyzu(1,ii) + sin(phi)*(vphiold - vadd)
       vxyzu(2,ii) = vxyzu(2,ii) - cos(phi)*(vphiold - vadd)

    enddo

    Bxyz(3,:) = 0.0

    print*,'Magnetic field added.'

! Set up poloidal magnetic field throughout (part of) the disc
!!  Bzero = sqrt(2.*polyk*rhosum/beta)
!!    do ii=1,npart
!!      if (abs(xyzh(3,ii)) < HonR) then  ! to only set the field up in a section of the disc
!!       theta=atan2(xyzh(2,ii),xyzh(1,ii))
!!       Bxyz(1,ii) = 0. !real(Bzero*sin(theta),kind=4)
!!       Bxyz(2,ii) = 0. !real(-Bzero*cos(theta),kind=4)
!!       Bxyz(3,ii) = 0.
!!      endif
!!    enddo
 endif



 return
end subroutine modify_dump

end module moddump

