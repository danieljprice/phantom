!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
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
 use setdisc, only:set_warp
 use physcon, only:pi
 use part, only:Bevol,mhd,maxBevol,rhoh,igas
 use setup_params, only:ihavesetupB
 integer, intent(in)    :: npartoftype(:)
 real,    intent(in)    :: massoftype(:)
 integer, intent(inout) :: npart
 real :: sininclination,rwarp,warp_smoothl,rsi,rso
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: npart_start_count,npart_tot,ii,i
 logical :: do_twist
 real    :: beta,rhosum,Bzero,pmassii,phi,inclination
 real    :: rhoc,r2,r2cyl,r,omega,cs,HonR,pressure,psimax
 real    :: vphiold2,vphiold,vadd,vphicorr2

! ihavesetupB=.true.
 HonR = 0.02       ! Must check that this is the same as in setup_rndisc

! Set warp parameters
 rwarp = 2.321!80.
 warp_smoothl = 0.0!20.
 inclination = 0.5 ! sine of inclination angle, 0->1

 rsi = rwarp - warp_smoothl
 rso = rwarp + warp_smoothl

! Similar to that in set_disc
 do_twist=.true.
 npart_start_count=1
 npart_tot=npart
!
!---------------------------------------------
! Call setwarp to actually calculate the warp
 call set_warp(npart_tot,npart_start_count,&
               xyzh,vxyzu,inclination,sininclination,&
               rwarp,psimax,rsi,rso,do_twist)
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
       Bevol(1,ii) = -Bzero*sin(phi)
       Bevol(2,ii) = Bzero*cos(phi)

       ! Calculate correction in v_phi due to B
       vphiold = (-xyzh(2,ii)*vxyzu(1,ii) + xyzh(1,ii)*vxyzu(2,ii))/r
       vphiold2 = vphiold**2
       vphicorr2 = -2.*cs**2
!    if (vphicorr2 > vphi
       vadd = sqrt(vphiold2 + vphicorr2)
       vxyzu(1,ii) = vxyzu(1,ii) + sin(phi)*(vphiold - vadd)
       vxyzu(2,ii) = vxyzu(2,ii) - cos(phi)*(vphiold - vadd)

    enddo

    Bevol(3,:) = 0.0

    print*,'Magnetic field added.'

! Set up poloidal magnetic field throughout (part of) the disc
!!  Bzero = sqrt(2.*polyk*rhosum/beta)
!!    do ii=1,npart
!!      if (abs(xyzh(3,ii)) < HonR) then  ! to only set the field up in a section of the disc
!!       theta=atan2(xyzh(2,ii),xyzh(1,ii))
!!       Bevol(1,ii) = 0. !real(Bzero*sin(theta),kind=4)
!!       Bevol(2,ii) = 0. !real(-Bzero*cos(theta),kind=4)
!!       Bevol(3,ii) = 0.
!!      endif
!!    enddo
 endif



 return
end subroutine modify_dump

end module moddump

