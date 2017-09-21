!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  this module does general accretion disc setups
!  Modified from an original routine by Giuseppe Lodato
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, extern_lensethirring, externalforces, io, options,
!    part, physcon, setdisc, setup_params, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,   only:set_disc
! use setbinary, only:set_binary
 use dim,          only:maxp,maxvxyzu
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihsoft,ihacc
 use part,         only:Bevol,mhd,maxBevol,rhoh
 use part,         only:iphase,iamtype,igas,maxphase
 use setup_params, only:ihavesetupB,rhozero
 use io,             only:master
 use externalforces, only:accradius1,iext_star
 use options,        only:iexternalforce,alpha
 use extern_lensethirring, only:blackhole_spin
 use units,          only:set_units,udist
 use physcon,        only:solarm,km
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: R_in,R_out,R_warp,ampl,Rs,theta!,xinc
 integer :: ii
 real :: Bzero,beta,HonR,rhosum,pmassii,r2,cs,r,r2cyl,omega,rhoc,pressure,phi
 real :: vzero

 call set_units(mass=1.*solarm,c=1.)

 !
 !  Set problem parameters
 !
 !--disc inner and outer radius
 Rs = 2.0
 R_in    = 4.233
 R_out   = 100.
 R_warp  = R_in
! npart   = size(xyzh(1,:))
 npart = 100000
 npartoftype(:) = 0
 npartoftype(1) = npart
 gamma   = 5./3. !Unless otherwise specified, 1.0
 time    = 0.
 hfact   = 1.0
 !xinc    = 45.*(3.141592654/180.0) ! Must be in radians
 ampl    = 0.00 !259 ! sine of inclination angle, 0->1
 alpha   = 0.1 !0.0096 !0443 !0206
 blackhole_spin = 0.5
 HonR = 0.1

 iexternalforce = 11
 accradius1 = R_in

 call set_disc(id,master=master,&
                nparttot  = npart,  &
                npart     = npart,  &
                rmin      = R_in,   &
                rmax      = R_out,  &
                rwarp     = R_warp, &
                p_index   = -1.0,   &
                q_index   = 0.5,    &
                HoverR    = HonR,   &
                disc_Q    = 168.,   &
                star_mass = 1.0,    &
                gamma     = gamma,  &
                particle_mass = massoftype(igas), &
                hfact     = hfact,  &
                xyzh      = xyzh,   &
                vxyzu     = vxyzu,  &
                polyk     = polyk,  &
                twist     = .false., &
                alpha     = alpha,  &
                ismooth   = .true., &
                sininclination = ampl, &
                warp_smoothl = 0.,  &
                bh_spin = blackhole_spin, &
                prefix = fileprefix)


! Use a point mass instead
! nptmass = 1
! xyzmh_ptmass(:,:) = 0.0
! xyzmh_ptmass(4,1) = 1.0
! xyzmh_ptmass(5,1) = 0.25*R_in
! vxyz_ptmass(:,:) = 0.0
! xyzmh_ptmass(ihsoft,1) = 0.0
! xyzmh_ptmass(ihacc,1) = 0.25*R_in

! polyk = 0.1717**2

! Add magnetic field
 if (mhd) then
    beta=1348.
    ! ! Calculate Bzero
    pmassii = massoftype(igas)
    rhosum = 0.0
    do ii=1,npart
       if (maxphase==maxp) pmassii = massoftype(iamtype(iphase(ii)))
       rhosum = rhosum + rhoh(xyzh(4,ii),pmassii)
    enddo
    rhosum = rhosum/npart
    rhoc = 1.0 !central density

    Bevol(:,:) = 0.0

    ! Set up magnetic field as in Gaberov and Vanaverbeke et al 2014
    do ii=1,npart
       r2 = xyzh(1,ii)**2 + xyzh(2,ii)**2 + xyzh(3,ii)**2
       phi = atan2(xyzh(2,ii),xyzh(1,ii))
       if ((r2 > 4.) .and. (r2 < 9)) then
          r2cyl = xyzh(1,ii)**2 + xyzh(2,ii)**2
          r = sqrt(r2cyl)
          omega = r**(-1.5)
          cs = HonR*r*omega
          vzero = 0.1*cs
          pressure = cs**2*rhoc
          Bzero = sqrt(2.*pressure/beta)
          Bevol(3,ii) = Bzero*sin(6.2832*(sqrt(r2)-2.))
          !Bevol(1,ii) = Bzero*(-sin(phi)) !sin(6.2832*(sqrt(r2)-2.))
          !Bevol(2,ii) = Bzero*cos(phi) !sin(6.2832*(sqrt(r2)-2.))
       else
          Bevol(3,ii) = 0.0
       endif
       Bevol(1,ii) = 0.0
       Bevol(2,ii) = 0.0

       ! Add a velocity perturbation as in PB13a
       !vxyzu(3,ii) = vxyzu(3,ii) + vzero*cos(2.5*r + 10.*phi + 5.*xyzh(3,ii))
    enddo

    print*,'Using this density to set B value: ',rhosum
    print*,'Original field in z direction has strength: ',Bzero
    print*,'And magnetic field is set up in disc.'

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

    ihavesetupB=.true.
 endif

 return
end subroutine setpart

end module setup
