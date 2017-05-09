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
!  This module sets up a single, warped accretion disc
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: externalforces, io, mpiutils, options, physcon, setdisc,
!    timestep, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up warped discs
! As in Lodato & Price (2010)
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,        only:set_disc
 use io,             only:master,fatal
 use externalforces, only:accradius1
 use options,        only:iexternalforce
 use timestep,       only:dtmax,tmax
 use physcon,        only:pi,au,solarm
 use mpiutils,       only:bcast_mpi
 use units,          only:set_units
 use prompting,      only:prompt
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: R_in,R_out,R_warp,ampl,dr_warp,psimax,term
 real    :: HonR,q_index,csfac,omegfac,discmass
 !
 !  Set problem parameters
 !
 call set_units(dist=au,mass=solarm,G=1.)

 !--disc inner and outer radius
 R_in    = 0.5
 R_out   = 10.
 R_warp  = 5.
 npart   = 100000. !size(xyzh(1,:))
 npartoftype(:) = 0
 npartoftype(1) = npart
 gamma   = 1.0
 time    = 0.
 dr_warp = 1.5
 term    = pi*R_warp/(4.*dr_warp)
 if (id==master) then
    psimax = 0.05
    call prompt(' Enter warp amplitude (maximum psi)',psimax,0.,1.)
 endif
 call bcast_mpi(psimax)
 if (psimax < 0. .or. psimax > 10.) call fatal('setup','error in psi-max')
 !psimax  = 0.2
 ampl    = psimax/(term*sqrt(1. + (psimax/(2.*term))**2))
 if (id==master) print*,' using warp amplitude ',psimax,' gives A = ',ampl
! ampl    = 0.5  ! sine of inclination angle, 0->1
 q_index = 0.75
!
! set_disc wants H/R at R_in, but we want to specify at R=1
!
 csfac   = R_in**(-q_index)
 omegfac = 1./R_in**(1.5)
 HonR    = 0.04*(csfac/omegfac)/R_in
 discmass = 1.525020711459903E-004  ! based on Q = 168 for H/R=0.02

 call set_disc(id,master=master,&
                npart   = npartoftype(1),&
                rmin    = R_in, &
                rmax    = R_out,&
                rwarp   = R_warp, &
                p_index = 1.5,    &
                q_index = q_index,   &
                HoverR  = HonR, &
                disc_mass = discmass,   &
                star_mass = 1.0,    &
                gamma   = gamma,  &
                particle_mass = massoftype(1), &
                hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,polyk=polyk,&
                ismooth=.true., &
                sininclination=ampl, &
                warp_smoothl=dr_warp, &
                prefix = fileprefix)
 !
 !--set default options for the input file
 !
 iexternalforce = 1
 accradius1 = R_in
 dtmax = 5.
 tmax  = 1000.

 return
end subroutine setpart

end module setup
