!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! this module does setup
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: inject, part, physcon, units
!
 use io, only:master

 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for gas sphere inside wind tunnel
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,        only:ihsoft,igas
 use eos,         only:ieos,gmw
 use setstar_utils,only:set_star_density
 use rho_profile, only:rho_polytrope
 use extern_densprofile, only:nrhotab
 use physcon,     only:solarm,solarr
 use units,       only:udist,umass,utime,set_units
 use inject,      only:init_inject,BHL_r_star,BHL_m_star,BHL_mach,BHL_pmass,BHL_closepacked,BHL_handled_layers,&
                       BHL_wind_cylinder_radius,BHL_wind_injection_x,BHL_wind_length,BHL_psep
 use mpidomain,   only:i_belong
 use timestep,    only:dtmax,tmax
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:),vxyzu(:,:),massoftype(:),polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real               :: m,hacc,rhocentre,rmin,tcrush,rho_inf,v_inf,pres_inf,cs_inf,&
                       pmass,element_volume,rho_star
 real, allocatable  :: r(:),den(:)
 integer            :: ierr,npts,nstar
 logical            :: use_exactN
 character(len=30)  :: lattice
 
 call set_units(mass=solarm,dist=solarr,G=1.d0)
 !
 !--general parameters
 !
 time  = 0.
 polyk = 1. ! not used but needs to be initialised to non-zero value
 gamma = 5./3.
 ieos  = 2
 gmw   = 0.6

 ! Wind parameters (see inject_BHL module)
 BHL_mach = 1.31
 rho_inf = 6.8e-5
 pres_inf = 5.9e-6
 cs_inf = sqrt(gamma*pres_inf/rho_inf)
 v_inf = BHL_mach*cs_inf

 ! Star parameters
 BHL_r_star = 0.1
 BHL_m_star = 1.e-3
 nstar = 10000
 pmass = BHL_m_star / real(nstar)
 massoftype(igas) = pmass
 lattice = 'closepacked'
 use_exactN = .true.

 ! Wind injection settings
 BHL_closepacked = 1.           ! do not change, this is hardwired at the moment
 BHL_handled_layers = 4.
 BHL_wind_cylinder_radius = 5.  ! in units of Rstar
 BHL_wind_injection_x = -5.     ! in units of Rstar
 BHL_wind_length = 20.          ! in units of Rstar

 ! Calculate particle separation between layers given rho_inf, depending on lattice type
 element_volume = pmass / rho_inf
 if (BHL_closepacked == 1.) then
    BHL_psep = (sqrt(2.)*element_volume)**(1./3.)
 else
    BHL_psep = element_volume**(1./3.) 
 endif
 BHL_psep = BHL_psep / BHL_r_star  ! need to provide in units of Rstar, separation between layers of wind particle

 ! Set default tmax and dtmax
 rho_star = BHL_m_star/BHL_r_star**3
 tcrush = 2.*BHL_r_star*sqrt(rho_star/rho_inf)/v_inf
 dtmax = 1.6*0.05*tcrush
 tmax  = 1.6*2.5*tcrush
 
 ! Set star
!  allocate(r(nrhotab),den(nrhotab))
!  call rho_polytrope(gamma,polyk,BHL_m_star,r,den,npts,rhocentre,set_polyk=.true.,Rstar=BHL_r_star)
!  rmin = r(1)
!  call set_star_density(lattice,id,master,rmin,BHL_r_star,BHL_m_star,hfact,&
!                        npts,den,r,npart,npartoftype,massoftype,xyzh,use_exactN,np,i_belong)
!  deallocate(r,den)


 call init_inject(ierr)
 npart = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 
 print *, "udist = ", udist, "; umass = ", umass, "; utime = ", utime

end subroutine setpart
    
end module setup
    