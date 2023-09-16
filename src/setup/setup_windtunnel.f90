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
!  setup for polytropic gas sphere inside wind tunnel
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
 use inject,      only:init_inject,nstar,Rstar,mach_inf,lattice_type,BHL_handled_layers,&
                       BHL_wind_cylinder_radius,BHL_wind_injection_x,BHL_wind_length,&
                       cs_inf,rho_inf
 use mpidomain,   only:i_belong
 use timestep,    only:dtmax,tmax
 use unifdis,     only:mask_prototype
 use kernel,      only:hfact_default
 use setup_params,only:rhozero,npart_total
 use mpidomain,   only:i_belong
 use table_utils, only:yinterp
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:),vxyzu(:,:),massoftype(:),polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real               :: rhocentre,rmin,tcrush,v_inf,pres_inf,pmass,rho_star,Mstar,densi,presi,ri
 real, allocatable  :: r(:),den(:),pres(:)
 integer            :: ierr,npts,np,i
 logical            :: use_exactN
 character(len=30)  :: lattice
 
 call set_units(mass=solarm,dist=solarr,G=1.)
 !
 !--general parameters
 !
 time  = 0.
 polyk = 1. ! not used but needs to be initialised to non-zero value
 gamma = 5./3.
 ieos  = 2
 gmw   = 0.6
 hfact = hfact_default

 ! Wind parameters (see inject_BHL module)
 mach_inf = 1.55
 rho_inf  = 0.0068
 pres_inf = 5.64e-4
 cs_inf   = sqrt(gamma*pres_inf/rho_inf)
 v_inf    = mach_inf*cs_inf

 ! Star parameters
 Rstar = 0.1
 Mstar = 1.e-3
 nstar = 10000000
 lattice = 'closepacked'
 use_exactN = .true.
 pmass = Mstar / real(nstar)
 massoftype(igas) = pmass

 ! Wind injection settings
 lattice_type = 1
 BHL_handled_layers = 4.
 BHL_wind_cylinder_radius = 10. ! in units of Rstar
 BHL_wind_injection_x = -5.     ! in units of Rstar
 BHL_wind_length = 20.          ! in units of Rstar

 ! Set default tmax and dtmax
 rho_star = Mstar/Rstar**3
 tcrush = 2.*Rstar*sqrt(rho_star/rho_inf)/v_inf
 dtmax = 0.1!1.6*0.05*tcrush
 tmax  = 45.2!1.6*2.5*tcrush
 
 ! Initialise particle injection
 call init_inject(ierr)
 npart = 0
 np = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.

 ! Set polytropic star
 allocate(r(nrhotab),den(nrhotab),pres(nrhotab))
 call rho_polytrope(gamma,polyk,Mstar,r,den,npts,rhocentre,set_polyk=.true.,Rstar=Rstar)
 pres = polyk*den**gamma
 rmin = r(1)
 call set_star_density(lattice,id,master,rmin,Rstar,Mstar,hfact,&
                       npts,den,r,npart,npartoftype,massoftype,xyzh,&
                       use_exactN,np,rhozero,npart_total,i_belong) ! Note: mass_is_set = .true., so np is not used
 ! Set thermal energy
 do i = 1,npart
    ri = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    densi = yinterp(den(1:npts),r(1:npts),ri)
    presi = yinterp(pres(1:npts),r(1:npts),ri)
    vxyzu(4,i) =  presi / ( (gamma-1.) * densi)
 enddo
 nstar = npart
 
 deallocate(r,den,pres)
 
 print*, "udist = ", udist, "; umass = ", umass, "; utime = ", utime

end subroutine setpart
    
end module setup
    