!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    T_wind              -- wind temperature (K)
!    primary_racc    -- accretion radius of the central star (au)
!    central_star_mass   -- mass of the central star (Msun)
!    secondary_star_mass -- mass of the secondary star (Msun)
!    secondary_star_r    -- radius of the secondary star (au)
!    eccentricity        -- eccentricity of the binary system
!    isecondary_star     -- set to 1 for a binary system
!    initial_wind_gamma  -- polytropic index
!    mass_of_particles   -- mass resolution (Msun)
!    semi_major_axis     -- semi-major axis of the binary system (au)
!    wind_gamma          -- polytropic index
!
!  DEPENDENCIES: eos, infile_utils, inject, io, part, physcon, prompting,
!    readwrite_infile, setbinary, units
!+
!--------------------------------------------------------------------------
module setup

 implicit none
 public :: setpart

 private
 real, public :: wind_gamma = 5./3.
#ifdef ISOTHERMAL
 real, public :: T_wind = 50000.
#else
 real, public :: T_wind = 3000.
#endif
 integer, public :: icompanion_star = 0
 real, public :: spacial_resolution
 real :: primary_Teff = 3000.
 real :: primary_Reff
 real :: primary_lum
 real :: primary_mass
 real :: primary_racc
 real :: secondary_Teff = 0.
 real :: secondary_Reff
 real :: secondary_lum
 real :: secondary_mass
 real :: secondary_racc
 real, public :: semi_major_axis
 real, public :: eccentricity
 real :: default_particle_mass = 1.e-11
 real :: spacial_resolution_au = 0.01
 real :: primary_lum_lsun = 5315.
 real :: primary_mass_msun = 1.2
 real :: primary_Reff_au = 0.
 real :: primary_racc_au = 1.
 real :: secondary_lum_lsun = 0.
 real :: secondary_mass_msun=0.6
 real :: secondary_Reff_au = 0.
 real :: secondary_racc_au
 real :: temp_exponent = 0.

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only: xyzmh_ptmass, vxyz_ptmass, nptmass, igas, iTeff, iLum, iReff
 use physcon,   only: au, solarm, mass_proton_cgs, kboltz
 use units,     only: umass, set_units,unit_velocity
 use inject,    only: init_inject
 use setbinary, only: set_binary
 use io,        only: master
 use eos,       only: gmw,ieos,isink,qfacdisc
 use spherical, only:set_sphere
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=*),  intent(in)    :: fileprefix
 character(len=len(fileprefix)+6) :: filename
 integer :: ierr
 logical :: iexist

 call set_units(dist=au,mass=solarm,G=1.)
!
!--general parameters
!
 time = 0.
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call setup_interactive()
       call write_setupfile(filename)
    endif
 endif

!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.

 if (icompanion_star > 0) then
    call set_binary(primary_mass, &
                    secondary_mass/primary_mass, &
                    semi_major_axis, &
                    eccentricity, &
                    primary_racc, &
                    secondary_racc, &
                    xyzmh_ptmass, vxyz_ptmass, nptmass)
    xyzmh_ptmass(iTeff,1) = primary_Teff
    xyzmh_ptmass(iReff,1) = primary_Reff
    xyzmh_ptmass(iLum,1)  = primary_lum
    xyzmh_ptmass(iTeff,2) = secondary_Teff
    xyzmh_ptmass(iReff,2) = secondary_Reff
    xyzmh_ptmass(iLum,2)  = secondary_lum
 else
    nptmass = 1
    xyzmh_ptmass(4,1) = primary_mass
    xyzmh_ptmass(5,1) = primary_racc
    xyzmh_ptmass(iTeff,1) = primary_Teff
    xyzmh_ptmass(iReff,1) = primary_Reff
    xyzmh_ptmass(iLum,1)  = primary_lum
 endif
 !
 ! for binary wind simulations this mass is IRRELEVANT
 ! since it will be over-written on the first call to init_inject
 !
 massoftype(igas) = default_particle_mass * (solarm / umass)

#ifdef ISOTHERMAL
 gamma = 1.
#else
 gamma = wind_gamma
#endif
 if (gamma < 1.0001) then
    if (temp_exponent > 0.) then
       ieos = 6
       qfacdisc = 0.5*temp_exponent
       isink = 1
       polyk = kboltz*primary_Teff/(mass_proton_cgs * gmw * unit_velocity**2)
    else
       ieos = 1
       polyk = kboltz*T_wind/(mass_proton_cgs * gmw * unit_velocity**2)
    endif
 endif

!
! add low density background medium
!
! sphere of 20^3 particles, density will be determined when mass is set
! call set_sphere('cubic',id,master,primary_racc,20.*primary_racc,primary_racc,hfact,npart,xyzh,rhofunc=rhor)
! npartoftype(igas) = npart

! contains
!  real function rhor(r)
!    real, intent(in) :: r
!
!    rhor = 1./r**2
!
!  end function rhor

end subroutine setpart

!----------------------------------------------------------------
!+
!  determine which problem to set up interactively
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 use physcon,   only:au,solarm
 use units,     only:umass,udist
 integer :: iproblem,iwind

 iproblem = 1
 iwind = 2
 call prompt('Type of wind: 1=isoT, 2=adia, 3=T(r)',iwind,1,3)
 if (iwind == 1 .or. iwind == 3) wind_gamma = 1.

 call prompt('Add binary?',icompanion_star,0,1)
 print "(a)",' 2: Mass: 1.2 Msun, accretion radius: 0.2568 au',&
      ' 1: Mass: 1.0 Msun, accretion radius: 1.2568 au', &
      ' 0: custom'
 call prompt('select mass and radius of primary',iproblem,0,2)
 select case(iproblem)
 case(2)
    primary_mass = 1.2 * (solarm / umass)
    primary_racc = 0.2568 * (au / udist)
 case(1)
    primary_mass = 1. * (solarm / umass)
    primary_racc = 1.2568 * (au / udist)
 case default
    primary_mass = 1.
    primary_racc = 1.
    call prompt('enter primary mass',primary_mass,0.,100.)
    call prompt('enter accretion radius in au ',primary_racc,0.)
    primary_mass = primary_mass * (solarm / umass)
    primary_racc = primary_racc * (au / udist)
 end select
 primary_racc_au = primary_racc*udist/au

end subroutine setup_interactive

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use physcon,      only:au,steboltz,solarl,solarm,pi
 use units,        only:utime,unit_energ
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for wind setup routine'
 if (primary_Teff == 0. .and. primary_lum_lsun > 0. .and. primary_Reff_au > 0.) &
      primary_Teff = (primary_lum_lsun*solarl/(4.*pi*steboltz*primary_Reff_au*au))**0.25
 if (primary_Reff_au == 0. .and. primary_lum_lsun > 0. .and. primary_Teff > 0.) &
      primary_Reff_au = sqrt(primary_lum_lsun*solarl/(4.*pi*steboltz*primary_Teff**4))/au
 primary_lum = primary_lum_lsun * (solarl * utime / unit_energ)
 call write_inopt(primary_mass_msun,'primary_mass','primary star mass (Msun)',iunit)
 call write_inopt(primary_racc_au,'primary_racc','primary star accretion radius (au)',iunit)
 call write_inopt(primary_lum_lsun,'primary_lum','primary star luminosity (Lsun)',iunit)
 call write_inopt(primary_Teff,'primary_Teff','primary star effective temperature (K)',iunit)
 call write_inopt(primary_Reff_au,'primary_Reff','primary star effective radius (au)',iunit)
 call write_inopt(icompanion_star,'icompanion_star','set to 1 for a binary system',iunit)
 if (icompanion_star > 0) then
    if (secondary_Teff == 0. .and. secondary_lum_lsun > 0. .and. secondary_Reff_au > 0.) &
         secondary_Teff = (secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Reff_au*au))**0.25
    if (secondary_Reff_au == 0. .and. secondary_lum_lsun > 0. .and. secondary_Teff > 0.) &
         secondary_Reff_au = sqrt(secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Teff**4))/au
    secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
    call write_inopt(secondary_mass_msun,'secondary_mass','secondary star mass (Msun)',iunit)
    call write_inopt(secondary_racc_au,'secondary_racc','secondary star accretion radius (au)',iunit)
    call write_inopt(secondary_lum_lsun,'secondary_lum','secondary star luminosity (Lsun)',iunit)
    call write_inopt(secondary_Teff,'secondary_Teff','secondary star effective temperature)',iunit)
    call write_inopt(secondary_Reff_au,'secondary_Reff','secondary star effective radius (au)',iunit)
    call write_inopt(semi_major_axis,'semi_major_axis','semi-major axis of the binary system (au)',iunit)
    call write_inopt(eccentricity,'eccentricity','eccentricity of the binary system',iunit)
 endif
 call write_inopt(default_particle_mass,'mass_of_particles','mass resolution (Msun)',iunit)

#ifdef ISOTHERMAL
 wind_gamma = 1.
#else
 call write_inopt(wind_gamma,'wind_gamma','adiabatic index (initial if Krome chemistry used)',iunit)
#endif
 if ( wind_gamma < 1.0001) then
    call write_inopt(temp_exponent,'temp_exponent','temperature profile T = R^-p (0 = isothermal)',iunit)
    if (abs(temp_exponent) < tiny(0.)) call write_inopt(T_wind,'T_wind','temperature (K)',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use physcon,      only:au,steboltz,solarl,solarm,pi
 use units,        only:udist,umass,utime,unit_energ
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)
 integer :: nerr,ichange

 nerr = 0
 ichange = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(primary_mass_msun,'primary_mass',db,min=0.,max=1000.,errcount=nerr)
 primary_mass = primary_mass_msun * (solarm / umass)
 call read_inopt(primary_lum_lsun,'primary_lum',db,min=0.,max=1e7,errcount=nerr)
 primary_lum = primary_lum_lsun * (solarl * utime / unit_energ)
 call read_inopt(primary_Teff,'primary_Teff',db,min=0.,errcount=nerr)
 call read_inopt(primary_Reff_au,'primary_Reff',db,min=0.,errcount=nerr)
 primary_Reff = primary_Reff_au * au / udist
 call read_inopt(primary_racc_au,'primary_racc',db,min=0.,errcount=nerr)
 primary_racc = primary_racc_au * au / udist
 call read_inopt(icompanion_star,'icompanion_star',db,min=0,errcount=nerr)
 if (icompanion_star > 0) then
    call read_inopt(secondary_mass_msun,'secondary_mass',db,min=0.,max=1000.,errcount=nerr)
    secondary_mass = secondary_mass_msun * (solarm / umass)
    call read_inopt(secondary_lum_lsun,'secondary_lum',db,min=0.,max=1000.,errcount=nerr)
    secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
    call read_inopt(secondary_Teff,'secondary_Teff',db,min=0.,max=1000.,errcount=nerr)
    call read_inopt(secondary_Reff_au,'secondary_Reff',db,min=0.,max=1000.,errcount=nerr)
    secondary_Reff = secondary_Reff_au * au / udist
    call read_inopt(secondary_racc_au,'secondary_racc',db,min=0.,max=1000.,errcount=nerr)
    secondary_racc = secondary_racc_au * au / udist
    call read_inopt(semi_major_axis,'semi_major_axis',db,min=0.,errcount=nerr)
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
 endif
 call read_inopt(default_particle_mass,'mass_of_particles',db,min=0.,errcount=nerr)
#ifdef ISOTHERMAL
 wind_gamma = 1.
#else
 call read_inopt(wind_gamma,'wind_gamma',db,min=1.,max=4.,errcount=nerr)
#endif
 if ( wind_gamma < 1.0001) then
    call read_inopt(temp_exponent,'temp_exponent',db,min=0.,max=5.,errcount=nerr)
    if (abs(temp_exponent) < tiny(0.)) call read_inopt(T_wind,'T_wind',db,min=0.,errcount=nerr)
 endif
 call close_db(db)
 if (primary_Teff == 0. .and. primary_lum_lsun > 0. .and. primary_Reff > 0.) then
    ichange = ichange+1
    primary_Teff = (primary_lum_lsun*solarl/(4.*pi*steboltz*primary_Reff*udist))**0.25
 endif
 if (primary_Reff == 0. .and. primary_lum_lsun > 0. .and. primary_Teff > 0.) then
    ichange = ichange+1
    primary_Reff = sqrt(primary_lum_lsun*solarl/(4.*pi*steboltz*primary_Teff**4))/udist
    primary_Reff_au = primary_Reff * udist / au
 endif
 if (icompanion_star > 0) then
    if (secondary_Teff == 0. .and. secondary_lum_lsun > 0. .and. secondary_Reff > 0.) then
       ichange = ichange+1
       secondary_Teff = (secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Reff*udist))**0.25
    endif
    if (secondary_Reff == 0. .and. secondary_lum_lsun > 0. .and. secondary_Teff > 0.) then
       ichange = ichange+1
       secondary_Reff = sqrt(secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Teff**4))/udist
       secondary_Reff_au = secondary_Reff * udist / au
    endif
 endif
 ierr = nerr
 if (ichange > 0 .and. ierr == 0) call write_setupfile(filename)

end subroutine read_setupfile
end module setup
