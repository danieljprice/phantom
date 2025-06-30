!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up a sphere-in-a-box: a cold, dense sphere placed in
!   a warm medium; the two media are in pressure-equilibrium.
!   This currently works for gas-only and one-fluid dust.
!   Set density_contrast=1 to simulate decaying turbulence in a box.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - BEfac            : *over-density factor of the BE sphere [code units]*
!   - BEmass           : *mass radius of the BE sphere [code units]*
!   - BErad_norm       : *normalised radius of the BE sphere*
!   - BErad_phys       : *physical radius of the BE sphere [code units]*
!   - BErho_cen        : *central density of the BE sphere [code units]*
!   - Bzero            : *Magnetic field strength in Gauss*
!   - ang_Bomega       : *Angle (degrees) between B and rotation axis*
!   - angvel           : *angular velocity in rad/s*
!   - beta_r           : *rotational-to-gravitational energy ratio*
!   - cs_sphere        : *sound speed in sphere (with units e.g. 0.2 km/s)*
!   - density_contrast : *density contrast in code units*
!   - form_binary      : *the intent is to form a central binary*
!   - h_acc            : *accretion radius (with units; e.g. au,pc,kpc,0.1pc)*
!   - iBE_options      : *The set of parameters to define the BE sphere*
!   - icreate_sinks    : *1: create sinks.  0: do not create sinks*
!   - lattice          : *particle lattice (random,cubic,closepacked,hcp,hexagonal)*
!   - lbox             : *length of a box side in terms of spherical radii*
!   - masstoflux       : *mass-to-magnetic flux ratio in units of critical value*
!   - np               : *requested number of particles in sphere*
!   - r_sphere         : *radius of sphere in code units*
!   - rho_final        : *final maximum density (<=0 to ignore) (cgs units)*
!   - rho_pert_amp     : *amplitude of density perturbation*
!   - rms_mach         : *turbulent rms mach number*
!   - shuffle_parts    : *relax particles by shuffling*
!   - totmass_sphere   : *mass of sphere in code units*
!   - use_BE_sphere    : *centrally condense as a BE sphere*
!
! :Dependencies: boundary, centreofmass, datafiles, dim, dust, eos,
!   eos_barotropic, infile_utils, io, kernel, mpidomain, options, part,
!   physcon, prompting, ptmass, rho_profile, set_dust_options, setunits,
!   setup_params, spherical, timestep, unifdis, units,
!   utils_shuffleparticles, velfield
!
 use part,     only:mhd,graindens,grainsize,ndusttypes,ndustsmall,ndustlarge
 use dim,      only:use_dust,maxvxyzu,periodic,maxdustsmall,gr,isothermal
 use options,  only:calc_erot,use_dustfrac
 use setunits, only:dist_unit,mass_unit
 implicit none

 public :: setpart

 private
 !--private module variables
 !  geometry and boundary parameters
 real :: r_sphere, lbox

 ! physical properties
 real :: density_contrast, totmass_sphere
 real :: angvel, beta_r, rms_mach
 real :: rho_pert_amp
 logical :: angvel_not_betar, shuffle_parts

 ! Bonnor-Ebert sphere parameters
 real :: BErho_cen, BErad_phys, BErad_norm, BEmass, BEfac
 logical :: BEsphere, binary
 integer :: iBEparam

 ! MHD parameters
 real :: Bzero_G, masstoflux, ang_Bomega
 logical :: mu_not_B

 ! sink particle parameters
 real :: rhofinal_setup
 integer :: icreate_sinks_setup
 character(len=20) :: h_acc_char

 ! particle setup parameters
 integer :: np
 character(len=20) :: lattice, cs_sphere_char

contains

!----------------------------------------------------------------
!+
!  setup for a sphere-in-a-box
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:npart_total
 use io,           only:master,fatal
 use eos,          only:polyk2
 use part,         only:Bxyz,igas,idust,set_particle_type
 use set_dust_options, only:dustbinfrac,set_dust_grain_distribution,dtg=>dust_to_gas,ilimitdustfluxinp,&
                            ndustsmallinp,ndustlargeinp,dust_method
 use options,      only:use_dustfrac
 use kernel,       only:hfact_default
 use infile_utils, only:get_options,infile_exists
 use dust,         only:ilimitdustflux
 use units,        only:umass,udist
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer            :: ierr,iBElast,npartsphere
 real               :: totmass,vol_box,vol_sphere,cs_sphere
 real               :: dens_sphere,dens_medium,cs_medium,angvel_code,przero
 real               :: totmass_box,t_ff,area,rmasstoflux_crit,Bzero,h_acc_setup
 real               :: central_density,edge_density
 real, allocatable  :: rtab(:),rhotab(:)

 print "(/,1x,63('-'),1(/,a),/,1x,63('-'),/)",&
   '  Sphere-in-box setup: Almost Archimedes'' greatest achievement.'

 ! set default parameters
 hfact = hfact_default
 time  = 0.0
 call set_defaults()

 ! read/write options to/from .setup file
 call get_options(trim(fileprefix)//'.setup', id==master, ierr, &
                  read_setupfile, write_setupfile, setup_interactive)
 if (ierr /= 0) stop 'please edit .setup file and rerun phantomsetup'

 ! set dust grid
 if (use_dust) then
    if (dust_method==1) use_dustfrac = .true.
    ndustsmall = ndustsmallinp
    ndustlarge = ndustlargeinp
    call set_dust_grain_distribution(ndusttypes,dustbinfrac,grainsize,graindens,udist,umass)
    ilimitdustflux = ilimitdustfluxinp
 endif

 ! setup geometry, boundaries, and physical properties
 call setup_geometry_and_physics(rtab,rhotab,dens_sphere,dens_medium,cs_sphere,cs_medium,&
                                 totmass,totmass_box,t_ff,angvel_code,Bzero,przero,&
                                 central_density,edge_density,area,rmasstoflux_crit,&
                                 polyk,polyk2,gamma,vol_box,vol_sphere,iBElast,h_acc_setup)

 ! setup particles (sphere and medium)
 call setup_particles(id,master,hfact,npart,npartoftype,npartsphere,npart_total,xyzh,vxyzu,massoftype,rtab,rhotab,&
                      dens_sphere,dens_medium,totmass,vol_box,vol_sphere,fileprefix,dtg,iBElast,dust_method,dustbinfrac)

 ! add turbulent velocity field
 if (rms_mach > 0.) call set_turbulent_velocity_field(npart,xyzh,vxyzu,cs_sphere,npartsphere)

 ! add uniform rotation (in velocity and magnetic field)
 call set_rotating_sphere(npart,xyzh,vxyzu,Bxyz,dens_sphere,cs_sphere,cs_medium,angvel_code,Bzero)

 ! amend .in file as necessary
 call setup_runtime_parameters(fileprefix,t_ff,h_acc_setup)

 ! print setup summary
 call print_setup_summary(dens_sphere,dens_medium,cs_sphere,cs_medium,t_ff,angvel_code,Bzero,przero,&
                          totmass,totmass_sphere,central_density,edge_density,area,rmasstoflux_crit,dtg)

end subroutine setpart

!----------------------------------------------------------------
!+
!  Initialize default parameters and read setup file
!+
!----------------------------------------------------------------
subroutine set_defaults()
 use set_dust_options, only:set_dust_default_options,dust_method,ndustsmallinp,ndustlargeinp

 ! default options
 dist_unit = '1.0d16 cm'
 mass_unit = 'solarm'
 lattice = 'closepacked'
 h_acc_char  = '1.0d-2'
 np = 1000000
 BEsphere = .false.
 binary = .false.
 density_contrast = 30.0
 r_sphere = 4.
 lbox = 4.
 totmass_sphere = 1.0
 cs_sphere_char = '21888.0 cm/s' ! cm/s ~ 8K assuming mu = 2.31 & gamma = 5/3
 angvel = 1.77d-13
 angvel_not_betar = .true.
 beta_r = 0.02
 rms_mach = 0.
 shuffle_parts = .false.
 rho_pert_amp = 0.1
 icreate_sinks_setup = 0
 rhofinal_setup = 0.15

 ! MHD defaults
 if (mhd) then
    Bzero_G    = 1.0d-4 ! G
    masstoflux = 5.0
    ang_Bomega = 180.0
    mu_not_B   = .true.
 endif

 !--Initialise dust distribution, if using dust
 if (use_dust) then
    call set_dust_default_options()
    ! override some defaults
    dust_method = 1
    ndustsmallinp = 1
    ndustlargeinp = 0
 endif

end subroutine set_defaults

!----------------------------------------------------------------
!+
!  Setup geometry, boundaries, and physical properties
!+
!----------------------------------------------------------------
subroutine setup_geometry_and_physics(rtab,rhotab,dens_sphere,dens_medium,cs_sphere,cs_medium,&
                                      totmass,totmass_box,t_ff,angvel_code,Bzero,przero,&
                                      central_density,edge_density,area,rmasstoflux_crit,&
                                      polyk,polyk2,gamma,vol_box,vol_sphere,iBElast,h_acc_setup)
 use physcon,        only:pi
 use units,          only:utime,unit_density,unit_Bfield,in_code_units
 use eos_barotropic, only:rhocrit0cgs,drhocrit0
 use eos,            only:gmw
 use io,             only:fatal
 use boundary,       only:set_boundary,dxbound,dybound,dzbound
 use setup_params,   only:rhozero,rmax,ihavesetupB
 use part,           only:Bextx,Bexty,Bextz
 use rho_profile,    only:rho_bonnorebert
 real, intent(out), allocatable :: rtab(:), rhotab(:)
 real, intent(out) :: vol_box,vol_sphere
 real, intent(out) :: dens_sphere, dens_medium, cs_sphere, cs_medium
 real, intent(out) :: totmass, totmass_box, t_ff, angvel_code, Bzero, przero
 real, intent(out) :: central_density, edge_density, area, rmasstoflux_crit
 real, intent(out) :: polyk, polyk2, gamma
 real, intent(out) :: h_acc_setup
 integer, intent(out) :: iBElast
 integer :: ierr, iBE
 real :: rhocritTcgs

 ! convert units of sound speed
 cs_sphere = in_code_units(cs_sphere_char,ierr,unit_type='velocity')
 if (ierr /= 0) call fatal('setup_sphereinbox','Error converting sound speed to code units')

 ! convert accretion radius to code units
 h_acc_setup = in_code_units(h_acc_char,ierr,unit_type='length')
 if (ierr /= 0) call fatal('setup_sphereinbox','Error converting accretion radius to code units')

 ! Bonnor-Ebert profile (if requested)
 if (BEsphere .and. density_contrast > 1.) then
    iBE = 8192
    allocate(rtab(iBE),rhotab(iBE))
    call rho_bonnorebert(iBEparam,BErho_cen,edge_density,BErad_phys,BErad_norm,BEmass,BEfac,cs_sphere, &
                         gmw,iBE,iBElast,rtab,rhotab,ierr)
    central_density = BErho_cen
    r_sphere        = BErad_phys
    totmass_sphere  = BEmass
    if (ierr > 0) call fatal('setup_sphereinbox','Error in calculating Bonnor-Ebert profile')
 endif

 ! boundaries
 if (density_contrast < 1.+epsilon(density_contrast)) lbox = 2.
 call set_boundary(l=lbox*r_sphere)

 ! general parameters
 gamma       = 5./3.
 if (isothermal) gamma = 1.
 rmax        = r_sphere
 if (angvel_not_betar) then
    angvel_code = angvel*utime
 else
    angvel_code = sqrt(3.0*totmass_sphere*beta_r/r_sphere**3)
    angvel      = angvel_code/utime
 endif

 vol_box     = dxbound*dybound*dzbound
 vol_sphere  = 4./3.*pi*r_sphere**3
 rhozero     = totmass_sphere / vol_sphere
 dens_sphere = rhozero

 if (BEsphere) then
    dens_medium = edge_density/density_contrast
    cs_medium   = cs_sphere*sqrt(edge_density/dens_medium)
 else
    dens_medium = dens_sphere/density_contrast
    cs_medium   = cs_sphere*sqrt(density_contrast)
 endif
 rhocrit0cgs = dens_medium*unit_density * density_contrast/7.5  ! end of transition region for ieos=8;for density_contrast=30, this yields a coefficient of 4
 rhocritTcgs = rhocrit0cgs*(1.0-drhocrit0)                      ! start of transition region for ieos=8
 totmass_box = (vol_box - vol_sphere)*dens_medium
 totmass     = totmass_box + totmass_sphere
 t_ff        = sqrt(3.*pi/(32.*dens_sphere))

 if (totmass_sphere > 0.9*totmass .and. density_contrast > 1.) then
    print*, 'resetting boundaries to increase the number of background particles'
    dxbound = (0.1*totmass_sphere/dens_medium)**(1./3.)
    call set_boundary(l=dxbound)
    vol_box     = dxbound*dybound*dzbound
    totmass_box = (vol_box - vol_sphere)*dens_medium
    totmass     = totmass_box + totmass_sphere
 endif

 if (dens_medium*unit_density > rhocritTcgs .and. density_contrast > 1.0) then
    print*, 'Medium density = ',dens_medium*unit_density,'g/cm^3'
    print*, 'Sphere density = ',dens_sphere*unit_density,'g/cm^3'
    print*, 'Density at start of EOS transition = ',rhocritTcgs,'g/cm^3'
    print*, 'Density at end   of EOS transition = ',rhocrit0cgs,'g/cm^3'
    call fatal('setup_sphereinbox','Error setting sound-speed transition region in EOS')
 endif
 !
 ! magnetic field
 !
 rmasstoflux_crit = 2./3.*0.53*sqrt(5./pi)
 if (mhd) then
    area = pi*r_sphere**2
    if (mu_not_B) then
       if (masstoflux > tiny(masstoflux)) then
          Bzero = totmass_sphere/(area*masstoflux*rmasstoflux_crit)
       else
          Bzero = 0.
       endif
    else
       Bzero      = Bzero_G/unit_Bfield
       masstoflux = totmass_sphere/(area*Bzero*rmasstoflux_crit)
    endif
    ihavesetupB = .true.
 else
    Bzero = 0.
 endif

 Bextx  = 0.
 Bexty  = 0.
 Bextz  = Bzero
 przero = cs_sphere**2*dens_sphere

 ! temperature set to give a pressure equilibrium
 polyk  = cs_sphere**2
 polyk2 = cs_medium**2

end subroutine setup_geometry_and_physics

!----------------------------------------------------------------
!+
!  Setup particles (sphere and medium)
!+
!----------------------------------------------------------------
subroutine setup_particles(id,master,hfact,npart,npartoftype,npartsphere,npart_total,xyzh,vxyzu,massoftype,rtab,rhotab,&
                           dens_sphere,dens_medium,totmass,vol_box,vol_sphere,fileprefix,dtg,iBElast,dust_method,&
                           dustbinfrac)
 use unifdis,                only:set_unifdis
 use spherical,              only:set_sphere
 use mpidomain,              only:i_belong
 use part,                   only:set_particle_type,igas,idust,dustfrac,ndusttypes
 use utils_shuffleparticles, only:shuffleparticles
 use centreofmass,           only:reset_centreofmass
 use io,                     only:iprint
 use boundary,               only:xmin,xmax,ymin,ymax,zmin,zmax,dxbound
 use setup_params,           only:rmax
 integer,           intent(in)   :: id,master
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:),npartsphere
 integer(kind=8),   intent(out)   :: npart_total
 integer,           intent(in)    :: iBElast,dust_method
 real,              intent(in)    :: dustbinfrac(:)
 real,              intent(in)    :: hfact
 real,              intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(in)    :: dens_sphere,dens_medium,totmass,vol_box,vol_sphere,dtg
 character(len=20), intent(in)    :: fileprefix
 real,              intent(inout), allocatable :: rtab(:), rhotab(:)
 integer :: i,np_in,ierr
 real :: psep,psep_box,pmass_dusttogas

 np_in = np

 ! setup particles in the sphere
 if (density_contrast > 1.) then
    if (BEsphere) then
       call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh, &
                       rhotab=rhotab(1:iBElast),rtab=rtab(1:iBElast),nptot=npart_total,&
                       exactN=.true.,np_requested=np,mask=i_belong)
    else
       call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh,nptot=npart_total,&
                       exactN=.true.,np_requested=np,mask=i_belong)
       if (trim(lattice)/='random') print "(a,es10.3)",' Particle separation in sphere = ',psep
    endif
    print "(a)",' Initialised sphere'
 else
    psep_box = dxbound/np**(1./3.)
    call set_unifdis(trim(lattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep_box, &
                      hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong,err=ierr)
 endif

 npartsphere = npart
 if (np_in /= npartsphere) np = npartsphere
 massoftype(igas) = totmass_sphere/npartsphere

 ! setup surrounding low density medium
 if (BEsphere .or. trim(lattice)=='random') then
    psep_box = dxbound/(vol_box*dens_medium/massoftype(igas))**(1./3.)
 else
    psep_box = psep*(density_contrast)**(1./3.)
 endif

 if (density_contrast > 1.0) then
    call set_unifdis(trim(lattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep_box, &
                      hfact,npart,xyzh,periodic,rmin=r_sphere,nptot=npart_total,mask=i_belong,err=ierr)
    print "(a,es10.3)",' Particle separation in low density medium = ',psep_box
    print "(a,i10,a)",' added ',npart-npartsphere,' particles in low-density medium'
    print*, ""
 endif

 ! set particle properties
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 dustfrac          = 0.
 if (massoftype(igas) < epsilon(massoftype(igas))) massoftype(igas) = totmass/npart_total

 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
    if (use_dust .and. use_dustfrac) then
       if (ndusttypes > 1) then
          dustfrac(1:ndusttypes,i) = dustbinfrac(1:ndusttypes)*dtg
       else
          dustfrac(1,i) = dtg/(1.+dtg)
       endif
    endif
 enddo
 !
 ! Set dust-as-particles (experimental)
 !
 if (use_dust .and. dust_method==2) then
    ! particle separation in dust sphere & sdjust for close-packed lattice
    pmass_dusttogas = 10.*dtg*massoftype(igas)
    psep = (vol_sphere/pmass_dusttogas/real(np))**(1./3.)
    psep = psep*sqrt(2.)**(1./3.)
    call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh,nptot=npart_total,&
                    exactN=.true.,np_requested=np/10,mask=i_belong)
    npartoftype(idust) = int(npart_total) - npartoftype(igas)
    massoftype(idust)  = totmass_sphere*dtg/npartoftype(idust)

    do i = npartoftype(igas)+1,npart
       call set_particle_type(i,idust)
    enddo

    print "(a,4(i10,1x))", ' particle numbers: (gas_total, gas_sphere, dust, total): ' &
                        , npartoftype(igas),npartsphere,npartoftype(idust),npart
    print "(a,2es10.3)"  , ' particle masses: (gas,dust): ',massoftype(igas),massoftype(idust)
 else
    print "(a,3(i10,1x))", ' particle numbers: (sphere, low-density medium, total): ' &
                        , npartsphere, npart-npartsphere,npart
    print "(a,es10.3)",' particle mass = ',massoftype(igas)
 endif
 !
 ! shuffle particles
 !
 if (shuffle_parts) then
    print*, "lets shuffle!"
    if (BEsphere) then
       call shuffleparticles(iprint,npart,xyzh,massoftype(igas),dmedium=dens_medium,ntab=iBElast, &
                             rtab=rtab,dtab=rhotab,dcontrast=density_contrast,is_setup=.true.,prefix=trim(fileprefix))
    else
       call shuffleparticles(iprint,npart,xyzh,massoftype(igas), &
                             rsphere=rmax,dsphere=dens_sphere,dmedium=dens_medium,is_setup=.true.,prefix=trim(fileprefix))
    endif
 endif

 ! reset to centre of mass
 if (trim(lattice)/='random' .and. .not.shuffle_parts) call reset_centreofmass(npart,xyzh,vxyzu)

 ! add binary perturbation if requested
 if (binary) call set_binary_perturbation(npart,xyzh,rho_pert_amp)

 if (allocated(rtab)) deallocate(rtab)
 if (allocated(rhotab)) deallocate(rhotab)

end subroutine setup_particles

!----------------------------------------------------------------
!+
!  add m=2 azimuthal perturbation to the density profile
!+
!----------------------------------------------------------------
subroutine set_binary_perturbation(npart,xyzh,amplitude)
 use physcon, only:pi
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: amplitude
 integer :: i
 real :: rxy2, rxyz2, phi, dphi

 !--Stretching the spatial distribution to perturb the density profile
 do i = 1,npart
    rxy2  = xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i)
    rxyz2 = rxy2 + xyzh(3,i)*xyzh(3,i)
    if (rxyz2 <= r_sphere**2) then
       phi       = atan(xyzh(2,i)/xyzh(1,i))
       if (xyzh(1,i) < 0.0) phi = phi + pi
       dphi      = 0.5*amplitude*sin(2.0*phi)
       phi       = phi - dphi
       xyzh(1,i) = sqrt(rxy2)*cos(phi)
       xyzh(2,i) = sqrt(rxy2)*sin(phi)
    endif
 enddo

end subroutine set_binary_perturbation

!----------------------------------------------------------------
!+
!  set turbulent velocity field by reading from cubes
!+
!----------------------------------------------------------------
subroutine set_turbulent_velocity_field(npart,xyzh,vxyzu,cs_sphere,npartsphere)
 use velfield,  only:set_velfield_from_cubes
 use datafiles, only:find_phantom_datafile
 use io,        only:fatal
 use boundary,  only:xmax
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(out)   :: vxyzu(:,:)
 real,    intent(in)    :: cs_sphere
 integer, intent(inout) :: npartsphere
 integer :: i,ierr
 real :: v2i,rmsmach,turbfac
 character(len=120) :: filex,filey,filez
 character(len=20), parameter :: filevx = 'cube_v1.dat'
 character(len=20), parameter :: filevy = 'cube_v2.dat'
 character(len=20), parameter :: filevz = 'cube_v3.dat'

 ! reset values if turbulence-in-a-box
 if (density_contrast < 1.+epsilon(density_contrast)) then
    r_sphere = xmax
    npartsphere = npart
 endif

 ! Velocity: Turbulent velocity field
 vxyzu = 0.
 filex  = find_phantom_datafile(filevx,'velfield')
 filey  = find_phantom_datafile(filevy,'velfield')
 filez  = find_phantom_datafile(filevz,'velfield')

 call set_velfield_from_cubes(xyzh(:,1:npartsphere),vxyzu(:,:npartsphere),npartsphere, &
                              filex,filey,filez,1.,r_sphere,.false.,ierr)
 if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

 rmsmach = 0.0
 print*, 'Turbulence being set by user'
 do i = 1,npartsphere
    v2i     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
    rmsmach = rmsmach + v2i/cs_sphere**2
 enddo
 rmsmach = sqrt(rmsmach/npartsphere)
 if (rmsmach > 0.) then
    turbfac = rms_mach/rmsmach ! normalise the energy to the desired mach number
 else
    turbfac = 0.
 endif
 do i = 1,npartsphere
    vxyzu(1:3,i) = turbfac*vxyzu(1:3,i)
 enddo

end subroutine set_turbulent_velocity_field

!----------------------------------------------------------------
!+
!  set uniform rotation (thermal energy & magnetic field too)
!+
!----------------------------------------------------------------
subroutine set_rotating_sphere(npart,xyzh,vxyzu,Bxyz,dens_sphere,cs_sphere,cs_medium,angvel_code,Bzero)
 use physcon, only:pi
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(inout) :: Bxyz(:,:)
 real,    intent(in)    :: dens_sphere,cs_sphere,cs_medium,Bzero
 real,    intent(in)    :: angvel_code
 integer :: i
 real :: r2

 do i=1,npart
    r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
    if (r2 < r_sphere**2) then
       vxyzu(1,i) = vxyzu(1,i) - angvel_code*xyzh(2,i)
       vxyzu(2,i) = vxyzu(2,i) + angvel_code*xyzh(1,i)
       if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*cs_sphere**2
    else
       if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*cs_medium**2
    endif
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(1,i) = Bzero*sin(ang_Bomega*pi/180.0)
       Bxyz(3,i) = Bzero*cos(ang_Bomega*pi/180.0)
    endif
 enddo

end subroutine set_rotating_sphere

!----------------------------------------------------------------
!+
!  Setup runtime parameters
!+
!----------------------------------------------------------------
subroutine setup_runtime_parameters(fileprefix,t_ff,h_acc_setup)
 use timestep,     only:dtmax,tmax,dtmax_dratio,dtmax_min
 use options,      only:nfulldump,rhofinal_cgs,icooling,calc_erot
 use ptmass,       only:icreate_sinks,h_acc,r_crit
 use eos,          only:ieos
 use infile_utils, only:infile_exists
 character(len=20), intent(in) :: fileprefix
 real, intent(in) :: t_ff
 real, intent(in) :: h_acc_setup
 ! set default runtime parameters if .in file does not exist
 !
 dtmax = t_ff/100.  ! Since this variable can change, always reset it if running phantomsetup
 if (.not. infile_exists(fileprefix)) then
    ! default values
    tmax          = 1.21*t_ff ! = 10.75 for default settings (Wurster, Price & Bate 2016)
    if (isothermal) ieos = 8
    nfulldump     = 1
    calc_erot     = .true.
    icreate_sinks = icreate_sinks_setup
    h_acc         = h_acc_setup
    r_crit        = 5.0*h_acc
    ! reset defaults based upon options
    if (density_contrast > 1.) dtmax_dratio = 1.258
    if (density_contrast < 1.+epsilon(density_contrast) .and. maxvxyzu>=4) then
       ieos     = 2
       icooling = 5
    endif
    if (binary) tmax = 1.50*t_ff ! = 13.33 for default settings (Wurster, Price & Bate 2017)
    if (icreate_sinks==1) then
       dtmax_min = dtmax/8.0
    else
       dtmax_min    = 0.0
       rhofinal_cgs = rhofinal_setup
    endif
 endif

end subroutine setup_runtime_parameters

!----------------------------------------------------------------
!+
!  Print setup summary
!+
!----------------------------------------------------------------
subroutine print_setup_summary(dens_sphere,dens_medium,cs_sphere,cs_medium,t_ff,angvel_code,Bzero,przero,&
                               totmass,totmass_sphere,central_density,edge_density,area,rmasstoflux_crit,dtg)
 use units, only:in_units,unit_Bfield
 real, intent(in) :: dens_sphere, dens_medium, cs_sphere, cs_medium
 real, intent(in) :: t_ff, angvel_code, Bzero, przero
 real, intent(in) :: totmass, totmass_sphere, central_density, edge_density, area, rmasstoflux_crit
 real, intent(in) :: dtg
 character(len=40) :: fmt

 !--Summarise the sphere
 print "(a,i10)",' Input npart_sphere = ',np
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass,in_units(totmass,'g'),' g'
 print fmt,' Mass in sphere   : ',totmass_sphere,in_units(totmass_sphere,'g'),' g'
 print fmt,' Radius of sphere : ',r_sphere,in_units(r_sphere,'cm'),' cm'
 if (BEsphere) then
    print fmt,' Mean rho sphere  : ',dens_sphere,in_units(dens_sphere,'g/cm^3'),' g/cm^3'
    print fmt,' central density  : ',central_density,in_units(central_density,'g/cm^3'),' g/cm^3'
    print fmt,' edge density     : ',edge_density,in_units(edge_density,'g/cm^3'),' g/cm^3'
    print fmt,' Mean rho medium  : ',dens_medium,in_units(dens_medium,'g/cm^3'),' g/cm^3'
 else
    print fmt,' Density sphere   : ',dens_sphere,in_units(dens_sphere,'g/cm^3'),' g/cm^3'
    print fmt,' Density medium   : ',dens_medium,in_units(dens_medium,'g/cm^3'),' g/cm^3'
 endif
 print fmt,' cs in sphere     : ',cs_sphere,in_units(cs_sphere,'cm/s'),' cm/s'
 print fmt,' cs in medium     : ',cs_medium,in_units(cs_medium,'cm/s'),' cm/s'
 print fmt,' Free fall time   : ',t_ff,in_units(t_ff,'yrs'),' yrs'
 print fmt,' Angular velocity : ',angvel_code,angvel,' rad/s'
 print fmt,' Turbulent Mach no: ',rms_mach
 print fmt,' Omega*t_ff       : ',angvel_code*t_ff
 if (mhd) then
    print fmt,' B field (z)      : ',Bzero,Bzero*unit_Bfield*1.d6,' micro-G'
    print fmt,' Alfven speed     : ',Bzero/sqrt(dens_sphere),in_units(Bzero/sqrt(dens_sphere),'cm/s'),' cm/s'
    if (Bzero > 0.) then
       print fmt,' plasma beta      : ',przero/(0.5*Bzero*Bzero)
       print fmt,' mass-to-flux     : ',totmass_sphere/(area*Bzero)/rmasstoflux_crit
    endif
 endif
 if (use_dust) then
    print fmt,' dust-to-gas ratio: ',dtg,dtg,' '
 endif
 print "(1x,50('-'))"

end subroutine print_setup_summary

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils,     only:write_inopt
 use setunits,         only:write_options_units
 use set_dust_options, only:write_dust_setup_options
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for sphere-in-box setup routines'

 ! units
 call write_options_units(iunit,gr)

 write(iunit,"(/,a)") '# resolution and particle placement'
 call write_inopt(np,'np','requested number of particles in sphere',iunit)
 call write_inopt(lattice,'lattice','particle lattice (random,cubic,closepacked,hcp,hexagonal)',iunit)
 call write_inopt(shuffle_parts,'shuffle_parts','relax particles by shuffling',iunit)

 write(iunit,"(/,a)") '# options for box'
 call write_inopt(lbox,'lbox','length of a box side in terms of spherical radii',iunit)

 write(iunit,"(/,a)") '# intended result'
 call write_inopt(binary,'form_binary','the intent is to form a central binary',iunit)

 write(iunit,"(/,a)") '# options for sphere'
 call write_inopt(BEsphere,'use_BE_sphere','centrally condense as a BE sphere',iunit)
 if (.not. BEsphere) then
    call write_inopt(r_sphere,'r_sphere','radius of sphere in code units',iunit)
    call write_inopt(totmass_sphere,'totmass_sphere','mass of sphere in code units',iunit)
 else
    call write_inopt(iBEparam,'iBE_options','The set of parameters to define the BE sphere',iunit)
    if (iBEparam==1 .or. iBEparam==2 .or. iBEparam==3) &
       call write_inopt(BErho_cen,'BErho_cen','central density of the BE sphere [code units]',iunit)
    if (iBEparam==1 .or. iBEparam==4 .or. iBEparam==6) &
        call write_inopt(BErad_phys,'BErad_phys','physical radius of the BE sphere [code units]',iunit)
    if (iBEparam==2 .or. iBEparam==4 .or. iBEparam==5) &
        call write_inopt(BErad_norm,'BErad_norm','normalised radius of the BE sphere',iunit)
    if (iBEparam==3 .or. iBEparam==5 .or. iBEparam==6) &
        call write_inopt(BEmass,'BEmass','mass radius of the BE sphere [code units]',iunit)
    if (iBEparam==4 .or. iBEparam==5)                  &
        call write_inopt(BEfac,'BEfac','over-density factor of the BE sphere [code units]',iunit)
 endif
 call write_inopt(density_contrast,'density_contrast','density contrast in code units',iunit)
 call write_inopt(cs_sphere_char,'cs_sphere','sound speed in sphere (with units e.g. 0.2 km/s)',iunit)
 if (angvel_not_betar) then
    call write_inopt(angvel,'angvel','angular velocity in rad/s',iunit)
 else
    call write_inopt(beta_r,'beta_r','rotational-to-gravitational energy ratio',iunit)
 endif
 call write_inopt(rms_mach,'rms_mach','turbulent rms mach number',iunit)
 if (mhd) then
    if (mu_not_B) then
       call write_inopt(masstoflux,'masstoflux','mass-to-magnetic flux ratio in units of critical value',iunit)
    else
       call write_inopt(Bzero_G,'Bzero','Magnetic field strength in Gauss',iunit)
    endif
    call write_inopt(ang_Bomega,'ang_Bomega','Angle (degrees) between B and rotation axis',iunit)
 endif
 if (use_dust) call write_dust_setup_options(iunit)
 if (binary) call write_inopt(rho_pert_amp,'rho_pert_amp','amplitude of density perturbation',iunit)

 write(iunit,"(/,a)") '# Sink properties (values in .in file, if present, will take precedence)'
 call write_inopt(icreate_sinks_setup,'icreate_sinks','1: create sinks.  0: do not create sinks',iunit)
 if (icreate_sinks_setup==1) then
    call write_inopt(h_acc_char,'h_acc','accretion radius (with units; e.g. au,pc,kpc,0.1pc)',iunit)
 else
    call write_inopt(rhofinal_setup,'rho_final','final maximum density (<=0 to ignore) (cgs units)',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils,     only:open_db_from_file,inopts,read_inopt,close_db
 use unifdis,          only:is_valid_lattice
 use setunits,         only:read_options_and_set_units
 use set_dust_options, only:read_dust_setup_options
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: nerr,jerr,kerr
 type(inopts), allocatable     :: db(:)

 !--Read values
 print "(a)",' reading setup options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr,gr)
 call read_inopt(BEsphere,'use_BE_sphere',db,errcount=nerr)
 call read_inopt(binary,'form_binary',db,errcount=nerr)
 call read_inopt(np,'np',db,errcount=nerr)
 call read_inopt(lattice,'lattice',db,ierr,errcount=nerr)
 if (ierr /= 0 .or. .not. is_valid_lattice(lattice)) then
    print*, ' invalid lattice.  Setting to closepacked'
    lattice = 'closepacked'
 endif
 call read_inopt(shuffle_parts,'shuffle_parts',db,errcount=nerr)
 call read_inopt(lbox,'lbox',db,errcount=nerr)

 if (.not. BEsphere) then
    call read_inopt(r_sphere,'r_sphere',db,errcount=nerr)
    call read_inopt(totmass_sphere,'totmass_sphere',db,errcount=nerr)
 else
    call read_inopt(iBEparam,'iBE_options',db,errcount=nerr)
    if (iBEparam==1 .or. iBEparam==2 .or. iBEparam==3) call read_inopt(BErho_cen,'BErho_cen',db,errcount=nerr)
    if (iBEparam==1 .or. iBEparam==4 .or. iBEparam==6) call read_inopt(BErad_phys,'BErad_phys',db,errcount=nerr)
    if (iBEparam==2 .or. iBEparam==4 .or. iBEparam==5) call read_inopt(BErad_norm,'BErad_norm',db,errcount=nerr)
    if (iBEparam==3 .or. iBEparam==5 .or. iBEparam==6) call read_inopt(BEmass,'BEmass',db,errcount=nerr)
    if (iBEparam==4 .or. iBEparam==5)                  call read_inopt(BEfac,'BEfac',db,errcount=nerr)
 endif

 call read_inopt(density_contrast,'density_contrast',db,errcount=nerr)
 call read_inopt(cs_sphere_char,'cs_sphere',db,errcount=nerr)
 call read_inopt(angvel,'angvel',db,jerr)
 call read_inopt(beta_r,'beta_r',db,kerr)
 angvel_not_betar = .true.
 if (jerr /= 0 .and. kerr == 0) then
    angvel_not_betar = .false.
 elseif (jerr == 0 .and. kerr /= 0) then
    angvel_not_betar = .true.
 else
    ierr = ierr + 1
 endif
 call read_inopt(rms_mach,'rms_mach',db,ierr)
 mu_not_B = .true.
 if (mhd) then
    call read_inopt(masstoflux,'masstoflux',db,jerr)
    call read_inopt(Bzero_G,   'Bzero',     db,kerr)
    call read_inopt(ang_Bomega,'ang_Bomega',db,ierr)
    if (jerr /= 0 .and. kerr == 0) then
       mu_not_B = .false.
    elseif (jerr == 0 .and. kerr /= 0) then
       mu_not_B = .true.
    else
       nerr = nerr + 1
    endif
 endif
 if (use_dust) call read_dust_setup_options(db,nerr)
 if (binary) call read_inopt(rho_pert_amp,'rho_pert_amp',db,errcount=nerr)
 call read_inopt(icreate_sinks_setup,'icreate_sinks',db,errcount=nerr)
 if (icreate_sinks_setup==1) then
    call read_inopt(h_acc_char,'h_acc',db,errcount=nerr)
 else
    call read_inopt(rhofinal_setup,'rho_final',db,errcount=nerr)
 endif
 call close_db(db)

 if (nerr > 0) then
    ierr = nerr
    print "(1x,a,i2,a)",'Setup_sphereinbox: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Interactive setup routine
!+
!-----------------------------------------------------------------------
subroutine setup_interactive()
 use prompting,        only:prompt
 use rho_profile,      only:prompt_BEparameters
 use part,             only:maxp
 use setunits,         only:set_units_interactive
 use unifdis,          only:ilattice_max,i_closepacked,get_latticetype,latticetype
 use infile_utils,     only:get_optstring
 use units,            only:udist,umass
 use physcon,          only:au,solarm
 use set_dust_options, only:set_dust_interactive
 integer :: ilattice,npmax
 character(len=100) :: string
 logical :: make_sinks

 ! Set default units

 call set_units_interactive(gr)

 ! Prompt user for settings
 npmax = int(2.0/3.0*maxp) ! approx max number allowed in sphere given size(xyzh(1,:))
 if (npmax < 300000) then
    np = npmax
 elseif (npmax < np) then
    np = 300000
 endif
 call prompt('Enter the approximate number of particles in the sphere',np,0,npmax)

 ilattice = i_closepacked
 call get_optstring(ilattice_max,latticetype,string,4)
 call prompt('Enter the type of particle lattice '//trim(string),ilattice,0,ilattice_max)
 lattice = get_latticetype(ilattice)

 shuffle_parts = .false.
 if (ilattice==1) shuffle_parts = .true.
 call prompt('Relax particles by shuffling?',shuffle_parts)

 call prompt('Centrally condense the sphere as a BE sphere?',BEsphere)

 if (.not. BEsphere) then
    call prompt('Enter radius of sphere in units of '//dist_unit,r_sphere,0.)
    call prompt('Enter the box size as a multiple of the sphere radius: ',lbox,1.)
    call prompt('Enter total mass in sphere in units of '//mass_unit,totmass_sphere,0.)
 else
    call prompt_BEparameters(iBEparam,BErho_cen,BErad_phys,BErad_norm,BEmass,BEfac,umass,udist,au,solarm)
    call prompt('Enter the box size as a multiple of the sphere radius: ',lbox,1.)
 endif

 call prompt('Enter density contrast between sphere and box ',density_contrast,1.)

 call prompt('Do you intend to form a binary system (i.e. add an m=2 perturbation)?',binary)

 if (binary) cs_sphere_char = '18696.96 cm/s' ! cm/s ~ 5K assuming mu = 2.31 & gamma = 5/3
 call prompt('Enter sound speed in sphere (with units e.g. 0.2 km/s)',cs_sphere_char)

 if (binary) angvel = 1.006d-12
 call prompt('Input angular velocity (true); else input ratio of rotational-to-potential energy ',angvel_not_betar)
 if (angvel_not_betar) then
    call prompt('Enter angular rotation speed in rad/s ',angvel,0.)
 else
    call prompt('Enter ratio of rotational-to-potential energy ',beta_r,0.)
 endif

 call prompt('Enter the Mach number of the cloud turbulence',rms_mach,0.)

 if (mhd) then
    call prompt('Input the mass-to-flux ratio (true); else input the magnetic field strength ',mu_not_B)
    if (mu_not_B) then
       call prompt('Enter mass-to-flux ratio in units of critical value ',masstoflux,0.)
    else
       call prompt('Enter magnetic field strength in Gauss ',Bzero_G,0.)
    endif
    call prompt('Enter the angle (degrees) between B and the rotation axis? ',ang_Bomega)
 endif

 if (use_dust) call set_dust_interactive()

 if (binary) call prompt('Enter the amplitute of the density perturbation ',rho_pert_amp,0.0,0.4)

 ! ask about sink particle details; these will not be saved to the .setup file since they exist in the .in file
 call prompt('Dynamically create sink particles? ',make_sinks)
 if (make_sinks) then
    if (binary) h_acc_char  = '3.35au'
    call prompt('Enter the accretion radius of the sink (with units; e.g. au,pc,kpc,0.1pc) ',h_acc_char)
    icreate_sinks_setup = 1
 else
    icreate_sinks_setup = 0
    rhofinal_setup = 0.15
    call prompt('Enter final maximum density in g/cm^3 (ignored for <= 0) ',rhofinal_setup)
 endif

end subroutine setup_interactive

end module setup
