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
! :Owner: James Wurster
!
! :Runtime parameters:
!   - BEfac             : *over-density factor of the BE sphere [code units]*
!   - BEmass            : *mass radius of the BE sphere [code units]*
!   - BErad_norm        : *normalised radius of the BE sphere*
!   - BErad_phys        : *physical radius of the BE sphere [code units]*
!   - BErho_cen         : *central density of the BE sphere [code units]*
!   - Bzero             : *Magnetic field strength in Gauss*
!   - ang_Bomega        : *Angle (degrees) between B and rotation axis*
!   - angvel            : *angular velocity in rad/s*
!   - beta_r            : *rotational-to-gravitational energy ratio*
!   - cs_sphere_cgs     : *sound speed in sphere in cm/s*
!   - density_contrast  : *density contrast in code units*
!   - dist_unit         : *distance unit (e.g. au)*
!   - dust_to_gas_ratio : *dust-to-gas ratio*
!   - form_binary       : *the intent is to form a central binary*
!   - graindenscgs      : *grain density [g/cm^3]*
!   - grainsizecgs      : *grain size in [cm]*
!   - h_acc             : *accretion radius (code units)*
!   - h_soft_sinksink   : *sink-sink softening radius (code units)*
!   - iBE_options       : *The set of parameters to define the BE sphere*
!   - icreate_sinks     : *1: create sinks.  0: do not create sinks*
!   - lattice           : *particle lattice (random,cubic,closepacked,hcp,hexagonal)*
!   - lbox              : *length of a box side in terms of spherical radii*
!   - mass_unit         : *mass unit (e.g. solarm)*
!   - masstoflux        : *mass-to-magnetic flux ratio in units of critical value*
!   - ndusttypes        : *number of grain sizes*
!   - np                : *requested number of particles in sphere*
!   - r_crit            : *critical radius (code units)*
!   - r_sphere          : *radius of sphere in code units*
!   - rho_final         : *final maximum density (<=0 to ignore) (cgs units)*
!   - rho_pert_amp      : *amplitude of density perturbation*
!   - rms_mach          : *turbulent rms mach number*
!   - shuffle_parts     : *relax particles by shuffling*
!   - sindex            : *power-law index, e.g. MRN*
!   - smaxcgs           : *maximum grain size [cm]*
!   - smincgs           : *minimum grain size [cm]*
!   - totmass_sphere    : *mass of sphere in code units*
!   - use_BE_sphere     : *centrally condense as a BE sphere*
!
! :Dependencies: boundary, centreofmass, datafiles, dim, dust, eos,
!   eos_barotropic, infile_utils, io, kernel, mpidomain, options, part,
!   physcon, prompting, ptmass, rho_profile, set_dust, set_dust_options,
!   setup_params, spherical, timestep, unifdis, units,
!   utils_shuffleparticles, velfield
!
 use part,             only:mhd,graindens,grainsize,ndusttypes,ndustsmall
 use dim,              only:use_dust,maxvxyzu,periodic,maxdustsmall
 use options,          only:calc_erot
 use dust,             only:grainsizecgs,graindenscgs
 use set_dust_options, only:grainsizeinp,graindensinp,igrainsize,igraindens,&
                            smincgs,smaxcgs,sindex,dustbinfrac
 implicit none

 public :: setpart

 private
 !--private module variables
 real :: xmini(3), xmaxi(3)
 real :: density_contrast,totmass_sphere,r_sphere,cs_sphere,cs_sphere_cgs
 real :: angvel,beta_r,Bzero_G,masstoflux,dtg,ang_Bomega,rms_mach
 real :: rho_pert_amp,lbox
 real :: BErho_cen,BErad_phys,BErad_norm,BEmass,BEfac
 real :: r_crit_setup,h_acc_setup,h_soft_sinksink_setup,rhofinal_setup
 real(kind=8)                 :: udist,umass
 integer                      :: np,iBEparam,icreate_sinks_setup
 logical                      :: BEsphere,binary,mu_not_B,cs_in_code,angvel_not_betar,shuffle_parts
 logical                      :: is_cube = .true. ! if false, then can set a rectangle if BEsphere=false; for backwards compatibility
 character(len=20)            :: dist_unit,mass_unit,lattice
 character(len= 1), parameter :: labelx(3) = (/'x','y','z'/)

contains

!----------------------------------------------------------------
!+
!  setup for a sphere-in-a-box
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact_out,time,fileprefix)
 use physcon,      only:pi,solarm,hours,years,au
 use dim,          only:maxdusttypes,use_dustgrowth,maxdustlarge
 use setup_params, only:rhozero,npart_total,rmax,ihavesetupB
 use io,           only:master,fatal,iprint
 use unifdis,      only:set_unifdis
 use spherical,    only:set_sphere
 use rho_profile,  only:rho_bonnorebert,prompt_BEparameters
 use boundary,     only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,unit_density,unit_Bfield,unit_velocity
 use eos,          only:polyk2,ieos,gmw
 use eos_barotropic, only:rhocrit0cgs,drhocrit0
 use part,         only:Bxyz,Bextx,Bexty,Bextz,igas,idust,set_particle_type,hfact,dustfrac
 use set_dust_options, only:dustbinfrac,set_dust_default_options,set_dust_interactively,dust_method
 use dust,         only:ilimitdustflux
 use timestep,     only:dtmax,tmax,dtmax_dratio,dtmax_min
 use centreofmass, only:reset_centreofmass
 use options,      only:nfulldump,rhofinal_cgs,use_dustfrac,icooling
 use kernel,       only:hfact_default
 use mpidomain,    only:i_belong
 use ptmass,       only:icreate_sinks,r_crit,h_acc,h_soft_sinksink
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use set_dust,     only:set_dustfrac,set_dustbinfrac
 use utils_shuffleparticles, only:shuffleparticles
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact_out
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=20), parameter     :: filevx = 'cube_v1.dat'
 character(len=20), parameter     :: filevy = 'cube_v2.dat'
 character(len=20), parameter     :: filevz = 'cube_v3.dat'
 real(kind=8)       :: h_acc_in
 integer            :: i,np_in,npartsphere,npmax,iBElast,ierr
 integer            :: iBE,ilattice
 real               :: totmass,vol_box,psep,psep_box,pmass_dusttogas
 real               :: vol_sphere,dens_sphere,dens_medium,cs_medium,angvel_code,przero
 real               :: totmass_box,t_ff,r2,area,Bzero,rmasstoflux_crit
 real               :: rxy2,rxyz2,phi,dphi,central_density,edge_density,rmsmach,v2i,turbfac,rhocritTcgs
 real, allocatable  :: rtab(:),rhotab(:)
 logical            :: iexist
 logical            :: make_sinks = .true.     ! the default prompt is to ask to make sinks
 character(len=120) :: filex,filey,filez
 character(len=100) :: filename,cwd
 character(len=40)  :: fmt
 character(len=10)  :: h_acc_char

 !--Initialise dust distribution, if using dust
 if (use_dust) call set_dust_default_options()

 filename = trim(fileprefix)//'.setup'
 print "(/,1x,63('-'),1(/,a),/,1x,63('-'),/)",&
   '  Sphere-in-box setup: Almost Archimedes'' greatest achievement.'

 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    np_in = np
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    dist_unit = '1.0d16 cm'
    mass_unit = 'solarm'
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
       call select_unit(mass_unit,umass,ierr)
       if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    enddo
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1 pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo
    !
    ! units
    !
    call set_units(dist=udist,mass=umass,G=1.d0)
    !
    ! prompt user for settings
    !
    npmax = int(2.0/3.0*size(xyzh(1,:))) ! approx max number allowed in sphere given size(xyzh(1,:))
    if (npmax < 300000) then
       np = npmax
    elseif (npmax < 1000000) then
       np = 300000
    else
       np = 1000000
    endif
    call prompt('Enter the approximate number of particles in the sphere',np,0,npmax)
    np_in    = np

    lattice  = 'closepacked'
    ilattice = 3
    call prompt('Enter the type of particle lattice (1=random,2=cubic,3=closepacked,4=hexagonal)',ilattice,0,4)
    if (ilattice==1) then
       lattice = 'random'
       shuffle_parts = .false.
    elseif (ilattice==2) then
       lattice = 'cubic'
    elseif (ilattice==4) then
       lattice = 'hexagonal'
    endif

    shuffle_parts = .false.
    if (ilattice==1) shuffle_parts = .true.
    call prompt('Relax particles by shuffling?',shuffle_parts)

    BEsphere = .false.
    call prompt('Centrally condense the sphere as a BE sphere?',BEsphere)

    if (.not. BEsphere) then
       r_sphere = 4.
       call prompt('Enter radius of sphere in units of '//dist_unit,r_sphere,0.)
       lbox     = 4.
       call prompt('Enter the box size in units of spherical radii: ',lbox,1.)
       if (.not. is_cube) then
          do i=1,3
             xmini(i) = -0.5*(lbox*r_sphere)
             xmaxi(i) = -xmini(i)
          enddo
       endif

       totmass_sphere = 1.0
       call prompt('Enter total mass in sphere in units of '//mass_unit,totmass_sphere,0.)
    else
       call prompt_BEparameters(iBEparam,BErho_cen,BErad_phys,BErad_norm,BEmass,BEfac,umass,udist,au,solarm)
       lbox     = 4.
       call prompt('Enter the box size in units of spherical radii: ',lbox,1.)
    endif

    density_contrast = 30.0
    call prompt('Enter density contrast between sphere and box ',density_contrast,1.)

    binary = .false.
    call prompt('Do you intend to form a binary system (i.e. add an m=2 perturbation)?',binary)

    if (binary) then
       cs_sphere_cgs = 18696.96 ! cm/s ~ 5K assuming mu = 2.31 & gamma = 5/3
    else
       cs_sphere_cgs = 21888.0  ! cm/s ~ 8K assuming mu = 2.31 & gamma = 5/3
    endif
    call prompt('Enter sound speed in sphere in units of cm/s',cs_sphere_cgs,0.)

    if (binary) then
       angvel = 1.006d-12
    else
       angvel = 1.77d-13
    endif
    angvel_not_betar = .true.
    beta_r = 0.02
    call prompt('Input angular velocity (true); else input ratio of rotational-to-potential energy ',angvel_not_betar)
    if (angvel_not_betar) then
       call prompt('Enter angular rotation speed in rad/s ',angvel,0.)
    else
       call prompt('Enter ratio of rotational-to-potential energy ',beta_r,0.)
    endif

    rms_mach = 0.
    call prompt('Enter the Mach number of the cloud turbulence',rms_mach,0.)

    if (mhd) then
       Bzero_G    = 1.0d-4 ! G
       masstoflux =   5.0
       ang_Bomega = 180.0
       mu_not_B   = .true.
       call prompt('Input the mass-to-flux ratio (true); else input the magnetic field strength ',mu_not_B)
       if (mu_not_B) then
          call prompt('Enter mass-to-flux ratio in units of critical value ',masstoflux,0.)
       else
          call prompt('Enter magnetic field strength in Gauss ',Bzero_G,0.)
       endif
       call prompt('Enter the angle (degrees) between B and the rotation axis? ',ang_Bomega)
    endif

    if (use_dust) then
       !--currently assume one fluid dust
       dtg          = 0.01
       grainsize    = 0.
       graindens    = 0.
       grainsizecgs = 0.1
       graindenscgs = 3.
       ndustsmall   = 1
       smincgs      = 1.e-5
       smaxcgs      = 1.
       sindex       = 3.5
       call prompt('Enter total dust to gas ratio',dtg,0.)
       call prompt('How many grain sizes do you want?',ndustsmall,1,maxdustsmall)
       ndusttypes = ndustsmall
       if (ndusttypes > 1) then
          !--grainsizes
          call prompt('Enter minimum grain size in cm',smincgs,0.)
          call prompt('Enter maximum grain size in cm',smaxcgs,0.)
          !--mass distribution
          call prompt('Enter power-law index, e.g. MRN',sindex)
          !--grain density
          call prompt('Enter grain density in g/cm^3',graindenscgs,0.)
       else
          call prompt('Enter grain size in cm',grainsizecgs,0.)
          call prompt('Enter grain density in g/cm^3',graindenscgs,0.)
       endif
    endif

    if (binary) then
       rho_pert_amp = 0.1
       call prompt('Enter the amplitute of the density perturbation ',rho_pert_amp,0.0,0.4)
    endif

    ! ask about sink particle details; these will not be saved to the .setup file since they exist in the .in file
    !
    call prompt('Do you wish to dynamically create sink particles? ',make_sinks)
    if (make_sinks) then
       if (binary) then
          h_acc_char  = '3.35au'
       else
          h_acc_char  = '1.0d-2'
       endif
       call prompt('Enter the accretion radius of the sink (with units; e.g. au,pc,kpc,0.1pc) ',h_acc_char)
       call select_unit(h_acc_char,h_acc_in,ierr)
       h_acc_setup = h_acc_in
       if (ierr==0 ) h_acc_setup = h_acc_setup/udist
       r_crit_setup        = 5.0*h_acc_setup
       icreate_sinks_setup = 1
       if (binary) h_soft_sinksink_setup = 0.4*h_acc_setup
    else
       icreate_sinks_setup = 0
       rhofinal_setup = 0.15
       call prompt('Enter final maximum density in g/cm^3 (ignored for <= 0) ',rhofinal_setup)
    endif
    if (id==master) call write_setupfile(filename)
    stop 'please edit .setup file and rerun phantomsetup'
 else
    stop ! MPI, stop on other threads, interactive on master
 endif
 !
 ! units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !
 ! set dust properties
 !
 if (use_dust) then
    use_dustfrac = .true.
    ndustsmall   = ndusttypes
    if (ndusttypes > 1) then
       call set_dustbinfrac(smincgs,smaxcgs,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
       grainsize(1:ndusttypes) = grainsize(1:ndusttypes)/udist
       graindens(1:ndusttypes) = graindenscgs/umass*udist**3
    else
       grainsize(1) = grainsizecgs/udist
       graindens(1) = graindenscgs/umass*udist**3
    endif
 endif
 !
 ! convert units of sound speed
 !
 if (cs_in_code) then
    cs_sphere_cgs = cs_sphere*unit_velocity
 else
    cs_sphere     = cs_sphere_cgs/unit_velocity
 endif
 !
 ! Bonnor-Ebert profile (if requested)
 !
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
 if (density_contrast < 1.+epsilon(density_contrast)) lbox = 2.
 do i = 1,3
    xmini(i) = -0.5*(lbox*r_sphere)
    xmaxi(i) = -xmini(i)
 enddo
 !
 ! boundaries
 !
 call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3))
 !
 ! general parameters
 !
 time        = 0.
 if (use_dust) dust_method = 1
 hfact       = hfact_default
 hfact_out   = hfact_default
 if (maxvxyzu >=4 ) then
    gamma    = 5./3.
 else
    gamma    = 1.
 endif
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
    do i = 1,3
       xmini(i) = -0.5*dxbound
       xmaxi(i) = -xmini(i)
    enddo
    call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3))
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
 rmasstoflux_crit = 2./3.*0.53*sqrt(5./pi) ! code units *see derivation at the end of the file*
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
 !
 ! setup particles in the sphere; use this routine to get N_sphere as close to np as possible
 !
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
    print*, psep_box, dxbound, vol_box,dens_sphere,massoftype(igas)
    call set_unifdis(trim(lattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep_box, &
                      hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong,err=ierr)
 endif
 npartsphere = npart
 if (np_in /= npartsphere) np = npartsphere
 massoftype(igas) = totmass_sphere/npartsphere  ! note: before 5 Oct 2021, this was based upon total mass & total number, not sphere numbers
 !
 ! setup surrounding low density medium
 !
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
 !
 ! set particle properties
 !
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 dustfrac          = 0.
 if (massoftype(igas) < epsilon(massoftype(igas))) massoftype(igas) = totmass/npart_total
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
    if (use_dust .and. dust_method==1) then
       if (ndusttypes > 1) then
          dustfrac(1:ndusttypes,i) = dustbinfrac(1:ndusttypes)*dtg
       else
          dustfrac(1,i) = dtg/(1.+dtg) ! call set_dustfrac(dtg,dustfrac(:,i))
       endif
    endif
 enddo
 !
 ! Set two-fluid dust
 ! (currently deactivated; will need to re-test before use to ensure it is fully compatible with the current dust algorithms)
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
 if (BEsphere) deallocate(rtab,rhotab)
 !
 ! reset to centre of mass
 ! (if random or shuffling, recentering may shift particles outside of the defined range)
 !
 if (trim(lattice)/='random' .and. .not.shuffle_parts) call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! temperature set to give a pressure equilibrium
 !
 polyk  = cs_sphere**2
 polyk2 = cs_medium**2
 !
 !--Stretching the spatial distribution to perturb the density profile, if requested
 !
 if (binary) then
    do i = 1,npart
       rxy2  = xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i)
       rxyz2 = rxy2 + xyzh(3,i)*xyzh(3,i)
       if (rxyz2 <= r_sphere**2) then
          phi       = atan(xyzh(2,i)/xyzh(1,i))
          if (xyzh(1,i) < 0.0) phi = phi + pi
          dphi      = 0.5*rho_pert_amp*sin(2.0*phi)
          phi       = phi - dphi
          xyzh(1,i) = sqrt(rxy2)*cos(phi)
          xyzh(2,i) = sqrt(rxy2)*sin(phi)
       endif
    enddo
 endif
 !
 ! reset values if turbulence-in-a-box
 !
 if (density_contrast < 1.+epsilon(density_contrast)) then
    r_sphere = xmax
    npartsphere = npart
 endif
 !
 ! Velocity: Turbulent velocity field
 !
 vxyzu = 0.
 if (rms_mach > 0.) then
    call getcwd(cwd)
    ! personal hack for J. Wurster since different computer clusters required different velocity fields
    if (index(cwd,'gpfs1/scratch/astro/jhw5') > 0 .or. index(cwd,'data/dp187/dc-wurs1') > 0 ) then
       ! Kennedy or Dial
       filex  = find_phantom_datafile(filevx,'velfield_sphng')
       filey  = find_phantom_datafile(filevy,'velfield_sphng')
       filez  = find_phantom_datafile(filevz,'velfield_sphng')
    else
       filex  = find_phantom_datafile(filevx,'velfield')
       filey  = find_phantom_datafile(filevy,'velfield')
       filez  = find_phantom_datafile(filevz,'velfield')
    endif

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
 endif
 !
 ! Velocity: uniform rotation (thermal energy & magnetic field too)
 !
 do i=1,npart
    r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
    if (r2 < r_sphere**2) then
       vxyzu(1,i) = vxyzu(1,i) - angvel_code*xyzh(2,i)
       vxyzu(2,i) = vxyzu(2,i) + angvel_code*xyzh(1,i)
       if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk
    else
       if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk2
    endif
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(1,i) = Bzero*sin(ang_Bomega*pi/180.0)
       Bxyz(3,i) = Bzero*cos(ang_Bomega*pi/180.0)
    endif
 enddo
 !
 ! set default runtime parameters if .in file does not exist
 !
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 dtmax = t_ff/100.  ! Since this variable can change, always reset it if running phantomsetup
 if (.not. iexist) then
    ! default values
    tmax          = 1.21*t_ff ! = 10.75 for default settings (Wurster, Price & Bate 2016)
    if (maxvxyzu < 4) ieos = 8
    nfulldump     = 1
    calc_erot     = .true.
    icreate_sinks = icreate_sinks_setup
    r_crit        = r_crit_setup
    h_acc         = h_acc_setup
    ! reset defaults based upon options
    if (density_contrast > 1.) dtmax_dratio = 1.258
    if (density_contrast < 1.+epsilon(density_contrast) .and. maxvxyzu>=4) then
       ieos     = 2
       icooling = 5
    endif
    if (binary) then
       h_soft_sinksink = h_soft_sinksink_setup
       tmax            = 1.50*t_ff ! = 13.33 for default settings (Wurster, Price & Bate 2017)
    endif
    if (icreate_sinks==1) then
       dtmax_min = dtmax/8.0
    else
       dtmax_min    = 0.0
       rhofinal_cgs = rhofinal_setup
    endif
    ilimitdustflux = .true.
 endif
 !
 !--Summarise the sphere
 !
 print "(a,i10)",' Input npart_sphere = ',np
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass,totmass*umass,' g'
 print fmt,' Mass in sphere   : ',totmass_sphere,totmass_sphere*umass,' g'
 print fmt,' Radius of sphere : ',r_sphere,r_sphere*udist,' cm'
 if (BEsphere) then
    print fmt,' Mean rho sphere  : ',dens_sphere,dens_sphere*unit_density,' g/cm^3'
    print fmt,' central density  : ',central_density,central_density*unit_density,' g/cm^3'
    print fmt,' edge density     : ',edge_density,edge_density*unit_density,' g/cm^3'
    print fmt,' Mean rho medium  : ',dens_medium,dens_medium*unit_density,' g/cm^3'
 else
    print fmt,' Density sphere   : ',dens_sphere,dens_sphere*unit_density,' g/cm^3'
    print fmt,' Density medium   : ',dens_medium,dens_medium*unit_density,' g/cm^3'
 endif
 print fmt,' cs in sphere     : ',cs_sphere,cs_sphere_cgs,' cm/s'
 print fmt,' cs in medium     : ',cs_medium,cs_medium*unit_velocity,' cm/s'
 print fmt,' Free fall time   : ',t_ff,t_ff*utime/years,' yrs'
 print fmt,' Angular velocity : ',angvel_code,angvel,' rad/s'
 print fmt,' Turbulent Mach no: ',rms_mach
 print fmt,' Omega*t_ff       : ',angvel_code*t_ff
 if (mhd) then
    print fmt,' B field (z)      : ',Bzero,Bzero*unit_Bfield*1.d6,' micro-G'
    print fmt,' Alfven speed     : ',Bzero/sqrt(dens_sphere),Bzero/sqrt(dens_sphere)*udist/utime,' cm/s'
    if (Bzero > 0.) then
       print fmt,' plasma beta      : ',przero/(0.5*Bzero*Bzero)
       print fmt,' mass-to-flux     : ',totmass_sphere/(area*Bzero)/rmasstoflux_crit
    endif
 endif
 if (use_dust) then
    print fmt,' dust-to-gas ratio: ',dtg,dtg,' '
 endif
 print "(1x,50('-'))"

end subroutine setpart

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20
 integer                      :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for sphere-in-box setup routines'
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)

 write(iunit,"(/,a)") '# particle resolution & placement'
 call write_inopt(np,'np','requested number of particles in sphere',iunit)
 call write_inopt(lattice,'lattice','particle lattice (random,cubic,closepacked,hcp,hexagonal)',iunit)
 call write_inopt(shuffle_parts,'shuffle_parts','relax particles by shuffling',iunit)

 write(iunit,"(/,a)") '# options for box'
 if (.not.BEsphere .and. .not.is_cube) then
    ! left here for backwards compatibility and for simplicity if the user requires a rectangle in the future
    do i=1,3
       call write_inopt(xmini(i),labelx(i)//'min',labelx(i)//' min',iunit)
       call write_inopt(xmaxi(i),labelx(i)//'max',labelx(i)//' max',iunit)
    enddo
 else
    call write_inopt(lbox,'lbox','length of a box side in terms of spherical radii',iunit)
 endif

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
 call write_inopt(cs_sphere_cgs,'cs_sphere_cgs','sound speed in sphere in cm/s',iunit)
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
 if (use_dust) then
    write(iunit,"(/,a)") '# Dust properties'
    call write_inopt(dtg,'dust_to_gas_ratio','dust-to-gas ratio',iunit)
    call write_inopt(ndusttypes,'ndusttypes','number of grain sizes',iunit)
    if (ndusttypes > 1) then
       call write_inopt(smincgs,'smincgs','minimum grain size [cm]',iunit)
       call write_inopt(smaxcgs,'smaxcgs','maximum grain size [cm]',iunit)
       call write_inopt(sindex, 'sindex', 'power-law index, e.g. MRN',iunit)
       call write_inopt(graindenscgs,'graindenscgs','grain density [g/cm^3]',iunit)
    else
       call write_inopt(grainsizecgs,'grainsizecgs','grain size in [cm]',iunit)
       call write_inopt(graindenscgs,'graindenscgs','grain density [g/cm^3]',iunit)
    endif
 endif
 if (binary) then
    call write_inopt(rho_pert_amp,'rho_pert_amp','amplitude of density perturbation',iunit)
 endif
 write(iunit,"(/,a)") '# Sink properties (values in .in file, if present, will take precedence)'
 call write_inopt(icreate_sinks_setup,'icreate_sinks','1: create sinks.  0: do not create sinks',iunit)
 if (icreate_sinks_setup==1) then
    call write_inopt(h_acc_setup,'h_acc','accretion radius (code units)',iunit)
    call write_inopt(r_crit_setup,'r_crit','critical radius (code units)',iunit)
    if (binary) then
       call write_inopt(h_soft_sinksink_setup,'h_soft_sinksink','sink-sink softening radius (code units)',iunit)
    endif
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
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use unifdis,      only: is_valid_lattice
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: i,nerr,jerr,kerr
 type(inopts), allocatable     :: db(:)

 !--Read values
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)
 call read_inopt(BEsphere,'use_BE_sphere',db,ierr)
 call read_inopt(binary,'form_binary',db,ierr)
 call read_inopt(np,'np',db,ierr)
 call read_inopt(lattice,'lattice',db,ierr)
 if (ierr/=0 .or. .not. is_valid_lattice(trim(lattice))) then
    print*, ' invalid lattice.  Setting to closepacked'
    lattice = 'closepacked'
 endif
 call read_inopt(shuffle_parts,'shuffle_parts',db,ierr)

 call read_inopt(lbox,'lbox',db,jerr)  ! for backwards compatibility
 if (jerr /= 0) then
    do i=1,3
       call read_inopt(xmini(i),labelx(i)//'min',db,ierr)
       call read_inopt(xmaxi(i),labelx(i)//'max',db,ierr)
    enddo
    lbox = -2.0*xmini(1)/r_sphere
 endif

 if (.not. BEsphere) then
    call read_inopt(r_sphere,'r_sphere',db,ierr)
    call read_inopt(totmass_sphere,'totmass_sphere',db,ierr)
 else
    call read_inopt(iBEparam,'iBE_options',db,ierr)
    if (iBEparam==1 .or. iBEparam==2 .or. iBEparam==3) call read_inopt(BErho_cen,'BErho_cen',db,ierr)
    if (iBEparam==1 .or. iBEparam==4 .or. iBEparam==6) call read_inopt(BErad_phys,'BErad_phys',db,ierr)
    if (iBEparam==2 .or. iBEparam==4 .or. iBEparam==5) call read_inopt(BErad_norm,'BErad_norm',db,ierr)
    if (iBEparam==3 .or. iBEparam==5 .or. iBEparam==6) call read_inopt(BEmass,'BEmass',db,ierr)
    if (iBEparam==4 .or. iBEparam==5)                  call read_inopt(BEfac,'BEfac',db,ierr)
 endif

 call read_inopt(density_contrast,'density_contrast',db,ierr)
 call read_inopt(cs_sphere,'cs_sphere',db,jerr)
 call read_inopt(cs_sphere_cgs,'cs_sphere_cgs',db,kerr)
 cs_in_code = .false.  ! for backwards compatibility
 if (jerr /= 0 .and. kerr == 0) then
    cs_in_code = .false.
 elseif (jerr == 0 .and. kerr /= 0) then
    cs_in_code = .true.
 else
    ierr = ierr + 1
 endif
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
       ierr = ierr + 1
    endif
 endif
 if (use_dust) then
    call read_inopt(dtg,'dust_to_gas_ratio',db,ierr)
    call read_inopt(ndusttypes,'ndusttypes',db,ierr)
    if (ndusttypes > 1) then
       call read_inopt(smincgs,'smincgs',db,ierr)
       call read_inopt(smaxcgs,'smaxcgs',db,ierr)
       call read_inopt(sindex,'cs_sphere',db,ierr)
       call read_inopt(graindenscgs,'graindenscgs',db,ierr)
    else
       call read_inopt(grainsizecgs,'grainsizecgs',db,ierr)
       call read_inopt(graindenscgs,'graindenscgs',db,ierr)
    endif
 endif
 if (binary) then
    call read_inopt(rho_pert_amp,'rho_pert_amp',db,ierr)
 endif
 call read_inopt(icreate_sinks_setup,'icreate_sinks',db,ierr)
 if (icreate_sinks_setup==1) then
    call read_inopt(h_acc_setup,'h_acc',db,ierr)
    call read_inopt(r_crit_setup,'r_crit',db,ierr)
    if (binary) then
       call read_inopt(h_soft_sinksink_setup,'h_soft_sinksink',db,ierr)
    endif
 else
    call read_inopt(rhofinal_setup,'rho_final',db,ierr)
 endif
 call close_db(db)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_sphereinbox','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_sphereinbox','length unit not recognised')
    ierr = ierr + 1
 endif

 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_sphereinbox: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile
!----------------------------------------------------------------
 !--Magnetic flux justification
 !  This shows how the critical mass-to-flux values translates from CGS to code units.
 !
 ! rmasstoflux_crit = 0.53/(3*pi)*sqrt(5./G)                                ! cgs units of g G^-1 cm^-2
 ! convert base units from cgs to code:
 ! rmasstoflux_crit = 0.53/(3*pi)*sqrt(5./G)    *unit_Bfield*udist**2/umass
 ! where
 ! unit_Bfield   = umass/(utime*sqrt(umass*udist/4*pi)) = sqrt(4.*pi*umass)/(utime*sqrt(udist))
 ! therefore
 ! rmasstoflux_crit = 0.53/(3*pi)*sqrt(5./G)    *sqrt(4.*pi*umass)*udist**2/(utime*sqrt(udist)*umass)
 ! rmasstoflux_crit = (2/3)*0.53*sqrt(5./(G*pi))*sqrt(umass)*udist**2/(utime*sqrt(udist)*umass)
 ! rmasstoflux_crit = (2/3)*0.53*sqrt(5./(G*pi))*udist**1.5/ (sqrt(umass)*utime)
 ! where
 ! G [cgs] = 1 * udist**3/(umass*utime**2)
 ! therefore
 ! rmasstoflux_crit = (2/3)*0.53*sqrt(5./pi)    *udist**1.5/ (sqrt(umass)*utime) / sqrt(udist**3/(umass*utime**2))
 ! rmasstoflux_crit = (2/3)*0.53*sqrt(5./pi)                                ! code units

!----------------------------------------------------------------
end module setup
