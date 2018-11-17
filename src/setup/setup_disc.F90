!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   This module sets up accretion discs. The central object(s) can be
!   modelled with sink particles or external potentials. Systems with two
!   sink particles:
!     (i)  in a bound binary can have circumbinary, circumprimary, and
!          circumsecondary discs,
!     (ii) in an unbound binary (i.e. a fly-by) can have circumprimary and
!          circumsecondary discs.
!   In addition to gas, each disc can contain dust, modelled with either the
!   one fluid or two fluid methods. The dust can only grow in the two-fluid method.
!   Embedded planets can be added to single or circumbinary discs.
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    Ratm_in       -- inner atmosphere radius (planet radii)
!    Ratm_out      -- outer atmosphere radius (planet radii)
!    accr1         -- central star accretion radius
!    accr2         -- perturber accretion radius
!    alphaSS       -- desired alphaSS
!    atm_type      -- atmosphere type (1:r**(-3); 2:r**(-1./(gamma-1.)))
!    bhspin        -- black hole spin
!    bhspinangle   -- black hole spin angle (deg)
!    binary_O      -- Omega, PA of ascending node (deg)
!    binary_a      -- binary semi-major axis
!    binary_e      -- binary eccentricity
!    binary_f      -- f, initial true anomaly (deg,180=apastron)
!    binary_i      -- i, inclination (deg)
!    binary_w      -- w, argument of periapsis (deg)
!    deltat        -- output interval as fraction of orbital period
!    dist_unit     -- distance unit (e.g. au,pc,kpc,0.1pc)
!    einst_prec    -- include Einstein precession
!    flyby_O       -- position angle of ascending node (deg)
!    flyby_a       -- distance of minimum approach
!    flyby_d       -- initial distance (units of dist. min. approach)
!    flyby_i       -- inclination (deg)
!    ibinary       -- binary orbit (0=bound,1=unbound [flyby])
!    ipotential    -- potential (1=central point mass,
!    m1            -- central star mass
!    m2            -- perturber mass
!    mass_unit     -- mass unit (e.g. solarm,jupiterm,earthm)
!    norbits       -- maximum number of orbits at outer disc
!    np            -- number of gas particles
!    nplanets      -- number of planets
!    nsinks        -- number of sinks
!    ramp          -- Do you want to ramp up the planet mass slowly?
!    rho_core      -- planet core density (cgs units)
!    setplanets    -- add planets? (0=no,1=yes)
!    surface_force -- model m1 as planet with surface
!    use_mcfost    -- use the mcfost library
!
!  DEPENDENCIES: centreofmass, dim, dust, eos, extern_binary,
!    extern_corotate, extern_lensethirring, externalforces, fileutils,
!    growth, infile_utils, io, kernel, options, part, physcon, prompting,
!    set_dust, set_dust_options, setbinary, setdisc, setflyby, spherical,
!    table_utils, timestep, units, vectorutils
!+
!--------------------------------------------------------------------------
module setup
 use dim,              only:maxp,use_dust,maxalpha,use_dustgrowth,maxdusttypes,&
                            maxdustlarge,maxdustsmall
 use externalforces,   only:iext_star,iext_binary,iext_lensethirring,&
                            iext_einsteinprec,iext_corot_binary,iext_corotate
 use options,          only:use_dustfrac,iexternalforce
#ifdef MCFOST
 use options,          only:use_mcfost,nfulldump
#endif
 use part,             only:ndusttypes,ndustsmall,ndustlarge,grainsize,graindens
 use physcon,          only:au,solarm
 use setdisc,          only:scaled_sigma
 use extern_binary,    only:ramp,surface_force
 use set_dust_options, only:dust_method,dust_to_gas,ndusttypesinp,&
                            grainsizeinp,graindensinp,igrainsize,igraindens,&
                            iprofile_dust,smincgs,smaxcgs,sindex,dustbinfrac,&
                            Kdrag,ilimitdustfluxinp
 use dust,             only:ilimitdustflux,grainsizecgs,graindenscgs
 use set_dust,         only:set_dustbinfrac
 use fileutils,        only:make_tags_unique

 implicit none
 public  :: setpart

 integer :: np,np_dust(maxdustlarge),norbits,i

 !--central objects
 real    :: m1,m2,accr1,accr2,bhspin,bhspinangle,flyby_a,flyby_d,flyby_O,flyby_i
 real    :: binary_a,binary_e,binary_i,binary_O,binary_w,binary_f,deltat
 integer :: icentral,ipotential,nsinks,ibinary
 logical :: einst_prec

 !--discs
 character(len=20) :: disclabel
 character(len=*), dimension(3), parameter :: disctype = &
    (/'binary   ', &
      'primary  ', &
      'secondary'/)
 logical :: iuse_disc(3),itapergas(3),itaperdust(3),iwarp(3)
 logical :: ismoothgas(3),ismoothdust(3),use_global_iso
 integer :: mass_set(3),ndiscs
 real    :: R_in(3),R_out(3),R_ref(3),R_c(3),R_warp(3),H_warp(3)
 real    :: pindex(3),qindex(3),H_R(3),posangl(3),incl(3)
 real    :: disc_m(3),sig_ref(3),sig_norm(3),annulus_m(3),R_inann(3),R_outann(3),Q_min(3)
 real    :: R_indust(3),R_indust_swap(3),R_outdust(3),R_outdust_swap(3),R_c_dust(3)
 real    :: pindex_dust(3),qindex_dust(3),H_R_dust(3)
 real    :: ldisc(3),lcentral(3)
 real    :: alphaSS,dist_bt_sinks

 !--planets
 integer, parameter :: maxplanets = 9
 integer :: nplanets,setplanets
 real    :: mplanet(maxplanets),rplanet(maxplanets),accrplanet(maxplanets),inclplan(maxplanets)
 character(len=*), dimension(maxplanets), parameter :: planets = &
    (/'1','2','3','4','5','6','7','8','9' /)

 !--planet atmosphere defaults
 integer :: atm_type     = 1
 real    :: rho_core_cgs = 5.
 real    :: Ratm_in      = 1.
 real    :: Ratm_out     = 3.
 real    :: Natmfrac     = 0.

 !--units
 character(len=20) :: dist_unit,mass_unit

contains

!----------------------------------------------------------------
!
! This subroutine sets up a one or more accretion disc(s)
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use centreofmass,         only:reset_centreofmass
 use eos,                  only:isink,qfacdisc
 use extern_binary,        only:accradius1,accradius2,binarymassr,eps_soft1
 use externalforces,       only:mass1,accradius1
 use extern_lensethirring, only:blackhole_spin,blackhole_spin_angle
 use io,                   only:master,warning,error,fatal
 use kernel,               only:hfact_default
 use options,              only:ieos,alpha,icooling
 use part,                 only:nptmass,xyzmh_ptmass,maxvxyzu,vxyz_ptmass,ihacc,&
                                ihsoft,igas,idust,iamtype,iphase,dustprop,dustfrac
 use physcon,              only:jupiterm,earthm,pi,years
 use setbinary,            only:set_binary,Rochelobe_estimate,get_mean_angmom_vector
 use setdisc,              only:set_disc,get_disc_mass
 use set_dust_options,     only:set_dust_default_options,check_dust_method
 use setflyby,             only:set_flyby,get_T_flyby
 use table_utils,          only:logspace
 use timestep,             only:tmax,dtmax
 use units,                only:set_units,select_unit,umass,udist,utime
 use vectorutils,          only:rotatevec
#ifdef MCFOST
 use options,              only:alphau,ipdv_heating,ishock_heating
#endif

 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 integer, parameter :: maxbins = 4096
 logical :: iexist,seq_exists,is_isothermal,ichange_method
 real    :: phi,vphi,sinphi,cosphi,omega,r2,disc_m_within_r,period_longest
 real    :: jdust_to_gas,Rj,period,Rochelobe,Hill(maxplanets)
 real    :: totmass_gas,totmass_dust,mcentral,R,Sigma,Sigmadust,Stokes(maxdusttypes)
 real    :: polyk_dust,xorigini(3),vorigini(3),alpha_returned(3)
 real    :: star_m(3),disc_mdust(3,maxdusttypes),sig_normdust(3),u(3)
 real    :: enc_m(maxbins),rad(maxbins),Q_mintmp,disc_mtmp(3),annulus_mtmp(3)
 integer :: maxdiscs
 integer :: ierr,j,idisc,nparttot,npingasdisc,npindustdisc,itype
 integer :: sigmaprofilegas(3),sigmaprofiledust(3),iprofilegas(3),iprofiledust(3)
 character(len=100) :: filename
 character(len=100) :: prefix
 character(len=20)  :: varstring(maxdusttypes)

 integer :: npart_planet_atm,npart_recentre
 integer :: npart_disc
 real, parameter :: a0 = 1.
 real    :: r_surface
 real    :: udens,rho_core

!--------------------------------------------------------------------------
!+
!  Initialization
!+
!--------------------------------------------------------------------------

 print "(/,65('-'),2(/,a),/,65('-'),/)"
 print "(a)",'     Welcome to the New Disc Setup'
 print "(/,65('-'),2(/,a),/,65('-'),/)"

 !
 !--set defaults
 !
 call set_dust_default_options()

 !
 !--get disc setup parameters from file or interactive setup
 !
 filename=trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then

    !--read from setup file
    call read_setupfile(filename,ierr)
    if (id==master) call write_setupfile(filename)
    if (ierr /= 0) then
       stop
    endif

    !--setup multiple discs each with different orbital parameters
    !  specified in orbits.dat
    inquire(file='orbits.dat',exist=iexist)
    inquire(file=trim(fileprefix)//'A.setup',exist=seq_exists)
    if (iexist .and. .not.seq_exists) then
       open(unit=23,file='orbits.dat',status='old',iostat=ierr)
       j = 0
       do while(ierr==0)
          read(23,*,iostat=ierr) binary_a,binary_e,binary_i,binary_O,binary_w,binary_f
          if (ierr==0) then
             j = j + 1
             write(filename,"(a)") trim(fileprefix)//achar(j+64)//'.setup'
             if (id==master) call write_setupfile(filename)
          endif
       enddo
       close(unit=23)
       stop
    endif

 elseif (id==master) then

    !--interactive setup
    print "(a,/)",' '//trim(filename)//' not found: using interactive setup'
    call setup_interactive()

    !--write setup file from interactive setup
    call write_setupfile(filename)
    print "(/,a)",' >>> please edit '//trim(filename)//' to set parameters for your problem then rerun phantomsetup <<<'
    stop

 else

    stop

 endif

!--------------------------------------------------------------------------
!+
!  Perform setup
!+
!--------------------------------------------------------------------------

 !
 !--set units
 !
 if (icentral==0 .and. ipotential==3) then
    !--black hole units
    !  note: distance unit not used but (currently) required in the .setup file
    call select_unit(mass_unit,umass,ierr)
    if (ierr /= 0) call error('setup_disc','mass unit not recognised')
    call set_units(mass=umass,c=1.d0)
 else
    !--stellar units
    call select_unit(mass_unit,umass,ierr)
    if (ierr /= 0) call error('setup_disc','mass unit not recognised')
    call select_unit(dist_unit,udist,ierr)
    if (ierr /= 0) call error('setup_disc','length unit not recognised')
    call set_units(dist=udist,mass=umass,G=1.d0)
 endif

 !
 !--how many discs?
 !
 maxdiscs = 1
 if ((icentral==1) .and. (nsinks==2)) maxdiscs = 3
 ndiscs = max(count(iuse_disc),1)
 !--index of disc (if only one)
 idisc = 0
 if (ndiscs==1) then
    do i=1,3
       if (iuse_disc(i)) idisc = i
    enddo
 endif

 is_isothermal = (maxvxyzu==3)
#ifdef MCFOST
 if (use_mcfost) then
    is_isothermal = .false.
    nfulldump = 1
 else ! We are in the isothermal case
    is_isothermal = .true.
 endif
#endif

 !
 !--equation of state
 !
 if (is_isothermal) then
    !--isothermal
    gamma = 1.0
    if (ndiscs /= 1) then
       !--multiple discs
       if (use_global_iso) then
          !--globally isothermal
          ieos = 1
          qindex = 0.
          qfacdisc = qindex(1)
          print "(/,a)",' setting ieos=1 for globally isothermal disc'
          if (iuse_disc(1)) then
             H_R(2) = sqrt(R_ref(2)/R_ref(1)*(m1+m2)/m1) * H_R(1)
             H_R(3) = sqrt(R_ref(3)/R_ref(1)*(m1+m2)/m2) * H_R(1)
             call warning('setup_disc','using circumbinary (H/R)_ref to set global temperature')
          elseif (iuse_disc(2)) then
             H_R(3) = sqrt(R_ref(3)/R_ref(2)*m1/m2) * H_R(2)
             call warning('setup_disc','using circumprimary (H/R)_ref to set global temperature')
          endif
       else
          !--locally isothermal prescription from Farris et al. (2014) for binary system
          ieos = 14
          print "(/,a)",' setting ieos=14 for locally isothermal from Farris et al. (2014)'
          if(iuse_disc(1)) then
             qfacdisc = qindex(1)
             call warning('setup_disc','using circumbinary (H/R)_ref to set global temperature')
          elseif(iuse_disc(2))then
             qfacdisc = qindex(2)
             call warning('setup_disc','using circumprimary (H/R)_ref to set global temperature')
          endif
       endif
    else
       !--single disc
       if (qindex(idisc) > 0.) then
          do i=1,3
             !--eos around sink
             if (iuse_disc(i)) isink = i-1
          enddo
          !--locally isothermal
          if (isink /= 0) then
             ieos = 6
             print "(/,a)",' setting ieos=6 for locally isothermal disc around sink'
          else
             ieos = 3
             print "(/,a)",' setting ieos=3 for locally isothermal disc around origin'
          endif
          qfacdisc = qindex(idisc)
       endif
    endif
 else
    !--adiabatic
    ieos = 2
    gamma = 5./3.
    icooling = 1

#ifdef MCFOST
    if (use_mcfost) then
       icooling = 0
       ipdv_heating = 0
       ishock_heating = 0
       alphau = 0
    endif
#endif
 endif

 !
 !--surface density profile
 !
 iprofilegas = 0
 sigmaprofilegas = 0
 do i=1,3
    if (itapergas(i)) then
       iprofilegas(i) = 1
       sigmaprofilegas(i) = 1
    endif
    if (ismoothgas(i)) sigmaprofilegas(i) = 2
    if (itapergas(i) .and. ismoothgas(i)) sigmaprofilegas(i) = 3
 enddo
 if (use_dust) then
    iprofiledust = 0
    sigmaprofiledust = 0
    do i=1,3
       if (itaperdust(i)) then
          iprofiledust(i) = 1
          sigmaprofiledust(i) = 1
       endif
       if (ismoothdust(i)) sigmaprofiledust(i) = 2
       if (itaperdust(i) .and. ismoothdust(i)) sigmaprofiledust(i) = 3
    enddo
 endif

 !
 !--set sink particle(s) or potential
 !
 mcentral = m1
 select case (icentral)
 case (0)
    select case (ipotential)
    case (1)
       print "(/,a)",' Central point mass represented by external force with accretion boundary'
       print "(a,g10.3,a)",'   Object mass:      ', m1,    trim(mass_unit)
       print "(a,g10.3,a)",'   Accretion Radius: ', accr1, trim(dist_unit)
       mass1      = m1
       accradius1 = accr1
       mcentral   = m1
    case (2)
       print "(/,a)",' Central binary represented by external force with accretion boundary'
       print "(a,g10.3,a)",'   Primary mass:       ', m2,    trim(mass_unit)
       print "(a,g10.3)",  '   Binary mass ratio:  ', m1/m2
       print "(a,g10.3,a)",'   Accretion Radius 1: ', accr1, trim(dist_unit)
       print "(a,g10.3,a)",'   Accretion Radius 2: ', accr2, trim(dist_unit)
       mass1       = m1
       binarymassr = m1/m2
       accradius1  = accr1
       accradius2  = accr2
       if (iexternalforce == iext_corot_binary) then
          mcentral = m2
       else
          mcentral = m1 + m2
       endif
    case (3)
       print "(/,a)",' Central black hole represented by external force with accretion boundary'
       print "(a,g10.3,a)",'   Black hole mass:        ', m1,    trim(mass_unit)
       print "(a,g10.3,a)",'   Accretion Radius:       ', accr1, trim(dist_unit)
       print "(a,g10.3)",  '   Black hole spin:        ', bhspin
       print "(a,g10.3,a)",'   Black hole spin angle:  ', bhspinangle, 'deg'
       mass1                = m1
       accradius1           = accr1
       blackhole_spin       = bhspin
       blackhole_spin_angle = bhspinangle*(pi/180.0)
       mcentral             = m1
    end select
 case (1)
    select case (nsinks)
    case (1)
       !--single star
       print "(/,a)",' Central object represented by a sink at the system origin'
       print "(a,g10.3,a)",'   Object mass:      ', m1,    trim(mass_unit)
       print "(a,g10.3,a)",'   Accretion Radius: ', accr1, trim(dist_unit)
       nptmass                      = 1
       xyzmh_ptmass(:,:)            = 0.
       xyzmh_ptmass(1:3,nptmass)    = 0.
       xyzmh_ptmass(4,nptmass)      = m1
       xyzmh_ptmass(ihacc,nptmass)  = accr1
       xyzmh_ptmass(ihsoft,nptmass) = accr1
       vxyz_ptmass                  = 0.
       mcentral                     = m1
    case (2)
       !--binary
       select case (ibinary)
       case (0)
          !--bound
          nptmass  = 0
          print "(/,a)",' Central objects represented by two sinks'
          print "(a,g10.3,a)",'   Primary mass:       ', m1,    trim(mass_unit)
          print "(a,g10.3)",  '   Binary mass ratio:  ', m2/m1
          print "(a,g10.3,a)",'   Accretion Radius 1: ', accr1, trim(dist_unit)
          print "(a,g10.3,a)",'   Accretion Radius 2: ', accr2, trim(dist_unit)
          call set_binary(m1,massratio=m2/m1,semimajoraxis=binary_a,eccentricity=binary_e, &
                          posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i, &
                          f=binary_f,accretion_radius1=accr1,accretion_radius2=accr2, &
                          xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass)
          mcentral = m1 + m2
       case (1)
          !--unbound (flyby)
          print "(/,a)",' Central object represented by a sink at the system origin with a perturber sink'
          print "(a,g10.3,a)",'   Primary mass:       ', m1,    trim(mass_unit)
          print "(a,g10.3,a)",'   Perturber mass:     ', m2,    trim(mass_unit)
          print "(a,g10.3,a)",'   Accretion Radius 1: ', accr1, trim(dist_unit)
          print "(a,g10.3,a)",'   Accretion Radius 2: ', accr2, trim(dist_unit)
          call set_flyby(mprimary=m1,massratio=m2/m1,minimum_approach=flyby_a, &
                         initial_dist=flyby_d,posang_ascnode=flyby_O,inclination=flyby_i, &
                         accretion_radius1=accr1,accretion_radius2=accr2, &
                         xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass)
          mcentral = m1
       end select
    end select
 end select
 !--set array of central object masses
 star_m = (/mcentral, m1, m2/)
 do i=1,3
    if (.not.iuse_disc(i)) star_m(i) = 0.
 enddo

 !
 !--dust
 !
 if (use_dust) then
    grainsize = 0.
    graindens = 0.
    if (ndusttypes > 1) then
       select case(igrainsize)
       case(0)
          call set_dustbinfrac(smincgs,smaxcgs,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
          grainsize(1:ndusttypes) = grainsize(1:ndusttypes)/udist
       case(1)
          grainsize(1:ndusttypes) = grainsizeinp(1:ndusttypes)/udist
       end select
       select case(igraindens)
       case(0)
          graindens(1:ndusttypes) = graindensinp(1)/umass*udist**3
       case(1)
          graindens(1:ndusttypes) = graindensinp(1:ndusttypes)/umass*udist**3
       end select
    else
       grainsize(1) = grainsizeinp(1)/udist
       graindens(1) = graindensinp(1)/umass*udist**3
       grainsizecgs = grainsizeinp(1)
       graindenscgs = graindensinp(1)
    endif
 endif

 !
 !--compute the disc mass for different mass_set values
 !  and calculate sigma normalisation
 !
 totmass_gas  = 0.
 totmass_dust = 0.
 do i=1,3
    if (iuse_disc(i)) then
       select case(mass_set(i))
       case (0)
          !--set sigma normalisation from disc mass
          sig_norm(i) = 1.d0
          call get_disc_mass(disc_mtmp(i),enc_m,rad,Q_mintmp,sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
          sig_norm(i) = sig_norm(i) * disc_m(i) / disc_mtmp(i)
       case (1)
          !--set disc mass from annulus mass
          sig_norm(i) = 1.d0
          call get_disc_mass(annulus_mtmp(i),enc_m,rad,Q_mintmp,sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_inann(i),R_outann(i),R_ref(i),R_c(i), &
                             H_R(i))
          sig_norm(i) = sig_norm(i) * annulus_m(i) / annulus_mtmp(i)
          call get_disc_mass(disc_m(i),enc_m,rad,Q_min(i),sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
       case (2)
          !--set disc mass from sigma normalisation
          call get_disc_mass(disc_m(i),enc_m,rad,Q_min(i),sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
       case (3)
          !--set disc mass from sigma at reference radius
          if (.not.(R_in(i) < R_ref(i)) .and. ismoothgas(i)) call fatal('setup_disc', &
             'if smoothing inner disc and setting disc mass by sigma(R_ref), require R_in < R_ref')
          sig_norm(i) = sig_ref(i) / scaled_sigma(R_ref(i),sigmaprofilegas(i),pindex(i),R_ref(i),R_in(i),R_c(i))
          call get_disc_mass(disc_m(i),enc_m,rad,Q_min(i),sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
       case (4)
          !--set disc mass from minimum Toomre Q
          sig_norm(i) = 1.d0
          call get_disc_mass(disc_mtmp(i),enc_m,rad,Q_mintmp,sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
          sig_norm(i) = sig_norm(i) * Q_mintmp / Q_min(i)
          !--recompute actual disc mass and Toomre Q
          call get_disc_mass(disc_m(i),enc_m,rad,Q_min(i),sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
       end select
       totmass_gas = totmass_gas + disc_m(i)
       if (use_dust) then
          disc_mdust(i,:) = 0.
          do j=1,ndusttypes
             disc_mdust(i,j) = disc_m(i) * dust_to_gas * dustbinfrac(j)
             totmass_dust    = totmass_dust + disc_mdust(i,j)
          enddo
       endif
    endif
 enddo

 !
 !--setup disc(s)
 !
 time  = 0.
 hfact = hfact_default
 incl    = incl*(pi/180.0)
 posangl = posangl*(pi/180.0)
 if (maxalpha==0) alpha = alphaSS
 nparttot = 0
 npartoftype(:) = 0
 do i=1,3
    if (iuse_disc(i)) then
       !--set disc origin
       if (ndiscs > 1) then
          print "(/,a)",'>>> Setting up circum'//trim(disctype(i))//' disc <<<'
          prefix = trim(fileprefix)//'-'//disctype(i)
       else
          prefix = fileprefix
       endif
       select case(i)
       case(3)
          !--circumsecondary
          xorigini  = xyzmh_ptmass(1:3,2)
          vorigini  = vxyz_ptmass(1:3,2)
          Rochelobe = Rochelobe_estimate(m1,m2,binary_a)
       case(2)
          !--circumprimary
          xorigini  = xyzmh_ptmass(1:3,1)
          vorigini  = vxyz_ptmass(1:3,1)
          Rochelobe = Rochelobe_estimate(m2,m1,binary_a)
       case default
          !--single disc or circumbinary
          !  centre of mass of binary defined to be zero (see set_binary)
          xorigini  = 0.
          vorigini  = 0.
          Rochelobe = huge(0.)
       end select
       if (maxdiscs > 1 .and. ibinary==0) then
          if (R_out(i) > Rochelobe) call warning('setup_disc', &
             'Outer disc radius for circum'//trim(disctype(i))//' > Roche lobe of ' &
             //trim(disctype(i)))
       endif

       !--set disc(s)
       if (use_dust .and. use_dustfrac) then
          !--gas and dust mixture disc
          npingasdisc = int(disc_m(i)/totmass_gas*np)
          call set_disc(id,master        = master,               &
                        mixture          = .true.,               &
                        npart            = npingasdisc,          &
                        npart_start      = nparttot + 1,         &
                        rref             = R_ref(i),             &
                        rmin             = R_in(i),              &
                        rmax             = R_out(i),             &
                        rmindust         = R_indust(i),          &
                        rmaxdust         = R_outdust(i),         &
                        indexprofile     = iprofilegas(i),       &
                        indexprofiledust = iprofiledust(i),      &
                        rc               = R_c(i),               &
                        rcdust           = R_c_dust(i),          &
                        p_index          = pindex(i),            &
                        p_indexdust      = pindex_dust(i),       &
                        q_index          = qindex(i),            &
                        q_indexdust      = qindex_dust(i),       &
                        HoverR           = H_R(i),               &
                        HoverRdust       = H_R_dust(i),          &
                        disc_mass        = disc_m(i),            &
                        disc_massdust    = sum(disc_mdust(i,:)), &
                        star_mass        = star_m(i),            &
                        gamma            = gamma,                &
                        particle_mass    = massoftype(igas),     &
                        xyz_origin       = xorigini,             &
                        vxyz_origin      = vorigini,             &
                        hfact            = hfact,                &
                        xyzh             = xyzh,                 &
                        vxyzu            = vxyzu,                &
                        polyk            = polyk,                &
                        alpha            = alpha,                &
                        ismooth          = ismoothgas(i),        &
                        position_angle   = posangl(i),           &
                        inclination      = incl(i),              &
                        rwarp            = R_warp(i),            &
                        warp_smoothl     = H_warp(i),            &
                        bh_spin          = bhspin,               &
                        prefix           = prefix)
          !--swap radii to keep dust profile the same as gas within [R_indust,R_outdust]
          if (iprofile_dust == 2) then
             R_indust_swap(i)  = R_in(i)
             R_outdust_swap(i) = R_out(i)
          else
             R_indust_swap(i)  = R_indust(i)
             R_outdust_swap(i) = R_outdust(i)
          endif

          !--set dustfrac
          sig_normdust(i) = 1.d0
          call get_disc_mass(disc_mtmp(i),enc_m,rad,Q_mintmp,sigmaprofiledust(i), &
                             sig_normdust(i),star_m(i),pindex_dust(i),qindex_dust(i), &
                             R_indust_swap(i),R_outdust_swap(i),R_ref(i),R_c_dust(i),H_R_dust(i))
          sig_normdust(i) = sig_normdust(i) * sum(disc_mdust(i,:)) / disc_mtmp(i)
          do j=nparttot+1,npingasdisc
             Rj = sqrt(dot_product(xyzh(1:2,j)-xorigini(1:2),xyzh(1:2,j)-xorigini(1:2)))
             if (iprofile_dust > 0 .and. (Rj<R_indust(i) .or. Rj>R_outdust(i))) then
                jdust_to_gas = tiny(jdust_to_gas)
             else
                call get_dust_to_gas(jdust_to_gas,Rj,sigmaprofilegas(i), &
                                     sigmaprofiledust(i),sig_norm(i),sig_normdust(i), &
                                     pindex(i),pindex_dust(i),R_in(i),R_indust_swap(i),R_ref(i),&
                                     xyzh(3,j)-xorigini(3),get_H(H_R(i)*R_ref(i),qindex(i),Rj/R_ref(i)),&
                                     get_H(H_R_dust(i)*R_ref(i),qindex_dust(i),Rj/R_ref(i)),R_c(i),R_c_dust(i))
             endif
             jdust_to_gas = max(jdust_to_gas,tiny(jdust_to_gas))
             dustfrac(:,j) = (jdust_to_gas/(1.+jdust_to_gas))*dustbinfrac(:)
          enddo
          nparttot = nparttot + npingasdisc
       else
          !--gas disc
          npingasdisc = int(disc_m(i)/totmass_gas*np)
          call set_disc(id,master       = master,             &
                        npart           = npingasdisc,        &
                        npart_start     = nparttot + 1,       &
                        particle_type   = igas,               &
                        rref            = R_ref(i),           &
                        rmin            = R_in(i),            &
                        rmax            = R_out(i),           &
                        indexprofile    = iprofilegas(i),     &
                        rc              = R_c(i),             &
                        p_index         = pindex(i),          &
                        q_index         = qindex(i),          &
                        HoverR          = H_R(i),             &
                        disc_mass       = disc_m(i),          &
                        star_mass       = star_m(i),          &
                        gamma           = gamma,              &
                        particle_mass   = massoftype(igas),   &
                        xyz_origin      = xorigini,           &
                        vxyz_origin     = vorigini,           &
                        hfact           = hfact,              &
                        xyzh            = xyzh,               &
                        vxyzu           = vxyzu,              &
                        polyk           = polyk,              &
                        alpha           = alpha,              &
                        ismooth         = ismoothgas(i),      &
                        position_angle  = posangl(i),         &
                        inclination     = incl(i),            &
                        rwarp           = R_warp(i),          &
                        warp_smoothl    = H_warp(i),          &
                        bh_spin         = bhspin,             &
                        prefix          = prefix)
          nparttot = nparttot + npingasdisc
          npartoftype(igas) = npartoftype(igas) + npingasdisc
          if (use_dust) then
             !--dust disc
             do j=1,ndustlarge
                npindustdisc = int(disc_mdust(i,j)/sum(disc_mdust(:,j))*np_dust(j))
                itype = idust + j - 1
                call set_disc(id,master      = master,             &
                              npart          = npindustdisc,       &
                              npart_start    = nparttot + 1,       &
                              particle_type  = itype,              &
                              rref           = R_ref(i),           &
                              rmin           = R_indust(i),        &
                              rmax           = R_outdust(i),       &
                              indexprofile   = iprofiledust(i),    &
                              rc             = R_c_dust(i),        &
                              p_index        = pindex_dust(i),     &
                              q_index        = qindex_dust(i),     &
                              HoverR         = H_R_dust(i),        &
                              disc_mass      = disc_mdust(i,j),    &
                              star_mass      = star_m(i),          &
                              gamma          = gamma,              &
                              particle_mass  = massoftype(itype),  &
                              xyz_origin     = xorigini,           &
                              vxyz_origin    = vorigini,           &
                              hfact          = hfact,              &
                              xyzh           = xyzh,               &
                              vxyzu          = vxyzu,              &
                              polyk          = polyk_dust,         &
                              ismooth        = ismoothdust(i),     &
                              position_angle = posangl(i),         &
                              inclination    = incl(i),            &
                              rwarp          = R_warp(i),          &
                              warp_smoothl   = H_warp(i),          &
                              bh_spin        = bhspin,             &
                              prefix         = prefix)
                npartoftype(itype) = npartoftype(itype) + npindustdisc
                nparttot  = nparttot  + npindustdisc
                sig_normdust(i) = 1.d0
                call get_disc_mass(disc_mtmp(i),enc_m,rad,Q_mintmp,sigmaprofiledust(i), &
                                   sig_normdust(i),star_m(i),pindex_dust(i),qindex_dust(i), &
                                   R_indust(i),R_outdust(i),R_ref(i),R_c_dust(i),H_R_dust(i))
                sig_normdust(i) = sig_normdust(i) * disc_mdust(i,j) / disc_mtmp(i)
             enddo
          endif
       endif
       !--reset alpha for each disc
       alpha_returned(i) = alpha
    endif
 enddo

 if (surface_force) then
    npart_planet_atm = floor(Natmfrac*np)
    npart_disc = nparttot - npart_planet_atm

    udens = umass/udist**3
    rho_core  = rho_core_cgs/udens
    r_surface = (3./(4.*pi)*m1/rho_core)**(1./3.)
    !--Note the surface of the planet is located at the Plumber softening length
    eps_soft1 = r_surface
    if (eps_soft1 <= 0.) then
       print*,'Something wrong in the surface radius: eps_soft1 =',eps_soft1
    endif
 else
    npart_disc = nparttot
 endif

 !--set up an atmosphere around one of the binary masses (i.e. planet)
 if (surface_force .and. npart_planet_atm > 0) then
    call set_planet_atm(id,xyzh,vxyzu,npartoftype,maxvxyzu,itype,a0,R_in(1), &
                        H_R(1),m2,qindex(1),gamma,Ratm_in,Ratm_out,r_surface, &
                        nparttot,npart_planet_atm,npart_disc,hfact)
 endif

 !--move into the corotating frame with the planet
 if (surface_force .or. iexternalforce == iext_corotate) then
    call make_corotate(xyzh,vxyzu,a0,m2,npart,npart_disc)
 endif

 !--number of particles
 npart = nparttot

 call check_dust_method(dust_method,ichange_method)
 if (ichange_method .and. id==master) then
    np_dust = npart/5
    call write_setupfile(filename)
    print "(/,a)",' >>> please rerun the setup routine <<<'
    stop
 endif

 !
 ! print information about the angular momenta
 !
 ldisc = get_mean_angmom_vector(npart,xyzh,vxyzu)
 print "(a,'(',3(es10.2,1x),')')",' Disc specific angular momentum = ',ldisc
 if (nptmass > 1) then
    lcentral = get_mean_angmom_vector(nptmass,xyzmh_ptmass,vxyz_ptmass)
    print "(a,'(',3(es10.2,1x),')')",' Binary specific angular momentum = ',lcentral
    ! make unit vectors
    lcentral = lcentral/sqrt(dot_product(lcentral,lcentral))
    ldisc    = ldisc/sqrt(dot_product(ldisc,ldisc))
    print "(a,f6.1,a)",' Angle between disc and binary = ',acos(dot_product(lcentral,ldisc))*180./pi,' deg'
 endif

 !--alpha viscosity
 if (ndiscs==1) then
    alpha = alpha_returned(idisc)
 else
    call warning('setup_disc', &
       'multiple discs: cannot use alpha_AV for alpha_SS, setting equal to 0.1')
    alpha = 0.1
 endif

 !
 !--dust
 !
 if (use_dust) then
    if (use_dustgrowth) then !-size and dens already in code units
       dustprop(1,:) = grainsize(1)
       dustprop(2,:) = graindens(1)
    endif
 endif
 if (maxdiscs > 1 .and. ibinary==1) then
    !--circumprimary in flyby
    i = 2
 else
    !--single disc or circumbinary
    i = idisc
 endif
 R = (R_in(i) + R_out(i))/2
 Sigma = sig_norm(i)*scaled_sigma(R,sigmaprofilegas(i),pindex(i),R_ref(i),R_in(i),R_c(i))
 if (use_dust) then
    Sigmadust = sig_normdust(i)*scaled_sigma(R,sigmaprofiledust(i),pindex_dust(i),R_ref(i),R_indust(i),R_c_dust(i))
    Stokes = 0.5*pi*graindens*grainsize/(Sigma+Sigmadust)
    print "(a,i2,a)",' -------------- added dust --------------'
    if (use_dustgrowth) then
       print "(a,g10.3,a)", ' initial grain size: ',grainsize(1)*udist,' cm'
    else
       varstring = 'grain size'
       call make_tags_unique(ndusttypes,varstring)
       do i=1,ndusttypes
          print*,adjustr(varstring(i))//' : ',grainsize(i)*udist,' cm'
       enddo
       varstring = 'grain density'
       call make_tags_unique(ndusttypes,varstring)
       do i=1,ndusttypes
          print*,adjustr(varstring(i))//' : ',graindens(i)*umass/udist**3,' g/cm^3'
       enddo
    endif
    varstring = 'approx. Stokes'
    call make_tags_unique(ndusttypes,varstring)
    do i=1,ndusttypes
       print*,'',adjustr(varstring(i))//' : ',Stokes(i)
    enddo
    print "(1x,40('-'),/)"
 else
    print "(/,a,/)",' There is no dust here!'
 endif

 !
 !--planets
 !
 period_longest = 0.
 if (setplanets==1) then
    print "(a,i2,a)",' --------- added ',nplanets,' planets ------------'
    do i=1,nplanets
       nptmass = nptmass + 1
       phi = 0.
       phi = phi*pi/180.
       cosphi = cos(phi)
       sinphi = sin(phi)
       disc_m_within_r = 0.
       !--disc mass correction
       do j=1,npart
          r2 = xyzh(1,j)**2 + xyzh(2,j)**2 + xyzh(3,j)**2
          if (r2 < rplanet(i)**2) then
             itype = iamtype(iphase(j))
             disc_m_within_r = disc_m_within_r + massoftype(itype)
          endif
       enddo
       !--inner planet mass correction
       if (nplanets>1) then
          do j=1,nplanets
             if (rplanet(j)<rplanet(i)) disc_m_within_r = disc_m_within_r + mplanet(j)*jupiterm/umass
          enddo
       endif
       !--set sink particles
       Hill(i) = (mplanet(i)*jupiterm/solarm/(3.*mcentral))**(1./3.) * rplanet(i)
       if (nsinks == 2) then
          dist_bt_sinks = sqrt(dot_product(xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,2)))
          if (rplanet(i) > dist_bt_sinks) Hill(i) = (mplanet(i)*jupiterm/solarm/(3.*m1))**(1./3.) * rplanet(i)
       endif
       xyzmh_ptmass(1:3,nptmass)    = (/rplanet(i)*cosphi,rplanet(i)*sinphi,0./)
       xyzmh_ptmass(4,nptmass)      = mplanet(i)*jupiterm/umass
       xyzmh_ptmass(ihacc,nptmass)  = accrplanet(i)*Hill(i)
       xyzmh_ptmass(ihsoft,nptmass) = accrplanet(i)*Hill(i)
       vphi                         = sqrt((mcentral + disc_m_within_r)/rplanet(i))
       if (nsinks == 2 .and. rplanet(i) < dist_bt_sinks) vphi = sqrt((m1 + disc_m_within_r)/rplanet(i))
       vxyz_ptmass(1:3,nptmass)     = (/-vphi*sinphi,vphi*cosphi,0./)
       if (nsinks == 2 .and. rplanet(i) < dist_bt_sinks) then
          vxyz_ptmass(1:3,nptmass) = vxyz_ptmass(1:3,nptmass) + vxyz_ptmass(1:3,1)
          xyzmh_ptmass(1:3,nptmass) = xyzmh_ptmass(1:3,nptmass) + xyzmh_ptmass(1:3,1)
       endif

       !--incline positions and velocities
       inclplan(i) = inclplan(i)*pi/180.
       u = (/-sin(phi),cos(phi),0./)
       call rotatevec(xyzmh_ptmass(1:3,nptmass),u,-inclplan(i))
       call rotatevec(vxyz_ptmass(1:3,nptmass), u,-inclplan(i))
       !--print planet information
       omega = vphi/rplanet(i)
       print "(a,i2,a)",             ' >>> planet ',i,' <<<'
       print "(a,g10.3,a)",          ' orbital radius: ',rplanet(i)*udist/au,' AU'
       print "(a,g10.3,a,2pf7.3,a)", '          M(<R): ',(disc_m_within_r + mcentral)*umass/solarm, &
                                     ' MSun, disc mass correction is ',disc_m_within_r/mcentral,'%'
       print "(a,g10.3,a)",          '    planet mass: ',mplanet(i),' MJup'
       print "(a,g10.3,a)",          '    planet mass: ',mplanet(i)*jupiterm/earthm,' MEarth'
       print "(a,g10.3,a)",          '    planet mass: ',mplanet(i)*jupiterm/solarm,' MSun'
       print "(a,2(g10.3,a))",       ' orbital period: ',2.*pi*rplanet(i)/vphi*utime/years,' years or ', &
                                                 2*pi*rplanet(i)/vphi,' in code units'
       print "(a,g10.3,a)",          '    Hill radius: ',Hill(i),' AU'
       print "(a,g10.3,a,i3,a)",     '   accr. radius: ',xyzmh_ptmass(ihacc,nptmass),' AU or ', &
                                                 int(100*xyzmh_ptmass(ihacc,nptmass)/Hill(i)), ' % of Hill radius'
       print "(a,g10.3,a)",          '     resonances:'
       print "(a,g10.3,a)",   '    3:1 : ',(sqrt(mcentral)/(3.*omega))**(2./3.)*udist/au,' AU'
       print "(a,g10.3,a)",   '    4:1 : ',(sqrt(mcentral)/(4.*omega))**(2./3.)*udist/au,' AU'
       print "(a,g10.3,a)",   '    5:1 : ',(sqrt(mcentral)/(5.*omega))**(2./3.)*udist/au,' AU'
       print "(a,g10.3,a)",   '    9:1 : ',(sqrt(mcentral)/(9.*omega))**(2./3.)*udist/au,' AU'
       !--check planet accretion radii
       if (accrplanet(i) < 0.05) then
          call warning('setup_disc','accretion radius of planet < 1/20 Hill radius: unnecessarily small')
       elseif (accrplanet(i) > 0.5) then
          call warning('setup_disc','accretion radius of planet > Hill radius: too large')
       elseif(accrplanet(i)*Hill(i) > accr1) then
          call warning('setup_disc','accretion radius of planet > accretion radius of primary star: this is unphysical')
       endif
       print *, ''
       !--determine longest period
       period_longest = max(period_longest, 2.*pi/omega)
    enddo
    print "(1x,45('-'))"
    print *, ''
 endif

 !
 !--reset centre of mass to the origin
 !
 if (iexternalforce == iext_corot_binary) then
    npart_recentre = npart-npart_planet_atm
 else
    npart_recentre = npart
 endif

 call reset_centreofmass(npart_recentre,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 !
 !--set tmax and dtmax
 !
 if (icentral==1 .and. nsinks==2 .and. ibinary==0) then
    !--bound binary
    period = sqrt(4.*pi**2*binary_a**3/mcentral)
 elseif (icentral==1 .and. nsinks==2 .and. ibinary==1) then
    !--unbound binary (flyby)
    period = get_T_flyby(m1,m2,flyby_a,flyby_d)
 elseif (setplanets==1) then
    !--outer planet set above
    period = period_longest
 elseif (iwarp(idisc)) then
    !--warp radius
    period = sqrt(4.*pi**2*R_warp(idisc)**3/mcentral)
 else
    !--outer disc
    period = sqrt(4.*pi**2*R_out(idisc)**3/mcentral)
 endif
 if (period > 0.) then
    if (deltat > 0.) dtmax = deltat*period
    if (norbits >= 0) tmax = norbits*period
 endif

 !
 !--remind user to check for warnings and errors
 !
 print "(/,a)",' + ----------------------------------------------------- +'
 print "(a)",  ' |                                                       |'
 print "(a)",  ' |   please check output above for WARNINGS and ERRORS   |'
 print "(a)",  ' |   before starting the calculation                     |'
 print "(a)",  ' |                                                       |'
 print "(a,/)",' + ----------------------------------------------------- +'

 return
end subroutine setpart

!------------------------------------------------------------------------
!
!  prompt user for desired setup options
!
!------------------------------------------------------------------------
subroutine setup_interactive()
 use growth,           only:ifrag,isnow,rsnow,Tsnow,vfragSI,vfraginSI,vfragoutSI,gsizemincgs
 use io,               only:warning
 use prompting,        only:prompt
 use set_dust_options, only:set_dust_interactively
 integer :: maxdiscs
 real    :: disc_mfac(3)
 logical :: questplanets

 !
 !--units
 !
 dist_unit = 'au'
 mass_unit = 'solarm'
 !
 !--set defaults for central object(s)
 !
 icentral = 1
 print "(a)",'==========================='
 print "(a)",'+++  CENTRAL OBJECT(S)  +++'
 print "(a)",'==========================='
 call prompt('Do you want to use sink particles or an external potential?'// &
             new_line('A')//' 0=potential'//new_line('A')//' 1=sinks'// &
             new_line('A'),icentral,0,1)
 select case (icentral)
 case (0)
    !--external potential
    ipotential = 1
    call prompt('Which potential?'//new_line('A')// &
               ' 1=central point mass'//new_line('A')// &
               ' 2=binary potential'//new_line('A')// &
               ' 3=spinning black hole'//new_line('A'),ipotential,1,3)
    select case (ipotential)
    case (1)
       !--point mass
       iexternalforce = iext_star
       m1       = 1.
       accr1    = 1.
    case (2)
       !--fixed binary
       surface_force = .false.
       call prompt('Do you want model a surface on one mass (i.e. planet)?',surface_force)
       if (surface_force) then
          iexternalforce = iext_corot_binary
          !--binary
          m1       = 0.001 ! Planet
          m2       = 1.    ! Central star
          accr1    = 0.    ! Surface force, so we don't need an accretion radius
          accr2    = 0.1

          call prompt('Enter orbital radius of the planet (e.g. 5.2au)',dist_unit)
          call prompt('Enter average core density for the planet (g/cm^3)',rho_core_cgs,0.)
          call prompt('Enter inner atmosphere radius in planet radii',Ratm_in,1.,10.)
          call prompt('Enter outer atmosphere radius in planet radii',Ratm_out,Ratm_in,10.)
          call prompt('Enter atmosphere type (1:r**(-3); 2:r**(-1./(gamma-1.)))',atm_type,1,2)
          call prompt('Enter fraction of particles to be used in planet atmosphere',Natmfrac,0.,1.)
          call prompt('Do you want to ramp up the planet mass slowly?',ramp)
       else
          iexternalforce = iext_binary
          m1       = 1.
          m2       = 1.
          accr1    = 1.
          accr2    = 1.
       endif
    case (3)
       !--spinning black hole (Lense-Thirring)
       iexternalforce = iext_lensethirring
       call prompt('Include Einstein precession?',einst_prec)
       if (einst_prec) iexternalforce = iext_einsteinprec
       m1          = 1.
       accr1       = 30.
       bhspin      = 1.
       bhspinangle = 0.
    end select
 case (1)
    !--sink particle(s)
    nsinks = 1
    call prompt('How many sinks?',nsinks,1,2)
    select case (nsinks)
    case (1)
       !--single star
       m1       = 1.
       accr1    = 1.
    case (2)
       !--binary
       ibinary = 0
       call prompt('Do you want the binary orbit to be bound (elliptic) or'// &
                  ' unbound (parabolic/hyperbolic) [flyby]?'//new_line('A')// &
                  ' 0=bound'//new_line('A')//' 1=unbound'//new_line('A'),ibinary,0,1)
       select case (ibinary)
       case (0)
          !--bound
          m1       = 1.
          m2       = 0.2
          binary_a = 10.
          binary_e = 0.
          binary_i = 0.
          binary_O = 0.
          binary_w = 270.
          binary_f = 180.
          accr1    = 1.
          accr2    = 0.5
       case (1)
          !--unbound (flyby)
          m1       = 1.
          m2       = 1.
          accr1    = 1.
          accr2    = 1.
          flyby_a  = 200.
          flyby_d  = 10.
          flyby_O  = 0.
          flyby_i  = 0.
       end select
    end select
 end select
!
!--multiple disc options
!
 print "(/,a)",'================='
 print "(a)",  '+++  DISC(S)  +++'
 print "(a)",  '================='
 iuse_disc = .false.
 iuse_disc(1) = .true.
 ndiscs = 1
 maxdiscs = 1
 if ((icentral==1) .and. (nsinks>=2)) then
    !--multiple discs possible
    if (ibinary==0) then
       !--bound binary: circum-binary, -primary, -secondary
       maxdiscs = 3
    elseif (ibinary==1) then
       !--unbound binary (flyby): circum-primary, -secondary
       maxdiscs = 2
       iuse_disc(2) = .true.
       iuse_disc(1) = .false.
    endif
    do i=4-maxdiscs,3
       call prompt('Do you want a circum'//trim(disctype(i))//' disc?',iuse_disc(i))
    enddo
    if (.not.any(iuse_disc)) iuse_disc(1) = .true.
    !--set number of discs
    ndiscs = count(iuse_disc)
    if (ndiscs > 1) then
       use_global_iso = .false.
    endif
 endif
!
!--set gas disc defaults
!
 R_in       = accr1
 R_out      = 150.
 R_ref      = R_in
 R_c        = R_out
 R_warp     = 0.
 H_warp     = 0.
 mass_set   = 0
 itapergas  = .false.
 ismoothgas = .true.
 iwarp      = .false.
 pindex     = 1.
 qindex     = 0.25
 if (ndiscs > 1) qindex = 0.
 if (maxalpha==0) alphaSS    = 0.005
 posangl    = 0.
 incl       = 0.
 H_R        = 0.05
 disc_mfac  = 1.
 if (surface_force) then
    R_in       = 0.1
    R_out      = 3.
    R_ref      = 1.
 endif
 if ((icentral==1 .and. nsinks>=2) .and. (ibinary==0)) then
    !--don't smooth circumbinary, by default
    ismoothgas(1) = .false.
    !--set appropriate disc radii for bound binary
    R_in      = (/2.5*binary_a, accr1, accr2/)
    R_out     = (/5.*R_in(1), 5.*accr1, 5.*accr2/)
    R_ref     = R_in
    R_c       = R_out
    disc_mfac = (/1., 0.1, 0.01/)
    if (ndiscs > 1) then
       !--set H/R so temperature is globally constant
       call prompt('Do you want a globally isothermal disc (if not Farris et al. 2014)?',use_global_iso)
       !--------------------------------
       ! N.B. The initializations of multiple discs is not done using the implementation of the eos
       ! a radial profile centred on CM, primary and secondary is used.
       ! The value of H_R used in setpart to set cs0 is the one of the circumbinary if cb disc is present,
       ! otherwise it uses the circumprimary.
       ! The values of H_R used for the other discs are set using the equations below, however changing them here
       ! is not enough. THey need to be changed also in the the setpart function.
       !--------------------------------
       if(.not. use_global_iso) then
          call prompt('Enter q_index',qindex(1))
          qindex=qindex(1)
          if(iuse_disc(1)) then
             call prompt('Enter H/R of circumbinary at R_ref',H_R(1))
             H_R(2) = (R_ref(2)/R_ref(1)*(m1+m2)/m1)**(0.5-qindex(1)) * H_R(1)
             H_R(3) = (R_ref(3)/R_ref(1)*(m1+m2)/m2)**(0.5-qindex(1)) * H_R(1)
          else
             if(iuse_disc(2))then
                call prompt('Enter H/R of circumprimary at R_ref',H_R(2))
                H_R(1) = (R_ref(1)/R_ref(2)*m1/(m1+m2))**(0.5-qindex(2)) * H_R(2)
                H_R(3) = (R_ref(3)/R_ref(2)*m2/m1)**(0.5-qindex(2)) * H_R(2)
             else
                call prompt('Enter H/R of circumsecondary at R_ref',H_R(3))
                H_R(1) = sqrt(R_ref(1)/R_ref(3)*m2/(m1+m2))**(0.5-qindex(3)) * H_R(3)
                H_R(2) = sqrt(R_ref(2)/R_ref(3)*m2/m1)**(0.5-qindex(3)) * H_R(3)
             endif
          endif
          H_R(2) = nint(H_R(2)*10000.)/10000.
          H_R(3) = nint(H_R(3)*10000.)/10000.
       else
          if (iuse_disc(1)) then
             H_R(2) = sqrt(R_ref(2)/R_ref(1)*(m1+m2)/m1) * H_R(1)
             H_R(3) = sqrt(R_ref(3)/R_ref(1)*(m1+m2)/m2) * H_R(1)
             qindex(2) = qindex(1)
             qindex(3) = qindex(1)
             call warning('setup_disc','using circumbinary (H/R)_ref to set global temperature')
          elseif (iuse_disc(2)) then
             H_R(3) = sqrt(R_ref(3)/R_ref(2)*m1/m2) * H_R(2)
             qindex(3) = qindex(2)
             call warning('setup_disc','using circumprimary (H/R)_ref to set global temperature')
          endif
          H_R(2) = nint(H_R(2)*10000.)/10000.
          H_R(3) = nint(H_R(3)*10000.)/10000.
       endif
    endif
 endif
 do i=1,3
    if (iuse_disc(i)) then
       if (ndiscs > 1) print "(/,a)",' >>>  circum'//trim(disctype(i))//' disc  <<<'
       call prompt('How do you want to set the gas disc mass?'//new_line('A')// &
                  ' 0=total disc mass'//new_line('A')// &
                  ' 1=mass within annulus'//new_line('A')// &
                  ' 2=surface density normalisation'//new_line('A')// &
                  ' 3=surface density at reference radius'//new_line('A')// &
                  ' 4=minimum Toomre Q'//new_line('A'),mass_set(i),0,4)
       call prompt('Do you want to exponentially taper the outer gas disc profile?',itapergas(i))
       call prompt('Do you want to warp the disc?',iwarp(i))
       select case (mass_set(i))
       case (0)
          disc_m(i)    = 0.05   * disc_mfac(i)
       case (1)
          annulus_m(i) = 0.05   * disc_mfac(i)
          R_inann(i)   = R_in(i)
          R_outann(i)  = R_out(i)
       case (2)
          sig_norm(i)  = 1.E-02 * disc_mfac(i)
       case (3)
          sig_ref(i)   = 1.E-02 * disc_mfac(i)
       case (4)
          Q_min(i)     = 1.0
       end select
       if (iwarp(i)) then
          R_warp = 0.5*(R_in + R_out)
          H_warp = 20.
          incl   = 30.
       endif
    endif
 enddo
!
!--set dust disc defaults
!
 if (use_dust) then
    R_indust       = R_in
    R_outdust      = R_out
    pindex_dust    = pindex
    qindex_dust    = qindex
    H_R_dust       = H_R
    itaperdust     = itapergas
    ismoothdust    = ismoothgas
    R_c_dust       = R_c
    print "(/,a)",'=============='
    print "(a)",  '+++  DUST  +++'
    print "(a)",  '=============='
    call set_dust_interactively()
    if (use_dustgrowth .and. dust_method == 2) then
       print "(/,a)",'================================'
       print "(a)",  '+++  GROWTH & FRAGMENTATION  +++'
       print "(a)",  '================================'
       !
       !--set growth parameters default
       !
       ifrag = 1
       isnow = 0
       rsnow = 100.
       Tsnow = 20.
       vfragSI = 15.
       vfraginSI = 5.
       vfragoutSI = 15.
       gsizemincgs = 1.e-3
       !
       !--growth parameters from user
       !
       call prompt('Enter fragmentation model (0=off,1=on,2=Kobayashi)',ifrag,-1,2)
       select case(ifrag)
       case(0)
          print "(a)",'-----------'
          print "(a)",'Pure growth'
          print "(a)",'-----------'
       case(1)
          print "(a)",'----------------------'
          print "(a)",'Growth + fragmentation'
          print "(a)",'----------------------'
       case(2)
          print "(a)",'----------------------------------------'
          print "(a)",'Growth + Kobayashi`s fragmentation model'
          print "(a)",'----------------------------------------'
       case default
       end select
       if (ifrag > 0) then
          call prompt('Enter minimum allowed grain size in cm',gsizemincgs)
          call prompt('Do you want a snow line ? (0=no,1=position based,2=temperature based)',isnow,0,2)
          if (isnow == 0) then
             call prompt('Enter uniform vfrag in m/s',vfragSI,1.)
          else
             if (isnow == 1) call prompt('How far from the star in AU ?',rsnow,0.)
             if (isnow == 2) call prompt('Enter snow line condensation temperature in K',Tsnow,0.)
             call prompt('Enter inward vfragin in m/s',vfraginSI,1.)
             call prompt('Enter outward vfragout in m/s',vfragoutSI,1.)
          endif
       endif
    elseif (use_dustgrowth .and. dust_method == 1) then
       print "(a)",'growth and fragmentation not available for one fluid method'
    endif
 endif
!
!--resolution
!
 np = 500000
 if (use_dust .and. .not.use_dustfrac) then
    np_dust = np/ndusttypesinp/5
 else
    np_dust = 0
 endif
!
!--add planets
!
 questplanets  = .false.
 setplanets    = 0
 nplanets      = 0
 mplanet       = 1.
 rplanet       = (/ (10.*i, i=1,maxplanets) /)
 accrplanet    = 0.25
 inclplan      = 0.
 print "(/,a)",'================='
 print "(a)",  '+++  PLANETS  +++'
 print "(a)",  '================='
 call prompt('Do you want to add planets?',questplanets)
 if (questplanets) then
    setplanets = 1
    nplanets   = 1
    call prompt('Enter the number of planets',nplanets,1,maxplanets)
 endif
!
!--determine simulation time
!
 print "(/,a)",'================'
 print "(a)",  '+++  OUTPUT  +++'
 print "(a)",  '================'
 deltat  = 0.1
 norbits = 100
 if (setplanets==1) then
    call prompt('Enter time between dumps as fraction of outer planet period',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 else if (icentral==1 .and. nsinks==2 .and. ibinary==0) then
    call prompt('Enter time between dumps as fraction of binary period',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 else if (icentral==1 .and. nsinks==2 .and. ibinary==1) then
    deltat  = 0.01
    norbits = 1
    call prompt('Enter time between dumps as fraction of flyby time',deltat,0.)
 elseif (any(iwarp)) then
    call prompt('Enter time between dumps as fraction of orbital time at warp',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 else
    call prompt('Enter time between dumps as fraction of outer disc orbital time',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 endif

end subroutine setup_interactive

!------------------------------------------------------------------------
!
! write setup file
!
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils,     only:write_inopt
 use set_dust_options, only:write_dust_setup_options
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 logical :: done_alpha
 integer :: i,maxdiscs
 character(len=20) :: varstring(maxdusttypes)

 done_alpha = .false.
 maxdiscs = 1
 if ((icentral==1) .and. (nsinks==2)) maxdiscs = 3

 print "(/,a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for disc setup routine'
 !--resolution
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','number of gas particles',iunit)
 if (use_dust .and. .not.use_dustfrac) then
    varstring = 'np_dust'
    call make_tags_unique(ndusttypesinp,varstring)
    do i=1,ndusttypesinp
       call write_inopt(np_dust(i),varstring(i),'number of dust particles',iunit)
    enddo
 endif
 !--units
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au,pc,kpc,0.1pc)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm,jupiterm,earthm)',iunit)
 !--central objects(s)/potential
 write(iunit,"(/,a)") '# central object(s)/potential'
 call write_inopt(icentral,'icentral', &
    'use sink particles or external potential (0=potential,1=sinks)',iunit)
 select case (icentral)
 case (0)
    !--external potential
    call write_inopt(ipotential,'ipotential','potential (1=central point mass,'// &
                     '2=binary potential,3=spinning black hole)',iunit)
    select case (ipotential)
    case (1)
       !--point mass
       call write_inopt(m1,'m1','star mass',iunit)
       call write_inopt(accr1,'accr1','star accretion radius',iunit)
    case (2)
       !--fixed binary
       call write_inopt(m1,'m1','primary mass',iunit)
       call write_inopt(m2,'m2','secondary mass',iunit)
       call write_inopt(accr1,'accr1','primary accretion radius',iunit)
       call write_inopt(accr2,'accr2','secondary accretion radius',iunit)

       !--options of planet surface/atmosphere
       write(iunit,"(/,a)") '# options for planet surface/atmosphere'
       call write_inopt(surface_force,'surface_force','model m1 as planet with surface',iunit)
       if (surface_force) then
          call write_inopt(rho_core_cgs,'rho_core','planet core density (cgs units)',iunit)
          call write_inopt(Ratm_in,'Ratm_in','inner atmosphere radius (planet radii)',iunit)
          call write_inopt(Ratm_out,'Ratm_out','outer atmosphere radius (planet radii)',iunit)
          call write_inopt(atm_type,'atm_type','atmosphere type (1:r**(-3); 2:r**(-1./(gamma-1.)))',iunit)
          call write_inopt(Natmfrac,'Natm/Npart','fraction of particles for planet atmosphere',iunit)
          call write_inopt(ramp,'ramp','Do you want to ramp up the planet mass slowly?',iunit)
          if (.not.ramp .and. Natmfrac == 0.) then
             print*,'Warning! Not ramping the mass or initialising an atmosphere will'// &
                     'likely cause explosive collisions between particles'
          elseif (ramp .and. Natmfrac /= 0.) then
             print*,'Warning! The atmosphere will be lost while ramping up the planet mass...'
             print*,'         ...try using one or the other'
          endif
       endif
    case (3)
       !--spinning black hole: Lense-Thirring (+ Einstein precession)
       call write_inopt(einst_prec,'einst_prec','include Einstein precession',iunit)
       call write_inopt(m1,'m1','black hole mass',iunit)
       call write_inopt(accr1,'accr1','black hole accretion radius',iunit)
       call write_inopt(bhspin,'bhspin','black hole spin',iunit)
       call write_inopt(bhspinangle,'bhspinangle','black hole spin angle (deg)',iunit)
    end select
 case (1)
    !--sink particle(s)
    call write_inopt(nsinks,'nsinks','number of sinks',iunit)
    select case (nsinks)
    case (1)
       !--single star
       write(iunit,"(/,a)") '# options for central star'
       call write_inopt(m1,'m1','star mass',iunit)
       call write_inopt(accr1,'accr1','star accretion radius',iunit)
    case (2)
       !--binary
       call write_inopt(ibinary,'ibinary','binary orbit (0=bound,1=unbound [flyby])',iunit)
       select case (ibinary)
       case (0)
          !--bound
          write(iunit,"(/,a)") '# options for binary'
          call write_inopt(m1,'m1','primary mass',iunit)
          call write_inopt(m2,'m2','secondary mass',iunit)
          call write_inopt(binary_a,'binary_a','binary semi-major axis',iunit)
          call write_inopt(binary_e,'binary_e','binary eccentricity',iunit)
          call write_inopt(binary_i,'binary_i','i, inclination (deg)',iunit)
          call write_inopt(binary_O,'binary_O','Omega, PA of ascending node (deg)',iunit)
          call write_inopt(binary_w,'binary_w','w, argument of periapsis (deg)',iunit)
          call write_inopt(binary_f,'binary_f','f, initial true anomaly (deg,180=apastron)',iunit)
          call write_inopt(accr1,'accr1','primary accretion radius',iunit)
          call write_inopt(accr2,'accr2','secondary accretion radius',iunit)
       case (1)
          !--unbound (flyby)
          !--central star
          write(iunit,"(/,a)") '# options for central star'
          call write_inopt(m1,'m1','central star mass',iunit)
          call write_inopt(accr1,'accr1','central star accretion radius',iunit)
          !--perturber
          write(iunit,"(/,a)") '# options for perturber'
          call write_inopt(m2,'m2','perturber mass',iunit)
          call write_inopt(accr2,'accr2','perturber accretion radius',iunit)
          call write_inopt(flyby_a,'flyby_a','distance of minimum approach',iunit)
          call write_inopt(flyby_d,'flyby_d','initial distance (units of dist. min. approach)',iunit)
          call write_inopt(flyby_O,'flyby_O','position angle of ascending node (deg)',iunit)
          call write_inopt(flyby_i,'flyby_i','inclination (deg)',iunit)
       end select
    end select
 end select
 !--multiple disc options
 if (maxdiscs > 1) then
    write(iunit,"(/,a)") '# options for multiple discs'
    do i=1,3
       call write_inopt(iuse_disc(i),'use_'//trim(disctype(i))//'disc','setup circum' &
                                     //trim(disctype(i))//' disc',iunit)
    enddo
    call write_inopt(use_global_iso,'use_global_iso',&
                     'globally isothermal or Farris et al. (2014)',iunit)
 endif
 !--individual disc(s)
 do i=1,3
    if (iuse_disc(i)) then
       if (maxdiscs > 1) then
          disclabel = disctype(i)
       else
          disclabel = ''
       endif
       !--gas disc
       if (maxdiscs > 1) then
          write(iunit,"(/,a)") '# options for circum'//trim(disclabel)//' gas disc'
       else
          write(iunit,"(/,a)") '# options for gas accretion disc'
       endif
       call write_inopt(mass_set(i),'mass_set'//trim(disclabel),'how to set gas density profile' // &
          ' (0=total disc mass,1=mass within annulus,2=surface density normalisation,' // &
          '3=surface density at reference radius,4=minimum Toomre Q)',iunit)
       call write_inopt(itapergas(i),'itapergas'//trim(disclabel), &
          'exponentially taper the outer disc profile',iunit)
       call write_inopt(ismoothgas(i),'ismoothgas'//trim(disclabel),'smooth inner disc',iunit)
       call write_inopt(iwarp(i),'iwarp'//trim(disclabel),'warp disc',iunit)
       call write_inopt(R_in(i),'R_in'//trim(disclabel),'inner radius',iunit)
       call write_inopt(R_ref(i),'R_ref'//trim(disclabel),'reference radius',iunit)
       call write_inopt(R_out(i),'R_out'//trim(disclabel),'outer radius',iunit)
       if (itapergas(i)) call write_inopt(R_c(i),'R_c'//trim(disclabel), &
          'characteristic radius of the exponential taper',iunit)
       select case (mass_set(i))
       case (0)
          call write_inopt(disc_m(i),'disc_m'//trim(disclabel),'disc mass',iunit)
       case (1)
          call write_inopt(annulus_m(i),'annulus_m'//trim(disclabel),'mass within annulus',iunit)
          call write_inopt(R_inann(i),'R_inann'//trim(disclabel),'inner annulus radius',iunit)
          call write_inopt(R_outann(i),'R_outann'//trim(disclabel),'outer annulus radius',iunit)
       case (2)
          if (itapergas(i)) then
             call write_inopt(sig_norm(i),'sig_norm'//trim(disclabel), &
                'sigma = sig_norm (R/R_ref)^-p exp[-(R/R_c)^(2-p)] (1-sqrt(R_in/R))',iunit)
          else
             call write_inopt(sig_norm(i),'sig_norm'//trim(disclabel), &
                'sigma = sig_norm (R/R_ref)^-p (1-sqrt(R_in/R))',iunit)
          endif
       case (3)
          call write_inopt(sig_ref(i),'sig_ref'//trim(disclabel),'sigma at reference radius',iunit)
       case (4)
          call write_inopt(Q_min(i),'Q_min'//trim(disclabel),'minimum Toomre Q',iunit)
       end select
       call write_inopt(pindex(i),'pindex'//trim(disclabel),'p index',iunit)
       call write_inopt(qindex(i),'qindex'//trim(disclabel),'q index',iunit)
       call write_inopt(posangl(i),'posangl'//trim(disclabel),'position angle (deg)',iunit)
       call write_inopt(incl(i),'incl'//trim(disclabel),'inclination (deg)',iunit)
       call write_inopt(H_R(i),'H_R'//trim(disclabel),'H/R at R=R_ref',iunit)
       if (iwarp(i)) then
          call write_inopt(R_warp(i),'R_warp'//trim(disclabel),'warp radius',iunit)
          call write_inopt(H_warp(i),'H_warp'//trim(disclabel),'warp smoothing length',iunit)
       endif
       if (.not.done_alpha) then
          if (maxalpha==0) call write_inopt(alphaSS,'alphaSS','desired alphaSS',iunit)
          done_alpha = .true.
       endif
       !--dust disc
       if (use_dust .and. (iprofile_dust == 1 .or. iprofile_dust == 2)) then
          if (maxdiscs > 1) then
             write(iunit,"(/,a)") '# options for circum'//trim(disclabel)//' dust disc'
          else
             write(iunit,"(/,a)") '# options for dust accretion disc'
          endif
          call write_inopt(itaperdust(i),'itaperdust'//trim(disclabel), &
             'exponentially taper the outer disc profile',iunit)
          call write_inopt(ismoothdust(i),'ismoothdust'//trim(disclabel),'smooth inner disc',iunit)
          call write_inopt(R_indust(i),'R_indust'//trim(disclabel),'inner radius',iunit)
          call write_inopt(R_outdust(i),'R_outdust'//trim(disclabel),'outer radius',iunit)
          if (itaperdust(i)) call write_inopt(R_c_dust(i),'R_c_dust'//trim(disclabel), &
             'characteristic radius of the exponential taper',iunit)
          call write_inopt(pindex_dust(i),'pindex_dust'//trim(disclabel),'p index',iunit)
          call write_inopt(qindex_dust(i),'qindex_dust'//trim(disclabel),'q index',iunit)
          call write_inopt(H_R_dust(i),'H_R_dust'//trim(disclabel),'H/R at R=R_ref',iunit)
       endif
    endif
 enddo
 !--dust & growth options
 if (use_dust) then
    call write_dust_setup_options(iunit)
 endif
 !--planets
 write(iunit,"(/,a)") '# set planets'
 call write_inopt(setplanets,'setplanets','add planets? (0=no,1=yes)',iunit)
 if (setplanets==1) then
    call write_inopt(nplanets,'nplanets','number of planets',iunit)
    do i=1,nplanets
       write(iunit,"(/,a)") '# planet:'//trim(planets(i))
       call write_inopt(mplanet(i),'mplanet'//trim(planets(i)),'planet mass (in Jupiter mass)',iunit)
       call write_inopt(rplanet(i),'rplanet'//trim(planets(i)),'planet distance from star',iunit)
       call write_inopt(inclplan(i),'inclplanet'//trim(planets(i)),'planet orbital inclination (deg)',iunit)
       call write_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),'planet accretion radius (in Hill radius)',iunit)
    enddo
 endif
 !--timestepping
 write(iunit,"(/,a)") '# timestepping'
 if (setplanets==1) then
    call write_inopt(norbits,'norbits','maximum number of outer planet orbits',iunit)
 else if (icentral==1 .and. nsinks==2 .and. ibinary==0) then
    call write_inopt(norbits,'norbits','maximum number of binary orbits',iunit)
 else
    call write_inopt(norbits,'norbits','maximum number of orbits at outer disc',iunit)
 endif
 call write_inopt(deltat,'deltat','output interval as fraction of orbital period',iunit)
#ifdef MCFOST
 !--mcfost
 write(iunit,"(/,a)") '# mcfost'
 call write_inopt(use_mcfost,'use_mcfost','use the mcfost library',iunit)
#endif

 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!
! read setup file
!
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils,     only:open_db_from_file,inopts,read_inopt,close_db
 use set_dust_options, only:read_dust_setup_options
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 character(len=20) :: varstring(maxdusttypes)
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)

 !--read old options for backwards compatibility
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(icentral,'icentral',db,err=ierr)
 if (ierr /= 0) then
    ! if old .setup file, read options and return error so it gets rewritten
    call read_obsolete_setup_options(db)
    ierr = 1
    print "(a)",' >>> re-writing obsolete .setup file: CHECK CAREFULLY <<<'
    return
 endif

 nerr = 0
 !--units
 call read_inopt(mass_unit,'mass_unit',db,errcount=nerr)
 call read_inopt(dist_unit,'dist_unit',db,errcount=nerr)
 !--central objects(s)/potential
 call read_inopt(icentral,'icentral',db,min=0,max=1,errcount=nerr)
 select case (icentral)
 case (0)
    !--external potential
    call read_inopt(ipotential,'ipotential',db,min=1,max=3,errcount=nerr)
    select case (ipotential)
    case (1)
       !--point mass
       iexternalforce = iext_star
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
    case (2)
       !--fixed binary
       iexternalforce = iext_binary
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
       call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)

       !--options of planet surface/atmosphere
       call read_inopt(surface_force,'surface_force',db,errcount=nerr)
       if (surface_force) then
          iexternalforce = iext_corot_binary
          call read_inopt(rho_core_cgs,'rho_core',db,min=0.,errcount=nerr)
          call read_inopt(Ratm_in,'Ratm_in',db,min=1.,errcount=nerr)
          call read_inopt(Ratm_out,'Ratm_out',db,min=1.,max=10.,errcount=nerr)
          call read_inopt(atm_type,'atm_type',db,min=1,max=2,errcount=nerr)
          call read_inopt(Natmfrac,'Natm/Npart',db,min=0.,max=1.,errcount=nerr)
          call read_inopt(ramp,'ramp',db,errcount=nerr)
       endif
    case (3)
       !--spinning black hole (Lense-Thirring)
       iexternalforce = iext_lensethirring
       call read_inopt(einst_prec,'einst_prec',db,errcount=nerr)
       if (einst_prec) iexternalforce = iext_einsteinprec
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
       call read_inopt(bhspin,'bhspin',db,min=0.,errcount=nerr)
       call read_inopt(bhspinangle,'bhspinangle',db,min=0.,errcount=nerr)
    end select
 case (1)
    !--sink particles
    call read_inopt(nsinks,'nsinks',db,min=1,max=2,errcount=nerr)
    select case (nsinks)
    case (1)
       !--single star
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
    case (2)
       !--binary
       call read_inopt(ibinary,'ibinary',db,min=0,max=1,errcount=nerr)
       select case (ibinary)
       case (0)
          !--bound
          call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
          call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
          call read_inopt(binary_a,'binary_a',db,errcount=nerr)
          call read_inopt(binary_e,'binary_e',db,min=0.,errcount=nerr)
          call read_inopt(binary_i,'binary_i',db,errcount=nerr)
          call read_inopt(binary_O,'binary_O',db,errcount=nerr)
          call read_inopt(binary_w,'binary_w',db,errcount=nerr)
          call read_inopt(binary_f,'binary_f',db,errcount=nerr)
          call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
          call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)
       case (1)
          !--unbound (flyby)
          !--central star
          call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
          call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
          !--perturber
          call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
          call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)
          call read_inopt(flyby_a,'flyby_a',db,min=0.,errcount=nerr)
          call read_inopt(flyby_d,'flyby_d',db,min=0.,errcount=nerr)
          call read_inopt(flyby_O,'flyby_O',db,min=0.,errcount=nerr)
          call read_inopt(flyby_i,'flyby_i',db,min=0.,errcount=nerr)
       end select
    end select
 end select

 !--dust
 if (use_dust) then
    call read_dust_setup_options(db,nerr)
    !--dust method
    select case(dust_method)
    case(1)
       use_dustfrac = .true.
       ilimitdustflux = ilimitdustfluxinp
       ndustsmall = ndusttypesinp
    case(2)
       use_dustfrac = .false.
       ndustlarge = ndusttypesinp
    end select
    ndusttypes = ndusttypesinp
 endif

 !--resolution
 call read_inopt(np,'np',db,min=0,errcount=nerr)
 if (use_dust .and. .not.use_dustfrac) then
    varstring = 'np_dust'
    call make_tags_unique(ndusttypesinp,varstring)
    do i=1,ndustlarge
       call read_inopt(np_dust(i),varstring(i),db,min=0,errcount=nerr)
    enddo
 endif

 !--multiple discs
 iuse_disc = .false.
 if ((icentral==1) .and. (nsinks==2)) then
    if (ibinary==0) then
       call read_inopt(iuse_disc(1),'use_'//trim(disctype(1))//'disc',db,errcount=nerr)
    endif
    do i=2,3
       call read_inopt(iuse_disc(i),'use_'//trim(disctype(i))//'disc',db,errcount=nerr)
    enddo
 else
    iuse_disc(1) = .true.
 endif
 ndiscs = count(iuse_disc)
 if (ndiscs > 1) then
    call read_inopt(use_global_iso,'use_global_iso',db,errcount=nerr)
 endif

 do i=1,3
    if (iuse_disc(i)) then
       if (nsinks == 2) then
          disclabel = disctype(i)
       else
          disclabel = ''
       endif
       !--gas disc
       call read_inopt(R_in(i),'R_in'//trim(disclabel),db,min=0.,errcount=nerr)
       call read_inopt(R_out(i),'R_out'//trim(disclabel),db,min=R_in(i),errcount=nerr)
       call read_inopt(R_ref(i),'R_ref'//trim(disclabel),db,min=R_in(i),errcount=nerr)
       call read_inopt(itapergas(i),'itapergas'//trim(disclabel),db,errcount=nerr)
       call read_inopt(ismoothgas(i),'ismoothgas'//trim(disclabel),db,errcount=nerr)
       call read_inopt(mass_set(i),'mass_set'//trim(disclabel),db,min=0,max=4,errcount=nerr)
       if (itapergas(i)) then
          call read_inopt(R_c(i),'R_c'//trim(disclabel),db,min=0.,errcount=nerr)
       endif
       select case (mass_set(i))
       case (0)
          call read_inopt(disc_m(i),'disc_m'//trim(disclabel),db,min=0.,errcount=nerr)
       case (1)
          call read_inopt(annulus_m(i),'annulus_m'//trim(disclabel),db,min=0.,errcount=nerr)
          call read_inopt(R_inann(i),'R_inann'//trim(disclabel),db,min=R_in(i),errcount=nerr)
          call read_inopt(R_outann(i),'R_outann'//trim(disclabel),db,min=R_in(i),errcount=nerr)
       case (2)
          call read_inopt(sig_norm(i),'sig_norm'//trim(disclabel),db,min=0.,errcount=nerr)
       case (3)
          call read_inopt(sig_ref(i),'sig_ref'//trim(disclabel),db,min=0.,errcount=nerr)
       case (4)
          call read_inopt(Q_min(i),'Q_min'//trim(disclabel),db,min=0.,errcount=nerr)
       end select
       call read_inopt(pindex(i),'pindex'//trim(disclabel),db,errcount=nerr)
       call read_inopt(qindex(i),'qindex'//trim(disclabel),db,errcount=nerr)
       call read_inopt(posangl(i),'posangl'//trim(disclabel),db,min=0.,max=360.,errcount=nerr)
       call read_inopt(incl(i),'incl'//trim(disclabel),db,min=0.,max=180.,errcount=nerr)
       call read_inopt(H_R(i),'H_R'//trim(disclabel),db,min=0.,errcount=nerr)
       call read_inopt(iwarp(i),'iwarp'//trim(disclabel),db,errcount=nerr)
       if (iwarp(i)) then
          call read_inopt(R_warp(i),'R_warp'//trim(disclabel),db,min=0.,errcount=nerr)
          call read_inopt(H_warp(i),'H_warp'//trim(disclabel),db,min=0.,errcount=nerr)
       endif
       !--dust disc
       if (use_dust) then
          select case (iprofile_dust)
          case (0)
             R_indust(i)    = R_in(i)
             R_outdust(i)   = R_out(i)
             pindex_dust(i) = pindex(i)
             qindex_dust(i) = qindex(i)
             H_R_dust(i)    = H_R(i)
             itaperdust(i)  = itapergas(i)
             ismoothdust(i) = ismoothgas(i)
             R_c_dust(i)    = R_c(i)
          case (1,2)
             call read_inopt(R_indust(i),'R_indust'//trim(disclabel),db,min=R_in(i),err=ierr,errcount=nerr)
             if (ierr /= 0) R_indust(i) = R_in(i)

             call read_inopt(R_outdust(i),'R_outdust'//trim(disclabel),db,min=R_indust(i),max=R_out(i),err=ierr,errcount=nerr)
             if (ierr /= 0) R_outdust(i) = R_out(i)
             call read_inopt(pindex_dust(i),'pindex_dust'//trim(disclabel),db,err=ierr,errcount=nerr)
             if (ierr /= 0) pindex_dust(i) = pindex(i)
             call read_inopt(itaperdust(i),'itaperdust'//trim(disclabel),db,err=ierr,errcount=nerr)
             call read_inopt(ismoothdust(i),'ismoothdust'//trim(disclabel),db,err=ierr,errcount=nerr)
             if (itaperdust(i)) then
                call read_inopt(R_c_dust(i),'R_c_dust'//trim(disclabel),db,min=0.,err=ierr,errcount=nerr)
                if (ierr /= 0) R_c_dust(i) = R_c(i)
             endif
             call read_inopt(qindex_dust(i),'qindex_dust'//trim(disclabel),db,min=qindex(i),err=ierr,errcount=nerr)
             if (ierr /= 0) qindex_dust(i) = qindex(i)
             call read_inopt(H_R_dust(i),'H_R_dust'//trim(disclabel),db,min=0.,max=H_R(i),err=ierr,errcount=nerr)
             if (ierr /= 0) H_R_dust(i) = H_R(i)
          end select
       endif
    endif
 enddo
 if (maxalpha==0) call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
 !--planets
 call read_inopt(setplanets,'setplanets',db,min=0,max=1,errcount=nerr)
 if (setplanets==1) then
    call read_inopt(nplanets,'nplanets',db,min=0,max=maxplanets,errcount=nerr)
    do i=1,nplanets
       call read_inopt(mplanet(i),'mplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
       call read_inopt(rplanet(i),'rplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
       call read_inopt(inclplan(i),'inclplanet'//trim(planets(i)),db,min=0.,max=180.,errcount=nerr)
       call read_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
    enddo
 endif
 !--timestepping
 !  following two are optional: not an error if not present
 call read_inopt(norbits,'norbits',db,err=ierr)
 call read_inopt(deltat,'deltat',db,err=ierr)
#ifdef MCFOST
 !--mcfost
 call read_inopt(use_mcfost,'use_mcfost',db,err=ierr)
 if (ierr /= 0) use_mcfost = .false. ! no mcfost by default
#endif

 call close_db(db)
 ierr = nerr
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
 endif

end subroutine read_setupfile

!------------------------------------------------------------------------
!
! subroutine for reading old options for backwards compatibility
!
!------------------------------------------------------------------------
subroutine read_obsolete_setup_options(db)
 use infile_utils, only:inopts,read_inopt
 type(inopts), allocatable, intent(inout) :: db(:)
 real, parameter :: tol = 0.01
 integer :: tmp_i,ierr
 logical :: tmp_l
 real    :: tmp_r

 ierr = 0
 tmp_r = 0.
 tmp_i = -1
 tmp_l = .false.

 call read_inopt(np,'npart',db,err=ierr)
 call read_inopt(tmp_i,'np_dust',db,err=ierr)
 dust_method = 2
 if (ierr /= 0) dust_method = 1
 call read_inopt(tmp_r,'udist',db,err=ierr)
 dist_unit = 'au'
 if (ierr==0 .and. abs(tmp_r - au) < au*tol) dist_unit = 'au'
 call read_inopt(tmp_r,'umass',db,err=ierr)
 mass_unit = 'solarm'
 if (ierr==0 .and. abs(tmp_r - solarm) < solarm*tol) mass_unit = 'solarm'
 call read_inopt(tmp_i,'icentralforce',db,err=ierr)
 icentral = 1
 nsinks = 1
 if (ierr==0 .and. tmp_i==1) then
    icentral = 0
    ipotential = 1
 endif
 call read_inopt(tmp_i,'binary_set',db,err=ierr)
 if (ierr==0 .and. tmp_i==1) then
    nsinks = 2
    ibinary = 0
    iuse_disc(1) = .true.
    call read_inopt(disc_m(1),'disc_m',db,err=ierr)
    call read_inopt(pindex(1),'pindex',db,err=ierr)
    call read_inopt(qindex(1),'qindex',db,err=ierr)
    call read_inopt(R_in(1),'R_in',db,err=ierr)
    call read_inopt(R_out(1),'R_out',db,err=ierr)
    call read_inopt(pindex_dust(1),'pindex_dust',db,err=ierr)
    call read_inopt(R_indust(1),'R_indust',db,err=ierr)
    call read_inopt(R_outdust(1),'R_outdust',db,err=ierr)
 endif
 call read_inopt(tmp_l,'use_binarydisc',db,err=ierr)
 if (ierr==0) then
    nsinks = 2
    ibinary = 0
    iuse_disc(1) = tmp_l
    call read_inopt(iuse_disc(2),'use_primarydisc',db,err=ierr)
    call read_inopt(iuse_disc(3),'use_secondardisc',db,err=ierr)
    do i=1,3
       if (iuse_disc(i)) then
          call read_inopt(disc_m(i),'discmass'//trim(disctype(i)),db,err=ierr)
          call read_inopt(pindex(i),'p_index'//trim(disctype(i)),db,err=ierr)
          call read_inopt(qindex(i),'q_index'//trim(disctype(i)),db,err=ierr)
          call read_inopt(H_R(i),'HoverR'//trim(disctype(i)),db,err=ierr)
       endif
    enddo
 endif
 call read_inopt(m1,'star_m',db,err=ierr)
 call read_inopt(m1,'object_mass',db,err=ierr)
 call read_inopt(tmp_r,'massratio',db,err=ierr)
 if (ierr==0) m2 = tmp_r*m1
 call read_inopt(accr1,'accradius',db,err=ierr)
 call read_inopt(accr2,'accretion_radius2',db,err=ierr)
 call read_inopt(binary_a,'binary_separation',db,err=ierr)
 call read_inopt(binary_a,'semimajoraxis',db,err=ierr)
 call read_inopt(binary_e,'ecc',db,err=ierr)
 call read_inopt(binary_e,'eccentricity',db,err=ierr)
 call read_inopt(tmp_i,'mass_set',db,err=ierr)
 mass_set = 0
 if (ierr==0 .and. tmp_i==1) mass_set(1) = 2
 call read_inopt(disc_m(1),'discmass',db,err=ierr)
 call read_inopt(disc_m(1),'disc_mass',db,err=ierr)
 call read_inopt(incl(1),'xinc',db,err=ierr)
 call read_inopt(incl(1),'inclination',db,err=ierr)
 call read_inopt(pindex(1),'p_index',db,err=ierr)
 call read_inopt(qindex(1),'q_index',db,err=ierr)
 call read_inopt(H_R(1),'HoverR',db,err=ierr)
 call read_inopt(tmp_i,'iprofilegas',db,err=ierr)
 itapergas = .false.
 if (ierr==0 .and. tmp_i==1) itapergas = .true.
 call read_inopt(sig_norm(1),'sigma_naught',db,err=ierr)
 call read_inopt(dust_to_gas,'dust_to_gas',db,err=ierr)
 call read_inopt(tmp_i,'profile_set_dust',db,err=ierr)
 if (ierr /= 0) then
    iprofile_dust = 1
    pindex_dust = pindex
 endif
 call read_inopt(tmp_r,'graindensinp',db,err=ierr)
 if (ierr /= 0) graindensinp = 3.
 nplanets = 0
 call read_inopt(nplanets,'nplanets',db,err=ierr)
 if (nplanets > 0) setplanets = 1
 call read_inopt(R_in(1),'R_in',db,err=ierr)
 R_ref = R_in(1)

end subroutine read_obsolete_setup_options

!------------------------------------------------------------------------
!
! calculates dust-to-gas ratio at a particular value of R
!
!------------------------------------------------------------------------
subroutine get_dust_to_gas(dust_to_gas,R,sigmaprofilegas,sigmaprofiledust, &
                                 sig_norm,sig_normdust,pindex,pindex_dust, &
                                 R_in,R_indust,R_ref,zi,hgas,hdust,R_c,R_c_dust)
 real,           intent(in)  :: R,pindex,pindex_dust,sig_norm,sig_normdust
 real,           intent(in)  :: R_in,R_indust,R_ref,zi,hgas,hdust
 real,           intent(in)  :: R_c,R_c_dust
 integer,        intent(in)  :: sigmaprofilegas,sigmaprofiledust
 real,           intent(out) :: dust_to_gas
 real :: sigma_gas,sigma_dust

 sigma_gas   = sig_norm     * scaled_sigma(R,sigmaprofilegas,pindex,R_ref,R_in,R_c)
 sigma_dust  = sig_normdust * scaled_sigma(R,sigmaprofiledust,pindex_dust,R_ref,R_indust,R_c_dust)
 dust_to_gas = (sigma_dust/sigma_gas) * (hgas/hdust) * exp(-0.5d0*(((zi/hdust)**2.)*(1.d0-(hdust/hgas)**2.)))

end subroutine get_dust_to_gas

!----------------------------------------------------------------
!
! height scale as a function of radius
!
!----------------------------------------------------------------
real function get_H(h0,qindex,r)
 real, intent(in) :: h0,qindex,r

 get_H = h0*(r**(-qindex+1.5))

end function get_H
!----------------------------------------------------------------
!
! spherical density profile as a function of radius
!
!----------------------------------------------------------------
real function atm_dens(r)
 use eos, only:gamma
 real, intent(in) :: r

 select case(atm_type)
 case(1)
    atm_dens = r**(-3)
 case(2)
    atm_dens = r**(-1./(gamma - 1.))
 case default
    !atm_dens = exp(-(r-r_planet)/scaleheight)
    stop 'atmosphere not yet implemented...stopping!'
 end select

end function atm_dens

!----------------------------------------------------------------
!
! function to return the sound speed given the radius
!
!----------------------------------------------------------------
pure real function cs_func(cs0,r,q_index)
 real, intent(in) :: cs0,r,q_index

 cs_func = cs0*r**(-q_index)

end function cs_func

!------------------------------------------------------------------------
!
! Set sphere of particles around one of the binary masses (i.e. planet)
!
!------------------------------------------------------------------------
subroutine set_planet_atm(id,xyzh,vxyzu,npartoftype,maxvxyzu,itype,a0,R_in, &
                          HoverR,Mstar,q_index,gamma,Ratm_in,Ratm_out,r_surface, &
                          npart,npart_planet_atm,npart_disc,hfact)
 use extern_binary, only:binarymassr,ramp
 use io,            only:master
 use part,          only:set_particle_type
 use spherical,     only:set_sphere
 use physcon,       only:pi
 integer, intent(in)    :: id,maxvxyzu,itype
 integer, intent(inout) :: npart,npart_planet_atm,npart_disc
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: Ratm_in,Ratm_out
 real,    intent(in)    :: a0,R_in,HoverR,Mstar,q_index,gamma,r_surface,hfact
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, parameter :: igas = 1
 integer(kind=8)    :: nptot
 integer            :: i,nx
 real :: xyz_orig(3)
 real :: a_orbit
 real :: psep,vol_sphere
 real :: cs0,cs
 !
 ! place particles in sphere
 !
 Ratm_in   = Ratm_in*r_surface
 Ratm_out  = Ratm_out*r_surface

 if (ramp) then
    xyz_orig(:) = (/a0,0.,0./)
 else
    a_orbit = a0 - binarymassr
    xyz_orig(:) = (/a_orbit,0.,0./)
 endif

 vol_sphere  = 4./3.*pi*Ratm_out**3
 nx          = int(npart_planet_atm**(1./3.))
 psep        = vol_sphere**(1./3.)/real(nx)
 nptot       = npart

 call set_sphere('closepacked',id,master,Ratm_in,Ratm_out,psep,hfact,npart,xyzh, &
                 rhofunc=atm_dens,nptot=nptot, &
                 np_requested=npart_planet_atm,xyz_origin=xyz_orig)

 npart_planet_atm = npart-npart_disc
 npartoftype(1) = npart
 do i=npart_disc+1,npart
    !--set the particle type for the atmosphere particles
    call set_particle_type(i,1)
    !-----------------------------------------
    !  Set thermal energy
    !  utherm generally should not be stored
    !  for an isothermal equation of state
    if (maxvxyzu >= 4) then
       if (itype==igas) then
          cs0 = sqrt(1.0/R_in)*(HoverR)*sqrt(Mstar)*(1.0/R_in)**(-q_index)
          cs = cs_func(cs0,a0,q_index)
          if (gamma > 1.) then
             vxyzu(4,i) = cs**2/(gamma - 1.)/gamma
          else
             vxyzu(4,i) = 1.5*cs**2
          endif
       else
          vxyzu(4,i) = 0.
       endif
    endif
 enddo

end subroutine set_planet_atm

!------------------------------------------------------------------------
!
! Convert to a corotating frame around one of the binary masses
!
!------------------------------------------------------------------------
subroutine make_corotate(xyzh,vxyzu,a0,Mstar,npart,npart_disc)
 use extern_corotate, only:omega_corotate
 integer, intent(inout) :: npart,npart_disc
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(in)    :: xyzh(:,:),a0,Mstar
 integer :: i
 real    :: phipart,r
 real    :: v_0(3),vmag,omega0

 !
 !--Change to corotating frame
 !
 ! Calculate velocity at planet
 vmag   = sqrt(Mstar/a0)
 omega0 = sqrt(Mstar/a0**3)

 ! v_phi = v_y at y=0
 ! Obtain the true v_phi at any point (r,phi) via rotation in z axis

 v_0 = (/0.0,vmag,0.0/)

 print *, 'Transforming to corotating frame: angular velocity ', omega0

 do i=1,npart_disc
    r          = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
    phipart    = atan2(xyzh(2,i),xyzh(1,i))
    vxyzu(1,i) = vxyzu(1,i) - r*(-omega0)*sin(phipart)
    vxyzu(2,i) = vxyzu(2,i) + r*(-omega0)*cos(phipart)
 enddo
 vxyzu(1:3,npart_disc+1:npart) = 0.

 omega_corotate = omega0

end subroutine make_corotate

end module setup
