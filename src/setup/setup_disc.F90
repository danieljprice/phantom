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
!   This module sets up an accretion disc with options for dust, embedded
!   planets, binaries, and flybys.
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    H_R               -- H/R at R=R_ref
!    H_R_dust          -- H/R at R=R_ref
!    R_c               -- characteristic radius of the exponential taper
!    R_c_dust          -- characteristic radius of the exponential taper
!    R_in              -- inner radius
!    R_indust          -- inner radius
!    R_out             -- outer radius
!    R_outdust         -- outer radius
!    R_ref             -- reference radius
!    accr1             -- central star accretion radius
!    accr2             -- perturber accretion radius
!    alphaSS           -- desired alphaSS
!    bhspin            -- black hole spin
!    bhspinangle       -- black hole spin angle
!    binary_O          -- Omega, PA of ascending node (deg)
!    binary_a          -- binary semi-major axis
!    binary_e          -- binary eccentricity
!    binary_f          -- f, initial true anomaly (deg,180=apastron)
!    binary_i          -- i, inclination (deg)
!    binary_w          -- w, argument of periapsis (deg)
!    deltat            -- output interval as fraction of total time
!    disc_m            -- disc mass
!    dist_unit         -- distance unit (e.g. au,pc,kpc,0.1pc)
!    dust_to_gas_ratio -- dust to gas ratio
!    flyby_a           -- flyby periastron distance
!    flyby_d           -- initial distance of flyby (in units of periastron distance)
!    flyby_r           -- roll angle of flyby
!    graindensinp      -- intrinsic grain density (in g/cm^3)
!    grainsizeinp      -- grain size (in cm)
!    ibinary           -- closed binary or flyby (0=binary,1=flyby)
!    ipotential        -- potential (1=central point mass,2=binary potential,3=spinning black hole)
!    ismoothdust       -- smooth the inner disc profile
!    ismoothgas        -- smooth the inner disc profile
!    itaperdust        -- exponentially taper the outer disc profile
!    itapergas         -- exponentially taper the outer disc profile
!    m1                -- central star mass
!    m2                -- perturber mass
!    mass_set          -- how to set gas density profile (0=disc mass,1=surface density)
!    mass_unit         -- mass unit (e.g. solarm,jupiterm,earthm)
!    norbits           -- maximum number of binary orbits
!    np                -- number of gas particles
!    np_dust           -- number of dust particles
!    nplanets          -- number of planets
!    nsinks            -- number of sinks
!    pindex            -- p index
!    pindex_dust       -- p index
!    profile_set_dust  -- how to set dust density profile (0=equal to gas, 1=custom)
!    qindex            -- q index
!    qindex_dust       -- q index
!    setplanets        -- add planets? (0=no,1=yes)
!    sigma_ref         -- sigma at R=R_ref
!    xinc              -- inclination angle
!
!  DEPENDENCIES: centreofmass, dim, dust, eos, extern_binary,
!    extern_lensethirring, externalforces, infile_utils, io, options, part,
!    physcon, prompting, setbinary, setdisc, setflyby, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 use dim, only:maxp,use_dust,use_dustfrac,maxalpha
 implicit none
 public  :: setpart
 real    :: R_in,R_out,R_ref,xinc,disc_m,pindex,qindex,pindex_dust,qindex_dust
 real    :: H_R,H_R_dust,sig_norm,sig_normdust,grainsizeinp,graindensinp
 real    :: disc_mdust,dust_to_gas_ratio,R_outdust,R_indust,R_c,R_c_dust
 real    :: m1,m2,mcentral,bhspin,bhspinangle,binary_a,binary_e,binary_i,binary_O,binary_w,binary_f
 real    :: accr1,accr2,alphaSS,deltat,flyby_a,flyby_d,flyby_r
 integer :: i,nplanets,np,np_dust,icentral,ipotential,nsinks,ibinary,setplanets,norbits
 integer :: itype,npart_inter,mass_set,profile_set_dust,iprofilegas,iprofiledust
 integer :: sigmaprofilegas,sigmaprofiledust
 logical :: ismoothgas,ismoothdust,itapergas,itaperdust
 integer, parameter :: maxplanets = 9
 real    :: mplanet(maxplanets),rplanet(maxplanets),accrplanet(maxplanets),inclplan(maxplanets)
 character(len=*), dimension(maxplanets), parameter :: planets = &
  (/'1','2','3','4','5','6','7','8','9' /)

 private
 character(len=20) :: dist_unit,mass_unit

contains

!----------------------------------------------------------------
!
! This subroutine sets up an accretion disc
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use centreofmass,         only:reset_centreofmass
 use dust,                 only:set_dustfrac,grainsizecgs,graindenscgs
 use eos,                  only:isink,qfacdisc
 use extern_binary,        only:accradius1,accradius2,binarymassr
 use externalforces,       only:mass1,accradius1,iext_star,iext_binary,iext_lensethirring
 use extern_lensethirring, only:blackhole_spin,blackhole_spin_angle
 use io,                   only:master,fatal,error
 use kernel,               only:hfact_default
 use options,              only:iexternalforce,ieos,alpha,icooling
 use part,                 only:nptmass,xyzmh_ptmass,maxvxyzu,vxyz_ptmass,ihacc,ihsoft,igas,idust,dustfrac
 use physcon,              only:au,solarm,jupiterm,pi,years
 use prompting,            only:prompt
 use setbinary,            only:set_binary
 use setdisc,              only:set_disc,scaled_discmass
 use setflyby,             only:set_flyby,get_T_flyby
 use timestep,             only:tmax,dtmax
 use units,                only:set_units,select_unit,umass,udist,utime
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=100)               :: filename
 logical :: iexist,questplanets,seq_exists
 real    :: phi,vphi,sinphi,cosphi,omega,r2,disc_m_within_r,period_longest
 real    :: idust_to_gas_ratio,ri,period
 real    :: sini,cosi,polyk_dust
 integer :: ierr,j

 print "(/,65('-'),2(/,a),/,65('-'),/)",&
   ' Welcome to the New Disc Setup'

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
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    !
    !--resolution
    !
    np = 1e6
    if (use_dust .and. .not. use_dustfrac) then
       np_dust = np/10
    else
       np_dust = 0
    endif
    !
    !--units
    !
    dist_unit = 'au'
    mass_unit = 'solarm'
    !
    !--set defaults for central object(s)
    !
    icentral = 1
    call prompt('Do you want to use sink particles or an external potential? (0=potential,1=sinks)',icentral,0,1)
    select case (icentral)
    case (0)
       !--external potential
       !--todo: check that potentials are implemented correctly; add more potentials
       ipotential = 1
       call prompt('Which potential? (1=central point mass,2=binary potential,3=spinning black hole)',ipotential,1,3)
       select case (ipotential)
       case (1)
          !--point mass
          iexternalforce = iext_star
          m1       = 1.
          accr1    = 1.
       case (2)
          !--fixed binary
          iexternalforce = iext_binary
          m1       = 1.
          m2       = 1.
          accr1    = 1.
          accr2    = 1.
       case (3)
          !--spinning black hole (Lense-Thirring)
          iexternalforce = iext_lensethirring
          m1          = 1.
          accr1       = 1.
          bhspin      = 1.
          bhspinangle = 0.
       end select
    case (1)
       !--sink particle(s)
       nsinks = 1
       call prompt('How many sinks? (1=single,2=binary)',nsinks,1,2)
       select case (nsinks)
       case (1)
          !--single star
          m1       = 1.
          accr1    = 1.
       case (2)
          ibinary = 0
          call prompt('Closed binary or flyby? (0=binary,1=flyby)',ibinary,0,1)
          select case (ibinary)
          case (0)
             !--binary
             m1       = 1.
             m2       = 1.
             binary_a = 1.
             binary_e = 0.
             binary_i = 0.
             binary_O = 0.
             binary_w = 0.
             binary_f = 180.
             accr1    = 0.25*binary_a
             accr2    = 0.25*binary_a
          case (1)
             !--flyby
             m1       = 1.
             m2       = 1.
             accr1    = 1.
             accr2    = 1.
             flyby_a  = 200.
             flyby_d  = 10.
             flyby_r  = 0.
          end select
       end select
    end select
    !
    !--set gas disc options
    !
    R_in            = 1.
    R_out           = 150.
    R_ref           = R_in
    R_c             = R_out
    mass_set        = 0
    itapergas       = .false.
    !--todo: add circumprimary & circumsecondary discs
    call prompt('How do you want to set the gas disc mass? (0=disc mass,1=surface density)',mass_set,0,1)
    call prompt('Do you want to exponentially taper the outer gas disc profile?',itapergas)
    select case (mass_set)
    case (0)
       disc_m = 0.05
    case (1)
       sig_norm = 1.0E-02
    end select
    pindex  = 1.
    qindex  = 0.25
    alphaSS = 0.005
    xinc    = 0.            !--todo: implement general rotations around a position angle
    H_R     = 0.05
    !
    !--set dust disc options
    !
    if (use_dust) then
       profile_set_dust  = 0
       dust_to_gas_ratio = 0.01
       R_indust          = R_in
       R_outdust         = R_out
       pindex_dust       = pindex
       qindex_dust       = qindex
       H_R_dust          = H_R
       itaperdust        = itapergas
       grainsizeinp      = 0.1
       graindensinp      = 3.
       R_c_dust          = R_c
       call prompt('How do you want to set the dust density profile? (0=equal to gas,1=custom)',profile_set_dust,0,1)
       if (profile_set_dust == 1) then
          call prompt('Do you want to exponentially taper the outer dust disc profile?',itaperdust)
       endif
       call prompt('Enter dust to gas ratio',dust_to_gas_ratio,0.)
       call prompt('Enter grain size in cm',grainsizeinp,0.)
    endif
    !
    !--add planets
    !
    questplanets  = .false.
    setplanets    = 0
    nplanets      = 0
    mplanet(:)    = 1.
    accrplanet(:) = 0.25
    do i=1,maxplanets
       rplanet(i) = 10.*i
       inclplan(i) = 0.
    enddo
    call prompt('Do you want to add planets?',questplanets)
    if (questplanets) then
       setplanets = 1
       call prompt('Enter the number of planets',nplanets,0,maxplanets)
    endif
    !
    !--determine simulation time
    !
    deltat  = 0.1
    norbits = 100
    if (nplanets > 0) then
       call prompt('Enter time between dumps as fraction of outer planet period',deltat,0.)
       call prompt('Enter number of orbits to simulate',norbits,0)
    else if (icentral==1 .and. nsinks==2) then
       if (ibinary==0) then
          call prompt('Enter time between dumps as fraction of binary period',deltat,0.)
          call prompt('Enter number of orbits to simulate',norbits,0)
       else if (ibinary==1) then
          deltat  = 0.01
          call prompt('Enter time between dumps as fraction of flyby time',deltat,0.)
       endif
    endif

    !
    !--write default input file
    !
    call write_setupfile(filename)

    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'

    stop
 else
    stop
 endif

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
 !--equation of state
 !
 if (maxvxyzu==3 .and. qindex > 0.) then
    !--use locally isothermal equation of state
    ieos = 3
    gamma = 1.0
 else
    ieos = 2
    gamma = 5./3.
    icooling = 1
 endif
#ifdef MCFOST
 icooling = 0
#endif

 !--sanity check on ieos = 6
 if (ieos==6 .and. isink==0) call fatal('setup_disc','something''s gone wrong with ieos & isink...')

 !
 !--surface density profile
 !
 ismoothgas      = .true.       !--smoothed at inner edge by default
 iprofilegas     = 0
 sigmaprofilegas = 0
 if (itapergas) then
    iprofilegas     = 1
    sigmaprofilegas = 1
 endif
 if (ismoothgas)                 sigmaprofilegas = 2
 if (itapergas .and. ismoothgas) sigmaprofilegas = 3
 if (use_dust) then
    ismoothdust      = .true.   !--smoothed at inner edge by default
    iprofiledust     = iprofilegas
    sigmaprofiledust = sigmaprofilegas
    if (itaperdust) then
       iprofiledust     = 1
       sigmaprofiledust = 1
    endif
    if (ismoothdust)                  sigmaprofiledust = 2
    if (itaperdust .and. ismoothdust) sigmaprofiledust = 3
 endif

 !
 !--set sink particle(s) or potential
 !
 select case (icentral)
 case (0)
    select case (ipotential)
    case (1)
       print "(a)", 'Central object represented by external force with accretion boundary'
       print*, ' Object mass:      ', m1
       print*, ' Accretion Radius: ', accr1
       mass1      = m1
       accradius1 = accr1
       mcentral   = m1
    case (2)
       print "(a)", 'Central binary represented by external force with accretion boundary'
       print*, ' Primary mass:       ', m1
       print*, ' Binary mass ratio:  ', m2/m1
       print*, ' Accretion Radius 1: ', accr1
       print*, ' Accretion Radius 2: ', accr2
       mass1       = m1
       binarymassr = m2/m1
       accradius1  = accr1
       accradius2  = accr2
       mcentral    = m1 + m2
    case (3)
       print "(a)", 'Central black hole represented by external force with accretion boundary'
       print*, ' Black hole mass:        ', m1
       print*, ' Accretion Radius:       ', accr1
       print*, ' Black hole spin:        ', bhspin
       print*, ' Black hole spin angle:  ', bhspinangle
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
       print*,'Disc around a single star '
       print "(a)", 'Central object represented by a sink at the system origin'
       isink                        = 1
       nptmass                      = 1
       xyzmh_ptmass(:,:)            = 0.
       xyzmh_ptmass(1:3,nptmass)    = 0.
       xyzmh_ptmass(4,nptmass)      = m1
       xyzmh_ptmass(ihacc,nptmass)  = accr1
       xyzmh_ptmass(ihsoft,nptmass) = accr1
       vxyz_ptmass                  = 0.
       mcentral                     = m1
    case (2)
       select case (ibinary)
       case (0)
          !--binary
          nptmass  = 0
          print*,'Setup the binary' !--the stars' barycentre is put on the origin
          call set_binary(m1,massratio=m2/m1,semimajoraxis=binary_a,eccentricity=binary_e, &
                          posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i,&
                          f=binary_f,accretion_radius1=accr1,accretion_radius2=accr2,&
                          xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass)
          mcentral = m1 + m2
       case (1)
          !--flyby
          print*,'Disc around a single star with flyby'
          print "(a)", 'Central object represented by a sink at the system origin'
          call set_flyby(mprimary=m1,massratio=m2/m1,dma=flyby_a,n0=flyby_d,roll=flyby_r, &
                         accretion_radius1=accr1,accretion_radius2=accr2, &
                         xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass)
          mcentral = m1
       end select
    end select
 end select

 !
 !--setup disc(s)
 !
 xinc = xinc*(pi/180.0)
 if (maxalpha==0) alpha = alphaSS
 time  = 0.
 hfact = hfact_default
 npartoftype(:)    = 0
 npartoftype(igas) = np
 !--compute the disc mass from sig_norm if required
 if (mass_set==1) disc_m = sig_norm*scaled_discmass(sigmaprofilegas,pindex,R_in,R_out,R_ref,R_c)
 if (use_dust) then
    grainsizecgs = grainsizeinp
    graindenscgs = graindensinp
    npartoftype(idust) = np_dust
    disc_mdust = disc_m*dust_to_gas_ratio
 endif
 npart = sum(npartoftype)
 npart_inter = 1
 itype       = igas
 if (use_dust .and. use_dustfrac) then
    !--gas and dust disc
    call set_disc(id,master        = master,             &
                  mixture          = .true.,             &
                  npart            = npartoftype(itype), &
                  npart_start      = npart_inter,        &
                  rref             = R_ref,              &
                  rmin             = R_in,               &
                  rmax             = R_out,              &
                  rmindust         = R_indust,           &
                  rmaxdust         = R_outdust,          &
                  indexprofile     = iprofilegas,        &
                  indexprofiledust = iprofiledust,       &
                  rc               = R_c,                &
                  rcdust           = R_c_dust,           &
                  p_index          = pindex,             &
                  p_indexdust      = pindex_dust,        &
                  q_index          = qindex,             &
                  HoverR           = H_R,                &
                  disc_mass        = disc_m,             &
                  disc_massdust    = disc_mdust,         &
                  star_mass        = mcentral,           &
                  gamma            = gamma,              &
                  particle_mass    = massoftype(itype),  &
                  hfact            = hfact,              &
                  xyzh             = xyzh,               &
                  vxyzu            = vxyzu,              &
                  polyk            = polyk,              &
                  alpha            = alpha,              &
                  ismooth          = ismoothgas,         &
                  inclination      = xinc,               &
                  twist            = .false.,            &
                  prefix           = fileprefix)
    !--set dustfrac
    sig_norm     = disc_m     / scaled_discmass(sigmaprofilegas,pindex,R_in,R_out,R_ref,R_c)
    sig_normdust = disc_mdust / scaled_discmass(sigmaprofiledust,pindex_dust,R_indust,R_outdust,R_ref,R_c_dust)
    do i=1,npart
       Ri = sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       if (Ri<R_indust .or. Ri>R_outdust) then
          idust_to_gas_ratio = 0.
       else
          call get_dust_to_gas_ratio(idust_to_gas_ratio,Ri,sigmaprofilegas,sigmaprofiledust,&
                                     sig_norm,sig_normdust,pindex,pindex_dust,&
                                     R_in,R_ref,R_c,R_indust,R_c_dust)
       endif
       call set_dustfrac(idust_to_gas_ratio,dustfrac(i))
    enddo
 else
    !--gas disc
    call set_disc(id,master       = master,             &
                  npart           = npartoftype(itype), &
                  npart_start     = npart_inter,        &
                  particle_type   = itype,              &
                  rref            = R_ref,              &
                  rmin            = R_in,               &
                  rmax            = R_out,              &
                  indexprofile    = iprofilegas,        &
                  rc              = R_c,                &
                  p_index         = pindex,             &
                  q_index         = qindex,             &
                  HoverR          = H_R,                &
                  disc_mass       = disc_m,             &
                  star_mass       = mcentral,           &
                  gamma           = gamma,              &
                  particle_mass   = massoftype(itype),  &
                  hfact           = hfact,              &
                  xyzh            = xyzh,               &
                  vxyzu           = vxyzu,              &
                  polyk           = polyk,              &
                  alpha           = alpha,              &
                  ismooth         = ismoothgas,         &
                  inclination     = xinc,               &
                  twist           = .false.,            &
                  prefix          = fileprefix)
    if (use_dust) then
       itype       = idust
       npart_inter = npartoftype(igas) + 1
       !--dust disc
       call set_disc(id,master     = master,             &
                     npart         = npartoftype(itype), &
                     npart_start   = npart_inter,        &
                     particle_type = itype,              &
                     rref          = R_ref,              &
                     rmin          = R_indust,           &
                     rmax          = R_outdust,          &
                     indexprofile  = iprofiledust,       &
                     rc            = R_c_dust,           &
                     p_index       = pindex_dust,        &
                     q_index       = qindex_dust,        &
                     HoverR        = H_R_dust,           &
                     disc_mass     = disc_mdust,         &
                     star_mass     = mcentral,           &
                     gamma         = gamma,              &
                     particle_mass = massoftype(itype),  &
                     hfact         = hfact,              &
                     xyzh          = xyzh,               &
                     vxyzu         = vxyzu,              &
                     polyk         = polyk_dust,         &
                     ismooth       = ismoothdust,        &
                     inclination   = xinc,               &
                     twist         = .false.,            &
                     prefix        = fileprefix)
       !--reset qfacdisc to gas disc value
       qfacdisc = qindex
    endif
 endif

 !
 !--planets
 !
 if (setplanets==1) then
    print "(a,i2,a)",' --------- added ',nplanets,' planets ------------'
    period_longest = 0.
    do i=1,nplanets
       nptmass = nptmass + 1
       phi = 0.*pi/180.
       cosphi = cos(phi)
       sinphi = sin(phi)
       disc_m_within_r = 0.
       do j=1,npart
          r2 = xyzh(1,j)**2 + xyzh(2,j)**2 + xyzh(3,j)**2
          if (r2 < rplanet(i)**2) then
             if (.not.use_dust .or. (use_dust .and. use_dustfrac)) then
                disc_m_within_r = disc_m_within_r + massoftype(igas)
             else
                if (j<=np) then
                   disc_m_within_r = disc_m_within_r + massoftype(igas)
                else
                   disc_m_within_r = disc_m_within_r + massoftype(idust)
                endif
             endif
          endif
       enddo
       if (nplanets>1) then
          do j=1,nplanets
             if (rplanet(j)<rplanet(i)) disc_m_within_r = disc_m_within_r + mplanet(j)*jupiterm/solarm
          enddo
       endif
       xyzmh_ptmass(1:3,nptmass)    = (/rplanet(i)*cosphi,rplanet(i)*sinphi,0./)
       xyzmh_ptmass(4,nptmass)      = mplanet(i)*jupiterm/umass
       xyzmh_ptmass(ihacc,nptmass)  = accrplanet(i)
       xyzmh_ptmass(ihsoft,nptmass) = accrplanet(i)
       vphi                         = sqrt((mcentral + disc_m_within_r)/rplanet(i))
       vxyz_ptmass(1:3,nptmass)     = (/-vphi*sinphi,vphi*cosphi,0./)
       !--rotate positions and velocities
       if (inclplan(i)  /=  0.) then
          cosi=cos(inclplan(i)*pi/180.)
          sini=sin(inclplan(i)*pi/180.)
          call rotate(xyzmh_ptmass(1:3,nptmass),cosi,sini)
          call rotate(vxyz_ptmass(1:3,nptmass),cosi,sini)
       endif
       print "(a,i2,a)",       ' planet ',i,':'
       print "(a,g10.3,a)",    ' radius: ',rplanet(i)*udist/au,' AU'
       print "(a,g10.3,a,2pf7.3,a)",    ' M(<R) : ',(disc_m_within_r + mcentral)*umass/solarm, &
             ' MSun, disc mass correction is ',disc_m_within_r/mcentral,'%'
       print "(a,2(g10.3,a))", ' mass  : ',mplanet(i),' MJup, or ',mplanet(i)*jupiterm/solarm,' MSun'
       print "(a,2(g10.3,a))", ' period: ',2.*pi*rplanet(i)/vphi*utime/years,' years or ',2*pi*rplanet(i)/vphi,' in code units'
       omega = vphi/rplanet(i)
       period_longest = max(period_longest, 2.*pi/omega)
       print "(a,g10.3,a)",   ' resonances: 3:1: ',(sqrt(mcentral)/(3.*omega))**(2./3.),' AU'
       print "(a,g10.3,a)",   '             4:1: ',(sqrt(mcentral)/(4.*omega))**(2./3.),' AU'
       print "(a,g10.3,a)",   '             5:1: ',(sqrt(mcentral)/(5.*omega))**(2./3.),' AU'
       print "(a,g10.3,a)",   '             9:1: ',(sqrt(mcentral)/(9.*omega))**(2./3.),' AU'
    enddo
    print "(1x,45('-'))"
    period = period_longest
 endif

 !
 !--reset centre of mass to the origin
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 !
 !--set tmax and dtmax
 !
 if (icentral==1 .and. nsinks==2 .and. ibinary==0) then
    !--closed binary
    period = sqrt(4.*pi**2*binary_a**3/mcentral)
 elseif (icentral==1 .and. nsinks==2 .and. ibinary==1) then
    !--flyby
    period = get_T_flyby(m1,m2,flyby_a,flyby_d)
    norbits = 1
 elseif (setplanets==1) then
    !--outer planet set above
 else
    !--outer disc
    period = sqrt(4.*pi**2*R_out**3/mcentral)
 endif
 if (period > 0.) then
    if (deltat > 0.)  dtmax = deltat*period
    if (norbits >= 0) tmax  = norbits*period
 endif

 return
end subroutine setpart

!------------------------------------------------------------------------
!
! write setup file
!
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for disc setup routine'
 !--resolution
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','number of gas particles',iunit)
 if (use_dust .and. .not.use_dustfrac) then
    call write_inopt(np_dust,'np_dust','number of dust particles',iunit)
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
    call write_inopt(ipotential,'ipotential','potential (1=central point mass,2=binary potential,3=spinning black hole)',iunit)
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
    case (3)
       !--spinning black hole (Lense-Thirring)
       call write_inopt(m1,'m1','black hole mass',iunit)
       call write_inopt(accr1,'accr1','black hole accretion radius',iunit)
       call write_inopt(bhspin,'bhspin','black hole spin',iunit)
       call write_inopt(bhspinangle,'bhspinangle','black hole spin angle',iunit)
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
       call write_inopt(ibinary,'ibinary','closed binary or flyby (0=binary,1=flyby)',iunit)
       select case (ibinary)
       case (0)
          !--closed binary
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
          !--flyby
          !--central star
          write(iunit,"(/,a)") '# options for central star'
          call write_inopt(m1,'m1','central star mass',iunit)
          call write_inopt(accr1,'accr1','central star accretion radius',iunit)
          !--perturber
          write(iunit,"(/,a)") '# options for perturber'
          call write_inopt(m2,'m2','perturber mass',iunit)
          call write_inopt(accr2,'accr2','perturber accretion radius',iunit)
          call write_inopt(flyby_a,'flyby_a','flyby periastron distance',iunit)
          call write_inopt(flyby_d,'flyby_d','initial distance of flyby (in units of periastron distance)',iunit)
          call write_inopt(flyby_r,'flyby_r','roll angle of flyby',iunit)
       end select
    end select
 end select
 !--gas disc
 write(iunit,"(/,a)") '# options for gas accretion disc'
 call write_inopt(mass_set,'mass_set','how to set gas density profile (0=disc mass,1=surface density)',iunit)
 call write_inopt(itapergas,'itapergas','exponentially taper the outer disc profile',iunit)
 call write_inopt(R_in,'R_in','inner radius',iunit)
 call write_inopt(R_ref,'R_ref','reference radius',iunit)
 call write_inopt(R_out,'R_out','outer radius',iunit)
 if (itapergas) call write_inopt(R_c,'R_c','characteristic radius of the exponential taper',iunit)
 select case (mass_set)
 case (0)
    call write_inopt(disc_m,'disc_m','disc mass',iunit)
 case (1)
    call write_inopt(sig_norm,'sig_norm','sigma normalisation',iunit)
 end select
 call write_inopt(pindex,'pindex','p index',iunit)
 call write_inopt(qindex,'qindex','q index',iunit)
 if (maxalpha==0) call write_inopt(alphaSS,'alphaSS','desired alphaSS',iunit)
 call write_inopt(xinc,'xinc','inclination angle',iunit)
 call write_inopt(H_R,'H_R','H/R at R=R_ref',iunit)
 !--dust disc
 if (use_dust) then
    write(iunit,"(/,a)") '# options for dust accretion disc'
    call write_inopt(dust_to_gas_ratio,'dust_to_gas_ratio','dust to gas ratio',iunit)
    call write_inopt(profile_set_dust,'profile_set_dust','how to set dust density profile (0=equal to gas, 1=custom)',iunit)
    if (profile_set_dust==1) then
       call write_inopt(itaperdust,'itaperdust','exponentially taper the outer disc profile',iunit)
       call write_inopt(R_indust,'R_indust','inner radius',iunit)
       call write_inopt(R_outdust,'R_outdust','outer radius',iunit)
       if (iprofiledust==1) call write_inopt(R_c_dust,'R_c_dust','characteristic radius of the exponential taper',iunit)
       call write_inopt(pindex_dust,'pindex_dust','p index',iunit)
       if (.not. use_dustfrac) then
          call write_inopt(qindex_dust,'qindex_dust','q index',iunit)
          call write_inopt(H_R_dust,'H_R_dust','H/R at R=R_ref',iunit)
       endif
    endif
    call write_inopt(grainsizeinp,'grainsizeinp','grain size (in cm)',iunit)
    call write_inopt(graindensinp,'graindensinp','intrinsic grain density (in g/cm^3)',iunit)
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
       call write_inopt(inclplan(i),'inclplanet'//trim(planets(i)),'planet inclination [deg] with respect to xy plane',iunit)
       call write_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),'planet radius',iunit)
    enddo
 endif
 if (nplanets > 0) then
    write(iunit,"(/,a)") '# timestepping'
    call write_inopt(norbits,'norbits','maximum number of outer planet orbits',iunit)
    call write_inopt(deltat,'deltat','output interval as fraction of orbital period',iunit)
 else if (icentral==1 .and. nsinks==2) then
    if (ibinary==0) then
       write(iunit,"(/,a)") '# timestepping'
       call write_inopt(norbits,'norbits','maximum number of binary orbits',iunit)
       call write_inopt(deltat,'deltat','output interval as fraction of binary orbital period',iunit)
    else if (ibinary==1) then
       write(iunit,"(/,a)") '# timestepping'
       call write_inopt(deltat,'deltat','output interval as fraction of total time',iunit)
    endif
 else
    write(iunit,"(/,a)") '# timestepping'
    call write_inopt(norbits,'norbits','maximum number of orbits at outer disc',iunit)
    call write_inopt(deltat,'deltat','output interval as fraction of orbital period',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!
! read setup file
!
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use externalforces, only:iext_star,iext_binary,iext_lensethirring
 use infile_utils,   only:open_db_from_file,inopts,read_inopt,close_db
 use options,        only:iexternalforce
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)

 !--read old options for backwards compatibility
 ! call read_obsolete_setup_options(db)

 !--resolution
 call read_inopt(np,'np',db,min=0,max=maxp,errcount=nerr)
 if (use_dust .and. .not.use_dustfrac) then
    call read_inopt(np_dust,'np_dust',db,min=0,max=maxp-np,errcount=nerr)
 endif
 !--units
 call read_inopt(mass_unit,'mass_unit',db,errcount=nerr)
 call read_inopt(dist_unit,'dist_unit',db,errcount=nerr)
 !--central objects(s)/potential
 call read_inopt(icentral,'icentral',db,min=0,max=1, errcount=nerr)
 select case (icentral)
 case (0)
    !--external potential
    call read_inopt(ipotential,'ipotential',db,min=1,max=3, errcount=nerr)
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
    case (3)
       !--spinning black hole (Lense-Thirring)
       iexternalforce = iext_lensethirring
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
       call read_inopt(bhspin,'bhspin',db,min=0.,errcount=nerr)
       call read_inopt(bhspinangle,'bhspinangle',db,min=0.,errcount=nerr)
    end select
 case (1)
    !--sink particles
    call read_inopt(nsinks,'nsinks',db,min=1,max=2, errcount=nerr)
    select case (nsinks)
    case (1)
       !--single star
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
    case (2)
       !--two sinks
       call read_inopt(ibinary,'ibinary',db,min=0,max=1, errcount=nerr)
       select case (ibinary)
       case (0)
          !--binary
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
          !--flyby
          !--central star
          call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
          call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
          !--perturber
          call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
          call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)
          call read_inopt(flyby_a,'flyby_a',db,min=0.,errcount=nerr)
          call read_inopt(flyby_d,'flyby_d',db,min=0.,errcount=nerr)
          call read_inopt(flyby_r,'flyby_r',db,min=0.,errcount=nerr)
       end select
    end select
 end select
 !--gas disc
 call read_inopt(R_in,'R_in',db,min=0.,errcount=nerr)
 call read_inopt(R_out,'R_out',db,min=R_in,errcount=nerr)
 call read_inopt(R_ref,'R_ref',db,min=R_in,errcount=nerr)
 call read_inopt(itapergas,'itapergas',db,errcount=nerr)
 call read_inopt(mass_set,'mass_set',db,min=0,max=1,errcount=nerr)
 if (itapergas) then
    call read_inopt(R_c,'R_c',db,min=0.,errcount=nerr)
 endif
 select case (mass_set)
 case (0)
    call read_inopt(disc_m,'disc_m',db,min=0.,errcount=nerr)
 case (1)
    call read_inopt(sig_norm,'sig_norm',db,min=0.,errcount=nerr)
 end select
 call read_inopt(pindex,'pindex',db,errcount=nerr)
 call read_inopt(qindex,'qindex',db,errcount=nerr)
 if (maxalpha==0) call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
 call read_inopt(xinc,'xinc',db,min=0.,max=180.,errcount=nerr)
 call read_inopt(H_R,'H_R',db,min=0.,errcount=nerr)
 !--dust disc
 if (use_dust) then
    call read_inopt(dust_to_gas_ratio,'dust_to_gas_ratio',db,min=0.,errcount=nerr)
    call read_inopt(profile_set_dust,'profile_set_dust',db,min=0,max=1,errcount=nerr)
    select case (profile_set_dust)
    case (0)
       R_indust         = R_in
       R_outdust        = R_out
       pindex_dust      = pindex
       qindex_dust      = qindex
       H_R_dust         = H_R
       itaperdust       = itapergas
       R_c_dust         = R_c
    case (1)
       call read_inopt(R_indust,'R_indust',db,min=R_in,errcount=nerr)
       call read_inopt(R_outdust,'R_outdust',db,min=R_indust,max=R_out,errcount=nerr)
       call read_inopt(pindex_dust,'pindex_dust',db,errcount=nerr)
       call read_inopt(itaperdust,'itaperdust',db,errcount=nerr)
       if (itaperdust) then
          call read_inopt(R_c_dust,'R_c_dust',db,min=0.,errcount=nerr)
       endif
       if (.not. use_dustfrac) then
          call read_inopt(qindex_dust,'qindex_dust',db,err=ierr,errcount=nerr)
          if (ierr /= 0) qindex_dust = qindex
          call read_inopt(H_R_dust,'H_R_dust',db,min=0.,err=ierr,errcount=nerr)
          if (ierr /= 0) H_R_dust = H_R
       endif
    end select
    call read_inopt(grainsizeinp,'grainsizeinp',db,min=0.,errcount=nerr)
    call read_inopt(graindensinp,'graindensinp',db,min=0.,errcount=nerr)
 endif
 !--planets
 setplanets = 0
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
 !--following two are optional: not an error if not present
 norbits = 0
 deltat = 0.
 call read_inopt(norbits,'norbits',db,err=ierr)
 call read_inopt(deltat,'deltat',db,err=ierr)
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
!subroutine read_obsolete_setup_options(db)
!end subroutine read_obsolete_setup_options

!------------------------------------------------------------------------
!
! calculates dust-to-gas ratio at a particular value of R
!
!------------------------------------------------------------------------
subroutine get_dust_to_gas_ratio(dust_to_gas,R,sigmaprofilegas,sigmaprofiledust,&
                                 sig_norm,sig_normdust,pindex,pindex_dust,&
                                 R_in,R_ref,R_c,R_indust,R_c_dust)
 use setdisc, only:scaled_sigma
 real,           intent(in)  :: R,pindex,pindex_dust,sig_norm,sig_normdust
 real,           intent(in)  :: R_ref,R_in,R_indust
 real, optional, intent(in)  :: R_c,R_c_dust
 integer,        intent(in)  :: sigmaprofilegas,sigmaprofiledust
 real,           intent(out) :: dust_to_gas
 real :: sigma_gas,sigma_dust

 sigma_gas   = sig_norm     * scaled_sigma(R,sigmaprofilegas,pindex,R_ref,R_in,R_c)
 sigma_dust  = sig_normdust * scaled_sigma(R,sigmaprofiledust,pindex_dust,R_ref,R_indust,R_c_dust)
 dust_to_gas = sigma_dust / sigma_gas

end subroutine get_dust_to_gas_ratio

!------------------------------------------------------------------------
!
! rotate
!
!------------------------------------------------------------------------
pure subroutine rotate(xyz,cosi,sini)
 real, intent(inout) :: xyz(3)
 real, intent(in)    :: cosi,sini
 real :: xi,yi,zi

 xi = xyz(1)
 yi = xyz(2)
 zi = xyz(3)
 xyz(1) =  xi*cosi + zi*sini
 xyz(2) =  yi
 xyz(3) = -xi*sini + zi*cosi

end subroutine rotate

end module setup
