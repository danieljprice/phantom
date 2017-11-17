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
!    R_inann           -- inner annulus radius
!    R_indust          -- inner radius
!    R_out             -- outer radius
!    R_outann          -- outer annulus radius
!    R_outdust         -- outer radius
!    R_ref             -- reference radius
!    accr1             -- central star accretion radius
!    accr2             -- perturber accretion radius
!    alphaSS           -- desired alphaSS
!    annulus_m         -- mass in annulus
!    bhspin            -- black hole spin
!    bhspinangle       -- black hole spin angle
!    binary_O          -- Omega, PA of ascending node (deg)
!    binary_a          -- binary semi-major axis
!    binary_e          -- binary eccentricity
!    binary_f          -- f, initial true anomaly (deg,180=apastron)
!    binary_i          -- i, inclination (deg)
!    binary_w          -- w, argument of periapsis (deg)
!    deltat            -- output interval as fraction of orbital period
!    disc_m            -- disc mass
!    dist_unit         -- distance unit (e.g. au,pc,kpc,0.1pc)
!    dust_method       -- dust method (1=one fluid,2=two fluid)
!    dust_to_gas_ratio -- dust to gas ratio
!    flyby_a           -- flyby periastron distance
!    flyby_d           -- initial distance of flyby (in units of periastron distance)
!    flyby_r           -- roll angle of flyby
!    graindensinp      -- intrinsic grain density (in g/cm^3)
!    grainsizeinp      -- grain size (in cm)
!    ibinary           -- closed binary or flyby (0=binary,1=flyby)
!    ipotential        -- potential (1=central point mass,2=binary potential,3=spinning black hole)
!    itaperdust        -- exponentially taper the outer disc profile
!    itapergas         -- exponentially taper the outer disc profile
!    m1                -- central star mass
!    m2                -- perturber mass
!    mass_set          -- how to set gas density profile
!    mass_unit         -- mass unit (e.g. solarm,jupiterm,earthm)
!    norbits           -- maximum number of orbits at outer disc
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
!    sig_norm          -- sigma = sig_norm (R/R_ref)^-p (1-sqrt(R_in/R))
!    sig_ref           -- sigma at reference radius
!    xinc              -- inclination angle
!
!  DEPENDENCIES: centreofmass, dim, dust, eos, extern_binary,
!    extern_lensethirring, externalforces, infile_utils, io, kernel,
!    options, part, physcon, prompting, setbinary, setdisc, setflyby,
!    timestep, units
!+
!--------------------------------------------------------------------------
module setup
 use dim,     only:maxp,use_dust,maxalpha
 use options, only:use_dustfrac
 implicit none
 public  :: setpart
 integer :: np,np_dust,i,itype,npart_inter,nparttot,npartdust,norbits
 real    :: m1,m2,mcentral,accr1,accr2,bhspin,bhspinangle
 real    :: binary_a,binary_e,binary_i,binary_O,binary_w,binary_f,Rochelobe
 real    :: alphaSS,deltat,flyby_a,flyby_d,flyby_r
 integer :: icentral,ipotential,nsinks,ibinary
 character(len=20) :: disclabel
 character(len=*), dimension(3), parameter :: disctype = &
    (/'binary   ', &
      'primary  ', &
      'secondary'/)
 logical :: iuse_disc(3),ismoothgas(3),ismoothdust(3),itapergas(3),itaperdust(3),multiple_disc_flag
 integer :: mass_set(3),profile_set_dust,iprofilegas(3),iprofiledust(3)
 integer :: sigmaprofilegas(3),sigmaprofiledust(3)
 real    :: R_in(3),R_out(3),R_ref(3),R_c(3),pindex(3),qindex(3),H_R(3),xinc(3)
 real    :: disc_m(3),disc_mfac(3),sig_ref(3),sig_norm(3),annulus_m(3),R_inann(3),R_outann(3)
 real    :: R_indust(3),R_outdust(3),R_c_dust(3),pindex_dust(3),qindex_dust(3)
 real    :: H_R_dust(3),disc_mdust(3),sig_normdust(3)
 real    :: grainsizeinp,graindensinp,dust_to_gas_ratio
 integer :: dust_method
 integer, parameter :: maxplanets = 9
 integer :: nplanets,setplanets
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
 use io,                   only:master,warning,error,fatal
 use kernel,               only:hfact_default
 use options,              only:iexternalforce,ieos,alpha,icooling
 use part,                 only:nptmass,xyzmh_ptmass,maxvxyzu,vxyz_ptmass,ihacc,ihsoft,igas,idust,dustfrac
 use physcon,              only:au,solarm,jupiterm,pi,years
 use prompting,            only:prompt
 use setbinary,            only:set_binary,Rochelobe_estimate
 use setdisc,              only:set_disc,scaled_sigma,scaled_discmass
 use setflyby,             only:set_flyby,get_T_flyby
 use timestep,             only:tmax,dtmax
 use units,                only:set_units,select_unit,umass,udist,utime
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
 character(len=100)               :: filename
 logical :: iexist,questplanets,seq_exists
 real    :: phi,vphi,sinphi,cosphi,omega,r2,disc_m_within_r,period_longest
 real    :: jdust_to_gas_ratio,Rj,period,Rochesizei
 real    :: totmass_gas,totmass_dust,starmass
 real    :: sini,cosi,polyk_dust,xorigini(3),vorigini(3),alpha_returned(3)
 integer :: ierr,j,ndiscs,idisc,npingasdisc,npindustdisc

 print "(/,65('-'),2(/,a),/,65('-'),/)", &
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
       !--todo: check that potentials are implemented correctly
       !        add more potentials
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
             m2       = 0.1
             binary_a = 10.
             binary_e = 0.
             binary_i = 0.
             binary_O = 0.
             binary_w = 0.
             binary_f = 180.
             accr1    = 0.25*binary_a
             accr2    = 0.25*binary_a
          case (1)
             !--unbound (flyby)
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
    !--multiple disc options
    !
    print*, ''
    print "(a)",'================='
    print "(a)",'+++  DISC(S)  +++'
    print "(a)",'================='
    iuse_disc = .false.
    if ((icentral==1) .and. (nsinks==2)) then
       !--multiple discs possible
       multiple_disc_flag = .true.
       if (ibinary==0) then
          !--bound binary: circum-binary, -primary, -secondary
          ndiscs = 3
          iuse_disc(1) = .true.
       elseif (ibinary==1) then
          !--unbound binary (flyby): circum-primary, -secondary
          ndiscs = 2
          iuse_disc(2) = .true.
       endif
       do i=4-ndiscs,3
          call prompt('Do you want a circum'//trim(disctype(i))//' disc?',iuse_disc(i))
       enddo
       if (.not.any(iuse_disc)) call fatal('setup','need to setup at least one disc!')
       !--set number of discs
       ndiscs = count(iuse_disc)
    else
       !--only a single disc possible
       multiple_disc_flag = .false.
       iuse_disc(1) = .true.
       ndiscs = 1
    endif
    !
    !--set gas disc defaults
    !
    R_in      = 1.
    R_out     = 150.
    R_ref     = R_in
    R_c       = R_out
    mass_set  = 0
    itapergas = .false.
    pindex    = 1.
    qindex    = 0.25
    if (ndiscs > 1) qindex = 0.
    alphaSS   = 0.005
    xinc      = 0.            !--todo: implement general rotations around a position angle
    H_R       = 0.05
    disc_mfac = 1.
    if (multiple_disc_flag .and. (ibinary==0)) then
       !--set appropriate disc radii for bound binary
       Rochelobe = Rochelobe_estimate(m1,m2,binary_a)
       R_in      = (/binary_a + Rochelobe + accr1, accr1, accr2/)
       R_out     = (/10.*R_in(1), &
                     binary_a - Rochelobe - accr1, &
                     binary_a - Rochelobe - accr2/)
       R_ref     = R_in
       R_c       = R_out
       disc_mfac = (/1., 0.1, 0.01/)
    endif
    do i=1,3
       if (iuse_disc(i)) then
          if (multiple_disc_flag) then
             print*, ''
             print "(a)",'>>>  circum'//trim(disctype(i))//' disc  <<<'
          endif
          call prompt('How do you want to set the gas disc mass?'//new_line('A')// &
                      ' 0=total disc mass'//new_line('A')// &
                      ' 1=mass within annulus'//new_line('A')// &
                      ' 2=surface density normalisation'//new_line('A')// &
                      ' 3=surface density at reference radius'//new_line('A'),mass_set(i),0,3)
          call prompt('Do you want to exponentially taper the outer gas disc profile?',itapergas(i))
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
          end select
       endif
    enddo
    !
    !--set dust disc defaults
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
       print*, ''
       print "(a)",'=============='
       print "(a)",'+++  DUST  +++'
       print "(a)",'=============='
       !
       !--dust method
       !
       dust_method  = 2
       use_dustfrac = .false.
       call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
       if (dust_method==1) use_dustfrac = .true.
       call prompt('How do you want to set the dust density profile? (0=equal to gas,1=custom)',profile_set_dust,0,1)
       call prompt('Enter dust to gas ratio',dust_to_gas_ratio,0.)
       call prompt('Enter grain size in cm',grainsizeinp,0.)
    endif
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
    print*, ''
    print "(a)",'================='
    print "(a)",'+++  PLANETS  +++'
    print "(a)",'================='
    call prompt('Do you want to add planets?',questplanets)
    if (questplanets) then
       setplanets = 1
       call prompt('Enter the number of planets',nplanets,0,maxplanets)
    endif
    !
    !--determine simulation time
    !
    print*, ''
    print "(a)",'================'
    print "(a)",'+++  OUTPUT  +++'
    print "(a)",'================'
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
    else
       call prompt('Enter time between dumps as fraction of outer disc orbital time',deltat,0.)
       call prompt('Enter number of orbits to simulate',norbits,0)
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
 !--index of disc (if only one)
 !
 idisc = 0
 if (ndiscs==1) then
    do i=1,3
       if (iuse_disc(i).eqv..true.) idisc = i
    enddo
 endif
 !
 !--equation of state
 !
 if (ndiscs > 1) then
    !--multiple discs
    if (maxvxyzu==3) then
       !--force globally isothermal
       if (sum(qindex) > maxval(qindex)) then
          call fatal('setup_disc','locally isothermal eos for more than one disc '// &
                     'requested, no ieos to handle this')
       else
          !--isothermal
          ieos = 1
          gamma = 1.0
       endif
    else
       !--adiabatic
       ieos = 2
       gamma = 5./3.
       icooling = 1
    endif
 else
    !--single disc
    if (maxvxyzu==3 .and. qindex(idisc) > 0.) then
       !--locally isothermal
       ieos = 3
       gamma = 1.0
    else
       !--adiabatic
       ieos = 2
       gamma = 5./3.
       icooling = 1
    endif
 endif
#ifdef MCFOST
 !--radiative equilibrium
 ieos = 2
 icooling = 0
 ipdv_heating = 0
 ishock_heating = 0
 alphau = 0
#endif
 !--sanity check on ieos = 6
 if (ieos==6 .and. isink==0) call fatal('setup_disc','something''s gone wrong with ieos & isink...')

 !
 !--surface density profile
 !
 ismoothgas      = .true.       !--smoothed at inner edge by default
 iprofilegas     = 0
 sigmaprofilegas = 0
 do i=1,3
    if (itapergas(i)) then
       iprofilegas(i)     = 1
       sigmaprofilegas(i) = 1
    endif
    if (ismoothgas(i))                    sigmaprofilegas(i) = 2
    if (itapergas(i) .and. ismoothgas(i)) sigmaprofilegas(i) = 3
 enddo
 if (use_dust) then
    ismoothdust      = .true.   !--smoothed at inner edge by default
    iprofiledust     = iprofilegas
    sigmaprofiledust = sigmaprofilegas
    do i=1,3
       if (itaperdust(i)) then
          iprofiledust(i)     = 1
          sigmaprofiledust(i) = 1
       endif
       if (ismoothdust(i))                     sigmaprofiledust(i) = 2
       if (itaperdust(i) .and. ismoothdust(i)) sigmaprofiledust(i) = 3
    enddo
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
       !--binary
       select case (ibinary)
       case (0)
          !--bound
          nptmass  = 0
          print*,'Setup the binary' !--the stars' barycentre is put on the origin
          print "(a)", 'Central objects represented by two sinks'
          call set_binary(m1,massratio=m2/m1,semimajoraxis=binary_a,eccentricity=binary_e, &
                          posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i, &
                          f=binary_f,accretion_radius1=accr1,accretion_radius2=accr2, &
                          xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass)
          mcentral = m1 + m2
       case (1)
          !--unbound (flyby)
          print*,'Disc around a single star with flyby'
          print "(a)", 'Central object represented by a sink at the system origin with a perturber sink'
          call set_flyby(mprimary=m1,massratio=m2/m1,dma=flyby_a,n0=flyby_d,roll=flyby_r, &
                         accretion_radius1=accr1,accretion_radius2=accr2, &
                         xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass)
          mcentral = m1
       end select
    end select
 end select

 !
 !--calculate total mass (for multiple discs)
 !
 totmass_gas  = 0.
 totmass_dust = 0.
 do i=1,3
    if (iuse_disc(i)) then
       !--compute the disc mass for different mass_set values
       if (mass_set(i)==1) then
          sig_norm(i) = annulus_m(i) / scaled_discmass(sigmaprofilegas(i),pindex(i),R_inann(i),R_outann(i),R_ref(i),R_c(i))
          disc_m(i)   = sig_norm(i)  * scaled_discmass(sigmaprofilegas(i),pindex(i),R_in(i),R_out(i),R_ref(i),R_c(i))
       endif
       if (mass_set(i)==2) disc_m(i) = sig_norm(i) * scaled_discmass(sigmaprofilegas(i),pindex(i),R_in(i),R_out(i),R_ref(i),R_c(i))
       if (mass_set(i)==3) then
          if (.not.(R_in(i) < R_ref(i)) .and. ismoothgas(i)) call fatal('set_disc', &
             'if smoothing at inner disc edge and setting disc mass by sigma(R_ref), must have R_in < R_ref')
          sig_norm(i) = sig_ref(i)  / scaled_sigma(R_ref(i),sigmaprofilegas(i),pindex(i),R_ref(i),R_in(i),R_c(i))
          disc_m(i)   = sig_norm(i) * scaled_discmass(sigmaprofilegas(i),pindex(i),R_in(i),R_out(i),R_ref(i),R_c(i))
       endif
       totmass_gas = totmass_gas + disc_m(i)
       if (use_dust) then
          disc_mdust(i) = disc_m(i) * dust_to_gas_ratio
          totmass_dust  = totmass_dust + disc_mdust(i)
       endif
    endif
 enddo
 print*,' Total gas mass of system =  ',totmass_gas
 if (use_dust) print*,' Total dust mass of system = ',totmass_dust

 !
 !--setup disc(s)
 !
 time  = 0.
 hfact = hfact_default
 xinc = xinc*(pi/180.0)
 if (maxalpha==0) alpha = alphaSS
 if (use_dust) then
    grainsizecgs = grainsizeinp
    graindenscgs = graindensinp
 endif
 nparttot  = 0
 npartdust = 0
 do i=1,3
    if (iuse_disc(i)) then

       !--set disc origin
       if (multiple_disc_flag) print "(/,a)",'>>> Setting up circum'//trim(disctype(i))//' disc <<<'
       select case(i)
       case (1)
          !--single disc or circumbinary
          !  centre of mass of binary defined to be zero (see set_binary)
          xorigini(:) = 0.
          vorigini(:) = 0.
          starmass    = mcentral
          Rochesizei  = huge(0.)
       case(2)
          !--circumprimary
          xorigini(:) = xyzmh_ptmass(1:3,1)
          vorigini(:) = vxyz_ptmass(1:3,1)
          starmass    = m1
          Rochesizei  = binary_a - Rochelobe
       case(3)
          !--circumsecondary
          xorigini(:) = xyzmh_ptmass(1:3,2)
          vorigini(:) = vxyz_ptmass(1:3,2)
          starmass    = m2
          Rochesizei  = Rochelobe
       end select
       if (multiple_disc_flag .and. ibinary==0) then
          if (R_out(i) > Rochesizei .and. R_out(i) < binary_a) then
             print "(/,a,/)",'*** WARNING: Outer disc radius for circum'//trim(disctype(i))// &
                              ' > Roche lobe of '//trim(disctype(i))//' ***'
          endif
       endif

       !--set disc(s)
       if (use_dust .and. use_dustfrac) then
          !--gas and dust mixture disc
          npingasdisc = int(disc_m(i)/totmass_gas*np)
          call set_disc(id,master        = master,             &
                        mixture          = .true.,             &
                        npart            = npingasdisc,        &
                        npart_start      = nparttot + 1,       &
                        rref             = R_ref(i),           &
                        rmin             = R_in(i),            &
                        rmax             = R_out(i),           &
                        rmindust         = R_indust(i),        &
                        rmaxdust         = R_outdust(i),       &
                        indexprofile     = iprofilegas(i),     &
                        indexprofiledust = iprofiledust(i),    &
                        rc               = R_c(i),             &
                        rcdust           = R_c_dust(i),        &
                        p_index          = pindex(i),          &
                        p_indexdust      = pindex_dust(i),     &
                        q_index          = qindex(i),          &
                        HoverR           = H_R(i),             &
                        disc_mass        = disc_m(i),          &
                        disc_massdust    = disc_mdust(i),      &
                        star_mass        = starmass,           &
                        gamma            = gamma,              &
                        particle_mass    = massoftype(igas),   &
                        xyz_origin       = xorigini(:),        &
                        vxyz_origin      = vorigini(:),        &
                        hfact            = hfact,              &
                        xyzh             = xyzh,               &
                        vxyzu            = vxyzu,              &
                        polyk            = polyk,              &
                        alpha            = alpha,              &
                        ismooth          = ismoothgas(i),      &
                        inclination      = xinc(i),            &
                        twist            = .false.,            &
                        prefix           = fileprefix)
          !--set dustfrac
          sig_norm(i)     = disc_m(i)     / scaled_discmass(sigmaprofilegas(i),pindex(i),R_in(i),R_out(i),R_ref(i),R_c(i))
          sig_normdust(i) = disc_mdust(i) / scaled_discmass(sigmaprofiledust(i),pindex_dust(i),R_indust(i),R_outdust(i),R_ref(i),R_c_dust(i))
          do j=nparttot+1,npingasdisc
             Rj = sqrt(dot_product(xyzh(1:2,j)-xorigini(1:2),xyzh(1:2,j)-xorigini(1:2)))
             if (Rj<R_indust(i) .or. Rj>R_outdust(i)) then
                jdust_to_gas_ratio = 0.
             else
                call get_dust_to_gas_ratio(jdust_to_gas_ratio,Rj,sigmaprofilegas(i),sigmaprofiledust(i), &
                                           sig_norm(i),sig_normdust(i),pindex(i),pindex_dust(i), &
                                           R_in(i),R_ref(i),R_c(i),R_indust(i),R_c_dust(i))
             endif
             call set_dustfrac(jdust_to_gas_ratio,dustfrac(j))
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
                        star_mass       = starmass,           &
                        gamma           = gamma,              &
                        particle_mass   = massoftype(igas),   &
                        xyz_origin      = xorigini(:),        &
                        vxyz_origin     = vorigini(:),        &
                        hfact           = hfact,              &
                        xyzh            = xyzh,               &
                        vxyzu           = vxyzu,              &
                        polyk           = polyk,              &
                        alpha           = alpha,              &
                        ismooth         = ismoothgas(i),      &
                        inclination     = xinc(i),            &
                        twist           = .false.,            &
                        prefix          = fileprefix)
          nparttot = nparttot + npingasdisc
          if (use_dust) then
             !--dust disc
             npindustdisc = int(disc_mdust(i)/totmass_dust*np_dust)
             call set_disc(id,master     = master,             &
                           npart         = npindustdisc,       &
                           npart_start   = nparttot + 1,       &
                           particle_type = idust,              &
                           rref          = R_ref(i),           &
                           rmin          = R_indust(i),        &
                           rmax          = R_outdust(i),       &
                           indexprofile  = iprofiledust(i),    &
                           rc            = R_c_dust(i),        &
                           p_index       = pindex_dust(i),     &
                           q_index       = qindex_dust(i),     &
                           HoverR        = H_R_dust(i),        &
                           disc_mass     = disc_mdust(i),      &
                           star_mass     = starmass,           &
                           gamma         = gamma,              &
                           particle_mass = massoftype(idust),  &
                           xyz_origin    = xorigini(:),        &
                           vxyz_origin   = vorigini(:),        &
                           hfact         = hfact,              &
                           xyzh          = xyzh,               &
                           vxyzu         = vxyzu,              &
                           polyk         = polyk_dust,         &
                           ismooth       = ismoothdust(i),     &
                           inclination   = xinc(i),            &
                           twist         = .false.,            &
                           prefix        = fileprefix)
             nparttot  = nparttot  + npindustdisc
             npartdust = npartdust + npindustdisc
             !--reset qfacdisc to gas disc value
             qfacdisc = qindex(i)
          endif
       endif
       !--reset alpha for each disc
       alpha_returned(i) = alpha
       alpha = alphaSS

    endif
 enddo

 !--number of particles
 npart = nparttot
 npartoftype(igas)  = nparttot - npartdust
 npartoftype(idust) = npartdust

 !--alpha viscosity
 if (ndiscs==1) then
    do i=1,3
       if (iuse_disc(i)) alpha = alpha_returned(i)
    enddo
 else
    call warning('setup_disc','multiple discs: cannot use alpha for alpha_SS, setting equal to 0.1 instead')
    alpha = 0.1
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
    !--bound binary
    period = sqrt(4.*pi**2*binary_a**3/mcentral)
 elseif (icentral==1 .and. nsinks==2 .and. ibinary==1) then
    !--unbound binary (flyby)
    period = get_T_flyby(m1,m2,flyby_a,flyby_d)
    norbits = 1
 elseif (setplanets==1) then
    !--outer planet set above
 else
    !--outer disc
    period = sqrt(4.*pi**2*R_out(idisc)**3/mcentral)
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
 logical :: done_alpha

 done_alpha = .false.

 print*, ''
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
       call write_inopt(ibinary,'ibinary','binary: bound or unbound [flyby] (0=bound,1=unbound)',iunit)
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
          call write_inopt(flyby_a,'flyby_a','flyby periastron distance',iunit)
          call write_inopt(flyby_d,'flyby_d','initial distance of flyby (in units of periastron distance)',iunit)
          call write_inopt(flyby_r,'flyby_r','roll angle of flyby',iunit)
       end select
    end select
 end select
 !--multiple disc options
 if (multiple_disc_flag) then
    write(iunit,"(/,a)") '# options for multiple discs'
    if (ibinary==0) then
       call write_inopt(iuse_disc(1),'use_'//trim(disctype(1))//'disc','setup circum'//trim(disctype(1))//' disc',iunit)
    endif
    do i=2,3
       call write_inopt(iuse_disc(i),'use_'//trim(disctype(i))//'disc','setup circum'//trim(disctype(i))//' disc',iunit)
    enddo
 endif
 !--individual disc(s)
 do i=1,3
    if (iuse_disc(i)) then
       if (multiple_disc_flag) then
          disclabel = disctype(i)
       else
          disclabel = ''
       endif
       !--gas disc
       if (multiple_disc_flag) then
          write(iunit,"(/,a)") '# options for circum'//trim(disclabel)//' gas disc'
       else
          write(iunit,"(/,a)") '# options for gas accretion disc'
       endif
       call write_inopt(mass_set(i),'mass_set'//trim(disclabel),'how to set gas density profile' // &
          ' (0=total disc mass,1=mass within annulus,2=surface density normalisation,' // &
          '3=surface density at reference radius)',iunit)
       call write_inopt(itapergas(i),'itapergas'//trim(disclabel),'exponentially taper the outer disc profile',iunit)
       call write_inopt(R_in(i),'R_in'//trim(disclabel),'inner radius',iunit)
       call write_inopt(R_ref(i),'R_ref'//trim(disclabel),'reference radius',iunit)
       call write_inopt(R_out(i),'R_out'//trim(disclabel),'outer radius',iunit)
       if (itapergas(i)) call write_inopt(R_c(i),'R_c'//trim(disclabel),'characteristic radius of the exponential taper',iunit)
       select case (mass_set(i))
       case (0)
          call write_inopt(disc_m(i),'disc_m'//trim(disclabel),'disc mass',iunit)
       case (1)
          call write_inopt(annulus_m(i),'annulus_m'//trim(disclabel),'mass within annulus',iunit)
          call write_inopt(R_inann(i),'R_inann'//trim(disclabel),'inner annulus radius',iunit)
          call write_inopt(R_outann(i),'R_outann'//trim(disclabel),'outer annulus radius',iunit)
       case (2)
          if (itapergas(i)) then
             call write_inopt(sig_norm(i),'sig_norm'//trim(disclabel),'sigma = sig_norm (R/R_ref)^-p exp[-(R/R_c)^(2-p)] (1-sqrt(R_in/R))',iunit)
          else
             call write_inopt(sig_norm(i),'sig_norm'//trim(disclabel),'sigma = sig_norm (R/R_ref)^-p (1-sqrt(R_in/R))',iunit)
          endif
       case (3)
          call write_inopt(sig_ref(i),'sig_ref'//trim(disclabel),'sigma at reference radius',iunit)
       end select
       call write_inopt(pindex(i),'pindex'//trim(disclabel),'p index',iunit)
       call write_inopt(qindex(i),'qindex'//trim(disclabel),'q index',iunit)
       call write_inopt(xinc(i),'xinc'//trim(disclabel),'inclination angle',iunit)
       call write_inopt(H_R(i),'H_R'//trim(disclabel),'H/R at R=R_ref',iunit)
       if (.not.done_alpha) then
          if (maxalpha==0) call write_inopt(alphaSS,'alphaSS','desired alphaSS',iunit)
          done_alpha = .true.
       endif
       !--dust disc
       if (profile_set_dust==1) then
          if (multiple_disc_flag) then
             write(iunit,"(/,a)") '# options for circum'//trim(disclabel)//' dust disc'
          else
             write(iunit,"(/,a)") '# options for dust accretion disc'
          endif
          call write_inopt(itaperdust(i),'itaperdust'//trim(disclabel),'exponentially taper the outer disc profile',iunit)
          call write_inopt(R_indust(i),'R_indust'//trim(disclabel),'inner radius',iunit)
          call write_inopt(R_outdust(i),'R_outdust'//trim(disclabel),'outer radius',iunit)
          if (iprofiledust(i)==1) call write_inopt(R_c_dust(i),'R_c_dust'//trim(disclabel),'characteristic radius of the exponential taper',iunit)
          call write_inopt(pindex_dust(i),'pindex_dust'//trim(disclabel),'p index',iunit)
          if (.not. use_dustfrac) then
             call write_inopt(qindex_dust(i),'qindex_dust'//trim(disclabel),'q index',iunit)
             call write_inopt(H_R_dust(i),'H_R_dust'//trim(disclabel),'H/R at R=R_ref',iunit)
          endif
       endif
    endif
 enddo
 !--dust options
 if (use_dust) then
    write(iunit,"(/,a)") '# options for dust'
    call write_inopt(dust_method,'dust_method','dust method (1=one fluid,2=two fluid)',iunit)
    call write_inopt(dust_to_gas_ratio,'dust_to_gas_ratio','dust to gas ratio',iunit)
    call write_inopt(profile_set_dust,'profile_set_dust','how to set dust density profile (0=equal to gas, 1=custom)',iunit)
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
 !--timestepping
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

 !--dust method
 if (use_dust) then
    call read_inopt(dust_method,'dust_method',db,min=1,max=2,errcount=nerr)
    if (dust_method==1) use_dustfrac = .true.
 endif
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
          call read_inopt(flyby_r,'flyby_r',db,min=0.,errcount=nerr)
       end select
    end select
 end select
 !--multiple discs
 multiple_disc_flag = .false.
 iuse_disc = .false.
 if ((icentral==1) .and. (nsinks==2)) then
    multiple_disc_flag = .true.
    if (ibinary==0) then
       call read_inopt(iuse_disc(1),'use_'//trim(disctype(1))//'disc',db,errcount=nerr)
    endif
    do i=2,3
       call read_inopt(iuse_disc(i),'use_'//trim(disctype(i))//'disc',db,errcount=nerr)
    enddo
 else
    iuse_disc(1) = .true.
 endif
 do i=1,3
    if (iuse_disc(i)) then
       if (multiple_disc_flag) then
          disclabel = disctype(i)
       else
          disclabel = ''
       endif
       !--gas disc
       call read_inopt(R_in(i),'R_in'//trim(disclabel),db,min=0.,errcount=nerr)
       call read_inopt(R_out(i),'R_out'//trim(disclabel),db,min=R_in(i),errcount=nerr)
       call read_inopt(R_ref(i),'R_ref'//trim(disclabel),db,min=R_in(i),errcount=nerr)
       call read_inopt(itapergas(i),'itapergas'//trim(disclabel),db,errcount=nerr)
       call read_inopt(mass_set(i),'mass_set'//trim(disclabel),db,min=0,max=3,errcount=nerr)
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
       end select
       call read_inopt(pindex(i),'pindex'//trim(disclabel),db,errcount=nerr)
       call read_inopt(qindex(i),'qindex'//trim(disclabel),db,errcount=nerr)
       call read_inopt(xinc(i),'xinc'//trim(disclabel),db,min=0.,max=180.,errcount=nerr)
       call read_inopt(H_R(i),'H_R'//trim(disclabel),db,min=0.,errcount=nerr)
       !--dust disc
       select case (profile_set_dust)
       case (0)
          R_indust(i)         = R_in(i)
          R_outdust(i)        = R_out(i)
          pindex_dust(i)      = pindex(i)
          qindex_dust(i)      = qindex(i)
          H_R_dust(i)         = H_R(i)
          itaperdust(i)       = itapergas(i)
          R_c_dust(i)         = R_c(i)
       case (1)
          call read_inopt(R_indust(i),'R_indust'//trim(disclabel),db,min=R_in(i),errcount=nerr)
          call read_inopt(R_outdust(i),'R_outdust'//trim(disclabel),db,min=R_indust(i),max=R_out(i),errcount=nerr)
          call read_inopt(pindex_dust(i),'pindex_dust'//trim(disclabel),db,errcount=nerr)
          call read_inopt(itaperdust(i),'itaperdust'//trim(disclabel),db,errcount=nerr)
          if (itaperdust(i)) then
             call read_inopt(R_c_dust(i),'R_c_dust'//trim(disclabel),db,min=0.,errcount=nerr)
          endif
          if (.not. use_dustfrac) then
             call read_inopt(qindex_dust(i),'qindex_dust'//trim(disclabel),db,err=ierr,errcount=nerr)
             if (ierr /= 0) qindex_dust(i) = qindex(i)
             call read_inopt(H_R_dust(i),'H_R_dust'//trim(disclabel),db,min=0.,err=ierr,errcount=nerr)
             if (ierr /= 0) H_R_dust(i) = H_R(i)
          endif
       end select
    endif
 enddo
 if (maxalpha==0) call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
 !--dust
 if (use_dust) then
    call read_inopt(dust_to_gas_ratio,'dust_to_gas_ratio',db,min=0.,errcount=nerr)
    call read_inopt(profile_set_dust,'profile_set_dust',db,min=0,max=1,errcount=nerr)
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
 !--timestepping
 !  following two are optional: not an error if not present
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
subroutine get_dust_to_gas_ratio(dust_to_gas,R,sigmaprofilegas,sigmaprofiledust, &
                                 sig_norm,sig_normdust,pindex,pindex_dust, &
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
