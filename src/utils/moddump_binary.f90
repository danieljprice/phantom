!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! take a snapshot containing a single star and convert it into a binary
! system, either by adding a sink particle companion, or by adding a
! second star from another dumpfile
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - a1                   : *1st companion semi-major axis (triple) [code units]*
!   - a2                   : *2nd companion semi-major axis (triple) [code units]*
!   - cmass1               : *1st companion mass (triple) [code units]*
!   - cmass2               : *2nd companion mass (triple) [code units]*
!   - comp_shift           : *code units to shift companion (+ve towards primary)*
!   - companion_hsoft      : *softening length for companion [code units]*
!   - companion_mass       : *companion mass [Msun]*
!   - densityfile          : *filename of the input stellar profile*
!   - ecc                  : *orbital eccentricity*
!   - gwinspiral           : *add gravitational radiation reaction*
!   - hacc                 : *accretion radius for the companion [code units]*
!   - hacc1                : *accretion radius for the primary [code units]*
!   - hacc2                : *accretion radius for the 1st companion [code units]*
!   - hacc3                : *accretion radius for the 2nd companion [code units]*
!   - hacc_sec             : *accretion radius of secondary [Rsun]*
!   - hsoft_core           : *softening length of the created point-mass core*
!   - hsoft_primary        : *softening length for primary [code units]*
!   - hsoft_sec            : *softening length of secondary [Rsun]*
!   - hsoft_secondary      : *softening length for secondary [code units]*
!   - hsoft_tertiary       : *softening length for tertiary [code units]*
!   - infile_name          : *name of the .in file to read omega_corotate from*
!   - iprimary_grav        : *replace primary core with a fixed potential (if 1 sink)*
!   - iproperty            : *sink property index to (re)set (0 = none)*
!   - iremove              : *which sink to remove (2 or 3)*
!   - iselect              : *which sink particle to (re)set properties for*
!   - mcomp                : *companion mass [code units]*
!   - mcut                 : *mass of the created point-mass core [code units]*
!   - nstar2               : *number of particles in the second dumpfile (0 = same as star 1)*
!   - operation            : *operation to perform (1-13)*
!   - propval              : *new value for the selected sink property [solar units]*
!   - second_dumpfile      : *name of the second dumpfile*
!   - separation           : *orbital separation [Rsun]*
!   - use_corotating_frame : *transform to a corotating frame*
!   - vel_shift            : *velocity to add in the direction of the primary [code units]*
!
! :Dependencies: centreofmass, dim, eos, extern_corotate,
!   extern_gwinspiral, externalforces, infile_utils, io, moddump_utils,
!   options, part, physcon, prompting, readwrite_dumps, readwrite_mesa,
!   setbinary, table_utils, timestep, units, vectorutils
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 !
 ! operation selector. The original moddump auto-selected a menu based on the
 ! number of sink particles present; here the operation is chosen explicitly
 ! and the sink-count requirement is validated inside modify_dump.
 !
 !  1) add a sink companion (binary)                    [nptmass <= 1]
 !  2) add a magnetic field in the star                 [nptmass <= 1]
 !  3) cut profile to create a sink in the core         [nptmass <= 1]
 !  4) manually create a sink in the core               [nptmass <= 1]
 !  5) set up a triple system                           [nptmass <= 1]
 !  6) set up star for relaxation in corotating frame   [nptmass <= 1]
 !  7) set up binary after relaxation in corotating fr. [nptmass <= 1]
 !  8) add a second star from another dumpfile          [nptmass <= 1]
 !  9) remove a sink from the simulation                [nptmass == 3]
 ! 10) transform from corotating to inertial frame      [nptmass == 2]
 ! 11) shift companion position in corotating frame     [nptmass == 2]
 ! 12) add velocity to companion                        [nptmass == 2]
 ! 13) (re)set sink properties                          [nptmass >= 1]
 !
 integer :: operation = 1

 ! operations 1 & 8 (add companion / binary from dumpfile)
 real    :: mcomp           = 0.6
 real    :: a1              = 100.
 real    :: ecc             = 0.
 real    :: hacc            = 0.
 real    :: companion_hsoft = 0.
 logical :: use_corotating_frame = .false.
 logical :: gwinspiral = .false.
 character(len=120) :: second_dumpfile = ''
 integer :: nstar2 = 0

 ! operations 3 & 4 (create core sink)
 real    :: mcut       = 0.
 real    :: hsoft_core = 3.
 character(len=120) :: densityfile = 'P12_Phantom_Profile.data'

 ! operation 5 (triple)
 real    :: cmass1 = 0.0095, cmass2 = 0.0095
 real    :: a2 = 336.8
 real    :: hacc1 = 0., hacc2 = 0., hacc3 = 0.
 real    :: hsoft_primary = 0., hsoft_secondary = 0., hsoft_tertiary = 0.

 ! operations 6 & 7 (corotating-frame relaxation / binary after relaxation)
 real    :: companion_mass = 1.26
 real    :: separation     = 865.24
 real    :: hacc_sec       = 0.
 real    :: hsoft_sec      = 0.
 logical :: iprimary_grav  = .false.
 character(len=120) :: infile_name = 'binary.in'

 ! operation 9 (remove sink)
 integer :: iremove = 2

 ! operation 11 / 12 (shift / add velocity to companion)
 real    :: comp_shift = 100.
 real    :: vel_shift  = 0.

 ! operation 13 ((re)set sink properties)
 integer :: iselect   = 1
 integer :: iproperty = 0
 real    :: propval   = 0.

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas, &
                             delete_dead_or_accreted_particles,mhd,rho,shuffle_part,&
                             kill_particle,copy_particle
 use setbinary,         only:set_binary
 use units,             only:utime
 use physcon,           only:pi
 use centreofmass,      only:reset_centreofmass
 use options,           only:iexternalforce
 use externalforces,    only:omega_corotate,iext_corotate,iext_gwinspiral
 use extern_corotate,   only:icompanion_grav,companion_xpos,companion_mass_ext=>companion_mass, &
                             primarycore_xpos,primarycore_mass,primarycore_hsoft,hsoft
 use extern_gwinspiral, only:Nstar_gw
 use infile_utils,      only:open_db_from_file,inopts,read_inopt,close_db
 use table_utils,       only:yinterp
 use readwrite_mesa,    only:read_mesa
 use dim,               only:maxptmass,maxp,nsinkproperties
 use io,                only:fatal,idisk1,iprint
 use timestep,          only:tmax,dtmax
 use readwrite_dumps,   only:read_dump
 use eos,               only:X_in,Z_in

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                   :: i,ierr,irhomax,n
 integer                   :: nstar1,nptmass1,nptmass2,iprim,isec
 real                      :: primary_mass,mass_ratio,m1,m2,pmass1,pmass2
 real                      :: mass_donor,newCoM,period,a,primarycore_xpos_old
 real                      :: mcore,sink_dist
 real                      :: rcut,Mstar,radi,rhopart,rhomax
 real                      :: time2,hfact2
 real                      :: xyzmh1_stash(nsinkproperties),xyzmh2_stash(nsinkproperties),vxyz1_stash(3),vxyz2_stash(3)
 real, allocatable         :: r(:),den(:),pres(:),temp(:),enitab(:),Xfrac(:),Yfrac(:),mu(:),m(:)
 character(len=120)        :: dumpname
 type(inopts), allocatable :: db(:)

 !
 ! prompt only when the driver did not find a parameter file
 !
 if (prompt_for_params) call read_interactive_moddumpfile()

 !
 ! validate that the requested operation is compatible with the sink count
 !
 if (nptmass > 3) call fatal('moddump_binary','Number of sink particles > 3')
 select case(operation)
 case(9)
    if (nptmass /= 3) call fatal('moddump_binary','operation 9 requires 3 sink particles')
 case(10,11,12)
    if (nptmass /= 2) call fatal('moddump_binary','operations 10-12 require 2 sink particles')
 case(13)
    if (nptmass < 1) call fatal('moddump_binary','operation 13 requires at least 1 sink particle')
 case(1:8)
    if (nptmass > 1) call fatal('moddump_binary','operations 1-8 require 1 or fewer sink particles')
 end select

 select case(operation)
    !
    !--------------------------------------------------------------------
 case(1,8)  ! add a sink companion, or add a star from another dumpfile
    !--------------------------------------------------------------------
    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
    call delete_dead_or_accreted_particles(npart,npartoftype)
    nptmass1 = nptmass
    nstar1 = npart
    pmass1 = massoftype(igas)
    m1 = nstar1 * pmass1
    if (nptmass1 > 1) then
       call fatal('moddump_binary', 'unexpected number of sink particles in dump file (nptmass > 1)')
    elseif (nptmass1 == 1) then  ! there is a sink stellar core
       m1 = m1 + xyzmh_ptmass(4,1)
       xyzmh1_stash(1:nsinkproperties) = xyzmh_ptmass(1:nsinkproperties,1)
       vxyz1_stash(1:3) = vxyz_ptmass(1:3,1)
       print*,'Dump contains one sink particle with m=',xyzmh1_stash(4),&
               ', hacc=',xyzmh1_stash(ihacc),', and hsoft=',xyzmh1_stash(ihsoft)
    endif
    m2 = mcomp

    print*, 'Current primary mass in code units is ',m1
    print*, 'Companion mass in code units is ',m2

    ! set the binary
    if (use_corotating_frame) then
       iexternalforce = iext_corotate  !turns on corotation
       call set_binary(m1,m2,a1,ecc,xyzmh1_stash(ihacc),hacc,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate)
       print "(/,a,es18.10,/)", ' The angular velocity in the corotating frame is: ', omega_corotate
       ! set all gas velocities in the corotating frame to 0 (corotating binary)
       do i=1,npart
          vxyzu(1:3,i) = 0.
       enddo
    else ! non corotating frame
       if (gwinspiral) iexternalforce = iext_gwinspiral
       call set_binary(m1,m2,a1,ecc,xyzmh1_stash(ihacc),hacc,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr)
    endif

    if (nptmass1 == 1) then
       iprim = 2
       isec = 3
    else ! nptmass1 = 0
       iprim = 1
       isec = 2
    endif

    !shifts star 1 gas to the primary sink
    do i=1,npart
       xyzh(1:3,i) = xyzh(1:3,i) + xyzmh_ptmass(1:3,iprim)
       vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz_ptmass(1:3,iprim)
    enddo

    ! Store sink velocity and position in binary orbit
    xyzmh1_stash(1:3) = xyzmh_ptmass(1:3,iprim)
    vxyz1_stash(1:3) = vxyz_ptmass(1:3,iprim)
    xyzmh2_stash(1:3) = xyzmh_ptmass(1:3,isec)
    vxyz2_stash(1:3) = vxyz_ptmass(1:3,isec)

    if (operation == 8) then
       dumpname = second_dumpfile
       if (nstar2 <= 0) nstar2 = nstar1

       ! Move star 1 particles to avoid getting overwritten when reading second dump file.
       if (2*nstar1 > maxp) then  ! Check if particle array is large enough to provide particle-copying buffer
          call fatal('moddump_binary','Two times number of particles in star 1 > array size. Run with --maxp=N '//&
                     'where N is desired number of particles')
       endif
       if (nstar1 > nstar2) then ! Move ith particle of star 1 to nstar1+i
          do i=1,nstar1
             call copy_particle(i,nstar1+i,.false.)
          enddo
       else ! Move ith particle of star 1 to nstar2+i
          do i=1,nstar1
             call copy_particle(i,nstar2+i,.false.)
          enddo
       endif

       ! read dump file containing star 2
       call read_dump(trim(dumpname),time2,hfact2,idisk1+1,iprint,0,1,ierr)
       if (ierr /= 0) call fatal('read_dump','error reading second dump file')
       nptmass2 = nptmass
       if (nptmass2 > 1) then
          call fatal('moddump_binary', 'unexpected number of sink particles in second dump file (nptmass > 1)')
       elseif (nptmass == 1) then
          xyzmh2_stash(4:nsinkproperties) = xyzmh_ptmass(4:nsinkproperties,1)
       endif

       pmass2 = massoftype(igas)
       if ( abs(1.-pmass2/pmass1) > 1.e-3) then
          call fatal('moddump_binary','unequal mass particles between dumps 1 and 2, pmass2 /= pmass1')
       endif
       print*,'Setting gas mass to be that from first dump,',pmass1
       massoftype(igas) = pmass1

       if (nstar1 > nstar2) then ! Move ith particle of star 1 to nstar2+i
          do i=1,nstar1
             call copy_particle(nstar1+i,nstar2+i,.false.)
          enddo
       endif

       npart = nstar1 + nstar2
       npartoftype(igas) = npart
       nptmass = nptmass1 + nptmass2

       ! shift star 2 gas to secondary sink
       do i=1,nstar2
          xyzh(1:3,i) = xyzh(1:3,i) + xyzmh2_stash(1:3)
          vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz2_stash(1:3)
       enddo

       if (nptmass2 == 1) then
          xyzmh_ptmass(1:nsinkproperties,nptmass1+nptmass2) = xyzmh2_stash(1:nsinkproperties)
          vxyz_ptmass(1:3,nptmass1+nptmass2) = vxyz2_stash(1:3)
       endif

       if (gwinspiral) then
          Nstar_gw(1) = nstar2
          Nstar_gw(2) = nstar1
       endif

    else
       nptmass = nptmass1 + 1
       xyzmh_ptmass(1:3,nptmass) = xyzmh2_stash(1:3)
       vxyz_ptmass(1:3,nptmass) = vxyz2_stash(1:3)
       xyzmh_ptmass(4,nptmass) = m2
       xyzmh_ptmass(ihacc,nptmass) = hacc
       xyzmh_ptmass(ihsoft,nptmass) = companion_hsoft

       if (gwinspiral) Nstar_gw(1) = npart
    endif

    if (nptmass1 == 1) then
       xyzmh_ptmass(1:nsinkproperties,1) = xyzmh1_stash(1:nsinkproperties)
       vxyz_ptmass(1:3,1) = vxyz1_stash(1:3)
    endif

    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    !--------------------------------------------------------------------
 case(2)  ! add a magnetic field in the star
    !--------------------------------------------------------------------
    if (mhd) then
       print "(/,a,/)", 'Automatic insertion of the magnetic field through the setBfield module'
    else
       print "(/,a,/)", 'Code not compiled with MHD=yes, no changes to the dump have been made'
    endif

    !--------------------------------------------------------------------
 case(3)  ! cut profile to create a sink in the core
    !--------------------------------------------------------------------
    call read_mesa(densityfile,den,r,pres,m,enitab,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr,cgsunits=.false.)
    rcut = yinterp(r,m,mcut)

    rhomax = 0.0
    irhomax = 1
    do i=1,npart
       rhopart = rho(i)
       if (rhopart > rhomax) then
          rhomax = rhopart
          irhomax = i
       endif
    enddo

    nptmass = nptmass + 1
    if (nptmass > maxptmass) call fatal('ptmass_create','nptmass > maxptmass')
    n = nptmass
    xyzmh_ptmass(:,n)   = 0.  ! zero all quantities by default
    xyzmh_ptmass(1:3,n) = xyzh(1:3,irhomax)
    xyzmh_ptmass(4,n)   = 0.  ! zero mass
    xyzmh_ptmass(ihsoft,n) = hsoft_core
    vxyz_ptmass(:,n) = 0.     ! zero velocity, get this by accreting

    do i=1,npart
       radi = sqrt((xyzh(1,i)-xyzh(1,irhomax))**2 + &
                (xyzh(2,i)-xyzh(2,irhomax))**2 + &
                (xyzh(3,i)-xyzh(3,irhomax))**2)
       if (radi < rcut) then
          xyzmh_ptmass(4,n) = xyzmh_ptmass(4,n) + massoftype(igas)
          npartoftype(igas) = npartoftype(igas) - 1
          call kill_particle(i)
       endif
    enddo

    call shuffle_part(npart)

    !--------------------------------------------------------------------
 case(4)  ! manually create a sink in the core
    !--------------------------------------------------------------------
    nptmass = nptmass + 1
    if (nptmass > maxptmass) call fatal('ptmass_create','nptmass > maxptmass')
    n = nptmass
    xyzmh_ptmass(:,n)      = 0.  ! zero all quantities by default
    xyzmh_ptmass(4,n)      = mcut
    xyzmh_ptmass(ihsoft,n) = hsoft_core
    vxyz_ptmass(:,n)       = 0.

    !--------------------------------------------------------------------
 case(5)  ! set up a triple system
    !--------------------------------------------------------------------
    !resets to (0,0,0) position and velocity of centre of mass for whole system before creating the binary
    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    !removes the dead or accreted particles for a correct total mass computation
    call delete_dead_or_accreted_particles(npart,npartoftype)
    print*,' Got ',npart,npartoftype(igas),' after deleting accreted particles'

    !sets up the binary system orbital parameters
    if (nptmass > 0) then
       mcore = xyzmh_ptmass(4,1)
    else
       mcore = 0.
    endif

    primary_mass = npartoftype(igas) * massoftype(igas) + mcore

    call set_triple(primary_mass,cmass1,cmass2,&
                     a1,a2,hacc1,hacc2,hacc3,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass)

    if (nptmass > 3) then
       xyzmh_ptmass(1:3,1) = xyzmh_ptmass(1:3,2)
       xyzmh_ptmass(:,2) = xyzmh_ptmass(:,3)
       xyzmh_ptmass(:,3) = xyzmh_ptmass(:,4)

       vxyz_ptmass(:,1) = vxyz_ptmass(:,2)
       vxyz_ptmass(:,2) = vxyz_ptmass(:,3)
       vxyz_ptmass(:,3) = vxyz_ptmass(:,4)
    endif

    xyzmh_ptmass(ihsoft,1) = hsoft_primary
    xyzmh_ptmass(ihsoft,2) = hsoft_secondary
    xyzmh_ptmass(ihsoft,3) = hsoft_tertiary

    !shifts gas to the primary point mass created in 'set_binary'
    do i=1,npart
       xyzh(1:3,i) = xyzh(1:3,i) + xyzmh_ptmass(1:3,1)
       vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz_ptmass(1:3,1)
    enddo

    !deletes third point mass
    nptmass = 3

    !resets to (0,0,0) position and velocity of centre of mass for whole system after creating the binary
    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    !--------------------------------------------------------------------
 case(6)  ! set up star for relaxation in corotating frame
    !--------------------------------------------------------------------
    iexternalforce = iext_corotate
    companion_mass_ext = companion_mass

    if (nptmass == 0) then ! Primary core already replaced with potential
       print*,'No sinks in dump. Using primary core properties from input file'
       call open_db_from_file(db,infile_name,20,ierr)
       call read_inopt(icompanion_grav,'icompanion_grav',db) ! Must have icompanion_grav = 2
       call read_inopt(primarycore_mass,'primarycore_mass',db)
       call read_inopt(primarycore_hsoft,'primarycore_hsoft',db)
       call close_db(db)
       print*,'Primary core mass is ',primarycore_mass,' Msun'
       print*,'Primary core softening length is ',primarycore_hsoft,' Rsun'
    elseif (nptmass == 1) then
       print*,'One sink found in dump. Assuming to be primary core.'
       print*,'Primary core mass is ',xyzmh_ptmass(4,1),' Msun'
       print*,'Primary core softening length is ',xyzmh_ptmass(ihsoft,1),' Rsun'
       primarycore_mass  = xyzmh_ptmass(4,1)
       primarycore_hsoft = xyzmh_ptmass(ihsoft,1)
       if (iprimary_grav) then
          icompanion_grav = 2
          nptmass = nptmass - 1
       else
          icompanion_grav = 1
       endif
    endif

    ! Centre to new CoM with the companion
    mass_donor = npartoftype(igas)*massoftype(igas) + primarycore_mass
    omega_corotate = sqrt((mass_donor + companion_mass)/separation**3)
    newCoM = companion_mass / (mass_donor + companion_mass) * separation
    companion_xpos = separation - newCoM
    if (icompanion_grav == 1) xyzmh_ptmass(1,1) = xyzmh_ptmass(1,1) - newCoM
    if (icompanion_grav == 2) then
       primarycore_xpos_old = primarycore_xpos
       primarycore_xpos = -newCoM
    endif
    do i = 1,npart
       ! Zero all particle velocities in the corotating frame, implying that the star is
       ! instantaneously spun up to the orbital frequency.
       vxyzu(1:3,i) = 0.
       ! Move star to new primary position
       xyzh(1,i) = xyzh(1,i) - newCoM
    enddo
    if (icompanion_grav == 1) vxyz_ptmass(1:3,1) = 0.

    ! Calculate softening length, hsoft, of companion gravity. Take hsoft to be 10% of
    ! the companion Roche radius, evaluated with Eggleton (1983)
    mass_ratio = companion_mass / mass_donor
    hsoft = 0.1 * 0.49 * mass_ratio**(2./3.) / (0.6*mass_ratio**(2./3.) + &
            log( 1 + mass_ratio**(1./3.) ) ) * separation
    print*,'Angular velocity of the corotating frame in code units is ',omega_corotate
    print*,'Orbital period is ',2*pi/omega_corotate * utime / 3.15E+07,' years'
    print*,'Softening radius of companion gravity is ',hsoft,' Rsun'

    !--------------------------------------------------------------------
 case(7)  ! set up binary after relaxation in corotating frame
    !--------------------------------------------------------------------
    call open_db_from_file(db,infile_name,20,ierr)
    call read_inopt(icompanion_grav,'icompanion_grav',db)
    call read_inopt(iexternalforce,'iexternalforce',db)
    call read_inopt(omega_corotate,'omega_corotate',db)
    call read_inopt(companion_mass,'companion_mass',db)
    call read_inopt(companion_xpos,'companion_xpos',db)
    if (icompanion_grav == 2) then
       call read_inopt(primarycore_mass,'primarycore_mass',db)
       call read_inopt(primarycore_xpos,'primarycore_xpos',db)
       call read_inopt(primarycore_hsoft,'primarycore_hsoft',db)
    elseif (icompanion_grav == 1) then
       primarycore_mass = xyzmh_ptmass(4,1)
    else
       call fatal('companion_gravity','icompanion_grav not equal to 1 or 2')
    endif
    call close_db(db)

    m1 = npartoftype(igas)*massoftype(igas) + primarycore_mass
    m2 = companion_mass
    a = abs(companion_xpos) + abs(primarycore_xpos)
    print*,' Primary mass from existing file = ',m1,' Msun'
    print*,' Secondary mass from existing file = ',companion_mass,' Msun'
    print*,' Mass of primary core = ',primarycore_mass,' Msun'
    print*,' Softening length of primary core = ',primarycore_hsoft,' Rsun'
    print*,' Orbital separation = ',a,' Rsun'
    if (icompanion_grav == 2) hacc1 = xyzmh_ptmass(ihacc,1)

    ! Add companion particle
    nptmass = nptmass + 1
    xyzmh_ptmass(1:3,2) = (/ companion_xpos, 0., 0. /)
    xyzmh_ptmass(4,2) = m2
    xyzmh_ptmass(ihsoft,2) = hsoft_sec
    xyzmh_ptmass(ihacc,2) = hacc_sec
    vxyz_ptmass(1:3,2) = 0.

    ! Add primary core
    if (icompanion_grav == 2) then
       nptmass = nptmass + 1
       xyzmh_ptmass(1:3,1) = (/ primarycore_xpos, 0., 0. /)
       xyzmh_ptmass(4,1) = primarycore_mass
       xyzmh_ptmass(ihsoft,1) = primarycore_hsoft
       vxyz_ptmass(1:3,1) = 0.
    endif

    call transform_from_corotating_to_inertial_frame(xyzh,vxyzu,npart,nptmass,&
         omega_corotate,xyzmh_ptmass,vxyz_ptmass)
    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    ! Set tmax and dtmax
    period = 2.*pi*sqrt(a**3/(m1 + m2))
    print*,' Orbital period = ',period
    tmax = 30.*period
    dtmax = 0.1*period

    !--------------------------------------------------------------------
 case(9)  ! remove a sink from the simulation (requires nptmass == 3)
    !--------------------------------------------------------------------
    do i=1,nptmass
       write(*,'(A,I2,A,ES10.3,A,ES10.3)') 'Point mass ',i,': M = ',xyzmh_ptmass(4,i),&
                                          ' and radial position = ',sqrt(dot_product(xyzmh_ptmass(1:3,i),xyzmh_ptmass(1:3,i)))
    enddo
    if (iremove == 3) then
       xyzmh_ptmass(:,iremove) = 0.
       vxyz_ptmass(:,iremove) = 0.
       nptmass = 2
    elseif (iremove == 2) then
       xyzmh_ptmass(:,2) = xyzmh_ptmass(:,3)
       vxyz_ptmass(:,2) = vxyz_ptmass(:,3)
       nptmass = 2
    endif

    !--------------------------------------------------------------------
 case(10)  ! transform from corotating to inertial frame (requires nptmass == 2)
    !--------------------------------------------------------------------
    call open_db_from_file(db,infile_name,20,ierr)
    call read_inopt(omega_corotate,'omega_corotate',db)
    call close_db(db)
    call transform_from_corotating_to_inertial_frame(xyzh,vxyzu,npart,nptmass,&
          omega_corotate,xyzmh_ptmass,vxyz_ptmass)

    !--------------------------------------------------------------------
 case(11)  ! shift companion position in the corotating frame (requires nptmass == 2)
    !--------------------------------------------------------------------
    sink_dist = sqrt((xyzmh_ptmass(1,1)-xyzmh_ptmass(1,2))**2 &
                   + (xyzmh_ptmass(2,1)-xyzmh_ptmass(2,2))**2 &
                   + (xyzmh_ptmass(3,1)-xyzmh_ptmass(3,2))**2)

    xyzmh_ptmass(1,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1)) - xyzmh_ptmass(1,2))
    xyzmh_ptmass(2,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(2,2)-xyzmh_ptmass(2,1)) - xyzmh_ptmass(2,2))
    xyzmh_ptmass(3,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(3,2)-xyzmh_ptmass(3,1)) - xyzmh_ptmass(3,2))

    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
    iexternalforce = iext_corotate
    omega_corotate = sqrt((sink_dist-comp_shift)* &
          (xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2)))/(sink_dist-comp_shift)**2

    do i=1,npart
       vxyzu(1:3,i) = 0.0
    enddo
    do i=1,nptmass
       vxyz_ptmass(1:3,i) = 0.0
    enddo

    !--------------------------------------------------------------------
 case(12)  ! add velocity to companion (requires nptmass == 2)
    !--------------------------------------------------------------------
    sink_dist = sqrt((xyzmh_ptmass(1,1)-xyzmh_ptmass(1,2))**2 &
              + (xyzmh_ptmass(2,1)-xyzmh_ptmass(2,2))**2 &
              + (xyzmh_ptmass(3,1)-xyzmh_ptmass(3,2))**2)

    vxyz_ptmass(1,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1)) - vxyz_ptmass(1,2))
    vxyz_ptmass(2,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(2,2)-xyzmh_ptmass(2,1)) - vxyz_ptmass(2,2))
    vxyz_ptmass(3,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(3,2)-xyzmh_ptmass(3,1)) - vxyz_ptmass(3,2))

    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    !--------------------------------------------------------------------
 case(13)  ! (re)set sink properties (requires nptmass >= 1)
    !--------------------------------------------------------------------
    call reset_sink_property(xyzmh_ptmass)

 end select

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  interactively set the moddump parameters
!+
!----------------------------------------------------------------
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt
 use part,      only:nptmass,xyzmh_ptmass,xyzmh_ptmass_label
 use dim,       only:nsinkproperties
 integer :: i,j

 print "(13(/,a))",'Operations:', &
    '  1) add a sink companion (binary)             [nptmass <= 1]', &
    '  2) add a magnetic field in the star          [nptmass <= 1]', &
    '  3) cut profile to create a sink in the core  [nptmass <= 1]', &
    '  4) manually create a sink in the core        [nptmass <= 1]', &
    '  5) set up a triple system                    [nptmass <= 1]', &
    '  6) relax in corotating frame w/ companion    [nptmass <= 1]', &
    '  7) binary after relaxation in corotating fr. [nptmass <= 1]', &
    '  8) add a second star from another dumpfile   [nptmass <= 1]', &
    '  9) remove a sink                             [nptmass == 3]', &
    ' 10/11/12) corotating->inertial / shift / add velocity [nptmass == 2]', &
    ' 13) (re)set sink properties                   [nptmass >= 1]'
 call prompt('Choose an operation ',operation,1,13)

 select case(operation)
 case(1,8)
    mcomp = 0.6
    a1 = 100.
    ecc = 0.
    call prompt('Enter companion mass in code units',mcomp, 0.)
    call prompt('Enter orbit semi-major axis in code units', a1, 0.)
    call prompt('Enter orbit eccentricity', ecc, 0., 1.)
    if (operation == 1) then
       call prompt('Enter accretion radius for the companion in code units', hacc, 0.)
       call prompt('Enter softening length for companion', companion_hsoft, 0.)
    endif
    call prompt('Do you want to transform to a corotating frame and simulate corotating binary?', use_corotating_frame)
    gwinspiral = .false.
    if (.not.use_corotating_frame) then
       call prompt('Do you want to add gravitational radiation reaction?',gwinspiral)
    endif
    if (operation == 8) then
       call prompt('Enter name of second dumpfile',second_dumpfile)
       call prompt('Enter no. of particles in second dumpfile (0 = same as star 1)',nstar2)
    endif
 case(3)
    call prompt('Enter filename of the input stellar profile', densityfile)
    call prompt('Enter mass of the created point mass core', mcut)
    call prompt('Enter softening length of the point mass', hsoft_core)
 case(4)
    call prompt('Enter mass of the created point mass core', mcut)
    call prompt('Enter softening length of the point mass', hsoft_core)
 case(5)
    cmass1 = 0.0095; cmass2 = 0.0095
    a1 = 166.5; a2 = 336.8
    call prompt('Enter 1st companion mass in code units',cmass1,0.)
    call prompt('Enter 2nd companion mass in code units',cmass2,0.)
    call prompt('Enter 1st companion orbit semi-major axis in code units', a1, 0.0)
    call prompt('Enter 2nd companion orbit semi-major axis in code units', a2, 0.0)
    call prompt('Enter accretion radius for the primary in code units', hacc1, 0.0)
    call prompt('Enter accretion radius for the 1st companion in code units', hacc2, 0.0)
    call prompt('Enter accretion radius for the 2nd companion in code units', hacc3, 0.0)
    call prompt('Enter softening length for primary',hsoft_primary,0.)
    call prompt('Enter softening length for secondary',hsoft_secondary,0.)
    call prompt('Enter softening length for tertiary',hsoft_tertiary,0.)
 case(6)
    companion_mass = 1.26
    separation = 865.24
    call prompt('Enter companion mass in Msun',companion_mass,0.)
    call prompt('Enter orbital separation in Rsun',separation,0.)
    call prompt('Name of the input .in file (used if no sink in dump)',infile_name)
    call prompt('Replace primary core with fixed gravitational potential? (if 1 sink present)',iprimary_grav)
 case(7)
    call prompt('Name of the input .in file',infile_name)
    call prompt('Enter eccentricity ',ecc,0.)
    call prompt('Enter accretion radius of secondary in Rsun: ',hacc_sec,0.)
    call prompt('Enter softening length of secondary in Rsun: ',hsoft_sec,0.)
 case(9)
    call prompt('Which sink would you like to remove (2 or 3) : ',iremove,2,3)
 case(10)
    call prompt('Name of the input .in file',infile_name)
 case(11)
    call prompt('How many code units to shift companion (+ve is towards primary)?',comp_shift)
 case(12)
    call prompt('Give velocity to add in direction of the primary : ',vel_shift, 0.0)
 case(13)
    if (nptmass >= 1) then
       do i = 1,nptmass
          print '("sink properties for #",i2," (in code units)")',i
          do j = 1,nsinkproperties
             print "(3x,i2,1x,a,es10.3)", j,xyzmh_ptmass_label(j),xyzmh_ptmass(j,i)
          enddo
       enddo
       call prompt('Select sink particle : ',iselect,1,nptmass)
    endif
    call prompt('Select sink property (0 to exit): ',iproperty,0,nsinkproperties)
    if (iproperty > 0) call prompt('What value (in solar units) ?',propval)
 end select

end subroutine read_interactive_moddumpfile

!----------------------------------------------------------------
!+
!  write options to .moddump file
!+
!----------------------------------------------------------------
subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(a)") '# operation: 1=add companion, 2=Bfield, 3/4=create core sink, 5=triple,'
 write(iunit,"(a)") '#            6=relax corotating, 7=binary after relax, 8=binary from dump,'
 write(iunit,"(a)") '#            9=remove sink, 10=corotating->inertial, 11=shift companion,'
 write(iunit,"(a)") '#            12=add velocity, 13=(re)set sink properties'
 call write_inopt(operation,'operation','operation to perform (1-13)',iunit)

 select case(operation)
 case(1,8)
    call write_inopt(mcomp,'mcomp','companion mass [code units]',iunit)
    call write_inopt(a1,'a1','semi-major axis [code units]',iunit)
    call write_inopt(ecc,'ecc','orbital eccentricity',iunit)
    if (operation == 1) then
       call write_inopt(hacc,'hacc','accretion radius for the companion [code units]',iunit)
       call write_inopt(companion_hsoft,'companion_hsoft','softening length for companion [code units]',iunit)
    endif
    call write_inopt(use_corotating_frame,'use_corotating_frame','transform to a corotating frame',iunit)
    if (.not.use_corotating_frame) then
       call write_inopt(gwinspiral,'gwinspiral','add gravitational radiation reaction',iunit)
    endif
    if (operation == 8) then
       call write_inopt(second_dumpfile,'second_dumpfile','name of the second dumpfile',iunit)
       call write_inopt(nstar2,'nstar2','number of particles in the second dumpfile (0 = same as star 1)',iunit)
    endif
 case(3)
    call write_inopt(densityfile,'densityfile','filename of the input stellar profile',iunit)
    call write_inopt(mcut,'mcut','mass of the created point-mass core [code units]',iunit)
    call write_inopt(hsoft_core,'hsoft_core','softening length of the created point-mass core',iunit)
 case(4)
    call write_inopt(mcut,'mcut','mass of the created point-mass core [code units]',iunit)
    call write_inopt(hsoft_core,'hsoft_core','softening length of the created point-mass core',iunit)
 case(5)
    call write_inopt(cmass1,'cmass1','1st companion mass (triple) [code units]',iunit)
    call write_inopt(cmass2,'cmass2','2nd companion mass (triple) [code units]',iunit)
    call write_inopt(a1,'a1','1st companion semi-major axis (triple) [code units]',iunit)
    call write_inopt(a2,'a2','2nd companion semi-major axis (triple) [code units]',iunit)
    call write_inopt(hacc1,'hacc1','accretion radius for the primary [code units]',iunit)
    call write_inopt(hacc2,'hacc2','accretion radius for the 1st companion [code units]',iunit)
    call write_inopt(hacc3,'hacc3','accretion radius for the 2nd companion [code units]',iunit)
    call write_inopt(hsoft_primary,'hsoft_primary','softening length for primary [code units]',iunit)
    call write_inopt(hsoft_secondary,'hsoft_secondary','softening length for secondary [code units]',iunit)
    call write_inopt(hsoft_tertiary,'hsoft_tertiary','softening length for tertiary [code units]',iunit)
 case(6)
    call write_inopt(companion_mass,'companion_mass','companion mass [Msun]',iunit)
    call write_inopt(separation,'separation','orbital separation [Rsun]',iunit)
    call write_inopt(infile_name,'infile_name','name of the .in file (used if no sink in dump)',iunit)
    call write_inopt(iprimary_grav,'iprimary_grav','replace primary core with a fixed potential (if 1 sink)',iunit)
 case(7)
    call write_inopt(infile_name,'infile_name','name of the .in file to read companion info from',iunit)
    call write_inopt(ecc,'ecc','orbital eccentricity',iunit)
    call write_inopt(hacc_sec,'hacc_sec','accretion radius of secondary [Rsun]',iunit)
    call write_inopt(hsoft_sec,'hsoft_sec','softening length of secondary [Rsun]',iunit)
 case(9)
    call write_inopt(iremove,'iremove','which sink to remove (2 or 3)',iunit)
 case(10)
    call write_inopt(infile_name,'infile_name','name of the .in file to read omega_corotate from',iunit)
 case(11)
    call write_inopt(comp_shift,'comp_shift','code units to shift companion (+ve towards primary)',iunit)
 case(12)
    call write_inopt(vel_shift,'vel_shift','velocity to add in the direction of the primary [code units]',iunit)
 case(13)
    call write_inopt(iselect,'iselect','which sink particle to (re)set properties for',iunit)
    call write_inopt(iproperty,'iproperty','sink property index to (re)set (0 = none)',iunit)
    call write_inopt(propval,'propval','new value for the selected sink property [solar units]',iunit)
 end select

 close(iunit)

end subroutine write_moddump

!----------------------------------------------------------------
!+
!  read options from .moddump file
!+
!----------------------------------------------------------------
subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use dim,          only:nsinkproperties
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(operation,'operation',db,min=1,max=13,errcount=nerr)

 select case(operation)
 case(1,8)
    call read_inopt(mcomp,'mcomp',db,min=0.,errcount=nerr)
    call read_inopt(a1,'a1',db,min=0.,errcount=nerr)
    call read_inopt(ecc,'ecc',db,min=0.,max=1.,errcount=nerr)
    if (operation == 1) then
       call read_inopt(hacc,'hacc',db,min=0.,errcount=nerr)
       call read_inopt(companion_hsoft,'companion_hsoft',db,min=0.,errcount=nerr)
    endif
    call read_inopt(use_corotating_frame,'use_corotating_frame',db,errcount=nerr)
    if (.not.use_corotating_frame) then
       call read_inopt(gwinspiral,'gwinspiral',db,errcount=nerr,default=.false.)
    else
       gwinspiral = .false.
    endif
    if (operation == 8) then
       call read_inopt(second_dumpfile,'second_dumpfile',db,errcount=nerr)
       call read_inopt(nstar2,'nstar2',db,min=0,errcount=nerr)
    endif
 case(3)
    call read_inopt(densityfile,'densityfile',db,errcount=nerr)
    call read_inopt(mcut,'mcut',db,min=0.,errcount=nerr)
    call read_inopt(hsoft_core,'hsoft_core',db,min=0.,errcount=nerr)
 case(4)
    call read_inopt(mcut,'mcut',db,min=0.,errcount=nerr)
    call read_inopt(hsoft_core,'hsoft_core',db,min=0.,errcount=nerr)
 case(5)
    call read_inopt(cmass1,'cmass1',db,min=0.,errcount=nerr)
    call read_inopt(cmass2,'cmass2',db,min=0.,errcount=nerr)
    call read_inopt(a1,'a1',db,min=0.,errcount=nerr)
    call read_inopt(a2,'a2',db,min=0.,errcount=nerr)
    call read_inopt(hacc1,'hacc1',db,min=0.,errcount=nerr)
    call read_inopt(hacc2,'hacc2',db,min=0.,errcount=nerr)
    call read_inopt(hacc3,'hacc3',db,min=0.,errcount=nerr)
    call read_inopt(hsoft_primary,'hsoft_primary',db,min=0.,errcount=nerr)
    call read_inopt(hsoft_secondary,'hsoft_secondary',db,min=0.,errcount=nerr)
    call read_inopt(hsoft_tertiary,'hsoft_tertiary',db,min=0.,errcount=nerr)
 case(6)
    call read_inopt(companion_mass,'companion_mass',db,min=0.,errcount=nerr)
    call read_inopt(separation,'separation',db,min=0.,errcount=nerr)
    call read_inopt(infile_name,'infile_name',db,errcount=nerr)
    call read_inopt(iprimary_grav,'iprimary_grav',db,errcount=nerr)
 case(7)
    call read_inopt(infile_name,'infile_name',db,errcount=nerr)
    call read_inopt(ecc,'ecc',db,min=0.,errcount=nerr)
    call read_inopt(hacc_sec,'hacc_sec',db,min=0.,errcount=nerr)
    call read_inopt(hsoft_sec,'hsoft_sec',db,min=0.,errcount=nerr)
 case(9)
    call read_inopt(iremove,'iremove',db,min=2,max=3,errcount=nerr)
 case(10)
    call read_inopt(infile_name,'infile_name',db,errcount=nerr)
 case(11)
    call read_inopt(comp_shift,'comp_shift',db,errcount=nerr)
 case(12)
    call read_inopt(vel_shift,'vel_shift',db,errcount=nerr)
 case(13)
    call read_inopt(iselect,'iselect',db,min=1,errcount=nerr)
    call read_inopt(iproperty,'iproperty',db,min=0,max=nsinkproperties,errcount=nerr)
    call read_inopt(propval,'propval',db,errcount=nerr)
 end select

 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

!----------------------------------------------------------------
!+
!  (re)set a single sink property, converting from solar units
!+
!----------------------------------------------------------------
subroutine reset_sink_property(xyzmh_ptmass)
 use part,    only:nptmass,ihacc,ihsoft,imacc,ilum,ireff,xyzmh_ptmass_label
 use units,   only:umass,udist,utime,unit_energ
 use physcon, only:solarm,solarr,solarl
 use dim,     only:nsinkproperties
 use io,      only:iprint
 real, intent(inout) :: xyzmh_ptmass(:,:)
 integer :: j
 real    :: fac

 do j = 1,nsinkproperties
    write(iprint,"(3x,i2,1x,a,es10.3)")  j,xyzmh_ptmass_label(j),xyzmh_ptmass(j,iselect)
 enddo

 if (iproperty < 1 .or. iproperty > nsinkproperties) return
 if (iselect < 1 .or. iselect > nptmass) stop 'wrong sink particle number'

 select case (iproperty)
 case (ihacc,ihsoft,iReff)
    fac =  solarr / udist
 case (ilum)
    fac =  solarl * utime / unit_energ
 case (imacc,4)
    fac = solarm / umass
 case default
    fac = 1.
 end select
 xyzmh_ptmass(iproperty,iselect) = propval*fac

 print *,'summary'
 do j = 1,nsinkproperties
    write(iprint,"(3x,i2,1x,a,es10.3)")  j,xyzmh_ptmass_label(j),xyzmh_ptmass(j,iselect)
 enddo

end subroutine reset_sink_property

subroutine transform_from_corotating_to_inertial_frame(xyzh,vxyzu,npart,nptmass,&
           omega_corotate,xyzmh_ptmass,vxyz_ptmass)
 use options,     only:iexternalforce
 use vectorutils, only:cross_product3D
 integer, intent(in)    :: npart,nptmass
 real,    intent(in)    :: omega_corotate,xyzh(:,:),xyzmh_ptmass(:,:)
 real,    intent(inout) :: vxyzu(:,:),vxyz_ptmass(:,:)
 real :: omega_vec(3),omegacrossr(3)
 integer :: i

 iexternalforce = 0
 omega_vec = (/ 0.,0.,omega_corotate /)
 do i=1,npart
    call cross_product3D(omega_vec,xyzh(1:3,i),omegacrossr)
    vxyzu(1,i) = vxyzu(1,i) + omegacrossr(1)
    vxyzu(2,i) = vxyzu(2,i) + omegacrossr(2)
    vxyzu(3,i) = vxyzu(3,i) + omegacrossr(3)
 enddo
 do i=1,nptmass
    call cross_product3D(omega_vec,xyzmh_ptmass(1:3,i),omegacrossr)
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + omegacrossr(1)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + omegacrossr(2)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + omegacrossr(3)
 enddo

end subroutine transform_from_corotating_to_inertial_frame

subroutine set_triple(mprimary,msecondary,mtertiary,semimajoraxis12,semimajoraxis13,&
                      accretion_radius1,accretion_radius2,accretion_radius3,&
                      xyzmh_ptmass,vxyz_ptmass,nptmass)
 real,    intent(in)    :: mprimary,msecondary,mtertiary
 real,    intent(in)    :: semimajoraxis12,semimajoraxis13
 real,    intent(in)    :: accretion_radius1,accretion_radius2,accretion_radius3
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass

 integer :: i1,i2,i3
 real    :: m1,m2,m3,mtot,dx12(3),dx13(3),dv12(3),dv13(3)
 real    :: x1(3),x2(3),x3(3),v1(3),v2(3),v3(3)

 i1 = nptmass + 1
 i2 = nptmass + 2
 i3 = nptmass + 3
 nptmass = nptmass + 3

 ! masses
 m1 = mprimary
 m2 = msecondary
 m3 = mtertiary
 mtot = m1 + m2 + m3

 dx12 = (/semimajoraxis12,0.,0./)
 dv12 = (/0.,sqrt((m1+m2)/dx12(1)),0./)

 dx13 = (/semimajoraxis13,0.,0./)
 dv13 = (/0.,sqrt(mtot/dx13(1)),0./)

 ! positions of each star so centre of mass is at zero
 x1 = -(dx12*m2 + dx13*m3)/mtot
 x2 = (dx12*m1 + dx12*m3 - dx13*m3)/mtot
 x3 = (dx13*m1 + dx13*m2 - dx12*m2)/mtot

 ! velocities
 v1 = -(dv12*m2 + dv13*m3)/mtot
 v2 = (dv12*m1 + dv12*m3 - dv13*m3)/mtot
 v3 = (dv13*m1 + dv13*m2 - dv12*m2)/mtot

 ! positions and accretion radii
 xyzmh_ptmass(:,i1:i3) = 0.
 xyzmh_ptmass(1:3,i1) = x1
 xyzmh_ptmass(1:3,i2) = x2
 xyzmh_ptmass(1:3,i3) = x3
 xyzmh_ptmass(4,i1) = m1
 xyzmh_ptmass(4,i2) = m2
 xyzmh_ptmass(4,i3) = m3
 xyzmh_ptmass(5,i1) = accretion_radius1
 xyzmh_ptmass(5,i2) = accretion_radius2
 xyzmh_ptmass(5,i3) = accretion_radius3
 xyzmh_ptmass(6,i1) = 0.0
 xyzmh_ptmass(6,i2) = 0.0
 xyzmh_ptmass(6,i3) = 0.0

 ! velocities
 vxyz_ptmass(:,i1) = v1
 vxyz_ptmass(:,i2) = v2
 vxyz_ptmass(:,i3) = v3

end subroutine set_triple

end module moddump
