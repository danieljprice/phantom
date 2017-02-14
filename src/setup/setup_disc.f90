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
!  this module does general accretion disc setups
!  Modified from an original routine by Giuseppe Lodato
!  Modified by Duncan Forgan to allow interactive setup and parameter file reads
!  plus the option to model central object via a sink or via external forces
!
!  REFERENCES: None
!
!  OWNER: Duncan Forgan
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    HoverR        -- disc aspect ratio at R=R_in
!    R_in          -- inner disc radius
!    R_out         -- outer disc radius
!    accradius     -- accretion radius
!    disc_mass     -- disc mass in code units
!    icentralforce -- central mass force choice (0=sink,1=external force)
!    npart         -- number of particles
!    object_mass   -- object mass in code units
!    p_index       -- surface density powerlaw index
!    q_index       -- sound speed powerlaw index
!    udist         -- distance unit in cm
!    umass         -- mass unit in g
!    xinc          -- disc inclination in degrees
!
!  DEPENDENCIES: dim, externalforces, infile_utils, io, options, part,
!    physcon, prompting, setdisc, setup_params, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

!--private module variables
 integer :: np, icentralforce
 real    :: R_in,R_out,xinc
 real    :: p_index,q_index,HoverR,disc_mass
 real    :: object_mass, accradius

 private
real(kind=8) :: udist,umass

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,        only:set_disc
 !use setbinary,     only:set_binary
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use io,             only:master
 use infile_utils,   only:read_next_inopt
 use externalforces, only:accradius1
 use options,        only:iexternalforce,icooling
 use physcon,        only:au,solarm,pi
 use prompting,      only:prompt
 use units,          only:set_units
 use physcon,        only:au,solarm
 use setup_params,   only: rhozero
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 integer                          :: ierr
 character(len=100)               :: filename
 logical                          :: iexist

 filename=trim(fileprefix)//".setup"

 print "(/,65('-'),2(/,a),/,65('-'),/)",&
  ' Welcome to the Disc Setup Routine^TM', &
  ' Brought to you by Duncan Forgan and Daniel Price'

 inquire(file=filename,exist=iexist)
 if (iexist) then

    ! Read disc parameters from file
    call read_discinputfile(filename,ierr)
    call set_units(dist=udist,mass=umass,G=1.)

    if(ierr/=0) then
        if(id==master) call write_discinputfile(filename)
        stop
    endif

 elseif (id==master) then

    print "(a,/)",trim(filename)//' not found: using interactive setup'

    ! Set units
    umass = solarm
    udist = au
    call prompt('Enter code units of mass in g ',umass,0.)
    call prompt('Enter code units of distance in cm ',udist,0.)

    call set_units(dist=au,mass=solarm,G=1.)

    ! Read in options (with preset defaults)

    np = 100000
    call prompt('How many SPH particles? ',np, 1000)

    R_in = 1.0
    call prompt('Enter inner disc radius (code units)',R_in,0.)

    R_out = 25.0
    call prompt('Enter outer disc radius (code units)',R_out,R_in)

    p_index = 1.0
    call prompt('Enter the surface density powerlaw index:',p_index)

    q_index = 0.75
    call prompt('Enter the sound speed powerlaw index:',q_index)

    HoverR = 0.1
    call prompt('Enter the disc aspect ratio (H/R), measured at R=R_in: ',HoverR, 0.)

    xinc = 0.0
    call prompt('What is the disc inclination (in degrees)? ', xinc, 0.0, 180.0)

    disc_mass = 0.1
    call prompt('What is the disc mass in code units?', disc_mass,0.)

    icentralforce=0
    call prompt('How is the central object defined? (0=sink particle, 1=external potential)',icentralforce,0,1)

    object_mass = 1.0
    call prompt('What is the central object mass?',object_mass,0.)

    accradius = R_in

    call prompt("What is the accretion radius? ",accradius,0.25,R_in)

    !
    !--Write these inputs to a new parameters file
    !
    call write_discinputfile(filename)

    print "(a)", '>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
    stop
 else
    stop
 endif

 xinc = xinc*(pi/180.0)

 npartoftype(:) = 0
 npartoftype(1) = np
 npart = np
 if (size(vxyzu(:,1)) >= 4) then
    gamma = 5./3.
 else
    gamma = 1.0
 endif
 time    = 0.

 call set_disc(id,master=master,&
                npart   = npartoftype(1),&
                rmin    = R_in, &
                rmax    = R_out,&
                p_index = p_index,    &
                q_index = q_index,   &
                HoverR  = HoverR,   &
                disc_mass = disc_mass,   &
                star_mass = object_mass,    &
                gamma   = gamma,  &
                particle_mass = massoftype(1), &
                hfact=hfact, &
                xyzh=xyzh, &
                vxyzu=vxyzu, &
                polyk=polyk, &
                inclination=xinc, &
                twist=.false., &
                prefix = fileprefix)


 iexternalforce = icentralforce
 icooling = 1 ! Switches on beta cooling

 !
 !--Define a typical density so that B-fields can be correctly initialised
 !

 rhozero = disc_mass/(2.0*pi*HoverR*R_out*R_out*R_out)

 ! Initialise either a sink or external forces at the origin

 if(iexternalforce==0) then
    print "(a)", ' Central object represented by a sink at the system origin'
    print*, ' Accretion Radius: ', accradius

    nptmass = 1
    xyzmh_ptmass(:,:) = 0.0
    vxyz_ptmass(:,:) = 0.0
    xyzmh_ptmass(4,nptmass) = object_mass
    xyzmh_ptmass(ihacc,nptmass) = accradius
    xyzmh_ptmass(ihsoft,nptmass) = 0.0
 endif

 if(iexternalforce==1) then
    print "(a)", ' Central object represented by external force with accretion boundary'
    print*, ' Accretion Radius: ', accradius

    accradius1 = accradius
 endif

 return
end subroutine setpart
!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_discinputfile(filename, ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use dim,          only:maxp
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer,          parameter   :: iunit = 21
 integer                       :: nerr
 type(inopts),     allocatable :: db(:)

 print "(a)", 'reading setup options from '//trim(filename)

 nerr = 0

 call open_db_from_file(db,filename, iunit,ierr)
 call read_inopt(np,'npart',db,min = 0, max=maxp, errcount=nerr)
 call read_inopt(umass,'umass',db,min = 0., errcount =nerr)
 call read_inopt(udist,'udist',db,min = 0., errcount =nerr)
 call read_inopt(R_in,'R_in',db,min = 0., errcount =nerr)
 call read_inopt(R_out,'R_out',db,min = 0., errcount =nerr)
 call read_inopt(p_index,'p_index',db, errcount =nerr)
 call read_inopt(q_index,'q_index',db, errcount =nerr)
 call read_inopt(HoverR,'HoverR',db,min = 0., errcount =nerr)
 call read_inopt(xinc,'xinc',db,min = 0.,max = 180.0, errcount =nerr)
 call read_inopt(disc_mass,'disc_mass',db,min = 0., errcount = nerr)
 call read_inopt(object_mass,'object_mass',db,min = 0., errcount = nerr)
 call read_inopt(icentralforce,'icentralforce',db,min=0,max=1, errcount=nerr)
 call read_inopt(accradius,'accradius',db,min = 0., errcount=nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

 return
end subroutine read_discinputfile
!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_discinputfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)

 open(unit=iunit,file=filename,status='replace',form='formatted')

 write(iunit,"(a)") '# input file for disc setup routines'
 write(iunit,"(a)") '# Generated by setup_disc.f90'

 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'npart','number of particles',iunit)

 write(iunit,"(/,a)") '# units'
 call write_inopt(udist,'udist','distance unit in cm',iunit)
 call write_inopt(umass,'umass','mass unit in g',iunit)

 write(iunit,"(/,a)") '# disc extent'
 call write_inopt(R_in,'R_in','inner disc radius',iunit)
 call write_inopt(R_out,'R_out','outer disc radius',iunit)

 write(iunit,"(/,a)") '# disc powerlaw indices'
 call write_inopt(p_index,'p_index','surface density powerlaw index',iunit)
 call write_inopt(q_index,'q_index','sound speed powerlaw index',iunit)

 write(iunit,"(/,a)") '# disc axial properties'
 call write_inopt(HoverR,'HoverR','disc aspect ratio at R=R_in',iunit)
 call write_inopt(xinc,'xinc','disc inclination in degrees',iunit)

 write(iunit,"(/,a)") '# disc and object masses (plus object mass options)'

 call write_inopt(disc_mass,'disc_mass','disc mass in code units',iunit)
 call write_inopt(object_mass,'object_mass','object mass in code units',iunit)

 call write_inopt(icentralforce,'icentralforce', 'central mass force choice (0=sink,1=external force)',iunit)
 call write_inopt(accradius, 'accradius','accretion radius',iunit)

close(iunit)

end subroutine write_discinputfile
!----------------------------------------------------------------
end module setup
