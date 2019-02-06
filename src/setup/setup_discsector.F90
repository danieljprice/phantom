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
!  this module does setups of a disc sector for local simulations
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    HoverR      -- disc aspect ratio at R=Rsect_in
!    R_in        -- inner total disc radius
!    R_out       -- outer total disc radius
!    Rsect_in    -- sector inner radius
!    Rsect_out   -- sector outer radius
!    accradius   --  accretion radius
!    disc_mass   -- Total Disc mass
!    dr_bound    -- Radial boundary thickness
!    npart       -- number of particles
!    object_mass -- object mass in code units
!    p_index     -- surface density powerlaw index
!    phi_inject  -- azimuthal range of injection zone
!    phimax      -- azimuthal extent (-phimax,phimax)
!    q_index     -- sound speed powerlaw index
!    udist       -- distance unit in cm
!    umass       -- mass unit in g
!
!  DEPENDENCIES: dim, extern_corotate, externalforces, infile_utils,
!    inject, io, options, part, physcon, prompting, setdisc, setup_params,
!    units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

!--private module variables
 integer :: np
 real    :: R_in,R_out,R_mid,Rsect_in,Rsect_out,phimax
 real    :: p_index,q_index,HoverR,disc_mass,dr_bound,phi_inject
 real    :: object_mass,accradius

 private
 real(kind=8) :: udist,umass

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,         only:set_disc
 use io,              only:master
 use infile_utils,    only:read_next_inopt
 use externalforces,  only:accradius1,iext_star,iext_corotate
 use extern_corotate, only:omega_corotate
#ifdef INJECT_PARTICLES
 use inject,          only:set_injection_parameters
#endif
 use options,         only:iexternalforce,icooling
 use part,            only:set_particle_type,iboundary
 use physcon,         only:au,solarm,pi
 use prompting,       only:prompt
 use units,           only:set_units
 use physcon,         only:au,solarm
 use setup_params,    only:rhozero
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 integer :: ierr, i, nboundary
 character(len=100) :: filename
 logical :: iexist
 real ::  phipart,sig0,sig_in,rpart,omega0,vmag,phimaxrad
 real :: v_0(3),v_subtract(3)
 real :: annulus_halfwidth, midsep

 filename=trim(fileprefix)//".setup"

 umass = solarm
 udist = au

 print "(/,65('-'),2(/,a),/,65('-'),/)",&
  ' Welcome to the Disc Sector Setup Routine^TM', &
  ' Brought to you by Duncan Forgan and Daniel Price'

 inquire(file=filename,exist=iexist)
 if (iexist) then

    ! Read disc parameters from file
    call read_setupfile(filename,ierr)
    call set_units(dist=udist,mass=umass,G=1.)

    if(ierr/=0) then
       if(id==master) call write_setupfile(filename)
       stop
    endif

 elseif (id==master) then

    print "(a,/)",trim(filename)//' not found: using interactive setup'

    ! Set units
    call prompt('Enter code units of mass in g ',umass,0.)
    call prompt('Enter code units of distance in cm ',udist,0.)

    call set_units(dist=au,mass=solarm,G=1.)

    ! Read in options (with preset defaults)

    np = 100000
    call prompt('How many SPH particles? ',np, 1000)

    R_in = 1.0
    call prompt('Enter inner (total) disc radius (code units)',R_in,0.)


    R_out = 25.0
    call prompt('Enter outer (total) disc radius (code units)',R_out,R_in)

    Rsect_in = 10.0
    call prompt("Enter disc sector's inner radius (code units)",Rsect_in,0.)

    Rsect_out = 15.0
    call prompt("Enter disc sector's outer radius (code units)",Rsect_out,Rsect_in,R_out )

    dr_bound = 0.1
    call prompt("How thick are the radial pressure boundaries?", dr_bound, 0.)

    phimax = 35.0
    call prompt('Enter azimuthal range of sector (-phimax,phimax)', phimax, 0.,45.0)

    p_index = 1.0
    call prompt('Enter the surface density powerlaw index:',p_index)

    q_index = 0.75
    call prompt('Enter the sound speed powerlaw index:',q_index)

    HoverR = 0.1
    call prompt('Enter the disc aspect ratio (H/R), measured at R=R_in: ',HoverR, 0.)

    disc_mass = 0.25
    call prompt('What is the total disc mass? ', disc_mass,0.)

    object_mass = 1.0
    call prompt('What is the central object mass?',object_mass,0.)

    accradius = 0.1*R_in
    call prompt("What is the accretion radius? ",accradius,0.0,0.5*R_in)



    phi_inject = 0.1*phimax

    !
    !--Write these inputs to a new parameters file
    !
    call write_setupfile(filename)

    print "(a)", '>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
    stop
 else
    stop
 endif

 phimaxrad = phimax*(pi/180.0)

 npartoftype(:) = 0
 npartoftype(1) = np
 npart = np
 if (size(vxyzu(:,1)) >= 4) then
    gamma = 5./3.
 else
    gamma = 1.0
 endif
 time    = 0.

 !--calculate sigma0 of entire disc
 sig0 = sigma0(disc_mass,R_in,R_out,p_index)

 !--set sig_in as required for set_disc (sig_norm at R=Rin)
 sig_in = sig0*R_in**p_index

 call set_disc(id,master     = master,             &
               npart         = npartoftype(1),     &
               rmin          = Rsect_in-dr_bound,  &
               rmax          = Rsect_out+dr_bound, &
               phimin        = -phimaxrad,         &
               phimax        = phimaxrad,          &
               p_index       = p_index,            &
               q_index       = q_index,            &
               HoverR        = HoverR,             &
               sig_norm      = sig_in,             &
               star_mass     = object_mass,        &
               gamma         = gamma,              &
               particle_mass = massoftype(1),      &
               hfact         = hfact,              &
               xyzh          = xyzh,               &
               vxyzu         = vxyzu,              &
               polyk         = polyk,              &
               prefix        = fileprefix)

 icooling = 1 ! Switches on beta cooling

 !
 !--Change to corotating frame

 ! Find middle of disc sector, calculate its velocity
 R_mid = (Rsect_out - Rsect_in)/2.0 + Rsect_in
 annulus_halfwidth = 0.5*(Rsect_out-Rsect_in)
 vmag = sqrt(object_mass/R_mid)
 omega0 = sqrt(object_mass/(R_mid)**3)

 ! v_phi = v_y at y=0
 ! Obtain the true v_phi at any point (r,phi) via rotation in z axis

 v_0 = (/0.0, vmag,0.0/)

 !print *, 'Transforming to corotating frame: angular velocity ', omega0

 do i=1,npart
    rpart = sqrt(xyzh(1,i)*xyzh(1,i)  + xyzh(2,i)*xyzh(2,i))
    phipart = atan2(xyzh(2,i),xyzh(1,i))
    call rotate_z(v_0, v_subtract,phipart)
    vxyzu(1:3,i) = vxyzu(1:3,i)-v_subtract(:)
 enddo

 !
 !--Define a typical density so that B-fields can be correctly initialised
 !--(here we use the midplane density at the inner edge of the sector)
 !

 rhozero = sig0/(2.0*HoverR*Rsect_in)

 ! Initialise external forces at the origin

 print "(a)", ' Central object represented by external force at the system origin'
 print*, ' Accretion Radius: ', accradius

 iexternalforce=iext_corotate
 omega_corotate = omega0
 accradius1=accradius

 ! Label boundary particles
 nboundary = 0
 do i=1,npart
    rpart = sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i))
    midsep = abs(rpart - R_mid)
    if(midsep > annulus_halfwidth) then
       call set_particle_type(i,iboundary)
       nboundary = nboundary+1
    endif
 enddo

 print "(a,i10,a,i2)",' Boundary particles initialised: total ', nboundary,' of type ', iboundary

#ifdef INJECT_PARTICLES
! call function to set injection parameters
 call set_injection_parameters(R_in, R_out, Rsect_in,Rsect_out,dr_bound,&
                               phimax,phi_inject, p_index,q_index,HoverR,&
                               disc_mass,object_mass)
#endif

 return
end subroutine setpart

subroutine read_setupfile(filename, ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use dim, only: maxp
 character(len=*), intent(in) :: filename
 integer, intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

!character (len=20) :: name
!character (len=120) :: valstring
!character (len=100) :: filename

 print "(a)", 'reading setup options from '//trim(filename)

 nerr = 0

 call open_db_from_file(db,filename, iunit,ierr)
 call read_inopt(np,'npart',db,min = 0, max=maxp, errcount=nerr)
 call read_inopt(umass,'umass',db,min = 0., errcount =nerr)
 call read_inopt(udist,'udist',db,min = 0., errcount =nerr)
 call read_inopt(R_in,'R_in',db,min = 0., errcount =nerr)
 call read_inopt(R_out,'R_out',db,min = 0., errcount =nerr)
 call read_inopt(Rsect_in,'Rsect_in',db,min = 0., errcount =nerr)
 call read_inopt(Rsect_out,'Rsect_out',db,min = 0., errcount =nerr)
 call read_inopt(dr_bound, 'dr_bound', db,min=0., errcount=nerr)
 call read_inopt(phimax, 'phimax',db,min=0., errcount=nerr)
 call read_inopt(phi_inject, 'phi_inject', db,min=0., errcount=nerr)
 call read_inopt(p_index,'p_index',db, errcount =nerr)
 call read_inopt(q_index,'q_index',db, errcount =nerr)
 call read_inopt(HoverR,'HoverR',db,min = 0., errcount =nerr)
 call read_inopt(disc_mass,'disc_mass',db,min = 0., errcount = nerr)
 call read_inopt(object_mass,'object_mass',db,min = 0., errcount = nerr)
 call read_inopt(accradius,'accradius',db,min = 0., errcount=nerr)

 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif


 return

end subroutine read_setupfile

subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)

 open(unit=iunit,file=filename,status='replace',form='formatted')

 write(iunit,"(a)") '# input file for disc setup routines'
 write(iunit,"(a)") '# Generated by setup_discsector.F90'

 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'npart','number of particles',iunit)

 write(iunit,"(/,a)") '# units'
 call write_inopt(udist,'udist','distance unit in cm',iunit)
 call write_inopt(umass,'umass','mass unit in g',iunit)

 write(iunit,"(/,a)") '# disc extent'
 call write_inopt(R_in,'R_in','inner total disc radius',iunit)
 call write_inopt(R_out,'R_out','outer total disc radius',iunit)
 call write_inopt(Rsect_in,'Rsect_in','sector inner radius',iunit)
 call write_inopt(Rsect_out,'Rsect_out','sector outer radius',iunit)
 call write_inopt(dr_bound,'dr_bound','Radial boundary thickness',iunit)

 call write_inopt(phimax, 'phimax', 'azimuthal extent (-phimax,phimax)',iunit)
 call write_inopt(phi_inject, 'phi_inject', 'azimuthal range of injection zone',iunit)

 write(iunit,"(/,a)") '# disc powerlaw indices'
 call write_inopt(p_index,'p_index','surface density powerlaw index',iunit)
 call write_inopt(q_index,'q_index','sound speed powerlaw index',iunit)

 write(iunit,"(/,a)") '# disc axial properties'
 call write_inopt(HoverR,'HoverR','disc aspect ratio at R=Rsect_in',iunit)

 write(iunit,"(/,a)") '# disc and object masses (plus object mass options)'

 call write_inopt(disc_mass,'disc_mass','Total Disc mass',iunit)
 call write_inopt(object_mass,'object_mass','object mass in code units',iunit)

 call write_inopt(accradius, 'accradius',' accretion radius',iunit)


 close(iunit)

end subroutine write_setupfile

! Rotates a vector in the z axis
subroutine rotate_z(oldvec,newvec,phi)
 real, intent(inout) :: oldvec(3), newvec(3)
 real, intent(in) :: phi

 newvec(1) = oldvec(1)*cos(phi) - oldvec(2)*sin(phi)
 newvec(2) = oldvec(1)*sin(phi) + oldvec(2)*cos(phi)

 return
end subroutine rotate_z

! Calculates sigma0 for a given set of disc parameters
real function sigma0(Mdisc, Rinner, Router, p_index)
 real, intent(in) :: Mdisc,Rinner, Router, p_index

 real :: exponent

 sigma0 = Mdisc/(2.0*3.141592654)
 exponent = 2.0-p_index

 if(p_index==2.0) then
    sigma0 = sigma0*log(Rinner/Router)
 else
    sigma0 = sigma0*exponent/(Router**exponent - Rinner**exponent)
 endif

end function sigma0

end module setup
