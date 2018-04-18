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
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, extern_lensethirring, externalforces, io, options,
!    part, physcon, setdisc, setup_params, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 real,    private :: mhole,mdisc,r_in,r_out,spin,honr,theta,accrad
 integer, private :: np

 private

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,        only:set_disc
 use part,           only:igas
 use io,             only:master
 use externalforces, only:accradius1,accradius1_hard
 use options,        only:iexternalforce,alpha,alphau,iexternalforce
 use units,          only:set_units,umass
 use physcon,        only:solarm,pi
#ifdef GR
 use metric,         only:mass1,a
#else
 use externalforces, only:mass1,iext_einsteinprec
 use extern_lensethirring, only:blackhole_spin
#endif
 use prompting,      only:prompt
 use timestep,       only:tmax,dtmax
 use eos,            only:ieos
 use kernel,         only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=120) :: filename
 integer :: ierr
 logical :: iexist

 time            = 0.
 alphau          = 0.0
 npartoftype(:)  = 0
 iexternalforce  = 1
 hfact           = hfact_default
#ifdef GR
 ieos  = 4
 gamma = 5./3.
#else
 ieos  = 2
 gamma = 1.
 iexternalforce = iext_einsteinprec
#endif

 tmax  = 1000.
 dtmax = 10.

!
! Set default problem parameters
!
 mhole  = 1.e6    ! (solarm)
 mdisc  = 10.     ! (solarm)
 r_in   = 40.     ! (GM/c^2)
 r_out  = 160.    ! (GM/c^2)
 spin   = 0.
 honr   = 0.02
 theta  = 0.      ! inclination angle (degrees)
 np     = 1e5
 accrad = 4.      ! (GM/c^2)

!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Tidal disruption in GR'
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
  if (id==master) then
     call write_setupfile(filename)
     print*,' Edit '//trim(filename)//' and rerun phantomsetup'
  endif
  stop
 endif
 accradius1 = accrad
 npart = np

!
! Convert to code units
!
 mhole = mhole*solarm
 call set_units(G=1.,c=1.,mass=mhole) ! Set central mass to M=1 in code units
 mdisc           = mdisc*solarm/umass
 accradius1_hard = accradius1

!
! Convert to radians
!
 theta = theta/180. * pi

 call set_disc(id,master,&
               npart         = npart,                &
               rmin          = r_in,                 &
               rmax          = r_out,                &
               p_index       = -1.0,                 &
               q_index       = 0.0,                  &
               HoverR        = honr,                 &
               gamma         = gamma,                &
               hfact         = hfact,                &
               xyzh          = xyzh,                 &
               vxyzu         = vxyzu,                &
               polyk         = polyk,                &
               particle_mass = massoftype(igas),     &
               ! star_mass     = 1.0,                &
               disc_mass     = mdisc,                &
               inclination   = theta,                &
               ! bh_spin       = spin,               &
               ! alpha         = alpha,              &
               prefix        = fileprefix)

#ifdef GR
 a     = spin
 polyk = vxyzu(4,1)
#else
 blackhole_spin = spin
 polyk          = (5./3. -1.)*7.2e-6
#endif

 npartoftype(1) = npart

 return
end subroutine setpart


!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'
 call write_inopt(mhole ,'mhole' ,'mass of black hole (solar mass)'        , iunit)
 call write_inopt(mdisc ,'mdisc' ,'mass of disc       (solar mass)'        , iunit)
 call write_inopt(r_in  ,'r_in'  ,'inner edge of disc (GM/c^2, code units)', iunit)
 call write_inopt(r_out ,'r_out' ,'outer edge of disc (GM/c^2, code units)', iunit)
 call write_inopt(spin  ,'spin'  ,'spin parameter of black hole |a|<1'     , iunit)
 call write_inopt(honr  ,'honr'  ,'scale height H/R for disc'              , iunit)
 call write_inopt(theta ,'theta' ,'inclination of disc (degrees)'          , iunit)
 call write_inopt(accrad,'accrad','accretion radius   (GM/c^2, code units)', iunit)
 call write_inopt(np    ,'np'    ,'number of particles in disc'            , iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mhole ,'mhole' ,db,min=0.,errcount=nerr)
 call read_inopt(mdisc ,'mdisc' ,db,min=0.,errcount=nerr)
 call read_inopt(r_in  ,'r_in'  ,db,min=0.,errcount=nerr)
 call read_inopt(r_out ,'r_out' ,db,min=0.,errcount=nerr)
 call read_inopt(spin  ,'spin'  ,db,min=0.,errcount=nerr)
 call read_inopt(honr  ,'honr'  ,db,min=0.,errcount=nerr)
 call read_inopt(theta ,'theta' ,db,min=0.,errcount=nerr)
 call read_inopt(accrad,'accrad',db,min=0.,errcount=nerr)
 call read_inopt(np    ,'np   ' ,db,min=0 ,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
