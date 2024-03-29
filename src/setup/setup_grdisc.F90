!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! this module does an accretion disc setup in general relativity
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - accrad  : *accretion radius   (GM/c^2, code units)*
!   - alpha   : *artificial viscosity*
!   - gamma   : *adiabatic gamma*
!   - honr    : *scale height H/R of disc (at r_ref)*
!   - ismooth : *smooth inner edge of the disc (logical)*
!   - mdisc   : *mass of disc       (solar mass)*
!   - mhole   : *mass of black hole (solar mass)*
!   - np      : *number of particles in disc*
!   - p_index : *power law index of surface density profile*
!   - q_index : *power law index of sound speed profile*
!   - r_in    : *inner edge of disc (GM/c^2, code units)*
!   - r_out   : *outer edge of disc (GM/c^2, code units)*
!   - r_ref   : *reference radius   (GM/c^2, code units)*
!   - spin    : *spin parameter of black hole |a|<1*
!   - theta   : *inclination of disc (degrees)*
!
! :Dependencies: eos, extern_lensethirring, externalforces, infile_utils,
!   io, kernel, metric, options, part, physcon, prompting, setdisc,
!   timestep, units
!
 use options, only:alpha
 implicit none
 public :: setpart

 real,    private :: mhole,mdisc,r_in,r_out,r_ref,spin,honr,theta,p_index,q_index,accrad,gamma_ad
 integer, private :: np
 logical, private :: ismooth

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
 use options,        only:iexternalforce,alphau,iexternalforce
 use units,          only:set_units,umass
 use physcon,        only:solarm,pi
#ifdef GR
 use metric,         only:a
#else
 use externalforces,       only:iext_einsteinprec
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
 real    :: cs2

 time            = 0.
 alphau          = 0.0
 npartoftype(:)  = 0
 iexternalforce  = 1
 hfact           = hfact_default

#ifndef GR
 iexternalforce = iext_einsteinprec
#endif

 tmax  = 2.e4
 dtmax = 100.

 ieos  = 2
!
! Set default problem parameters
!

 mhole  = 1.e6    ! (solarm)
 mdisc  = 10.     ! (solarm)
 r_in   = 4.      ! (GM/c^2)
 r_out  = 160.    ! (GM/c^2)
 r_ref  = r_in
 ismooth= .true.
 spin   = 0.9
 honr   = 0.05
 alpha  = 0.368205218
 theta  = 3.      ! inclination angle (degrees)
 p_index= 1.5
 q_index= 0.75
 gamma_ad= 1.001
 np     = 1e6
 accrad = 4.      ! (GM/c^2)

!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",'Disc setup'
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

 !-- Set gamma from the option read from .setup file
 gamma = gamma_ad

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
               rref          = r_ref,                &
               p_index       = p_index,              &
               q_index       = q_index,              &
               HoverR        = honr,                 &
               ismooth       = ismooth,              &
               gamma         = gamma,                &
               hfact         = hfact,                &
               xyzh          = xyzh,                 &
               vxyzu         = vxyzu,                &
               polyk         = cs2,                  &
               particle_mass = massoftype(igas),     &
 ! star_mass     = 1.0,                &
               disc_mass     = mdisc,                &
               inclination   = theta,                &
               bh_spin       = spin,                 &
               prefix        = fileprefix)

#ifdef GR
 a     = spin
 ! Overwrite thermal energies to be correct for GR
 ! And use polyk to store the constant thermal energy
 ! polyk = cs2/(gamma-1.-cs2)
 ! vxyzu(4,:) = polyk

 polyk = 0. !?

#else
 blackhole_spin = spin
 polyk = cs2
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
 write(iunit,"(a)") '# input file for grdisc setup'
 call write_inopt(mhole  ,'mhole'  ,'mass of black hole (solar mass)'           , iunit)
 call write_inopt(mdisc  ,'mdisc'  ,'mass of disc       (solar mass)'           , iunit)
 call write_inopt(r_in   ,'r_in'   ,'inner edge of disc (GM/c^2, code units)'   , iunit)
 call write_inopt(r_out  ,'r_out'  ,'outer edge of disc (GM/c^2, code units)'   , iunit)
 call write_inopt(r_ref  ,'r_ref'  ,'reference radius   (GM/c^2, code units)'   , iunit)
 call write_inopt(ismooth,'ismooth','smooth inner edge of the disc (logical)'   , iunit)
 call write_inopt(spin   ,'spin'   ,'spin parameter of black hole |a|<1'        , iunit)
 call write_inopt(honr   ,'honr'   ,'scale height H/R of disc (at r_ref)'       , iunit)
 call write_inopt(alpha  ,'alpha'  ,'artificial viscosity'                      , iunit)
 call write_inopt(theta  ,'theta'  ,'inclination of disc (degrees)'             , iunit)
 call write_inopt(p_index,'p_index','power law index of surface density profile', iunit)
 call write_inopt(q_index,'q_index','power law index of sound speed profile'    , iunit)
 call write_inopt(gamma_ad,'gamma' ,'adiabatic gamma'                           , iunit)
 call write_inopt(accrad ,'accrad' ,'accretion radius   (GM/c^2, code units)'   , iunit)
 call write_inopt(np     ,'np'     ,'number of particles in disc'               , iunit)
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
 call read_inopt(mhole  ,'mhole'  ,db,min=0.,errcount=nerr)
 call read_inopt(mdisc  ,'mdisc'  ,db,min=0.,errcount=nerr)
 call read_inopt(r_in   ,'r_in'   ,db,min=0.,errcount=nerr)
 call read_inopt(r_out  ,'r_out'  ,db,min=0.,errcount=nerr)
 call read_inopt(r_ref  ,'r_ref'  ,db,min=0.,errcount=nerr)
 call read_inopt(ismooth,'ismooth',db,errcount=nerr)
 call read_inopt(spin   ,'spin'   ,db,min=-1.,max=1.,errcount=nerr)
 call read_inopt(honr   ,'honr'   ,db,min=0.,errcount=nerr)
 call read_inopt(alpha  ,'alpha'  ,db,min=0.,errcount=nerr)
 call read_inopt(theta  ,'theta'  ,db,min=0.,max=90.,errcount=nerr)
 call read_inopt(p_index,'p_index',db,errcount=nerr)
 call read_inopt(q_index,'q_index',db,errcount=nerr)
 call read_inopt(gamma_ad,'gamma' ,db,min=1.,errcount=nerr)
 call read_inopt(accrad ,'accrad' ,db,min=0.,errcount=nerr)
 call read_inopt(np     ,'np   '  ,db,min=0 ,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
