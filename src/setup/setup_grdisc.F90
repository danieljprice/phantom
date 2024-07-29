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
! :Dependencies: dim, eos, extern_lensethirring, externalforces,
!   infile_utils, io, kernel, metric, mpidomain, options, part, physcon,
!   prompting, setdisc, setorbit, setstar, setunits, setup_params,
!   timestep, units
!
 use options,  only:alpha
 use setstar,  only:star_t
 use setorbit, only:orbit_t
 implicit none
 public :: setpart

 real,    private :: mhole,mdisc,r_in,r_out,r_ref,spin,honr,theta,p_index,q_index,accrad,gamma_ad
 integer, private :: np,nstars
 logical, private :: ismooth,relax
 integer, parameter :: max_stars = 10
 type(star_t), private :: star(max_stars)
 type(orbit_t),private :: orbit(max_stars)

 private

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,        only:set_disc
 use part,           only:igas,nsinkproperties,eos_vars,rad,xyzmh_ptmass,vxyz_ptmass,nptmass
 use io,             only:master
 use externalforces, only:accradius1,accradius1_hard
 use options,        only:iexternalforce,alphau,iexternalforce,ipdv_heating,ishock_heating
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
 use eos,            only:ieos,use_var_comp,X_in,Z_in
 use kernel,         only:hfact_default
 use setstar,        only:shift_star,set_stars
 use setorbit,       only:set_defaults_orbit,set_orbit
 use setunits,       only:mass_unit
 use mpidomain,      only:i_belong
 use setup_params,   only:rhozero
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
 integer :: ierr,nptmass_in,i
 integer(kind=8) :: npart_total
 logical :: iexist,write_profile
 real    :: cs2
 real :: xyzmh_ptmass_in(nsinkproperties,2),vxyz_ptmass_in(3,2)

 time            = 0.
 alphau          = 0.0
 npartoftype(:)  = 0
 nptmass         = 0
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
 mass_unit = '1e6*solarm'
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
 gamma_ad= 5./3.
 np     = 1e6
 accrad = 4.      ! (GM/c^2)
 accradius1 = accrad
 gamma = gamma_ad
 ! default units
 call set_units(G=1.,c=1.,mass=mhole*solarm) ! Set central mass to M=1 in code units

 ! stars
 nstars = 0
 do i=1,size(orbit)
    call set_defaults_orbit(orbit(i))
 enddo
 relax = .true.

!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),(/,1x,a),/,65('-'),/)",'General relativistic disc setup'
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

 !
 ! add stars on desired orbits around the black hole, these could be
 ! either sink particles or balls of gas
 !
 if (nstars > 0) then
    write_profile = .false.
    call set_stars(id,master,nstars,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                   massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
                   X_in,Z_in,relax,use_var_comp,write_profile,&
                   rhozero,npart_total,i_belong,ierr)
    do i=1,nstars
       nptmass_in = 0
       call set_orbit(orbit(i),mhole/umass,star(i)%mstar,r_in,star(i)%rstar, &
                      xyzmh_ptmass_in,vxyz_ptmass_in,nptmass_in,(id==master),ierr)

       ! shift the star to the position of the second body
       if (star(i)%iprofile > 0) then
          call shift_star(npart,xyzh,vxyzu,x0=xyzmh_ptmass_in(:,2),v0=vxyz_ptmass_in(:,2),itype=i)
       else
          nptmass = nptmass + 1
          xyzmh_ptmass(:,nptmass) = xyzmh_ptmass_in(:,2)
          vxyz_ptmass(:,nptmass) = vxyz_ptmass_in(:,2)
       endif
    enddo
 endif

 ipdv_heating = 0
 ishock_heating = 0
 if (id==master) print "(/,a,/)",' ** SETTING ipdv_heating=0 and ishock_heating=0 for grdisc setup **'

end subroutine setpart


!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setstar,      only:write_options_stars
 use setorbit,     only:write_options_orbit
 use setunits,     only:write_options_units
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_options_units(iunit,gr=.true.)

 write(iunit,"(/,a)") '# disc parameters'
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

 write(iunit,"(/,a)") '# stars'
 call write_options_stars(star,relax,iunit,nstar=nstars)
 do i=1,nstars
    call write_options_orbit(orbit(i),iunit,label=achar(i+48))
 enddo
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 use setstar,      only:read_options_stars
 use setorbit,     only:read_options_orbit
 use eos,          only:ieos,polyk
 use setunits,     only:read_options_and_set_units
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr,need_iso,i
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr,gr=.true.)
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
 call read_options_stars(star,need_iso,ieos,polyk,relax,db,nerr,nstars)
 do i=1,nstars
    call read_options_orbit(orbit(i),db,nerr,label=achar(i+48))
 enddo
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
