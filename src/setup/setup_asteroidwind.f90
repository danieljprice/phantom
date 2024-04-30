!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Wind from an evaporating/disintegrating asteroid, as used
! in Trevascus et al. (2021)
!
! :References:
!   Trevascus et al. (2021), MNRAS 505, L21-L25
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - dumpsperorbit : *number of dumps per orbit*
!   - eccentricity  : *eccentricity*
!   - gastemp       : *gas temperature in K*
!   - hacc1         : *white dwarf (sink) accretion radius (solar radii)*
!   - ipot          : *wd modelled by 0=sink or 1=externalforce*
!   - m1            : *mass of white dwarf (solar mass)*
!   - m2            : *mass of asteroid (ceres mass)*
!   - mdot          : *mass injection rate (g/s)*
!   - norbits       : *number of orbits*
!   - npart_at_end  : *number of particles injected after norbits*
!   - rasteroid     : *radius of asteroid (km)*
!   - semia         : *semi-major axis (solar radii)*
!
! :Dependencies: eos, extern_lensethirring, externalforces, infile_utils,
!   inject, io, kernel, options, part, physcon, setbinary, spherical,
!   timestep, units
!
 use inject, only:mdot
 implicit none
 public :: setpart

 real :: m1,m2,ecc,semia,hacc1,rasteroid,norbits,gastemp
 integer :: npart_at_end,dumpsperorbit,ipot

 private

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,idust,set_particle_type,igas
 use setbinary, only:set_binary,get_a_from_period
 use spherical, only:set_sphere
 use units,     only:set_units,umass,udist,utime,unit_velocity
 use physcon,   only:solarm,au,pi,solarr,ceresm,km,kboltz,mass_proton_cgs
 use externalforces,   only:iext_binary, iext_einsteinprec, update_externalforce, &
                            mass1,accradius1
 use io,        only:master,fatal
 use timestep,  only:tmax,dtmax
 !use inject,    only:inject_particles
 use eos,       only:gmw
 use options,   only:iexternalforce
 use extern_lensethirring, only:blackhole_spin
 use kernel,    only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(inout)   :: massoftype(:)
 real,              intent(inout)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename
 integer :: ierr
 logical :: iexist
 real    :: period,hacc2,temperature_coef
 real    :: rp
!
!--Default runtime parameters (values for SDSS J1228+1040)
!
 ipot          = 1         ! (0=sink or 1=externalforce)
 m1            = 0.705     ! (solar masses)
 m2            = 0.1       ! (ceres masses)
 ecc           = 0.54      ! (eccentricity)
 semia         = 0.73      ! (solar radii)
 hacc1         = 0.1679    ! (solar radii)
 rasteroid     = 2338.3      ! (km)
 gastemp       = 5000.     ! (K)
 norbits       = 1000.
 mdot          = 5.e8      ! Mass injection rate (g/s)
 npart_at_end  = 1.0e6       ! Number of particles after norbits
 dumpsperorbit = 1

!
!--Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Asteroid wind'
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

!
!-- Set units
!
 m1 = m1*solarm
 if (ipot == 0) then
    call set_units(mass=solarm,dist=solarr,G=1.d0)
 else
    call set_units(c=1.0,G=1.0,mass=m1)
    print*,'Code units changed to c=G=m1=1.0 for LT'
 endif

!
!--Convert to code units
!
 m1    = m1/umass
 m2    = m2*ceresm/umass

 semia = semia*solarr/udist
 hacc1 = hacc1*solarr/udist
 rasteroid = rasteroid*km/udist

!
!--general parameters
!
 time             = 0.

 temperature_coef = mass_proton_cgs/kboltz * unit_velocity**2
 polyk            = gastemp/(temperature_coef*gmw)
 gamma            = 1.

!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 period = sqrt(4.*pi**2*semia**3/(m1+m2))
 hacc2  = 0.                                 ! asteroid should not accrete
 tmax   = norbits*period
 dtmax  = period/dumpsperorbit

! If using sink particle for central mass, use binary to setup orbits
! and have two point masses. Otherwise use iexternalforce and only one
! point mass for the asteroid

 if (ipot == 0) then
    !
    !--Set a binary orbit given the desired orbital parameters
    !
    call set_binary(m1,m2,semia,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr)
    xyzmh_ptmass(ihsoft,2) = rasteroid ! Asteroid radius softening

 else

    !
    !--Set the asteroid on orbit around the fixed potential
    !
    mass1                = m1
    accradius1           = hacc1
    iexternalforce       = 11
    blackhole_spin       = 0.
    call update_externalforce(iexternalforce,time,0.)

    ! Orbit and position
    nptmass = 1
    xyzmh_ptmass(1:3,1) = (/semia*(1. + ecc),0.,0./)
    vxyz_ptmass(1:3,1)  = (/0.,sqrt(semia*(1.-ecc**2)*(m1+m2))/xyzmh_ptmass(1,1),0./)

    xyzmh_ptmass(4,1)      = m2
    xyzmh_ptmass(ihacc,1)  = hacc2        ! asteroid should not accrete
    xyzmh_ptmass(ihsoft,1) = rasteroid    ! asteroid radius softening
 endif

 ! we use the estimated injection rate and the final time to set the particle mass
 massoftype(igas) = tmax*mdot/(umass/utime)/npart_at_end
 hfact = hfact_default
 !npart_old = npart
 !call inject_particles(time,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinj)

!
!-- check for silly parameter choices
!
 rp = semia*(1. - ecc)
 if (rp < hacc1)   call fatal('setup','periapsis is within racc of central sink')
 if (ipot > 1)     call fatal('setup','choice of potential not recognised, try 1')
 if (nptmass == 0) call fatal('setup','no sink particles setup')
 !if (npart == 0)   call fatal('setup','no hydro particles setup')
 if (ierr /= 0)    call fatal('setup','ERROR during setup')

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
 call write_inopt(ipot,         'ipot',         'wd modelled by 0=sink or 1=externalforce',         iunit)
 call write_inopt(m1,           'm1',           'mass of white dwarf (solar mass)',                 iunit)
 call write_inopt(m2,           'm2',           'mass of asteroid (ceres mass)',                    iunit)
 call write_inopt(ecc,          'ecc',          'eccentricity',                                     iunit)
 call write_inopt(semia,        'semia',        'semi-major axis (solar radii)',                    iunit)
 call write_inopt(hacc1,        'hacc1',        'white dwarf (sink) accretion radius (solar radii)',iunit)
 call write_inopt(rasteroid,    'rasteroid',    'radius of asteroid (km)',                          iunit)
 call write_inopt(gastemp,      'gastemp',      'gas temperature in K',                             iunit)
 call write_inopt(norbits,      'norbits',      'number of orbits',                                 iunit)
 call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',                        iunit)
 call write_inopt(npart_at_end,'npart_at_end','number of particles injected after norbits',iunit)
 call write_inopt(mdot,'mdot','mass injection rate (g/s)',iunit)
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
 call read_inopt(ipot,         'ipot',         db,min=0 ,errcount=nerr)
 call read_inopt(m1,           'm1',           db,min=0.,errcount=nerr)
 call read_inopt(m2,           'm2',           db,min=0.,errcount=nerr)
 call read_inopt(ecc,          'ecc',          db,min=0.,errcount=nerr)
 call read_inopt(semia,        'semia',        db,min=0.,errcount=nerr)
 call read_inopt(hacc1,        'hacc1',        db,min=0.,errcount=nerr)
 call read_inopt(rasteroid,    'rasteroid',    db,min=0.,errcount=nerr)
 call read_inopt(gastemp,      'gastemp',      db,min=0.,errcount=nerr)
 call read_inopt(norbits,      'norbits',      db,min=0.,errcount=nerr)
 call read_inopt(dumpsperorbit,'dumpsperorbit',db,min=0 ,errcount=nerr)
 call read_inopt(npart_at_end, 'npart_at_end', db,min=0 ,errcount=nerr)
 call read_inopt(mdot,         'mdot',         db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
