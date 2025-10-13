!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup a polytrope in flat space, Minkowski metric (special relativity)
!
! :References: Liptai & Price (2019)
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - nr : *resolution (number of radial particles)*
!
! :Dependencies: checksetup, cons2prim, deriv, eos, infile_utils, io,
!   kernel, memory, metric_tools, part, physcon, prompting, rho_profile,
!   setup_params, spherical, timestep, units
!
 implicit none

 public :: setpart

 private

 integer :: nr = 25 ! Default number of particles in star radius

contains

!----------------------------------------------------------------
!+
!  Setup for SR polytrope
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,         only:igas,set_particle_type,rhoh,maxp
 use spherical,    only:set_sphere
 use units,        only:set_units,umass,udist
 use physcon,      only:solarm,solarr
 use io,           only:master,fatal
 use timestep,     only:tmax,dtmax
 use eos,          only:ieos
 use rho_profile,  only:rho_polytrope
 use prompting,    only:prompt
 use setup_params, only:npart_total
 use infile_utils, only:get_options,infile_exists
 use kernel,       only:hfact_default
 use metric_tools, only:init_metric
 use cons2prim,    only:prim2consall
 use deriv,        only:get_derivs_global
 use checksetup,   only:check_setup
 use memory,       only:allocate_memory
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer, parameter :: ntab=5000
 integer :: i,npts,ierr,nerror,nwarn
 real    :: psep
 real    :: rtab(ntab),rhotab(ntab)
 real    :: densi,mstar,rstar

 ! general parameters
 time  = 0.
 polyk = 1.e-10
 gamma = 5./3.
 ieos  = 2

 ! set tmax and dtmax if no infile found, otherwise we use whatever values it had
 if (.not.infile_exists(fileprefix)) then
    tmax  = 20000.
    dtmax = 100.
 endif

 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.
 hfact          = hfact_default
 ! set units
 call set_units(mass=1.e6*solarm,c=1.,G=1.)
 mstar = 1.*solarm/umass
 rstar = 1.*solarr/udist

 !
 !--Read runtime parameters from setup file
 !
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' SR polytrope'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

!-- resolution
 psep  = rstar/nr

!-- polytrope
 npart_total = 0
 call rho_polytrope(gamma,polyk,mstar,rtab,rhotab,npts,set_polyk=.true.,Rstar=rstar)
 call set_sphere('cubic',id,master,0.,rstar,psep,hfact,npart,xyzh,nptot=npart_total,&
                  xyz_origin=(/0.,0.,0./),rhotab=rhotab(1:npts),rtab=rtab(1:npts))

!-- mass and number of gas particles
 npartoftype(igas) = npart
 massoftype(igas)  = mstar/npart

 ! actually compute density so we can correctly set the entropy
 call check_setup(nerror,nwarn)
 call allocate_memory(int(maxp,kind=8)) ! allocate memory for tree
 call get_derivs_global()

 !-- set thermal energy from density
 do i=1,npart
    call set_particle_type(i,igas)
    densi        = rhoh(xyzh(4,i),massoftype(igas))
    vxyzu(4,i)   = polyk*densi**(gamma-1.) / (gamma-1.)
 enddo

 if (id==master) print "(/,a,i10,/)",' Number of particles setup = ',npart
 if (id==master) print*,' polyk = ',polyk
 if (id==master) print*,' mstar = ',mstar
 if (npart == 0) call fatal('setup','no particles setup')

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for SR polytrope setup'
 call write_inopt(nr,'nr','resolution (number of radial particles)',iunit)
 close(iunit)

end subroutine write_setupfile

!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
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
 call read_inopt(nr,'nr',db,min=0,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Interactive setup
!+
!-----------------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt

 call prompt('Resolution -- number of radial particles',nr,2)

end subroutine setup_interactive

end module setup
