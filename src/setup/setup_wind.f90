!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    central_star_mass   -- mass of the central star (Msun)
!    central_star_radius -- radius of the central star (au)
!    companion_star_mass -- mass of the companion star (Msun)
!    companion_star_r    -- radius of the companion star (au)
!    eccentricity        -- eccentricity of the binary system
!    icompanion_star     -- set to 1 for a binary system
!    mass_of_particles   -- mass resolution (Msun)
!    semi_major_axis     -- semi-major axis of the binary system (au)
!    wind_gamma          -- polytropic index
!
!  DEPENDENCIES: infile_utils, inject, io, part, physcon, prompting,
!    setbinary, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private
 real, public :: wind_gamma = 5./3.
 real, public :: central_star_mass = 1.
 real, public :: central_star_radius = 1.
 integer, public :: icompanion_star = 0
 real, public :: companion_star_mass
 real, public :: companion_star_r
 real, public :: semi_major_axis
 real, public :: eccentricity
 real, public :: mass_of_particles = 1.e-6

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only: xyzmh_ptmass, vxyz_ptmass, nptmass
 use physcon,   only: au, solarm
 use units,     only: umass, set_units
 use inject,    only: wind_init, mass_of_particles
 use setbinary, only: set_binary
 use io,        only: master
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=*), intent(in)    :: fileprefix
 character(len=len(fileprefix)+6) :: filename
 integer :: ierr
 logical :: iexist

 call set_units(dist=au,mass=solarm,G=1.)
!
!--general parameters
!
 time = 0.
 polyk = 0.
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call setup_interactive()
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
 endif

 gamma = wind_gamma
 call wind_init(.true.)
!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0

 massoftype = mass_of_particles / umass

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.

 if (icompanion_star > 0) then
    call set_binary(central_star_mass, &
                    companion_star_mass/central_star_mass, &
                    semi_major_axis, &
                    eccentricity, &
                    central_star_radius, &
                    companion_star_r, &
                    xyzmh_ptmass, vxyz_ptmass, nptmass)
 else
    nptmass = 1
    xyzmh_ptmass(4,1) = central_star_mass
    xyzmh_ptmass(5,1) = central_star_radius
 endif

end subroutine setpart

!----------------------------------------------------------------
!+
!  determine which problem to set up interactively
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 use physcon,   only:au,solarm
 use units,     only:umass,udist
 integer :: iproblem

 iproblem = 1
 call prompt('Which defaults to use? (1=adiabatic wind 2=Bowen)',iproblem,0,1)
 call prompt('Add binary?',icompanion_star,0,1)
 select case(iproblem)
 case(2)
    central_star_mass = 1.2 * (solarm / umass)
    central_star_radius = 1.2568 * (au / udist)
 case default
    central_star_mass = 1. * (solarm / umass)
    central_star_radius = 1. * (au / udist)
 end select

end subroutine setup_interactive

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for wind setup routine'
 call write_inopt(central_star_mass,'central_star_mass','mass of the central star (Msun)',iunit)
 call write_inopt(central_star_radius,'central_star_radius','radius of the central star (au)',iunit)
 call write_inopt(icompanion_star,'icompanion_star','set to 1 for a binary system',iunit)
 if (icompanion_star > 0) then
    call write_inopt(companion_star_mass,'companion_star_mass','mass of the companion star (Msun)',iunit)
    call write_inopt(companion_star_r,'companion_star_r','radius of the companion star (au)',iunit)
    call write_inopt(semi_major_axis,'semi_major_axis','semi-major axis of the binary system (au)',iunit)
    call write_inopt(eccentricity,'eccentricity','eccentricity of the binary system',iunit)
 endif
 call write_inopt(mass_of_particles,'mass_of_particles','mass resolution (Msun)',iunit)
 call write_inopt(wind_gamma,'wind_gamma','polytropic index',iunit)
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)
 integer :: nerr

 nerr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(central_star_mass,'central_star_mass',db,min=0.,max=1000.,errcount=nerr)
 call read_inopt(central_star_radius,'central_star_radius',db,min=0.,errcount=nerr)
 call read_inopt(icompanion_star,'icompanion_star',db,min=0,errcount=nerr)
 if (icompanion_star > 0) then
    call read_inopt(companion_star_mass,'companion_star_mass',db,min=0.,max=1000.,errcount=nerr)
    call read_inopt(companion_star_r,'companion_star_r',db,min=0.,max=1000.,errcount=nerr)
    call read_inopt(semi_major_axis,'semi_major_axis',db,min=0.,errcount=nerr)
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
    call read_inopt(wind_gamma,'wind_gamma',db,min=1.,max=4.,errcount=nerr)
 endif
 call read_inopt(mass_of_particles,'mass_of_particles',db,min=0.,errcount=nerr)
 call close_db(db)
 ierr = nerr

end subroutine read_setupfile
end module setup
