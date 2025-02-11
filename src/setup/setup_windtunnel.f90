!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! this module does setup
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters:
!   - Mstar            : *sphere mass in code units*
!   - Rstar            : *sphere radius in code units*
!   - gamma            : *adiabatic index*
!   - handled_layers   : *number of handled layers*
!   - lattice_type     : *0: cubic, 1: close-packed cubic*
!   - pres_inf         : *wind pressure / dyn cm^2*
!   - rho_inf          : *wind density / g cm^-3*
!   - v_inf            : *wind speed / km s^-1*
!   - wind_injection_x : *injection x in units of Rstar*
!   - wind_length      : *wind length in units of Rstar*
!   - wind_radius      : *injection radius in units of Rstar*
!
! :Dependencies: dim, eos, extern_densprofile, infile_utils, inject, io,
!   kernel, mpidomain, part, physcon, relaxstar, rho_profile,
!   setstar_utils, setunits, setup_params, table_utils, timestep, unifdis,
!   units
!
 use io,     only:master,fatal
 use inject, only:init_inject,nstar,Rstar,lattice_type,handled_layers,&
                  wind_radius,wind_injection_x,wind_length,&
                  rho_inf,v_inf,mach

 implicit none
 public :: setpart

 real    :: Mstar,pmass
 integer :: nstar_in  ! guess for no. of star prticles
 logical :: add_star = .false.

 private

contains

!----------------------------------------------------------------
!+
!  setup for polytropic gas sphere inside wind tunnel
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,        only:ihsoft,igas
 use eos,         only:ieos,gmw
 use setstar_utils,only:set_star_density
 use rho_profile, only:rho_polytrope
 use relaxstar,   only:relax_star
 use extern_densprofile, only:nrhotab
 use physcon,     only:solarm,solarr,pi
 use units,       only:udist,umass,utime,set_units,unit_velocity,unit_density
 use setunits,    only:mass_unit,dist_unit
 use mpidomain,   only:i_belong
 use timestep,    only:dtmax,tmax
 use unifdis,     only:mask_prototype
 use kernel,      only:hfact_default
 use setup_params,only:rhozero,npart_total
 use table_utils, only:yinterp
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:),vxyzu(:,:),massoftype(:),polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real               :: rhocentre,rmin,densi,presi,ri
 real, allocatable  :: r(:),den(:),pres(:),Xfrac(:),Yfrac(:),mu(:)
 integer            :: ierr,ierr_relax,npts,np,i
 logical            :: use_exactN,setexists,use_var_comp
 character(len=30)  :: lattice
 character(len=120) :: setupfile

 call set_units(mass=solarm,dist=solarr,G=1.)
 !
 ! Initialise parameters, including those that will not be included in *.setup
 !
 ! units
 mass_unit = 'solarm'
 dist_unit = 'solarr'

 time  = 0.
 polyk = 1. ! not used but needs to be initialised to non-zero value
 gamma = 5./3.
 ieos  = 2
 gmw   = 0.6
 hfact = hfact_default

 ! Wind parameters (see inject_windtunnel module)
 v_inf    = 230e5 / unit_velocity
 rho_inf  = 4.e-2 / unit_density
 mach     = 13.86

 ! Star parameters
 add_star = .false.
 Rstar = 0.1
 Mstar = 1.e-3
 nstar_in = 1000
 lattice = 'closepacked'
 use_exactN = .true.

 ! Wind injection settings
 lattice_type = 1
 handled_layers = 4
 wind_radius = 1.         ! in code units
 wind_injection_x = -0.2  ! in code units
 wind_length = 1.5        ! in code units

 ! default setting for particle mass based on
 ! filling the injection box with 10*nstar_in particles
 pmass = rho_inf*wind_length*pi*wind_radius**2/(10.*nstar_in)

 ! Set default tmax and dtmax
 dtmax = 0.1
 tmax  = 6.8

 ! determine if the .setup file exists
 setupfile = trim(fileprefix)//'.setup'
 inquire(file=setupfile,exist=setexists)
 if (setexists) then
    call read_setupfile(setupfile,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(setupfile)
       stop 'please rerun phantomsetup with revised .setup file'
    endif
    !--Prompt to get inputs and write to file
 elseif (id==master) then
    print "(a,/)",trim(setupfile)//' not found: using default parameters'
    call write_setupfile(setupfile)
    stop 'please check and edit .setup file and rerun phantomsetup'
 endif
 !
 ! override the particle mass if we are setting up a star based
 ! on the number of particles desired inside the star
 !
 if (add_star) then
    pmass = Mstar / real(nstar_in)
    call check_setup(pmass,ierr)
    if (ierr /= 0) call fatal('windtunnel','errors in setup parameters')
 endif

 massoftype(igas) = pmass

 ! Initialise particle injection
 call init_inject(ierr)
 npart = 0
 np = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.

 ! Set polytropic star
 if (add_star) then
    allocate(r(nrhotab),den(nrhotab),pres(nrhotab))
    call rho_polytrope(gamma,polyk,Mstar,r,den,npts,rhocentre,set_polyk=.true.,Rstar=Rstar)
    pres = polyk*den**gamma
    rmin = r(1)
    call set_star_density(lattice,id,master,rmin,Rstar,Mstar,hfact,&
                       npts,den,r,npart,npartoftype,massoftype,xyzh,&
                       use_exactN,np,rhozero,npart_total,i_belong) ! Note: mass_is_set = .true., so np is not used
    nstar = npart
    use_var_comp = .false.
    call relax_star(npts,den,pres,r,npart,xyzh,use_var_comp,Xfrac,Yfrac,mu,ierr_relax)

    ! Set thermal energy
    do i = 1,npart
       ri = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
       densi = yinterp(den(1:npts),r(1:npts),ri)
       presi = yinterp(pres(1:npts),r(1:npts),ri)
       vxyzu(4,i) =  presi / ( (gamma-1.) * densi)
    enddo

    deallocate(r,den,pres)
 endif

 print*, "udist = ", udist, "; umass = ", umass, "; utime = ", utime

end subroutine setpart



!-----------------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils,  only:write_inopt
 use dim,           only:tagline
 use eos,           only:gamma
 use setunits,      only:write_options_units
 use units,         only:unit_density,unit_velocity
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20

 write(*,"(a)") ' Writing '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# '//trim(tagline)
 write(iunit,"(a)") '# input file for Phantom wind tunnel setup'

 call write_options_units(iunit)

 write(iunit,"(/,a)") '# sphere settings'
 call write_inopt(add_star,'add_star','add polytropic star to the wind tunnel',iunit)
 if (add_star) then
    call write_inopt(nstar_in,'nstar','number of particles resolving gas sphere',iunit)  ! note: this is an estimate, actual no. of particles is npart outputted from set_sphere
    call write_inopt(Mstar,'Mstar','sphere mass in code units',iunit)
    call write_inopt(Rstar,'Rstar','sphere radius in code units',iunit)
 endif

 write(iunit,"(/,a)") '# wind settings'
 call write_inopt(v_inf*unit_velocity/1.e5,'v_inf','wind speed / km s^-1',iunit)
 call write_inopt(rho_inf*unit_density,'rho_inf','wind density / g cm^-3',iunit)
 call write_inopt(mach,'mach','wind mach number',iunit)
 call write_inopt(gamma,'gamma','adiabatic index',iunit)
 if (.not.add_star) then
    call write_inopt(pmass,'pmass','particle mass in code units (i.e. the mass resolution)',iunit)
 endif

 write(iunit,"(/,a)") '# wind injection settings'
 call write_inopt(lattice_type,'lattice_type','0: cubic, 1: close-packed cubic',iunit)
 call write_inopt(handled_layers,'handled_layers','number of handled layers',iunit)
 call write_inopt(wind_radius,'wind_radius','injection radius in code units',iunit)
 call write_inopt(wind_injection_x,'wind_injection_x','injection x in code units',iunit)
 call write_inopt(wind_length,'wind_length','wind length in code units',iunit)

 close(iunit)

end subroutine write_setupfile


!-----------------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils,  only:open_db_from_file,inopts,close_db,read_inopt
 use io,            only:error
 use units,         only:select_unit,unit_density,unit_velocity
 use setunits,      only:read_options_and_set_units
 use eos,           only:gamma
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer,          parameter   :: lu = 21
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return

 nerr = 0

 call read_options_and_set_units(db,nerr)

 call read_inopt(add_star,'add_star',db,errcount=nerr)
 if (add_star) then
    call read_inopt(nstar_in,'nstar',db,errcount=nerr)
    call read_inopt(Mstar,'Mstar',db,errcount=nerr)
    call read_inopt(Rstar,'Rstar',db,errcount=nerr)
 else
    call read_inopt(pmass,'pmass',db,errcount=nerr)
 endif

 call read_inopt(v_inf,'v_inf',db,errcount=nerr)
 call read_inopt(rho_inf,'rho_inf',db,errcount=nerr)
 call read_inopt(mach,'mach',db,errcount=nerr)
 call read_inopt(gamma,'gamma',db,errcount=nerr)

 ! Convert wind quantities to code units
 v_inf = v_inf / unit_velocity * 1.e5
 rho_inf = rho_inf / unit_density

 call read_inopt(lattice_type,'lattice_type',db,errcount=nerr)
 call read_inopt(handled_layers,'handled_layers',db,errcount=nerr)
 call read_inopt(wind_radius,'wind_radius',db,errcount=nerr)
 call read_inopt(wind_injection_x,'wind_injection_x',db,errcount=nerr)
 call read_inopt(wind_length,'wind_length',db,errcount=nerr)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_windtunnel: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile


!-----------------------------------------------------------------------
!+
!  Check that setup is sensible
!+
!-----------------------------------------------------------------------
subroutine check_setup(pmass,ierr)
 real, intent(in)     :: pmass
 integer, intent(out) :: ierr
 real                 :: min_layer_sep

 ierr = 0

 min_layer_sep = (pmass / rho_inf)**(1./3.)

 if ( abs(wind_injection_x - Rstar) < real(handled_layers)*min_layer_sep ) then
    print*,'error: Handled layers overlap with sphere. Try decreasing wind_injection_x or handled_layers'
    ierr = 1
 endif
 if (wind_radius < 1.) then
    print*,'error: Wind cross-section should not be smaller than the sphere'
    ierr = 1
 endif
 if ( wind_injection_x + wind_length < 1. ) then
    print*,'error: Wind not long enough to cover initial sphere position. Try increasing wind_injection_x or wind_length'
    ierr = 1
 endif

end subroutine check_setup

end module setup

