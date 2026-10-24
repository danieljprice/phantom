!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testsetstar
!
! Unit tests of set_star and relax_star modules
!
! :References: None
!
! :Owner: Francesco Radica
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, checksetup, damping, datafiles, deriv, dim,
!   energies, eos, gravwaveutils, infile_utils, io, kernel, mpidomain, options, part,
!   physcon, ptmass, radiation_utils, readwrite_mesa, relaxstar, setstar,
!   setstar_utils, sortutils, stellar_hcp, step_lf_global, table_utils,
!   testutils, timestep, units, vectorutils
!
 use testutils, only:checkval,update_test_scores
 implicit none
 public :: test_setstar

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of the set_star routine and associated utilities
!+
!-----------------------------------------------------------------------
subroutine test_setstar(ntests,npass,string)
 use checksetup, only:check_setup
 use io,         only:id,master
 use dim,        only:gravity,do_radiation,mpi,use_apr,ind_timesteps,mhd,use_dust,gr,maxvxyzu
 integer, intent(inout) :: ntests,npass
 character(len=*), intent(in), optional :: string

 if (.not.gravity) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING STAR SETUP TESTS'
    return
 else
    if (id==master) write(*,"(/,a,/)") '--> TESTING STAR SETUP'
 endif

 ! test the get_mass_coordinate routine in relax_star
 call test_get_mass_coord(ntests,npass)
 call test_spin_options(ntests,npass)
 call test_weighted_hcp(ntests,npass)
 call test_hcp_surface_buffer(ntests,npass)
 call test_hcp_options(ntests,npass)
 call test_composition_without_mu(ntests,npass)
 if (.not.(do_radiation.or.mpi.or.use_apr.or.ind_timesteps.or.mhd.or.use_dust.or.gr).and.maxvxyzu>=4) &
    call test_analytic_hcp(ntests,npass)
 if (present(string)) then
    if (string=='hcp'.or.string=='weighted_hcp') return
 endif

 ! test that a star moved onto a circular orbit is set corotating
 call test_shift_star_corotate(ntests,npass)

 ! test polytrope
 if (.not. do_radiation) call test_polytrope(ntests,npass)

 ! test red supergiant
 call test_redsupergiant(ntests,npass)

 ! test white dwarf
 call test_whitedwarf(ntests,npass)

end subroutine test_setstar

subroutine test_weighted_hcp(ntests,npass)
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 use stellar_hcp, only:generate_stellar_hcp,natural_spline,spline_value
 real, parameter :: tol = max(1.e-12,64.*epsilon(1.))
 integer, intent(inout) :: ntests,npass
 real, parameter :: pi = acos(-1.d0)
 real :: r(101),den(101),menc(101),rho2(101),redge,a1,ri,h_expected
 real, allocatable :: xyz(:,:),h(:),pmass(:),rhoi(:)
 integer :: i,ierr,nfailed(9),nfailed_origin(3)

 do i=1,size(r)
    r(i) = real(i)/size(r)
    den(i) = 1. - 0.5*r(i)**2
    menc(i) = 4.*pi*(r(i)**3/3. - r(i)**5/10.)
 enddo
 call generate_stellar_hcp(r,den,menc,1.,menc(101),10000,0.,23,1.3,0.,&
                           xyz,h,pmass,rhoi,redge,a1,ierr)
 nfailed = 0
 call checkval(ierr,0,0,nfailed(1),'weighted HCP construction')
 call update_test_scores(ntests,nfailed(1:1),npass)
 if (ierr /= 0) return
 call natural_spline(r,den,rho2)
 ri = sqrt(sum(xyz(:,1)**2))
 h_expected = (3./32./pi*1.9*23*pmass(1)/rhoi(1))**(1./3.)
 call checkval(mod(size(h),2),0,0,nfailed(2),'antipodal particle count')
 call checkval(maxval(abs(sum(xyz,dim=2))),0.,tol,nfailed(3),'symmetric HCP sphere')
 call checkval(sum(pmass),menc(101),tol,nfailed(4),'normalised stellar mass')
 call checkval(maxval(pmass/rhoi)-minval(pmass/rhoi),0.,tol,nfailed(5),'uniform HCP number density')
 call checkval(rhoi(1),spline_value(r,den,rho2,ri),tol,nfailed(6),'cubic density interpolation')
 call checkval(h(1),h_expected,tol,nfailed(7),'initial HCP smoothing length')
 call checkval(maxval(abs(xyz(:,1)+xyz(:,2))),0.,0.,nfailed(8),'antipodal pair coordinates')
 call checkval(maxval(sqrt(sum(xyz**2,dim=1))) < redge,.true.,nfailed(9),'support boundary')
 call update_test_scores(ntests,nfailed(2:),npass)

 ! Analytic profiles may include the origin, unlike the MESA reference grid.
 do i=1,size(r)
    r(i) = real(i-1)/(size(r)-1)
    den(i) = 1. - 0.5*r(i)**2
    menc(i) = 4.*pi*(r(i)**3/3. - r(i)**5/10.)
 enddo
 call generate_stellar_hcp(r,den,menc,1.,menc(101),10000,1.,23,1.3,0.,&
                           xyz,h,pmass,rhoi,redge,a1,ierr)
 nfailed_origin = 0
 call checkval(ierr,0,0,nfailed_origin(1),'mapped HCP profile including the origin')
 if (ierr == 0) then
    call checkval(all(ieee_is_finite(xyz)).and.all(ieee_is_finite(h)),.true.,&
                  nfailed_origin(2),'finite central HCP mapping')
    call checkval(maxval(pmass)-minval(pmass),0.,tol,nfailed_origin(3),'mapped HCP equal-mass limit')
 else
    nfailed_origin(2:) = 1
 endif
 call update_test_scores(ntests,nfailed_origin,npass)

end subroutine test_weighted_hcp

subroutine test_hcp_surface_buffer(ntests,npass)
 use, intrinsic :: ieee_arithmetic, only:ieee_value,ieee_quiet_nan
 use stellar_hcp, only:generate_stellar_hcp
 use testutils, only:checkval,update_test_scores
 real, parameter :: tol = max(1.e-12,64.*epsilon(1.))
 integer, intent(inout) :: ntests,npass
 real :: r(1001),den(1001),menc(1001),redge,a1,hc,automatic_edge
 real, allocatable :: xyz(:,:),h(:),pmass(:),rhoi(:)
 integer :: i,ierr,nfailed(7)

 do i=1,size(r)
    r(i) = real(i)/size(r)
    den(i) = (1.+(r(i)/0.2)**2)**(-2.5)
    menc(i) = r(i)**3/(r(i)**2+0.2**2)**1.5
 enddo
 hc = 0.5*(1.3*23/10000.)**(1./3.)
 nfailed = 0
 call generate_stellar_hcp(r,den,menc,1.,menc(size(r)),10000,0.,23,1.3,0.,&
                          xyz,h,pmass,rhoi,redge,a1,ierr)
 automatic_edge = redge
 call checkval(redge,1.-2.*hc,tol,nfailed(1),'default single-zone surface buffer')
 call generate_stellar_hcp(r,den,menc,1.,menc(size(r)),10000,0.,23,1.3,0.,&
                          xyz,h,pmass,rhoi,redge,a1,ierr,surface_buffer=3.)
 call checkval(ierr,0,0,nfailed(2),'explicit surface buffer construction')
 call checkval(redge,1.-3.*hc,tol,nfailed(3),'three h_c surface buffer')
 if (ierr == 0) then
    call checkval(sum(pmass),menc(size(r)),tol,nfailed(4),'buffered stellar mass normalisation')
 else
    nfailed(4) = 1
 endif
 call generate_stellar_hcp(r,den,menc,1.,menc(size(r)),10000,0.,23,1.3,0.,&
                          xyz,h,pmass,rhoi,redge,a1,ierr,surface_buffer=-1.)
 call checkval(redge,automatic_edge,0.,nfailed(5),'automatic buffer preserves default')
 call generate_stellar_hcp(r,den,menc,1.,menc(size(r)),10000,0.,23,1.3,0.,&
                          xyz,h,pmass,rhoi,redge,a1,ierr,surface_buffer=0.)
 call checkval(ierr,1,0,nfailed(6),'zero surface buffer rejected')
 call generate_stellar_hcp(r,den,menc,1.,menc(size(r)),10000,0.,23,1.3,0.,&
                          xyz,h,pmass,rhoi,redge,a1,ierr,surface_buffer=ieee_value(0.,ieee_quiet_nan))
 call checkval(ierr,1,0,nfailed(7),'nonfinite surface buffer rejected')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_hcp_surface_buffer

subroutine test_hcp_options(ntests,npass)
 use infile_utils, only:inopts,open_db_from_file,close_db
 use setstar, only:star_t,set_defaults_star,read_options_star,write_options_star
 use units, only:set_units
 use physcon, only:solarm,solarr
 use io, only:id
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star,defaults
 type(inopts) :: db(6)
 type(inopts), allocatable :: written(:)
 integer :: i,nerr,nfailed(25),iu
 character(len=48) :: filename
 character(len=32), parameter :: tags(*) = [character(len=32) :: &
    'hcp_density_exponent1','hcp_density_exponent1','hcp_density_exponent1', &
    'hcp_nnopt1','hcp_surface_factor1','hcp_hfloor1','hcp_surface_buffer1', &
    'hcp_surface_buffer1','hcp_profile_mass_unit1','hcp_profile_radius_unit1']
 character(len=16), parameter :: values(*) = [character(len=16) :: &
    '1.5','-0.5','not_a_number','0','0','-1','0','-2','0','-1']

 call set_units(dist=solarr,mass=solarm,G=1.d0)
 call set_defaults_star(defaults)
 db(1:5) = [inopts('iprofile1','2',.false.),inopts('np1','2000',.false.), &
            inopts('Mstar1','1*msun',.false.),inopts('Rstar1','1*rsun',.false.), &
            inopts('particle_lattice1','weighted_hcp',.false.)]
 nfailed = 0
 star = defaults
 nerr = 0
 call read_options_star(star,db(1:5),nerr,label='1')
 call checkval(nerr,0,0,nfailed(1),'missing optional HCP parameters are accepted')
 call checkval(star%hcp_density_exponent,defaults%hcp_density_exponent,0.,nfailed(2),'missing exponent keeps default')
 call checkval(star%hcp_nnopt,defaults%hcp_nnopt,0,nfailed(3),'missing neighbour target keeps default')
 call checkval(star%hcp_surface_buffer,defaults%hcp_surface_buffer,0.,nfailed(4),'missing surface buffer keeps default')

 do i=1,size(tags)
    star = defaults
    db(6) = inopts(tags(i),values(i),.false.)
    nerr = 0
    call read_options_star(star,db,nerr,label='1')
    call checkval(nerr,1,0,nfailed(4+i),'invalid '//trim(tags(i))//' is a setup error')
 enddo
 do i=0,2
    star = defaults
    db(6) = inopts('hcp_density_exponent1','',.false.)
    write(db(6)%val,*) 0.5*i
    nerr = 0
    call read_options_star(star,db,nerr,label='1')
    call checkval(nerr,0,0,nfailed(15+2*i),'valid HCP exponent is accepted')
    call checkval(star%hcp_density_exponent,0.5*i,0.,nfailed(16+2*i),'valid exponent is retained')
 enddo
 star = defaults
 nerr = 0
 call read_options_star(star,db(1:4),nerr,label='1')
 call checkval(nerr,0,0,nfailed(21),'legacy star options need no HCP parameters')
 call checkval(trim(star%particle_lattice)=='closepacked',.true.,nfailed(22),'legacy lattice default is retained')

 ! Profile-unit constants need a round trip, not the default three decimals.
 star%particle_lattice = 'weighted_hcp'
 star%hcp_profile_mass_unit = 1.9892e33
 star%hcp_profile_radius_unit = 6.9598e10
 write(filename,'("test_hcp_options_",i0,".setup")') id
 open(newunit=iu,file=trim(filename),status='replace')
 call write_options_star(star,iu,label='1')
 close(iu)
 call open_db_from_file(written,trim(filename),23,nerr)
 star = defaults
 if (nerr==0) call read_options_star(star,written,nerr,label='1')
 call checkval(nerr,0,0,nfailed(23),'written HCP options can be read')
 call checkval(star%hcp_profile_mass_unit,1.9892e33,0.,nfailed(24),'MESA mass unit survives setup rewrite')
 call checkval(star%hcp_profile_radius_unit,6.9598e10,0.,nfailed(25),'MESA radius unit survives setup rewrite')
 call close_db(written)
 open(newunit=iu,file=trim(filename),status='old')
 close(iu,status='delete')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_hcp_options

subroutine test_composition_without_mu(ntests,npass)
 use part, only:init_part,xyzh,eos_vars,particle_mass,variable_mass,imu,iX,iZ
 use setstar_utils, only:set_star_composition
 use eos, only:get_mean_molecular_weight
 real, parameter :: tol = max(1.e-12,64.*epsilon(1.))
 integer, intent(inout) :: ntests,npass
 integer :: i,nfailed(4)
 real :: mtab(3),Xfrac(3),Yfrac(3),mu(3),expected

 call init_part()
 variable_mass = .true.
 xyzh(:,1:3) = 0.
 do i=1,3
    xyzh(1,i) = 0.1*i
 enddo
 particle_mass(1:3) = 1./3.
 mtab = [0.,0.5,1.]
 Xfrac = 0.7
 Yfrac = 0.28
 mu = 0.
 nfailed = 0
 call set_star_composition(.false.,3,xyzh,Xfrac,Yfrac,mu,mtab,eos_vars)
 expected = get_mean_molecular_weight(0.7,0.02)
 call checkval(maxval(abs(eos_vars(imu,1:3)-expected)),0.,tol,nfailed(1),'missing mu uses ionized composition')
 call checkval(maxval(abs(eos_vars(iX,1:3)-0.7)),0.,tol,nfailed(2),'hydrogen abundance retained')
 call checkval(maxval(abs(eos_vars(iZ,1:3)-0.02)),0.,tol,nfailed(3),'metallicity retained')
 mu = 0.61
 call set_star_composition(.false.,3,xyzh,Xfrac,Yfrac,mu,mtab,eos_vars)
 call checkval(maxval(abs(eos_vars(imu,1:3)-0.61)),0.,tol,nfailed(4),'provided mu takes precedence')
 call update_test_scores(ntests,nfailed,npass)
 call init_part()
end subroutine test_composition_without_mu

subroutine test_analytic_hcp(ntests,npass)
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 use io, only:id,master
 use part, only:init_part,npart,npartoftype,massoftype,xyzh,vxyzu,eos_vars,rad,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,particle_mass,variable_mass,hfact,rho
 use setstar, only:star_t,set_defaults_star,set_star,shift_star,ipoly
 use setstar_utils, only:set_star_thermalenergy
 use stellar_hcp, only:natural_spline,spline_value
 use table_utils, only:yinterp
 use eos, only:gamma,polyk,init_eos
 use options, only:ieos,iexternalforce
 use ptmass, only:icreate_sinks
 use gravwaveutils, only:calc_gravitwaves
 use mpidomain, only:i_belong
 use units, only:set_units
 use physcon, only:solarr,solarm
 use relaxstar, only:maxits,tol_ekin
 use units, only:utime
 use step_lf_global, only:init_step,step
 use timestep, only:dtmax,dtcourant,dtforce
 use damping, only:idamp,damp
 use deriv, only:get_derivs_global
 use energies, only:compute_energies,etot
 use centreofmass, only:get_total_angular_momentum
 use checksetup, only:check_setup
 use kernel, only:hfact_default,radkern,kernelname
 real, parameter :: tol = max(1.e-12,64.*epsilon(1.))
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star,companion
 integer :: ierr,nfailed(8),ieosprev,nfailed_pressure(4),maxits_saved,nfailed_relax(5),nsteps,nfailed_free(7),maxpenetrated
 integer(kind=8) :: ntot
 integer :: nfirst,nerror,nwarn,nfailed_pair(7),sinkcreate_saved,external_saved
 logical :: gravitwaves_saved
 real :: rhozero,gammaprev,polykprev,hfactprev
 real, parameter :: rtest(4) = [0.,1./3.,2./3.,1.],ptest(4) = [1.,0.1,0.01,0.001]
 real :: p2test(4),expected(2)
 real :: residual,omega,velocity_error
 real :: t,dt,dtnew,dtext,energy0,energy_error,radius0,radius_error,L0(3),L(3)
 real, allocatable :: masses(:),initial_radius(:),initial_support(:)

 ieosprev = ieos
 gammaprev = gamma
 polykprev = polyk
 hfactprev = hfact
 sinkcreate_saved = icreate_sinks
 external_saved = iexternalforce
 gravitwaves_saved = calc_gravitwaves
 icreate_sinks = 0
 iexternalforce = 0
 calc_gravitwaves = .false.
 call init_part()
 call set_units(dist=solarr,mass=solarm,G=1.d0)
 npart = 0
 npartoftype = 0
 nptmass = 0
 gamma = 5./3.
 ieos = 2
 call init_eos(ieos,ierr)
 hfact = hfact_default
 call set_defaults_star(star)
 star%iprofile = ipoly
 star%m = '1*msun'
 star%r = '1*rsun'
 star%np = 2000
 star%particle_lattice = 'weighted_hcp'
 if (kernelname/='Wendland 2/3D C^4') star%particle_lattice = 'closepacked'
 stellar_models: block
 call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,2,gamma,0.7,0.02,.false.,.false.,.false.,&
               rhozero,ntot,i_belong,ierr,write_files=.false.)
 if (kernelname/='Wendland 2/3D C^4') then
    nfailed = 0
    call checkval(ierr,0,0,nfailed(1),'unsupported individual-mass kernel still allows legacy stellar setup')
    call checkval(npart>0,.true.,nfailed(2),'legacy stellar setup creates particles')
    call checkval(variable_mass,.false.,nfailed(3),'legacy stellar setup keeps equal-mass mode')
    call check_setup(nerror,nwarn)
    call checkval(nerror,0,0,nfailed(4),'legacy equal-mass setup keeps its kernel')
    particle_mass(1:npart) = massoftype(1)
    variable_mass = .true.
    call check_setup(nerror,nwarn)
    call checkval(nerror,1,0,nfailed(5),'unsupported kernel blocks individual-mass run')
    call check_setup(nerror,nwarn,restart=.true.)
    call checkval(nerror,1,0,nfailed(6),'unsupported kernel blocks individual-mass restart')
    call update_test_scores(ntests,nfailed,npass)
    exit stellar_models
 endif
 nfailed = 0
 call checkval(ierr,0,0,nfailed(1),'analytic weighted HCP setup')
 call checkval(variable_mass,.true.,nfailed(2),'analytic HCP enables individual masses')
 call checkval(npart > 0,.true.,nfailed(3),'analytic HCP contains particles')
 if (npart > 0) then
    call checkval(sum(particle_mass(1:npart)),1.,tol,nfailed(4),'analytic HCP stellar mass')
    call checkval(all(ieee_is_finite(xyzh(:,1:npart))),.true.,nfailed(5),'finite analytic HCP geometry')
    call checkval(all(particle_mass(1:npart) > 0.),.true.,nfailed(6),'positive analytic HCP masses')
    call checkval(all(ieee_is_finite(vxyzu(:,1:npart))),.true.,nfailed(7),'finite analytic HCP thermal state')
    call checkval(all(vxyzu(4,1:npart) > 0.),.true.,nfailed(8),'positive analytic HCP internal energy')
 endif
 call update_test_scores(ntests,nfailed,npass)

 ! Exercise the upstream asynchronous relaxer with fixed unequal masses and spin.
 masses = particle_mass(1:npart)
 icreate_sinks = 1
 call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,2,gamma,0.7,0.02,.true.,.false.,.false.,&
               rhozero,ntot,i_belong,ierr,write_files=.false.,energy_error=residual)
 nfailed = 0
 call checkval(ierr/=0,.true.,nfailed(1),'hard relaxation failure reaches the setup caller')
 call checkval(residual>tol_ekin,.true.,nfailed(2),'failed relaxation does not report convergence')
 call checkval(all(xyzh(4,1:size(masses))>0.),.true.,nfailed(3),'failed relaxation unmasks the existing star')
 call checkval(all(particle_mass(1:size(masses))==masses),.true.,nfailed(4),'failed relaxation preserves the existing masses')
 call update_test_scores(ntests,nfailed,npass)
 icreate_sinks = 0
 call init_part()
 maxits_saved = maxits
 maxits = 4000
 star%spin_period = '12 hr'
 call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,2,gamma,0.7,0.02,.true.,.false.,.false.,&
               rhozero,ntot,i_belong,ierr,write_files=.false.,energy_error=residual)
 nfailed_relax = 0
 omega = 2.*acos(-1.)*utime/(12.*3600.)
 velocity_error = max(maxval(abs(vxyzu(1,1:npart)+omega*xyzh(2,1:npart))), &
                      maxval(abs(vxyzu(2,1:npart)-omega*xyzh(1,1:npart))))
 call checkval(ierr,0,0,nfailed_relax(1),'rotating unequal-mass star setup')
 call checkval(residual < tol_ekin,.true.,nfailed_relax(2),'unequal-mass asynchronous relaxation converges')
 call checkval(sum(particle_mass(1:npart)),1.,tol,nfailed_relax(3),'relaxation preserves stellar mass')
 call checkval(velocity_error,0.,tol,nfailed_relax(4),'relaxed unequal-mass star has physical spin')
 call checkval(all(particle_mass(1:npart)==masses),.true.,nfailed_relax(5),'relaxation preserves individual masses')
 call update_test_scores(ntests,nfailed_relax,npass)
 maxits = maxits_saved

 ! Release without damping, using the production integrator and adaptive steps.
 t = 0.
 idamp = 0
 damp = 0.
 dtmax = 0.02
 call init_step(npart,t,dtmax)
 call get_derivs_global()
 call compute_energies(t)
 energy0 = etot
 radius0 = sqrt(sum(masses*sum(xyzh(1:3,1:npart)**2,dim=1)))
 initial_radius = sqrt(sum(xyzh(1:3,1:npart)**2,dim=1))
 initial_support = radkern*xyzh(4,1:npart)
 maxpenetrated = 0
 call get_total_angular_momentum(xyzh,vxyzu,npart,L0)
 energy_error = 0.
 radius_error = 0.
 nsteps = 0
 do while (t < 10. .and. nsteps < 20000)
    dt = min(dtmax,dtcourant,dtforce,10.-t)
    dtext = dt
    call step(npart,npart,t,dt,dtext,dtnew)
    t = t + dt
    nsteps = nsteps + 1
    call compute_energies(t)
    energy_error = max(energy_error,abs(etot/energy0-1.))
    radius_error = max(radius_error,abs(sqrt(sum(masses*sum(xyzh(1:3,1:npart)**2,dim=1)))/radius0-1.))
    maxpenetrated = max(maxpenetrated,count(initial_radius-sqrt(sum(xyzh(1:3,1:npart)**2,dim=1))>initial_support))
 enddo
 call get_total_angular_momentum(xyzh,vxyzu,npart,L)
 nfailed_free = 0
 call checkval(t,10.,tol,nfailed_free(1),'released unequal-mass star reaches requested time')
 call checkval(energy_error,0.,1.e-3,nfailed_free(2),'released unequal-mass stellar energy drift')
 call checkval(radius_error,0.,0.1,nfailed_free(3),'released unequal-mass stellar RMS radius variation')
 call checkval(sqrt(sum((L-L0)**2)/sum(L0**2)),0.,1.e-6,nfailed_free(4),'released stellar angular momentum drift')
 call checkval(maxpenetrated,0,0,nfailed_free(5),'no deep inward penetration in released unequal-mass star')
 call checkval(npart,size(masses),0,nfailed_free(6),'released stellar particle count stays fixed')
 call checkval(all(particle_mass(1:npart)==masses),.true.,nfailed_free(7),'released stellar particle masses stay fixed')
 call update_test_scores(ntests,nfailed_free,npass)

 ! The earlier star is temporarily masked and tagged while its companion relaxes.
 call init_part()
 maxits = 4000
 residual = huge(residual)
 call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,2,gamma,0.7,0.02,.true.,.false.,.false.,&
               rhozero,ntot,i_belong,ierr,itype=1,write_files=.false.,energy_error=residual)
 nfailed_pair = 0
 call checkval(residual < tol_ekin,.true.,nfailed_pair(1),'primary in unequal-mass pair relaxes')
 nfirst = npart
 masses = particle_mass(1:nfirst)
 companion = star
 companion%m = '0.3*msun'
 companion%r = '0.3*rsun'
 companion%np = 1000
 companion%spin_period = '0 hr'
 residual = huge(residual)
 call set_star(id,master,companion,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,2,gamma,0.7,0.02,.true.,.false.,.false.,&
               rhozero,ntot,i_belong,ierr,itype=2,write_files=.false.,energy_error=residual)
 call checkval(residual < tol_ekin,.true.,nfailed_pair(2),'companion relaxes with tagged primary masked')
 call checkval(all(particle_mass(1:nfirst)==masses),.true.,nfailed_pair(3),'companion relaxation preserves primary masses')
 call checkval(sum(particle_mass(nfirst+1:npart)),0.3,tol,nfailed_pair(4),'companion retains its independent stellar mass')
 call shift_star(npart,npartoftype,xyzh,vxyzu,x0=[-2.,0.,0.],v0=[0.,0.,0.],itype=1)
 call shift_star(npart,npartoftype,xyzh,vxyzu,x0=[2.,0.,0.],v0=[0.,0.,0.],itype=2)
 call check_setup(nerror,nwarn)
 call checkval(nerror,0,0,nfailed_pair(5),'released unequal-mass pair passes setup checks')
 call checkval(all(xyzh(4,1:npart)>0.),.true.,nfailed_pair(6),'both stars unmasked after relaxation')
 call checkval(maxval(abs(vxyzu(1:3,nfirst+1:npart))),0.,tol,nfailed_pair(7),'static companion has no inherited primary spin')
 call update_test_scores(ntests,nfailed_pair,npass)
 maxits = maxits_saved

 ! Positive, steeply falling pressure samples can produce a negative cubic interpolant.
 nfailed_pressure = 0
 call natural_spline(rtest,ptest,p2test)
 expected = [spline_value(rtest,ptest,p2test,1./6.),spline_value(rtest,ptest,p2test,0.5)]
 call checkval(expected(1)>0. .and. expected(2)<0.,.true.,nfailed_pressure(1),'pressure spline overshoot fixture')
 expected(2) = yinterp(ptest,rtest,0.5)
 variable_mass = .true.
 xyzh(:,1:2) = 0.
 xyzh(1,1) = 1./6.
 xyzh(1,2) = 0.5
 rho(1:2) = 1.
 call set_star_thermalenergy(2,[1.,1.,1.,1.],ptest,[0.,0.,0.,0.],rtest,4,2,xyzh,vxyzu,rad,eos_vars,&
                            .false.,.false.,0.,0.)
 call checkval(vxyzu(4,1),expected(1)/(gamma-1.),tol,nfailed_pressure(2),'positive pressure keeps cubic interpolation')
 call checkval(vxyzu(4,2),expected(2)/(gamma-1.),tol,nfailed_pressure(3),'negative pressure uses linear interpolation')
 call checkval(all(ieee_is_finite(vxyzu(4,1:2))).and.all(vxyzu(4,1:2)>0.),.true.,&
               nfailed_pressure(4),'finite positive thermal energy after pressure fallback')
 call update_test_scores(ntests,nfailed_pressure,npass)
 end block stellar_models
 call init_part()
 ieos = ieosprev
 gamma = gammaprev
 polyk = polykprev
 hfact = hfactprev
 icreate_sinks = sinkcreate_saved
 iexternalforce = external_saved
 calc_gravitwaves = gravitwaves_saved

end subroutine test_analytic_hcp

!-----------------------------------------------------------------------
!+
!   Check that shift_star with corotate puts a star into rigid rotation
!   with its orbit: on a circular orbit v = omega x r for every particle,
!   so the star spins the same way as it orbits, and on an eccentric orbit
!   the star keeps its full orbital velocity (including the radial part),
!   v = v0 + omega x (r - x0)
!+
!-----------------------------------------------------------------------
subroutine test_shift_star_corotate(ntests,npass)
 use setstar,     only:shift_star
 use vectorutils, only:cross_product3D
 integer, intent(inout) :: ntests,npass
 integer, parameter :: np = 3
 integer :: i,npartoftype(1),nfail(2)
 real    :: xyz0(3,np),xyz(3,np),vxyz(3,np),vexact(3,np),x0(3),v0(3),omega(3),dv(3)

 ! particles at the centre and around it, initially at rest
 xyz0(:,1) = (/0.,0.,0./)
 xyz0(:,2) = (/0.1,0.,0./)
 xyz0(:,3) = (/0.,0.2,-0.1/)
 npartoftype = np
 xyz  = xyz0
 vxyz = 0.

 ! circular orbit at radius 2 with speed 1, so omega = (0,0,1/2)
 x0 = (/2.,0.,0./)
 v0 = (/0.,1.,0./)
 omega = (/0.,0.,0.5/)
 call shift_star(np,npartoftype,xyz,vxyz,x0,v0,corotate=.true.)
 do i=1,np
    call cross_product3D(omega,xyz(:,i),vexact(:,i))
 enddo
 call checkval(3*np,reshape(vxyz,(/3*np/)),reshape(vexact,(/3*np/)),epsilon(0.),nfail(1),'v = omega x r for corotating star')

 ! eccentric orbit: same position and angular momentum, so the same omega, plus a radial velocity
 xyz  = xyz0
 vxyz = 0.
 v0 = (/0.3,1.,0./)
 call shift_star(np,npartoftype,xyz,vxyz,x0,v0,corotate=.true.)
 do i=1,np
    call cross_product3D(omega,xyz(:,i)-x0,dv)
    vexact(:,i) = v0 + dv
 enddo
 call checkval(3*np,reshape(vxyz,(/3*np/)),reshape(vexact,(/3*np/)),epsilon(0.),nfail(2),&
               'v = v0 + omega x (r - x0) on eccentric orbit')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_shift_star_corotate

!-----------------------------------------------------------------------
!+
!   Check that the procedure to get the Lagrangian mass coordinate
!   from the radius works
!+
!-----------------------------------------------------------------------
subroutine test_get_mass_coord(ntests,npass)
 use sortutils,     only:r2func,find_rank
 use part,          only:xyzh,massoftype,igas
 use setstar_utils, only:get_mass_coord
 use table_utils,   only:linspace
 use testutils,     only:checkvalbuf,checkvalbuf_end
 integer, intent(inout) :: ntests,npass
 real,    allocatable   :: mass_enclosed_r(:)
 integer, allocatable :: iorder(:)
 integer :: i,i1,itest,np,nfail(1),ncheck
 real    :: massri,errmax,x0(3)
 real, parameter :: tol = 1.e-14

 i1 = 10  ! start at non-zero offset to check this
 np = 100
 x0 = 0.

 ! place the particles in a line between x=[0,1]
 xyzh(:,1:np) = 0.
 call linspace(xyzh(1,i1+1:np),0.,1.)

 ! give them equal masses
 massoftype(igas) = 3.e-6

 ! call the routine we are trying to test
 call get_mass_coord(i1,np,xyzh,mass_enclosed_r,x0)

 ! check memory was allocated correctly
 call checkval(size(mass_enclosed_r),np-i1,0,nfail(1),'size of mass_enclosed_r array')
 call update_test_scores(ntests,nfail,npass)

 ! check that the mass ranking is correct
 ncheck = 0; nfail = 0; errmax = 0.
 do i=i1+1,np
    call checkvalbuf(mass_enclosed_r(i-i1),(i-i1-1)*massoftype(igas),&
                     tol,'m(x) for particles in x=[0,1]',nfail(1),ncheck,errmax)
 enddo
 call checkvalbuf_end('m(x) for particles in x=[0,1]',ncheck,nfail(1),errmax,tol)
 call update_test_scores(ntests,nfail,npass)

 do itest=1,2
    if (itest==2) then
       ! now check pathological case where half the particles are at y=0, and half are at y=1
       xyzh(:,1:np) = 0.
       xyzh(2,i1+np/2:np) = 1.
       call get_mass_coord(i1,np,xyzh,mass_enclosed_r,x0)
    endif
    !
    ! check that this agrees with our previous method that
    ! works only for equal mass particles
    !
    mass_enclosed_r = mass_enclosed_r / ((np-i1)*massoftype(igas))  ! make it a mass fraction
    call find_rank(np-i1,r2func,xyzh(:,i1+1:np),iorder)

    ncheck = 0; nfail = 0; errmax = 0.
    do i=i1+1,np
       massri = real(iorder(i-i1)-1) / real(np-i1)
       call checkvalbuf(mass_enclosed_r(i-i1),massri,tol,'m(<r) agrees with previous method',nfail(1),ncheck,errmax)
    enddo
    call checkvalbuf_end('m(<r) same as previous method',ncheck,nfail(1),errmax,tol)
    call update_test_scores(ntests,nfail,npass)
 enddo

end subroutine test_get_mass_coord

!-----------------------------------------------------------------------
!+
!   Check signed spin periods and backwards-compatible setup files
!+
!-----------------------------------------------------------------------
subroutine test_spin_options(ntests,npass)
 use infile_utils, only:inopts
 use setstar, only:star_t,set_defaults_star,read_options_star
 use units, only:set_units,in_code_units
 use physcon, only:solarm,solarr
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 type(inopts) :: db(5)
 integer :: i,nerr,nfailed(2)
 character(len=20), parameter :: periods(*) = [character(len=20) :: &
    '24 hr','-24 hr','0 hr','24 rsun','invalid','NaN hr','Inf hr']

 call set_units(dist=solarr,mass=solarm,G=1.d0)
 db(1:4) = [inopts('iprofile1','2',.false.),inopts('np1','1000',.false.), &
            inopts('Mstar1','1*msun',.false.),inopts('Rstar1','1*rsun',.false.)]
 call set_defaults_star(star)
 nerr = 0
 call read_options_star(star,db(1:4),nerr,label='1')
 call checkval(nerr,0,0,nfailed(1),'legacy setup needs no spin period')
 call checkval(in_code_units(star%spin_period,nerr),0.,0.,nfailed(2),'default star has no spin')
 call update_test_scores(ntests,nfailed,npass)
 do i=1,size(periods)
    call set_defaults_star(star)
    db(5) = inopts('spin_period1',periods(i),.false.)
    nerr = 0
    call read_options_star(star,db,nerr,label='1')
    call checkval(nerr,merge(0,1,i<=3),0,nfailed(1),'spin period validation: '//trim(periods(i)))
    call update_test_scores(ntests,nfailed(1:1),npass)
 enddo

end subroutine test_spin_options

!-----------------------------------------------------------------------
!+
!   test that we can successfully relax a polytrope
!   and that the profile matches the analytic solution
!+
!-----------------------------------------------------------------------
subroutine test_polytrope(ntests,npass)
 use io,        only:id,master,iverbose
 use part,      only:npart,npartoftype,xyzh,vxyzu,eos_vars,rad,massoftype,hfact,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass
 use mpidomain, only:i_belong
 use options,   only:ieos
 use physcon,   only:solarr,solarm,pi
 use eos,       only:gamma,X_in,Z_in,polyk
 use setstar,   only:star_t,set_star,set_defaults_star,ipoly
 use units,     only:set_units,in_code_units
 use checksetup, only:check_setup
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 real :: rhozero,rmserr,ekin,x0(3),v0(3),omega,velocity_error,jz
 integer(kind=8) :: ntot
 integer :: ierr,nfail(1),i,j,i1,nerror,nwarn

 npart = 0
 npartoftype = 0
 massoftype = 0.
 iverbose = 0
 call set_units(dist=10.*solarr,mass=10.*solarm,G=1.d0)
 ieos = 2
 gamma = 5./3.
 call set_defaults_star(star)
 star%m = '1*msun'
 star%r = '1*rsun'
 star%iprofile = ipoly  ! a polytrope
 star%np = 1000
 x0 = 0.
 v0 = 0.
 ! Also relax prograde and retrograde stars after an existing star.
 do i=1,3
    if (i > 1) then
       x0 = [3.*(i-1),0.,0.]
       v0 = [0.2,-0.1,0.3]
       star%spin_period = '24 hr'
       if (i==3) star%spin_period = '-24 hr'
    endif
    omega = 0.
    if (i > 1) omega = 2.*pi/in_code_units(star%spin_period,ierr)
    i1 = npart

    call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
               npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
               relax=.true.,use_var_comp=.false.,write_rho_to_file=.false.,&
               rhozero=rhozero,npart_total=ntot,mask=i_belong,ierr=ierr,&
               write_files=.false.,density_error=rmserr,energy_error=ekin,x0=x0,v0=v0)

    call checkval(ierr,0,0,nfail(1),'set_star runs with ierr = 0')
    call update_test_scores(ntests,nfail,npass)

    call check_setup(nerror,nwarn,restart=.false.)
    call checkval(nerror+nwarn,0,0,nfail(1),'no errors or warnings')
    call update_test_scores(ntests,nfail,npass)

    call checkval(rhozero,100./(4./3.*pi),1e-6,nfail(1),'mean density')
    call update_test_scores(ntests,nfail,npass)

    call checkval(star%polyk,1.9694457e-2,1e-6,nfail(1),'polyk value for M=0.1,R=0.1')
    call update_test_scores(ntests,nfail,npass)

    call checkval(polyk,1.9694457e-2,1e-6,nfail(1),'polyk value for M=0.1,R=0.1')
    call update_test_scores(ntests,nfail,npass)

    call checkval(rmserr,0.0,0.04,nfail(1),'error in density profile')
    call update_test_scores(ntests,nfail,npass)

    call checkval(ekin,0.,1e-7,nfail(1),'ekin/epot < 1.e-7')
    call update_test_scores(ntests,nfail,npass)
    velocity_error = 0.
    jz = 0.
    do j=i1+1,npart
       velocity_error = max(velocity_error,abs(vxyzu(1,j)-v0(1)+omega*(xyzh(2,j)-x0(2))), &
                            abs(vxyzu(2,j)-v0(2)-omega*(xyzh(1,j)-x0(1))),abs(vxyzu(3,j)-v0(3)))
       jz = jz + (xyzh(1,j)-x0(1))*(vxyzu(2,j)-v0(2)) - (xyzh(2,j)-x0(2))*(vxyzu(1,j)-v0(1))
    enddo
    call checkval(velocity_error,0.,1.e-12,nfail(1),'physical spin and bulk velocity after relaxation')
    call update_test_scores(ntests,nfail,npass)
    if (i > 1) then
       call checkval(jz*omega > 0.,.true.,nfail(1),'spin angular momentum has requested sign')
       call update_test_scores(ntests,nfail,npass)
    endif
 enddo

 call checkval(npart,3*star%np,0,nfail(1),'np = 3000')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_polytrope

!-----------------------------------------------------------------------
!+
!   test that we can successfully relax a red supergiant star from
!   Lau et al. (2022a,b)
!+
!-----------------------------------------------------------------------
subroutine test_redsupergiant(ntests,npass)
 use eos,     only:gmw,X_in,Z_in,irecomb
 use setstar, only:star_t,set_defaults_star
 integer, intent(inout) :: ntests,npass
 integer, parameter :: ieos_list(4) = (/2,12,10,20/)  ! ideal gas, gas+rad, MESA, gas+rad+recombination
 type(star_t) :: star

 gmw     = 0.61821
 X_in    = 0.698
 Z_in    = 1.-X_in-0.287
 irecomb = 0
 call set_defaults_star(star)
 star%isinkcore = .true.
 star%hsoft     = '9.25'
 star%hacc      = '0.0'
 star%mcore     = '3.8405'
 call test_relax_mesa_star(ntests,npass,'lau22.data','red supergiant',star,ieos_list,100,0.04)

end subroutine test_redsupergiant

!-----------------------------------------------------------------------
!+
!   test that we can successfully relax a white dwarf with the
!   Helmholtz and zero temperature equations of state (Ali Pourmand)
!+
!-----------------------------------------------------------------------
subroutine test_whitedwarf(ntests,npass)
 use eos,     only:X_in,Z_in
 use setstar, only:star_t,set_defaults_star
 integer, intent(inout) :: ntests,npass
 integer, parameter :: ieos_list(2) = (/15,25/)  ! Helmholtz, zero temperature
 type(star_t) :: star

 X_in = 0.
 Z_in = 1.
 call set_defaults_star(star)
 call test_relax_mesa_star(ntests,npass,'pourmandwd.data','white dwarf',star,ieos_list,1500,0.05)

end subroutine test_whitedwarf

!-----------------------------------------------------------------------
!+
!   set up and relax a star from a MESA profile with each of the
!   requested equations of state, and check the relaxed star matches
!   the profile. Skipped if the profile is not found.
!+
!-----------------------------------------------------------------------
subroutine test_relax_mesa_star(ntests,npass,profile,label,star,ieos_list,maxits_in,tol_relax)
 use io,        only:id,master,iverbose,warning
 use dim,       only:do_radiation
 use datafiles, only:find_phantom_datafile
 use part,      only:init_part,npart,npartoftype,xyzh,vxyzu,eos_vars,rad,massoftype,hfact,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass,rho,iradxi,itemp
 use mpidomain, only:i_belong
 use options,   only:ieos
 use physcon,   only:solarr,solarm
 use eos,       only:init_eos,gamma,X_in,Z_in,eos_works_with_radiation
 use setstar,   only:star_t,set_star,imesa,ierr_radiation_conflict
 use units,     only:set_units
 use relaxstar, only:maxits
 use checksetup,  only:check_setup
 use testutils,   only:checkvalbuf,checkvalbuf_end
 use table_utils, only:yinterp
 use readwrite_mesa,  only:read_mesa
 use radiation_utils, only:Trad_from_radxi
 integer,          intent(inout) :: ntests,npass
 character(len=*), intent(in)    :: profile,label
 type(star_t),     intent(inout) :: star
 integer,          intent(in)    :: ieos_list(:),maxits_in
 real,             intent(in)    :: tol_relax
 real, parameter :: tolprof = 0.05
 character(len=500) :: filepath
 real :: rhozero,rmserr,ekin,errmax(2),Mstar,rhoj,rhoj_mesa,rhomax,rj,Tgas,Trad
 integer(kind=8) :: ntot
 integer :: ierr,nfail(2),ncheck(2),i,j,nerror,nwarn,iunit,expected_error
 real, allocatable :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),mu(:)

 filepath = find_phantom_datafile(profile,'star_data_files')
 open(newunit=iunit,file=filepath,status="old",iostat=ierr)
 if (ierr /= 0) then
    call warning('test_setstar',trim(profile)//' not found, skipping '//trim(label)//' test')
    return
 else
    close(unit=iunit)
 endif

 iverbose = 0
 call set_units(dist=solarr,mass=solarm,G=1.d0)
 gamma = 5./3.
 maxits = maxits_in
 star%np            = 10000
 star%input_profile = trim(filepath)
 star%iprofile      = imesa

 nfail = 0; ncheck = 0; errmax = 0.
 call read_mesa(filepath,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)

 do i=1,size(ieos_list)
    call init_part()
    ieos = ieos_list(i)
    if (id==master) write(*,"(/,a,i2/)") '--> Testing '//trim(label)//' with ieos = ',ieos
    call init_eos(ieos,ierr)
    if (ierr /= 0) then  ! e.g. eos tables not available
       call warning('test_setstar','could not initialise eos, skipping '//trim(label)//' test',var='ieos',ival=ieos)
       cycle
    endif

    rmserr = 0.
    call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
               npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
               relax=.true.,use_var_comp=.false.,write_rho_to_file=.false.,&
               rhozero=rhozero,npart_total=ntot,mask=i_belong,ierr=ierr,&
               write_files=.false.,density_error=rmserr,energy_error=ekin)

    if (do_radiation .and. (.not. eos_works_with_radiation(ieos))) then
       expected_error = ierr_radiation_conflict
    else
       expected_error = 0
    endif
    call checkval(ierr,expected_error,0,nfail(1),'set_star runs with expected ierr')
    call update_test_scores(ntests,nfail,npass)
    if (ierr == ierr_radiation_conflict) cycle  ! skip the rest of the checks if we get a radiation conflict error

    call check_setup(nerror,nwarn,restart=.false.)
    call checkval(nerror+nwarn,0,0,nfail(1),'no errors or warnings')
    call update_test_scores(ntests,nfail,npass)

    call checkval(rmserr,0.0,tol_relax,nfail(1),'error in density profile')
    call update_test_scores(ntests,nfail,npass)

    call checkval(ekin,0.,5e-6,nfail(1),'ekin/epot < 1.e-7')
    call update_test_scores(ntests,nfail,npass)

    rmserr = 0.
    rhomax = 0.
    do j=1,npart
       rj = sqrt(dot_product(xyzh(1:3,j),xyzh(1:3,j)))
       rhoj = rho(j)
       rhoj_mesa = yinterp(den,r,rj)
       rmserr = rmserr + (rhoj - rhoj_mesa)**2
       rhomax = max(rhomax,rhoj_mesa)

       if (do_radiation) then
          Trad = Trad_from_radxi(rhoj, rad(iradxi,j))
          Tgas = eos_vars(itemp,j)
          call checkvalbuf(Trad/Tgas,1.,1e-15,'Trad/Tgas = 1 for radiative star',nfail(2),ncheck(2),errmax(2),use_rel_tol=.true.)
       endif
    enddo
    ! L2 error norm normalised by the maximum density, as in splash (a relative error
    ! per particle is dominated by the unresolved low density surface layers)
    rmserr = sqrt(rmserr/npart)/rhomax
    call checkval(rmserr,0.,tolprof,nfail(1),'density matches MESA profile')
    call update_test_scores(ntests,nfail,npass)

    if (do_radiation) then
       call checkvalbuf_end('Trad/Tgas ~ 1 for radiative star',ncheck(2),nfail(2),errmax(2),1.e-15)
       call update_test_scores(ntests,nfail,npass)
    endif
 enddo

end subroutine test_relax_mesa_star

end module testsetstar
