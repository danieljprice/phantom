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
! :Dependencies: checksetup, datafiles, dim, eos, io, mpidomain, options,
!   part, physcon, radiation_utils, readwrite_mesa, relaxstar, setstar,
!   setstar_utils, sortutils, stellar_hcp, table_utils, testutils, units
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
 use dim,        only:gravity,do_radiation,mpi,use_apr,ind_timesteps
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
 call test_weighted_hcp(ntests,npass)
 call test_hcp_surface_buffer(ntests,npass)
 call test_composition_without_mu(ntests,npass)
 if (.not.(do_radiation.or.mpi.or.use_apr.or.ind_timesteps)) call test_analytic_hcp(ntests,npass)
 if (present(string)) then
    if (string=='hcp'.or.string=='weighted_hcp') return
 endif

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
               xyzmh_ptmass,vxyz_ptmass,nptmass,particle_mass,variable_mass,hfact
 use setstar, only:star_t,set_defaults_star,set_star,ipoly
 use eos, only:gamma,polyk,init_eos
 use options, only:ieos
 use mpidomain, only:i_belong
 use units, only:set_units
 use physcon, only:solarr,solarm
 real, parameter :: tol = max(1.e-12,64.*epsilon(1.))
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 integer :: ierr,nfailed(8),ieosprev
 integer(kind=8) :: ntot
 real :: rhozero,gammaprev,polykprev,hfactprev

 ieosprev = ieos
 gammaprev = gamma
 polykprev = polyk
 hfactprev = hfact
 call init_part()
 call set_units(dist=solarr,mass=solarm,G=1.d0)
 npart = 0
 npartoftype = 0
 nptmass = 0
 gamma = 5./3.
 ieos = 2
 call init_eos(ieos,ierr)
 hfact = 1.2
 call set_defaults_star(star)
 star%iprofile = ipoly
 star%m = '1*msun'
 star%r = '1*rsun'
 star%np = 2000
 star%particle_lattice = 'weighted_hcp'
 call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,2,gamma,0.7,0.02,.false.,.false.,.false.,&
               rhozero,ntot,i_belong,ierr,write_files=.false.)
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
 call init_part()
 ieos = ieosprev
 gamma = gammaprev
 polyk = polykprev
 hfact = hfactprev

end subroutine test_analytic_hcp

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
 use units,     only:set_units
 use checksetup, only:check_setup
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 real :: rhozero,rmserr,ekin,x0(3)
 integer(kind=8) :: ntot
 integer :: ierr,nfail(1),i,nerror,nwarn

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
 ! do this test twice, to check the second star relaxes...
 do i=1,2
    if (i==2) x0 = [3.,0.,0.]

    call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
               npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
               relax=.true.,use_var_comp=.false.,write_rho_to_file=.false.,&
               rhozero=rhozero,npart_total=ntot,mask=i_belong,ierr=ierr,&
               write_files=.false.,density_error=rmserr,energy_error=ekin,x0=x0)

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
 enddo

 call checkval(npart,2*star%np,0,nfail(1),'np = 2000')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_polytrope

!-----------------------------------------------------------------------
!+
!   test that we can successfully relax a red supergiant star from
!   Lau et al. (2022a,b)
!+
!-----------------------------------------------------------------------
subroutine test_redsupergiant(ntests,npass)
 use io,        only:id,master,iverbose,warning
 use dim,       only:do_radiation
 use datafiles, only:find_phantom_datafile
 use part,      only:init_part,npart,npartoftype,xyzh,vxyzu,eos_vars,rad,massoftype,hfact,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass,rho,igas,igasP,imu,iX,iZ,iradxi,&
                     eos_vars,itemp
 use mpidomain, only:i_belong
 use options,   only:ieos
 use physcon,   only:solarr,solarm
 use eos,       only:init_eos,gamma,gmw,X_in,Z_in,irecomb,eos_works_with_radiation
 use setstar,   only:star_t,set_star,set_defaults_star,imesa,ierr_radiation_conflict
 use units,     only:set_units
 use relaxstar, only:maxits
 use checksetup,  only:check_setup
 use testutils,   only:checkvalbuf,checkvalbuf_end
 use table_utils, only:yinterp
 use readwrite_mesa,  only:read_mesa
 use radiation_utils, only:Trad_from_radxi
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 character(len=500) :: filepath
 real :: rhozero,rmserr,rmserr_mu,rmserr_X,rmserr_Z,ekin,x0(3),errmax(2)
 real :: Mstar,tolprof,tolcomposition,rhoj,rhoj_mesa,rj,Tgas,Trad
 integer(kind=8) :: ntot
 integer :: ierr,nfail(2),ncheck(2),i,j,nerror,nwarn,iunit,expected_error
 logical :: relax_star,var_comp
 real, allocatable :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),Zfrac(:),mu(:)

 filepath = find_phantom_datafile('lau22.data','star_data_files')
 open(newunit=iunit,file=filepath,status="old",iostat=ierr)
 if (ierr /= 0) then
    call warning('test_setstar','lau22.data not found, skipping red supergiant test')
    return
 else
    close(unit=iunit)
 endif

 iverbose = 0
 call set_units(dist=solarr,mass=solarm,G=1.d0)
 gamma = 5./3.
 gmw = 0.61821
 X_in = 0.698
 Z_in = 1.-X_in-0.287
 irecomb = 0
 relax_star = .true.
 maxits = 100
 var_comp = .false.
 call set_defaults_star(star)
 star%isinkcore   = .true.
 star%hsoft          = '9.25'
 star%hacc           = '0.0'
 star%mcore          = '3.8405'
 star%np             = 10000
 star%input_profile  = trim(filepath)
 star%iprofile = imesa
 x0 = 0.

 nfail = 0; ncheck = 0; errmax = 0.
 tolprof = 0.08
 tolcomposition = 0.08

 call read_mesa(filepath,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)
 Zfrac = 1.-Xfrac-Yfrac

 do i=1,4
    call init_part()
    if (i==1) then  ! ideal gas EoS
       ieos = 2
    elseif (i==2) then  ! ideal gas + radiation EoS
       ieos = 12
    elseif (i==3) then  ! MESA EoS
       ieos = 10
    elseif (i==4) then  ! Gas + radiation + recombination energy EoS
       ieos = 20
    elseif (i==5) then  ! ideal gas + radiation EoS, with variable composition
       ieos = 12
       var_comp = .true.
    endif

    if (id==master) write(*,"(/,a,i2/)") '--> Testing ieos = ',ieos
    call init_eos(ieos,ierr)

    rmserr = 0.
    call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
               npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
               relax=relax_star,use_var_comp=var_comp,write_rho_to_file=.false.,&
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

    call checkval(rmserr,0.0,0.04,nfail(1),'error in density profile')
    call update_test_scores(ntests,nfail,npass)

    call checkval(ekin,0.,5e-6,nfail(1),'ekin/epot < 1.e-7')
    call update_test_scores(ntests,nfail,npass)

    rmserr = 0.
    rmserr_mu = 0.
    rmserr_X = 0.
    rmserr_Z = 0.
    do j=1,npart
       rj = sqrt(dot_product(xyzh(1:3,j),xyzh(1:3,j)))
       rhoj = rho(j)
       rhoj_mesa = yinterp(den,r,rj)
       rmserr = rmserr + (1-rhoj/rhoj_mesa)**2

       if (do_radiation) then
          Trad = Trad_from_radxi(rhoj, rad(iradxi,j))
          Tgas = eos_vars(itemp,j)
          call checkvalbuf(Trad/Tgas,1.,1e-15,'Trad/Tgas = 1 for radiative star',nfail(2),ncheck(2),errmax(2),use_rel_tol=.true.)
       endif

       if (i==5) then  ! calculate rms error in mu, X, Z, not implemented yet
          rmserr_mu = rmserr_mu + (1-abs(eos_vars(imu,j)/yinterp(mu,r,rj)))**2
          rmserr_X = rmserr_X + (1-abs(eos_vars(iX,j)/yinterp(Xfrac,r,rj)))**2
          rmserr_Z = rmserr_Z + (1-abs(eos_vars(iZ,j)/yinterp(Zfrac,r,rj)))**2
       endif
    enddo
    rmserr = sqrt(rmserr/npart)  ! root mean square relative error
    call checkval(rmserr,0.,tolprof,nfail(1),'density matches MESA profile')
    call update_test_scores(ntests,nfail,npass)

    call checkvalbuf_end('Trad/Tgas ~ 1 for radiative star',ncheck(2),nfail(2),errmax(2),1.e-15)
    call update_test_scores(ntests,nfail,npass)

    if (i==5) then  ! not implemented yet
       rmserr_mu = sqrt(rmserr_mu/npart)
       rmserr_X = sqrt(rmserr_X/npart)
       rmserr_Z = sqrt(rmserr_Z/npart)
       call checkval(rmserr_mu,0.,tolcomposition,nfail(1),'mu matches MESA profile')
       call update_test_scores(ntests,nfail,npass)
       call checkval(rmserr_X,0.,tolcomposition,nfail(1),'X matches MESA profile')
       call update_test_scores(ntests,nfail,npass)
       call checkval(rmserr_Z,0.,tolcomposition,nfail(1),'Z matches MESA profile')
       call update_test_scores(ntests,nfail,npass)
    endif

 enddo

end subroutine test_redsupergiant

!-----------------------------------------------------------------------
!+
!   test that we can successfully relax a white dwarf with the
!   zero temperature EoS and Helmholtz eos(Ali Pourmand, similar to test_redsupergiant)
!+
!-----------------------------------------------------------------------
subroutine test_whitedwarf(ntests,npass)
 use io,        only:id,master,iverbose,warning
 use dim,       only:do_radiation
 use datafiles, only:find_phantom_datafile
 use part,      only:init_part,npart,npartoftype,xyzh,vxyzu,eos_vars,rad,massoftype,hfact,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass,rho,igas,&
                     eos_vars
 use mpidomain, only:i_belong
 use options,   only:ieos
 use physcon,   only:solarr,solarm
 use eos,       only:init_eos,gamma,X_in,Z_in
 use setstar,   only:star_t,set_star,set_defaults_star,imesa
 use units,     only:set_units
 use relaxstar, only:maxits
 use checksetup,  only:check_setup
 use testutils,   only:checkvalbuf,checkvalbuf_end
 use table_utils, only:yinterp
 use readwrite_mesa,  only:read_mesa
 integer, intent(inout) :: ntests,npass
 type(star_t) :: star
 character(len=500) :: filepath
 real :: rhozero,rmserr,ekin,x0(3),errmax(2)
 real :: Mstar,tolprof,rhoj,rhoj_mesa,rj
 integer(kind=8) :: ntot
 integer :: ierr,nfail(2),ncheck(2),i,j,nerror,nwarn,iunit,expected_error
 logical :: relax_star,var_comp
 real, allocatable :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),Zfrac(:),mu(:)

 filepath = find_phantom_datafile('pourmandwd.data','star_data_files')
 open(newunit=iunit,file=filepath,status="old",iostat=ierr)
 if (ierr /= 0) then
    call warning('test_setstar','pourmandwd.data not found, skipping white dwarf test')
    return
 else
    close(unit=iunit)
 endif

 iverbose = 0
 call set_units(dist=solarr,mass=solarm,G=1.d0)
 gamma = 5./3.
 X_in = 0.
 Z_in = 1.
 relax_star = .true.
 maxits = 1500
 var_comp = .false.
 call set_defaults_star(star)
 star%np             = 10000
 star%input_profile  = trim(filepath)
 star%iprofile = imesa
 x0 = 0.

 nfail = 0; ncheck = 0; errmax = 0.
 tolprof = 0.1

 call read_mesa(filepath,den,r,pres,mtab,en,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr)
 Zfrac = 1.-Xfrac-Yfrac

 do i=1,2
    call init_part()
    if (i==1) then  ! EOS Helmholtz
       ieos = 15
    elseif (i==2) then  ! zero Temperature EOS
       ieos = 25
    endif

    if (id==master) write(*,"(/,a,i2/)") '--> Testing ieos = ',ieos
    call init_eos(ieos,ierr)

    rmserr = 0.
    call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
               npart,npartoftype,massoftype,hfact,&
               xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
               relax=relax_star,use_var_comp=var_comp,write_rho_to_file=.false.,&
               rhozero=rhozero,npart_total=ntot,mask=i_belong,ierr=ierr,&
               write_files=.false.,density_error=rmserr,energy_error=ekin)

    call checkval(ierr,expected_error,0,nfail(1),'set_star runs with expected ierr')
    call update_test_scores(ntests,nfail,npass)

    call check_setup(nerror,nwarn,restart=.false.)
    call checkval(nerror+nwarn,0,0,nfail(1),'no errors or warnings')
    call update_test_scores(ntests,nfail,npass)

    call checkval(rmserr,0.0,0.05,nfail(1),'error in density profile')
    call update_test_scores(ntests,nfail,npass)

    call checkval(ekin,0.,5e-6,nfail(1),'ekin/epot < 1.e-7')
    call update_test_scores(ntests,nfail,npass)

    rmserr = 0.
    do j=1,npart
       rj = sqrt(dot_product(xyzh(1:3,j),xyzh(1:3,j)))
       rhoj = rho(j)
       rhoj_mesa = yinterp(den,r,rj)
       rmserr = rmserr + (1-rhoj/rhoj_mesa)**2

    enddo
    rmserr = sqrt(rmserr/npart)  ! root mean square relative error
    call checkval(rmserr,0.,tolprof,nfail(1),'density matches MESA profile')
    call update_test_scores(ntests,nfail,npass)

 enddo

end subroutine test_whitedwarf

end module testsetstar
