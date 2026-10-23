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
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: checksetup, datafiles, dim, eos, infile_utils, io,
!   mpidomain, options, part, physcon, radiation_utils, readwrite_mesa,
!   relaxstar, setstar, setstar_utils, sortutils, table_utils, testutils,
!   units, vectorutils
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
subroutine test_setstar(ntests,npass)
 use checksetup, only:check_setup
 use io,         only:id,master
 use dim,        only:gravity,do_radiation
 integer, intent(inout) :: ntests,npass

 if (.not.gravity) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING STAR SETUP TESTS'
    return
 else
    if (id==master) write(*,"(/,a,/)") '--> TESTING STAR SETUP'
 endif

 ! test the get_mass_coordinate routine in relax_star
 call test_get_mass_coord(ntests,npass)
 call test_spin_options(ntests,npass)

 ! test that a star moved onto a circular orbit is set corotating
 call test_shift_star_corotate(ntests,npass)

 ! test polytrope
 if (.not. do_radiation) call test_polytrope(ntests,npass)

 ! test red supergiant
 call test_redsupergiant(ntests,npass)

 ! test white dwarf
 call test_whitedwarf(ntests,npass)

end subroutine test_setstar

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
