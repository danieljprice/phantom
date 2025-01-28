!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testeos
!
! Unit tests of the equation of state module
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, eos_barotropic, eos_gasradrec, eos_helmholtz,
!   eos_idealplusrad, eos_tillotson, io, ionization_mod, mpiutils, part,
!   physcon, table_utils, testeos_stratified, testutils, units
!
 implicit none
 public :: test_eos
 public :: test_helmholtz ! to avoid compiler warning for unused routine

 private
 logical :: use_rel_tol = .true.

contains
!----------------------------------------------------------
!+
!  unit tests of equation of state module
!+
!----------------------------------------------------------
subroutine test_eos(ntests,npass)
 use io,            only:id,master,stdout
 use physcon,       only:solarm
 use units,         only:set_units
 use eos_gasradrec, only:irecomb
 use testeos_stratified, only:test_eos_stratified
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING EQUATION OF STATE MODULE'

 call set_units(mass=solarm,dist=1.d16,G=1.d0)

 !
 ! perform tests that can be applied to most equations of state
 !
 call test_all(ntests,npass)

 !
 ! unit tests for particular equations of state below
 !
 !  call test_helmholtz(ntests, npass)
 call test_idealplusrad(ntests, npass)

 do irecomb = 0,3
    call test_hormone(ntests,npass)
 enddo

 call test_eos_stratified(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- EQUATION OF STATE TEST COMPLETE'

end subroutine test_eos


!----------------------------------------------------------
!+
!  test that the initialisation of all eos works correctly
!+
!----------------------------------------------------------
subroutine test_all(ntests, npass)
 use eos,       only:maxeos,init_eos,isink,polyk,polyk2,qfacdisc,&
                     ierr_file_not_found,ierr_option_conflict,&
                     eos_is_not_implemented
 use io,        only:id,master
 use testutils, only:checkval,update_test_scores
 use dim,       only:do_radiation
 use part,      only:xyzmh_ptmass,vxyz_ptmass,nptmass
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(1)
 integer :: ierr,ieos,correct_answer
 character(len=20) :: pdir
 logical :: got_phantom_dir
 integer, parameter :: eos_to_test_for_u_from_Prho(4)=(/2,5,12,17/)

 nfailed = 0

 ! ieos=6 is for an isothermal disc around a sink particle, use isink=1
 isink = 1
 nptmass = 1
 xyzmh_ptmass(:,1) = 0.
 xyzmh_ptmass(4,1) = 1.
 vxyz_ptmass(:,1) = 0.

 ! ieos=8, barotropic eos, requires polyk to be set to avoid undefined
 polyk = 0.1
 polyk2 = 0.1
 qfacdisc = 0.5

 ! ieos=10 and 15, MESA and Helmholtz eos, require table read from files
 call get_environment_variable('PHANTOM_DIR',pdir)
 got_phantom_dir = (len_trim(pdir) > 0)

 do ieos=1,maxeos
    ! skip equations of state that are not implemented
    if (eos_is_not_implemented(ieos)) cycle

    if (id==master) write(*,"(/,a,i2)") '--> testing equation of state ',ieos
    call init_eos(ieos,ierr)
    correct_answer = 0
    if (ieos==10 .and. ierr /= 0 .and. .not. got_phantom_dir) cycle ! skip mesa
    if (ieos==15 .and. ierr /= 0 .and. .not. got_phantom_dir) cycle ! skip helmholtz
    if (ieos==16 .and. ierr /= 0 .and. .not. got_phantom_dir) cycle ! skip Shen
    if (ieos==24 .and. ierr /= 0 .and. .not. got_phantom_dir) cycle ! skip Stamatellos
    if (do_radiation .and. (ieos==10 .or. ieos==12 .or. ieos==20)) correct_answer = ierr_option_conflict
    call checkval(ierr,correct_answer,0,nfailed(1),'eos initialisation')
    call update_test_scores(ntests,nfailed,npass)
    !
    !--check u(P,rho) is inverse of P(u,rho) where implemented
    !
    if (any(eos_to_test_for_u_from_Prho==ieos)) then
       call test_u_from_Prho(ntests,npass,ieos)
    endif
    !
    !--check pressure is a continuous function
    !
    call test_p_is_continuous(ntests,npass,ieos)
 enddo

end subroutine test_all

!----------------------------------------------------------------------------
!+
!  test that the routine to solve for u from pressure and density works
!+
!----------------------------------------------------------------------------
subroutine test_u_from_Prho(ntests,npass,ieos)
 use io,               only:id,master,stdout
 use eos,              only:equationofstate,calc_temp_and_ene,gamma
 use testutils,        only:checkval,checkvalbuf_start,checkvalbuf,checkvalbuf_end,update_test_scores
 use units,            only:unit_density,unit_ergg
 use table_utils,      only:logspace
 integer, intent(inout) :: ntests,npass
 integer, intent(in)    :: ieos
 integer                :: npts,ierr,i,j,nfail(1),ncheck(1)
 real                   :: rhoi,eni,pri,ponrhoi,spsoundi,dum,tempi,tol,en_back
 real                   :: errmax(1)
 real, allocatable      :: rhogrid(:),ugrid(:)

 gamma = 5./3.
 npts = 50
 allocate(rhogrid(npts),ugrid(npts))
 if (ieos == 24) then
    call logspace(rhogrid,1e-24,1e0) ! cgs
    call logspace(ugrid,53020000.0,1.877E+14)  ! cgs
 else
    call logspace(rhogrid,1e-30,1e1) ! cgs
    call logspace(ugrid,1e-20,1e-4)  ! cgs
 endif

 dum = 0.
 tol = 1.e-15
 nfail = 0; ncheck = 0; errmax = 0.

 over_grid: do i=1,npts
    do j=1,npts
       ! get u from P, rho
       rhoi  = rhogrid(i)/unit_density
       eni   = ugrid(i)/unit_ergg
       tempi = -1. ! no initial guess
       call equationofstate(ieos,ponrhoi,spsoundi,rhoi,dum,dum,dum,tempi,eni)
       pri = ponrhoi * rhoi

       call calc_temp_and_ene(ieos,rhoi,pri,en_back,tempi,ierr) ! out = energy and temp
       if (ierr /= 0) exit over_grid

       call checkvalbuf(ugrid(i),real(en_back*unit_ergg),tol,'recovery of u(rho,P)',nfail(1),ncheck(1),errmax(1),use_rel_tol)
    enddo
 enddo over_grid

 if (ierr == 0) then
    call checkvalbuf_end('recovery of u(rho,P)',ncheck(1),nfail(1),errmax(1),tol)
    call update_test_scores(ntests,nfail,npass)
 else
    if (id==master) write(*,"(a)") ' skipped: not implemented for this eos'
 endif

end subroutine test_u_from_Prho

!----------------------------------------------------------------------------
!+
!  test ideal gas plus radiation eos: Check that P, T calculated from rho, u gives back
!  rho, u (assume fixed composition)
!+
!----------------------------------------------------------------------------
subroutine test_idealplusrad(ntests, npass)
 use io,               only:id,master,stdout
 use eos,              only:init_eos,equationofstate
 use eos_idealplusrad, only:get_idealplusrad_enfromtemp,get_idealplusrad_pres
 use testutils,        only:checkval,checkvalbuf_start,checkvalbuf,checkvalbuf_end,update_test_scores
 use units,            only:unit_density,unit_pressure,unit_ergg
 use physcon,          only:Rg
 integer, intent(inout) :: ntests,npass
 integer                :: npts,ieos,ierr,i,j,nfail(2),ncheck(2)
 real                   :: rhocodei,gamma,presi,dum,csound,eni,temp,ponrhoi,mu,tol,errmax(2),pres2,code_eni
 real, allocatable      :: rhogrid(:),Tgrid(:)

 if (id==master) write(*,"(/,a)") '--> testing ideal gas + radiation equation of state'

 ieos = 12
 mu = 0.6

 call get_rhoT_grid(npts,rhogrid,Tgrid)
 dum = 0.
 tol = 1.e-15
 nfail = 0; ncheck = 0; errmax = 0.
 call init_eos(ieos,ierr)
 do i=1,npts
    do j=1,npts
       ! Get u, P from rho, T
       call get_idealplusrad_enfromtemp(rhogrid(i),Tgrid(j),mu,eni)
       call get_idealplusrad_pres(rhogrid(i),Tgrid(j),mu,presi)

       ! Recalculate T, P, from rho, u
       code_eni = eni/unit_ergg
       temp = eni*mu/Rg ! guess
       rhocodei = rhogrid(i)/unit_density
       call equationofstate(ieos,ponrhoi,csound,rhocodei,dum,dum,dum,temp,code_eni,mu_local=mu,gamma_local=gamma)
       pres2 = ponrhoi * rhocodei * unit_pressure

       call checkvalbuf(temp,Tgrid(j),tol,'T from rho, u',nfail(1),ncheck(1),errmax(1),use_rel_tol)
       call checkvalbuf(pres2,presi,tol,'P from rho, u',nfail(2),ncheck(2),errmax(2),use_rel_tol)
    enddo
 enddo
 call checkvalbuf_end('T from rho, u',ncheck(1),nfail(1),errmax(1),tol)
 call checkvalbuf_end('P from rho, u',ncheck(2),nfail(2),errmax(2),tol)
 call update_test_scores(ntests,nfail,npass)

end subroutine test_idealplusrad

!----------------------------------------------------------------------------
!+
!  test HORMONE eos's: Check that P, T calculated from rho, u gives back
!  rho, u
!+
!----------------------------------------------------------------------------
subroutine test_hormone(ntests, npass)
 use io,        only:id,master,stdout
 use eos,       only:init_eos,equationofstate
 use eos_idealplusrad, only:get_idealplusrad_enfromtemp,get_idealplusrad_pres
 use eos_gasradrec, only:calc_uT_from_rhoP_gasradrec
 use ionization_mod, only:get_erec,get_imurec
 use testutils, only:checkval,checkvalbuf_start,checkvalbuf,checkvalbuf_end,update_test_scores
 use units,     only:unit_density,unit_pressure,unit_ergg
 integer, intent(inout) :: ntests,npass
 integer                :: npts,ieos,ierr,i,j,nfail(6),ncheck(6)
 real                   :: imurec,mu,eni_code,presi,pres2,dum,csound,eni,tempi,gamma_eff
 real                   :: ponrhoi,X,Z,tol,errmax(6),gasrad_eni,eni2,rhocodei,gamma,mu2
 real, allocatable      :: rhogrid(:),Tgrid(:)

 if (id==master) write(*,"(/,a)") '--> testing HORMONE equation of states'

 ieos = 20
 X = 0.69843
 Z = 0.01426

 call get_rhoT_grid(npts,rhogrid,Tgrid)

 ! Testing
 dum = 0.
 tol = 1.e-14
 tempi = -1.
 nfail = 0; ncheck = 0; errmax = 0.
 call init_eos(ieos,ierr)
 do i=1,npts
    do j=1,npts
       gamma = 5./3.
       ! Get mu, u, P from rho, T
       call get_imurec(log10(rhogrid(i)),Tgrid(j),X,1.-X-Z,imurec)
       mu = 1./imurec
       call get_idealplusrad_enfromtemp(rhogrid(i),Tgrid(j),mu,gasrad_eni)
       eni = gasrad_eni + get_erec(log10(rhogrid(i)),Tgrid(j),X,1.-X-Z)
       call get_idealplusrad_pres(rhogrid(i),Tgrid(j),mu,presi)

       ! Recalculate P, T from rho, u
       tempi = 1.
       eni_code = eni/unit_ergg
       rhocodei = rhogrid(i)/unit_density
       call equationofstate(ieos,ponrhoi,csound,rhocodei,0.,0.,0.,tempi,eni_code,&
                            mu_local=mu2,Xlocal=X,Zlocal=Z,gamma_local=gamma_eff)  ! mu and gamma_eff are outputs
       pres2 = ponrhoi * rhocodei * unit_pressure
       call checkvalbuf(mu2,mu,tol,'mu from rho, u',nfail(1),ncheck(1),errmax(1),use_rel_tol)
       call checkvalbuf(tempi,Tgrid(j),tol,'T from rho, u',nfail(2),ncheck(2),errmax(2),use_rel_tol)
       call checkvalbuf(pres2,presi,tol,'P from rho, u',nfail(3),ncheck(3),errmax(3),use_rel_tol)

       ! Recalculate u, T, mu from rho, P
       call calc_uT_from_rhoP_gasradrec(rhogrid(i),presi,X,1.-X-Z,tempi,eni2,mu2,ierr)
       call checkvalbuf(mu2,mu,tol,'mu from rho, P',nfail(4),ncheck(4),errmax(4),use_rel_tol)
       call checkvalbuf(tempi,Tgrid(j),tol,'T from rho, P',nfail(5),ncheck(5),errmax(5),use_rel_tol)
       call checkvalbuf(eni2,eni,tol,'u from rho, P',nfail(6),ncheck(6),errmax(6),use_rel_tol)
    enddo
 enddo
 call checkvalbuf_end('mu from rho, u',ncheck(1),nfail(1),errmax(1),tol)
 call checkvalbuf_end('T from rho, u',ncheck(2),nfail(2),errmax(2),tol)
 call checkvalbuf_end('P from rho, u',ncheck(3),nfail(3),errmax(3),tol)
 call checkvalbuf_end('mu from rho, P',ncheck(4),nfail(4),errmax(4),tol)
 call checkvalbuf_end('T from rho, P',ncheck(5),nfail(5),errmax(5),tol)
 call checkvalbuf_end('u from rho, P',ncheck(6),nfail(6),errmax(6),tol)
 call update_test_scores(ntests,nfail,npass)

end subroutine test_hormone

!----------------------------------------------------------------------------
!+
!  Helper routine to allocate density and temperature grids
!+
!----------------------------------------------------------------------------
subroutine get_rhoT_grid(npts,rhogrid,Tgrid)
 integer, intent(out) :: npts
 real, allocatable, intent(out) :: rhogrid(:),Tgrid(:)
 integer :: i
 real :: logQmin,logQmax,logTmin,logTmax
 real :: delta_logQ,delta_logT,logQi,logTi

 ! Initialise grids in Q and T (cgs units)
 npts = 30
 logQmin = -10.
 logQmax = -2.
 logTmin = 1.
 logTmax = 8.

 ! Note: logQ = logrho - 2logT + 12 in cgs units
 delta_logQ = (logQmax-logQmin)/real(npts-1)
 delta_logT = (logTmax-logTmin)/real(npts-1)

 allocate(rhogrid(npts),Tgrid(npts))
 do i=1,npts
    logQi = logQmin + (i-1)*delta_logQ
    logTi = logTmin + (i-1)*delta_logT
    rhogrid(i) = 10.**( logQi + 2.*logTi - 12. )
    Tgrid(i) = 10.**logTi
 enddo

end subroutine get_rhoT_grid

!----------------------------------------------------------------------------
!+
!  test piecewise barotropic eos has continuous pressure over density range
!+
!----------------------------------------------------------------------------
subroutine test_p_is_continuous(ntests, npass,ieos)
 use eos,            only:equationofstate,eos_requires_isothermal
 use eos_barotropic, only:rhocrit1cgs
 use eos_helmholtz,  only:eos_helmholtz_get_minrho
 use eos_tillotson,  only:rho_0,u_iv
 use testutils,      only:checkvalbuf,checkvalbuf_start,checkvalbuf_end,update_test_scores
 use units,          only:unit_density,unit_ergg!,unit_pressure,unit_velocity
 use mpiutils,       only:barrier_mpi
 integer, intent(inout) :: ntests,npass
 integer, intent(in)    :: ieos
 integer :: nfailed(2),ncheck(2)
 integer :: i,maxpts,ierrmax,itest
 real    :: rhoi,eni,xi,yi,zi,tempi,ponrhoi,spsoundi,ponrhoprev,spsoundprev
 real    :: errmax,rho_test
 character(len=3) :: var

 call barrier_mpi
 maxpts = 5001
 errmax = 0.
 select case(ieos)
 case(23)
    rhoi = 0.01*rho_0/unit_density
    eni = 0.01*u_iv/unit_ergg
    print*,' rho_0 = ',rho_0,' g/cm^3'
    rho_test = 0.5*rho_0/unit_density
 case(16,21,22,24)
    return
 case(15)
    rhoi = eos_helmholtz_get_minrho()
    eni = 1.e20/unit_ergg
    tempi = 1e4
    rho_test = 1000.*rhoi
 case default
    rhoi = 1.e-6*rhocrit1cgs/unit_density
    eni  = 1.e-2/unit_ergg
    tempi  = -1.
    rho_test = rhoi ! initial guess to avoid compiler warning
 end select
 xi = 3.
 yi = 2.
 zi = 1.

 over_tests: do itest=1,2
    nfailed = 0
    ncheck  = 0
    ! first test, fix u and vary rho
    var = 'rho'
    ! second test, fix rho and vary u
    if (itest==2) then
       var = 'u'
       rhoi = rho_test
       if (eos_requires_isothermal(ieos) .or. ieos==20) cycle over_tests
    endif
    !if (ieos==23 .and. itest==2) write(1,"(a)") '# rho,u,pressure,spsound'
    do i=1,maxpts
       if (itest==2) then
          eni = 1.002*eni
          tempi = 1.002*tempi
       else
          rhoi = 1.001*rhoi
       endif
       if (eos_requires_isothermal(ieos)) then
          call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi)
       else
          call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi,eni)
       endif
       !if (ieos==23 .and. itest==2) write(1,*) rhoi*unit_density,eni*unit_ergg,ponrhoi*rhoi*unit_pressure,spsoundi*unit_velocity
       if (i > 1) call checkvalbuf(ponrhoi,ponrhoprev,1.e-2,'p/rho continuous with '//trim(var),nfailed(1),ncheck(1),errmax)
       !if (i > 1) call checkvalbuf(spsoundi,spsoundprev,1.e-2,'cs is continuous',nfailed(2),ncheck(2),errmax)
       ponrhoprev = ponrhoi
       spsoundprev = spsoundi
    enddo
    ierrmax = 0
    call checkvalbuf_end('p/rho continuous with '//trim(var),ncheck(1),nfailed(1),ierrmax,0,maxpts-1)
 enddo over_tests

!call checkvalbuf_end('cs is continuous',ncheck(2),nfailed(2),0,0,maxpts)

 call update_test_scores(ntests,nfailed(1:1),npass)

end subroutine test_p_is_continuous


!----------------------------------------------------------------------------
!+
!  test helmholtz eos has continuous pressure
!+
!----------------------------------------------------------------------------
subroutine test_helmholtz(ntests, npass)
 use eos,           only:maxeos,equationofstate,eosinfo,init_eos
 use eos_helmholtz, only:eos_helmholtz_get_minrho, eos_helmholtz_get_maxrho, &
                         eos_helmholtz_get_mintemp, eos_helmholtz_get_maxtemp
 use io,            only:id,master,stdout
 use testutils,     only:checkval,checkvalbuf,checkvalbuf_start,checkvalbuf_end
 use units,         only:unit_density
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2),ncheck(2)
 integer :: i,j,ierr,maxpts,ierrmax,ieos
 real    :: rhoi,eni,tempi,xi,yi,zi,ponrhoi,spsoundi,ponrhoprev,spsoundprev
 real    :: errmax
 real    :: rhomin, rhomax, tempmin, tempmax, logdtemp, logdrho, logrhomin

 if (id==master) write(*,"(/,a)") '--> testing Helmholtz free energy equation of state'

 ieos = 15

 call init_eos(ieos, ierr)
 if (ierr /= 0) then
    write(*,"(/,a)") '--> skipping Helmholtz eos test due to init_eos() fail'
    return
 endif

 ntests  = ntests + 1
 nfailed = 0
 ncheck  = 0

 call eosinfo(ieos,stdout)
 call checkvalbuf_start('equation of state is continuous')

 rhomin  = eos_helmholtz_get_minrho()
 rhomax  = eos_helmholtz_get_maxrho()
 tempmin = eos_helmholtz_get_mintemp()
 tempmax = eos_helmholtz_get_maxtemp()

 maxpts = 7500
 errmax = 0.
 rhoi   = rhomin
 tempi  = tempmin

 logdtemp = log10(tempmax - tempmin) / (maxpts)
 logdrho  = (log10(rhomax) - log10(rhomin)) / (maxpts)
 logrhomin = log10(rhomin)

 ! run through temperature from 10^n with n=[5,13]
 do i=3,13
    tempi = 10.0**(i)
    ! run through density in log space
    do j=1,maxpts
       rhoi = 10**(logrhomin + j * logdrho)
       call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi,eni)
       write(1,*) rhoi*unit_density,ponrhoi,ponrhoi*rhoi,spsoundi
       if (j > 1) call checkvalbuf(ponrhoi, ponrhoprev, 5.e-2,'p/rho is continuous',nfailed(1),ncheck(1),errmax)
       !if (j > 1) call checkvalbuf(spsoundi,spsoundprev,1.e-2,'cs is continuous',   nfailed(2),ncheck(2),errmax)
       ponrhoprev  = ponrhoi
       spsoundprev = spsoundi
    enddo
 enddo

 ierrmax = 0
 call checkvalbuf_end('p/rho is continuous',ncheck(1),nfailed(1),ierrmax,0,maxpts)
 !call checkvalbuf_end('cs is continuous',   ncheck(2),nfailed(2),ierrmax,0,maxpts)

 if (nfailed(1)==0) npass = npass + 1

end subroutine test_helmholtz

end module testeos
