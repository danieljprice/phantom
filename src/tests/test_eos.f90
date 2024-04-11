!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
!   eos_idealplusrad, io, ionization_mod, mpiutils, physcon,
!   testeos_stratified, testutils, units
!
 implicit none
 public :: test_eos
 public :: test_helmholtz ! to avoid compiler warning for unused routine

 private

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

 call test_init(ntests, npass)
 call test_barotropic(ntests, npass)
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
subroutine test_init(ntests, npass)
 use eos,       only:maxeos,init_eos,isink,polyk,polyk2,&
                     ierr_file_not_found,ierr_option_conflict
 use io,        only:id,master
 use testutils, only:checkval,update_test_scores
 use dim,       only:do_radiation
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(maxeos)
 integer :: ierr,ieos,correct_answer
 character(len=20) :: pdir
 logical :: got_phantom_dir

 if (id==master) write(*,"(/,a)") '--> testing equation of state initialisation'

 nfailed = 0

 ! ieos=6 is for an isothermal disc around a sink particle, use isink=1
 isink = 1

 ! ieos=8, barotropic eos, requires polyk to be set to avoid undefined
 polyk = 0.1
 polyk2 = 0.1

 ! ieos=10 and 15, MESA and Helmholtz eos, require table read from files
 call get_environment_variable('PHANTOM_DIR',pdir)
 got_phantom_dir = (len_trim(pdir) > 0)

 do ieos=1,maxeos
    call init_eos(ieos,ierr)
    correct_answer = 0
    if (ieos==10 .and. ierr /= 0 .and. .not. got_phantom_dir) cycle ! skip mesa
    if (ieos==15 .and. ierr /= 0 .and. .not. got_phantom_dir) cycle ! skip helmholtz
    if (ieos==16 .and. ierr /= 0 .and. .not. got_phantom_dir) cycle ! skip Shen
    if (do_radiation .and. (ieos==10 .or. ieos==12 .or. ieos==20)) correct_answer = ierr_option_conflict
    call checkval(ierr,correct_answer,0,nfailed(ieos),'eos initialisation')
 enddo
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_init


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
 use physcon,          only:kb_on_mh
 integer, intent(inout) :: ntests,npass
 integer                :: npts,ieos,ierr,i,j,nfail(2),ncheck(2)
 real                   :: rhocodei,gamma,presi,dum,csound,eni,temp,ponrhoi,mu,tol,errmax(2),pres2,code_eni
 real, allocatable      :: rhogrid(:),Tgrid(:)

 if (id==master) write(*,"(/,a)") '--> testing ideal gas + radiation equation of state'

 ieos = 12
 mu = 0.6

 call get_rhoT_grid(npts,rhogrid,Tgrid)
 dum = 0.
 tol = 1.e-12
 nfail = 0; ncheck = 0; errmax = 0.
 call init_eos(ieos,ierr)
 do i=1,npts
    do j=1,npts
       ! Get u, P from rho, T
       call get_idealplusrad_enfromtemp(rhogrid(i),Tgrid(j),mu,eni)
       call get_idealplusrad_pres(rhogrid(i),Tgrid(j),mu,presi)

       ! Recalculate T, P, from rho, u
       code_eni = eni/unit_ergg
       temp = eni*mu/kb_on_mh
       rhocodei = rhogrid(i)/unit_density
       call equationofstate(ieos,ponrhoi,csound,rhocodei,dum,dum,dum,temp,code_eni,mu_local=mu,gamma_local=gamma)
       pres2 = ponrhoi * rhocodei * unit_pressure

       call checkvalbuf(temp,Tgrid(j),tol,'Check recovery of T from rho, u',nfail(1),ncheck(1),errmax(1))
       call checkvalbuf(pres2,presi,tol,'Check recovery of P from rho, u',nfail(2),ncheck(2),errmax(2))
    enddo
 enddo
 call checkvalbuf_end('Check recovery of T from rho, u',ncheck(1),nfail(1),errmax(1),tol)
 call checkvalbuf_end('Check recovery of P from rho, u',ncheck(2),nfail(2),errmax(2),tol)
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
 integer                :: npts,ieos,ierr,i,j,nfail(4),ncheck(4)
 real                   :: imurec,mu,eni_code,presi,pres2,dum,csound,eni,tempi
 real                   :: ponrhoi,X,Z,tol,errmax(4),gasrad_eni,eni2,rhocodei,gamma
 real, allocatable      :: rhogrid(:),Tgrid(:)

 if (id==master) write(*,"(/,a)") '--> testing HORMONE equation of states'

 ieos = 20
 X = 0.69843
 Z = 0.01426

 call get_rhoT_grid(npts,rhogrid,Tgrid)

 ! Testing
 dum = 0.
 tol = 1.e-12
 tempi = -1.
 nfail = 0; ncheck = 0; errmax = 0.
 call init_eos(ieos,ierr)
 tempi = 1.
 eni_code =  764437650.64783347/unit_ergg
 rhocodei = 3.2276168501594796E-015/unit_density
 call equationofstate(ieos,ponrhoi,csound,rhocodei,0.,0.,0.,tempi,eni_code,mu_local=mu,Xlocal=X,Zlocal=Z,gamma_local=gamma)
 do i=1,npts
    do j=1,npts
       gamma = 5./3.
       ! Get mu from rho, T
       call get_imurec(log10(rhogrid(i)),Tgrid(j),X,1.-X-Z,imurec)
       mu = 1./imurec

       ! Get u, P from rho, T, mu
       call get_idealplusrad_enfromtemp(rhogrid(i),Tgrid(j),mu,gasrad_eni)
       eni = gasrad_eni + get_erec(log10(rhogrid(i)),Tgrid(j),X,1.-X-Z)
       call get_idealplusrad_pres(rhogrid(i),Tgrid(j),mu,presi)

       ! Recalculate P, T from rho, u, mu
       tempi = 1.
       eni_code = eni/unit_ergg
       rhocodei = rhogrid(i)/unit_density
       call equationofstate(ieos,ponrhoi,csound,rhocodei,0.,0.,0.,tempi,eni_code,mu_local=mu,Xlocal=X,Zlocal=Z,gamma_local=gamma)
       pres2 = ponrhoi * rhocodei * unit_pressure
       call checkvalbuf(tempi,Tgrid(j),tol,'Check recovery of T from rho, u',nfail(1),ncheck(1),errmax(1))
       call checkvalbuf(pres2,presi,tol,'Check recovery of P from rho, u',nfail(2),ncheck(2),errmax(2))

       ! Recalculate u, T, mu from rho, P
       call calc_uT_from_rhoP_gasradrec(rhogrid(i),presi,X,1.-X-Z,tempi,eni2,mu,ierr)
       call checkvalbuf(tempi,Tgrid(j),tol,'Check recovery of T from rho, P',nfail(3),ncheck(3),errmax(3))
       call checkvalbuf(eni2,eni,tol,'Check recovery of u from rho, P',nfail(4),ncheck(4),errmax(4))
    enddo
 enddo
 call checkvalbuf_end('Check recovery of T from rho, u',ncheck(1),nfail(1),errmax(1),tol)
 call checkvalbuf_end('Check recovery of P from rho, u',ncheck(2),nfail(2),errmax(2),tol)
 call checkvalbuf_end('Check recovery of T from rho, P',ncheck(3),nfail(3),errmax(3),tol)
 call checkvalbuf_end('Check recovery of u from rho, P',ncheck(4),nfail(4),errmax(4),tol)
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
subroutine test_barotropic(ntests, npass)
 use eos,            only:equationofstate,polyk,polyk2,eosinfo,init_eos
 use eos_barotropic, only:rhocrit1cgs
 use io,             only:id,master,stdout
 use testutils,      only:checkvalbuf,checkvalbuf_start,checkvalbuf_end,update_test_scores
 use units,          only:unit_density
 use mpiutils,       only:barrier_mpi
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2),ncheck(2)
 integer :: i,ierr,maxpts,ierrmax,ieos
 real    :: rhoi,xi,yi,zi,tempi,ponrhoi,spsoundi,ponrhoprev,spsoundprev
 real    :: errmax

 if (id==master) write(*,"(/,a)") '--> testing barotropic equation of state'

 ieos = 8

 ! polyk has to be specified here before we call init_eos()
 polyk  = 0.1
 polyk2 = 0.1

 call init_eos(ieos, ierr)
 if (ierr /= 0) then
    if (id==master) write(*,"(/,a)") '--> skipping barotropic eos test due to init_eos() fail'
    return
 endif

 nfailed = 0
 ncheck  = 0

 if (id==master) call eosinfo(ieos,stdout)
 call barrier_mpi
 call checkvalbuf_start('equation of state is continuous')

 maxpts = 5000
 errmax = 0.
 rhoi   = 1.e-6*rhocrit1cgs/unit_density
 tempi  = -1. ! initial guess to avoid compiler warning

 do i=1,maxpts
    rhoi = 1.01*rhoi
    call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi)
    !write(1,*) rhoi*unit_density,ponrhoi,ponrhoi*rhoi,spsoundi
    if (i > 1) call checkvalbuf(ponrhoi,ponrhoprev,1.e-2,'p/rho is continuous',nfailed(1),ncheck(1),errmax)
    !if (i > 1) call checkvalbuf(spsoundi,spsoundprev,1.e-2,'cs is continuous',nfailed(2),ncheck(2),errmax)
    ponrhoprev = ponrhoi
    spsoundprev = spsoundi
 enddo

 ierrmax = 0
 call checkvalbuf_end('p/rho is continuous',ncheck(1),nfailed(1),ierrmax,0,maxpts)
 !call checkvalbuf_end('cs is continuous',ncheck(2),nfailed(2),0,0,maxpts)

 call update_test_scores(ntests,nfailed(1:1),npass)

end subroutine test_barotropic


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
