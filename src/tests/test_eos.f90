!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testeos
!
!  DESCRIPTION:
!  Unit tests of the equation of state module
!
!  REFERENCES: None
!
!  OWNER: Terrence Tricco
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, eos_helmholtz, io, mpiutils, physcon, testutils,
!    units
!+
!--------------------------------------------------------------------------
module testeos
 implicit none
 public :: test_eos

 private

contains

!----------------------------------------------------------
!+
!  run all unit tests of equation of state module
!+
!----------------------------------------------------------
subroutine test_eos(ntests,npass)
 use io,        only:id,master,stdout
 use physcon,   only:solarm
 use units,     only:set_units
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING EQUATION OF STATE MODULE'

 call set_units(mass=solarm,dist=1.d16,G=1.d0)

 call test_init(ntests, npass)
 call test_barotropic(ntests, npass)
 !call test_helmholtz(ntests, npass)

 if (id==master) write(*,"(/,a)") '<-- EQUATION OF STATE TEST COMPLETE'

end subroutine test_eos


!----------------------------------------------------------
!+
!  test that the initialisation of all eos works correctly
!+
!----------------------------------------------------------
subroutine test_init(ntests, npass)
 use eos,       only:maxeos,init_eos,isink,polyk,polyk2
 use io,        only:id,master
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(maxeos)
 integer :: ierr,ieos
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
    if (ieos==10 .and. .not. got_phantom_dir) cycle ! skip mesa
    if (ieos==15 .and. .not. got_phantom_dir) cycle ! skip helmholtz
    call checkval(ierr,0,0,nfailed(ieos),'eos initialisation')
 enddo
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_init


!----------------------------------------------------------------------------
!+
!  test piecewise barotropic eos has continuous pressure over density range
!+
!----------------------------------------------------------------------------
subroutine test_barotropic(ntests, npass)
 use eos,       only:equationofstate,rhocrit1cgs,polyk,polyk2,eosinfo,init_eos
 use io,        only:id,master,stdout
 use testutils, only:checkvalbuf,checkvalbuf_start,checkvalbuf_end,update_test_scores
 use units,     only:unit_density
 use mpiutils,  only:barrier_mpi
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2),ncheck(2)
 integer :: i,ierr,maxpts,ierrmax,ieos
 real    :: rhoi,xi,yi,zi,ponrhoi,spsoundi,ponrhoprev,spsoundprev
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

 do i=1,maxpts
    rhoi = 1.01*rhoi
    call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi)
    write(1,*) rhoi*unit_density,ponrhoi,ponrhoi*rhoi,spsoundi
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
                         eos_helmholtz_get_mintemp, eos_helmholtz_get_maxtemp, eos_helmholtz_set_relaxflag
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
       call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,eni,tempi)
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
