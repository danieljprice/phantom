!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module testcooling
!
! Unit tests of the cooling module
!
! :References:
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: chem, cooling, h2cooling, io, part, physcon, testutils,
!   units
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_cooling, test_cooling_rate

 private

contains

!--------------------------------------------
!+
!  Various tests of the cooling module
!+
!--------------------------------------------
subroutine test_cooling(ntests,npass)
 use physcon, only:solarm,kpc
 use units,   only:set_units
 integer, intent(inout) :: ntests,npass
 !integer :: nfailed(10),ierr,iregime

 if (id==master) write(*,"(/,a)") '--> TESTING COOLING MODULE'

 call set_units(mass=1e7*solarm,dist=kpc,G=1.d0)

 !call test_cooling_rate(ntests,npass)

 call test_coolfunc(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- COOLING TEST COMPLETE'

end subroutine test_cooling

!--------------------------------------------
!+
!  Check that cooling function is continuous
!+
!--------------------------------------------
subroutine test_cooling_rate(ntests,npass)
 use h2cooling, only:nrates,dphot0,init_h2cooling,energ_h2cooling,dphotflag,abundsi,abundo,abunde,abundc,nabn
 use chem,      only:update_abundances,init_chem,get_dphot
 use part,      only:nabundances,iHI
 use physcon,   only:Rg,mass_proton_cgs
 use units,     only:unit_ergg,unit_density,udist,utime
 real :: abundance(nabundances)
 !real :: ratesq(nrates)
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nt = 400
 real :: logtmin,logtmax,logt,dlogt,t,crate
 real :: tempiso,ndens,xi,yi,zi,gmwvar,rhoi,ui,dt
 real :: h2ratio,dudti,dphot
 real(kind=4) :: divv_cgs
 integer :: i,ichem,iunit
 real    :: abundi(nabn)

 if (id==master) write(*,"(/,a)") '--> testing h2cooling rate'

 logtmax = log10(1.e8)
 logtmin = log10(1.e0)
 abundance(:)   = 0.
 abundance(iHI) = 1.  ! assume all atomic hydrogen initially
 ichem = 0 ! don't do chemistry
 xi = 0.
 yi = 0.
 zi = 0.
 dt = 0.1
 rhoi = 2.3e-24/unit_density
 h2ratio = 0.
 gmwvar=1.4/1.1
 ndens = rhoi*unit_density/(gmwvar*mass_proton_cgs)
 print*,' rho = ',rhoi, ' ndens = ',ndens
 call init_chem()
 call init_h2cooling()

 open(newunit=iunit,file='cooltable.txt',status='replace')
 write(iunit,"(a)") '#   T   \Lambda_E(T) erg s^{-1} cm^3   \Lambda erg s^{-1} cm^{-3}'
 dlogt = (logtmax - logtmin)/real(nt)
 divv_cgs = 0.
 do i=1,nt
    dudti = 0.
    logt = logtmin + (i-1)*dlogt
    t = 10**logt
    ui = 1.5*t*(Rg/gmwvar)/unit_ergg
    dphot = get_dphot(dphotflag,dphot0,xi,yi,zi)
    call update_abundances(ui,rhoi,abundance,nabundances,dphot,dt,abundi,nabn,gmwvar,abundc,abunde,abundo,abundsi)
    call energ_h2cooling(ui,rhoi,divv_cgs,gmwvar,abundi,dudti)
!call cool_func(tempiso,ndens,dlq,divv_cgs,abund,crate,ratesq)
    ndens = (rhoi*unit_density/mass_proton_cgs)*5.d0/7.d0
    crate = dudti*udist**2/utime**3*(rhoi*unit_density)
    if (abs(logt-7.) < dlogt) then
       print "(a,9(es10.3,1x))",' abundances= ',abundance(:)
       print*,'T = ',t,' Tiso = ',tempiso,' n = ',ndens,' cm^-3',&
       ' Lam = ',crate/ndens**2
    endif
    write(iunit,*) t,crate/ndens**2,crate
 enddo
 close(iunit)

end subroutine test_cooling_rate

!--------------------------------------------
!+
!  Check that cooling function is continuous
!+
!--------------------------------------------
subroutine test_coolfunc(ntests,npass)
 use cooling,  only:find_in_table
 use testutils, only:checkvalbuf,checkvalbuf_start
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nt = 100
 integer :: i,k,ndiff(1),ncheck
 real    :: table(nt),val
 logical :: my_test

 if (id==master) write(*,"(/,a)") '--> testing find_in_table routine'
 !
 ! set up table
 !
 do i=1,nt
    table(i) = i
 enddo

 ndiff = 0
 ncheck = 0
 call checkvalbuf_start('table(i) < val < table(i+1)')
 do i=1,nt-1
    val = i+0.5
    k = find_in_table(nt,table,val)
    !print*,k,table(k),val,table(k+1)
    if (k < nt) then
       my_test = (val > table(k) .and. val < table(k+1))
       call checkvalbuf(my_test,.true.,'table(i) < val < table(i+1)',ndiff(1),ncheck)
    endif
 enddo
 call update_test_scores(ntests,ndiff,npass)

end subroutine test_coolfunc

end module testcooling
