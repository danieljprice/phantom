!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! various tests of the cooling solver module
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: cooling, cooling_functions, cooling_solver, dim,
!   dust_formation, options, physcon, prompting, units
!

 use cooling
 use cooling_functions
 use cooling_solver
 use physcon,          only:mass_proton_cgs,kboltz,atomic_mass_unit,patm
 use dust_formation,   only:init_muGamma,set_abundances,kappa_gas,calc_kappa_bowen,&
                              chemical_equilibrium_light,init_nucleation, eps
 use dim,              only:nElements, nabn_AGB
 use io,               only:id,master

 implicit none

 character(len=20), parameter, public :: analysistype = 'cooling'
 public :: do_analysis

 private
 integer :: analysis_to_perform
 real    :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867]
! Indices for cooling species:
 integer, parameter :: icoolH=1, icoolC=2, icoolO=3, icoolSi=4, icoolH2=5, icoolCO=6, &
                       icoolH2O=7, icoolOH=8, icoolC2=9, icoolC2H=10, icoolC2H2=11, &
                       icoolHe=12, icoolSiO=13, icoolCH4=14

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 use prompting,  only:prompt

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real(kind=8),     intent(in) :: xyzh(:,:),vxyzu(:,:)
 real(kind=8),     intent(in) :: particlemass,time


 print "(29(a,/))", &
 ' 1) check only for temperatute', &
 ' 2) check for temperature and density (to be checked)'

analysis_to_perform = 1

call prompt('Choose analysis type ',analysis_to_perform,1,2)
print *,''

!analysis
select case(analysis_to_perform)
case(1) !test rate
call test_cooling()
case(2)
call cooling_temp_dens()
!  case(3)
!     call test_cooling_solvers
end select

end subroutine do_analysis

!--------------------------------------------
!+
!  Various tests of the cooling module
!+
!--------------------------------------------
subroutine test_cooling()
  use physcon, only:solarm,kpc
  use units,   only:set_units
  !integer :: nfailed(10),ierr,iregime
  
  if (id==master) write(*,"(/,a)") '--> TESTING COOLING MODULE'
  
  call set_units(mass=1e7*solarm,dist=kpc,G=1.d0)
  
  call test_cooling_rate()
  
  if (id==master) write(*,"(/,a)") '<-- COOLING TEST COMPLETE'
  
  end subroutine test_cooling

!--------------------------------------------
!+
!  Cooling rates on temperature grid
!+
!--------------------------------------------
subroutine test_cooling_rate()
  use cooling_AGBwinds, only:nrates,dphot0,init_cooling_AGB,energ_cooling_AGB,dphotflag
  !use cooling,     only:energ_cooling
  ! use chem,           only:init_chem,get_dphot
  use dust_formation, only:chemical_equilibrium_light,wind_CO_ratio,init_muGamma,mass_per_H
  use physcon,        only:Rg,mass_proton_cgs,kboltz,patm
  use units,          only:unit_ergg,unit_density,udist,utime
  use dim,            only:nElements

implicit none

  integer, parameter :: nt = 1000
  real :: logtmin,logtmax,logt,dlogt,t,crate
  real :: tempiso,ndens,xi,yi,zi,mu,rho_cgs,rhoi,ui,dt
  real :: h2ratio,dudti,dphot
  real(kind=4) :: divv_cgs
  integer :: i,ichem,iunit
  real    :: h2_H, o_H, oh_H, h20_H, co_H, cI_H, CII_H, siI_H, siII_H, e_H, &
             hp_H, hI_H, hd_H, heI_H, heII_H, heIII_H
  real    :: ratesq(nrates)
  real    :: abundi(nabn_AGB)
  integer, parameter :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10
  real    :: ndens_H, epsC
  real    :: gamma
  
  if (id==master) write(*,"(/,a)") '--> testing cooling_AGB rate'
  
  logtmax = log10(2.2d4)
  logtmin = log10(2.d0)
 

  call set_abundances
 

  epsC = eps(3)

  xi = 0.
  yi = 0.
  zi = 0.
  dt = 1.0d0
  rho_cgs = 2.0d-14
  rhoi = rho_cgs/unit_density 
  ndens_H = rhoi*unit_density / mass_per_H


  call init_cooling_AGB()
  
  ! icooling = 1      ! use cooling solver
  ! excitation_HI = 1 ! H1 cooling
  ! icool_method = 1  ! explicit
  
  open(newunit=iunit,file='cooltable.txt',status='replace')
  write(iunit,'(A, E12.4)') '#   T   \Lambda_E(T) erg s^{-1} cm^3   N dens H: ', ndens_H
  dlogt = (logtmax - logtmin)/real(nt)
  divv_cgs = 0.
  ! T = 1.5d3
  ! ui = 1.5*T*(Rg/mu)/unit_ergg
  ! call energ_cooling_AGB(ui,rhoi,divv_cgs,mu,abundance,dudti,ratesq)
  ! ndens = (rhoi*unit_density/mass_proton_cgs)*5.d0/7.d0
  ! crate = dudti*udist**2/utime**3*(rhoi*unit_density)
  ! write(iunit,*) t,crate/ndens**2,ratesq(:)/ndens**2

  do i=1,nt
    dudti = 0.
    logt = logtmin + (i-1)*dlogt
    t = 10**logt
    call init_muGamma(rho_cgs, t, mu, gamma)
    
    ! dphot = get_dphot(dphotflag,dphot0,xi,yi,zi)

    call chemical_equilibrium_light(rho_cgs, t,  mu, gamma, abundi)
    
    abundi = abundi / ndens_H
    
    call energ_cooling_AGB(t,rhoi,divv_cgs,mu,abundi,dudti,ratesq)

    ndens = rhoi*unit_density/(mu*mass_proton_cgs)
    crate = dudti*(rhoi*unit_density)
    write(iunit,*)  t,crate/ndens**2,                     &
                    abundi(icoolH2), abundi(icoolOH), abundi(icoolH2O), &
                    abundi(icoolCO), abundi(icoolH), abundi(icoolHe), &
                    abundi(icoolO), abundi(icoolC2), abundi(icoolC), &
                    abundi(icoolC2H2), abundi(icoolSi), abundi(icoolSiO), &
                    abundi(icoolCH4), &
                    ratesq(:)/ndens**2
  enddo
  close(iunit)
  
  end subroutine test_cooling_rate


!--------------------------------------------
!+
!  Various tests of the cooling module for temperature and density
!+
!--------------------------------------------
subroutine cooling_temp_dens()
  use physcon, only:solarm,kpc
  use units,   only:set_units
  !integer :: nfailed(10),ierr,iregime

  if (id==master) write(*,"(/,a)") '--> TESTING COOLING MODULE'

  call set_units(mass=1e7*solarm,dist=kpc,G=1.d0)

  call cooling_rate_temp_dens()

  if (id==master) write(*,"(/,a)") '<-- COOLING TEST COMPLETE'

end subroutine cooling_temp_dens

!------------------------------------------------
!+
!  Test cooling rates for temperature and density
!+
!------------------------------------------------
subroutine cooling_rate_temp_dens()
  use cooling_AGBwinds, only:nrates,dphot0,init_cooling_AGB,energ_cooling_AGB,dphotflag
  use dust_formation, only:chemical_equilibrium_light,eps,wind_CO_ratio,init_muGamma,set_abundances,mass_per_H
  use physcon,        only:Rg,mass_proton_cgs,kboltz,patm
  use units,          only:unit_ergg,unit_density,udist,utime
  integer, parameter :: nt_grid = 15
  integer, parameter :: nd_grid = 15
  real :: logtmin,logtmax,logrhomin,logrhomax,logt,logrho,dlogt,dlogrho,t,crate
  real :: tempiso,ndens,xi,yi,zi,mu,rho_cgs,rhoi,ui,dt
  real :: h2ratio,dudti,dphot
  real(kind=4) :: divv_cgs
  integer :: i,l,iunit
  real    :: h2_H, o_H, oh_H, h20_H, co_H, cI_H, CII_H, siI_H, siII_H, e_H, &
             hp_H, hI_H, hd_H, heI_H, heII_H, heIII_H
  real    :: ratesq(nrates)
  real    :: abundi(nabn_AGB)
  integer, parameter :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10
  real    :: ndens_H, epsC
  real    :: gamma
  real    :: pC_old, pO_old, pSi_old, pS_old, pTi_old


  if (id==master) write(*,"(/,a)") '--> testing cooling_AGB rate'

  logtmax = log10(2.2d4)
  logtmin = log10(5.d2)
  logrhomin = 4.0d0
  logrhomax = 1.45d1

  !set atomic abundances

  call set_abundances

  epsC = eps(3)

  xi = 0.
  yi = 0.
  zi = 0.
  dt = 1.0d0
  ! rho_cgs = 2.0d-12
  ! rhoi = rho_cgs/unit_density 
  ! ndens_H = rhoi*unit_density / mass_per_H

  call init_cooling_AGB()


  open(newunit=iunit,file='cooltable_dens.txt',status='replace')
  write(iunit,"(a)") '#   T   \Lambda_E(T) erg s^{-1} cm^3   \Lambda erg s^{-1} cm^{-3}'

  dlogt = (logtmax - logtmin)/real(nt_grid)
  dlogrho = (logrhomax - logrhomin)/real(nd_grid)
  divv_cgs = 0.

  do l=1,nd_grid
    logrho = logrhomin + (l-1)*dlogrho
    ndens_H = 10**logrho
    rho_cgs = ndens_H * mass_per_H
    rhoi = rho_cgs/unit_density 
    do i=1,nt_grid
      logt = logtmin + (i-1)*dlogt
      t = 10**logt
      call init_muGamma(rho_cgs, t, mu, gamma)
      

      call chemical_equilibrium_light(rho_cgs, t,  mu, gamma, abundi)

      abundi = abundi / ndens_H

      call energ_cooling_AGB(t,rhoi,divv_cgs,mu,abundi,dudti,ratesq)

      ndens = (rhoi*unit_density/mass_proton_cgs)*5.d0/7.d0
      crate = dudti*(rhoi*unit_density)
      write(iunit,*)  t,ndens_H,crate/ndens**2,                     &
                      abundi(icoolH2), abundi(icoolOH), abundi(icoolH2O), &
                      abundi(icoolCO), abundi(icoolH), abundi(icoolH), &
                      abundi(icoolO), abundi(icoolC2), abundi(icoolC), &
                      abundi(icoolC2H2), abundi(icoolSi), abundi(icoolSiO), &
                      abundi(icoolCH4), &
                      ratesq(:)/ndens**2
    enddo
  enddo
  close(iunit)

end subroutine cooling_rate_temp_dens

end module analysis
