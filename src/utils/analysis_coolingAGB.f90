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
                              chemical_equilibrium_light
 use dim,              only:nElements
 use io,               only:id,master

 implicit none

 character(len=20), parameter, public :: analysistype = 'cooling'
 public :: do_analysis

 private
 integer :: analysis_to_perform
 real    :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867]
 real    :: eps(nElements) = [1.d0, 1.04d-1, 0.0,  6.d-4, 2.52d-4, 1.17d-4, 3.58d-5, 1.85d-5, 3.24d-5, 8.6d-8]

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
  use cooling_AGBwinds, only:nrates,dphot0,init_cooling_AGB,energ_cooling_AGB,dphotflag,nabn
  !use cooling,     only:energ_cooling
  ! use chem,           only:init_chem,get_dphot
  use dust_formation, only:chemical_equilibrium_light,eps,wind_CO_ratio,init_muGamma,set_abundances,mass_per_H
  use physcon,        only:Rg,mass_proton_cgs,kboltz,patm
  use units,          only:unit_ergg,unit_density,udist,utime
  integer, parameter :: nabundances = 16
  real :: abundance(nabundances)
  integer, parameter :: nt = 1000
  real :: logtmin,logtmax,logt,dlogt,t,crate
  real :: tempiso,ndens,xi,yi,zi,mu,rho_cgs,rhoi,ui,dt
  real :: h2ratio,dudti,dphot
  real(kind=4) :: divv_cgs
  integer :: i,ichem,iunit
  real    :: h2_H, o_H, oh_H, h20_H, co_H, cI_H, CII_H, siI_H, siII_H, e_H, &
             hp_H, hI_H, hd_H, heI_H, heII_H, heIII_H
  real    :: ratesq(nrates)
  real    :: pC, pC2, pC2H, pC2H2
  real    :: nH, nH2, nHe, nCO, nH2O, nOH, nO, nSi, nC2, nC, nC2H2, nSiO, nCH4
  integer, parameter :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10
  real    :: ndens_H, epsC
  real    :: gamma
  
  if (id==master) write(*,"(/,a)") '--> testing cooling_AGB rate'
  
  logtmax = log10(1.d5)
  logtmin = log10(5.d0)
 ! Abundances relative to H for C-rich outflows
  ! h2_H    = 0.38d0
  ! o_H     = 3.d-5   ! Atomic oxygen
  ! oh_H    = 1.d-7   ! Hydroxyl radical
  ! h20_H   = 2.55d-6 ! Water
  ! co_H    = 8.d-4   ! Carbon monoxide
  ! cI_H    = 2.d-5   ! Neutral carbon
  ! CII_H   = 1.d-6   ! Ionized carbon
  ! siI_H   = 1.d-7   ! Neutral silicon
  ! siII_H  = 5.d-6   ! Ionized silicon
  ! e_H     = 1.d-4   ! Free electrons (depends on ionization fraction)
  ! hp_H    = 1.d-8   ! Ionized hydrogen
  ! hI_H    = 0.d0    ! Neutral hydrogen
  ! hd_H    = 3.d-5   ! Hydrogen deuteride
  ! heI_H   = 2.4d-1  ! Neutral helium
  ! heII_H  = 1.d-9   ! Singly ionized helium
  ! heIII_H = 1.d-15  ! Doubly ionized helium (very rare in AGB outflows)

  !set atomic abundances

  call set_abundances
 
  ! abundance(1) = h2_H
  ! abundance(2) = o_H
  ! abundance(3) = oh_H
  ! abundance(4) = h20_H
  ! abundance(5) = co_H
  ! abundance(6) = cI_H
  ! abundance(7) = CII_H
  ! abundance(8) = siI_H
  ! abundance(9) = siII_H
  ! abundance(10) = e_H
  ! abundance(11) = hp_H
  ! abundance(12) = eps(1)
  ! abundance(13) = hd_H
  ! abundance(14) = eps(iHe)
  ! abundance(15) = heII_H
  ! abundance(16) = heIII_H

  epsC = eps(3)

  xi = 0.
  yi = 0.
  zi = 0.
  dt = 1.0d0
  rho_cgs = 2.0d-14
  rhoi = rho_cgs/unit_density 
  ndens_H = rhoi*unit_density / mass_per_H

  ! nH2  = ndens_H * h2_H
  ! nOH  = ndens_H * oh_H
  ! nH2O = ndens_H * h20_H
  ! nCO  = ndens_H * co_H                 
  ! nH   = ndens_H * hI_H
  ! nHe  = ndens_H * heI_H

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
    ui = 1.5*t*(Rg/mu)/unit_ergg  ! Comment by Davide: this formula is used in test_cooling, but I don't understand it
    
    ! dphot = get_dphot(dphotflag,dphot0,xi,yi,zi)

    call chemical_equilibrium_light(rho_cgs, t, epsC, pC, pC2, pC2H, pC2H2, mu, gamma, &
                                    nH, nH2, nHe, nCO, nH2O, nOH, nO, nSi, nC2, nC,    &
                                    nC2H2, nSiO, nCH4)
    
    abundance(1) = nH2 / ndens_H
    abundance(3) = nOH / ndens_H
    abundance(4) = nH2O / ndens_H
    abundance(5) = nCO / ndens_H
    abundance(6) = nC / ndens_H
    abundance(12) = nH / ndens_H
    abundance(14) = nHe / ndens_H
    abundance(2) = nO / ndens_H
    abundance(8) = nSi / ndens_H
    
                          
    h2_H    = abundance(1)
    oh_H    = abundance(3)
    h20_H   = abundance(4)
    co_H    = abundance(5)
    hI_H    = abundance(12)
    heI_H   = abundance(14)
    ! mu = 2.5d0
    call energ_cooling_AGB(ui,rhoi,divv_cgs,mu,abundance,dudti,ratesq)

    ndens = rhoi*unit_density/(mu*mass_proton_cgs)
    crate = dudti*udist**2/utime**3*(rhoi*unit_density)
    write(iunit,*)  t,crate/ndens**2,                     &
                    h2_H, oh_H, h20_H, co_H, hI_H, heI_H, &
                    abundance(2), nC2/ndens_H, nC/ndens_H, &
                    nC2H2/ndens_H, nSi/ndens_H, nSiO/ndens_H, &
                    nCH4/ndens_H, &
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
  use cooling_AGBwinds, only:nrates,dphot0,init_cooling_AGB,energ_cooling_AGB,dphotflag,nabn
  use dust_formation, only:chemical_equilibrium_light,eps,wind_CO_ratio,init_muGamma,set_abundances,mass_per_H
  use physcon,        only:Rg,mass_proton_cgs,kboltz,patm
  use units,          only:unit_ergg,unit_density,udist,utime
  integer, parameter :: nabundances = 16
  real :: abundance(nabundances)
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
  real    :: pC, pC2, pC2H, pC2H2
  real    :: nH, nH2, nHe, nCO, nH2O, nOH, nO, nSi, nC2, nC, nC2H2, nSiO, nCH4
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
      ui = 1.5*t*(Rg/mu)/unit_ergg  ! Comment by Davide: this formula is used in test_cooling, but I don't understand it
      

      ! print '(A19, F8.0, A7, ES12.2)', "Next temp will be :", t, "  yn: ", yn

      call chemical_equilibrium_light(rho_cgs, t, epsC, pC, pC2, pC2H, pC2H2, mu, gamma, &
      nH, nH2, nHe, nCO, nH2O, nOH, nO, nSi, nC2, nC, nC2H2, nSiO, nCH4)

      abundance(1) = nH2 / ndens_H
      abundance(3) = nOH / ndens_H
      abundance(4) = nH2O / ndens_H
      abundance(5) = nCO / ndens_H
      abundance(6) = nC / ndens_H
      abundance(12) = nH / ndens_H
      abundance(14) = nHe / ndens_H
      abundance(2) = nO / ndens_H
      abundance(8) = nSi / ndens_H


      h2_H    = abundance(1)
      oh_H    = abundance(3)
      h20_H   = abundance(4)
      co_H    = abundance(5)
      hI_H    = abundance(12)
      heI_H   = abundance(14)
      call energ_cooling_AGB(ui,rhoi,divv_cgs,mu,abundance,dudti,ratesq)

      ndens = (rhoi*unit_density/mass_proton_cgs)*5.d0/7.d0
      crate = dudti*udist**2/utime**3*(rhoi*unit_density)
      write(iunit,*)  t,ndens_H,crate/ndens**2,                     &
                      h2_H, oh_H, h20_H, co_H, hI_H, heI_H, &
                      abundance(2), nC2/ndens_H, nC/ndens_H,&
                      nC2H2/ndens_H, nSi/ndens_H, nSiO/ndens_H, &
                      nCH4/ndens_H, &
                      ratesq(:)/ndens**2
    enddo
  enddo
  close(iunit)

end subroutine cooling_rate_temp_dens

end module analysis