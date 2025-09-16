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
 use part,             only:isdead_or_accreted

 implicit none

 character(len=20), parameter, public :: analysistype = 'cooling'
 public :: do_analysis

 private
 integer :: analysis_to_perform
 real    :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867]
! Indices for cooling species:
 integer, parameter :: icoolH=1, icoolC=2, icoolO=3, icoolSi=4, icoolH2=5, icoolCO=6, &
                       icoolH2O=7, icoolOH=8, icoolC2=9, icoolC2H=10, icoolC2H2=11, &
                       icoolHe=12, icoolSiO=13, icoolCH4=14, icoolS=15, icoolTi=16

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 use prompting,  only:prompt

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real(kind=8),     intent(in) :: xyzh(:,:),vxyzu(:,:)
 real(kind=8),     intent(in) :: particlemass,time


 print "(29(a,/))", &
 ' 1) check cooling only for temperature', &
 ' 2) check cooling for temperature and density (to be checked)', &
 ' 3) check dust formation and destruction', &
 ' 4) total dust mass'

analysis_to_perform = 1

call prompt('Choose analysis type ',analysis_to_perform,1,4)
print *,''

!analysis
select case(analysis_to_perform)
case(1) !test rate
  call test_cooling()
case(2)
  call cooling_temp_dens()
case(3)
  call compute_dust_formation()
case(4)
  print*, 'particle mass: ', particlemass
  call total_dust_mass(time,npart,particlemass,xyzh)
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
  use dust_formation, only:chemical_equilibrium_light,wind_CO_ratio,init_muGamma,mass_per_H, &
                            chemical_equilibrium_light_fixed_mu_gamma_broyden, &
                            chemical_equilibrium_light_fixed_mu_gamma
  use physcon,        only:Rg,mass_proton_cgs,kboltz,patm
  use units,          only:unit_ergg,unit_density,udist,utime
  use dim,            only:nElements

implicit none

  integer, parameter :: nt = 1000
  real :: logtmin,logtmax,logt,dlogt,t,crate,Tdust
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
  real    :: start, finish
  
  if (id==master) write(*,"(/,a)") '--> testing cooling_AGB rate'
  
  logtmax = log10(1.5d4)
  logtmin = log10(5.d0)
 

  call set_abundances
 

  epsC = eps(3) ! ignoring nucleation

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

  call cpu_time(start)

  do i=1,nt
    dudti = 0.
    logt = logtmin + (i-1)*dlogt
    t = 10**logt
    call init_muGamma(rho_cgs, t, mu, gamma)
    
    ! dphot = get_dphot(dphotflag,dphot0,xi,yi,zi)

    ! call chemical_equilibrium_light_fixed_mu_gamma(rho_cgs, t, epsC, mu, gamma, abundi)
    ! call chemical_equilibrium_light_fixed_mu_gamma_broyden(rho_cgs, t, epsC, mu, gamma, abundi)
    call chemical_equilibrium_light(rho_cgs, t, epsC, mu, gamma, abundi)
    
    abundi = abundi / ndens_H
    Tdust = t
    
    call energ_cooling_AGB(t,Tdust,rhoi,divv_cgs,mu,abundi,dudti,ratesq)

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
  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
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
  use dust_formation, only:chemical_equilibrium_light_fixed_mu_gamma_broyden,eps, &
                           wind_CO_ratio,init_muGamma,set_abundances,mass_per_H
  use physcon,        only:Rg,mass_proton_cgs,kboltz,patm
  use units,          only:unit_ergg,unit_density,udist,utime
  integer, parameter :: nt_grid = 15
  integer, parameter :: nd_grid = 15
  real :: logtmin,logtmax,logrhomin,logrhomax,logt,logrho,dlogt,dlogrho,t,crate,Tdust
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

  logtmax = log10(2.d3)
  logtmin = log10(8.d2)
  logrhomin = 4.0d0
  logrhomax = 1.45d1

  !set atomic abundances

  call set_abundances

  epsC = eps(3)  ! ignoring nucleation

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
      

      ! call chemical_equilibrium_light_fixed_mu_gamma_broyden(rho_cgs, t, epsC, mu, gamma, abundi)
      call chemical_equilibrium_light(rho_cgs, t, epsC, mu, gamma, abundi)

      abundi = abundi / ndens_H
      Tdust = t

      call energ_cooling_AGB(t,Tdust,rhoi,divv_cgs,mu,abundi,dudti,ratesq)

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

! Subroutines to compute dust formation and destruction

subroutine compute_dust_formation()
  use dust_formation, only: set_abundances, evolve_chem,   &
                            wind_CO_ratio, eps, mass_per_H
  use physcon,        only: pi
  use part,           only: n_nucleation, idJstar, idK0, idK1, idK2, idK3, &
                            idmu, idgamma, idsat, idkappa, idalpha
  implicit none

  real :: rho_avg, rho_cgs, T_avg, delta_T, gamma, a0_radius, P, Scrit,    &
          dt, epsC, ndens_H
  integer, parameter :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10
  integer, parameter :: nTg_points = 10000
  real :: nTg(nTg_points)
  real :: T_range(nTg_points)
  integer :: idx, i, iunit
  logical :: id_exist
  real :: JKmuS(n_nucleation)

  rho_avg = 1.0d-13   ! in cgs
  T_avg = 1600
  delta_T = 400
  gamma = 5./3.  ! Adiabatic index
  a0_radius = 1.28d-8
  P = 2.6d7 ! in s
  Scrit = 1
  ! wind_CO_ratio = 3
  JKmuS(idgamma) = gamma
  JKmuS(idK0:idK3) = 0.
  
  do idx = 1, nTg_points
    nTg(idx) = real(idx - 1) / real(nTg_points - 1)
  end do

  T_range = T_avg + delta_T * cos(2*pi*(1-nTg))
  dt = P / nTg_points  ! Time step in seconds

  call set_abundances
  eps(iOx) = 9e-4  ! Adjusted O abundance to better match Gauger 1990
  eps(iC) = eps(iOx) * wind_CO_ratio
  ! eps(iH) = 1.0
  epsC = eps(iC)

  inquire(file='abundances_output.txt', exist=id_exist)
  if (id_exist) then
    call execute_command_line('rm -f abundances_output.txt')
  end if

  open(newunit=iunit,file='dustFormation_table.txt',status='replace')
  ! write(iunit,"(a)") '#   T   \Lambda_E(T) erg s^{-1} cm^3   \Lambda erg s^{-1} cm^{-3}'

  do i = 1, nTg_points
    rho_cgs = rho_avg * (T_range(i) / T_avg)**(1/(gamma - 1))
    ndens_H = rho_cgs / mass_per_H
    ! eps(iC) = abundi(icoolC) / ndens_H
    ! Call the dust formation routine for each temperature grid point
    ! if (i > 1) then
    !   JKmuS(i, :) = JKmuS(i-1, :)
    ! end if
    call evolve_chem(dt, T_range(i), rho_cgs, JKmuS(:))

    write(iunit,*)  T_range(i), JKmuS(idJstar)*ndens_H, JKmuS(idK0:), &
                   JKmuS(idK3) / (eps(iC) - eps(iOx))
    ! Jstar = 1, K0 = 2, K1 = 3, K2 = 4, K3 = 5, mu = 6,
    ! gamma = 7, sat = 8, kappa = 9, alpha = 10
  end do
  close(iunit)

end subroutine compute_dust_formation


subroutine total_dust_mass(time,npart,particlemass,xyzh)
 use part,           only:nucleation,idK3,idK0,idK1, idJstar
 use dust_formation, only:set_abundances, mass_per_H
 use physcon, only:atomic_mass_unit
 real, intent(in)               :: time,particlemass,xyzh(:,:)
 integer, intent(in)            :: npart
 integer                        :: i,ncols,j
 integer                        :: dump_number = 0
 real, dimension(2)             :: dust_mass
 character(len=17), allocatable :: columns(:)
 real, allocatable              :: temp(:) !npart
 real                           :: median,mass_factor,grain_size
 real, parameter :: a0 = 1.28e-4 !radius of a carbon atom in micron

 call set_abundances !initialize mass_per_H
 dust_mass = 0.
 ncols = 2
 print *,'size(nucleation,1) = ',size(nucleation,1)
 print *,'size(nucleation,2) = ',size(nucleation,2)
 allocate(columns(ncols),temp(npart))
 columns = (/'Dust mass [Msun]', &
             'median size [um]'/)
 j=0
 mass_factor = 12.*atomic_mass_unit*particlemass/mass_per_H
 do i = 1,npart
    if (.not. isdead_or_accreted(xyzh(4,i))) then
       dust_mass(1) = dust_mass(1) + nucleation(idK3,i) *mass_factor
       grain_size = a0*nucleation(idK1,i)/(nucleation(idK0,i)+1.0E-99) !in micron
       if (grain_size > a0) then
          j = j+1
          temp(j) = grain_size
       endif
    endif
 enddo

 call sort(temp,j)
 if (mod(j,2)==0) then !npart
    median = (temp(j/2)+temp(j/2+1))/2.0 !(temp(npart/2)+temp(npart/2+1))/2.0
 else
    median = (temp(j/2)+temp(j/2+1))/2.0 !temp(npart/2+1)
 endif

 dust_mass(2) = median

 call write_time_file('total_dust_mass_vs_time', columns, time, dust_mass, ncols, dump_number)
 !after execution of the analysis routine, a file named "total_dust_mass_vs_time.ev" appears
 deallocate(columns,temp)

end subroutine total_dust_mass

! ------------------------------------------------------
! The following subroutines are used by total_dust_mass
! ------------------------------------------------------

subroutine write_time_file(name_in, cols, time, data_in, ncols, num)
 !outputs a file over a series of dumps
 character(len=*), intent(in) :: name_in
 integer, intent(in)          :: ncols, num
 character(len=*), dimension(ncols), intent(in) :: cols
 character(len=20), dimension(ncols) :: columns
 character(len=40)             :: data_formatter, column_formatter
 character(len(name_in)+9)    :: file_name
 real, intent(in)             :: time
 real, dimension(ncols), intent(in) :: data_in
 integer                      :: i, unitnum

 write(column_formatter, "(a,I2.2,a)") "('#',2x,", ncols+1, "('[',a15,']',3x))"
 write(data_formatter, "(a,I2.2,a)") "(", ncols+1, "(2x,es18.11e2))"
 write(file_name,"(2a,i3.3,a)") name_in, '.ev'

 if (num == 0) then
    unitnum = 1000

    open(unit=unitnum,file=file_name,status='replace')
    do i=1,ncols
       write(columns(i), "(I2,a)") i+1, cols(i)
    enddo

    !set column headings
    write(unitnum, column_formatter) '1         time', columns(:)
    close(unit=unitnum)
 endif

 unitnum=1001+num

 open(unit=unitnum,file=file_name, position='append')

 write(unitnum,data_formatter) time, data_in(:ncols)

 close(unit=unitnum)

end subroutine write_time_file

! --------------------------------------------------------------------
! subroutine  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

subroutine  Sort(x, longitud)
 implicit  none
 integer, intent(in)                   :: longitud
 real, dimension(longitud), intent(inout) :: x
 integer                               :: i
 integer                               :: location

 do i = 1, longitud-1             ! except for the last
    location = findminimum(x, i, longitud)  ! find min from this to last
    call swap(x(i), x(location))  ! swap this and the minimum
 enddo
end subroutine Sort

! --------------------------------------------------------------------
! integer function  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

integer function  FindMinimum(x, Start, Fin)
 implicit  none
 integer, intent(in)                   :: start, fin
 real, dimension(Fin), intent(in) :: x
 real                            :: minimum
 integer                            :: location
 integer                            :: i

 minimum  = x(start)          ! assume the first is the min
 location = start             ! record its position
 do i = start+1, fin          ! start with next elements
    if (x(i) < minimum) then  !   if x(i) less than the min?
       minimum  = x(i)        !      yes, a new minimum found
       location = i                !      record its position
    endif
 enddo
 findminimum = location            ! return the position
end function FindMinimum

subroutine swap(a,b)
 real, intent(inout) :: a,b
 real                :: c

 c = a
 a = b
 b = c

end subroutine swap


end module analysis
