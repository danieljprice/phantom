!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_koyamainutsuka
!
! Koyama & Inutsuka (2002) cooling curve
!
! :References:
!   Koyama & Inutsuka (2002), ApJL 564, 97-100
!   Vazquez-Semadeni, et.al (2007), ApJ 657, 870-883
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, physcon, units
!
 implicit none
 public :: init_cooling_KI02, write_options_cooling_KI02, read_options_cooling_KI02
 public :: cooling_KoyamaInutsuka_explicit, cooling_KoyamaInutsuka_implicit

 integer, parameter :: maxt = 1000
 real :: rhov4_KI02(2,maxt)
 real :: LambdaKI_coef,GammaKI
 real :: KI02_rho_min_cgs = 1.0d-30              ! minimum density of the KI02 cooling curve
 real :: KI02_rho_max_cgs = 1.0d-14              ! maximum density of the KI02 cooling curve
 real :: KI02_rho_min,KI02_rho_max
 real, public :: GammaKI_cgs = 2.d-26        ! [erg/s] heating rate for Koyama & Inutuska cooling

 private

contains

!-----------------------------------------------------------------------
!+
!   Koyama & Inutsuka (2002) cooling
!+
!-----------------------------------------------------------------------
subroutine init_cooling_KI02(ierr)
 use io,      only:error
 use physcon, only:mass_proton_cgs
 use units,   only:utime,umass,udist
 integer, intent(out) :: ierr

 LambdaKI_coef = GammaKI_cgs*umass*utime**3/(mass_proton_cgs**2 * udist**5)
 GammaKI       = GammaKI_cgs*utime**3/(mass_proton_cgs*udist**2)
 call init_hv4table(ierr)
 if (ierr > 0) call error('init_cooling','Failed to create KI02 cooling table')

end subroutine init_cooling_KI02

!-----------------------------------------------------------------------
!+
!  create a h-v4 table based upon the cooling curve of KI02
!+
!-----------------------------------------------------------------------
subroutine init_hv4table(ierr)
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:unit_density,unit_velocity
 use eos,     only:gmw,gamma

 integer, intent(out) :: ierr

 integer              :: i,ctr
 real                 :: nrho0_min,nrho0_max,nrho,dnrho,dGammaKI,Lambda,dLambda
 real                 :: T,Tnew,Trat,fatT,faTdT
 logical              :: iterate
 logical              :: print_cc = .false. ! Print the cooling curve (for testing)

 !--Initialise densities
 KI02_rho_min = KI02_rho_min_cgs/unit_density
 KI02_rho_max = KI02_rho_max_cgs/unit_density
 nrho0_min    = KI02_rho_min_cgs/mass_proton_cgs
 nrho0_max    = KI02_rho_max_cgs/mass_proton_cgs
 dnrho        = (log10(nrho0_max) - log10(nrho0_min))/maxt
 !--Initialise additional variables
 dGammaKI     = 0.0
 ierr         = 0

 if (print_cc) open(unit=1031,file='coolingcurve.dat')

 !--Iterate (in cgs units)!
 T = 20000.
 do i = 1,maxt
    ctr     = 0
    iterate = .true.
    nrho    = 10**(log10(nrho0_min) + (i-1)*dnrho)
    do while ( iterate )
       Lambda  = 1.d7*exp(-1.184d5/(T+1.d3)) + 0.014*sqrt(T)*exp(-92./T) ! This is actually Lamda / Gamma
       dLambda = 0.007*exp(-92./T)*(T+184.)*T**(-1.5) + 1.184d12*exp(-1.184d5/(T+1.d3))*(T+1.d3)**(-2)
       fatT    =  Lambda*GammaKI_cgs*nrho -  GammaKI_cgs
       faTdT   = dLambda*GammaKI_cgs*nrho - dGammaKI
       Tnew    = abs(T - fatT/faTdT)
       Trat    = abs( 1.0 - T/Tnew )
       T       = Tnew
       ctr     = ctr + 1
       !--converged
       if (Trat < 1.0d-6) iterate = .false.
       !--failed to converge
       if (T < 0. .or. ctr > 2000) then
          iterate = .false.
          ierr    = 1
       endif
    enddo
    if (print_cc) write(1031,*) nrho,nrho*mass_proton_cgs,T,T*nrho,Lambda*GammaKI_cgs
    rhov4_KI02(1,i) = nrho
    rhov4_KI02(2,i) = T
 enddo
 if (print_cc) close(1031)

 !--Convert to useful values
 do i = 1,maxt
    rhov4_KI02(1,i) = rhov4_KI02(1,i)*mass_proton_cgs/unit_density                               ! number density (cm^-3) -> mass density (code units)
    rhov4_KI02(2,i) = kboltz*rhov4_KI02(2,i)/(gmw*mass_proton_cgs*(gamma-1.0))/unit_velocity**2  ! T -> internal energy (code units)
 enddo

end subroutine init_hv4table

!-----------------------------------------------------------------------
!+
!   Cooling rate as per Koyama & Inutuska (2002; eqns 4 & 5);
!   typos corrected as per Vazquez-Semadeni+ (2007)
!   This is for the explicit calculation
!   In equilibrium, n*LambdaKI = (rho/mp)*LambdaKI = GammaKI
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutsuka_explicit(rhoi,Tgas,dudti)
 real, intent(in)    :: rhoi,Tgas
 real, intent(inout) :: dudti
 real                :: LambdaKI

 ! Derivation to obtain correct units; used Koyama & Inutuska (2002) as the reference
 !LambdaKI = GammaKI_cgs * (1.d7*exp(-118400./(Tgas+1000))+0.014*sqrt(Tgas)*exp(-92./Tgas)) ! The cooling rate in erg cm^3/s = g cm^5/s^3
 !LambdaKI = LambdaKI/mass_proton_cgs**2                                                    ! units are now cm^5/(g s^3) ! since [u] = erg/g = cm^2/s^2
 !LambdaKI = LambdaKI*umass*utime**3/udist**5                                               ! convert to from cm^5/(g s^3) to code units
 !dudti    = dudti - LambdaKI*rhoi*fac                                                      ! multiply by rho (code) to get l^5/(m t^3) * m/l^3 = l^2/s^3 = [u]
 !
 !GammaKI = GammaKI_cgs                                                                     ! The heating rate in erg /s = g cm^2/s^3
 !GammaKI = GammaKI/mass_proton_cgs                                                         ! divide by proton mass.  Units are now g cm^2 / s^3 / g = cm^2/s^3
 !GammaKI = GammaKI*utime**3/udist**2                                                       ! convert from cm^2/s^3 to code units
 !dudti   = dudti + GammaKI                                                                 ! units and dependencies are correct

 LambdaKI = LambdaKI_coef*(1.d7*exp(-118400./(Tgas+1000.))+0.014*sqrt(Tgas)*exp(-92./Tgas))
 dudti    = dudti - LambdaKI*rhoi + GammaKI

end subroutine cooling_KoyamaInutsuka_explicit

!-----------------------------------------------------------------------
!+
!   Cooling rate as per Koyama & Inutuska (2002; eqns 4 & 5);
!   typos corrected as per Vazquez-Semadeni+ (2007)
!   This is the implicit method given by (5)-(6) in Vazquez-Semadeni+ (2007)
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutsuka_implicit(eni,rhoi,dt,dudti)
 use io,  only:fatal
 use eos, only:gamma,temperature_coef,gmw

 real, intent(in)    :: rhoi,eni,dt
 real, intent(out)   :: dudti

 integer             :: i,j,jm1
 real                :: ponrhoi,tempi,eni_equil,eni_final,deni,tau1,LambdaKI

 !--Determine the indicies surrounding the input h
 i = minloc(abs(rhov4_KI02(1,1:maxt)-rhoi), 1)
 if (i==1) then
    !print*, 'min density too large! extrapolating using two smallest densities'
    j = 2
 elseif (i==maxt) then
    !print*, 'max density too small! extrapolating using two largest densities'
    j = maxt
 elseif (rhov4_KI02(1,i-1) <= rhoi .and. rhoi <= rhov4_KI02(1,i  )) then
    j = i
 elseif (rhov4_KI02(1,i  ) <= rhoi .and. rhoi <= rhov4_KI02(1,i+1)) then
    j = i+1
 else
    j = 0 ! avoid compiler warning
    print*, rhoi,rhov4_KI02(1,i-1:i+1)
    call fatal('cooling_koyama','this should not happen')
 endif

 !--Calculate the equilibrium energy by linear interpolation
 jm1       = j - 1
 eni_equil = rhov4_KI02(2,j) + (rhov4_KI02(2,jm1)-rhov4_KI02(2,j))/(rhov4_KI02(1,jm1)-rhov4_KI02(1,j))*(rhoi-rhov4_KI02(1,j))

 !--Determine the inverse time require to radiate/acquire excess/deficit energy & Update energy
 ponrhoi  = (gamma-1.)*eni
 tempi    = temperature_coef*gmw*ponrhoi
 LambdaKI = LambdaKI_coef*(1.d7*exp(-118400./(tempi+1000.))+0.014*sqrt(tempi)*exp(-92./tempi))
 dudti    = LambdaKI*rhoi - GammaKI
 deni     = eni - eni_equil

 if (abs(deni) > 0.) then
    ! in both limits, this will approach the correct value
    tau1      = abs(dudti/deni)
    eni_final = eni_equil + deni*exp(-dt*tau1)
    dudti     = -(eni - eni_final)/dt
 else
    ! in the unlikly chance deni = 0
    dudti = -dudti
 endif

end subroutine cooling_KoyamaInutsuka_implicit

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling_KI02(iunit)
 !use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

end subroutine write_options_cooling_KI02

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling_KI02(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 imatch  = .true.
 igotall = .true. ! nothing to read

end subroutine read_options_cooling_KI02

end module cooling_koyamainutsuka
