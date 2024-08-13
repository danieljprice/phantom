!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_tillotson
!
! EoS from Tillotson et al. (1962)
!
! :References: https://apps.dtic.mil/sti/pdfs/AD0486711.pdf
!
! Implementation from Brundage, A. 2013
! :DOI: 10.1016/j.proeng.2013.05.053
!
! :Owner:
!
! :Runtime parameters: None
!
! :Dependencies: 
!
 implicit none

contains
!-----------------------------------------------------------------------
!+
!  EoS from Tillotson et al. (1962) ; Implementation from Brundage (2013) and Ahrens and O'Keefe (1977)
!+
!-----------------------------------------------------------------------
subroutine equationofstate_tillotson(rho,energy,temperature,mu,X,Z,pressure,spsound,gamma)
 real, intent(in) :: rho,energy,temperature,mu,X,Z
 real, intent(out) :: pressure,spsound,gamma
! 
! Define: rho_0, E_0, V_0 (STP (std T and P) specific volume)
! Calc: rho (g/cm3), A (bulk modulus in kbar)
!
! Internal Energy, E (erg/g)
! E(rho,T) = E_C(rho) + E_T(rho,T) + E_e(rho,T)
! E_T = C_V * T  , where E_e ~ 0 (small relative to the other terms)
!
! --- Compressed state ---
! if rho >= rho_0 and E >= 0.
! Define: 
! a (Tillotson param), a = 0.5 (polytropic constant - 1, at high Temp)
! b (Tillotson param in kbar), b: defined sth (a+b) = STP Gruneisen param
! b, B, E_0 obtained by fitting to Thomas-Fermi calcs of EOS in 10^2 Mbar range and Hugoniot data.
! gamma = V(del P/ del E)_v , gamma_e = 2/3 (free electron gas), 1/2 (real gas) 
! eta (compression, rho/rho_0, and/or V_0/V), mu (strain, eta-1)
! 
! pressure(rho,E) = [ a + [ b / ( E / ( E_0 * eta**2 ) + 1. ) ] ] * rho*E + A*mu + B*mu**2
! 
! Cold curve: compression (4th order Runge-Kutta)
! d E_C / d rho = P1(rho,E_C) / rho**2 , (IC: rho = rho_0, E_C = 0) , P1 = pressure @ compressed state (C)
! 
! --- Expanded state --- [4x regions]
! Cold curve: expansion same method to solve as compression except RHS piecewise function of eqs 1-4 below
!
! 1: Cold expanded state:
! Define: rho_IV, E_IV (IV = incipient vaporisation) 
! if rho_0 > rho >= rho_IV and E <= E_IV
! pressure(rho,E) = (same eq as compressed)
! 
! 2: Hot expanded state: 
! Define: E_CV (CV = completely vaporised)
! if rho_0 > rho and E >= E_CV
! Define: alpha, beta, a, b, A, B
! pressure(rho,E) = a*rho*E + [ (b*rho*E)/(E/(E_0 * eta**2)+1) + A*mu*exp( -beta*[ (rho_0/rho) -1 ] )] * exp( -alpha*[ (rho_0/rho) -1 ]**2 )
! As density approaches 0, second term in expression disappears and the material approaches classical Thomas-Fermi limit (high pressures and energy during shock compression)
! P_e = rho*gamma_e *E_e [rho_0 >> rho, E >> E_CV], gamma_e = 2/3 (free electron gas), 1/2 (real gas)
! Pressure in 2. chosen to match Pressure in 1. at rho_0 (init)
! 
! 3: Mixing region between gas and liquid (solid):
! Define: P3 = pressure @ 2. , P2 = pressure @ 1.
! if rho_0 > rho > rho_IV and E_IV < E < E_CV
! pressure(rho,E) = [(E-E_IV)*P3 + (E_CV - E)*P2] / (E_CV - E_IV)
! 
! 4: Low energy expansion state:
! if rho < rho_IV and E < E_CV
! Define: a, b, A
! pressure(rho,E) = [ a + [ b / ( E / ( E_0 * eta**2 ) + 1 ) ] ] * rho*E + A*mu  (same eq as compressed but without B*mu**2 term)
! 
! ---Tillotson EOS params [new subroutine?]---
! From Benz and Asphaug 1999 [Table 2] DOI: 10.1006/icar.1999.6204
! ---Basalt---
! rho_0 = 2.7 [g/cm3]; A = 2.67e11 [erg/cm3]; B = 2.67e11[erg/cm3]; 
! E_0 = 4.87e12 [erg/g]; E_IV = 4.72e10 [erg/g]; E_CV = 1.82e11 [erg/g]; 
! a = 0.5; b = 1.5; alpha = 5.0; beta = 5.0
!
! ---Ice---
! rho_0 = 0.917 [g/cm3]; A = 9.47e10 [erg/cm3]; B = 9.47e10[erg/cm3]; 
! E_0 = 1.00e11 [erg/g]; E_IV = 7.73e9 [erg/g]; E_CV = 3.04e10 [erg/g]; 
! a = 0.3; b = 0.1; alpha = 10.0; beta = 5.0
! 
end subroutine equationofstate_tillotson

!----------------------------------------------------------------
!+
!  print eos information
!+
!----------------------------------------------------------------
subroutine eos_info_tillotson(iprint)
 integer, intent(in) :: iprint

 write(iprint,"(/,a)") ' Tillotson EoS'

end subroutine eos_info_tillotson

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
 subroutine read_options_eos_tillotson(name,valstring,imatch,igotall,ierr)
  use io, only:fatal
  character(len=*),  intent(in)  :: name,valstring
  logical,           intent(out) :: imatch,igotall
  integer,           intent(out) :: ierr
  integer,           save        :: ngot  = 0
  character(len=30), parameter   :: label = 'eos_tillotson'
 
  imatch  = .true.
  select case(trim(name))
!   case('irecomb')
!      read(valstring,*,iostat=ierr) irecomb
!      if ((irecomb < 0) .or. (irecomb > 3)) call fatal(label,'irecomb = 0,1,2,3')
!      ngot = ngot + 1
  case default
     imatch = .false.
  end select
 
  igotall = (ngot >= 1)
 
 end subroutine read_options_eos_tillotson

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_tillotson(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

!  call write_inopt(irecomb,'irecomb','recombination energy to include. 0=H2+H+He, 1=H+He, 2=He, 3=none',iunit)

end subroutine write_options_eos_tillotson

end module eos_tillotson