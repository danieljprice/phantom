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
! Implementation from Benz and Asphaug (1999)
!
! :Owner:
!
! :Runtime parameters: None
!
! :Dependencies: 
!
 implicit none
 real :: rho0 = 2.7 ! g/cm^3 zero-pressure density (Basalt) from Benz & Asphaug 1999
! aparam, bparam, A, B, energy0 material-dependent Tillotson parameters (Basalt)
 real :: aparam = 0.5 , bparam = 1.5 , alpha = 5. , beta = 5. 
 real :: A = 2.67e11 , B = 2.67e11 ! erg/cm^3
 real :: u0 = 4.87e12 ! erg/g
 private :: rho0, aparam, bparam, A, B, alpha, beta, u0

contains
!-----------------------------------------------------------------------
!+
!  EoS from Tillotson et al. (1962) ; Implementation from Benz and Asphaug (1999)
!+
!-----------------------------------------------------------------------
subroutine init_eos_tillotson(ierr)
 integer, intent(out) :: ierr

 ierr = 0

end subroutine init_eos_tillotson

subroutine equationofstate_tillotson(rho,u,pressure,spsound,gamma)
 real, intent(in) :: rho,u
 real, intent(out) :: pressure,spsound,gamma
 real :: eta, mu, neu, omega, spsoundmin

eta = rho / rho0
mu = eta - 1.
neu = (1. / eta) - 1.  ! Kegerreis et al. 2020
omega = (u / (u0*eta**2) + 1.)

if (rho >= rho0 .and. u >= 0.) then ! Condensed state
  ! pressure = ( aparam + ( bparam / ( u / ( u0 * eta**2 ) + 1. ) ) ) * rho*u + A*mu + B*mu**2
  pressure = ( aparam + ( bparam / omega ) ) * rho*u + A*mu + B*mu**2
  if (pressure <= 0.) then
    pressure = 0.
  endif
  ! sound speed from ! Kegereis et al. 2020
  spsound = sqrt(pressure / rho * (1. + aparam + (bparam / omega )) + &
                  (bparam*(omega - 1.) / omega**2) * (2.*u - pressure/rho)) + &
                  ((1./rho)*(A + B*((eta**2) -1.)))
  spsoundmin = sqrt(A / rho) ! wave speed estimate
  if (spsound < spsoundmin) then
    spsound = spsoundmin
  endif
else ! Expanded state
  ! pressure = aparam*rho*u + ( (bparam*rho*u)/(u/(u0 * eta**2) + 1.) + &
            !  A*mu*exp(-beta*( (rho0/rho) - 1.))) * exp(-alpha*( (rho0/rho) - 1.)**2)
  pressure = aparam*rho*u + ( ((bparam*rho*u)/omega) + &
             A*mu*exp(-beta*neu)) * exp(-alpha*(neu**2))
  if (pressure <= 0.) then
    pressure = 0.
  endif
  spsound = sqrt( (pressure/rho)*( 1. + aparam + ((bparam / omega )*exp(-alpha*neu**2))) + &
                 ( ((bparam*rho*u)/((omega**2) * (eta**2))) * ((1/(u0*rho))*(2.*u - pressure/rho) &
                  + ((2.*alpha*eta*omega)/rho0)) + ((A/rho0)*(1.+ (mu/eta**2)*(beta*2.*alpha - eta))) &
                  * exp(-beta*neu) )*exp(-alpha*neu**2) )
  spsoundmin = sqrt(A / rho) ! wave speed estimate
  if (spsound < spsoundmin) then
    spsound = spsoundmin
  endif
endif
! Hybrid Pressure (inbetween)
! uiv = 
! ucv =
! frac = ((u / uiv) - 1.) / ((ucv/uiv) - 1.)
! pressure = 
! ce = spsound @ expanded
! cc = spsound @ condensed
! spsound = sqrt(( ((u - uiv)* ce**2) + ((ucv - u)*cc**2))/(ucv - uiv))

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
  case('rho0')
      read(valstring,*,iostat=ierr) rho0
      if ((rho0 < 0.)) call fatal(label,'rho0 < 0')
      ngot = ngot + 1
  case('aparam')
      read(valstring,*,iostat=ierr) aparam
      if ((aparam < 0.)) call fatal(label,'aparam < 0')
      ngot = ngot + 1
  case('bparam')
      read(valstring,*,iostat=ierr) bparam
      if ((bparam < 0.)) call fatal(label,'bparam < 0')
      ngot = ngot + 1
  case('A')
      read(valstring,*,iostat=ierr) A
      if ((A < 0.)) call fatal(label,'A < 0')
      ngot = ngot + 1
  case('B')
      read(valstring,*,iostat=ierr) B
      if ((B < 0.)) call fatal(label,'B < 0')
      ngot = ngot + 1
  case('alpha_t')
      read(valstring,*,iostat=ierr) alpha
      if ((alpha < 0.)) call fatal(label,'alpha < 0')
      ngot = ngot + 1
  case('beta_t')
      read(valstring,*,iostat=ierr) beta
      if ((beta < 0.)) call fatal(label,'beta < 0')
      ngot = ngot + 1
  case('u0')
      read(valstring,*,iostat=ierr) u0
      if ((u0 < 0.)) call fatal(label,'u0 < 0')
      ngot = ngot + 1
  case default
     imatch = .false.
  end select
 
  igotall = (ngot >= 8)
 
 end subroutine read_options_eos_tillotson

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_tillotson(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit
 
 call write_inopt(rho0,'rho0','reference density g/cm^3',iunit)
 call write_inopt(aparam,'aparam','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(bparam,'bparam','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(A,'A','material-dependent Tillotson parameter, erg/cm^3',iunit)
 call write_inopt(B,'B','material-dependent Tillotson parameter, erg/cm^3',iunit)
 call write_inopt(alpha,'alpha_t','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(beta,'beta_t','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(u0,'u0','material-dependent Tillotson parameter, erg/g',iunit)

end subroutine write_options_eos_tillotson

end module eos_tillotson