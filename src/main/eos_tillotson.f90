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
 real :: energy0 = 4.87e12 ! erg/g
 real :: eta, mu
 private :: rho0, aparam, bparam, A, B, alpha, beta, energy0, eta, mu

contains
!-----------------------------------------------------------------------
!+
!  EoS from Tillotson et al. (1962) ; Implementation from Benz and Asphaug (1999)
!+
!-----------------------------------------------------------------------
subroutine equationofstate_tillotson(rho,energy,pressure,spsound,gamma)
 real, intent(in) :: rho,energy
 real, intent(out) :: pressure,spsound,gamma
!
eta = rho / rho0
mu = eta - 1.

if (rho >= rho0 .and. energy >= 0.) then
  pressure = ( aparam + ( bparam / ( energy / ( energy0 * eta**2 ) + 1. ) ) ) * rho*energy + A*mu + B*mu**2
else
  pressure = aparam*rho*energy + ( (bparam*rho*energy)/(energy/(energy0 * eta**2) + 1.) + &
  A*mu*exp(-beta*( (rho0/rho) - 1.))) * exp(-alpha*( (rho0/rho) - 1.)**2)
endif

gamma = 0.5 ! gamma = V(del P/ del E)_v , gamma_e = 2/3 (free electron gas), 1/2 (real gas) 
spsound = sqrt(A / rho) ! wave speed estimate

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

 call write_inopt(rho0,'rho0','reference density',iunit)

end subroutine write_options_eos_tillotson

end module eos_tillotson