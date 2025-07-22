!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_barotropic_simple
!
! Implements simple barotropic equation of state, e.g. for star formation
!
! :References: Originally from Hennebelle papers
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - rho0cgs : *critical density 0 in g/cm^3 (barotropic eos)*
!   - rhocritcgs : *critical density 1 in g/cm^3 (barotropic eos)*
!   - gamma1 : *adiabatic index 1 (barotropic eos)*
!
! :Dependencies: infile_utils, io, units
!
 use physcon, only:solarm,pi
 implicit none
 !--Default initial parameters for Barotropic Eos
 real :: rho0cgs    = solarm/(4./3.*pi*4e16**3)
 real :: rhocritcgs = 1.e-13
 real :: gamma1     = 7./5.

 public :: init_eos_barotropic_simple
 public :: get_eos_barotropic_simple,eos_info_barotropic_simple
 public :: write_options_eos_barotropic_simple,read_options_eos_barotropic_simple
 public :: gamma_barotropic_simple

 real :: rho0,rhocrit

 private

contains

!-----------------------------------------------------------------------
!+
!  Initialise the equation of state
!+
!-----------------------------------------------------------------------
subroutine init_eos_barotropic_simple(ierr)
 use units, only:unit_density
 integer, intent(out) :: ierr
 !
 !--calculate initial variables for the barotropic equation of state
 !
 if (unit_density <= 0.) then
    ierr = 3
    return
 endif

 ! Convert to code units, and calculate constants
 rho0 = real(rho0cgs/unit_density)
 rhocrit = real(rhocritcgs/unit_density)

end subroutine init_eos_barotropic_simple

!-----------------------------------------------------------------------
!+
!  Main eos routine: calculates pressure at a given density
!+
!-----------------------------------------------------------------------
subroutine get_eos_barotropic_simple(rhoi,polyk,polyk2,ponrhoi,spsoundi,gammai)
 real, intent(in)  :: rhoi,polyk,polyk2
 real, intent(out) :: ponrhoi,spsoundi,gammai

 if (rhoi < rho0) then
    gammai  = 1.0
    ponrhoi = polyk * rho0/rhoi ! pressure = polyk * rho0
 else
    ponrhoi = polyk * (1. + (rhoi/rhocrit)**(gamma1-1.))
    gammai = 1. + gamma1*(rhoi/rhocrit)**(gamma1-1.)
 endif
 spsoundi = sqrt(gammai*ponrhoi)

end subroutine get_eos_barotropic_simple

!-----------------------------------------------------------------------
!+
!  Get gamma for thermal energy calculations when using the
!  piecewise polytrope
!+
!-----------------------------------------------------------------------
real function gamma_barotropic_simple(rhoi) result(gammai)
 real, intent(in) :: rhoi

 gammai = 1.0
 !call get_eos_barotropic_simple(rhoi,polyk,polyk2,ponrhoi,spsoundi,gammai)

end function gamma_barotropic_simple

!-----------------------------------------------------------------------
!+
!  print information about the equation of state parameters
!+
!-----------------------------------------------------------------------
subroutine eos_info_barotropic_simple(polyk,polyk2,iprint)
 use units, only:unit_velocity,unit_density
 real,    intent(in) :: polyk,polyk2
 integer, intent(in) :: iprint

 write(iprint,"(/,a)") ' Barotropic equation of state with continuous sound speed profile'
 write(iprint,"(a)")   '    P = K (rho + (rho/rho_c)^gamma)'
 write(iprint,"(a)")   ' where rho_c is the critical density and gamma is the adiabatic index.'
 write(iprint,"(a,es10.3,a)") ' rho0 = ',rho0*unit_density,' g/cm^3'
 write(iprint,"(a,es10.3,a)") ' rhocrit = ',rhocrit*unit_density,' g/cm^3'
 write(iprint,"(a,es10.3)")   ' gamma1 = ',gamma1
 write(iprint,"(a,es10.3,a)") ' isothermal sound speed = ',sqrt(polyk)*unit_velocity,' cm/s'
 write(iprint,"(a,es10.3,a,/)") ' sound speed in low density region = ',sqrt(polyk2)*unit_velocity,' cm/s'

end subroutine eos_info_barotropic_simple

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_barotropic_simple(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(rho0cgs,'rho0','critical density 0 in g/cm^3 (barotropic eos)',iunit)
 call write_inopt(rhocritcgs,'rhocrit','critical density 1 in g/cm^3 (barotropic eos)',iunit)
 call write_inopt(gamma1,'gamma1','adiabatic index 1 (barotropic eos)',iunit)

end subroutine write_options_eos_barotropic_simple

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_eos_barotropic_simple(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot  = 0
 character(len=30), parameter  :: label = 'eos_barotropic_simple'

 imatch  = .true.
 select case(trim(name))
 case('rho0')
    read(valstring,*,iostat=ierr) rho0cgs
    ngot = ngot + 1
 case('rhocrit')
    read(valstring,*,iostat=ierr) rhocritcgs
    ! if (rhocrit0cgs <= 0.) call fatal(label,'rhocrit0 <= 0')  ! This region can be 0 if the warm medium is undefined
    ngot = ngot + 1
 case('gamma1')
    read(valstring,*,iostat=ierr) gamma1
    if (gamma1 < 1.) call fatal(label,'gamma1 < 1.0')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 2)

end subroutine read_options_eos_barotropic_simple  

end module eos_barotropic_simple
