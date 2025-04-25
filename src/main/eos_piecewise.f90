!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_piecewise
!
! Implements piecewise polytropic equation of state, e.g. for neutron stars
!
! :References: Read et al. (2009)
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - gamma0pwp   : *adiabatic index 0 (piecewise polytropic eos)*
!   - gamma1pwp   : *adiabatic index 1 (piecewise polytropic eos)*
!   - gamma2pwp   : *adiabatic index 2 (piecewise polytropic eos)*
!   - gamma3pwp   : *adiabatic index 3 (piecewise polytropic eos)*
!   - p1pwp       : *pressure at cutoff density rhocrit1pwp (piecewise polytropic eos)*
!   - rhocrit0pwp : *critical density 0 in g/cm^3 (piecewise polytropic eos)*
!   - rhocrit1pwp : *critical density 1 in g/cm^3 (piecewise polytropic eos)*
!   - rhocrit2pwp : *critical density 2 in g/cm^3 (piecewise polytropic eos)*
!
! :Dependencies: infile_utils, io, units
!
 use units, only:unit_density,unit_pressure
 implicit none

 !--Default initial parameters for piecewise polytrope Eos
 integer, parameter :: maxEOSopt =  4 ! maximum number of piecewise polytrope defaults
 real(kind=8) :: rhocrit0pwpcgs = 2.62780d12
 real(kind=8) :: rhocrit1pwpcgs = 5.01187d14
 real(kind=8) :: rhocrit2pwpcgs = 1.0d15
 real(kind=8) :: p1pwpcgs       = 2.46604d34
 real :: gamma0pwp      = 5./3.
 real :: gamma1pwp      = 3.166
 real :: gamma2pwp      = 3.573
 real :: gamma3pwp      = 3.281
 real :: rhocrit0pwp,rhocrit1pwp,rhocrit2pwp,p0pwp,p1pwp,p2pwp,k0pwp,k1pwp,k2pwp,k3pwp

 public :: init_eos_piecewise,init_eos_piecewise_preset
 public :: get_eos_piecewise,eos_info_piecewise
 public :: write_options_eos_piecewise,read_options_eos_piecewise
 public :: gamma_pwp,get_dPdrho_piecewise

 private

contains

!-----------------------------------------------------------------------
!+
!  The default piecewise polytrope options, as per Read et al (2009)
!  The unlisted values are common to all options; all values are in cgs
!  The array is
!  pw(i,:) = (/ prescrit,gamma1,gamma2,gamma3 /)
!  pw(:,j) = (/ ARP3,SLy,MS1,ENG/)
!+
!-----------------------------------------------------------------------
subroutine init_eos_piecewise_preset(EOSopt)
 integer, parameter :: numparam =  4 ! number of parameters governing the piecewise polytrope
 integer, intent(in) :: EOSopt
 real :: pw(maxEOSopt,numparam)
 !
 ! Define the default options
 !
 pw(1,:)  = (/ 10**34.392, 3.166, 3.573, 3.281 /)
 pw(2,:)  = (/ 10**34.384, 3.005, 2.988, 2.851 /)
 pw(3,:)  = (/ 10**34.858, 3.224, 3.033, 1.325 /)
 pw(4,:)  = (/ 10**34.437, 3.514, 3.130, 3.168 /)
 !
 ! Choose the default option
 !
 p1pwpcgs  = pw(EOSopt,1)
 gamma1pwp = pw(EOSopt,2)
 gamma2pwp = pw(EOSopt,3)
 gamma3pwp = pw(EOSopt,4)

end subroutine init_eos_piecewise_preset

!-----------------------------------------------------------------------
!+
!  Initialise the equation of state
!+
!-----------------------------------------------------------------------
subroutine init_eos_piecewise(ierr)
 integer, intent(out) :: ierr
 !
 !--calculate initial variables for the piecewise polytrope equation of state
 !
 if (unit_density <= 0.0 .or. unit_pressure <= 0.0) then
    ierr = 3
    return
 endif
 rhocrit0pwp = real(rhocrit0pwpcgs/unit_density)
 rhocrit1pwp = real(rhocrit1pwpcgs/unit_density)
 rhocrit2pwp = real(rhocrit2pwpcgs/unit_density)
 p1pwp       = real(p1pwpcgs/unit_pressure)
 k1pwp       = p1pwp/rhocrit1pwp**gamma1pwp
 k2pwp       = p1pwp/rhocrit1pwp**gamma2pwp
 p2pwp       = k2pwp*rhocrit2pwp**gamma2pwp
 k3pwp       = p2pwp/rhocrit2pwp**gamma3pwp
 k0pwp       = k1pwp/(rhocrit0pwp**(gamma0pwp-gamma1pwp))
 p0pwp       = k0pwp*rhocrit0pwp**gamma0pwp

end subroutine init_eos_piecewise

!-----------------------------------------------------------------------
!+
!  Main eos routine: calculates pressure at a given density
!+
!-----------------------------------------------------------------------
subroutine get_eos_piecewise(rhoi,ponrhoi,spsoundi,gammai)
 real, intent(in)  :: rhoi
 real, intent(out) :: ponrhoi,spsoundi,gammai

 if (rhoi < rhocrit0pwp) then
    gammai  = gamma0pwp
    ponrhoi = k0pwp*rhoi**(gamma0pwp-1.)
 elseif (rhoi < rhocrit1pwp) then
    gammai  = gamma1pwp
    ponrhoi = k1pwp*rhoi**(gamma1pwp-1.)
 elseif (rhoi < rhocrit2pwp) then
    gammai  = gamma2pwp
    ponrhoi = k2pwp*rhoi**(gamma2pwp-1.)
 else
    gammai  = gamma3pwp
    ponrhoi = k3pwp*rhoi**(gamma3pwp-1.)
 endif
 spsoundi = sqrt(gammai*ponrhoi)

end subroutine get_eos_piecewise

!-----------------------------------------------------------------------
!+
!  Get gamma for thermal energy calculations when using the
!  piecewise polytrope
!+
!-----------------------------------------------------------------------
real function gamma_pwp(rhoi)
 real, intent(in) :: rhoi

 if (rhoi < rhocrit0pwp) then
    gamma_pwp = gamma0pwp
 elseif (rhoi < rhocrit1pwp) then
    gamma_pwp = gamma1pwp
 elseif (rhoi < rhocrit2pwp) then
    gamma_pwp = gamma2pwp
 else
    gamma_pwp = gamma3pwp
 endif

end function gamma_pwp

!-----------------------------------------------------------------------
!+
!  Calculates derivative of pressure with respect to density
!  at a given density
!+
!-----------------------------------------------------------------------
real function get_dPdrho_piecewise(rho) result(get_dPdrho)
 real, intent(in)  :: rho
 real              :: rhocrit0pwp,rhocrit1pwp,rhocrit2pwp,presscrit
 real              :: polyk0,polyk1,polyk2,polyk3
 real              :: gamma,polyk

 rhocrit0pwp = real(rhocrit0pwpcgs/unit_density)
 rhocrit1pwp = real(rhocrit1pwpcgs/unit_density)
 rhocrit2pwp = real(rhocrit2pwpcgs/unit_density)
 presscrit   = real(p1pwpcgs/unit_pressure)
 polyk1      = presscrit/rhocrit1pwp**gamma1pwp
 polyk2      = presscrit/rhocrit1pwp**gamma2pwp
 polyk3      = polyk2*rhocrit2pwp**(gamma2pwp-gamma3pwp)
 polyk0      = polyk1*rhocrit0pwp**(gamma1pwp-gamma0pwp)

 if (rho < rhocrit0pwp) then
    gamma = 5./3.
    polyk = polyk0
 elseif (rho < rhocrit1pwp) then
    gamma = gamma1pwp
    polyk = polyk1
 elseif (rho < rhocrit2pwp) then
    gamma = gamma2pwp
    polyk = polyk2
 else
    gamma = gamma3pwp
    polyk = polyk3
 endif
 get_dPdrho = gamma * polyk * rho**(gamma-1.0)

end function get_dPdrho_piecewise

!-----------------------------------------------------------------------
!+
!  print information about the equation of state parameters
!+
!-----------------------------------------------------------------------
subroutine eos_info_piecewise(iprint)
 integer, intent(in) :: iprint

 write(iprint,"(/,a,3(es10.3),a,4(es10.3))") ' Piecewise polytropic eq of state (code units) : rhocrit = '&
                                                 ,rhocrit0pwp,rhocrit1pwp,rhocrit2pwp, '; K = ',k0pwp,k1pwp,k2pwp,k3pwp
 write(iprint,"(  a,3(es10.3)            )") ' Piecewise polytropic eq of state (g/cm^3)     : rhocrit = '&
                                                 ,rhocrit0pwp*unit_density,rhocrit1pwp*unit_density,rhocrit2pwp*unit_density

end subroutine eos_info_piecewise

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_piecewise(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(rhocrit0pwpcgs,'rhocrit0pwp','critical density 0 in g/cm^3 (piecewise polytropic eos)',iunit)
 call write_inopt(rhocrit1pwpcgs,'rhocrit1pwp','critical density 1 in g/cm^3 (piecewise polytropic eos)',iunit)
 call write_inopt(rhocrit2pwpcgs,'rhocrit2pwp','critical density 2 in g/cm^3 (piecewise polytropic eos)',iunit,exp=.true.)
 call write_inopt(gamma0pwp,'gamma0pwp','adiabatic index 0 (piecewise polytropic eos)',iunit)
 call write_inopt(gamma1pwp,'gamma1pwp','adiabatic index 1 (piecewise polytropic eos)',iunit)
 call write_inopt(gamma2pwp,'gamma2pwp','adiabatic index 2 (piecewise polytropic eos)',iunit)
 call write_inopt(gamma3pwp,'gamma3pwp','adiabatic index 3 (piecewise polytropic eos)',iunit)
 call write_inopt(p1pwpcgs,'p1pwp','pressure at cutoff density rhocrit1pwp (piecewise polytropic eos)',iunit)

end subroutine write_options_eos_piecewise

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_eos_piecewise(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot  = 0
 character(len=30), parameter  :: label = 'eos_piecewise'

 imatch  = .true.
 select case(trim(name))
 case('rhocrit0pwp')
    read(valstring,*,iostat=ierr) rhocrit0pwpcgs
    if (rhocrit0pwpcgs <= 0.) call fatal(label,'rhocrit0pwp <= 0')
    ngot = ngot + 1
 case('rhocrit1pwp')
    read(valstring,*,iostat=ierr) rhocrit1pwpcgs
    if (rhocrit1pwpcgs <= 0.) call fatal(label,'rhocrit1pwp <= 0')
    ngot = ngot + 1
 case('rhocrit2pwp')
    read(valstring,*,iostat=ierr) rhocrit2pwpcgs
    if (rhocrit2pwpcgs <= 0.) call fatal(label,'rhocrit2pwp <= 0')
    ngot = ngot + 1
 case('gamma0pwp')
    read(valstring,*,iostat=ierr) gamma0pwp
    if (gamma0pwp <= 0.) call fatal(label,'gamma0pwp < 1.0')
    ngot = ngot + 1
 case('gamma1pwp')
    read(valstring,*,iostat=ierr) gamma1pwp
    if (gamma1pwp < 1.) call fatal(label,'gamma1pwp < 1.0')
    ngot = ngot + 1
 case('gamma2pwp')
    read(valstring,*,iostat=ierr) gamma2pwp
    if (gamma2pwp < 1.) call fatal(label,'gamma2pwp < 1.0')
    ngot = ngot + 1
 case('gamma3pwp')
    read(valstring,*,iostat=ierr) gamma3pwp
    if (gamma3pwp < 1.) call fatal(label,'gamma3pwp < 1.0')
    ngot = ngot + 1
 case('p1pwp')
    read(valstring,*,iostat=ierr) p1pwpcgs
    if (p1pwpcgs <= 0.) call fatal(label,'p1pwp <= 0.0')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 8)

end subroutine read_options_eos_piecewise

end module eos_piecewise
