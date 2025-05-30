!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module physcon
!
! Physical and mathematical constants (as in sphNG)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
!
!--Mathematical constants (keep these in variable precision to avoid warnings)
!
 real, parameter :: pi       =  3.1415926536d0
 real, parameter :: twopi    =  6.2831853072d0
 real, parameter :: fourpi   = 12.5663706144d0
 real, parameter :: piontwo  =  1.5707963268d0
 real, parameter :: rpiontwo =  1.2533141373d0          !square root of (Pi/2)
 real, parameter :: roottwo  =  1.4142135624d0
 real, parameter :: deg_to_rad = pi/180.
!
!--Physical constants
!
 real(kind=8), parameter :: c = 2.997924d10                     !Speed of light            cm/s
 real(kind=8), parameter :: gg = 6.672041d-8                    !Gravitational constant    dyn cm^2 g^-2
 !                          cm^3 s^-2 g^-1

 real(kind=8), parameter :: Rg = 8.31446261815324d7             !Gas constant              erg/K/g
 real(kind=8), parameter :: cgsmu0 = 4.*pi
 real(kind=8), parameter :: mass_electron_cgs = 9.10938291d-28  !Electron mass             g
 real(kind=8), parameter :: mass_proton_cgs = 1.67262158d-24    !Proton mass               g
 real(kind=8), parameter :: atomic_mass_unit = 1.660538921d-24  !Atomic mass unit          g
 real(kind=8), parameter :: cross_section_H2_cgs = 2.367d-15    !Hydrogen molecule cs      cm^-2
 real(kind=8), parameter :: radconst = 7.5646d-15               !Radiation constant        erg cm^-3 K^-4
 real(kind=8), parameter :: kboltz = 1.38066d-16                !Boltzmann constant        erg/K
 real(kind=8), parameter :: kb_on_mh = kboltz/mass_proton_cgs   !kB/m_H                    erg/K/g
 real(kind=8), parameter :: eV     = 1.60219d-12                !electron volt             erg
 real(kind=8), parameter :: qe     = 4.8032068d-10              !charge on electron        esu
 real(kind=8), parameter :: planckh  =   6.6260755d-27          !Planck's Constant         erg.s
 real(kind=8), parameter :: planckhbar = 1.05457266d-27         !Planck's Constant/(2pi)   erg.s
 real(kind=8), parameter :: thomcs     = 6.6525d-25             !Thomson cross section     cm^2
 real(kind=8), parameter :: finestr    = 7.2974d-3              !Fine structure constant   unitless
 real(kind=8), parameter :: steboltz   = 5.67051d-5             !Stefan-Boltzmann constant erg cm^-2K^-4 s^-1
 real(kind=8), parameter :: avogadro   = 6.0221408577d23        !Avogadro's number         mole^-1
 real(kind=8), parameter :: Ro         = 3.00000000             !Rossby number without dimension
 real(kind=8), parameter :: patm       = 1.013250d6             !Standard atmospheric pressure in cgs

!--Astronomical constants (cgs units)
!
!--Solar mass and radius
!
 real(kind=8), parameter :: solarm = 1.9891d33                  !Mass of the Sun           g
 real(kind=8), parameter :: solarr = 6.959500d10                !Radius of the Sun         cm
 real(kind=8), parameter :: solarl = 3.9d33                     !Luminosity of the Sun     erg/s
!
!--Earth mass and radius
!
 real(kind=8), parameter :: earthm = 5.979d27                   !Mass of the Earth         g
 real(kind=8), parameter :: earthr = 6.371315d8                 !Radius of the Earth       cm
 real(kind=8), parameter :: jupiterm = 1.89813d30               !Mass of Jupiter           g
 real(kind=8), parameter :: jupiterr = 7.1492e9                 !Equatorial radius Jupiter cm
 real(kind=8), parameter :: ceresm = 8.958d23                   !Mass of Ceres             g
 real(kind=8), parameter :: kg = 1.d3
 real(kind=8), parameter :: gram = 1.d0
!
!--Distance scale
!
 real(kind=8), parameter :: au = 1.496d13                       !Astronomical unit         cm
 real(kind=8), parameter :: ly = 9.4605d17                      !Light year                cm
 real(kind=8), parameter :: pc = 3.086d18                       !Parsec                    cm
 real(kind=8), parameter :: kpc = 3.086d21                      !Kiloparsec                cm
 real(kind=8), parameter :: Mpc = 3.086d24                      !Megaparsec                cm
 real(kind=8), parameter :: km = 1.d5                           !Kilometer                 cm
 real(kind=8), parameter :: metre = 1.d2                        !Metre                     cm
 real(kind=8), parameter :: cm = 1.d0                           !Centimetre                cm
 real(kind=8), parameter :: mm = 0.1d0                          !Millimetre                cm
 real(kind=8), parameter :: micron = 1.d-4                      !Micron                    cm
 real(kind=8), parameter :: nm = 1.d-7                          !Nanometre                 cm
 real(kind=8), parameter :: angstrom = 1.d-8                    !Angstrom                  cm
!
!--Time scale
!
 real(kind=8), parameter :: seconds = 1.d0
 real(kind=8), parameter :: minutes = 6.0d1
 real(kind=8), parameter :: hours = 3.6d3
 real(kind=8), parameter :: days = 8.64d4
 real(kind=8), parameter :: years = 3.1556926d7
 real(kind=8), parameter :: myr   = 3.1556926d13
!
!--Energy conversion
!
 real(kind=8), parameter :: eVtoK = 1.1604519d4                 !Degrees kelvin per eV     K/eV

end module physcon
