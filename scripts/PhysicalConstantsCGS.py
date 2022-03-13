#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#----------------------------------
# Python translation of physcon.f90
#----------------------------------

# Physical and mathematical constants (as in sphNG)


#
#--Mathematical constants
#
pi       =  3.1415926536e0
twopi    =  6.2831853072e0
fourpi   = 12.5663706144e0
piontwo  =  1.5707963268e0
rpiontwo =  1.2533141373e0          #square root of (Pi/2)
roottwo  =  1.4142135624e0
#
#--Physical constants
#
c = 2.997924e10                     #Speed of light            cm/s
gg = 6.672041e-8                    #Gravitational constant    dyn cm^2 g^-2, cm^3 s^-2 g^-1

Rg = 8.31446261815324e7             #Gas constant              erg/K/g
cgsmu0 = 4.*pi
mass_electron_cgs = 9.10938291e-28  #Electron mass             g
mass_proton_cgs = 1.67262158e-24    #Proton mass               g
atomic_mass_unit = 1.660538921e-24  #Atomic mass unit          g
cross_section_H2_cgs = 2.367e-15    #Hydrogen molecule cs      cm^-2
radconst = 7.5646e-15               #Radiation constant        erg cm^-3 K^-4
kboltz = 1.38066e-16                #Boltzmann constant        erg/K
kb_on_mh = kboltz/mass_proton_cgs   #kB/m_H                    erg/K/g
eV     = 1.60219e-12                #electron volt             erg
qe     = 4.8032068e-10              #charge on electron        esu
planckh  =   6.6260755e-27          #Planck's Constant         erg.s
planckhbar = 1.05457266e-27         #Planck's Constant/(2pi)   erg.s
thomcs     = 6.6525e-25             #Thomson cross section     cm^2
finestr    = 7.2974e-3              #Fine structure constant   unitless
steboltz   = 5.67051e-5             #Stefan-Boltzmann constant erg cm^-2K^-4 s^-1
avogadro   = 6.0221408577e23        #Avogadro's number         mole^-1
Ro         = 3.00000000             #Rossby number without dimension
#
#--Astronomical constants (cgs units)
#
#--Solar mass and radius
#
solarm = 1.9891e33                  #Mass of the Sun           g
solarr = 6.959500e10                #Radius of the Sun         cm
solarl = 3.9e33                     #Luminosity of the Sun     erg/s
#
#--Earth mass and radius
#
earthm = 5.979e27                   #Mass of the Earth         g
earthr = 6.371315e8                 #Radius of the Earth       cm
jupiterm = 1.89813e30               #Mass of Jupiter           g
ceresm = 8.958e23                   #Mass of Ceres             g
gram = 1.e0
#
#--Distance scale
#
au = 1.496e13                       #Astronomical unit         cm
ly = 9.4605e17                      #Light year                cm
pc = 3.086e18                       #Parsec                    cm
kpc = 3.086e21                      #Kiloparsec                cm
Mpc = 3.086e24                      #Megaparsec                cm
km = 1.e5                           #Kilometer                 cm
cm = 1.e0                           #Centimetre                cm
mm = 0.1e0                          #Millimetre                cm
micron = 1.e-4                      #Micron                    cm
nm = 1.e-7                          #Nanometre                 cm
angstrom = 1.e-8                    #Angstrom                  cm
#
#--Time scale
#
seconds = 1.e0
minutes = 6.0e1
hours = 3.6e3
days = 8.64e4
years = 3.1556926e7
#
#--Energy conversion
#
eVtoK = 1.1604519e4                 #Degrees kelvin per eV     K/eV


# Additionnal Constants not in physcon
patm = 1.013250e6                   # 1 atmosphere in cgs units
bar = 1.e6                          # dyne/cm²
golden_number = 1.618033988749895   #(1+sqrt(5))/2
