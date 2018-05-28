!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: eos
!
!  DESCRIPTION:
!  This module contains stuff to do with the equation of state
!  Curent options:
!     1 = isothermal eos
!     2 = adiabatic/polytropic eos
!     3 = eos for a locally isothermal disc as in Lodato & Pringle (2007)
!     6 = eos for a locally isothermal disc as in Lodato & Pringle (2007),
!         centered on a sink particle
!     7 = z-dependent locally isothermal eos
!     8 = Barotropic eos
!     9 = Piecewise polytrope
!    10 = MESA EoS
!    11 = isothermal eos with zero pressure
!    14 = locally isothermal prescription from Farris et al. (2014) for binary system
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    X           -- hydrogen mass fraction
!    drhocrit    -- transition size between rhocrit0 & 1 (fraction of rhocrit0; barotropic eos)
!    gamma0pwp   -- adiabatic index 0 (piecewise polytropic eos)
!    gamma1      -- adiabatic index 1 (barotropic eos)
!    gamma1pwp   -- adiabatic index 1 (piecewise polytropic eos)
!    gamma2      -- adiabatic index 2 (barotropic eos)
!    gamma2pwp   -- adiabatic index 2 (piecewise polytropic eos)
!    gamma3      -- adiabatic index 3 (barotropic eos)
!    gamma3pwp   -- adiabatic index 3 (piecewise polytropic eos)
!    ieos        -- eqn of state (1=isoth;2=adiab;3=locally iso;8=barotropic)
!    metallicity -- metallicity
!    mu          -- mean molecular weight
!    p1pwp       -- pressure at cutoff density rhocrit1pwp (piecewise polytropic eos)
!    rhocrit0    -- critical density 0 in g/cm^3 (barotropic eos)
!    rhocrit0pwp -- critical density 0 in g/cm^3 (piecewise polytropic eos)
!    rhocrit1    -- critical density 1 in g/cm^3 (barotropic eos)
!    rhocrit1pwp -- critical density 1 in g/cm^3 (piecewise polytropic eos)
!    rhocrit2    -- critical density 2 in g/cm^3 (barotropic eos)
!    rhocrit2pwp -- critical density 2 in g/cm^3 (piecewise polytropic eos)
!    rhocrit3    -- critical density 3 in g/cm^3 (barotropic eos)
!
!  DEPENDENCIES: dim, eos_helmholtz, eos_mesa, infile_utils, io, part,
!    physcon, units
!+
!--------------------------------------------------------------------------
module eos
 implicit none
 integer, parameter, public :: maxeos = 16
 real,               public :: polyk, polyk2, gamma
 real,               public :: qfacdisc
 logical, parameter, public :: use_entropy = .false.
 logical,            public :: extract_eos_from_hdr = .false.
 integer,            public :: isink = 0

 data qfacdisc /0.75/

 public  :: equationofstate,setpolyk,eosinfo,utherm,en_from_utherm
 public  :: get_spsound,get_temperature,get_temperature_from_ponrho
 public  :: gamma_pwp
 public  :: init_eos, finish_eos, write_options_eos, read_options_eos
 public  :: print_eos_to_file

 private

 integer, public :: ieos        = 1
 !--Default initial parameters for Barotropic Eos
 real,    public :: drhocrit0   = 0.50
 real,    public :: rhocrit0cgs = 1.e-18
 real,    public :: rhocrit1cgs = 1.e-14
 real,    public :: rhocrit2cgs = 1.e-10
 real,    public :: rhocrit3cgs = 1.e-3
 real,    public :: gamma1      = 1.4
 real,    public :: gamma2      = 1.1
 real,    public :: gamma3      = 5./3.
 !--Default initial parameters for piecewise polytrope Eos
 real,    public :: rhocrit0pwpcgs = 2.62780d12
 real,    public :: rhocrit1pwpcgs = 5.01187d14
 real,    public :: rhocrit2pwpcgs = 1.0d15
 real,    public :: p1pwpcgs       = 2.46604d34
 real,    public :: gamma0pwp      = 5./3.
 real,    public :: gamma1pwp      = 3.166
 real,    public :: gamma2pwp      = 3.573
 real,    public :: gamma3pwp      = 3.281
 !--Mean molecular weight if temperature required
 real,    public :: gmw            = 2.381
 real,    public :: X_in = 0.74, Z_in = 0.02
 !
 real            :: rhocritT,rhocrit0,rhocrit1,rhocrit2,rhocrit3
 real            :: fac2,fac3,log10polyk2,log10rhocritT,rhocritT0slope
 real            :: rhocrit0pwp,rhocrit1pwp,rhocrit2pwp,p0pwp,p1pwp,p2pwp,k0pwp,k1pwp,k2pwp,k3pwp
 real, public    :: temperature_coef
 !
 logical, public :: done_init_eos = .false.

contains

!----------------------------------------------------------------
!+
!  subroutine returns pressure/density as a function of density
!  (and position in the case of the isothermal disc)
!+
!----------------------------------------------------------------
subroutine equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xi,yi,zi,eni,tempi)
 use io,    only:fatal,error
 use part,  only:xyzmh_ptmass
 use units,   only:unit_density,unit_pressure,unit_ergg
 use eos_mesa, only:get_eos_pressure_gamma1_mesa
 use eos_helmholtz, only:eos_helmholtz_pres_sound

 integer, intent(in)  :: eos_type
 real,    intent(in)  :: rhoi,xi,yi,zi
 real,    intent(out) :: ponrhoi,spsoundi
 real,    intent(inout), optional :: eni
 real,    intent(inout), optional :: tempi
 real :: r,omega,bigH,polyk_new,r1,r2
 real :: gammai
 real :: cgsrhoi, cgseni, cgspgas, pgas, gam1

 select case(eos_type)
 case(1)
!
!--isothermal eos
!
    ponrhoi = polyk
    spsoundi = sqrt(ponrhoi)

 case(2)
!
!--adiabatic/polytropic eos
!  (polytropic using polyk if energy not stored, adiabatic if utherm stored)
!
!   check value of gamma
    if (gamma < tiny(gamma)) call fatal('eos','gamma not set for adiabatic eos',var='gamma',val=gamma)

    if (present(eni)) then
       if (eni < 0.) call fatal('eos','utherm < 0',var='u',val=eni)

       if (use_entropy) then
          ponrhoi = eni*rhoi**(gamma-1.)  ! use this if en is entropy
       elseif (gamma > 1.0001) then
          ponrhoi = (gamma-1.)*eni   ! use this if en is thermal energy
       else
          ponrhoi = 2./3.*eni ! en is thermal energy and gamma = 1
       endif
    else
       ponrhoi = polyk*rhoi**(gamma-1.)
    endif
    spsoundi = sqrt(gamma*ponrhoi)

 case(3)
!
!--this is for a locally isothermal disc as in Lodato & Pringle (2007)
!   cs = cs_0*R^(-q) -- polyk is cs^2, so this is (R^2)^(-q)
!
    ponrhoi = polyk*(xi**2 + yi**2 + zi**2)**(-qfacdisc)
    spsoundi = sqrt(ponrhoi)

 case(6)
!
!--this is for a locally isothermal disc as in Lodato & Pringle (2007), centered on a sink particle
!   cs = cs_0*R^(-q) -- polyk is cs^2, so this is (R^2)^(-q)
    ponrhoi = polyk*((xi-xyzmh_ptmass(1,isink))**2 + (yi-xyzmh_ptmass(2,isink))**2 + &
                     (zi-xyzmh_ptmass(3,isink))**2)**(-qfacdisc)
    spsoundi = sqrt(ponrhoi)

 case(7)
!
!-- z-dependent locally isothermal eos
!
    r = sqrt(xi**2 + yi**2 + zi**2)
    omega = r**(-1.5)
    bigH = (sqrt(polyk)*r**(-qfacdisc))/omega

    if (abs(zi) <= bigH) then
       polyk_new = polyk
    elseif (abs(zi) <= 2.0*bigH) then
       polyk_new = polyk*(1.0+99.0*(abs(zi)-bigH))
    else
       polyk_new = 100.0*polyk
    endif

    ponrhoi = polyk_new*(xi**2 + yi**2 + zi**2)**(-qfacdisc)
    spsoundi = sqrt(ponrhoi)

 case(8)
!
!--Barotropic equation of state
!
    ! variables calculated in the eos initialisation routine:
    !    fac2 = polyk*(rhocrit2/rhocrit1)**(gamma1-1.)
    !    fac3 =  fac2*(rhocrit3/rhocrit2)**(gamma2-1.)
    !    rhocritT0slope = (log10(polyk)-log10(polyk2)) &
    !                   /(log10(rhocritT)-log10(rhocrit0)))
    !
    if (rhoi < rhocritT) then
       gammai  = 1.0
       ponrhoi = polyk2
    elseif (rhoi < rhocrit0) then
       gammai  = 1.0
       ponrhoi = 10**(log10polyk2 + rhocritT0slope*(log10rhocritT-log10(rhoi))  )
    elseif (rhoi < rhocrit1) then
       gammai  = 1.0
       ponrhoi = polyk
    elseif (rhoi < rhocrit2) then
       gammai  = gamma1
       ponrhoi = polyk*(rhoi/rhocrit1)**(gamma1-1.)
    elseif (rhoi < rhocrit3) then
       gammai  = gamma2
       ponrhoi = fac2*(rhoi/rhocrit2)**(gamma2-1.)
    else
       gammai  = gamma3
       ponrhoi = fac3*(rhoi/rhocrit3)**(gamma3-1.)
    endif
    spsoundi = sqrt(gammai*ponrhoi)

 case(9)
!
!--Piecewise Polytropic equation of state
!
    if (rhoi < rhocrit0pwp) then
       gammai  = gamma0pwp
       ponrhoi = k0pwp*rhoi**(gamma0pwp-1.)
    else if (rhoi < rhocrit1pwp) then
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

 case(10)
!
!--MESA eos
!
    cgsrhoi = rhoi * unit_density
    cgseni = eni * unit_ergg

    call get_eos_pressure_gamma1_mesa(cgsrhoi,cgseni,cgspgas,gam1)
    pgas = cgspgas / unit_pressure

    ponrhoi = pgas / rhoi
    spsoundi = sqrt(gam1*ponrhoi)

 case(11)
!
!--isothermal eos with zero pressure
!
    ponrhoi = 0.
    spsoundi = sqrt(polyk)

 case(14)
!
!--locally isothermal prescription from Farris et al. (2014) for binary system
!
    r1=sqrt((xi-xyzmh_ptmass(1,1))**2+(yi-xyzmh_ptmass(2,1))**2 + (zi-xyzmh_ptmass(3,1))**2)
    r2=sqrt((xi-xyzmh_ptmass(1,2))**2+(yi-xyzmh_ptmass(2,2))**2 + (zi-xyzmh_ptmass(3,2))**2 )
    ponrhoi=polyk*((xyzmh_ptmass(4,1)/r1+xyzmh_ptmass(4,2)/r2)/(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2)))**(2*qfacdisc)
    spsoundi=sqrt(ponrhoi)

 case(15)
!
!--helmholtz free energy eos
!
    if (present(tempi)) then
       call eos_helmholtz_pres_sound(tempi, rhoi, ponrhoi, spsoundi, eni)
    else
       call fatal('eos','tried to call Helmholtz free energy eos without passing temperature')
    endif


 case default
    spsoundi = 0. ! avoids compiler warnings
    ponrhoi  = 0.
    call fatal('eos','unknown equation of state')
 end select

 return
end subroutine equationofstate

!----------------------------------------------------------------
!+
!  query function to return the sound speed
!  (called from step for decay timescale in alpha switches)
!+
!----------------------------------------------------------------
real function get_spsound(eos_type,xyzi,rhoi,vxyzui,tempi)
 use dim, only:maxvxyzu
 integer,      intent(in) :: eos_type
 real,         intent(in) :: xyzi(:),rhoi
 real,         intent(inout) :: vxyzui(maxvxyzu)
 real, intent(inout), optional :: tempi
 real :: spsoundi,ponrhoi

 if (maxvxyzu==4) then
    if (present(tempi)) then
       call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),vxyzui(4),tempi)
    else
       call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),vxyzui(4))
    endif
 else
    call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3))
 endif
 get_spsound = spsoundi

end function get_spsound

!-----------------------------------------------------------------------
!+
!  query function to return the temperature
!  (currently only required for non-ideal MHD)
!+
!-----------------------------------------------------------------------
real function get_temperature(eos_type,xyzi,rhoi,vxyzui)
 use dim, only:maxvxyzu
 integer,      intent(in) :: eos_type
 real,         intent(in) :: xyzi(:),rhoi
 real,         intent(inout) :: vxyzui(maxvxyzu)
 real                     :: spsoundi,ponrhoi
 !
 if (maxvxyzu==4) then
    call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),vxyzui(4))
 else
    call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3))
 endif
 get_temperature = temperature_coef*gmw*ponrhoi
 !
end function get_temperature
!-----------------------------------------------------------------------
real function get_temperature_from_ponrho(ponrho)
 real,         intent(in) :: ponrho
 !
 get_temperature_from_ponrho = temperature_coef*gmw*ponrho
 !
end function get_temperature_from_ponrho
!-----------------------------------------------------------------------
!+
!  Get gamma for thermal energy calculations when using the
!  piecewise polytrope
!+
!-----------------------------------------------------------------------
real function gamma_pwp(rhoi)
 real, intent(in) :: rhoi
 !
 if (rhoi < rhocrit0pwp) then
    gamma_pwp = gamma0pwp
 else if (rhoi < rhocrit1pwp) then
    gamma_pwp = gamma1pwp
 elseif (rhoi < rhocrit2pwp) then
    gamma_pwp = gamma2pwp
 else
    gamma_pwp = gamma3pwp
 endif
 !
end function gamma_pwp
!-----------------------------------------------------------------------
!+
!  initialise equation of state (read tables etc.)
!+
!-----------------------------------------------------------------------
subroutine init_eos(eos_type,ierr)
 use units,    only:unit_density,unit_velocity,unit_pressure
 use physcon,  only:mass_proton_cgs,kboltz
 use io,       only:error
 use eos_mesa, only:init_eos_mesa
 use eos_helmholtz, only:eos_helmholtz_init

 integer, intent(in)  :: eos_type
 integer, intent(out) :: ierr
 real                 :: logrhomin,logrhomax
 !
 ierr = 0
 logrhomin = -22.  ! for printing the EoS to file [cgs]; value is for ieos=8
 logrhomax =  -8.  ! for printing the EoS to file [cgs]; value is for ieos=8
 !
 !--Set coefficient to convert P/rho into temperature
 !  calculation will be in cgs; the mean molecular weight, gmw, will be
 !  included in the function call rather than here
 !  c_s^2 = gamma*P/rho = gamma*kT/(gmw*m_p) -> T = P/rho * (gmw*m_p)/k
 !
 temperature_coef = mass_proton_cgs/kboltz * unit_velocity**2
 !
 select case(eos_type)
 case(6)
    !
    !--Check that if using ieos=6, then isink is set properly
    !
    if (isink==0) then
       call error('eos','ieos=6, but isink is not set')
       ierr = 1
       return
    endif

 case(8)
    !
    !--calculate initial variables for the barotropic equation of state
    !
    if (unit_density <= 0.) then
       ierr = 1
       return
    endif

    ! Convert to code units, and calculate constants
    rhocrit0 = rhocrit0cgs/unit_density
    rhocrit1 = rhocrit1cgs/unit_density
    rhocrit2 = rhocrit2cgs/unit_density
    rhocrit3 = rhocrit3cgs/unit_density
    fac2     = polyk*(rhocrit2/rhocrit1)**(gamma1-1.)
    fac3     =  fac2*(rhocrit3/rhocrit2)**(gamma2-1.)

    ! verify that the rhocrit's are in the correct order
    call verify_less_than(ierr,rhocrit0,rhocrit1)
    call verify_less_than(ierr,rhocrit1,rhocrit2)
    call verify_less_than(ierr,rhocrit2,rhocrit3)
    ! Calculate values for the first transition region (no transition if drhocrit0=0)
    if (polyk < tiny(polyk) .or. polyk2 < tiny(polyk2)) drhocrit0 = 0.0

    if (drhocrit0 > 0.0) then
       rhocritT       = rhocrit0*(1.0-drhocrit0)
       log10polyk2    = log10(polyk2)
       log10rhocritT  = log10(rhocritT)
       rhocritT0slope = (log10(polyk)-log10(polyk2)) /(log10(rhocritT)-log10(rhocrit0))
    else
       rhocritT       = rhocrit0  ! moving the transition boundary to rhocrit0
       rhocrit0       = 0.0       ! removing the valid threshhold to enter the transition region
       log10polyk2    = 0.0
       log10rhocritT  = 0.0
       rhocritT0slope = 0.0
    endif

 case(9)
    !
    !--calculate initial variables for the piecewise polytrope equation of state
    !
    if (unit_density <= 0.0 .or. unit_pressure<=0.0) then
       ierr = 1
       return
    endif
    rhocrit0pwp = rhocrit0pwpcgs/unit_density
    rhocrit1pwp = rhocrit1pwpcgs/unit_density
    rhocrit2pwp = rhocrit2pwpcgs/unit_density
    p1pwp       = p1pwpcgs/unit_pressure
    k1pwp       = p1pwp/rhocrit1pwp**gamma1pwp
    k2pwp       = p1pwp/rhocrit1pwp**gamma2pwp
    p2pwp       = k2pwp*rhocrit2pwp**gamma2pwp
    k3pwp       = p2pwp/rhocrit2pwp**gamma3pwp
    k0pwp       = k1pwp/(rhocrit0pwp**(gamma0pwp-gamma1pwp))
    p0pwp       = k0pwp*rhocrit0pwp**gamma0pwp
    !
    ! for testing the EoS
    logrhomin = 10.  ! for testing the EoS [cgs]
    logrhomax = 20.  ! for testing the EoS [cgs]

 case(10)
    !
    !--MESA EoS initialisation
    !
    call init_eos_mesa(X_in,Z_in,ierr)

 case(15)

    call eos_helmholtz_init(ierr)

 end select
 done_init_eos = .true.

end subroutine init_eos

!-----------------------------------------------------------------------
!+
!  finish equation of state
!+
!-----------------------------------------------------------------------

subroutine finish_eos(eos_type,ierr)
 use eos_mesa, only: finish_eos_mesa

 integer, intent(in)  :: eos_type
 integer, intent(out) :: ierr

 ierr = 0

 select case(eos_type)

 case(10)
    !
    !--MESA EoS deallocation
    !
    call finish_eos_mesa

 end select
 done_init_eos=.false.
end subroutine finish_eos

!-----------------------------------------------------------------------
!+
!  verify that val1 < val2
!+
!-----------------------------------------------------------------------
subroutine verify_less_than(ierr,val1,val2)
 use io, only: error
 integer, intent(inout) :: ierr
 real,    intent(in)    :: val1,val2
 !
 if (val1 > val2) then
    ierr = ierr + 1
    call error('eos','incorrect ordering of rhocrit')
 endif
 !
end subroutine verify_less_than
!-----------------------------------------------------------------------
!+
!  allow the user to print the eos to file
!+
!-----------------------------------------------------------------------
subroutine print_eos_to_file(logrhomin,logrhomax,unit_density,unit_velocity)
 use io,               only: iuniteos
 real,         intent(in) :: logrhomin,logrhomax
 real(kind=8), intent(in) :: unit_density,unit_velocity
 integer,      parameter  :: nlogrho   = 1000
 real                     :: rho,drho,ponrhoi,spsoundi,dummy,temperaturei
 integer                  :: i
 !
 !--Open file
 open(unit=iuniteos,file="EOS.dat",form='formatted',status='replace')
 write(iuniteos,'("# Equation of state properties; all values in cgs")')
 write(iuniteos,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'rho', &
       2,'P', &
       3,'P/rho', &
       4,'c_s', &
       5,'T'
 !
 dummy = 0.0  ! initialise to avoind compiler warning
 drho  = (logrhomax - logrhomin)/float(nlogrho)
 do i = 1,nlogrho
    rho = 10**(logrhomin +(i-1)*drho)/unit_density
    call equationofstate(ieos,ponrhoi,spsoundi,rho,dummy,dummy,dummy)
    temperaturei = get_temperature_from_ponrho(ponrhoi)
    write(iuniteos,'(5(1pe18.10,1x))') &
       rho*unit_density, ponrhoi*unit_velocity**2*rho*unit_density, &
       ponrhoi*unit_velocity**2, spsoundi*unit_velocity**2, temperaturei
 enddo
 close(iuniteos)
 !
end subroutine print_eos_to_file
!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos(iunit)
 use infile_utils, only:write_inopt
 use eos_helmholtz, only:eos_helmholtz_write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling equation of state'
 call write_inopt(ieos,'ieos','eqn of state (1=isoth;2=adiab;3=locally iso;8=barotropic)',iunit)
 call write_inopt(gmw,'mu','mean molecular weight',iunit)
 select case(ieos)
 case(8)
    call write_inopt(drhocrit0,  'drhocrit','transition size between rhocrit0 & 1 (fraction of rhocrit0; barotropic eos)',iunit)
    call write_inopt(rhocrit0cgs,'rhocrit0','critical density 0 in g/cm^3 (barotropic eos)',iunit)
    call write_inopt(rhocrit1cgs,'rhocrit1','critical density 1 in g/cm^3 (barotropic eos)',iunit)
    call write_inopt(rhocrit2cgs,'rhocrit2','critical density 2 in g/cm^3 (barotropic eos)',iunit)
    call write_inopt(rhocrit3cgs,'rhocrit3','critical density 3 in g/cm^3 (barotropic eos)',iunit,exp=.true.)
    call write_inopt(gamma1,'gamma1','adiabatic index 1 (barotropic eos)',iunit)
    call write_inopt(gamma2,'gamma2','adiabatic index 2 (barotropic eos)',iunit)
    call write_inopt(gamma3,'gamma3','adiabatic index 3 (barotropic eos)',iunit)

 case(9)
    call write_inopt(rhocrit0pwpcgs,'rhocrit0pwp','critical density 0 in g/cm^3 (piecewise polytropic eos)',iunit)
    call write_inopt(rhocrit1pwpcgs,'rhocrit1pwp','critical density 1 in g/cm^3 (piecewise polytropic eos)',iunit)
    call write_inopt(rhocrit2pwpcgs,'rhocrit2pwp','critical density 2 in g/cm^3 (piecewise polytropic eos)',iunit,exp=.true.)
    call write_inopt(gamma0pwp,'gamma0pwp','adiabatic index 0 (piecewise polytropic eos)',iunit)
    call write_inopt(gamma1pwp,'gamma1pwp','adiabatic index 1 (piecewise polytropic eos)',iunit)
    call write_inopt(gamma2pwp,'gamma2pwp','adiabatic index 2 (piecewise polytropic eos)',iunit)
    call write_inopt(gamma3pwp,'gamma3pwp','adiabatic index 3 (piecewise polytropic eos)',iunit)
    call write_inopt(p1pwpcgs,'p1pwp','pressure at cutoff density rhocrit1pwp (piecewise polytropic eos)',iunit)

 case(10)
    call write_inopt(X_in,'X','hydrogen mass fraction',iunit)
    call write_inopt(Z_in,'Z','metallicity',iunit)

 case(15) ! helmholtz eos
    call eos_helmholtz_write_inopt(iunit)

 end select

end subroutine write_options_eos

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_eos(name,valstring,imatch,igotall,ierr)
 use io,            only:fatal
 use eos_helmholtz, only:eos_helmholtz_set_relaxflag
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot  = 0
 character(len=30), parameter  :: label = 'read_options_eos'
 integer :: tmp

 imatch  = .true.
 select case(trim(name))
 case('ieos')
    read(valstring,*,iostat=ierr) ieos
    ngot = ngot + 1
    if (ieos <= 0 .or. ieos > maxeos) call fatal(label,'equation of state choice out of range')
 case('mu')
    read(valstring,*,iostat=ierr) gmw
    ! not compulsory to read in
    if (gmw <= 0.)  call fatal(label,'mu <= 0')

 case('drhocrit')
    read(valstring,*,iostat=ierr) drhocrit0
    if (drhocrit0 < 0.)  call fatal(label,'drhocrit0 < 0: Negative transition region is nonsense')
    if (drhocrit0 > 1.)  call fatal(label,'drhocrit0 > 1: Too large of transition region')
    ngot = ngot + 1
 case('rhocrit0')
    read(valstring,*,iostat=ierr) rhocrit0cgs
    if (rhocrit0cgs <= 0.) call fatal(label,'rhocrit0 <= 0')
    ngot = ngot + 1
 case('rhocrit1')
    read(valstring,*,iostat=ierr) rhocrit1cgs
    if (rhocrit1cgs <= 0.) call fatal(label,'rhocrit1 <= 0')
    ngot = ngot + 1
 case('rhocrit2')
    read(valstring,*,iostat=ierr) rhocrit2cgs
    if (rhocrit2cgs <= 0.) call fatal(label,'rhocrit2 <= 0')
    ngot = ngot + 1
 case('rhocrit3')
    read(valstring,*,iostat=ierr) rhocrit3cgs
    if (rhocrit3cgs <= 0.) call fatal(label,'rhocrit3 <= 0')
    ngot = ngot + 1
 case('gamma1')
    read(valstring,*,iostat=ierr) gamma1
    if (gamma1 < 1.) call fatal(label,'gamma1 < 1.0')
    ngot = ngot + 1
 case('gamma2')
    read(valstring,*,iostat=ierr) gamma2
    if (gamma2 < 1.) call fatal(label,'gamma2 < 1.0')
    ngot = ngot + 1
 case('gamma3')
    read(valstring,*,iostat=ierr) gamma3
    if (gamma3 < 1.) call fatal(label,'gamma3 < 1.0')
    ngot = ngot + 1
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
 case('X')
    read(valstring,*,iostat=ierr) X_in
    if (X_in <= 0.) call fatal(label,'X <= 0.0')
    ngot = ngot + 1
 case('Z')
    read(valstring,*,iostat=ierr) Z_in
    if (Z_in <= 0.) call fatal(label,'Z <= 0.0')
    ngot = ngot + 1
 case('relaxflag')
    ! ideally would like this to be self-contained within eos_helmholtz,
    ! but it's a bit of a pain and this is easy
    read(valstring,*,iostat=ierr) tmp
    call eos_helmholtz_set_relaxflag(tmp)
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 !--make sure we have got all compulsory options (otherwise, rewrite input file)
 if (ieos==8) then
    igotall = (ngot >= 9)
 elseif (ieos==9) then
    igotall = (ngot >= 9)
 else
    igotall = (ngot >= 1)
 endif

end subroutine read_options_eos

!----------------------------------------------------------------
!+
!  subroutine sets polyk based on utherm/positions
!  read from an sphNG dump file
!+
!----------------------------------------------------------------
subroutine setpolyk(eos_type,iprint,utherm,xyzhi,npart)
 use part, only:xyzmh_ptmass
 integer, intent(in) :: eos_type,iprint
 real,    intent(in) :: utherm(:)
 real,    intent(in) :: xyzhi(:,:)
 integer, intent(in) :: npart
 integer :: ipart
 real :: r2,polykalt

 !-- pick a random particle from which to extract polyk
 ipart = npart/2

 select case(eos_type)
 case(1,8)
!
!--isothermal eos
!
    polykalt = 2./3.*utherm(ipart)
    !--check all other utherms identical
    if (any(utherm(1:npart) /= utherm(ipart))) then
       write(iprint,*) 'WARNING! different utherms but run is isothermal'
    endif

 case(2)
!
!--adiabatic/polytropic eos
!  this routine is ONLY called if utherm is NOT stored, so polyk matters
!
    write(iprint,*) 'Using polytropic equation of state, gamma = ',gamma
    polykalt = 2./3.*utherm(ipart)
    if (gamma <= 1.00000001) then
       stop 'silly to use gamma==1 without using isothermal eos'
    endif

 case(3)
!
!--locally isothermal disc as in Lodato & Pringle (2007)
!   cs = cs_0*R^(-q) -- polyk is cs^2, so this is (R^2)^(-q)
!
    r2 = xyzhi(1,ipart)*xyzhi(1,ipart) + xyzhi(2,ipart)*xyzhi(2,ipart) &
       + xyzhi(3,ipart)*xyzhi(3,ipart)
    polykalt = 2./3.*utherm(ipart)*r2**qfacdisc

 case(6)
!
!--locally isothermal disc as in Lodato & Pringle (2007), centered on specified sink particle
!   cs = cs_0*R^(-q) -- polyk is cs^2, so this is (R^2)^(-q)
!
    r2 = (xyzhi(1,ipart)-xyzmh_ptmass(1,isink))**2 + &
         (xyzhi(2,ipart)-xyzmh_ptmass(2,isink))**2 + &
         (xyzhi(3,ipart)-xyzmh_ptmass(3,isink))**2

    polykalt = 2./3.*utherm(ipart)*r2**qfacdisc
 case default
!
!--don't die in this routine as it can be called from readdump
!  (ie. not necessarily as part of a run)
!
    write(iprint,*) ' WARNING! unknown equation of state in setpolyk'
    polykalt = polyk

 end select

 if (diff(polykalt,polyk)) then
    write(iprint,*) 'WARNING! polyk set using RK2 in dump differs from that set using thermal energy'
    write(iprint,*) 'using polyk = ',polykalt, ' (from RK2 = ',polyk,')'
 endif
 polyk = polykalt
!
!--warn if polyk is zero, die if negative
!
 if (polyk < 0.) then
    write(iprint,*) 'ERROR: polyk < 0 in setting equation of state'
    stop
 elseif (polyk < tiny(polyk)) then
    write(iprint,*) 'WARNING: polyk = 0 in equation of state'
 endif

 return
end subroutine setpolyk
!----------------------------------------------------------------
!+
!  small utility returns whether two real numbers differ
!+
!----------------------------------------------------------------
logical pure function diff(r1,r2)
 real, intent(in) :: r1,r2

 diff = abs(r1-r2) > tiny(r1)

end function diff

!----------------------------------------------------------------
!+
!  prints equation of state info in the run header
!+
!----------------------------------------------------------------

subroutine eosinfo(eos_type,iprint)
 use dim,           only:maxvxyzu
 use io,            only:fatal
 use units,         only:unit_density, unit_velocity
 use eos_helmholtz, only:eos_helmholtz_eosinfo
 integer, intent(in) :: eos_type,iprint
 real, parameter :: uthermcheck = 3.14159, rhocheck = 23.456

 select case(eos_type)
    !
 case(1,11)
    write(iprint,"(/,a,f10.6)") ' Isothermal equation of state:     cs^2 = ',polyk
    if (eos_type==11) write(iprint,*) ' (ZERO PRESSURE) '
    !
 case(2)
    if (use_entropy) then
       write(iprint,"(/,a,f10.6,a,f10.6)") ' Adiabatic equation of state (evolving ENTROPY): polyk = ',polyk,' gamma = ',gamma
!
!--run a unit test on the en-> utherm and utherm-> en conversion utilities
!
       write(iprint,"(a)",ADVANCE='NO') ' checking utherm -> entropy -> utherm conversion ...'
       if (abs(utherm(en_from_utherm(uthermcheck,rhocheck),rhocheck)-uthermcheck) > epsilon(uthermcheck)) then
          call fatal('eosinfo','failed consistency check in eos: utherm0 -> entropy -> utherm  /=  utherm0')
       elseif (abs(en_from_utherm(utherm(uthermcheck,rhocheck),rhocheck)-uthermcheck) > epsilon(uthermcheck)) then
          call fatal('eosinfo','failed consistency check in eos: entropy0 -> utherm -> entropy  /=  entropy0')
       else
          write(iprint,*) 'OK'
       endif
    elseif (maxvxyzu >= 4) then
       write(iprint,"(/,a,f10.6)") ' Adiabatic equation of state (evolving UTHERM): P = (gamma-1)*rho*u, gamma = ',gamma
    else
       write(iprint,"(/,a,f10.6,a,f10.6)") ' Polytropic equation of state: P = ',polyk,'*rho^',gamma
    endif
    !
 case(3)
    write(iprint,"(/,a,f10.6,a,f10.6)") ' Locally isothermal eq of state (R_sph): cs^2_0 = ',polyk,' qfac = ',qfacdisc
    !
 case(6)
    write(iprint,"(/,a,i2,a,f10.6,a,f10.6)") ' Locally (on sink ',isink, &
          ') isothermal eos (R_sph): cs^2_0 = ',polyk,' qfac = ',qfacdisc
    !
 case(8)
    write(iprint,"(/,a,2(es10.3,a))")    ' Barotropic eq of state: cs_ld            = ',sqrt(polyk2),' code units = '&
                                         ,sqrt(polyk2)*unit_velocity,' cm/s'
    write(iprint,"(  a,2(es10.3,a))")    ' Barotropic eq of state: cs               = ',sqrt(polyk), ' code units = '&
                                         ,sqrt(polyk)*unit_velocity, ' cm/s'
    if (drhocrit0 > 0.0) then
       write(iprint,"(  a,2(es10.3,a))") ' Barotropic eq of state: rhocritT == rhoT = ',rhocritT,    ' code units = '&
                                         ,rhocritT*unit_density,     ' g/cm^3'
       write(iprint,"(  a,2(es10.3,a))") ' Barotropic eq of state: rhocrit0 == rho0 = ',rhocrit0,    ' code units = '&
                                         ,rhocrit0*unit_density,     ' g/cm^3'
    else
       write(iprint,"(  a,2(es10.3,a))") ' Barotropic eq of state: rhocrit0 == rho0 = ',rhocritT,    ' code units = '&
                                         ,rhocritT*unit_density,     ' g/cm^3'
    endif

    write(iprint,"(  a,2(es10.3,a))")    ' Barotropic eq of state: rhocrit1 == rho1 = ',rhocrit1,    ' code units = '&
                                         ,rhocrit1*unit_density,     ' g/cm^3'
    write(iprint,"(  a,2(es10.3,a))")    ' Barotropic eq of state: rhocrit2 == rho2 = ',rhocrit2,    ' code units = '&
                                         ,rhocrit2*unit_density,     ' g/cm^3'
    write(iprint,"(  a,2(es10.3,a))")    ' Barotropic eq of state: rhocrit3 == rho3 = ',rhocrit3,    ' code units = '&
                                         ,rhocrit3*unit_density,     ' g/cm^3'
    write(iprint,"(a)")                  ' Barotropic eq of state:'
    if (drhocrit0 > 0.0) then
       write(iprint,"(a,56x,a)")         ' Barotropic eq of state: P = cs_ld*rho','for         rho/(g/cm^3) < rhoT'
       write(iprint,"(a,14x,a)")         ' Barotropic eq of state: P = 10**(log10(cs_bg**2) + M*(log10(rhoT)-log10(rho)))' &
                                         ,' for rhoT <= rho/(g/cm^3) < rho0'
    else
       write(iprint,"(a,56x,a)")         ' Barotropic eq of state: P = cs_ld*rho','for         rho/(g/cm^3) < rho0'
    endif

    write(iprint,"(a,59x,a)")            ' Barotropic eq of state: P = cs*rho','for rho0 <= rho/(g/cm^3) < rho1'
    write(iprint,"(a,f6.3,39x,a)")       ' Barotropic eq of state: P = cs*rho1*(rho /rho1)^',gamma1 &
                                         ,'for rho1 <= rho/(g/cm^3) < rho2'
    write(iprint,"(2(a,f6.3),a)")        ' Barotropic eq of state: P = cs*rho1*(rho2/rho1)^',gamma1 &
                                         ,'*(rho /rho2)^',gamma2,'                    for rho2 <= rho/(g/cm^3) < rho3'
    write(iprint,"(3(a,f6.3),a)")        ' Barotropic eq of state: P = cs*rho1*(rho2/rho1)^',gamma1 &
                                        ,'*(rho3/rho2)^',gamma2,'*(rho /rho3)^',gamma3,' for rho3 <= rho/(g/cm^3)'
    !
 case(9)
    write(iprint,"(/,a,3(es10.3),a,4(es10.3))") ' Piecewise polytropic eq of state (code units) : rhocrit = '&
                                                 ,rhocrit0pwp,rhocrit1pwp,rhocrit2pwp, '; K = ',k0pwp,k1pwp,k2pwp,k3pwp
    write(iprint,"(  a,3(es10.3)            )") ' Piecewise polytropic eq of state (g/cm^3)     : rhocrit = '&
                                                 ,rhocrit0pwp*unit_density,rhocrit1pwp*unit_density,rhocrit2pwp*unit_density
 case(15)
    call eos_helmholtz_eosinfo(iprint)

 end select
 write(iprint,*)

 return
end subroutine eosinfo

!----------------------------------------------------------------
!+
!  the following two functions transparently handle evolution
!  of the entropy instead of the thermal energy
!+
!----------------------------------------------------------------
real function utherm(en,rho)
 real, intent(in) :: en, rho
 real :: gamm1

 if (use_entropy) then
    gamm1 = gamma - 1.
    if (gamm1 > tiny(gamm1)) then
       utherm = (en/gamm1)*rho**gamm1
    else
       stop 'gamma=1 using entropy evolution'
    endif
 else
    utherm = en
 endif

 return
end function utherm

!----------------------------------------------------------------
!+
!  function to transparently handle evolution of the entropy
!  instead of the thermal energy
!+
!----------------------------------------------------------------
real function en_from_utherm(utherm,rho)
 real, intent(in) :: utherm, rho
 real :: gamm1

 if (use_entropy) then
    gamm1 = gamma - 1.
    if (gamm1 > tiny(gamm1)) then
       en_from_utherm = gamm1*utherm*rho**(1.-gamma)
    else
       stop 'gamma=1 using entropy evolution'
    endif
 else
    en_from_utherm = utherm
 endif

 return
end function en_from_utherm

end module eos
