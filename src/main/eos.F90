!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module eos
!
! This module contains stuff to do with the equation of state
!  Current options:
!     1 = isothermal eos
!     2 = adiabatic/polytropic eos
!     3 = eos for a locally isothermal disc as in Lodato & Pringle (2007)
!     4 = GR isothermal
!     6 = eos for a locally isothermal disc as in Lodato & Pringle (2007),
!         centered on a sink particle
!     7 = z-dependent locally isothermal eos
!     8 = Barotropic eos
!     9 = Piecewise polytrope
!    10 = MESA EoS
!    11 = isothermal eos with zero pressure
!    12 = ideal gas with radiation pressure
!    14 = locally isothermal prescription from Farris et al. (2014) for binary system
!    15 = Helmholtz free energy eos
!    16 = Shen eos
!    20 = Ideal gas + radiation + various forms of recombination energy from HORMONE (Hirai et al., 2020)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - X           : *H mass fraction (ignored if variable composition)*
!   - Z           : *metallicity (ignored if variable composition)*
!   - ieos        : *eqn of state (1=isoth;2=adiab;3=locally iso;8=barotropic)*
!   - metallicity : *metallicity*
!   - mu          : *mean molecular weight*
!
! :Dependencies: dim, eos_barotropic, eos_gasradrec, eos_helmholtz,
!   eos_idealplusrad, eos_mesa, eos_piecewise, eos_shen, infile_utils, io,
!   mesa_microphysics, options, part, physcon, units
!
 implicit none
 integer, parameter, public :: maxeos = 20
 real,               public :: polyk, polyk2, gamma
 real,               public :: qfacdisc = 0.75, qfacdisc2 = 0.75
 logical, parameter, public :: use_entropy = .false.
 logical,            public :: extract_eos_from_hdr = .false.
 integer,            public :: isink = 0.


 public  :: equationofstate,setpolyk,eosinfo,utherm,en_from_utherm,get_mean_molecular_weight
 public  :: get_spsound,get_temperature,get_temperature_from_ponrho,eos_is_non_ideal,eos_outputs_mu
#ifdef KROME
 public  :: get_local_u_internal
#endif
 public  :: calc_rec_ene,calc_temp_and_ene,entropy,get_rho_from_p_s
 public  :: init_eos,finish_eos,write_options_eos,read_options_eos
 public  :: print_eos_to_file
 public  :: write_headeropts_eos, read_headeropts_eos

 private

 integer, public :: ieos          = 1
 integer, public :: iopacity_type = 0 ! used for radiation
 logical, public :: use_var_comp = .false. ! use variable composition
 !--Mean molecular weight if temperature required
 real,    public :: gmw           = 2.381
 real,    public :: X_in = 0.74, Z_in = 0.02
 !--Minimum temperature (failsafe to prevent u < 0)
 real,    public :: Tfloor = 0. ![K]
 real,    public :: temperature_coef

 logical, public :: done_init_eos = .false.
 !
 ! error codes for calls to init_eos
 !
 integer, public, parameter :: &
    ierr_file_not_found  = 1, &
    ierr_option_conflict = 2, &
    ierr_units_not_set   = 3, &
    ierr_isink_not_set   = 4

!
! 2D temperature structure fit parameters for HD 163296
!
real, public :: z0      = 1.
real, public :: alpha_z = 3.01
real, public :: beta_z  = 0.42


contains

!----------------------------------------------------------------
!+
!  subroutine returns pressure/density as a function of density
!  (and position in the case of the isothermal disc)
!+
!----------------------------------------------------------------
subroutine equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xi,yi,zi,eni,tempi,gamma_local,mu_local,Xlocal,Zlocal)
 use io,            only:fatal,error,warning
 use part,          only:xyzmh_ptmass
 use units,         only:unit_density,unit_pressure,unit_ergg,unit_velocity
 use physcon,       only:kb_on_mh
 use eos_mesa,      only:get_eos_pressure_temp_gamma1_mesa
 use eos_helmholtz, only:eos_helmholtz_pres_sound
 use eos_shen,      only:eos_shen_NL3
 use eos_idealplusrad
 use eos_gasradrec, only:equationofstate_gasradrec
 use eos_barotropic, only:get_eos_barotropic
 use eos_piecewise,  only:get_eos_piecewise
 integer, intent(in)  :: eos_type
 real,    intent(in)  :: rhoi,xi,yi,zi
 real,    intent(out) :: ponrhoi,spsoundi
 real,    intent(inout), optional :: eni
 real,    intent(inout), optional :: tempi,mu_local
 real,    intent(in)   , optional :: gamma_local,Xlocal,Zlocal
 real :: r1,r2
 real :: gammai,temperaturei,mui,imui,X_i,Z_i
 real :: cgsrhoi,cgseni,cgspresi,presi,gam1,cgsspsoundi
 integer :: ierr
 real :: uthermconst
 real :: zq,cs2atm,cs2mid,cs2
#ifdef GR
 real :: enthi,pondensi
! Check to see if adiabatic equation of state is being used.
 if (eos_type /= 2 .and. eos_type /= 4 .and. eos_type /= 11 .and. eos_type /= 12) &
 call fatal('eos','GR is only compatible with an adiabatic equation of state (ieos=2), for the time being.',&
 var='eos_type',val=real(eos_type))
#endif

 gammai = gamma
 mui = gmw
 X_i = X_in
 Z_i = Z_in
 if (present(gamma_local)) gammai = gamma_local
 if (present(mu_local)) mui = mu_local
 if (present(Xlocal)) X_i = Xlocal
 if (present(Zlocal)) Z_i = Zlocal

 select case(eos_type)
 case(1)
!
!--isothermal eos
!
    ponrhoi  = polyk
    spsoundi = sqrt(ponrhoi)
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(2)
!
!--adiabatic/polytropic eos
!  (polytropic using polyk if energy not stored, adiabatic if utherm stored)
!
!   check value of gamma
    if (gammai < tiny(gammai)) call fatal('eos','gamma not set for adiabatic eos',var='gamma',val=gammai)

#ifdef GR
    if (.not. present(eni)) call fatal('eos','GR call to equationofstate requires thermal energy as input!')
    if (eni < 0.) call fatal('eos','utherm < 0',var='u',val=eni)
    if (gammai == 1.) then
       call fatal('eos','GR not compatible with isothermal equation of state, yet...',var='gamma',val=gammai)
    elseif (gammai > 1.0001) then
       pondensi = (gammai-1.)*eni   ! eni is the thermal energy
       enthi = 1. + eni + pondensi    ! enthalpy
       spsoundi = sqrt(gammai*pondensi/enthi)
       ponrhoi = pondensi ! With GR this routine actually outputs pondensi (i.e. pressure on primitive density, not conserved.)
    endif
#else
    if (present(eni)) then
       if (eni < 0.) then
          !write(iprint,'(a,Es18.4,a,4Es18.4)')'Warning: eos: u = ',eni,' < 0 at {x,y,z,rho} = ',xi,yi,zi,rhoi
          call fatal('eos','utherm < 0',var='u',val=eni)
       endif
       if (use_entropy) then
          ponrhoi = eni*rhoi**(gammai-1.)  ! use this if en is entropy
       elseif (gammai > 1.0001) then
          ponrhoi = (gammai-1.)*eni   ! use this if en is thermal energy
       else
          ponrhoi = 2./3.*eni ! en is thermal energy and gamma = 1
       endif
    else
       ponrhoi = polyk*rhoi**(gammai-1.)
    endif
    spsoundi = sqrt(gammai*ponrhoi)
#endif

    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(3)
!
!--this is for a locally isothermal disc as in Lodato & Pringle (2007)
!   cs = cs_0*R^(-q) -- polyk is cs^2, so this is (R^2)^(-q)
!
    ponrhoi  = polyk*(xi**2 + yi**2 + zi**2)**(-qfacdisc)
    spsoundi = sqrt(ponrhoi)
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

!
!--GR isothermal
!
 case(4)
    uthermconst = polyk
    ponrhoi = (gammai-1.)*uthermconst
    spsoundi = sqrt(ponrhoi/(1.+uthermconst))
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(6)
!
!--this is for a locally isothermal disc as in Lodato & Pringle (2007), centered on a sink particle
!   cs = cs_0*R^(-q) -- polyk is cs^2, so this is (R^2)^(-q)
    ponrhoi  = polyk*((xi-xyzmh_ptmass(1,isink))**2 + (yi-xyzmh_ptmass(2,isink))**2 + &
                      (zi-xyzmh_ptmass(3,isink))**2)**(-qfacdisc)
    spsoundi = sqrt(ponrhoi)
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(7)
!
!-- z-dependent locally isothermal eos
    r2 = xi**2 + yi**2
    cs2mid = polyk * r2**(-qfacdisc)
    cs2atm = polyk2 * r2**(-qfacdisc2)
    zq = z0 * r2**(0.5*beta_z)

    ! modified equation 6 from Law et al. (2021)
    cs2 = (cs2mid**4 + 0.5*(1 + tanh((zi - alpha_z*zq)/zq))*cs2atm**4)**(1./4.)
    ponrhoi = cs2
    spsoundi = sqrt(cs2)
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(8)
!
!--Barotropic equation of state
!
    call get_eos_barotropic(rhoi,polyk,polyk2,ponrhoi,spsoundi,gammai)
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(9)
!
!--Piecewise Polytropic equation of state
!
    call get_eos_piecewise(rhoi,ponrhoi,spsoundi,gammai)
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(10)
!
!--MESA eos
!
    cgsrhoi = rhoi * unit_density
    cgseni  = eni * unit_ergg
    call get_eos_pressure_temp_gamma1_mesa(cgsrhoi,cgseni,cgspresi,temperaturei,gam1,ierr)
    presi = cgspresi / unit_pressure

    ponrhoi  = presi / rhoi
    spsoundi = sqrt(gam1*ponrhoi)
    if (present(tempi)) tempi = temperaturei
    if (ierr /= 0) call warning('eos_mesa','extrapolating off tables')

 case(11)
!
!--isothermal eos with zero pressure
!
    ponrhoi  = 0.
    spsoundi = sqrt(polyk)
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(12)
!
!--ideal gas plus radiation pressure
!
    if (present(tempi)) then
       temperaturei = tempi
    else
       temperaturei = -1. ! Use gas temperature as initial guess
    endif
    cgsrhoi = rhoi * unit_density
    cgseni  = eni * unit_ergg
    call get_idealplusrad_temp(cgsrhoi,cgseni,mui,gammai,temperaturei,ierr)
    call get_idealplusrad_pres(cgsrhoi,temperaturei,mui,cgspresi)
    call get_idealplusrad_spsoundi(cgsrhoi,cgspresi,cgseni,spsoundi)
    spsoundi = spsoundi / unit_velocity
    presi = cgspresi / unit_pressure
    ponrhoi = presi / rhoi
    if (present(tempi)) tempi = temperaturei
    if (ierr /= 0) call warning('eos_idealplusrad','temperature iteration did not converge')

 case(14)
!
!--locally isothermal prescription from Farris et al. (2014) for binary system
!
    r1=sqrt((xi-xyzmh_ptmass(1,1))**2+(yi-xyzmh_ptmass(2,1))**2 + (zi-xyzmh_ptmass(3,1))**2)
    r2=sqrt((xi-xyzmh_ptmass(1,2))**2+(yi-xyzmh_ptmass(2,2))**2 + (zi-xyzmh_ptmass(3,2))**2)
!  ponrhoi=polyk*(xyzmh_ptmass(4,1)/r1+xyzmh_ptmass(4,2)/r2)**(2*qfacdisc)/(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2))**(2*qfacdisc)
    ponrhoi=polyk*(xyzmh_ptmass(4,1)/r1+xyzmh_ptmass(4,2)/r2)**(2*qfacdisc)/(xyzmh_ptmass(4,1))**(2*qfacdisc)
    spsoundi=sqrt(ponrhoi)
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi

 case(15)
!
!--helmholtz free energy eos
!
    if (present(tempi)) then
       call eos_helmholtz_pres_sound(tempi, rhoi, ponrhoi, spsoundi, eni)
    else
       ponrhoi  = 0.
       spsoundi = 0.
       call fatal('eos','tried to call Helmholtz free energy eos without passing temperature')
    endif

 case(16)
!
!--shen eos
!
!    if (present(enei)) then
    cgsrhoi = rhoi * unit_density
    !note eni is actually tempi
    call eos_shen_NL3(cgsrhoi,eni,0.05,cgspresi,cgsspsoundi)
    spsoundi=cgsspsoundi / unit_velocity
    presi = cgspresi / unit_pressure
    ponrhoi = presi / rhoi
    if (present(tempi)) tempi = eni
!    else
!       call fatal('eos','tried to call NL3 eos without passing temperature')
!    endif

 case(20)
!
!--gas + radiation + various forms of recombination (from HORMONE, Hirai+20)
!
    cgsrhoi = rhoi * unit_density
    cgseni  = eni * unit_ergg
    imui = 1./mui
    if (present(tempi) .and. (tempi > 0.)) then
       temperaturei = tempi
    else
       temperaturei = 0.67 * cgseni * mui / kb_on_mh
    endif
    call equationofstate_gasradrec(cgsrhoi,cgseni*cgsrhoi,temperaturei,imui,X_i,1.-X_i-Z_i,cgspresi,cgsspsoundi)
    ponrhoi = cgspresi / (unit_pressure * rhoi)
    spsoundi = cgsspsoundi / unit_velocity
    if (present(mu_local)) mu_local = 1./imui
    if (present(tempi)) tempi = temperaturei

 case default
    spsoundi = 0. ! avoids compiler warnings
    ponrhoi  = 0.
    if (present(tempi)) tempi = temperature_coef*mui*ponrhoi
    call fatal('eos','unknown equation of state')
 end select

end subroutine equationofstate

!----------------------------------------------------------------
!+
!  Query function to return whether an EoS is non-ideal
!  Mainly used to decide whether it is necessary to write
!  things like pressure and temperature in the dump file or not
!+
!----------------------------------------------------------------
logical function eos_is_non_ideal(ieos)
 integer, intent(in) :: ieos

 select case(ieos)
 case(10,12,15,20)
    eos_is_non_ideal = .true.
 case default
    eos_is_non_ideal = .false.
 end select

end function eos_is_non_ideal

!----------------------------------------------------------------
!+
!  Query function to return whether an EoS outputs mean molecular weight
!+
!----------------------------------------------------------------
logical function eos_outputs_mu(ieos)
 integer, intent(in) :: ieos

 select case(ieos)
 case(20)
    eos_outputs_mu = .true.
 case default
    eos_outputs_mu = .false.
 end select

end function eos_outputs_mu

!----------------------------------------------------------------
!+
!  query function to return the sound speed
!  (called from step for decay timescale in alpha switches)
!+
!----------------------------------------------------------------
real function get_spsound(eos_type,xyzi,rhoi,vxyzui,tempi,gammai,mui,Xi,Zi)
 use dim, only:maxvxyzu
 integer,      intent(in)      :: eos_type
 real,         intent(in)      :: xyzi(:),rhoi
 real,         intent(inout)   :: vxyzui(:)
 real, intent(inout)   , optional    :: tempi
 real, intent(in)      , optional    :: gammai,mui,Xi,Zi
 real :: spsoundi,ponrhoi,mu,X,Z

 mu = gmw
 X = X_in
 Z = Z_in
 if (present(mui)) mu = mui
 if (present(Xi)) X = Xi
 if (present(Zi)) Z = Zi

 if (maxvxyzu==4) then
    if (present(gammai)) then
       call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),vxyzui(4),&
                            gamma_local=gammai,mu_local=mu,Xlocal=X,Zlocal=Z)
    elseif (present(tempi)) then
       call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),vxyzui(4),&
                            tempi=tempi,mu_local=mu,Xlocal=X,Zlocal=Z)
    else
       call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),vxyzui(4),mu_local=mu,Xlocal=X,Zlocal=Z)
    endif
 elseif (present(tempi)) then
    call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),tempi=tempi,mu_local=mu)
 else
    call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),mu_local=mu,Xlocal=X,Zlocal=Z)
 endif
 get_spsound = spsoundi

end function get_spsound

!-----------------------------------------------------------------------
!+
!  query function to return the temperature given density,
!  position and/or thermal energy
!+
!-----------------------------------------------------------------------
real function get_temperature(eos_type,xyzi,rhoi,vxyzui,gammai,mui,Xi,Zi) result(tempi)
 use dim, only:maxvxyzu
 integer,      intent(in)    :: eos_type
 real,         intent(in)    :: xyzi(:),rhoi
 real,         intent(inout) :: vxyzui(:)
 real, intent(in), optional  :: gammai,mui,Xi,Zi
 real :: spsoundi,ponrhoi,mu,X,Z

 mu = gmw
 X = X_in
 Z = Z_in
 tempi = -1.  ! needed because temperature is an in/out to some equations of state, -ve == use own guess
 if (present(mui)) mu = mui
 if (present(Xi)) X = Xi
 if (present(Zi)) Z = Zi
 if (maxvxyzu==4) then
    if (present(gammai)) then
       call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),vxyzui(4),&
                            gamma_local=gammai,mu_local=mu,Xlocal=X,Zlocal=Z,tempi=tempi)
    else
       call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),vxyzui(4),&
                            mu_local=mu,Xlocal=X,Zlocal=Z,tempi=tempi)
    endif
 else
    call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xyzi(1),xyzi(2),xyzi(3),mu_local=mu,Xlocal=X,Zlocal=Z,tempi=tempi)
 endif

end function get_temperature

!----------------------------------------------------------------------------
!+
!  query function to return the internal energyfor calculations with a local
!  mean molecular weight and local adiabatic index
!+
!----------------------------------------------------------------------------
real function get_local_u_internal(gammai, gmwi, gas_temp_local)
 real,         intent(in)    :: gammai, gmwi, gas_temp_local
 real :: ponrhoi

 ponrhoi              = gas_temp_local/(gmwi*temperature_coef)
 get_local_u_internal = ponrhoi/(gammai-1.)

end function get_local_u_internal

!-----------------------------------------------------------------------
!
!  query function to get (gas) temperature given P/rho, assuming fixed
!  mean molecular weight (gmw)
!
!-----------------------------------------------------------------------
real function get_temperature_from_ponrho(ponrho)
 real, intent(in) :: ponrho

 get_temperature_from_ponrho = temperature_coef*gmw*ponrho

end function get_temperature_from_ponrho

!-----------------------------------------------------------------------
!+
!  initialise equation of state (read tables etc.)
!+
!-----------------------------------------------------------------------
subroutine init_eos(eos_type,ierr)
 use units,          only:unit_velocity
 use physcon,        only:mass_proton_cgs,kboltz
 use io,             only:error,warning
 use eos_mesa,       only:init_eos_mesa
 use eos_helmholtz,  only:eos_helmholtz_init
 use eos_piecewise,  only:init_eos_piecewise
 use eos_barotropic, only:init_eos_barotropic
 use eos_shen,       only:init_eos_shen_NL3
 use eos_gasradrec,  only:init_eos_gasradrec
 use dim,            only:maxvxyzu,do_radiation
 integer, intent(in)  :: eos_type
 integer, intent(out) :: ierr

 ierr = 0
 !
 !--Set coefficient to convert P/rho into temperature
 !  calculation will be in cgs; the mean molecular weight, gmw, will be
 !  included in the function call rather than here
 !  c_s^2 = gamma*P/rho = gamma*kT/(gmw*m_p) -> T = P/rho * (gmw*m_p)/k
 !
 temperature_coef = mass_proton_cgs/kboltz * unit_velocity**2

 select case(eos_type)
 case(6)
    !
    !--Check that if using ieos=6, then isink is set properly
    !
    if (isink==0) then
       call error('eos','ieos=6, but isink is not set')
       ierr = ierr_isink_not_set
       return
    endif

 case(8)
    !
    ! barotropic equation of state
    !
    call init_eos_barotropic(polyk,polyk2,ierr)

 case(9)
    !
    ! piecewise polytropic equation of state (similar to barotropic)
    !
    call init_eos_piecewise(ierr)

 case(10)
    !
    !--MESA EoS initialisation
    !
    write(*,'(1x,a,f7.5,a,f7.5)') 'Initialising MESA EoS with X = ',X_in,', Z = ',Z_in
    call init_eos_mesa(X_in,Z_in,ierr)
    if (do_radiation .and. ierr==0) then
       call error('eos','ieos=10, cannot use eos with radiation, will double count radiation pressure')
       ierr=ierr_option_conflict !return error if using radiation and mesa EOS, shouldn't use mesa eos, as it will double count rad pres
    endif

 case(12)
    !
    ! ideal plus radiation
    !
    write(*,'(1x,a,f7.5)') 'Using ideal plus radiation EoS with mu = ',gmw
    if (do_radiation) then
       call error('eos','ieos=12, cannot use eos with radiation, will double count radiation pressure')
       ierr = ierr_option_conflict
    endif

 case(15)

    call eos_helmholtz_init(ierr)

 case(16)

    call init_eos_shen_NL3(ierr)

 case(20)

    call init_eos_gasradrec(ierr)
    if (.not. use_var_comp) then
       write(*,'(a,f7.5,a,f7.5)') 'Assuming fixed composition X = ',X_in,', Z = ',Z_in
    endif
    if (do_radiation) then
       call error('eos','ieos=20, cannot use eos with radiation, will double count radiation pressure')
       ierr = ierr_option_conflict
    endif

 end select
 done_init_eos = .true.

 if (do_radiation .and. iopacity_type==1) call init_eos_mesa(X_in,Z_in,ierr)

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
 !
 open(unit=iuniteos,file="EOS.dat",form='formatted',status='replace')
 write(iuniteos,'("# Equation of state properties; all values in cgs")')
 write(iuniteos,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'rho', &
       2,'P', &
       3,'P/rho', &
       4,'c_s', &
       5,'T'

 dummy = 0.0  ! initialise to avoid compiler warning
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

end subroutine print_eos_to_file
!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos(iunit)
 use infile_utils,   only:write_inopt
 use eos_helmholtz,  only:eos_helmholtz_write_inopt
 use eos_barotropic, only:write_options_eos_barotropic
 use eos_piecewise,  only:write_options_eos_piecewise
 use eos_gasradrec,  only:write_options_eos_gasradrec
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling equation of state'
 call write_inopt(ieos,'ieos','eqn of state (1=isoth;2=adiab;3=locally iso;8=barotropic)',iunit)
#ifndef KROME
 if (.not. eos_outputs_mu(ieos)) call write_inopt(gmw,'mu','mean molecular weight',iunit)
#endif
 select case(ieos)
 case(8)
    call write_options_eos_barotropic(iunit)
 case(9)
    call write_options_eos_piecewise(iunit)
 case(10)
    call write_inopt(X_in,'X','hydrogen mass fraction',iunit)
    call write_inopt(Z_in,'Z','metallicity',iunit)
 case(15) ! helmholtz eos
    call eos_helmholtz_write_inopt(iunit)
 case(20)
    call write_options_eos_gasradrec(iunit)
    if (.not. use_var_comp) then
       call write_inopt(X_in,'X','H mass fraction (ignored if variable composition)',iunit)
       call write_inopt(Z_in,'Z','metallicity (ignored if variable composition)',iunit)
    endif
 end select

end subroutine write_options_eos

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_eos(name,valstring,imatch,igotall,ierr)
 use io,             only:fatal
 use eos_helmholtz,  only:eos_helmholtz_set_relaxflag
 use eos_barotropic, only:read_options_eos_barotropic
 use eos_piecewise,  only:read_options_eos_piecewise
 use eos_gasradrec,  only:read_options_eos_gasradrec
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot  = 0
 character(len=30), parameter  :: label = 'read_options_eos'
 integer :: tmp
 logical :: igotall_barotropic,igotall_piecewise,igotall_gasradrec

 imatch  = .true.
 igotall_barotropic = .true.
 igotall_piecewise = .true.
 igotall_gasradrec = .true.

 select case(trim(name))
 case('ieos')
    read(valstring,*,iostat=ierr) ieos
    ngot = ngot + 1
    if (ieos <= 0 .or. ieos > maxeos) call fatal(label,'equation of state choice out of range')
 case('mu')
    read(valstring,*,iostat=ierr) gmw
    ! not compulsory to read in
    if (gmw <= 0.)  call fatal(label,'mu <= 0')
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
 if (.not.imatch .and. ieos==8) call read_options_eos_barotropic(name,valstring,imatch,igotall_barotropic,ierr)
 if (.not.imatch .and. ieos==9) call read_options_eos_piecewise(name,valstring,imatch,igotall_piecewise,ierr)
 if (.not.imatch .and. ieos==20) call read_options_eos_gasradrec(name,valstring,imatch,igotall_gasradrec,ierr)

 !--make sure we have got all compulsory options (otherwise, rewrite input file)
 igotall = (ngot >= 1) .and. igotall_piecewise .and. igotall_barotropic .and. igotall_gasradrec

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
 use dim,            only:maxvxyzu,gr
 use io,             only:fatal
 use eos_helmholtz,  only:eos_helmholtz_eosinfo
 use eos_barotropic, only:eos_info_barotropic
 use eos_piecewise,  only:eos_info_piecewise
 use eos_gasradrec,  only:eos_info_gasradrec
 integer, intent(in) :: eos_type,iprint
 real, parameter     :: uthermcheck = 3.14159, rhocheck = 23.456

 select case(eos_type)
 case(1,11)
    if (1.0d-5 < polyk .and. polyk < 1.0d3) then
       write(iprint,"(/,a,f10.6)")  ' Isothermal equation of state:     cs^2 = ',polyk
    else
       write(iprint,"(/,a,Es13.6)") ' Isothermal equation of state:     cs^2 = ',polyk
    endif
    if (eos_type==11) write(iprint,*) ' (ZERO PRESSURE) '
 case(2)
    if (use_entropy) then
       write(iprint,"(/,a,f10.6,a,f10.6,a,f10.6)") ' Adiabatic equation of state (evolving ENTROPY): polyk = ',polyk,&
                                                   ' gamma = ',gamma,' gmw = ',gmw
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
       if (gr) then
          write(iprint,"(/,a,f10.6,a,f10.6)") ' Adiabatic equation of state with gamma = ',gamma,' gmw = ',gmw
       else
          write(iprint,"(/,a,f10.6,a,f10.6)") ' Adiabatic equation of state (evolving UTHERM): P = (gamma-1)*rho*u, gamma = ',&
                                              gamma,' gmw = ',gmw
       endif
    else
       write(iprint,"(/,a,f10.6,a,f10.6,a,f10.6)") ' Polytropic equation of state: P = ',polyk,'*rho^',gamma,' gmw = ',gmw
    endif
 case(3)
    write(iprint,"(/,a,f10.6,a,f10.6)") ' Locally isothermal eq of state (R_sph): cs^2_0 = ',polyk,' qfac = ',qfacdisc
 case(6)
    write(iprint,"(/,a,i2,a,f10.6,a,f10.6)") ' Locally (on sink ',isink, &
          ') isothermal eos (R_sph): cs^2_0 = ',polyk,' qfac = ',qfacdisc
 case(8)
    call eos_info_barotropic(polyk,polyk2,iprint)
 case(9)
    call eos_info_piecewise(iprint)
 case(10)
    write(iprint,"(/,a,f10.6,a,f10.6,a,f10.6,a)") ' MESA EoS: X = ',X_in,' Z = ',Z_in,' (1-X-Z = ',1.-X_in-Z_in,')'
 case(12)
    write(iprint,"(/,a,f10.6,a,f10.6)") ' Gas + radiation equation of state: gmw = ',gmw,' gamma = ',gamma
 case(15)
    call eos_helmholtz_eosinfo(iprint)
 case(20)
    call eos_info_gasradrec(iprint)
    if (use_var_comp) then
       write(*,'(1x,a,i1,a)') 'Using variable composition'
    else
       write(*,'(1x,a,f10.6,a,f10.6)') 'Using fixed composition X = ',X_in,", Z = ",Z_in
    endif
 end select
 write(iprint,*)

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

end function en_from_utherm

!----------------------------------------------------------------
!+
!  Get recombination energy (per unit mass) assumming complete
!  ionisation
!+
!----------------------------------------------------------------
subroutine calc_rec_ene(XX,YY,e_rec)
 real, intent(in)  :: XX, YY
 real, intent(out) :: e_rec
 real              :: e_H2,e_HI,e_HeI,e_HeII
 real, parameter   :: e_ion_H2   = 1.312d13, & ! ionisation energies in erg/mol
                      e_ion_HI   = 4.36d12, &
                      e_ion_HeI  = 2.3723d13, &
                      e_ion_HeII = 5.2505d13

 ! XX     : Hydrogen mass fraction
 ! YY     : Helium mass fraction
 ! e_rec  : Total ionisation energy due to H2, HI, HeI, and HeII

 e_H2   = 0.5 * XX * e_ion_H2
 e_HI   = XX * e_ion_HI
 e_HeI  = 0.25 * YY * e_ion_HeI
 e_HeII = 0.25 * YY * e_ion_HeII
 e_rec  = e_H2 + e_HI + e_HeI + e_HeII

end subroutine calc_rec_ene

!----------------------------------------------------------------
!+
!  Calculate temperature and specific internal energy from
!  pressure and density. Inputs and outputs are in cgs units.
!
!  Note on composition:
!  For ieos=2 and 12, mu_local is an input, X & Z are not used
!  For ieos=10, mu_local is not used
!  For ieos=20, mu_local is not used but available as an output
!+
!----------------------------------------------------------------
subroutine calc_temp_and_ene(eos_type,rho,pres,ene,temp,ierr,guesseint,mu_local,X_local,Z_local)
 use physcon,          only:kb_on_mh
 use eos_idealplusrad, only:get_idealgasplusrad_tempfrompres,get_idealplusrad_enfromtemp
 use eos_mesa,         only:get_eos_eT_from_rhop_mesa
 use eos_gasradrec,    only:calc_uT_from_rhoP_gasradrec
 integer, intent(in)        :: eos_type
 real, intent(in)           :: rho,pres
 real, intent(inout)        :: ene,temp
 real, intent(in), optional :: guesseint,X_local,Z_local
 real, intent(inout), optional :: mu_local
 integer, intent(out)       :: ierr
 real                       :: mu,X,Z

 ierr = 0
 mu = gmw
 X = X_in
 Z = Z_in
 if (present(mu_local)) mu = mu_local
 if (present(X_local)) X = X_local
 if (present(Z_local)) Z = Z_local
 select case(eos_type)
 case(2) ! Ideal gas
    temp = pres / (rho * kb_on_mh) * mu
    ene = pres / ( (gamma-1.) * rho)
 case(12) ! Ideal gas + radiation
    call get_idealgasplusrad_tempfrompres(pres,rho,mu,temp)
    call get_idealplusrad_enfromtemp(rho,temp,mu,gamma,ene)
 case(10) ! MESA EoS
    call get_eos_eT_from_rhop_mesa(rho,pres,ene,temp,guesseint)
 case(20) ! Ideal gas + radiation + recombination (from HORMONE, Hirai et al., 2020)
    call calc_uT_from_rhoP_gasradrec(rho,pres,X,1.-X-Z,temp,ene,mu,ierr)
    if (present(mu_local)) mu_local = mu
 case default
    ierr = 1
 end select

end subroutine calc_temp_and_ene

!-----------------------------------------------------------------------
!+
!  Calculates specific entropy (gas + radiation + recombination)
!  up to an additive integration constant, from density and pressure.
!+
!-----------------------------------------------------------------------
function entropy(rho,pres,mu_in,ientropy,eint_in,ierr)
 use io,                only:fatal
 use physcon,           only:radconst,kb_on_mh
 use eos_idealplusrad,  only:get_idealgasplusrad_tempfrompres
 use eos_mesa,          only:get_eos_eT_from_rhop_mesa
 use mesa_microphysics, only:getvalue_mesa
 real, intent(in)               :: rho,pres
 real, intent(in), optional     :: mu_in,eint_in
 integer, intent(in)            :: ientropy
 integer, intent(out), optional :: ierr
 real                           :: mu,entropy,logentropy,temp,eint

 if (present(ierr)) ierr=0

 mu = mu_in
 select case(ientropy)
 case(1) ! Include only gas entropy (up to additive constants)
    temp = pres * mu / (rho * kb_on_mh)
    entropy = kb_on_mh / mu * log(temp**1.5/rho)

 case(2) ! Include both gas and radiation entropy (up to additive constants)
    temp = pres * mu / (rho * kb_on_mh) ! Guess for temp
    call get_idealgasplusrad_tempfrompres(pres,rho,mu,temp) ! First solve for temp from rho and pres
    entropy = kb_on_mh / mu * log(temp**1.5/rho) + 4.*radconst*temp**3 / (3.*rho)

 case(3) ! Get entropy from MESA tables if using MESA EoS
    if (ieos /= 10 .and. ieos /= 20) call fatal('eos','Using MESA tables to calculate S from rho and pres, but not using MESA EoS')

    if (present(eint_in)) then
       eint = eint_in
    else
       call get_eos_eT_from_rhop_mesa(rho,pres,eint,temp)
    endif

    ! Get entropy from rho and eint from MESA tables
    if (present(ierr)) then
       call getvalue_mesa(rho,eint,9,logentropy,ierr)
    else
       call getvalue_mesa(rho,eint,9,logentropy)
    endif
    entropy = 10.**logentropy

 case default
    entropy = 0.
    call fatal('eos','Unknown ientropy (can only be 1, 2, or 3)')
 end select

end function entropy

!-----------------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using Newton-Raphson
!  method
!+
!-----------------------------------------------------------------------
subroutine get_rho_from_p_s(pres,S,rho,mu,rhoguess,ientropy)
 use physcon, only:kb_on_mh
 real, intent(in)    :: pres,S,mu,rhoguess
 real, intent(inout) :: rho
 real                :: srho,srho_plus_dsrho,S_plus_dS,dSdsrho
 real(kind=8)        :: corr
 real, parameter     :: eoserr=1d-9,dfac=1d-12
 integer, intent(in) :: ientropy

 ! We apply the Newton-Raphson method directly to rho^1/2 ("srho") instead
 ! of rho since S(rho) cannot take a negative argument.
 srho = sqrt(rhoguess) ! Initial guess
 corr = huge(corr);
 do while (abs(corr) > eoserr*abs(srho))
    ! First calculate dS/dsrho
    srho_plus_dsrho = srho * (1. + dfac)
    S_plus_dS = entropy(srho_plus_dsrho**2,pres,mu,ientropy)
    dSdsrho = (S_plus_dS - entropy(srho**2,pres,mu,ientropy)) / (srho_plus_dsrho - srho)
    corr = ( entropy(srho**2,pres,mu,ientropy) - S ) / dSdsrho
    srho = srho - corr
 enddo
 rho = srho**2

end subroutine get_rho_from_p_s

!-----------------------------------------------------------------------
!+
!  Calculate mean molecular weight from X and Z, assumming complete
!  ionisation
!+
!-----------------------------------------------------------------------
real function get_mean_molecular_weight(XX,ZZ) result(mu)
 real, intent(in) :: XX,ZZ
 real :: YY

 YY = 1.-XX-ZZ
 mu = 1./(2.*XX + 0.75*YY + 0.5*ZZ)

end function get_mean_molecular_weight

!-----------------------------------------------------------------------
!+
!  write relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_eos(ieos,hdr,ierr)
 use dump_utils,        only:dump_h,add_to_rheader,add_to_iheader
 integer,      intent(in)    :: ieos
 type(dump_h), intent(inout) :: hdr
 integer,      intent(out)   :: ierr

 call add_to_iheader(isink,'isink',hdr,ierr)
 call add_to_rheader(gamma,'gamma',hdr,ierr)
 call add_to_rheader(1.5*polyk,'RK2',hdr,ierr)
 call add_to_rheader(polyk2,'polyk2',hdr,ierr)
 call add_to_rheader(qfacdisc,'qfacdisc',hdr,ierr)
 call add_to_rheader(qfacdisc2,'qfacdisc2',hdr,ierr)

 if (ieos==7) then
   call add_to_rheader(alpha_z,'alpha_z',hdr,ierr)
   call add_to_rheader(beta_z,'beta_z',hdr,ierr)
   call add_to_rheader(z0,'z0',hdr,ierr)

 endif

end subroutine write_headeropts_eos

!-----------------------------------------------------------------------
!+
!  read relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_eos(ieos,hdr,ierr)
 use dump_utils,        only:dump_h, extract
 use io,                only:iprint,id,master
 use dim,               only:use_krome,maxvxyzu
 integer,      intent(in)  :: ieos
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 real :: RK2


 call extract('gamma',gamma,hdr,ierr)
 call extract('RK2',rk2,hdr,ierr)
 polyk = 2./3.*rk2
 if (id==master) then
    if (maxvxyzu >= 4) then
       if (use_krome) then
          write(iprint,*) 'KROME eos: initial gamma = 1.666667'
       else
          write(iprint,*) 'adiabatic eos: gamma = ',gamma
       endif
    else
       write(iprint,*) 'setting isothermal sound speed^2 (polyk) = ',polyk,' gamma = ',gamma
       if (polyk <= tiny(polyk)) write(iprint,*) 'WARNING! sound speed zero in dump!, polyk = ',polyk
    endif
 endif
 call extract('polyk2',polyk2,hdr,ierr)
 call extract('qfacdisc',qfacdisc,hdr,ierr)
 call extract('qfacdisc2',qfacdisc2,hdr,ierr)
 call extract('isink',isink,hdr,ierr)

 if (abs(gamma-1.) > tiny(gamma) .and. maxvxyzu < 4) then
    write(*,*) 'WARNING! compiled for isothermal equation of state but gamma /= 1, gamma=',gamma
 endif

 ierr = 0
 if (ieos==3 .or. ieos==6 .or. ieos==7) then
    if (qfacdisc <= tiny(qfacdisc)) then
       write(iprint,*) 'ERROR: qfacdisc <= 0'
       ierr = 2
    else
       write(iprint,*) 'qfacdisc = ',qfacdisc
    endif
  endif

  if (ieos==7) then
    call extract('alpha_z',alpha_z,hdr,ierr)
    call extract('beta_z', beta_z, hdr,ierr)
    call extract('z0',z0,hdr,ierr)
    if (qfacdisc2 <= tiny(qfacdisc2)) then
       write(iprint,*) 'ERROR: qfacdisc2 <= 0'
       ierr = 2
    else
       write(iprint,*) 'qfacdisc2 = ',qfacdisc2
    endif
  endif

end subroutine read_headeropts_eos


end module eos
