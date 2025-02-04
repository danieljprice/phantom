!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos
!
! This module contains stuff to do with the equation of state
!  Current options:
!     1 = isothermal eos
!     2 = adiabatic/polytropic eos
!     3 = eos for a locally isothermal disc as in Lodato & Pringle (2007)
!     4 = GR isothermal
!     5 = polytropic EOS with varying mu and gamma depending on H2 formation
!     6 = eos for a locally isothermal disc as in Lodato & Pringle (2007),
!         centered on a sink particle
!     7 = z-dependent locally isothermal eos
!     8 = Barotropic eos
!     9 = Piecewise polytrope
!    10 = MESA EoS
!    11 = isothermal eos with zero pressure
!    12 = ideal gas with radiation pressure
!    13 = locally isothermal prescription from Farris et al. (2014) generalised for generic hierarchical systems
!    14 = locally isothermal prescription from Farris et al. (2014) for binary system
!    15 = Helmholtz free energy eos
!    16 = Shen eos
!    17 = polytropic EOS with varying mu (depending on H2 formation)
!    20 = Ideal gas + radiation + various forms of recombination energy from HORMONE (Hirai et al., 2020)
!    23 = Hypervelocity Impact of solids-fluids from Tillotson EOS (Tillotson 1962 - implemented by Brundage A. 2013
!    24 = read tabulated eos (for use with icooling == 9)
!
! :References:
!    Lodato & Pringle (2007)
!    Hirai et al. (2020)
!
! :Owner: Alison Young
!
! :Runtime parameters:
!   - X           : *H mass fraction (ignored if variable composition)*
!   - Z           : *metallicity (ignored if variable composition)*
!   - ieos        : *eqn of state (1=isoth;2=adiab;3=locally iso;8=barotropic)*
!   - metallicity : *metallicity*
!   - mu          : *mean molecular weight*
!
! :Dependencies: dim, dump_utils, eos_HIIR, eos_barotropic, eos_gasradrec,
!   eos_helmholtz, eos_idealplusrad, eos_mesa, eos_piecewise, eos_shen,
!   eos_stamatellos, eos_stratified, eos_tillotson, infile_utils, io,
!   mesa_microphysics, part, physcon, units
!
 use part,          only:ien_etotal,ien_entropy,ien_type
 use dim,           only:gr
 use eos_gasradrec, only:irecomb
 implicit none
 integer, parameter, public :: maxeos = 24
 real,               public :: polyk, polyk2, gamma
 real,               public :: qfacdisc = 0.75, qfacdisc2 = 0.75
 logical,            public :: extract_eos_from_hdr = .false.
 integer,            public :: isink = 0.

 public  :: equationofstate,setpolyk,eosinfo,get_mean_molecular_weight
 public  :: get_TempPresCs,get_spsound,get_temperature,get_pressure,get_cv
 public  :: eos_is_non_ideal,eos_outputs_mu,eos_outputs_gasP
 public  :: get_local_u_internal,get_temperature_from_u
 public  :: calc_rec_ene,calc_temp_and_ene,entropy,get_rho_from_p_s,get_u_from_rhoT
 public  :: calc_rho_from_PT,get_entropy,get_p_from_rho_s
 public  :: init_eos,finish_eos,write_options_eos,read_options_eos
 public  :: write_headeropts_eos,read_headeropts_eos
 public  :: eos_requires_isothermal,eos_requires_polyk
 public  :: eos_is_not_implemented

 public :: irecomb  ! propagated from eos_gasradrec

 private

 integer, public :: ieos          = 1
 integer, public :: iopacity_type = 0      ! used for radiation
 real,    public :: gmw           = 2.381  ! default mean molecular weight
 real,    public :: X_in          = 0.74   ! default metallicities
 real,    public :: Z_in          = 0.02   ! default metallicities
 logical, public :: use_var_comp  = .false. ! use variable composition
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
! Default temperature prescription for vertical stratification (0=MAPS, 1=Dartois)
!
 integer, public:: istrat = 0.
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
subroutine equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi,eni,gamma_local,mu_local,Xlocal,Zlocal,isionised)
 use io,            only:fatal,error,warning
 use part,          only:xyzmh_ptmass, nptmass
 use units,         only:unit_density,unit_pressure,unit_ergg,unit_velocity
 use physcon,       only:Rg,radconst,kb_on_mh
 use eos_mesa,      only:get_eos_pressure_temp_gamma1_mesa,get_eos_1overmu_mesa
 use eos_helmholtz, only:eos_helmholtz_pres_sound
 use eos_shen,      only:eos_shen_NL3
 use eos_idealplusrad, only:get_idealplusrad_pres,get_idealplusrad_temp,get_idealplusrad_spsoundi
 use eos_gasradrec,    only:equationofstate_gasradrec
 use eos_stratified,   only:get_eos_stratified
 use eos_barotropic,   only:get_eos_barotropic
 use eos_piecewise,    only:get_eos_piecewise
 use eos_tillotson,    only:equationofstate_tillotson
 use eos_stamatellos
 use eos_HIIR,         only:get_eos_HIIR_iso,get_eos_HIIR_adiab
 integer, intent(in)    :: eos_type
 real,    intent(in)    :: rhoi,xi,yi,zi
 real,    intent(out)   :: ponrhoi,spsoundi
 real,    intent(inout) :: tempi
 real,    intent(in),    optional :: eni
 real,    intent(inout), optional :: mu_local,gamma_local
 real,    intent(in)   , optional :: Xlocal,Zlocal
 logical, intent(in),    optional :: isionised
 integer :: ierr, i
 real    :: r1,r2
 real    :: mass_r, mass ! defined for generalised Farris prescription
 real    :: gammai,temperaturei,mui,imui,X_i,Z_i
 real    :: cgsrhoi,cgseni,cgspresi,presi,gam1,cgsspsoundi
 real    :: uthermconst,kappaBar,kappaPart
 real    :: enthi,pondensi
 logical :: isionisedi
 !
 ! Check to see if equation of state is compatible with GR cons2prim routines
 !
 if (gr .and. .not.any((/2,4,11,12/)==eos_type)) then
    ponrhoi = 0.; spsoundi = 0. ! avoid compiler warning
    call fatal('eos','GR currently only works for ieos=2,12 or 11',&
         var='eos_type',val=real(eos_type))
 endif

 gammai = gamma
 mui    = gmw
 X_i    = X_in
 Z_i    = Z_in
 if (present(gamma_local)) gammai = gamma_local
 if (present(mu_local)) mui = mu_local
 if (present(Xlocal)) X_i = Xlocal
 if (present(Zlocal)) Z_i = Zlocal
 if (present(isionised)) isionisedi = isionised

 select case(eos_type)
 case(1)
!
!--Isothermal eos
!
!  :math:`P = c_s^2 \rho`
!
!  where :math:`c_s^2 \equiv K` is a constant stored in the dump file header
!
    ponrhoi  = polyk
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*mui*ponrhoi

 case(2,5,17)
!
!--Adiabatic equation of state (code default)
!
!  :math:`P = (\gamma - 1) \rho u`
!
!  if the code is compiled with ISOTHERMAL=yes, ieos=2 gives a polytropic eos:
!
!  :math:`P = K \rho^\gamma`
!
!  where K is a global constant specified in the dump header
!
    if (gammai < tiny(gammai)) call fatal('eos','gamma not set for adiabatic eos',var='gamma',val=gammai)

    if (gr) then
       if (.not. present(eni)) call fatal('eos','GR call to equationofstate requires thermal energy as input!')
       if (eni < 0.) call fatal('eos','utherm < 0',var='u',val=eni)
       if (gammai <= 1.) then
          spsoundi = 0.; ponrhoi = 0. ! avoid compiler warning
          call fatal('eos','GR not compatible with isothermal equation of state, yet...',var='gamma',val=gammai)
       elseif (gammai > 1.0001) then
          pondensi = (gammai-1.)*eni   ! eni is the thermal energy
          enthi = 1. + eni + pondensi    ! enthalpy
          spsoundi = sqrt(gammai*pondensi/enthi)
          ponrhoi = pondensi ! With GR this routine actually outputs pondensi (i.e. pressure on primitive density, not conserved.)
       endif
    else
       if (present(eni)) then
          if (eni < 0.) then
             !write(iprint,'(a,Es18.4,a,4Es18.4)')'Warning: eos: u = ',eni,' < 0 at {x,y,z,rho} = ',xi,yi,zi,rhoi
             call fatal('eos','utherm < 0',var='u',val=eni)
          endif
          if (gammai > 1.0001) then
             ponrhoi = (gammai-1.)*eni   ! use this if en is thermal energy
          else
             ponrhoi = 2./3.*eni ! en is thermal energy and gamma = 1
          endif
       else
          ponrhoi = polyk*rhoi**(gammai-1.)
       endif
       spsoundi = sqrt(gammai*ponrhoi)
    endif

    tempi = temperature_coef*mui*ponrhoi

 case(3)
!
!--Locally isothermal disc as in Lodato & Pringle (2007) where
!
!  :math:`P = c_s^2 (r) \rho`
!
!  sound speed (temperature) is prescribed as a function of radius using:
!
!  :math:`c_s = c_{s,0} r^{-q}` where :math:`r = \sqrt{x^2 + y^2 + z^2}`
!
    ponrhoi  = polyk*(xi**2 + yi**2 + zi**2)**(-qfacdisc) ! polyk is cs^2, so this is (R^2)^(-q)
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*mui*ponrhoi

 case(4)
!
!--Isothermal equation of state for GR, enforcing cs = constant
!
!  .. WARNING:: this is experimental: use with caution
!
    uthermconst = polyk
    ponrhoi  = (gammai-1.)*uthermconst
    spsoundi = sqrt(ponrhoi/(1.+uthermconst))
    tempi    = temperature_coef*mui*ponrhoi

 case(6)
!
!--Locally isothermal disc centred on sink particle
!
!  As in ieos=3 but in this version radius is taken with respect to a designated
!  sink particle (by default the first sink particle in the simulation)
!
    ponrhoi  = polyk*((xi-xyzmh_ptmass(1,isink))**2 + (yi-xyzmh_ptmass(2,isink))**2 + &
                      (zi-xyzmh_ptmass(3,isink))**2)**(-qfacdisc) ! polyk is cs^2, so this is (R^2)^(-q)
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*mui*ponrhoi

 case(7)
!
!--Vertically stratified equation of state
!
!  sound speed is prescribed as a function of (cylindrical) radius R and
!  height z above the x-y plane
!
!  .. WARNING:: should not be used for misaligned discs
!
    call get_eos_stratified(istrat,xi,yi,zi,polyk,polyk2,qfacdisc,qfacdisc2,alpha_z,beta_z,z0,ponrhoi,spsoundi)
    tempi = temperature_coef*mui*ponrhoi

 case(8)
!
!--Barotropic equation of state
!
!  :math:`P = K \rho^\gamma`
!
!  where the value of gamma (and K) are prescribed functions of density
!
    call get_eos_barotropic(rhoi,polyk,polyk2,ponrhoi,spsoundi,gammai)
    tempi = temperature_coef*mui*ponrhoi

 case(9)
!
!--Piecewise Polytropic equation of state
!
!  :math:`P = K \rho^\gamma`
!
!  where the value of gamma (and K) are a prescribed function of density.
!  Similar to ieos=8 but with different defaults and slightly different
!  functional form
!
    call get_eos_piecewise(rhoi,ponrhoi,spsoundi,gammai)
    tempi = temperature_coef*mui*ponrhoi

 case(10)
!
!--MESA equation of state
!
!  a tabulated equation of state including gas, radiation pressure
!  and ionisation/dissociation. MESA is a stellar evolution code, so
!  this equation of state is designed for matter inside stars
!
    cgsrhoi = rhoi * unit_density
    cgseni  = eni * unit_ergg
    call get_eos_pressure_temp_gamma1_mesa(cgsrhoi,cgseni,cgspresi,temperaturei,gam1,ierr)
    presi = cgspresi / unit_pressure

    ponrhoi  = presi / rhoi
    spsoundi = sqrt(gam1*ponrhoi)
    tempi    = temperaturei
    if (present(gamma_local)) gamma_local = gam1 ! gamma is an output
    if (present(mu_local)) mu_local = 1./get_eos_1overmu_mesa(cgsrhoi,cgseni)
    if (ierr /= 0) call warning('eos_mesa','extrapolating off tables')

 case(11)
!
!--Equation of state with pressure and temperature equal to zero
!
!  :math:`P = 0`
!
!  useful for simulating test particle dynamics using SPH particles
!
    ponrhoi  = 0.
    spsoundi = sqrt(polyk)
    tempi    = 0.

 case(12)
!
!--Ideal gas plus radiation pressure
!
!  :math:`P = (\gamma - 1) \rho u`
!
!  but solved by first solving the quartic equation:
!
!  :math:`u = \frac32 \frac{k_b T}{\mu m_H} + \frac{a T^4}{\rho}`
!
!  for temperature (given u), then solving for pressure using
!
!  :math:`P = \frac{k_b T}{\mu m_H} + \frac13 a T^4`
!
!  hence in this equation of state gamma (and temperature) are an output
!
    temperaturei = tempi  ! Required as initial guess
    cgsrhoi      = rhoi * unit_density
    cgseni       = eni * unit_ergg
    call get_idealplusrad_temp(cgsrhoi,cgseni,mui,temperaturei,ierr)
    call get_idealplusrad_pres(cgsrhoi,temperaturei,mui,cgspresi)
    call get_idealplusrad_spsoundi(cgsrhoi,cgspresi,cgseni,spsoundi,gammai)
    if (present(gamma_local)) gamma_local = gammai ! gamma is an output
    spsoundi = spsoundi / unit_velocity
    presi    = cgspresi / unit_pressure
    ponrhoi  = presi / rhoi
    tempi    = temperaturei
    if (ierr /= 0) call warning('eos_idealplusrad','temperature iteration did not converge')


 case(13)
!
!--Locally isothermal eos for generic hierarchical system
!
!  Assuming all sink particles are stars.
!  Generalisation of Farris et al. (2014; for binaries) to N stars.
!  For two sink particles this is identical to ieos=14
!
    mass_r = 0
    mass = 0

    do i=1,nptmass
       mass_r = mass_r+xyzmh_ptmass(4,i)/sqrt((xi-xyzmh_ptmass(1,i))**2 + (yi-xyzmh_ptmass(2,i))**2 + (zi-xyzmh_ptmass(3,i))**2)
       mass = mass + xyzmh_ptmass(4,i)
    enddo
    ponrhoi=polyk*(mass_r)**(2*qfacdisc)/mass**(2*qfacdisc)
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*mui*ponrhoi


 case(14)
!
!--Locally isothermal eos from Farris et al. (2014) for binary system
!
!  uses the locations of the first two sink particles
!
    r1 = sqrt((xi-xyzmh_ptmass(1,1))**2+(yi-xyzmh_ptmass(2,1))**2 + (zi-xyzmh_ptmass(3,1))**2)
    r2 = sqrt((xi-xyzmh_ptmass(1,2))**2+(yi-xyzmh_ptmass(2,2))**2 + (zi-xyzmh_ptmass(3,2))**2)
    ponrhoi=polyk*(xyzmh_ptmass(4,1)/r1+xyzmh_ptmass(4,2)/r2)**(2*qfacdisc)/(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2))**(2*qfacdisc)
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*mui*ponrhoi

 case(15)
!
!--Helmholtz equation of state for fully ionized gas
!
!  .. WARNING:: not widely tested in phantom, better to use ieos=10
!
    call eos_helmholtz_pres_sound(tempi, rhoi, ponrhoi, spsoundi, eni)

 case(16)
!
!--Shen (2012) equation of state for neutron stars
!
!  this equation of state requires evolving temperature as the energy variable
!
!  .. WARNING:: not tested: use with caution
!
    if (present(eni)) then
       cgsrhoi = rhoi * unit_density
       !note eni is actually tempi
       call eos_shen_NL3(cgsrhoi,eni,0.05,cgspresi,cgsspsoundi)
       spsoundi=cgsspsoundi / unit_velocity
       presi = cgspresi / unit_pressure
       ponrhoi = presi / rhoi
       tempi = eni
       call warning('eos','Not sure if this is correct now that temperature is always passed into eos')
    else
       spsoundi = 0.; presi = 0.; ponrhoi = 0.; tempi = 0. ! to avoid compiler warnings
       call fatal('eos','tried to call NL3 eos without passing temperature')
    endif

 case(20)
!
!--Gas + radiation + various forms of recombination
!
!  from HORMONE, Hirai+2020, as used in Lau+2022b
!
    cgsrhoi = rhoi * unit_density
    cgseni  = eni * unit_ergg
    imui    = 1./mui
    if (tempi > 0.) then
       temperaturei = tempi
    else
       temperaturei = min(0.67 * cgseni * mui / Rg, (cgseni*cgsrhoi/radconst)**0.25)
    endif
    call equationofstate_gasradrec(cgsrhoi,cgseni*cgsrhoi,temperaturei,imui,X_i,1.-X_i-Z_i,cgspresi,cgsspsoundi,gammai)
    ponrhoi  = real(cgspresi / (unit_pressure * rhoi))
    spsoundi = real(cgsspsoundi / unit_velocity)
    tempi    = temperaturei
    if (present(mu_local)) mu_local = 1./imui
    if (present(gamma_local)) gamma_local = gammai

 case(21)
!
!--HII region two temperature "equation of state"
!
!  flips the temperature depending on whether a particle is ionised or not,
!  use with ISOTHERMAL=yes
!
    call get_eos_HIIR_iso(polyk,temperature_coef,mui,tempi,ponrhoi,spsoundi,isionisedi)
 case(22)
!
!--Same as ieos=21 but sets the thermal energy
!
!  for use when u is stored (ISOTHERMAL=no)
!
    call get_eos_HIIR_adiab(polyk,temperature_coef,mui,tempi,ponrhoi,rhoi,eni,gammai,spsoundi,isionisedi)
 case(23)
!
!--Tillotson (1962) equation of state for solids (basalt, granite, ice, etc.)
!
!  Implementation from Benz et al. (1986), Asphaug & Melosh (1993) and Kegerreis et al. (2019)
!
!  In the compressed (:math:`\rho > \rho_0`) or cold (:math:`u < u_{\rm iv}`) state gives
!
!  :math:`P_c = \left[a + \frac{b}{(u/(u_0 \eta^2) + 1}\right]\rho u + A \mu + B\mu^2`
!
!  where :math:`\eta = rho/rho_0`, :math:`\mu = \eta - 1`, u is the specific internal energy
!  and a,b,A,B and :math:`u_0` are input parameters chosen for a particular material
!
!  In the hot, expanded state (:math:`\rho < \rho_0` and :math:`u > u_{\rm iv}`) gives
!
!  :math:`P_e = a\rho u + \left[\frac{b\rho u}{u/(u_0 \eta^2) + 1} + A\mu \exp{-\beta \nu} \right] \exp(-\alpha \nu^2)`
!
!  where :math:`\nu = \rho/\rho_0 - 1`. In the intermediate state pressure is interpolated using
!
!  :math:`P = \frac{(u - u_{\rm iv}) P_e + (u_{\rm cv} - u) P_c}{u_{\rm cv} - u_{\rm iv}}.`
!
!  When using this equation of state bodies should be set up with uniform density equal, or close to the
!  reference density :math:`\rho_0`, e.g. 2.7 g/cm^3 for basalt
!
    cgsrhoi = rhoi * unit_density
    cgseni  = eni * unit_ergg
    call equationofstate_tillotson(cgsrhoi,cgseni,cgspresi,cgsspsoundi,gammai)
    ponrhoi  = real(cgspresi / (unit_pressure * rhoi))
    spsoundi = real(cgsspsoundi / unit_velocity)
    !  tempi    = 0. !temperaturei
case (24)
!--Interpolate tabulated EoS from Stamatellos et al. (2007).
!   
!  Tabulated equation of state with opacities from Lombardi et al. 2015. For use
!  with icooling = 9, the radiative cooling approximation (Young et al. 2024).
!
    if (eni < 0.) then
       call fatal('eos (stamatellos)','utherm < 0',var='u',val=eni)
    endif
    cgsrhoi = rhoi * unit_density
    cgseni = eni * unit_ergg
    call getopac_opdep(cgseni,cgsrhoi,kappaBar,kappaPart,tempi,mui)
    cgspresi = kb_on_mh*cgsrhoi*tempi/mui
    presi = cgspresi/unit_pressure
    ponrhoi = presi/rhoi
    gammai = 1.d0 + presi/(eni*rhoi)
    spsoundi = sqrt(gammai*ponrhoi)

 case default
    spsoundi = 0. ! avoids compiler warnings
    ponrhoi  = 0.
    tempi    = 0.
    call fatal('eos','unknown equation of state')
 end select

end subroutine equationofstate

!-----------------------------------------------------------------------
!+
!  initialise equation of state (read tables etc.)
!+
!-----------------------------------------------------------------------
subroutine init_eos(eos_type,ierr)
 use units,          only:unit_velocity
 use physcon,        only:Rg
 use io,             only:error,warning
 use eos_mesa,       only:init_eos_mesa
 use eos_helmholtz,  only:eos_helmholtz_init
 use eos_piecewise,  only:init_eos_piecewise
 use eos_barotropic, only:init_eos_barotropic
 use eos_shen,       only:init_eos_shen_NL3
 use eos_gasradrec,  only:init_eos_gasradrec
 use eos_stamatellos,only:read_optab,init_coolra,eos_file
 use eos_HIIR,       only:init_eos_HIIR
 use dim,            only:maxvxyzu,do_radiation
 use eos_tillotson,  only:init_eos_tillotson
 integer, intent(in)  :: eos_type
 integer, intent(out) :: ierr
 integer              :: ierr_mesakapp,ierr_ra

 ierr = 0
 !
 !--Set coefficient to convert P/rho into temperature
 !  calculation will be in cgs; the mean molecular weight, gmw, will be
 !  included in the function call rather than here
 !  c_s^2 = gamma*P/rho = gamma*kT/(gmw*m_p) -> T = P/rho * (gmw*m_p)/k
 !
 temperature_coef = unit_velocity**2  / Rg

 select case(eos_type)
 case(6)
    !
    !--Check that if using ieos=6, then isink is set properly
    !
    if (isink==0) then
       call error('eos','ieos=6, but isink is not set, setting to 1')
       isink = 1
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

 case(21,22)

    call init_eos_HIIR()

 case(23)
    call init_eos_tillotson(ierr)

 case(24)
    call read_optab(eos_file,ierr_ra)
    if (ierr_ra > 0) call warning('init_eos','Failed to read EOS file')
    call init_coolra
 end select
 done_init_eos = .true.

 if (do_radiation .and. iopacity_type==1) then
    write(*,'(1x,a,f7.5,a,f7.5)') 'Using radiation with MESA opacities. Initialising MESA EoS with X = ',X_in,', Z = ',Z_in
    call init_eos_mesa(X_in,Z_in,ierr_mesakapp)
    ierr = max(ierr,ierr_mesakapp)
 endif

end subroutine init_eos

!-----------------------------------------------------------------------
!+
!  finish equation of state
!+
!-----------------------------------------------------------------------
subroutine finish_eos(eos_type,ierr)
 use eos_mesa,       only:finish_eos_mesa
 use eos_helmholtz,  only:eos_helmholtz_finish
 use eos_stamatellos,only:finish_coolra
 integer, intent(in)  :: eos_type
 integer, intent(out) :: ierr

 ierr = 0

 select case(eos_type)
 case(15)
    !
    !--Helmholtz eos deallocation
    !
    call eos_helmholtz_finish(ierr)
 case(10)
    !
    !--MESA EoS deallocation
    !
    call finish_eos_mesa
 case(24)
    ! Stamatellos deallocation
    call finish_coolra
 end select
 done_init_eos=.false.

end subroutine finish_eos

!-----------------------------------------------------------------------
!+
!  Calculate gas temperature, sound speed, and pressure.
!  This will be required for various analysis routines if eos_vars
!  is not saved in the dump files
!+
!-----------------------------------------------------------------------
subroutine get_TempPresCs(eos_type,xyzi,vxyzui,rhoi,tempi,presi,spsoundi,gammai,mui,Xi,Zi)
  use dim, only:maxvxyzu
  use io,  only:warning
 integer, intent(in)              :: eos_type
 real,    intent(in)              :: vxyzui(:),xyzi(:),rhoi
 real,    intent(inout)           :: tempi
 real,    intent(out),   optional :: presi,spsoundi
 real,    intent(inout), optional :: gammai,mui
 real,    intent(in),    optional :: Xi,Zi
 real                             :: csi,ponrhoi,mu,X,Z
 logical                          :: use_gamma

 mu = gmw
 X  = X_in
 Z  = Z_in
 if (present(mui)) mu = mui
 if (present(Xi))  X = Xi
 if (present(Zi))  Z = Zi
 use_gamma = .false.
 if (present(gammai)) then
    if (gammai > 0.) use_gamma = .true.
 endif

 if (maxvxyzu==4) then
    if (vxyzui(4) < 0.) call warning('eos','ui negative in eos')
    if (use_gamma) then
       call equationofstate(eos_type,ponrhoi,csi,rhoi,xyzi(1),xyzi(2),xyzi(3),tempi,vxyzui(4),&
                            gamma_local=gammai,mu_local=mu,Xlocal=X,Zlocal=Z)
    else
       call equationofstate(eos_type,ponrhoi,csi,rhoi,xyzi(1),xyzi(2),xyzi(3),tempi,vxyzui(4),&
                            mu_local=mu,Xlocal=X,Zlocal=Z)
    endif
 else
    call equationofstate(eos_type,ponrhoi,csi,rhoi,xyzi(1),xyzi(2),xyzi(3),tempi,mu_local=mu)
 endif

 if (present(presi))    presi    = ponrhoi*rhoi
 if (present(spsoundi)) spsoundi = csi
 if (present(mui))     mui = mu
 if (present(gammai)) gammai = gamma

end subroutine get_TempPresCs

!-----------------------------------------------------------------------
!+
!  Wrapper function to calculate sound speed
!+
!-----------------------------------------------------------------------
real function get_spsound(eos_type,xyzi,rhoi,vxyzui,gammai,mui,Xi,Zi)
 integer, intent(in)             :: eos_type
 real,    intent(in)             :: xyzi(:),rhoi
 real,    intent(in)             :: vxyzui(:)
 real,    intent(in),    optional :: Xi,Zi
 real,    intent(inout), optional :: gammai,mui
 real                            :: spsoundi,tempi,gam,mu,X,Z

 !set defaults for variables not passed in
 mu    = gmw
 X     = X_in
 Z     = Z_in
 tempi = -1.  ! needed because temperature is an in/out to some equations of state, -ve == use own guess
 gam   = -1.  ! to indicate gamma is not being passed in
 if (present(Xi))  X  = Xi
 if (present(Zi))  Z  = Zi
 if (present(gammai)) gam = gammai
 if (present(mui))    mu = mui

 call get_TempPresCs(eos_type,xyzi,vxyzui,rhoi,tempi,spsoundi=spsoundi,gammai=gam,mui=mu,Xi=X,Zi=Z)

 get_spsound = spsoundi

 if (present(mui))    mui = mu
 if (present(gammai)) gammai = gam

end function get_spsound

!-----------------------------------------------------------------------
!+
!  Wrapper function to calculate temperature
!+
!-----------------------------------------------------------------------
real function get_temperature(eos_type,xyzi,rhoi,vxyzui,gammai,mui,Xi,Zi)
 integer, intent(in)             :: eos_type
 real,    intent(in)             :: xyzi(:),rhoi
 real,    intent(in)             :: vxyzui(:)
 real,    intent(in),   optional :: Xi,Zi
 real,    intent(inout),optional :: gammai,mui
 real                            :: tempi,gam,mu,X,Z

 !set defaults for variables not passed in
 mu    = gmw
 X     = X_in
 Z     = Z_in
 tempi = -1.  ! needed because temperature is an in/out to some equations of state, -ve == use own guess
 gam   = -1.  ! to indicate gamma is not being passed in
 if (present(Xi))  X  = Xi
 if (present(Zi))  Z  = Zi
 if (present(gammai)) gam = gammai
 if (present(mui)) mu = mui

 call get_TempPresCs(eos_type,xyzi,vxyzui,rhoi,tempi,gammai=gam,mui=mu,Xi=X,Zi=Z)

 get_temperature = tempi

 if (present(mui))    mui = mu
 if (present(gammai)) gammai = gam

end function get_temperature


!-----------------------------------------------------------------------
!+
!  Wrapper function to calculate temperature
!+
!-----------------------------------------------------------------------
real function get_temperature_from_u(eos_type,xpi,ypi,zpi,rhoi,ui,gammai,mui,Xi,Zi)
 integer, intent(in)             :: eos_type
 real,    intent(in)             :: xpi,ypi,zpi,rhoi
 real,    intent(in)             :: ui
 real,    intent(in),   optional :: Xi,Zi
 real,    intent(inout),optional :: gammai,mui
 real                            :: tempi,gam,mu,X,Z
 real :: vxyzui(4),xyzi(3)

 !set defaults for variables not passed in
 mu    = gmw
 X     = X_in
 Z     = Z_in
 tempi = -1.  ! needed because temperature is an in/out to some equations of state, -ve == use own guess
 gam   = -1.  ! to indicate gamma is not being passed in
 if (present(Xi))  X  = Xi
 if (present(Zi))  Z  = Zi
 if (present(gammai)) gam = gammai
 if (present(mui)) mu = mui

 vxyzui = (/0.,0.,0.,ui/)
 xyzi = (/xpi,ypi,zpi/)
 call get_TempPresCs(eos_type,xyzi,vxyzui,rhoi,tempi,gammai=gam,mui=mu,Xi=X,Zi=Z)

 get_temperature_from_u = tempi

 if (present(mui))    mui = mu
 if (present(gammai)) gammai = gam


end function get_temperature_from_u
!-----------------------------------------------------------------------
!+
!  Wrapper function to calculate pressure
!+
!-----------------------------------------------------------------------
real function get_pressure(eos_type,xyzi,rhoi,vxyzui,gammai,mui,Xi,Zi)
 integer, intent(in)             :: eos_type
 real,    intent(in)             :: xyzi(:),rhoi,vxyzui(:)
 real,    intent(in),   optional :: Xi,Zi
 real,    intent(inout),optional :: gammai,mui
 real                            :: presi,tempi,gam,mu,X,Z

 !set defaults for variables not passed in
 mu    = gmw
 X     = X_in
 Z     = Z_in
 tempi = -1.  ! needed because temperature is an in/out to some equations of state, -ve == use own guess
 gam   = -1.  ! to indicate gamma is not being passed in
 if (present(mui)) mu = mui
 if (present(Xi))  X  = Xi
 if (present(Zi))  Z  = Zi
 if (present(gammai)) gam = gammai
 if (present(mui)) mu = mui

 call get_TempPresCs(eos_type,xyzi,vxyzui,rhoi,tempi,presi=presi,gammai=gam,mui=mu,Xi=X,Zi=Z)

 get_pressure = presi

 if (present(mui))    mui = mu
 if (present(gammai)) gammai = gam

end function get_pressure

!-----------------------------------------------------------------------
!+
!  query function to return the internal energy for calculations with a
!  local mean molecular weight and local adiabatic index
!+
!-----------------------------------------------------------------------
real function get_local_u_internal(gammai, gmwi, gas_temp_local)
 real, intent(in) :: gammai, gmwi, gas_temp_local
 real             :: ponrhoi

 ponrhoi              = gas_temp_local/(gmwi*temperature_coef)
 get_local_u_internal = ponrhoi/(gammai-1.)

end function get_local_u_internal

!-----------------------------------------------------------------------
!+
!  get u from rho, T
!+
!-----------------------------------------------------------------------
real function get_u_from_rhoT(rho,temp,eos_type,uguess) result(u)
 use eos_mesa, only:get_eos_u_from_rhoT_mesa
 integer, intent(in)        :: eos_type
 real, intent(in)           :: rho,temp
 real, intent(in), optional :: uguess

 select case (eos_type)
 case(10) ! MESA EoS
    if (present(uguess)) then
       call get_eos_u_from_rhoT_mesa(rho,temp,u,uguess)
    else
       call get_eos_u_from_rhoT_mesa(rho,temp,u)
    endif

 case default
    u = temp/(gmw*temperature_coef*(gamma-1.))
 end select

end function get_u_from_rhoT

!-----------------------------------------------------------------------
!+
!  Get recombination energy (per unit mass) assumming complete
!  ionisation
!+
!-----------------------------------------------------------------------
subroutine calc_rec_ene(XX,YY,e_rec)
 real, intent(in)  :: XX, YY
 real, intent(out) :: e_rec
 real              :: e_H2,e_HI,e_HeI,e_HeII
 real, parameter   :: e_ion_H2   = 1.312e13, & ! ionisation energies in erg/mol
                      e_ion_HI   = 4.36e12, &
                      e_ion_HeI  = 2.3723e13, &
                      e_ion_HeII = 5.2505e13

 ! XX     : Hydrogen mass fraction
 ! YY     : Helium mass fraction
 ! e_rec  : Total ionisation energy due to H2, HI, HeI, and HeII

 e_H2   = 0.5 * XX * e_ion_H2
 e_HI   = XX * e_ion_HI
 e_HeI  = 0.25 * YY * e_ion_HeI
 e_HeII = 0.25 * YY * e_ion_HeII
 e_rec  = e_H2 + e_HI + e_HeI + e_HeII

end subroutine calc_rec_ene

!-----------------------------------------------------------------------
!+
!  Calculate temperature and specific internal energy from
!  pressure and density. Inputs and outputs are in cgs units.
!
!  Note on composition:
!  For ieos=2, 5, 12 and 17, mu_local is an input, X & Z are not used
!  For ieos=10, mu_local is not used
!  For ieos=20, mu_local is not used but available as an output
!+
!-----------------------------------------------------------------------
subroutine calc_temp_and_ene(eos_type,rho,pres,ene,temp,ierr,guesseint,mu_local,X_local,Z_local)
 use physcon,          only:Rg
 use eos_idealplusrad, only:get_idealgasplusrad_tempfrompres,get_idealplusrad_enfromtemp
 use eos_mesa,         only:get_eos_eT_from_rhop_mesa
 use eos_gasradrec,    only:calc_uT_from_rhoP_gasradrec
 use eos_stamatellos,  only:getintenerg_opdep
 integer, intent(in)              :: eos_type
 real,    intent(in)              :: rho,pres
 real,    intent(inout)           :: ene,temp
 real,    intent(in),    optional :: guesseint,X_local,Z_local
 real,    intent(inout), optional :: mu_local
 integer, intent(out)             :: ierr
 real                             :: mu,X,Z

 ierr = 0
 mu   = gmw
 X    = X_in
 Z    = Z_in
 if (present(mu_local)) mu = mu_local
 if (present(X_local))  X  = X_local
 if (present(Z_local))  Z  = Z_local
 select case(eos_type)
 case(2,5,17) ! Ideal gas
    temp = pres / (rho * Rg) * mu
    ene = pres / ( (gamma-1.) * rho)
 case(12) ! Ideal gas + radiation
    call get_idealgasplusrad_tempfrompres(pres,rho,mu,temp)
    call get_idealplusrad_enfromtemp(rho,temp,mu,ene)
 case(10) ! MESA EoS
    call get_eos_eT_from_rhop_mesa(rho,pres,ene,temp,guesseint)
 case(20) ! Ideal gas + radiation + recombination (from HORMONE, Hirai et al., 2020)
    call calc_uT_from_rhoP_gasradrec(rho,pres,X,1.-X-Z,temp,ene,mu,ierr)
    if (present(mu_local)) mu_local = mu
 case(24) ! Stamatellos
    temp = pres /(rho * Rg) * mu
    call getintenerg_opdep(temp, rho, ene)
 case default
    ierr = 1
 end select

end subroutine calc_temp_and_ene

!-----------------------------------------------------------------------
!+
!  Calculate density from pressure and temperature. Inputs and outputs
!  are in cgs units.
!
!  Note on composition:
!  For ieos=2, 5, 12 and 17, mu_local is an input, X & Z are not used
!  For ieos=10, mu_local is not used
!  For ieos=20, mu_local is not used but available as an output
!+
!-----------------------------------------------------------------------
subroutine calc_rho_from_PT(eos_type,pres,temp,rho,ierr,mu_local,X_local,Z_local)
 use physcon,          only:Rg
 use eos_idealplusrad, only:get_idealplusrad_rhofrompresT
 use eos_mesa,         only:get_eos_eT_from_rhop_mesa
 use eos_gasradrec,    only:calc_uT_from_rhoP_gasradrec
 integer, intent(in)              :: eos_type
 real,    intent(in)              :: pres,temp
 real,    intent(inout)           :: rho
 real,    intent(in),    optional :: X_local,Z_local
 real,    intent(inout), optional :: mu_local
 integer, intent(out)             :: ierr
 real                             :: mu,X,Z

 ierr = 0
 mu   = gmw
 X    = X_in
 Z    = Z_in
 if (present(mu_local)) mu = mu_local
 if (present(X_local))  X  = X_local
 if (present(Z_local))  Z  = Z_local
 select case(eos_type)
 case(2) ! Ideal gas
    rho = pres / (temp * Rg) * mu
 case(12) ! Ideal gas + radiation
    call get_idealplusrad_rhofrompresT(pres,temp,mu,rho)
 case default
    ierr = 1
 end select

end subroutine calc_rho_from_PT

!-----------------------------------------------------------------------
!+
!  Calculates specific entropy (gas + radiation + recombination)
!  up to an additive integration constant, from density and pressure.
!+
!-----------------------------------------------------------------------
function entropy(rho,pres,mu_in,ientropy,eint_in,ierr,T_in,Trad_in)
 use io,                only:fatal,warning
 use physcon,           only:radconst,kb_on_mh,Rg
 use eos_idealplusrad,  only:get_idealgasplusrad_tempfrompres
 use eos_mesa,          only:get_eos_eT_from_rhop_mesa
 use mesa_microphysics, only:getvalue_mesa
 real,    intent(in)            :: rho,pres,mu_in
 real,    intent(in),  optional :: eint_in,T_in,Trad_in
 integer, intent(in)            :: ientropy
 integer, intent(out), optional :: ierr
 real                           :: mu,entropy,logentropy,temp,Trad,eint

 if (present(ierr)) ierr=0

 mu = mu_in
 if (present(T_in)) then  ! is gas temperature specified?
    temp = T_in
 else
    temp = pres * mu / (rho * Rg)  ! used as initial guess for case 2
 endif

 select case(ientropy)
 case(1) ! Only include gas contribution
    ! check temp
    if (temp < tiny(0.)) call warning('entropy','temperature = 0 will give minus infinity with s entropy')
    entropy = kb_on_mh / mu * log(temp**1.5/rho)

 case(2) ! Include both gas and radiation contributions (up to additive constants)
    if (present(T_in)) then
       temp = T_in
       if (present(Trad_in)) then
          Trad = Trad_in
       else
          Trad = temp
       endif
    else
       call get_idealgasplusrad_tempfrompres(pres,rho,mu,temp) ! First solve for temp from rho and pres
       Trad = temp
    endif
    ! check temp
    if (temp < tiny(0.)) call warning('entropy','temperature = 0 will give minus infinity with s entropy')
    entropy = kb_on_mh / mu * log(temp**1.5/rho) + 4.*radconst*Trad**3 / (3.*rho)

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

real function get_entropy(rho,pres,mu_in,ieos)
 use units,   only:unit_density,unit_pressure,unit_ergg
 use physcon, only:kboltz
 integer, intent(in) :: ieos
 real, intent(in)    :: rho,pres,mu_in
 real                :: cgsrho,cgspres,cgss

 cgsrho = rho * unit_density
 cgspres = pres * unit_pressure
 select case (ieos)
 case (12)
    cgss = entropy(cgsrho,cgspres,mu_in,2)
 case (10, 20)
    cgss = entropy(cgsrho,cgspres,mu_in,3)
 case default
    cgss = entropy(cgsrho,cgspres,mu_in,1)
 end select
 cgss = cgss/kboltz ! s/kb
 get_entropy = cgss/unit_ergg

end function get_entropy

!-----------------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using Newton-Raphson
!  method
!+
!-----------------------------------------------------------------------
subroutine get_rho_from_p_s(pres,S,rho,mu,rhoguess,ientropy)
 real, intent(in)    :: pres,S,mu,rhoguess
 real, intent(inout) :: rho
 real                :: srho,srho_plus_dsrho,S_plus_dS,dSdsrho
 real(kind=8)        :: corr
 real,    parameter  :: eoserr=1e-9,dfac=1e-12
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
!  Calculate temperature given density and entropy using Newton-Raphson
!  method
!+
!-----------------------------------------------------------------------
subroutine get_p_from_rho_s(ieos,S,rho,mu,P,temp)
 use physcon, only:radconst,Rg,mass_proton_cgs,kboltz
 use io,      only:fatal
 use eos_idealplusrad, only:get_idealgasplusrad_tempfrompres,get_idealplusrad_pres
 use units,   only:unit_density,unit_pressure,unit_ergg
 real, intent(in)    :: S,mu,rho
 real, intent(inout) :: temp
 real, intent(out)   :: P
 integer, intent(in) :: ieos
 real                :: corr,df,f,temp_new,cgsrho,cgsp,cgss
 real,    parameter  :: eoserr=1e-12
 integer             :: niter
 integer, parameter  :: nitermax = 1000

 ! change to cgs unit
 cgsrho = rho*unit_density
 cgss   = s*unit_ergg

 niter = 0
 select case (ieos)
 case (2,5,17)
    temp = (cgsrho * exp(mu*cgss*mass_proton_cgs))**(2./3.)
    cgsP = cgsrho*Rg*temp / mu
 case (12)
    corr = huge(corr)
    do while (abs(corr) > eoserr .and. niter < nitermax)
       f = 1. / (mu*mass_proton_cgs) * log(temp**1.5/cgsrho) + 4.*radconst*temp**3 / (3.*cgsrho*kboltz) - cgss
       df = 1.5 / (mu*temp*mass_proton_cgs) + 4.*radconst*temp**2 / (cgsrho*kboltz)
       corr = f/df
       temp_new = temp - corr
       if (temp_new > 1.2 * temp) then
          temp = 1.2 * temp
       elseif (temp_new < 0.8 * temp) then
          temp = 0.8 * temp
       else
          temp = temp_new
       endif
       niter = niter + 1
    enddo
    call get_idealplusrad_pres(cgsrho,temp,mu,cgsP)
 case default
    cgsP = 0.
    call fatal('eos','[get_p_from_rho_s] only implemented for eos 2 and 12')
 end select

 ! check temp
 if (temp > huge(0.)) call fatal('entropy','entropy too large gives infinite temperature:'// &
                                 ' suggest to reduce C_ent for one dump')

 ! change back to code unit
 P = cgsP / unit_pressure

end subroutine get_p_from_rho_s

!-----------------------------------------------------------------------
!+
!  Calculate mean molecular weight from X and Z, assuming complete
!  ionisation
!+
!-----------------------------------------------------------------------
real function get_mean_molecular_weight(XX,ZZ) result(mu)
 real, intent(in) :: XX,ZZ
 real :: YY

 YY = 1.-XX-ZZ
 mu = 1./(2.*XX + 0.75*YY + 0.5*ZZ)

end function get_mean_molecular_weight

!---------------------------------------------------------
!+
!  return cv from rho, u in code units
!+
!---------------------------------------------------------
real function get_cv(rho,u,cv_type) result(cv)
 use mesa_microphysics, only:getvalue_mesa
 use units,             only:unit_ergg,unit_density
 use physcon,           only:Rg
 real, intent(in)    :: rho,u
 integer, intent(in) :: cv_type
 real                :: rho_cgs,u_cgs,temp

 select case (cv_type)

 case(1)  ! MESA EoS
    rho_cgs = rho*unit_density
    u_cgs = u*unit_ergg
    call getvalue_mesa(rho_cgs,u_cgs,4,temp)
    cv = u_cgs/temp / unit_ergg
 case default  ! constant cv
    cv = Rg/((gamma-1.)*gmw*unit_ergg)
 end select

end function get_cv

!-----------------------------------------------------------------------
!+
!  subroutine sets polyk based on utherm/positions
!  read from an sphNG dump file
!+
!-----------------------------------------------------------------------
subroutine setpolyk(eos_type,iprint,utherm,xyzhi,npart)
 use part, only:xyzmh_ptmass
 use io,   only:id,master
 integer, intent(in) :: eos_type,iprint
 real,    intent(in) :: utherm(:)
 real,    intent(in) :: xyzhi(:,:)
 integer, intent(in) :: npart
 integer :: ipart
 real    :: r2,polykalt

 !-- pick a random particle from which to extract polyk
 ipart = npart/2

 select case(eos_type)
 case(1,8)
!
!--isothermal eos
!
    polykalt = 2./3.*utherm(ipart)
    !--check all other utherms identical
    if (any(utherm(1:npart) /= utherm(ipart)) .and. id==master) then
       write(iprint,*) 'WARNING! different utherms but run is isothermal'
    endif

 case(2,5,17)
!
!--adiabatic/polytropic eos
!  this routine is ONLY called if utherm is NOT stored, so polyk matters
!
    if (id==master) write(iprint,*) 'Using polytropic equation of state, gamma = ',gamma
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
    if (id==master) write(iprint,*) ' WARNING! unknown equation of state in setpolyk'
    polykalt = polyk

 end select

 if (diff(polykalt,polyk) .and. id==master) then
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
 elseif (polyk < tiny(polyk) .and. id==master) then
    write(iprint,*) 'WARNING: polyk = 0 in equation of state'
 endif

end subroutine setpolyk
!-----------------------------------------------------------------------
!+
!  small utility returns whether two real numbers differ
!+
!-----------------------------------------------------------------------
logical pure function diff(r1,r2)
 real, intent(in) :: r1,r2

 diff = abs(r1-r2) > tiny(r1)

end function diff

!-----------------------------------------------------------------------
!+
!  Query function to return whether an EoS is non-ideal
!  Mainly used to decide whether it is necessary to write
!  things like pressure and temperature in the dump file or not
!+
!-----------------------------------------------------------------------
logical function eos_is_non_ideal(ieos)
 integer, intent(in) :: ieos

 select case(ieos)
 case(10,12,15,20)
    eos_is_non_ideal = .true.
 case default
    eos_is_non_ideal = .false.
 end select

end function eos_is_non_ideal

!-----------------------------------------------------------------------
!+
!  Query function to return whether an EoS outputs mean molecular weight
!+
!-----------------------------------------------------------------------
logical function eos_outputs_mu(ieos)
 integer, intent(in) :: ieos

 select case(ieos)
 case(20)
    eos_outputs_mu = .true.
 case(24)
    eos_outputs_mu = .true.
 case default
    eos_outputs_mu = .false.
 end select

end function eos_outputs_mu

!-----------------------------------------------------------------------
!+
!  Query function to whether to print pressure to dump file
!+
!-----------------------------------------------------------------------
logical function eos_outputs_gasP(ieos)
 integer, intent(in) :: ieos

 select case(ieos)
 case(8,9,10,15,23)
    eos_outputs_gasP = .true.
 case default
    eos_outputs_gasP = .false.
 end select

end function eos_outputs_gasP

!-----------------------------------------------------------------------
!+
!  Query function for whether the equation of state requires
!  the code to be compiled without a thermal energy variable
!+
!-----------------------------------------------------------------------
logical function eos_requires_isothermal(ieos)
 integer, intent(in) :: ieos

 select case(ieos)
 case(1,3,6,7,8,13,14,21)
    eos_requires_isothermal = .true.
 case default
    !case(2,5,4,10,11,12,15,16,17,20,22,23,24,9)
    eos_requires_isothermal = .false.
 end select

end function eos_requires_isothermal

!-----------------------------------------------------------------------
!+
!  Query function for whether the equation of state requires
!  the polyk variable to be set
!+
!-----------------------------------------------------------------------
logical function eos_requires_polyk(ieos)
 integer, intent(in) :: ieos

 eos_requires_polyk = eos_requires_isothermal(ieos) .or. ieos==9

end function eos_requires_polyk

!-----------------------------------------------------------------------
!+
!  function to skip unimplemented eos choices
!+
!-----------------------------------------------------------------------
logical function eos_is_not_implemented(ieos)
 integer, intent(in) :: ieos

 select case(ieos)
 case(18,19)
    eos_is_not_implemented = .true.
 case default
    eos_is_not_implemented = .false.
 end select

end function eos_is_not_implemented

!-----------------------------------------------------------------------
!+
!  prints equation of state info in the run header
!+
!-----------------------------------------------------------------------
subroutine eosinfo(eos_type,iprint)
 use dim,            only:maxvxyzu
 use io,             only:fatal,id,master
 use eos_helmholtz,  only:eos_helmholtz_eosinfo
 use eos_barotropic, only:eos_info_barotropic
 use eos_piecewise,  only:eos_info_piecewise
 use eos_gasradrec,  only:eos_info_gasradrec
 use eos_stamatellos, only:eos_file
 integer, intent(in) :: eos_type,iprint

 if (id/=master) return

 select case(eos_type)
 case(1,11)
    if (1.0d-5 < polyk .and. polyk < 1.0d3) then
       write(iprint,"(/,a,f10.6)")  ' Isothermal equation of state:     cs^2 = ',polyk
    else
       write(iprint,"(/,a,Es13.6)") ' Isothermal equation of state:     cs^2 = ',polyk
    endif
    if (eos_type==11) write(iprint,*) ' (ZERO PRESSURE) '
 case(2)
    if (maxvxyzu >= 4) then
       write(iprint,"(/,a,f10.6,a,f10.6)") ' Adiabatic equation of state: P = (gamma-1)*rho*u, gamma = ',&
                                              gamma,' gmw = ',gmw
    else
       write(iprint,"(/,a,f10.6,a,f10.6,a,f10.6)") ' Polytropic equation of state: P = ',polyk,'*rho^',gamma,' gmw = ',gmw
    endif
 case(3)
    write(iprint,"(/,a,f10.6,a,f10.6)") ' Locally isothermal eq of state (R_sph): cs^2_0 = ',polyk,' qfac = ',qfacdisc
 case(5,17)
    if (maxvxyzu >= 4) then
       write(iprint,"(' Adiabatic equation of state: P = (gamma-1)*rho*u, where gamma & mu depend on the formation of H2')")
    else
       write(iprint,*) 'ERROR: eos = 5,17 cannot assume isothermal conditions'
    endif
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

 case(24)
    write(iprint,"(/,a,a)") 'Using tabulated Eos from file:', eos_file, 'and calculated gamma.'
 end select
 write(iprint,*)

end subroutine eosinfo

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
    call add_to_iheader(istrat,'istrat',hdr,ierr)
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
 use io,                only:iprint,id,master,iverbose
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
         if (iverbose >= 0) write(iprint,*) 'KROME eos: initial gamma = 1.666667'
       else
          if (iverbose >= 0) write(iprint,*) 'adiabatic eos: gamma = ',gamma
       endif
    else
      if (iverbose >= 0) write(iprint,*) 'setting isothermal sound speed^2 (polyk) = ',polyk,' gamma = ',gamma
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
       if (id==master) write(iprint,*) 'ERROR: qfacdisc <= 0'
       ierr = 2
    else
       if (id==master .and. iverbose >= 0) write(iprint,*) 'qfacdisc = ',qfacdisc
    endif
 endif

 if (ieos==7) then
    call extract('istrat',istrat,hdr,ierr)
    call extract('alpha_z',alpha_z,hdr,ierr)
    call extract('beta_z', beta_z, hdr,ierr)
    call extract('z0',z0,hdr,ierr)
    if (abs(qfacdisc2) <= tiny(qfacdisc2)) then
       if (id==master) write(iprint,*) 'ERROR: qfacdisc2 == 0'
       ierr = 2
    else
       if (id==master .and. iverbose >= 0) write(iprint,*) 'qfacdisc2 = ',qfacdisc2
    endif
 endif

end subroutine read_headeropts_eos

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos(iunit)
 use dim,            only:use_krome
 use infile_utils,   only:write_inopt
 use eos_helmholtz,  only:eos_helmholtz_write_inopt
 use eos_barotropic, only:write_options_eos_barotropic
 use eos_piecewise,  only:write_options_eos_piecewise
 use eos_gasradrec,  only:write_options_eos_gasradrec
 use eos_tillotson,  only:write_options_eos_tillotson
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling equation of state'
 call write_inopt(ieos,'ieos','eqn of state (1=isoth;2=adiab;3=locally iso;8=barotropic)',iunit)

 if (.not.use_krome .or. .not.eos_outputs_mu(ieos)) then
    call write_inopt(gmw,'mu','mean molecular weight',iunit)
 endif

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
 case(23)
    call write_options_eos_tillotson(iunit)
 end select

end subroutine write_options_eos

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_eos(name,valstring,imatch,igotall,ierr)
 use dim,            only:store_dust_temperature,update_muGamma
 use io,             only:fatal
 use eos_barotropic, only:read_options_eos_barotropic
 use eos_piecewise,  only:read_options_eos_piecewise
 use eos_gasradrec,  only:read_options_eos_gasradrec
 use eos_tillotson,  only:read_options_eos_tillotson
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot  = 0
 character(len=30), parameter  :: label = 'read_options_eos'
 logical :: igotall_barotropic,igotall_piecewise,igotall_gasradrec
 logical :: igotall_tillotson

 imatch  = .true.
 igotall_barotropic = .true.
 igotall_piecewise  = .true.
 igotall_gasradrec =  .true.
 igotall_tillotson =  .true.

 select case(trim(name))
 case('ieos')
    read(valstring,*,iostat=ierr) ieos
    ngot = ngot + 1
    if (ieos <= 0 .or. ieos > maxeos) call fatal(label,'equation of state choice out of range')
    if (ieos == 5 .or. ieos == 17) then
       store_dust_temperature = .true.
       update_muGamma = .true.
    endif
 case('mu')
    read(valstring,*,iostat=ierr) gmw
    ! not compulsory to read in
    if (gmw <= 0.)  call fatal(label,'mu <= 0')
 case('X')
    read(valstring,*,iostat=ierr) X_in
    if (X_in <= 0. .or. X_in >= 1.) call fatal(label,'X must be between 0 and 1')
    ngot = ngot + 1
 case('Z')
    read(valstring,*,iostat=ierr) Z_in
    if (Z_in <= 0. .or. Z_in > 1.) call fatal(label,'Z must be between 0 and 1')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if (.not.imatch .and. ieos== 8) call read_options_eos_barotropic(name,valstring,imatch,igotall_barotropic,ierr)
 if (.not.imatch .and. ieos== 9) call read_options_eos_piecewise( name,valstring,imatch,igotall_piecewise, ierr)
 if (.not.imatch .and. ieos==20) call read_options_eos_gasradrec( name,valstring,imatch,igotall_gasradrec, ierr)
 if (.not.imatch .and. ieos==23) call read_options_eos_tillotson(name,valstring,imatch,igotall_tillotson,ierr)

 !--make sure we have got all compulsory options (otherwise, rewrite input file)
 igotall = (ngot >= 1) .and. igotall_piecewise .and. igotall_barotropic .and. igotall_gasradrec &
                       .and. igotall_tillotson

end subroutine read_options_eos


!-----------------------------------------------------------------------

end module eos
