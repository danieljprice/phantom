!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_potekhin
!
! Equation of state for fully ionized electron-ion plasmas (EOS EIP)
!   by A. Y. Potekhin & G. Chabrier (2013), with 2022 updates (eos22.f / MELANGE9).
!   Evaluated on-the-fly with smooth ideal-gas blending at low densities.
!
! :References:
!   Potekhin & Chabrier (2013), A&A, 550, A43
!   Baiko & Chugunov (2022), MNRAS, 510, 2628
!
! :Owner: Marat
!
! :Runtime parameters:
!   - xc  : *Carbon mass fraction*
!   - xh  : *Hydrogen mass fraction*
!   - xhe : *Helium mass fraction*
!   - xmg : *Magnesium mass fraction*
!   - xne : *Neon mass fraction*
!   - xo  : *Oxygen mass fraction*
!
! :Dependencies: infile_utils, io, physcon, units
!
 implicit none

 public :: eos_potekhin_init
 public :: read_options_eos_potekhin
 public :: write_options_eos_potekhin
 public :: eos_potekhin_pres_sound
 public :: eos_potekhin_compute_pres_sound
 public :: eos_potekhin_energy_from_rhoT
 public :: eos_potekhin_cv_dpresdt
 public :: eos_potekhin_get_minrho
 public :: eos_potekhin_get_maxrho
 public :: eos_potekhin_get_mintemp
 public :: eos_potekhin_get_maxtemp
 public :: eos_potekhin_eosinfo
 public :: eos_potekhin_finish

 private

 ! Mass fractions set by user at runtime (consistent with eos_helmholtz)
 real :: xh  = 0.0
 real :: xhe = 0.0
 real :: xc  = 0.5
 real :: xo  = 0.5
 real :: xne = 0.0
 real :: xmg = 0.0

 integer, parameter :: speciesmax = 15
 character(len=10) :: speciesname(speciesmax)
 real :: xmass(speciesmax) ! mass fraction of species
 real :: Aion(speciesmax)  ! number of nucleons
 real :: Zion(speciesmax)  ! number of protons
 real :: abar, zbar, mu_ideal

 ! Precomputed mixture arrays for MELANGE9 (sorted in ascending Z)
 integer :: nmix_pot = 0
 real    :: ay_pot(speciesmax), azion_pot(speciesmax), acmi_pot(speciesmax)

 ! Physical / operational validity limits (CGS & K)
 real, parameter :: rhomincgs = 1.0e-12
 real, parameter :: rhomaxcgs = 1.0e14
 real, parameter :: tempmin   = 1.0e3
 real, parameter :: tempmax   = 1.0e11

 ! Low-density blending thresholds (g/cm^3):
 ! Below rho_blend_lo: pure ideal gas + radiation
 ! Above rho_blend_hi: pure Potekhin MELANGE9 + radiation
 real, parameter :: rho_blend_lo = 10.0
 real, parameter :: rho_blend_hi = 100.0

 ! Conversion factor for Potekhin temperature in atomic units: T[au] = T6 / UN_T6
 real, parameter :: un_t6_pot = 0.315775

contains

!----------------------------------------------------------------
!+
!  Initialise Potekhin EOS composition and warm up legacy routines
!+
!----------------------------------------------------------------
subroutine eos_potekhin_init(ierr)
 use io, only:warning,id,master
 integer, intent(out) :: ierr
 integer :: i
 real    :: dummy_pres, dummy_sound, dummy_ener, dummy_denerdt

 ierr = 0

 speciesname(1)  = "hydrogen"
 speciesname(2)  = "helium"
 speciesname(3)  = "carbon"
 speciesname(4)  = "oxygen"
 speciesname(5)  = "neon"
 speciesname(6)  = "magnesium"
 speciesname(7)  = "silicon"
 speciesname(8)  = "sulphur"
 speciesname(9)  = "argon"
 speciesname(10) = "calcium"
 speciesname(11) = "titanium"
 speciesname(12) = "chromium"
 speciesname(13) = "iron"
 speciesname(14) = "nickel"
 speciesname(15) = "zinc"

 Aion(1)  = 1.0  ; Zion(1)  = 1.0
 Aion(2)  = 4.0  ; Zion(2)  = 2.0
 Aion(3)  = 12.0 ; Zion(3)  = 6.0
 Aion(4)  = 16.0 ; Zion(4)  = 8.0
 Aion(5)  = 20.0 ; Zion(5)  = 10.0
 Aion(6)  = 24.0 ; Zion(6)  = 12.0
 Aion(7)  = 28.0 ; Zion(7)  = 14.0
 Aion(8)  = 32.0 ; Zion(8)  = 16.0
 Aion(9)  = 36.0 ; Zion(9)  = 18.0
 Aion(10) = 40.0 ; Zion(10) = 20.0
 Aion(11) = 44.0 ; Zion(11) = 22.0
 Aion(12) = 48.0 ; Zion(12) = 24.0
 Aion(13) = 52.0 ; Zion(13) = 26.0
 Aion(14) = 56.0 ; Zion(14) = 28.0
 Aion(15) = 60.0 ; Zion(15) = 30.0

 xmass(:) = 0.0
 xmass(1) = xh
 xmass(2) = xhe
 xmass(3) = xc
 xmass(4) = xo
 xmass(5) = xne
 xmass(6) = xmg

 if (sum(xmass(:)) > 1.0 + 1.0e-6 .or. sum(xmass(:)) < 1.0 - 1.0e-6) then
    if (id == master) call warning('eos_potekhin','mass fractions total != 1')
    ierr = 1
    return
 endif

 abar = 1.0 / sum(xmass(:) / Aion(:))
 zbar = abar * sum(xmass(:) * Zion(:) / Aion(:))
 mu_ideal = 1.0 / sum(xmass(:) * (1.0 + Zion(:)) / Aion(:))

 ! Build active species list (Aion/Zion are already sorted in ascending Z)
 nmix_pot = 0
 ay_pot(:) = 0.0
 azion_pot(:) = 0.0
 acmi_pot(:) = 0.0
 do i = 1, speciesmax
    if (xmass(i) > 1.0e-12) then
       nmix_pot = nmix_pot + 1
       azion_pot(nmix_pot) = Zion(i)
       acmi_pot(nmix_pot)  = Aion(i)
       ay_pot(nmix_pot)    = (xmass(i) / Aion(i)) * abar
    endif
 enddo

 if (nmix_pot == 0) then
    if (id == master) call warning('eos_potekhin','no active chemical species')
    ierr = 1
    return
 endif

 ! Single-threaded warm-up calls so legacy SAVE initialization inside eos22.f
 ! (e.g. KINIT in ANHBC and KRUN in BILN9) runs once before any OpenMP region.
 call eos_potekhin_compute_pres_sound(1.0e7, 1.0e6, dummy_pres, dummy_sound, dummy_ener, dummy_denerdt)
 call eos_potekhin_compute_pres_sound(1.0e5, 1.0e6, dummy_pres, dummy_sound, dummy_ener, dummy_denerdt)

end subroutine eos_potekhin_init

!----------------------------------------------------------------
!+
!  Clean up (nothing dynamically allocated for on-the-fly evaluation)
!+
!----------------------------------------------------------------
subroutine eos_potekhin_finish(ierr)
 integer, intent(out) :: ierr
 ierr = 0
end subroutine eos_potekhin_finish

!----------------------------------------------------------------
!+
!  Write options to the input file
!+
!----------------------------------------------------------------
subroutine write_options_eos_potekhin(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(xh ,'xh' ,'Hydrogen mass fraction',iunit)
 call write_inopt(xhe,'xhe','Helium mass fraction',iunit)
 call write_inopt(xc ,'xc' ,'Carbon mass fraction',iunit)
 call write_inopt(xo ,'xo' ,'Oxygen mass fraction',iunit)
 call write_inopt(xne,'xne','Neon mass fraction',iunit)
 call write_inopt(xmg,'xmg','Magnesium mass fraction',iunit)

end subroutine write_options_eos_potekhin

!----------------------------------------------------------------
!+
!  Read options from the input file
!+
!----------------------------------------------------------------
subroutine read_options_eos_potekhin(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(xh ,'xh' ,db,errcount=nerr,min=0.,max=1.,default=xh)
 call read_inopt(xhe,'xhe',db,errcount=nerr,min=0.,max=1.,default=xhe)
 call read_inopt(xc ,'xc' ,db,errcount=nerr,min=0.,max=1.,default=xc)
 call read_inopt(xo ,'xo' ,db,errcount=nerr,min=0.,max=1.,default=xo)
 call read_inopt(xne,'xne',db,errcount=nerr,min=0.,max=1.,default=xne)
 call read_inopt(xmg,'xmg',db,errcount=nerr,min=0.,max=1.,default=xmg)

end subroutine read_options_eos_potekhin

real function eos_potekhin_get_minrho()
 use units, only:unit_density
 eos_potekhin_get_minrho = real(rhomincgs / unit_density)
end function eos_potekhin_get_minrho

real function eos_potekhin_get_maxrho()
 use units, only:unit_density
 eos_potekhin_get_maxrho = real(rhomaxcgs / unit_density)
end function eos_potekhin_get_maxrho

real function eos_potekhin_get_mintemp()
 eos_potekhin_get_mintemp = tempmin
end function eos_potekhin_get_mintemp

real function eos_potekhin_get_maxtemp()
 eos_potekhin_get_maxtemp = tempmax
end function eos_potekhin_get_maxtemp

!----------------------------------------------------------------
!+
!  Print EOS information
!+
!----------------------------------------------------------------
subroutine eos_potekhin_eosinfo(iprint)
 integer, intent(in) :: iprint
 integer :: i

 write(iprint,"(/,a)") ' Potekhin & Chabrier (2013, 2022) EIP equation of state (on-the-fly)'
 write(iprint,"(a)")   '   mass fractions of each species:'
 do i = 1, speciesmax
    if (xmass(i) > 0.0) then
       write(iprint,"(a,a,a,f5.3)") '     ', speciesname(i), ': ', xmass(i)
    endif
 enddo

end subroutine eos_potekhin_eosinfo

!----------------------------------------------------------------
!+
!  Called by eos.f90: Newton-Raphson inversion to get T, P/rho, cs from (rho, u)
!  Inputs/outputs are in Phantom code units (except tempi in K).
!+
!----------------------------------------------------------------
subroutine eos_potekhin_pres_sound(tempi,rhoi,ponrhoi,spsoundi,eni)
 use units, only:unit_density,unit_pressure,unit_ergg,unit_velocity
 use io,    only:warning
 real, intent(inout) :: tempi
 real, intent(in)    :: rhoi
 real, intent(out)   :: ponrhoi
 real, intent(out)   :: spsoundi
 real, intent(in)    :: eni
 integer, parameter  :: maxiter = 30
 real,    parameter  :: tol = 1.0e-4
 logical :: done
 integer :: itercount
 real :: cgsrhoi, cgspresi, cgsspsoundi, cgseni, cgseni_eos, cgsdendti
 real :: tnew, tprev, dt_step

 cgsrhoi = rhoi * unit_density
 cgseni  = eni * unit_ergg

 if (tempi <= 0.0) tempi = 1.0e6
 if (tempi < tempmin) tempi = tempmin
 if (tempi > tempmax) tempi = tempmax

 call eos_potekhin_compute_pres_sound(tempi, cgsrhoi, cgspresi, cgsspsoundi, cgseni_eos, cgsdendti)

 tprev = tempi
 dt_step = (cgseni_eos - cgseni) / max(cgsdendti, 1.0e-30)
 tnew = tempi - dt_step
 if (tnew > 2.0 * tempi) tnew = 2.0 * tempi
 if (tnew < 0.5 * tempi) tnew = 0.5 * tempi
 if (tnew > tempmax) tnew = tempmax
 if (tnew < tempmin) tnew = tempmin

 itercount = 0
 done = (abs(tnew - tprev) < tempi * tol)

 do while (.not. done)
    itercount = itercount + 1
    tprev = tnew
    call eos_potekhin_compute_pres_sound(tnew, cgsrhoi, cgspresi, cgsspsoundi, cgseni_eos, cgsdendti)

    dt_step = (cgseni_eos - cgseni) / max(cgsdendti, 1.0e-30)
    tnew = tnew - dt_step

    if (tnew > 2.0 * tprev) tnew = 2.0 * tprev
    if (tnew < 0.5 * tprev) tnew = 0.5 * tprev

    if (abs(tnew - tprev) < max(tprev, tempmin) * tol) then
       done = .true.
    endif

    if (tnew >= tempmax) then
       tnew = tempmax
       done = .true.
    elseif (tnew <= tempmin) then
       tnew = tempmin
       done = .true.
    endif

    if (itercount >= maxiter .and. .not. done) then
       call warning('eos','Potekhin eos failed to converge on temperature')
       done = .true.
    endif
 enddo

 ! Recompute final thermodynamics at converged tnew if it was clamped or updated
 if (abs(tnew - tprev) > 0.0) then
    call eos_potekhin_compute_pres_sound(tnew, cgsrhoi, cgspresi, cgsspsoundi, cgseni_eos, cgsdendti)
 endif

 tempi    = tnew
 ponrhoi  = cgspresi / (unit_pressure * rhoi)
 spsoundi = cgsspsoundi / unit_velocity

end subroutine eos_potekhin_pres_sound

!----------------------------------------------------------------
!+
!  Core thermodynamic calculation in CGS units from (temp [K], den [g/cm^3]).
!+
!----------------------------------------------------------------
subroutine eos_potekhin_compute_pres_sound(temp,den,pres,sound,ener,denerdt)
 use physcon, only:kboltz,c,steboltz,avogadro
 use io,      only:fatal
 real, intent(in)  :: temp, den
 real, intent(out) :: pres, sound, ener, denerdt

 real, parameter :: asol   = 4.0 * steboltz / c
 real, parameter :: asoli3 = asol / 3.0
 real, parameter :: light2 = c * c

 real :: t_safe, d_safe, prad, erad, dpraddt, deraddt
 real :: p_id, u_id, cv_id, cs_id
 real :: p_pot, u_pot, cv_pot, cs_pot
 real :: x, w

 if (temp <= 0.0) call fatal('eos_potekhin','temperature <= 0')
 if (den  <= 0.0) call fatal('eos_potekhin','density <= 0')

 t_safe = min(max(temp, tempmin), tempmax)
 d_safe = min(max(den, rhomincgs), rhomaxcgs)

 ! Radiation contribution (common to both regimes)
 prad    = asoli3 * t_safe**4
 dpraddt = 4.0 * prad / t_safe
 erad    = 3.0 * prad / d_safe
 deraddt = 3.0 * dpraddt / d_safe

 if (d_safe <= rho_blend_lo) then
    call compute_ideal_plus_rad(t_safe, d_safe, prad, erad, deraddt, pres, sound, ener, denerdt)
 elseif (d_safe >= rho_blend_hi) then
    call compute_potekhin_plus_rad(t_safe, d_safe, prad, erad, deraddt, pres, sound, ener, denerdt)
 else
    call compute_ideal_plus_rad(t_safe, d_safe, prad, erad, deraddt, p_id, cs_id, u_id, cv_id)
    call compute_potekhin_plus_rad(t_safe, d_safe, prad, erad, deraddt, p_pot, cs_pot, u_pot, cv_pot)

    x = (log10(d_safe) - log10(rho_blend_lo)) / (log10(rho_blend_hi) - log10(rho_blend_lo))
    x = min(max(x, 0.0), 1.0)
    w = x * x * (3.0 - 2.0 * x)

    pres    = (1.0 - w) * p_id   + w * p_pot
    ener    = (1.0 - w) * u_id   + w * u_pot
    denerdt = (1.0 - w) * cv_id  + w * cv_pot
    sound   = (1.0 - w) * cs_id  + w * cs_pot
 endif

end subroutine eos_potekhin_compute_pres_sound

!----------------------------------------------------------------
!+
!  Ideal fully-ionized gas + Planck radiation (low-density regime)
!+
!----------------------------------------------------------------
subroutine compute_ideal_plus_rad(temp,den,prad,erad,deraddt,pres,sound,ener,denerdt)
 use physcon, only:kboltz,avogadro,c
 real, intent(in)  :: temp, den, prad, erad, deraddt
 real, intent(out) :: pres, sound, ener, denerdt

 real :: pgas, egas, cvgas, chit, chid, gam1, z
 real, parameter :: light2 = c * c

 pgas  = den * (avogadro / mu_ideal) * kboltz * temp
 egas  = 1.5 * pgas / den
 cvgas = 1.5 * (avogadro / mu_ideal) * kboltz

 pres    = pgas + prad
 ener    = egas + erad
 denerdt = cvgas + deraddt

 chid  = pgas / pres
 chit  = (pgas + 4.0 * prad) / pres
 gam1  = chid + (pres / (den * temp * denerdt)) * chit * chit
 z     = 1.0 + (ener + light2) * (den / pres)
 sound = c * sqrt(max(gam1 / z, 1.0e-30))

end subroutine compute_ideal_plus_rad

!----------------------------------------------------------------
!+
!  Potekhin MELANGE9 plasma + Planck radiation (high-density regime)
!+
!----------------------------------------------------------------
subroutine compute_potekhin_plus_rad(temp,den,prad,erad,deraddt,pres,sound,ener,denerdt)
 use physcon, only:kboltz,avogadro,c
 real, intent(in)  :: temp, den, prad, erad, deraddt
 real, intent(out) :: pres, sound, ener, denerdt

 integer :: liqsol, nmix_loc
 real(kind=8) :: ay_loc(speciesmax), az_loc(speciesmax), ac_loc(speciesmax)
 real(kind=8) :: rho_d, temp_au, pradnkt, dens_au, zmean, cmimean, z2mean
 real(kind=8) :: gamimean, chi, tpt, pnkt, unkt, snk, cv_nk, chir, chit_pl
 real :: nkt_vol, nkt_mass, p_plasma, u_plasma, cv_plasma
 real :: chit_eff, chir_eff, gam1, z
 real, parameter :: light2 = c * c

 nmix_loc = nmix_pot
 ay_loc(1:nmix_loc) = ay_pot(1:nmix_loc)
 az_loc(1:nmix_loc) = azion_pot(1:nmix_loc)
 ac_loc(1:nmix_loc) = acmi_pot(1:nmix_loc)

 rho_d   = real(den, kind=8)
 temp_au = real((temp * 1.0e-6) / un_t6_pot, kind=8)

 ! Thread-safety lock: eos22.f uses explicit SAVE in subroutines.
 !$omp critical(eos_potekhin_melange)
 call melange9(nmix_loc, ay_loc, az_loc, ac_loc, rho_d, temp_au, &
               pradnkt, dens_au, zmean, cmimean, z2mean, gamimean, chi, tpt, liqsol, &
               pnkt, unkt, snk, cv_nk, chir, chit_pl)
 !$omp end critical(eos_potekhin_melange)

 nkt_mass = (avogadro / real(cmimean)) * kboltz * temp
 nkt_vol  = den * nkt_mass

 p_plasma  = real(pnkt) * nkt_vol
 u_plasma  = real(unkt) * nkt_mass
 cv_plasma = real(cv_nk) * (avogadro / real(cmimean)) * kboltz

 pres    = max(p_plasma + prad, 1.0e-30)
 ener    = u_plasma + erad
 denerdt = max(cv_plasma + deraddt, 1.0e-30)

 chir_eff = (p_plasma / pres) * real(chir)
 chit_eff = (p_plasma / pres) * real(chit_pl) + 4.0 * (prad / pres)

 gam1  = chir_eff + (pres / (den * temp * denerdt)) * chit_eff * chit_eff
 gam1  = max(gam1, 1.0e-6)
 z     = 1.0 + (ener + light2) * (den / pres)
 sound = c * sqrt(gam1 / z)

end subroutine compute_potekhin_plus_rad

!----------------------------------------------------------------
!+
!  Get specific internal energy [erg/g] from rho [g/cm^3] and T [K]
!+
!----------------------------------------------------------------
subroutine eos_potekhin_energy_from_rhoT(rho,temp,ener)
 real, intent(in)  :: rho, temp
 real, intent(out) :: ener
 real :: dummy_pres, dummy_sound, dummy_denerdt

 call eos_potekhin_compute_pres_sound(temp, rho, dummy_pres, dummy_sound, ener, dummy_denerdt)

end subroutine eos_potekhin_energy_from_rhoT

!----------------------------------------------------------------
!+
!  Get c_v [erg/(g*K)] and dP/dT [dyn/(cm^2*K)] from rho [g/cm^3] and T [K]
!+
!----------------------------------------------------------------
subroutine eos_potekhin_cv_dpresdt(temp,den,cv,dpresdt)
 real, intent(in)  :: temp, den
 real, intent(out) :: cv, dpresdt
 real :: pres, sound, ener, t1, t2, p1, p2, s_dum, e_dum, cv_dum

 call eos_potekhin_compute_pres_sound(temp, den, pres, sound, ener, cv)

 ! Numerical dP/dT across small relative step
 t1 = temp * 0.999
 t2 = temp * 1.001
 call eos_potekhin_compute_pres_sound(t1, den, p1, s_dum, e_dum, cv_dum)
 call eos_potekhin_compute_pres_sound(t2, den, p2, s_dum, e_dum, cv_dum)
 dpresdt = (p2 - p1) / (t2 - t1)

end subroutine eos_potekhin_cv_dpresdt

end module eos_potekhin
