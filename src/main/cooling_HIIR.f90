!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_HIIR
!
! Two phase thermal model of HII regions with an adiabatic eos (icooling=10)
!
!  the gas is either ionised (Tion) or cold (Tcold), the ionised fraction phi
!  of each particle being read from mu, set by the HII routines. Each phase
!  has a target T_target = phi*Tion + (1-phi)*Tcold:
!
!   - fully ionised gas is held at uIon by the stars
!   - neutral gas away from the ionisation front is held at Tcold, as in
!     ieos=21 (in force, its hydro du/dt is ignored)
!   - the gas of the front (partially ionised, or neutral with an ionised
!     neighbour) is heated at once to its target by the stars and, above it,
!     cools over t_relax = Crelax*h/c_s(u): the heat conducted from the
!     ionised gas and of the shocks spreads the front over a few h
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters:
!   - Crelax : *relaxation time of the gas of the front in units of h/c_s*
!
! :Dependencies: eos, eos_HIIR, infile_utils
!
 implicit none

 real, public :: Crelax = 1.

 public :: init_cooling_HIIR,relax_HIIR,heated_u_HIIR,held_HIIR,is_ionised
 public :: write_options_cooling_HIIR,read_options_cooling_HIIR

 private

 real :: mucold,tempcoef

contains

!-----------------------------------------------------------------------
!+
!  initialise the two phases, from the cold gas of the eos
!+
!-----------------------------------------------------------------------
subroutine init_cooling_HIIR(ierr)
 use eos,      only:gamma,polyk,gmw,temperature_coef
 use eos_HIIR, only:init_eos_HIIR
 integer, intent(out) :: ierr

 call init_eos_HIIR(gamma,polyk,gmw,temperature_coef,ierr)
 mucold   = gmw
 tempcoef = temperature_coef

end subroutine init_cooling_HIIR

!-----------------------------------------------------------------------
!+
!  thermal energy update in the step:
!    u_new = u_target + (u - u_target)*exp(-dt/t_relax)
!  exact for any dt; uIon for the fully ionised gas and u_target at once
!  below the target. The neutral gas held at Tcold has no du/dt (force),
!  so it stays there
!+
!-----------------------------------------------------------------------
subroutine relax_HIIR(ui,hi,mui,gammai,dt,dudti)
 use eos_HIIR, only:uIon
 real, intent(in)  :: ui,hi,mui,gammai,dt
 real, intent(out) :: dudti
 real :: utarget,csi,trelax,unew

 dudti = 0.
 if (dt <= 0. .or. gammai <= 1.) return

 if (is_fully_ionised(mui)) then
    unew = uIon
 else
    utarget = get_utarget(mui,gammai)
    if (ui <= utarget) then
       unew = utarget  ! heated by the stars
    else
       csi    = sqrt(gammai*(gammai - 1.)*ui)
       trelax = Crelax*hi/csi
       if (trelax > 0.) then
          unew = utarget + (ui - utarget)*exp(-dt/trelax)
       else
          unew = utarget  ! Crelax = 0: two isothermal phases, as ieos=21
       endif
    endif
 endif
 dudti = (unew - ui)/dt

end subroutine relax_HIIR

!-----------------------------------------------------------------------
!+
!  thermal energy given to the forces: uIon for the fully ionised gas,
!  at least the target otherwise. Used in cons2prim on the predicted
!  thermal energy of derivs, so the forces see the gas just ionised
!+
!-----------------------------------------------------------------------
real function heated_u_HIIR(ui,mui,gammai) result(uheated)
 use eos_HIIR, only:uIon
 real, intent(in) :: ui,mui,gammai

 if (is_fully_ionised(mui)) then
    uheated = uIon
 elseif (gammai > 1.) then
    uheated = max(ui,get_utarget(mui,gammai))
 else
    uheated = ui
 endif

end function heated_u_HIIR

!-----------------------------------------------------------------------
!+
!  phases held at their temperature during the step (their hydro du/dt is
!  ignored in force): the fully ionised gas, and the neutral gas with no
!  ionised neighbour
!+
!-----------------------------------------------------------------------
logical function held_HIIR(mui,nionised_neigh)
 real, intent(in) :: mui,nionised_neigh

 if (is_fully_ionised(mui)) then
    held_HIIR = .true.
 else
    held_HIIR = (.not.is_ionised(mui) .and. nionised_neigh <= 0.)
 endif

end function held_HIIR

!-----------------------------------------------------------------------
!+
!  gas (partially) ionised by the HII routines; mu not set yet is neutral
!+
!-----------------------------------------------------------------------
logical function is_ionised(mui)
 real, intent(in) :: mui

 is_ionised = (mui > 0. .and. mui < mucold)

end function is_ionised

logical function is_fully_ionised(mui)
 use eos_HIIR, only:muion
 real, intent(in) :: mui

 is_fully_ionised = (mui > 0. .and. mui <= muion)

end function is_fully_ionised

!-----------------------------------------------------------------------
!+
!  thermal energy of the phase, T_target = phi*Tion + (1-phi)*Tcold
!+
!-----------------------------------------------------------------------
real function get_utarget(mui,gammai) result(utarget)
 use eos_HIIR, only:Tion,Tcold,muion
 real, intent(in) :: mui,gammai
 real :: mu,phi

 mu = mui
 if (mu <= 0.) mu = mucold
 phi = min(max((mucold - mu)/(mucold - muion),0.),1.)
 utarget = (phi*Tion + (1. - phi)*Tcold)/(tempcoef*mu*(gammai - 1.))

end function get_utarget

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling_HIIR(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(Crelax,'Crelax','relaxation time of the gas of the HII fronts in units of h/c_s',iunit)

end subroutine write_options_cooling_HIIR

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling_HIIR(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(Crelax,'Crelax',db,errcount=nerr,min=0.,default=Crelax)

end subroutine read_options_cooling_HIIR

end module cooling_HIIR
