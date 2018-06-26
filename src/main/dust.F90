!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: dust
!
!  DESCRIPTION:
!  Contains routine for gas-dust drag term
!
!  REFERENCES:
!    Laibe & Price (2012a,b)
!    Kwok (1975), Draine et al. (2006)
!
!  OWNER: Mark Hutchison
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, eos, io, physcon, units
!+
!--------------------------------------------------------------------------

module dust
 use dim, only: ndusttypes
 implicit none
 !--Default values for the dust in the infile
 real, public     :: K_code                   = 1.
 real, public     :: grainsizecgs(ndusttypes) = 0.1
 real, public     :: smincgs                  = 1.e-5
 real, public     :: smaxcgs                  = 0.1
 real, public     :: sindex                   = 3.5
 real, public     :: graindenscgs(ndusttypes) = 3.

 integer, public  :: idrag             = 1
 integer, public  :: icut_backreaction = 0
 logical, public  :: ilimitdustflux    = .false. ! to limit spurious dust generation in outer disc
 real, public     :: grainsize(ndusttypes),graindens(ndusttypes)
 real, private    :: grainmass(ndusttypes)
 public           :: get_ts
 public           :: init_drag
 public           :: print_dustinfo
 ! generic interface to set_dustfrac
 interface set_dustfrac
  module procedure set_dustfrac_single, set_dustfrac_power_law
 end interface set_dustfrac
 public           :: set_dustfrac
 public           :: set_grainsize

 real, private    :: cste_mu,coeff_gei_1,seff
 private

contains

!-----------------------------------------------------------------------
!+
!  initialize the drag: compute the quantities that are used once
!+
!-----------------------------------------------------------------------
subroutine init_drag(ierr)
 use physcon,  only:pi
 use io,       only:error
 use units,    only:udist,umass
 use physcon,  only:mass_proton_cgs,cross_section_H2_cgs
 use eos,      only:gamma
 integer, intent(out) :: ierr
 integer :: i
 real    :: cste_seff
 real    :: mass_mol_gas, cross_section_gas

 ierr = 0
 !--compute constants which are used in the ts calculation
 if (gamma < 1.) then
    call error('init_drag','gamma < 1',var='gamma',val=gamma)
    ierr = 1
 endif
 cste_mu           = sqrt(2./(pi*gamma))
 coeff_gei_1       = sqrt(8./(pi*gamma))

 select case(idrag)
 case(1)
    !--compute the grain mass (spherical compact grains of radius s)
    call set_grainsize(smincgs,smaxcgs)
    do i = 1,ndusttypes
       if (grainmass(i) <= 0. .and. idrag == 1) then
          call error('init_drag','grain size/density <= 0',var='grainmass',val=grainmass(i))
          ierr = 2
       endif
       if (grainsize(i) <= 0.) then
          call error('init_drag','grain size <= 0',var='grainsize',val=grainsize(i))
          ierr = 2
       endif
       if (graindens(i) <= 0.) then
          call error('init_drag','grain density <= 0',var='graindens',val=graindens(i))
          ierr = 2
       endif
    enddo
 case(2,3)
    !--check the value of K_code
    if (K_code < 0.) then
       call error('init_drag','K_code < 0',var='K_code',val=K_code)
       ierr = 4
    endif
 case default
 end select

 !--compute the effective surface density used to calculate the mean free path
 cste_seff         = pi/sqrt(2.)*5./64.
 mass_mol_gas      = (2.*mass_proton_cgs)/umass
 cross_section_gas = cross_section_H2_cgs/(udist*udist)
 seff              = cste_seff*mass_mol_gas/cross_section_gas
 if (seff <= 0.) then
    call error('init_drag','effective surface density <= 0',var='seff',val=seff)
    ierr = 3
 endif

end subroutine init_drag

!--------------------------------------------
!+
!  print information about the dust physics
!+
!--------------------------------------------
subroutine print_dustinfo(iprint)
 use units, only:unit_density,umass,udist
 use physcon,  only:pi
 use dim, only:use_dustgrowth
 integer, intent(in) :: iprint
 integer :: i
 real    :: rhocrit

 select case(idrag)
 case(1)
    if (use_dustgrowth) then
       write(iprint,"(a)") ' Using Epstein/Stokes drag with variable grain size. '
    else
       write(iprint,"(a)") ' Using Epstein/Stokes drag with constant grain size: '
       do i = 1,ndusttypes
          write(iprint,"(2(a,1pg10.3),a)") '        Grain size = ',grainsize(i)*udist,      &
                                           ' cm     = ',grainsize(i),' (code units)'
          write(iprint,"(2(a,1pg10.3),a)") '        Grain mass = ',grainmass(i)*umass,      &
                                           ' g      = ',grainmass(i),' (code units)'
          write(iprint,"(2(a,1pg10.3),a)") '     Grain density = ',graindens(i)*unit_density,  &
                                           ' g/cm^3 = ',graindens(i),' (code units)'
          write(iprint,"(2(a,1pg10.3),a)") '  Gas mfp at rho=1 = ',seff*udist/unit_density, &
                                           ' cm     = ',seff,' (code units)'
          rhocrit = 9.*seff/(4.*grainsize(i))
       enddo
       write(iprint,"(/,a)") ' Density above which Stokes drag is used:'
       write(iprint,"(2(a,1pg10.3),a)")    '           rhocrit = ',rhocrit*unit_density,    &
                                           ' g/cm^3 = ',rhocrit,' (code units)'
    endif
 case(2)
    write(iprint,"(/,a,1pg12.5)") ' Using K=const drag with K = ',K_code
 case(3)
    write(iprint,"(/,a,1pg12.5)") ' Using ts=const drag with ts = ',K_code
 case default
    write(iprint,"(/,a)") ' Drag regime not set'
 end select

end subroutine print_dustinfo


!----------------------------------------------------------------
!+
!  utility function to set the dust fraction given the
!  dust-to-gas ratio. Equation (57) in Price & Laibe (2015)
!+
!----------------------------------------------------------------
subroutine set_dustfrac_single(dust_to_gas,dustfrac)
 real, intent(in)  :: dust_to_gas
 real, intent(out) :: dustfrac(:)

 dustfrac = dust_to_gas/(1. + dust_to_gas)

end subroutine set_dustfrac_single

!----------------------------------------------------------------
!+
!  utility function to set the dust fraction given the
!  dust-to-gas ratio.
!+
!----------------------------------------------------------------
subroutine set_dustfrac_power_law(dust_to_gas_tot,dustfrac,smin,smax,sind)
 use io,  only: fatal
 real, intent(in)  :: dust_to_gas_tot,smin,smax,sind
 real, intent(out) :: dustfrac(:)
 integer :: i
 real :: dustfrac_tot
 real :: norm
 real :: rhodtot
 real :: grid(ndusttypes+1) = 0.
 real :: rhodusti(ndusttypes)
 real :: exact
 real :: power = 0.
 real, parameter :: tol = 1.e-10

 !--reset global power-law index
 sindex = sind

 if (smax==smin .or. ndusttypes==1) then
    !--If all the same grain size, then just scale the dust fraction
    dustfrac = dust_to_gas_tot/(1.+dust_to_gas_tot)*1./real(ndusttypes)
 else
    call set_grainsize(smin,smax,grid)

    !--Dust density is computed from drhodust ∝ dn*mdust where dn ∝ s**(-p)*ds
    !  and mdust ∝ s**(3). This is then integrated across each cell to account
    !  for mass contributions from unrepresented grain sizes
    do i = 1,ndusttypes
       if (sindex == 4.) then
          rhodusti(i) = log(grid(i+1)/grid(i))
       else
          power = 4. - sindex
          rhodusti(i) = 1./power*(grid(i+1)**power - grid(i)**power)
       endif
    enddo

    !--Sum the contributions from each cell to get total relative dust content
    rhodtot = sum(rhodusti)

    !--Calculate the total dust fraction from the dust-to-gas ratio
    dustfrac_tot = dust_to_gas_tot/(1.+dust_to_gas_tot)

    !--Calculate the normalisation factor (∝ 1/rhotot) and scale the dust fractions
    !  Note: dust density and dust fraction have the same power-law dependence on s.
    norm         = dustfrac_tot/rhodtot
    dustfrac(:)  = norm*rhodusti(:)

    !--Check to make sure the integral determining the contributions is correct
    if (sindex == 4.) then
       exact = log(grid(ndusttypes+1)/grid(1))
    else
       exact = 1./power*(grid(ndusttypes+1)**power - grid(1)**power)
    endif
    if (abs(rhodtot-exact)/exact>tol) &
       call fatal('dust','Piecewise integration of MRN distribution not matching the exact solution!')
 endif

end subroutine set_dustfrac_power_law

!-----------------------------------------------------------------------------
!+
!  utility function to set the grain size
!
!  if spread of sizes, optionally returns an array of size bins in log space
!+
!-----------------------------------------------------------------------------
subroutine set_grainsize(smin,smax,grid)
 use physcon, only:pi
 use io,      only: fatal
 use units,   only: udist,unit_density
 real, intent(in)  :: smin,smax
 real, optional, intent(out) :: grid(:)
 integer :: i
 real :: log_ds
 real :: log_grid(ndusttypes+1)

 smincgs = smin
 smaxcgs = smax

 ! check whether grid is passed in, and if so, that is large enough
 if (present(grid)) then
    if (size(grid) < ndusttypes+1) then
       print *, 'error trying to pass grid of insufficient size to set_grainsize()'
    endif
 endif

 if (ndusttypes==1) then
    !--Grain size is set in the input file
 elseif (smax==smin .and. ndusttypes>1) then
    !--If all the same grain size, then just scale the dust fraction
    grainsizecgs(:) = smax
 else
    !--Create a uniform grid with N+1 points between smax and smin (inclusive)
    log_ds = log10(smax/smin)/real(ndusttypes)
    do i = 1,ndusttypes+1
       log_grid(i) = log10(smin) + (i-1)*log_ds
    enddo

    !--Convert grid coordinates back to real space
    log_grid = 10.**log_grid

    !--Find representative s for each cell
    !  (skewed towards small grains because there are more small grains than large grains)
    do i = 1,ndusttypes
       grainsizecgs(i) = sqrt(log_grid(i)*log_grid(i+1))
    enddo

    ! if we supplied grid, then return it
    if (present(grid)) then
       grid = log_grid ! this is no longer log of the grid at this point
    endif
 endif

 !--Set the grain properties relating to grain size
 grainsize(:) = grainsizecgs(:)/udist
 graindens(:) = graindenscgs(:)/unit_density
 grainmass(:) = 4./3.*pi*graindens(:)*grainsize(:)**3

end subroutine set_grainsize


!----------------------------------------------------------------------------
!+
!  get the stopping time (rhoi*rhoj)/(K*(rhoi+rhoj)) for a pair of particles
!
!  idrag = 1 : Epstein/Stokes with automatic switching
!  idrag = 2 : const K
!+
!----------------------------------------------------------------------------
subroutine get_ts(idrag,sgrain,densgrain,rhogas,rhodust,spsoundgas,dv2, &
                  ts,iregime)
 use physcon,     only:pi
 integer, intent(in)  :: idrag
 integer, intent(out) :: iregime
 real,    intent(in)  :: sgrain,densgrain,rhogas,rhodust,spsoundgas,dv2
 real,    intent(out) :: ts

 real :: tol_super
 real :: rhosum,abs_dv,kwok
 real :: lambda,kn_eff,viscmol_nu,Re_dust
 real :: dragcoeff,f,ts1

 ! initialise variables
 tol_super  = 0.1
 dragcoeff  = 0.
 f          = 0.
 ts1        = 0.
 ts         = 0.
 rhosum     = rhogas + rhodust
 ! rhosum     = rhogas ! this is an approx. to allow calculations to proceed
 ! efficiently in case of numerical dust trapping
 ! should be rhosum =  rhogas + rhodust
 ! however, without full dust grain population rhodust is an
 ! approx. anyway

 ! compute quantities specific to the drag regime
 select case(idrag)
 case(1)
    !
    ! physical drag (Epstein or Stokes regime)
    ! check if the regime is Epstein or Stokes
    !
    lambda = seff/rhogas
    if (sgrain > 0.) then
       kn_eff = 9.*lambda/(4.*sgrain)
    else
       kn_eff = huge(kn_eff)
    endif

    if (kn_eff >= 1.) then
       !
       ! Epstein regime
       !
       if (densgrain > tiny(densgrain)) then
          dragcoeff = coeff_gei_1*spsoundgas/(densgrain*sgrain)
          !if (dragcoeff > 1.e10) print*,dragcoeff
       else
          dragcoeff = huge(dragcoeff) ! so get ts=0 in this case
       endif
       if (spsoundgas > 0. .and. dv2 > 0.) then
          kwok = 9.*pi/128.*dv2/(spsoundgas*spsoundgas)
          f = sqrt(1.+kwok)
       else
          kwok = 0.
          f = 1. ! not important
       endif
       iregime   = 1
       ! count where Kwok (1975) correction for supersonic drag is important
       if (kwok > tol_super) iregime = 2
    else
       !
       ! Stokes regime
       !
       viscmol_nu = cste_mu*lambda*spsoundgas  ! kinematic viscosity
       !--compute the local Stokes number
       abs_dv  = sqrt(dv2)
       Re_dust = 2.*sgrain*abs_dv/viscmol_nu
       if (Re_dust  <=  1.) then
          dragcoeff = 4.5*viscmol_nu/(densgrain*sgrain*sgrain)
          f         = 1.
          iregime   = 3
       elseif (Re_dust  <=  800.) then
          dragcoeff = 9./(densgrain*sgrain*Re_dust**0.6)
          f         = abs_dv
          iregime   = 4
       else
          dragcoeff = 0.163075/(densgrain*sgrain)  ! coeff is (3/8)*24/800**0.6
          f         = abs_dv
          iregime   = 5
       endif
    endif
    if (dragcoeff == huge(dragcoeff)) then
       ts1 = huge(ts1)
    else
       ts1 = dragcoeff*f*rhosum
    endif
    if (ts1 > 0.) then
       ts  = 1./ts1
    else
       ts = huge(ts)
    endif

 case(2)
    !
    ! constant drag coefficient
    !
    if (K_code > 0.) then
       ! WARNING! When ndusttypes > 1, K_code ONLY makes sense
       ! if all of the grains are identical to one another.
       ts = rhogas*rhodust/(K_code*rhosum)
    else
       ts = huge(ts)
    endif
    iregime = 0

 case(3)
    !
    ! constant ts
    !
    ts = K_code
    iregime = 0

 case default
    ts = 0.
    iregime = 0 ! unknown
 end select

end subroutine get_ts

end module dust
