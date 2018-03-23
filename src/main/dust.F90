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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    K_code            -- drag constant when constant drag is used
!    graindens         -- Intrinsic grain density in g/cm^3
!    grainsize         -- Initial grain size in cm
!    icut_backreaction -- cut the drag on the gas phase (0=no, 1=yes)
!    idrag             -- gas/dust drag (0=off,1=Epstein/Stokes,2=const K,3=const ts)
!
!  DEPENDENCIES: dim, eos, infile_utils, io, physcon, units
!+
!--------------------------------------------------------------------------

module dust
 implicit none
 !--Default values for the dust in the infile
 real, public     :: K_code          = 1.
 real, public     :: grainsizecgs    = 0.1
 real, public     :: graindenscgs    = 3.

 integer, public  :: idrag             = 1
 integer, public  :: icut_backreaction = 0
 real, public     :: grainsize,graindens
 public           :: write_options_dust,read_options_dust
 public           :: get_ts, init_drag
 public           :: print_dustinfo

 ! generic interface to set_dustfrac
 interface set_dustfrac
  module procedure set_dustfrac_single, set_dustfrac_power_law
 end interface set_dustfrac
 public           :: set_dustfrac

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
 use units,    only:udist,umass,unit_density
 use physcon,  only:mass_proton_cgs,cross_section_H2_cgs
 use eos,      only:gamma
 integer, intent(out) :: ierr
 real :: cste_seff
 real :: mass_mol_gas, cross_section_gas

 ierr = 0
 !--compute constants which are used in the ts calculation
 if (gamma < 1.) then
    call error('init_drag','gamma < 1',var='gamma',val=gamma)
    ierr = 1
 endif
 cste_mu           = sqrt(2./(pi*gamma))
 coeff_gei_1       = sqrt(8./(pi*gamma))

 !--compute the grain mass (spherical compact grains of radius s)
 !--change this line for fractal grains or grains of different sizes
 grainsize         = grainsizecgs/udist
 graindens         = graindenscgs/unit_density
 if (grainsize <= 0.) then
    call error('init_drag','grain size <= 0',var='grainsize',val=grainsize)
    ierr = 2
 endif
 if (graindens <= 0.) then
    call error('init_drag','grain density <= 0',var='graindens',val=graindens)
    ierr = 2
 endif
 !--compute the effective surface density used to calculate the mean free path
 cste_seff         = pi/sqrt(2.)*5./64.
 mass_mol_gas      = (2.*mass_proton_cgs)/umass
 cross_section_gas = cross_section_H2_cgs/(udist*udist)
 seff              = cste_seff*mass_mol_gas/cross_section_gas
 if (seff <= 0.) then
    call error('init_drag','effective surface density <= 0',var='seff',val=seff)
    ierr = 3
 endif
 !--check the value of K_code
 if ((idrag == 2 .or. idrag == 3) .and. K_code < 0.) then
    call error('init_drag','K_code < 0',var='K_code',val=K_code)
    ierr = 4
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
 real :: rhocrit, grainmass

 select case(idrag)
 case(1)
    if (use_dustgrowth) then
       write(iprint,"(a)")              ' Using Epstein/Stokes drag with variable grain size. '
    else
       grainmass = 4./3.*pi*graindens*grainsize**3
       write(iprint,"(a)")              ' Using Epstein/Stokes drag with constant grain size: '
       write(iprint,"(2(a,1pg10.3),a)") '        Grain size = ',grainsize*udist,        ' cm     = ',grainsize,' (code units)'
       write(iprint,"(2(a,1pg10.3),a)") '        Grain mass = ',grainmass*umass,        ' g      = ',grainmass,' (code units)'
       write(iprint,"(2(a,1pg10.3),a)") '     Grain density = ',graindens*unit_density, ' g/cm^3 = ',graindens,' (code units)'
       write(iprint,"(2(a,1pg10.3),a)") '  Gas mfp at rho=1 = ',seff*udist/unit_density,' cm     = ',seff,' (code units)'
       rhocrit = 9.*seff/(4.*grainsize)
       write(iprint,"(/,a)")              ' Density above which Stokes drag is used:'
       write(iprint,"(2(a,1pg10.3),a)") '           rhocrit = ',rhocrit*unit_density,   ' g/cm^3 = ',rhocrit,' (code units)'
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
 real, intent(out) :: dustfrac

 dustfrac = dust_to_gas/(1. + dust_to_gas)

end subroutine set_dustfrac_single

!----------------------------------------------------------------
!+
!  utility function to set the dust fraction given the
!  dust-to-gas ratio. Equation (57) in Price & Laibe (2015)
!+
!----------------------------------------------------------------
subroutine set_dustfrac_power_law(dust_to_gas_tot,dustfrac,s,smin,smax,sindex)
 use io, only: fatal
 real, intent(in)  :: dust_to_gas_tot, s, smin, smax, sindex
 real, intent(out) :: dustfrac
 real :: dust_to_gas_ind

 !dn/ds=const*(s/smax)^-sindex  where dn is the number of particles per cm^3 per interval ds in size (Draine et al. 2006)
 if (sindex /= 4.) then
    dust_to_gas_ind = (dust_to_gas_tot*(4.-sindex)*s**(3.-sindex))/(smax**(4.-sindex)-smin**(4.-sindex))
 elseif (sindex == 4.) then
    dust_to_gas_ind = (dust_to_gas_tot*s**(3.-sindex))/(log(smax/smin))
 else
    call fatal('dust','congratulations on discovering a new number (please check your value for sindex)')
 endif
 dustfrac = dust_to_gas_ind/(1. + dust_to_gas_ind)

end subroutine set_dustfrac_power_law

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
    if (ts1 >= 0.) then
       ts  = 1./ts1
    else
       ts = huge(ts)
    endif

 case(2)
    !
    ! constant drag coefficient
    !
    if (K_code > 0.) then
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

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_dust(iunit)
 use infile_utils, only:write_inopt
 use dim, only:use_dustgrowth
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling dust'
 call write_inopt(idrag,'idrag','gas/dust drag (0=off,1=Epstein/Stokes,2=const K,3=const ts)',iunit)
 if (.not.use_dustgrowth) then
    call write_inopt(grainsizecgs,'grainsize','Grain size in cm',iunit)
 else
    call write_inopt(grainsizecgs,'grainsize','Initial grain size in cm',iunit)
 endif
 call write_inopt(graindenscgs,'graindens','Intrinsic grain density in g/cm^3',iunit)
 call write_inopt(K_code,'K_code','drag constant when constant drag is used',iunit)
 call write_inopt(icut_backreaction,'icut_backreaction','cut the drag on the gas phase (0=no, 1=yes)',iunit)

end subroutine write_options_dust

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_dust(name,valstring,imatch,igotall,ierr)
 use units, only:udist,umass
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 real(kind=8)  :: udens

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('idrag')
    read(valstring,*,iostat=ierr) idrag
    ngot = ngot + 1
 case('grainsize')
    read(valstring,*,iostat=ierr) grainsizecgs
    grainsize = grainsizecgs/udist
    ngot = ngot + 1
 case('graindens')
    read(valstring,*,iostat=ierr) graindenscgs
    udens = umass/udist**3
    graindens = graindenscgs/udens
 case('K_code')
    read(valstring,*,iostat=ierr) K_code
    ngot = ngot + 1
 case('icut_backreaction')
    read(valstring,*,iostat=ierr) icut_backreaction
 case default
    imatch = .false.
 end select

 if (idrag == 2 .or. idrag==3) then
    if (ngot >= 3) igotall = .true.
 else
    if (ngot >= 2) igotall = .true.
 endif

end subroutine read_options_dust

end module dust
