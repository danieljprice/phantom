!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dust
!
! Contains routine for gas-dust drag term
!
! :References:
!    Laibe & Price (2012a,b)
!    Kwok (1975), Draine et al. (2006)
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - K_code            : *drag constant when constant drag is used*
!   - graindens         : *Intrinsic grain density in g/cm^3*
!   - grainsize         : *Grain size in cm*
!   - icut_backreaction : *cut the drag on the gas phase (0=no, 1=yes)*
!   - idrag             : *gas/dust drag (0=off,1=Epstein/Stokes,2=const K,3=const ts)*
!   - ilimitdustflux    : *limit the dust flux using Ballabio et al. (2018)*
!   - irecon            : *use reconstruction in gas/dust drag (-1=off,0=no slope limiter,1=van leer MC)*
!
! :Dependencies: dim, eos, fileutils, infile_utils, io, options, part,
!   physcon, units
!

 use dim,     only:use_dustgrowth,maxdusttypes
 use part,    only:ndusttypes,grainsize,graindens
 use physcon, only:pi
 use units,   only:umass,udist
 implicit none
 !--Default values for the dust in the infile
 real,    public  :: K_code(maxdusttypes) = 1.
 real,    public  :: grainsizecgs         = 0.1
 real,    public  :: graindenscgs         = 3.
 integer, public  :: idrag                = 1
 integer, public  :: icut_backreaction    = 0
 integer, public  :: irecon               = 1
 logical, public  :: ilimitdustflux       = .false. ! to limit spurious dust generation in outer disc
 logical, public  :: drag_implicit        = .false.  ! use implicit scheme for 2-fluids drag forces

 public :: get_ts
 public :: init_drag
 public :: print_dustinfo
 public :: write_options_dust
 public :: read_options_dust
 public :: get_viscmol_nu

 real, private :: cste_mu,coeff_gei_1,seff
 private

contains

!--------------------------------------------------------------------------
!+
!  initialize the drag: compute the quantities that are used once
!+
!--------------------------------------------------------------------------
subroutine init_drag(ierr)
 use eos,      only:gamma
 use io,       only:error
 use physcon,  only:mass_proton_cgs,cross_section_H2_cgs
 integer, intent(out) :: ierr
 real    :: cste_seff
 real    :: mass_mol_gas, cross_section_gas
 integer :: i

 ierr = 0
 !--compute constants which are used in the ts calculation
 if (gamma < 1.) then
    call error('init_drag','gamma < 1',var='gamma',val=gamma)
    ierr = 1
 endif
 cste_mu           = sqrt(2./(pi*gamma))
 coeff_gei_1       = sqrt(8./(pi*gamma))

 select case(idrag)
 case(2,3)
    !--check the value of K_code
    do i=1,maxdusttypes
       if (K_code(i) < 0.) then
          call error('init_drag','K_code < 0',var='K_code',val=K_code(i))
          ierr = 4
       endif
    enddo
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

!--------------------------------------------------------------------------
!+
!  print information about the dust physics
!+
!--------------------------------------------------------------------------
subroutine print_dustinfo(iprint)
 use dim,   only:maxdusttypes
 use units, only:unit_density
 integer, intent(in) :: iprint
 real    :: rhocrit,grainmass(maxdusttypes)
 integer :: i

 select case(idrag)
 case(1)
    if (use_dustgrowth) then
       write(iprint,"(a)") ' Using Epstein/Stokes drag with variable grain size. '
    else
       write(iprint,"(a)") ' Using Epstein/Stokes drag with constant grain size: '
       !--compute the grain mass (spherical compact grains of radius s)
       grainmass(:) = 4./3.*pi*graindens(:)*grainsize(:)**3
       do i=1,ndusttypes
          write(iprint,"(2(a,1pg10.3),a)") '        Grain size = ',grainsize(i)*udist,      &
                                           ' cm     = ',grainsize(i),' (code units)'
          write(iprint,"(2(a,1pg10.3),a)") '        Grain mass = ',grainmass(i)*umass,      &
                                           ' g      = ',grainmass(i),' (code units)'
          write(iprint,"(2(a,1pg10.3),a)") '     Grain density = ',graindens(i)*unit_density,  &
                                           ' g/cm^3 = ',graindens(i),' (code units)'
          write(iprint,"(2(a,1pg10.3),a)") '  Gas mfp at rho=1 = ',seff*udist/unit_density, &
                                           ' cm     = ',seff,' (code units)'
          rhocrit = 9.*seff/(4.*grainsize(i))
          write(iprint,"(/,a)") ' Density above which Stokes drag is used:'
          write(iprint,"(2(a,1pg10.3),a)")    '           rhocrit = ',rhocrit*unit_density,    &
                                              ' g/cm^3 = ',rhocrit,' (code units)'
       enddo
    endif
 case(2)
    do i=1,ndusttypes
       write(iprint,"(/,a,1pg12.5)") ' Using K=const drag with K = ',K_code(i)
    enddo
 case(3)
    do i=1,ndusttypes
       write(iprint,"(/,a,1pg12.5)") ' Using ts=const drag with ts = ',K_code(i)
    enddo
 case default
    write(iprint,"(/,a)") ' Drag regime not set'
 end select

end subroutine print_dustinfo

!--------------------------------------------------------------------------
!+
!  get the stopping time (rhoi*rhoj)/(K*(rhoi+rhoj)) for a pair of
!  particles
!
!  idrag = 1 : Epstein/Stokes with automatic switching
!  idrag = 2 : const K
!+
!--------------------------------------------------------------------------
subroutine get_ts(idrag,idust,sgrain,densgrain,rhogas,rhodust,spsoundgas,dv2, &
                  ts,iregime)
 integer, intent(in)  :: idrag,idust
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
    if (K_code(idust) > 0.) then
       ts = rhogas*rhodust/(K_code(idust)*rhosum)
    else
       ts = huge(ts)
    endif
    iregime = 0

 case(3)
    !
    ! constant ts
    !
    ts = K_code(idust)
    iregime = 0

 case default
    ts = 0.
    iregime = 0 ! unknown
 end select

end subroutine get_ts

!--------------------------------------------------------------------------
!+
!  writes input dust options to the input file
!  Note: ndustypes & use_dustfract are read from the dump file, so will
!  not be correctly printed in the header, where iunit=iprint
!+
!--------------------------------------------------------------------------
subroutine write_options_dust(iunit)
 use fileutils,    only:make_tags_unique
 use infile_utils, only:write_inopt
 use options,      only:use_dustfrac
 integer, intent(in) :: iunit
 character(len=10)   :: numdust
 character(len=20)   :: duststring(maxdusttypes)
 integer             :: i

 write(iunit,"(/,a)") '# options controlling dust'

 call write_inopt(idrag,'idrag','gas/dust drag (0=off,1=Epstein/Stokes,2=const K,3=const ts)',iunit)

 select case(idrag)
 case(1)
    if (ndusttypes <= 1) then
       if (use_dustgrowth) then
          call write_inopt(grainsizecgs,'grainsize','Initial grain size in cm',iunit)
       else
          call write_inopt(grainsizecgs,'grainsize','Grain size in cm',iunit)
       endif
       call write_inopt(graindenscgs,'graindens','Intrinsic grain density in g/cm^3',iunit)
    endif
 case(2,3)
    if (ndusttypes > 1) then
       write(numdust,'(i10)') ndusttypes
       duststring='K_code'
       call make_tags_unique(ndusttypes,duststring)
       do i=1,ndusttypes
          call write_inopt(K_code(i),duststring(i),'drag constant when constant drag is used',iunit)
       enddo
    else
       call write_inopt(K_code(1),'K_code','drag constant when constant drag is used',iunit)
    endif
 end select

 if (use_dustfrac) then
    call write_inopt(ilimitdustflux,'ilimitdustflux','limit the dust flux using Ballabio et al. (2018)',iunit)
 else
    call write_inopt(irecon,'irecon','use reconstruction in gas/dust drag (-1=off,0=no slope limiter,1=van leer MC)',iunit)
    call write_inopt(drag_implicit,'drag_implicit','gas/dust drag implicit scheme (works only with IND_TIMESTEPS=no)',iunit)
 endif

 call write_inopt(icut_backreaction,'icut_backreaction','cut the drag on the gas phase (0=no, 1=yes)',iunit)

end subroutine write_options_dust

!--------------------------------------------------------------------------
!+
!  reads input dust options from the input file
!+
!--------------------------------------------------------------------------
subroutine read_options_dust(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 real(kind=8)  :: udens
 integer, parameter :: nvalues = 5
 integer, parameter :: iidrag        = 1, &
                       ibackreact    = 2, &
                       igrainsize    = 3, &
                       igraindens    = 4, &
                       iKcode        = 5
 integer, save :: igot(nvalues) = 0
 integer       :: ineed(nvalues)
 integer       :: int
 character(len=10) :: str

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('idrag')
    read(valstring,*,iostat=ierr) idrag
    igot(iidrag) = 1
 case('grainsize')
    read(valstring,*,iostat=ierr) grainsizecgs
    grainsize(1) = grainsizecgs/udist
    !--no longer a compulsory parameter
 case('graindens')
    read(valstring,*,iostat=ierr) graindenscgs
    udens = umass/udist**3
    graindens(1) = graindenscgs/udens
    !--no longer a compulsory parameter
 case('K_code')
    read(valstring,*,iostat=ierr) K_code(1)
    igot(iKcode) = 1
 case('icut_backreaction')
    read(valstring,*,iostat=ierr) icut_backreaction
    igot(ibackreact) = 1
 case('ilimitdustflux')
    read(valstring,*,iostat=ierr) ilimitdustflux
    !--no longer a compulsory parameter
 case('irecon')
    read(valstring,*,iostat=ierr) irecon
 case('drag_implicit')
    read(valstring,*,iostat=ierr) drag_implicit
 case default
    imatch = .false.
 end select

 if (name(1:6) == 'K_code') then
    str = trim(name(7:len(name)))
    read(str,*,iostat=ierr) int
    if (ierr /= 0) int = 1
    if (int > 0) read(valstring,*,iostat=ierr) K_code(int)
    igot(iKcode) = 1
    imatch = .true.
 endif

 ineed = 0

 !--Parameters needed by all combinations
 ineed(iidrag)     = 1
 ineed(ibackreact) = 1

 !--Parameters specific to particular setups
 select case(idrag)
 case(0,1)
    ineed(iKcode) = 0
 case(2,3)
    ineed(iKcode) = 0 !1
 case default
    call fatal('read_dust_infile_options','Invalid option',var='idrag',ival=idrag)
 end select

 !--Check that we have just the *necessary* parameters
 if (all(igot >= ineed)) igotall = .true.

end subroutine read_options_dust

real function get_viscmol_nu(spsoundgas,rhogas)
 real,intent(in)  :: spsoundgas,rhogas

 get_viscmol_nu = cste_mu*seff*spsoundgas/rhogas

end function get_viscmol_nu

end module dust
