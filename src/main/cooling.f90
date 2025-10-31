!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling
!
! Gas cooling
!  Current options:
!     0 = none
!     1 = cooling library                 [implicit/exact]
!     2 = cooling library                 [explicit]
!     3 = Gammie cooling                  [explicit]
!     5 = Koyama & Inutuska (2002)        [explicit]
!     6 = Koyama & Inutuska (2002)        [implicit]
!     7 = Gammie cooling power law        [explicit]
!     9 = Young et al. (2024)             [implicit]
!
! :References:
!   Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - C_cool   : *factor controlling cooling timestep*
!   - Tfloor   : *temperature floor (K); on if > 0*
!   - icooling : *cooling function (0=off, 1=library (step), 2=library (force),*
!
! :Dependencies: chem, cooling_gammie, cooling_gammie_PL, cooling_ism,
!   cooling_koyamainutsuka, cooling_molecular, cooling_radapprox,
!   cooling_solver, dim, eos, infile_utils, io, part, physcon, timestep,
!   units, viscosity
!

 use eos,      only:icooling
 use timestep, only:C_cool
 use cooling_solver, only:T0_value,lambda_shock_cgs ! expose to other routines

 implicit none
 character(len=*), parameter :: label = 'cooling'

 public :: init_cooling,energ_cooling
 public :: write_options_cooling,read_options_cooling

 logical, public :: cooling_in_step  = .false.

 !--Minimum temperature (failsafe to prevent u < 0); optional for ALL cooling options
 real,    public :: Tfloor = 0.                     ! [K]; set in .in file.  On if Tfloor > 0.
 real,    public :: ufloor = 0.                     ! [code units]; set in init_cooling
 public :: T0_value,lambda_shock_cgs ! expose to public

 private

contains

!-----------------------------------------------------------------------
!+
!  Initialise cooling
!+
!-----------------------------------------------------------------------
subroutine init_cooling(id,master,iprint,ierr)
 use dim,               only:maxvxyzu
 use units,             only:unit_ergg
 use physcon,           only:mass_proton_cgs,kboltz
 use io,                only:error,fatal,warning
 use eos,               only:gamma,gmw,ieos
 use part,              only:iHI
 use cooling_ism,       only:init_cooling_ism,abund_default
 use cooling_koyamainutsuka, only:init_cooling_KI02
 use cooling_solver,         only:init_cooling_solver
 use cooling_radapprox, only:init_star
 use viscosity,         only:irealvisc

 integer, intent(in)  :: id,master,iprint
 integer, intent(out) :: ierr

 cooling_in_step = .true.
 ierr = 0
 select case(icooling)
 case(4,8)
    if (id==master) write(iprint,*) 'initialising ISM cooling functions...'
    abund_default(iHI) = 1.
    call init_cooling_ism()
    if (icooling==8) cooling_in_step = .false.
 case(9)
    if (ieos /= 24 )  call fatal('cooling','icooling=9 requires ieos=24',&
         var='ieos',ival=ieos)
    if (irealvisc > 0) call warning('cooling',&
         'Using real viscosity will affect optical depth estimate',var='irealvisc',ival=irealvisc)
    call init_star()
 case(6)
    call init_cooling_KI02(ierr)
 case(5)
    call init_cooling_KI02(ierr)
    cooling_in_step = .false.
 case(3)
    ! Gammie
    cooling_in_step = .false.
 case(7)
    ! Gammie PL
    cooling_in_step = .false.
 case(2)
    cooling_in_step = .false.
 case default
    call init_cooling_solver(ierr)
 end select

 !--calculate the energy floor in code units
 if (icooling == 9) then
    ufloor = 0. ! because we calculate & use umin separately
 elseif (Tfloor > 0.) then
    if (gamma > 1.) then
       ufloor = kboltz*Tfloor/((gamma-1.)*gmw*mass_proton_cgs)/unit_ergg
    else
       ufloor = 3.0*kboltz*Tfloor/(2.0*gmw*mass_proton_cgs)/unit_ergg
    endif
    if (maxvxyzu < 4) ierr = 1
 else
    ufloor = 0.
 endif

end subroutine init_cooling

!-----------------------------------------------------------------------
!
!   this routine returns the effective cooling rate du/dt
!
!-----------------------------------------------------------------------

subroutine energ_cooling(xi,yi,zi,ui,rho,dt,divv,dudt,Tdust_in,mu_in,gamma_in,K2_in,kappa_in,abund_in,duhydro,ipart)
 use io,      only:fatal
 use dim,     only:nabundances
 use eos,     only:gmw,gamma,ieos,get_temperature_from_u
 use chem,    only:get_extra_abundances
 use cooling_ism,            only:nabn,energ_cooling_ism,abund_default,abundc,abunde,abundo,abundsi
 use cooling_gammie,         only:cooling_Gammie_explicit
 use cooling_gammie_PL,      only:cooling_Gammie_PL_explicit
 use cooling_solver,         only:energ_cooling_solver
 use cooling_koyamainutsuka, only:cooling_KoyamaInutsuka_explicit,&
                                  cooling_KoyamaInutsuka_implicit
 use cooling_radapprox,      only:radcool_update_du

 real(kind=4), intent(in)   :: divv               ! in code units
 real, intent(in)           :: xi,yi,zi,ui,rho,dt                      ! in code units
 real, intent(in), optional :: Tdust_in,mu_in,gamma_in,K2_in,kappa_in   ! in cgs
 real, intent(in), optional :: abund_in(nabn),duhydro
 integer, intent(in), optional :: ipart
 real, intent(out)          :: dudt                                ! in code units
 real                       :: mui,gammai,Tgas,Tdust,K2,kappa
 real :: abundi(nabn)

 dudt   = 0.
 mui    = gmw
 gammai = gamma
 kappa  = 0.
 K2     = 0.
 if (present(gamma_in)) gammai = gamma_in
 if (present(mu_in))    mui        = mu_in
 if (present(K2_in))    K2        = K2_in
 if (present(kappa_in)) kappa     = kappa_in
 if (gammai < 1.) call fatal('energ_cooling','gamma < 1')
 if (present(abund_in)) then
    abundi = abund_in
 elseif (icooling==4 .or. icooling==8) then
    call get_extra_abundances(abund_default,nabundances,abundi,nabn,mui,&
         abundc,abunde,abundo,abundsi)
 endif

 Tgas  = get_temperature_from_u(ieos,xi,yi,zi,rho,ui,gammai,mui)
 Tdust = Tgas
 if (present(Tdust_in)) Tdust = Tdust_in

 select case (icooling)
 case (6)
    call cooling_KoyamaInutsuka_implicit(ui,rho,dt,dudt)
 case (5)
    call cooling_KoyamaInutsuka_explicit(rho,Tgas,dudt)
 case (4,8)
    call energ_cooling_ism(ui,rho,divv,mui,abundi,dudt)
 case (3)
    call cooling_Gammie_explicit(xi,yi,zi,ui,dudt)
 case (7)
    call cooling_Gammie_PL_explicit(xi,yi,zi,ui,dudt)
 case (9)
    call radcool_update_du(ipart,xi,yi,zi,rho,ui,duhydro,Tfloor)
 case default
    call energ_cooling_solver(ui,dudt,rho,dt,mui,gammai,Tdust,K2,kappa)
 end select

end subroutine energ_cooling

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling(iunit)
 use infile_utils,      only:write_inopt
 use cooling_ism,       only:write_options_cooling_ism
 use cooling_gammie,    only:write_options_cooling_gammie
 use cooling_gammie_PL, only:write_options_cooling_gammie_PL
 use cooling_molecular, only:write_options_molecularcooling
 use cooling_solver,    only:write_options_cooling_solver
 use cooling_radapprox, only:write_options_cooling_radapprox
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling cooling'
 call write_inopt(icooling,'icooling','cooling function (0=off, 1=library (step), 2=library (force),'// &
                     '3=Gammie, 4=ISM, 5,6=KI02, 7=powerlaw, 9=radiative approx)',iunit)
 if (icooling > 0) call write_inopt(C_cool,'C_cool','factor controlling cooling timestep',iunit)
 select case(icooling)
 case(0,5,6)
    ! do nothing
 case(4,8)
    call write_options_cooling_ism(iunit)
 case(3)
    call write_options_cooling_gammie(iunit)
 case(7)
    call write_options_cooling_gammie_PL(iunit)
 case(9)
    call write_options_cooling_radapprox(iunit)
 case default
    call write_options_cooling_solver(iunit)
 end select
 if (icooling > 0) call write_inopt(Tfloor,'Tfloor','temperature floor (K); on if > 0',iunit)

end subroutine write_options_cooling

!-----------------------------------------------------------------------
!+
!  reads options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling(db,nerr)
 use io,                only:fatal
 use eos,               only:ipdv_heating,ishock_heating,eos_allows_shock_and_work,ieos
 use cooling_gammie,    only:read_options_cooling_gammie
 use cooling_gammie_PL, only:read_options_cooling_gammie_PL
 use cooling_ism,       only:read_options_cooling_ism
 use cooling_molecular, only:read_options_molecular_cooling
 use cooling_solver,    only:read_options_cooling_solver
 use cooling_radapprox, only:read_options_cooling_radapprox
 use infile_utils,      only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 character(len=*), parameter  :: label = 'read_infile'

 call read_inopt(icooling,'icooling',db,errcount=nerr,min=0,max=9,default=icooling)
 if (icooling > 0 .and. .not. eos_allows_shock_and_work(ieos)) &
    call fatal(label,'cooling requires adiabatic eos (e.g. ieos=2)')
 if (icooling > 0 .and. (ipdv_heating <= 0 .or. ishock_heating <= 0)) &
    call fatal(label,'cooling requires shock and work contributions')
 call read_inopt(C_cool,'C_cool',db,errcount=nerr,min=0.,default=C_cool)
 call read_inopt(Tfloor,'Tfloor',db,errcount=nerr,min=0.,default=Tfloor)

 select case(icooling)
 case(0,5,6)
    ! do nothing
 case(4,8)
    call read_options_cooling_ism(db,nerr)
 case(3)
    call read_options_cooling_gammie(db,nerr)
 case(7)
    call read_options_cooling_gammie_PL(db,nerr)
 case(9)
    call read_options_cooling_radapprox(db,nerr)
 case default
    call read_options_cooling_solver(db,nerr)
 end select

end subroutine read_options_cooling

end module cooling
