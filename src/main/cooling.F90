!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!
! :References:
!   Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - C_cool   : *factor controlling cooling timestep*
!   - Tfloor   : *temperature floor (K); on if > 0*
!   - icooling : *cooling function (0=off, 1=cooling library (step), 2=cooling library (force),*
!
! :Dependencies: chem, cooling_gammie, cooling_ism, cooling_koyamainutsuka,
!   cooling_molecular, cooling_solver, dim, eos, infile_utils, io, options,
!   part, physcon, timestep, units
!

 use options,  only:icooling
 use timestep, only:C_cool
 use cooling_solver, only:T0_value ! expose to other routines

 implicit none
 character(len=*), parameter :: label = 'cooling'

 public :: init_cooling,energ_cooling
 public :: write_options_cooling,read_options_cooling

 logical, public :: cooling_in_step  = .false.

 !--Minimum temperature (failsafe to prevent u < 0); optional for ALL cooling options
 real,    public :: Tfloor = 0.                     ! [K]; set in .in file.  On if Tfloor > 0.
 real,    public :: ufloor = 0.                     ! [code units]; set in init_cooling
 public :: T0_value ! expose to public

 private

contains

!-----------------------------------------------------------------------
!+
!  Initialise cooling
!+
!-----------------------------------------------------------------------
subroutine init_cooling(id,master,iprint,ierr)
 use dim,               only:maxvxyzu,h2chemistry
 use units,             only:unit_ergg
 use physcon,           only:mass_proton_cgs,kboltz
 use io,                only:error
 use eos,               only:gamma,gmw
 use cooling_ism,       only:init_cooling_ism
 use chem,              only:init_chem
 use cooling_molecular,      only:init_cooling_molec
 use cooling_koyamainutsuka, only:init_cooling_KI02
 use cooling_solver,         only:init_cooling_solver

 integer, intent(in)  :: id,master,iprint
 integer, intent(out) :: ierr

 cooling_in_step = .true.
 ierr = 0
 if (h2chemistry) then
    if (id==master) write(iprint,*) 'initialising cooling function...'
    call init_chem()
    call init_cooling_ism()
 else
    select case(icooling)
    case(6)
       call init_cooling_KI02(ierr)
    case(5)
       call init_cooling_KI02(ierr)
       cooling_in_step = .false.
    case(4)
       ! Initialise molecular cooling
       call init_cooling_molec
    case(3)
       ! Gammie
       cooling_in_step = .false.
    case default
       call init_cooling_solver(ierr)
    end select
 endif

 !--calculate the energy floor in code units
 if (Tfloor > 0.) then
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
subroutine energ_cooling(xi,yi,zi,ui,dudt,rho,dt,Tdust_in,mu_in,gamma_in,K2_in,kappa_in)
 use io,      only:fatal
 use eos,     only:gmw,gamma
 use physcon, only:Rg
 use units,   only:unit_ergg
 use cooling_gammie,         only:cooling_Gammie_explicit
 use cooling_solver,         only:energ_cooling_solver
 use cooling_koyamainutsuka, only:cooling_KoyamaInutsuka_explicit,&
                                  cooling_KoyamaInutsuka_implicit

 real, intent(in)           :: xi,yi,zi,ui,rho,dt                  ! in code units
 real, intent(in), optional :: Tdust_in,mu_in,gamma_in,K2_in,kappa_in   ! in cgs
 real, intent(out)          :: dudt                                ! in code units
 real                       :: mu,polyIndex,T_on_u,Tgas,Tdust,K2,kappa

 dudt       = 0.
 mu         = gmw
 polyIndex  = gamma
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 Tgas   = T_on_u*ui
 Tdust  = Tgas
 kappa  = 0.
 K2     = 0.
 if (present(gamma_in)) polyIndex = gamma_in
 if (present(mu_in))    mu        = mu_in
 if (present(Tdust_in)) Tdust     = Tdust_in
 if (present(K2_in))    K2        = K2_in
 if (present(kappa_in)) kappa     = kappa_in

 select case (icooling)
 case (6)
    call cooling_KoyamaInutsuka_implicit(ui,rho,dt,dudt)
 case (5)
    call cooling_KoyamaInutsuka_explicit(rho,Tgas,dudt)
 case (4)
    !call cooling_molecular
 case (3)
    call cooling_Gammie_explicit(xi,yi,zi,ui,dudt)
 case default
    call energ_cooling_solver(ui,dudt,rho,dt,mu,polyIndex,Tdust,K2,kappa)
 end select

end subroutine energ_cooling

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling(iunit)
 use infile_utils,      only:write_inopt
 use part,              only:h2chemistry
 use cooling_ism,       only:write_options_cooling_ism
 use cooling_gammie,    only:write_options_cooling_gammie
 use cooling_molecular, only:write_options_molecularcooling
 use cooling_solver,    only:write_options_cooling_solver
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling cooling'
 call write_inopt(C_cool,'C_cool','factor controlling cooling timestep',iunit)
 if (h2chemistry) then
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=on)',iunit)
    if (icooling > 0) call write_options_cooling_ism(iunit)
 else
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=cooling library (step), 2=cooling library (force),'// &
                     '3=Gammie, 5,6=KI02)',iunit)
    select case(icooling)
    case(0,4,5,6)
       ! do nothing
    case(3)
       call write_options_cooling_gammie(iunit)
    case default
       call write_options_cooling_solver(iunit)
    end select
 endif
 if (icooling > 0) call write_inopt(Tfloor,'Tfloor','temperature floor (K); on if > 0',iunit)

end subroutine write_options_cooling

!-----------------------------------------------------------------------
!+
!  reads options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling(name,valstring,imatch,igotall,ierr)
 use part,              only:h2chemistry
 use io,                only:fatal
 use cooling_gammie,    only:read_options_cooling_gammie
 use cooling_ism,       only:read_options_cooling_ism
 use cooling_molecular, only:read_options_molecular_cooling
 use cooling_solver,    only:read_options_cooling_solver
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 logical :: igotallism,igotallmol,igotallgammie,igotallfunc

 imatch        = .true.
 igotall       = .false.  ! cooling options are compulsory
 igotallism    = .true.
 igotallmol    = .true.
 igotallgammie = .true.
 igotallfunc   = .true.

 select case(trim(name))
 case('icooling')
    read(valstring,*,iostat=ierr) icooling
    ngot = ngot + 1
 case('C_cool')
    read(valstring,*,iostat=ierr) C_cool
    ngot = ngot + 1
 case('Tfloor')
    ! not compulsory to read in
    read(valstring,*,iostat=ierr) Tfloor
 case default
    imatch = .false.
    if (h2chemistry) then
       call read_options_cooling_ism(name,valstring,imatch,igotallism,ierr)
    else
       select case(icooling)
       case(0,4,5,6)
          ! do nothing
       case(3)
          call read_options_cooling_gammie(name,valstring,imatch,igotallgammie,ierr)
       case default
          call read_options_cooling_solver(name,valstring,imatch,igotallfunc,ierr)
       end select
    endif
 end select
 ierr = 0
 if (h2chemistry .and. igotallism .and. ngot >= 2) then
    igotall = .true.
 elseif (icooling >= 0 .and. ngot >= 2 .and. igotallgammie .and. igotallfunc) then
    igotall = .true.
 else
    igotall = .false.
 endif

end subroutine read_options_cooling

end module cooling
