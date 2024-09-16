!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_tillotson
!
! EoS from Tillotson et al. (1962)
!
! :References: https://apps.dtic.mil/sti/pdfs/AD0486711.pdf
!
! Implementation from Benz and Asphaug (1999)
!
! :Owner:
!
! :Runtime parameters: None
!
! :Dependencies: 
!
 implicit none
 real :: rho0 = 2.7 ! g/cm^3 zero-pressure density (Basalt) from Benz & Asphaug 1999
! aparam, bparam, A, B, energy0 material-dependent Tillotson parameters (Basalt)
 real :: aparam = 0.5 , bparam = 1.5 , alpha = 5. , beta = 5. 
 real :: A = 2.67e11 , B = 2.67e11 ! erg/cm^3
 real :: u0 = 4.87e12 ! erg/g
 real :: u_iv = 4.72e10 ! erg/g
 real :: u_cv = 1.82e11 ! erg/g
 real :: rho_iv = 2.57 ! g/cm^3 incipient vaporisation density of gabbroic anorthosite lpp from Ahrens and Okeefe (1977) table 1
!  real :: rho_cv = 8.47 ! g/cm^3 complete vaporisation density of gabbroic anorthosite lpp from Ahrens and Okeefe (1977) table 2
 private :: rho0, aparam, bparam, A, B, alpha, beta, u0, rho_iv, u_iv, u_cv
!  private :: spsoundmin, spsound2, spsound3, pressure2, pressure3
!  private :: expterm, expterm2, sqrtterm, sqrtterm1, sqrtterm2

contains
!-----------------------------------------------------------------------
!+
!  EoS from Tillotson et al. (1962) ; Implementation from Benz and Asphaug (1999) and Brundage (2013)
!+
!-----------------------------------------------------------------------
subroutine init_eos_tillotson(ierr)
 integer, intent(out) :: ierr

 ierr = 0

end subroutine init_eos_tillotson

subroutine equationofstate_tillotson(rho,u,pressure,spsound,gamma)
 real, intent(in) :: rho,u
 real, intent(out) :: pressure, spsound, gamma
 real :: eta, mu, neu, omega, spsoundmin, spsound2, spsound2h, spsound3, pressure2, pressure3

 eta = rho / rho0
 mu = eta - 1.
 neu = (1. / eta) - 1.  ! Kegerreis et al. 2020
 omega = (u / (u0*eta**2) + 1.)
 spsoundmin = sqrt(A / rho) ! wave speed estimate
!  gamma = 2./3.

! Brundage 2013 eq (2) cold expanded and compressed states
 pressure2 = (( aparam + ( bparam / omega ) ) * rho*u) + A*mu + B*mu**2
 if (pressure2 <= 0.) then
  pressure2 = 0.
 endif

! Brundage 2013 eq (3) completely vaporised
 pressure3 = (aparam*rho*u) + ( ((bparam*rho*u)/omega) + &
              A*mu*exp(-beta*neu)) * exp(-alpha*(neu**2))
 if (pressure3 <= 0.) then
  pressure3 = 0.
 endif

! sound speed from Kegereis et al. 2020
 spsound2 = sqrt( ((pressure2 / rho) * (1. + aparam + (bparam / omega ))) + &
            (((bparam*(omega - 1.)) / (omega**2)) * (2.*u - pressure2/rho)) + &
            ((1./rho)*(A + B*((eta**2) - 1.))))
 if (spsound2 < spsoundmin) then
  spsound2 = spsoundmin
 endif

 spsound3 = sqrt( (pressure3/rho)*( 1. + aparam + ((bparam / omega )*exp(-alpha*(neu**2)))) + &
            ( ((bparam*rho*u)/((omega**2)*(eta**2))) * ((1. /(u0*rho))*((2.*u) - (pressure3/rho)) &
            + ((2.*alpha*neu*omega)/rho0)) + ((A/rho0)*(1. + (mu/eta**2)*((beta*2.*alpha*neu) - eta))) &
            * exp(-beta*neu) )*exp(-alpha*(neu**2)) )
 if (spsound3 < spsoundmin) then
  spsound3 = spsoundmin
 endif

 spsound2h = sqrt( ((pressure2 / rho) * (1. + aparam + (bparam / omega ))) + &
            (((bparam*(omega - 1.)) / (omega**2)) * (2.*u - pressure2/rho)) + &
            ((1./rho)*(A + B*((eta**2) - 1.))) - ((((2.*B)/rho0**2)*rho) - (2.*B/rho0)))
 if (spsound2h < spsoundmin) then
  spsound2h = spsoundmin
 endif
 
 if (rho >= rho0) then !           compressed state
    pressure = pressure2
    spsound = spsound2
 else ! rho < rho0  
    if (u >= u_cv) then !          vaporised expanded state
        pressure = pressure3
        spsound = spsound3
    ! need to check rho ~0 and u >> u_cv (zero density becomes pressure = rho*gamma*u fermi, gamma = 2/3)
    else ! u < u_cv
        if (rho >= rho_iv) then 
            if (u > u_iv) then !   Hybrid state
                pressure = ( (u - u_iv)*pressure3 + (u_cv - u)*pressure2 ) / (u_cv - u_iv)
                spsound = sqrt(( ((u-u_iv)*(spsound3)**2) + ((u_cv-u)*(spsound2)**2) ) / (u_cv - u_iv))
            else ! u <= u_iv       cold expanded 
            pressure = pressure2
            spsound = spsound2
            endif
        else ! rho < rho_iv        Low energy state
            pressure = pressure2 - (B*mu**2)
            spsound = spsound2h
        endif
    endif
 endif

end subroutine equationofstate_tillotson

!----------------------------------------------------------------
!+
!  print eos information
!+
!----------------------------------------------------------------
subroutine eos_info_tillotson(iprint)
 integer, intent(in) :: iprint

 write(iprint,"(/,a)") ' Tillotson EoS'

end subroutine eos_info_tillotson

!-----------------------------------------------------------------------
!+
!  calc internal energy from pressure and density
!+
!-----------------------------------------------------------------------
subroutine calc_uT_from_rhoP_tillotson(rho,pres,temp,u,ierr)
 real, intent(in) :: rho,pres
 real, intent(out) :: temp,u
 integer, intent(out) :: ierr
 real :: eta, mu, neu, omega, expterm, expterm2, sqrtterm1, sqrtterm2, sqrtterm
 
 eta = rho / rho0
 mu = eta - 1.
 neu = (1. / eta) - 1.  ! Kegerreis et al. 2020
 omega = (u / (u0*eta**2) + 1.)
 expterm2 = exp(-alpha*(neu**2))
 expterm = exp(-beta*neu)

 if (rho >= rho0) then !           compressed state
    u = (mu*(eta**2)*u0*(aparam-bparam)*(A+(B*mu))) / &
        (aparam*(A*mu+bparam*(eta**2)*u0*rho+B*(mu**2)-pres))
 elseif (rho >= rho_iv) then !      cold expanded 
     u = (mu*(eta**2)*u0*(aparam-bparam)*(A+(B*mu))) / &
        (aparam*(A*mu+bparam*(eta**2)*u0*rho+B*(mu**2)-pres))
 elseif (rho < rho_iv) then ! rho < rho_iv        Low energy state
     u = (A*mu*(eta**2)*u0*(aparam-bparam)) / &
         (aparam*(A*mu+bparam*(eta**2)*u0*rho-pres))
 else ! rho < rho0  vaporised expanded state
     u = ((eta**2)*u0*(bparam-aparam)*(pres-A*mu*expterm*expterm2) ) / &
         (aparam*rho*(A*mu*expterm + bparam*(eta**2)*u0 - pres))
        ! need to check rho ~0 and u >> u_cv (zero density becomes u = pres / (rho*gamma) fermi, gamma = 2/3)
 endif

 temp = 300.
 ierr = 0
 if (u < 0.) then
  u = 0.
  ierr = 1
 endif

end subroutine calc_uT_from_rhoP_tillotson  
!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
 subroutine read_options_eos_tillotson(name,valstring,imatch,igotall,ierr)
  use io, only:fatal
  character(len=*),  intent(in)  :: name,valstring
  logical,           intent(out) :: imatch,igotall
  integer,           intent(out) :: ierr
  integer,           save        :: ngot  = 0
  character(len=30), parameter   :: label = 'eos_tillotson'
 
  imatch  = .true.
  select case(trim(name))
  case('rho0')
      read(valstring,*,iostat=ierr) rho0
      if ((rho0 < 0.)) call fatal(label,'rho0 < 0')
      ngot = ngot + 1
  case('aparam')
      read(valstring,*,iostat=ierr) aparam
      if ((aparam < 0.)) call fatal(label,'aparam < 0')
      ngot = ngot + 1
  case('bparam')
      read(valstring,*,iostat=ierr) bparam
      if ((bparam < 0.)) call fatal(label,'bparam < 0')
      ngot = ngot + 1
  case('A')
      read(valstring,*,iostat=ierr) A
      if ((A < 0.)) call fatal(label,'A < 0')
      ngot = ngot + 1
  case('B')
      read(valstring,*,iostat=ierr) B
      if ((B < 0.)) call fatal(label,'B < 0')
      ngot = ngot + 1
  case('alpha_t')
      read(valstring,*,iostat=ierr) alpha
      if ((alpha < 0.)) call fatal(label,'alpha < 0')
      ngot = ngot + 1
  case('beta_t')
      read(valstring,*,iostat=ierr) beta
      if ((beta < 0.)) call fatal(label,'beta < 0')
      ngot = ngot + 1
  case('u0')
      read(valstring,*,iostat=ierr) u0
      if ((u0 < 0.)) call fatal(label,'u0 < 0')
      ngot = ngot + 1
  case default
     imatch = .false.
  end select
 
  igotall = (ngot >= 8)
 
 end subroutine read_options_eos_tillotson

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_tillotson(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit
 
 call write_inopt(rho0,'rho0','reference density g/cm^3',iunit)
 call write_inopt(aparam,'aparam','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(bparam,'bparam','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(A,'A','material-dependent Tillotson parameter, erg/cm^3',iunit)
 call write_inopt(B,'B','material-dependent Tillotson parameter, erg/cm^3',iunit)
 call write_inopt(alpha,'alpha_t','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(beta,'beta_t','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(u0,'u0','material-dependent Tillotson parameter, erg/g',iunit)

end subroutine write_options_eos_tillotson

end module eos_tillotson