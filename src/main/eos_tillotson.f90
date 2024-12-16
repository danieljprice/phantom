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
! :References:
!   Benz, Slattery & Cameron (1986), Icarus 66, 515-535
!   Benz & Asphaug (1999), I
!   Kegerreis et al. (2019), MNRAS 487, 5029-5040
!   https://apps.dtic.mil/sti/pdfs/AD0486711.pdf
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
 real :: rho_0 = 2.7 ! g/cm^3 zero-pressure density (Basalt) from Benz & Asphaug 1999
 ! aparam, bparam, A, B, u_0 material-dependent Tillotson parameters
 real :: aparam = 0.5 , bparam = 1.5 , alpha = 5. , beta = 5.
 real :: A = 2.67e11 , B = 2.67e11 ! erg/cm^3
 real :: u_0 = 4.87e12 ! erg/g
 real :: u_iv = 4.72e10 ! erg/g
 real :: u_cv = 1.82e11 ! erg/g

 public :: init_eos_tillotson, equationofstate_tillotson
 public :: eos_info_tillotson, read_options_eos_tillotson, write_options_eos_tillotson

 private

contains
!-----------------------------------------------------------------------
!+
!  initialise the equation of state
!+
!-----------------------------------------------------------------------
subroutine init_eos_tillotson(ierr)
 integer, intent(out) :: ierr

 ierr = 0

end subroutine init_eos_tillotson

!-----------------------------------------------------------------------
!+
!  EoS from Tillotson et al. (1962) ; Implementation from
!  Benz et al. (1986) and Kegerreis et al. (2019)
!  notes:
!   - u_iv (incipient vaporisation) is called u_s  in Benz et al. 1986
!   - u_cv (complete vaporisation) is called u_s' in Benz et al. 1986
!+
!-----------------------------------------------------------------------
subroutine equationofstate_tillotson(rho,u,pressure,spsound,gamma)
 real, intent(inout) :: rho,u
 real, intent(out)   :: pressure, spsound, gamma
 real :: eta,mu_t,nu_t,omega,spsoundmin,pc,cc2,pe,ce2,denom

 eta = rho/rho_0
 mu_t = eta - 1.
 nu_t = (1./eta) - 1.  ! Kegerreis et al. 2020
 omega = (u / (u_0*eta**2)) + 1.
 spsoundmin = sqrt(A / rho_0) ! wave speed estimate Benz and Asphaug

 if (rho >= rho_0 .or. u < u_iv) then
    ! compressed or cold state
    call get_compressed_state(pc,cc2,rho,u,mu_t,omega,eta,aparam,bparam,A,B)
    pressure = pc
    spsound  = sqrt(cc2)
 elseif (rho >= rho_0 .and. u < u_cv) then
    ! hybrid state
    call get_compressed_state(pc,cc2,rho,u,mu_t,omega,eta,aparam,bparam,A,B)
    call get_expanded_state(pe,ce2,rho,u,mu_t,omega,eta,aparam,bparam,A,nu_t,alpha,beta)
    denom = 1./(u_cv - u_iv)
    pressure = ((u - u_iv)*pe + (u_cv - u)*pc)*denom
    spsound  = sqrt((u-u_iv)*ce2 + (u_cv - u)*cc2)*denom
 else
    ! expanded state
    call get_expanded_state(pe,ce2,rho,u,mu_t,omega,eta,aparam,bparam,A,nu_t,alpha,beta)
    pressure = pe
    spsound  = sqrt(ce2)
 endif

end subroutine equationofstate_tillotson

!----------------------------------------------------------------
!+
!  functional form in compressed or cold state
!+
!----------------------------------------------------------------
pure subroutine get_compressed_state(pr,cs2,rho,u,mu_t,omega,eta,aa,bb,A,B)
 real, intent(out) :: pr,cs2
 real, intent(in)  :: rho,u,mu_t,omega,eta,aa,bb,A,B

 pr  = (aa + bb / omega)*rho*u + A*mu_t + B*mu_t**2
 cs2 = pr/rho*(1. + aa + bb/omega) &
     + bb*(omega-1.)/omega**2*(2.*u - pr/rho) &
     + 1./rho*(A + B*(eta**2 - 1))

end subroutine get_compressed_state

!----------------------------------------------------------------
!+
!  functional form in expanded state
!+
!----------------------------------------------------------------
pure subroutine get_expanded_state(pr,cs2,rho,u,mu_t,omega,eta,aa,bb,A,nu_t,alpha,beta)
 real, intent(out) :: pr,cs2
 real, intent(in)  :: rho,u,mu_t,omega,eta,aa,bb,A,nu_t,alpha,beta
 real :: exp_alpha,exp_beta

 exp_alpha = exp(-alpha*nu_t**2)
 exp_beta  = exp(-beta*nu_t)
 pr  = aa*rho*u + (bb*rho*u/omega + A*mu_t*exp_beta)*exp_alpha
 cs2 = pr/rho*(1. + aa + bb/omega*exp_alpha) &
     + (bb*rho*u/(omega**2*eta**2)*(1./(u_0*rho)*(2.*u - pr/rho) + 2*alpha*nu_t*omega/rho_0) &
     +  A/rho_0*(1. + mu_t/eta**2*(beta + 2.*alpha*nu_t - eta))*exp_beta) * exp_alpha

end subroutine get_expanded_state

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
 case('rho_0')
    read(valstring,*,iostat=ierr) rho_0
    if ((rho_0 < 0.)) call fatal(label,'rho_0 < 0')
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
 case('u_0')
    read(valstring,*,iostat=ierr) u_0
    if ((u_0 < 0.)) call fatal(label,'u_0 < 0')
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

 call write_inopt(rho_0,'rho_0','reference density g/cm^3',iunit)
 call write_inopt(aparam,'aparam','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(bparam,'bparam','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(A,'A','material-dependent Tillotson parameter, erg/cm^3',iunit)
 call write_inopt(B,'B','material-dependent Tillotson parameter, erg/cm^3',iunit)
 call write_inopt(alpha,'alpha_t','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(beta,'beta_t','material-dependent Tillotson parameter, unitless',iunit)
 call write_inopt(u_0,'u_0','material-dependent Tillotson parameter, erg/g',iunit)

end subroutine write_options_eos_tillotson

end module eos_tillotson