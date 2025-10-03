!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_gasradrec
!
! EoS from HORMONE that includes internal energy from ideal gas,
!  radiation, and recombination (H2, H, He) (Hirai et al., 2020)
!
! :References: Appendix C, https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.1154H/abstract
!
! :Owner: Mike Lau
!
! :Runtime parameters:
!   - irecomb : *recombination energy to include. 0=H2+H+He, 1=H+He, 2=He, 3=none*
!
! :Dependencies: infile_utils, io, ionization_mod, physcon
!
 implicit none
 integer, public :: irecomb = 0 ! types of recombination energy to include for ieos=20
 public :: equationofstate_gasradrec,calc_uT_from_rhoP_gasradrec,calc_uP_from_rhoT_gasradrec,&
           read_options_eos_gasradrec,write_options_eos_gasradrec,eos_info_gasradrec,init_eos_gasradrec
 private
 real, parameter :: eoserr=1.e-15,W4err=1.e-2

contains
!-----------------------------------------------------------------------
!+
!  EoS from HORMONE (Hirai et al., 2020).
!  Note: eint is internal energy per unit volume, xi must be included
!  to include Trad term in sound speed
!+
!-----------------------------------------------------------------------
subroutine equationofstate_gasradrec(d,eint,T,imu,X,Y,p,cf,gamma_eff,cveff_out,do_radiation,xi)
 use ionization_mod, only:get_erec_cveff,get_imurec
 use physcon,        only:radconst,Rg
 use io,             only:fatal
 real, intent(in)    :: d,eint,X,Y
 real, intent(in), optional :: xi
 logical, intent(in), optional :: do_radiation
 real, intent(inout) :: T,imu ! imu is 1/mu, an output
 real, intent(out)   :: p,cf,gamma_eff
 real, intent(out), optional :: cveff_out
 real                :: corr,erec,derecdT,Tdot,logd,dt,Tguess,cveff,dcveffdlnT,cs2
 logical             :: do_radiation_local
 integer, parameter  :: nmax = 500
 integer n

 if (.not. present(do_radiation)) then
    do_radiation_local = .false.
 else
    do_radiation_local = do_radiation
 endif
 corr=huge(0.); Tdot=0.; logd=log10(d); dt=0.9; Tguess=T

 do n = 1,nmax
    call get_erec_cveff(logd,T,X,Y,erec,cveff,derecdT,dcveffdlnT)
    if (d*erec>=eint) then ! avoid negative thermal energy
       T = 0.9*T; Tdot=0.;cycle
    endif
    if (do_radiation_local) then
       corr = (eint-d*(Rg*cveff*T+erec)) &
              / ( -d*(Rg*(cveff+dcveffdlnT)+derecdT) )
    else
       corr = (eint-(radconst*T**3+Rg*d*cveff)*T-d*erec) &
              / ( -4.*radconst*T**3-d*(Rg*(cveff+dcveffdlnT)+derecdT) )
    endif
    if (-corr > 10.*T) then  ! do not let temperature guess increase by more than a factor of 1000
       T = 10.*T
       Tdot = 0.
    elseif (abs(corr) > W4err*T) then
       T = T + Tdot*dt
       Tdot = (1.-2.*dt)*Tdot - dt*corr
    else
       T = T-corr*dt
       Tdot = 0.
       dt = 1.
    endif
    if (abs(corr)<max(eoserr*T,eoserr)) exit
    if (n>50) dt=0.5
    if (n>100) dt=0.25
 enddo
 call get_imurec(logd,T,X,Y,imu)
 if (n > nmax) then
    print*,'d=',d,'eint=',eint/d,'Tguess=',Tguess,'mu=',1./imu,'T=',T,'erec=',erec
    print*,'n = ',n,' nmax = ',n,' correction is ',abs(corr),' needs to be < ',eoserr*T
    call fatal('eos_gasradrec','Failed to converge on temperature in equationofstate_gasradrec')
 endif
 if (do_radiation_local) then
    p = Rg*imu*d*T
 else
    p = ( Rg*imu*d + radconst*T**3/3. )*T
 endif
 if (present(xi)) then
    cs2 = get_cs2(d,T,X,Y,do_radiation_local,xi)  ! not used at the moment
 else
    cs2 = get_cs2(d,T,X,Y,do_radiation_local)
 endif
 gamma_eff=cs2*d/p
 cf = sqrt(cs2)
 if (present(cveff_out)) cveff_out = cveff

end subroutine equationofstate_gasradrec

!-----------------------------------------------------------------------
!+
!  To compute sound speed squared from d and T
!+
!-----------------------------------------------------------------------
function get_cs2(d,T,X,Y,do_radiation,xi) result(cs2)
 use ionization_mod, only:get_erec_cveff,get_imurec
 use physcon,        only:radconst,Rg
 real, intent(in)              :: d,T,X,Y
 real, intent(in), optional    :: xi
 logical, intent(in), optional :: do_radiation
 logical :: do_radiation_local
 real    :: cs2
 real    :: erec,cveff,derecdT,dcveffdlnT,imurec,dimurecdlnT,dimurecdlnd
 real    :: logd,deraddT

 if (.not. present(do_radiation)) then
    do_radiation_local = .false.
 else
    do_radiation_local = do_radiation
 endif

 logd=log10(d)
 call get_erec_cveff(logd,T,X,Y,erec,cveff,derecdT,dcveffdlnT)
 call get_imurec(logd,T,X,Y,imurec,dimurecdlnT,dimurecdlnd)
 if (do_radiation_local) then
    cs2 = Rg*(imurec+dimurecdlnd)*T &
          + ( Rg*(imurec+dimurecdlnT))**2*T &
          / (Rg*(cveff+dcveffdlnT)+derecdT)
    if (present(xi)) cs2 = cs2 + 4.*xi/9.  ! need xi for Tgas\=Trad
 else
    deraddT = 4.*radconst*T**3/d
    cs2 = Rg*(imurec+dimurecdlnd)*T &
          + ( Rg*(imurec+dimurecdlnT)+deraddT/3.)**2*T &
          / (Rg*(cveff+dcveffdlnT)+deraddT+derecdT)
 endif

end function get_cs2

!-----------------------------------------------------------------------
!+
!  Compute pressure, internal energy, and mean molecular weight from
!  density and temperature
!+
!-----------------------------------------------------------------------
subroutine calc_uP_from_rhoT_gasradrec(rho,T,X,Y,eint,p,imu,do_radiation)
 use ionization_mod, only:get_erec_cveff,get_imurec
 use physcon,        only:Rg,radconst
 real, intent(in)  :: rho,T,X,Y
 logical, intent(in), optional :: do_radiation
 real, intent(out) :: eint,p,imu
 real              :: logrho,erec,cveff
 logical           :: do_radiation_local

 if (.not. present(do_radiation)) then
    do_radiation_local = .false.
 else
    do_radiation_local = do_radiation
 endif

 logrho = log10(rho)
 call get_erec_cveff(logrho,T,X,Y,erec,cveff)
 call get_imurec(logrho,T,X,Y,imu)

 if (do_radiation_local) then
    p = Rg*imu*rho*T
    eint = Rg*cveff*T + erec
 else
    p = ( Rg*imu*rho + radconst*T**3/3. )*T
    eint = (Rg*cveff + radconst*T**3/rho )*T + erec
 endif

end subroutine calc_uP_from_rhoT_gasradrec

!-----------------------------------------------------------------------
!+
!  Solve for u and T (and mu) from rho and P (ideal gas + radiation + recombination)
!  Inputs and outputs in cgs units
!+
!-----------------------------------------------------------------------
subroutine calc_uT_from_rhoP_gasradrec(rhoi,presi,X,Y,T,u,mui,ierr,do_radiation)
 use ionization_mod, only:get_imurec
 use physcon,        only: radconst,Rg
 real, intent(in)              :: rhoi,presi,X,Y
 logical, intent(in), optional :: do_radiation
 real, intent(inout)           :: T
 real, intent(out)             :: u
 real, optional, intent(out)   :: mui
 integer, intent(out)          :: ierr
 integer                       :: n
 real                          :: logrhoi,imu,dimurecdlnT,dT,Tdot,corr,p1
 logical                       :: do_radiation_local
 integer, parameter            :: nmax = 500

 if (.not. present(do_radiation)) then
    do_radiation_local = .false.
 else
    do_radiation_local = do_radiation
 endif

 if (T <= 0.) T = min((3.*presi/radconst)**0.25, presi/(rhoi*Rg))  ! initial guess for temperature
 ierr = 0
 corr = huge(0.)
 Tdot = 0.
 dT = 0.9
 logrhoi = log10(rhoi)
 do n = 1,nmax
    call get_imurec(logrhoi,T,X,Y,imu,dimurecdlnT)
    if (do_radiation_local) then
       corr = ( presi - rhoi*Rg*imu*T ) / &
              ( -rhoi*Rg * ( imu + dimurecdlnT ) )
    else
       corr = ( presi - (rhoi*Rg*imu + radconst*T**3/3.)*T ) / &
              ( -rhoi*Rg * ( imu + dimurecdlnT ) - 4.*radconst*T**3/3. )
    endif
    if (abs(corr)>W4err*T) then
       T = T + Tdot*dT
       Tdot = (1.-2.*dT)*Tdot - dT*corr
    else
       T = T - corr
       Tdot = 0.
    endif
    if (abs(corr) < eoserr*T) exit
    if (n > 50) dT = 0.5
 enddo
 if (n > nmax) then
    print*,'Error in calc_uT_from_rhoP_hormone'
    print*,'rho=',rhoi,'P=',presi,'mu=',1./imu
    ierr = 1
    return
 endif
 call calc_uP_from_rhoT_gasradrec(rhoi,T,X,Y,u,p1,imu,do_radiation_local)
 mui = 1./imu

end subroutine calc_uT_from_rhoP_gasradrec

!-----------------------------------------------------------------------
!+
!  Initialise eos by setting ionisation energy array according to
!  irecomb option
!+
!-----------------------------------------------------------------------
subroutine init_eos_gasradrec(ierr)
 use ionization_mod, only:ionization_setup,eion
 integer, intent(out) :: ierr

 ierr = 0
 call ionization_setup
 if (irecomb == 1) then
    eion(1) = 0.  ! H and He recombination only (no recombination to H2)
 elseif (irecomb == 2) then
    eion(1:2) = 0.  ! He recombination only
 elseif (irecomb == 3) then
    eion(1:4) = 0.  ! No recombination energy
 elseif (irecomb < 0 .or. irecomb > 3) then
    ierr = 1
 endif
 write(*,'(1x,a,i1,a)') 'Initialising gas+rad+rec EoS (irecomb = ',irecomb,')'

end subroutine init_eos_gasradrec

!----------------------------------------------------------------
!+
!  print eos information
!+
!----------------------------------------------------------------
subroutine eos_info_gasradrec(iprint)
 integer, intent(in) :: iprint

 write(iprint,"(/,a,i1)") ' Gas+rad+rec EoS: irecomb = ',irecomb

end subroutine eos_info_gasradrec

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_eos_gasradrec(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(irecomb,'irecomb',db,errcount=nerr,min=0,max=3)

end subroutine read_options_eos_gasradrec

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_gasradrec(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(irecomb,'irecomb','recombination energy to include. 0=H2+H+He, 1=H+He, 2=He, 3=none',iunit)

end subroutine write_options_eos_gasradrec

end module eos_gasradrec
