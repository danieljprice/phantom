!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
 public :: equationofstate_gasradrec,calc_uT_from_rhoP_gasradrec,read_options_eos_gasradrec,&
           write_options_eos_gasradrec,eos_info_gasradrec,init_eos_gasradrec
 private
 real, parameter :: eoserr=1.e-15,W4err=1.e-2

contains
!-----------------------------------------------------------------------
!+
!  EoS from HORMONE (Hirai et al., 2020). Note eint is internal energy per unit volume
!+
!-----------------------------------------------------------------------
subroutine equationofstate_gasradrec(d,eint,T,imu,X,Y,p,cf,gamma_eff)
 use ionization_mod, only:get_erec_imurec
 use physcon,        only:radconst,Rg
 use io,             only:fatal
 real, intent(in)    :: d,eint
 real, intent(inout) :: T,imu ! imu is 1/mu, an output
 real, intent(in)    :: X,Y
 real, intent(out)   :: p,cf,gamma_eff
 real                :: corr,erec,derecdT,dimurecdT,Tdot,logd,dt,Tguess
 integer, parameter  :: nmax = 500
 integer n

 corr=huge(0.); Tdot=0.; logd=log10(d); dt=0.9; Tguess=T

 do n = 1,nmax
    call get_erec_imurec(logd,T,X,Y,erec,imu,derecdT,dimurecdT)
    if (d*erec>=eint) then ! avoid negative thermal energy
       T = 0.9*T; Tdot=0.;cycle
    endif
    corr = (eint-(radconst*T**3+1.5*Rg*d*imu)*T-d*erec) &
           / ( -4.*radconst*T**3-d*(1.5*Rg*(imu+dimurecdT*T)+derecdT) )
    if (-corr > 10.*T) then  ! do not let temperature guess increase by more than a factor of 1000
       T = 10.*T
       Tdot = 0.
    elseif (abs(corr) > W4err*T) then
       T = T + Tdot*dt
       Tdot = (1.-2.*dt)*Tdot - dt*corr
    else
       T = T-corr
       Tdot = 0.
    endif
    if (abs(corr)<max(eoserr*T,eoserr)) exit
    if (n>50) dt=0.5
 enddo
 if (n > nmax) then
    print*,'d=',d,'eint=',eint/d,'Tguess=',Tguess,'mu=',1./imu,'T=',T,'erec=',erec
    print*,'n = ',n,' nmax = ',n,' correction is ',abs(corr),' needs to be < ',eoserr*T
    call fatal('eos_gasradrec','Failed to converge on temperature in equationofstate_gasradrec')
 endif
 call get_erec_imurec(logd,T,X,Y,erec,imu)
 p = ( Rg*imu*d + radconst*T**3/3. )*T
 gamma_eff = 1.+p/(eint-d*erec)
 cf = sqrt(gamma_eff*p/d)

end subroutine equationofstate_gasradrec


!-----------------------------------------------------------------------
!+
!  Solve for u and T (and mu) from rho and P (ideal gas + radiation + recombination)
!  Inputs and outputs in cgs units
!+
!-----------------------------------------------------------------------
subroutine calc_uT_from_rhoP_gasradrec(rhoi,presi,X,Y,T,eni,mui,ierr)
 use ionization_mod, only: get_imurec,get_erec
 use physcon,        only: radconst,Rg
 real, intent(in)            :: rhoi,presi,X,Y
 real, intent(inout)         :: T
 real, intent(out)           :: eni
 real, optional, intent(out) :: mui
 integer, intent(out)        :: ierr
 integer                     :: n
 real                        :: logrhoi,imu,dimurecdT,dT,Tdot,corr

 if (T <= 0.) T = min((3.*presi/radconst)**0.25, presi/(rhoi*Rg))  ! initial guess for temperature
 ierr = 0
 corr = huge(0.)
 Tdot = 0.
 dT = 0.9
 logrhoi = log10(rhoi)
 do n = 1,500
    call get_imurec(logrhoi,T,X,Y,imu,dimurecdT)
    corr = ( presi - (rhoi*Rg*imu + radconst*T**3/3.)*T ) / &
           ( -rhoi*Rg * ( imu + T*dimurecdT ) - 4.*radconst*T**3/3. )
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
 if (n > 500) then
    print*,'Error in calc_uT_from_rhoP_hormone'
    print*,'rho=',rhoi,'P=',presi,'mu=',1./imu
    ierr = 1
    return
 endif
 eni = ( 1.5*Rg*imu + radconst*T**3/rhoi )*T + get_erec(logrhoi,T,X,Y)
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
subroutine read_options_eos_gasradrec(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*),  intent(in)  :: name,valstring
 logical,           intent(out) :: imatch,igotall
 integer,           intent(out) :: ierr
 integer,           save        :: ngot  = 0
 character(len=30), parameter   :: label = 'eos_gasradrec'

 imatch  = .true.
 select case(trim(name))
 case('irecomb')
    read(valstring,*,iostat=ierr) irecomb
    if ((irecomb < 0) .or. (irecomb > 3)) call fatal(label,'irecomb = 0,1,2,3')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

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
