!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_gasradrec
! ALi: I am changing this to to degenerate eos plus gasrad for C/O WDs based on Ryo's code
! in ryo's solver e is internal energy per volume and u is internal energy per gram and they convert with e=u*rho

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
           read_options_eos_gasradrec,write_options_eos_gasradrec,eos_info_gasradrec,init_eos_gasradrec,&
           get_pt_from_de_ryo, get_cs_from_dT_ryo, get_eT_from_dp_ryo, init_eos_ryo
 private
 real, parameter :: eoserr=1.e-15,W4err=1.e-2


 !!! constants from constants.f90, I might have to change this because the ones in eos_degen are based on the ones below:
 real(8),parameter:: pi=acos(-1d0), avo=6.02214076d23
 real(8),parameter:: clight=2.99792458d10, kbol=1.380649d-16
 real(8),parameter:: hplanck=6.62607015d-27,  hbar = hplanck/(2d0*pi)
 real(8),parameter:: m_p = 1.67262192369d-24, m_e=9.1093837015d-28, amu = 1d0/avo
 real(8),parameter:: G=6.67430d-8, sigma = (pi*pi*kbol*kbol*kbol*kbol)/(60*hbar*hbar*hbar*clight*clight) !5.67051d-5
 real(8),parameter:: arad = 4d0*sigma/clight
 real(8),parameter:: msun=1.3271244d26/6.67430d-8, rsun = 6.96d10, lsun=3.9d33

 real(8),parameter:: Knr=hplanck**2/(20d0*m_e*amu**(5d0/3d0))*(3d0/pi)**(2d0/3d0)
 real(8),parameter:: Ker=hplanck*clight/(8d0*amu**(4d0/3d0))*(3d0/pi)**(1d0/3d0)
 
!  real(8),parameter:: sigma = 5.67051d-5
! I removed Knr and Ker from constants because it seems like eos_degen wasn't using it.
 !!! from fermi_dirac_approx.f90
 real(8),parameter:: exp20=exp(20d0)
 real(8),parameter:: gamma12=gamma(0.5d0), gamma32=gamma(1.5d0), &
                     gamma52=gamma(2.5d0), gamma72=gamma(3.5d0), &
                     smooth32 = 2.5d0/sqrt(2d0), smooth52 = 2.5d0/sqrt(3d0)
 !! from eos_degenerate.f90
 real(8),parameter:: err = 1d-12
 real(8),parameter:: fac_ne = 8d0*pi*sqrt(2d0)*(m_e*clight /hplanck)**3
 real(8),parameter:: fac_pe = 16d0*pi*sqrt(2d0)/3d0/hplanck**3*m_e**4*clight**5
 real(8),parameter:: fac_ue = 8d0*pi*sqrt(2d0)/hplanck**3*m_e**4*clight**5
 real(8),parameter:: fac_theta = kbol/(m_e*clight**2)
 real(8),parameter:: fac_pf = 3d0*hplanck**3/(8d0*pi*(m_e*clight)**3)
!  real(8),parameter:: Knr=hplanck**2/(20d0*m_e)*(3d0/pi)**(2d0/3d0)
!  real(8),parameter:: Ker=hplanck*clight/8d0*(3d0/pi)**(1d0/3d0)

 real(8),parameter:: etalim=30d0, etamax=1d5, Tfloor=1d0/etamax
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

! ------------------------------------------------------------------------------

! analytical degenerate EOS by Ryo Hirai (which handles gasrad as well as fermi dirac integrals)

! -------------------------------------------------------------------------------

subroutine init_eos_ryo(ierr)
!   use eos, only: X_in, Z_in. ! this is bad architecture, so I am uncommenting it
!   real, intent(in) :: x, z
  integer, intent(out) :: ierr

  ierr = 0

!   ! sanity check (copy MESA style)
!   if ((x + z > 1.) .or. (x < 0.) .or. (z < 0.)) then
!      ierr = -1
!      return
!   endif

!   ! store them somewhere if needed (module variables)
!   X_in = x
!   Z_in = z

  write(*,*) 'Init degenerate EOS'
!   X=', x, ' Z=', z
end subroutine init_eos_ryo


 function get_eta_fully_ionized(d,mue,theta,eta_guess) result(eta)
!   use fermi_dirac_approx
  real(8),intent(in):: d,mue,theta
  real(8),intent(in),optional:: eta_guess
  real(8):: eta, Nele_matter, Nele, Npos, f, dfde, corr, dnede, dnpde, etapos
  real(8):: fe12, fe32, fp12, fp32, dfe12de, dfe32de, dfp12de, dfp32de
  real(8):: coeff, maxcorr = 2d0, thetalim, logne
  integer:: n, nmax=10000

  Nele_matter = d/(amu*mue)
  logne = log(Nele_matter)
  coeff = fac_ne*theta**1.5d0

  if(present(eta_guess))then
   eta = eta_guess
  else
   eta = Efermi(Nele_matter)*fac_theta/(kbol*theta)
  end if

  do n = 1, nmax

   call get_nele(eta,theta,Nele,dnede)
   if(theta>1d-5)then
    call get_nele(-eta-2d0/theta,theta,Npos,dnpde)
   else
    Npos = 0d0
    dnpde = 0d0
   end if

   f = 1d0-log(Nele)/logne
   dfde = - Nele/(dnede*logne)

   corr = - f/dfde
   if(eta<10d0)&
    corr = sign(min(abs(corr),max(abs(eta),1d-3)*0.5d0),corr)

   eta = eta + corr

   if(abs(f)<err*theta*1d4)exit

  end do

 end function get_eta_fully_ionized


 subroutine etaT_from_dp(d,p,mue,mui,eta,T)
  real(8),intent(in):: d,p,mue,mui
  real(8),intent(inout):: eta,T
  real(8):: f,dfdt,corr
  integer:: i,imax=10000
  real(8):: pp,dpdt,ne

  ne = d/(mue*amu)
  if(p <= Pdeg(ne)*(1.d0+epsilon(1d0)))then
   eta = etamax
   T = Efermi(ne)/(kbol*eta)
   return
  end if

! Get good initial guess if guess is not provided
  if(T<=0d0)call getT_from_dp(d,p,get_mu(mue,mui),T)

  do i = 1, imax

   call get_p_from_dT(d,T,mue,mui,pp,eta=eta,dpdt=dpdt)
   if(abs(pp/p-1d0)<err)exit

! Switch to bisection if it gets stuck in a limit cycle
   if(i>20)then
    call bisection_p(d,p,T,mue,mui,err)
    exit
   end if

   f = pp/p-1d0
   dfdt = dpdt/p

   corr = -f/dfdt
   corr = max(corr,-0.9d0*T)
   T = T + corr

  end do

  if(i>imax)then
   T = -1d0
   eta = 1d99
  end if

 end subroutine etaT_from_dp

 subroutine bisection_p(d,p,T,mue,mui,tol)
  real(8),intent(inout):: T
  real(8),intent(in):: d,p,mue,mui,tol
  real(8):: t1,t2,p1,p2,newt,newp,fac=0.5d0
  integer:: i

  call get_p_from_dT(d,T,mue,mui,p1);p1=p1-p

  if(p1 > 0)then
   t1 = T
   do i = 1, 1000
    t1 = fac*t1
    p2 = p1
    call get_p_from_dT(d,t1,mue,mui,p1);p1=p1-p
    if(p1<0)exit
   end do
   t2 = t1/fac
  else
   t2 = T
   do i = 1, 1000
    t2 = t2/fac
    p1 = p2
    call get_p_from_dT(d,t2,mue,mui,p2);p2=p2-p
    if(p2>0)exit
   end do
   t1 = t2*fac
  end if

  do i = 1, 1000
   newt = 0.5d0*(t1+t2)
   call get_p_from_dT(d,newt,mue,mui,newp);newp=newp-p
   if(p1*newp>0d0)then
    t1 = newt
    p1 = newp
   else
    t2 = newt
    p2 = newp
   end if
   if(abs(newp/p)<tol)exit
   if(t2/t1-1d0<err)exit
  end do

  T = newt

 end subroutine bisection_p

 subroutine etaT_from_de(d,e,mue,mui,eta,T)
  real(8),intent(in):: d,e,mue,mui
  real(8),intent(inout):: eta,T
  real(8):: f,dfdt,corr
  integer:: i,imax=1000
  real(8):: u,dudt,ne

  ne = d/(mue*amu)
  if(e <= Edeg(ne)*(1.d0+epsilon(1d0)))then
   eta = etamax
   T = Efermi(ne)/(kbol*eta)
   return
  end if

! Get good initial guess if guess is not provided
  if(T<=0d0)call getT_from_de(d,e,get_mu(mue,mui),T)

  do i = 1, imax

   call get_u_from_dT(d,T,mue,mui,u,dudt=dudt)
   if(abs(u*d/e-1d0)<err)exit

! Switch to bisection if it gets stuck in a limit cycle
   if(i>30)then
    call bisection_u(d,e/d,T,mue,mui,err)
    exit
   end if
  
   f = u*d/e-1d0
   dfdt = dudt*d/e

   corr = max(-f/dfdt,-0.9d0*T)
   T = T + corr

   if(abs(u*d/e-1d0)<err)exit
  end do

  call get_u_from_dT(d,T,mue,mui,u,eta=eta)

  if(i>imax)then
   T = -1d0
   eta = 1d99
  end if

 end subroutine etaT_from_de


 subroutine bisection_u(d,u,T,mue,mui,tol)
  real(8),intent(inout):: T
  real(8),intent(in):: d,u,mue,mui,tol
  real(8):: t1,t2,u1,u2,newt,newu,fac=0.5d0,eta=0d0
  integer:: i

  call get_u_from_dT(d,T,mue,mui,u1,eta=eta);u1=u1-u

  if(u1 > 0)then
   t1 = T
   do i = 1, 1000
    t1 = fac*t1
    u2 = u1
    call get_u_from_dT(d,t1,mue,mui,u1,eta=eta);u1=u1-u
    if(u1<0)exit
   end do
   t2 = t1/fac
  else
   t2 = T
   do i = 1, 1000
    t2 = t2/fac
    u1 = u2
    call get_u_from_dT(d,t2,mue,mui,u2,eta=eta);u2=u2-u
    if(u2>0)exit
   end do
   t1 = t2*fac
  end if

  do i = 1, 1000
   newt = 0.5d0*(t1+t2)
   call get_u_from_dT(d,newt,mue,mui,newu,eta=eta);newu=newu-u
   if(u1*newu>0d0)then
    t1 = newt
    u1 = newu
   else
    t2 = newt
    u2 = newu
   end if
   if(abs(newu/u)<tol)exit
   if(t2/t1-1d0<err)exit
  end do

  T = newt

 end subroutine bisection_u

 subroutine get_pT_from_de_ryo(d,e,mue,mui,p,T,eta)
  real(8),intent(in):: d,e,mue,mui
  real(8),intent(inout):: p,T
  real(8),intent(inout),optional:: eta
  real(8):: eta1

  call etaT_from_de(d,e,mue,mui,eta1,T)
  call get_p_from_dT(d,T,mue,mui,p,eta1)
  if(present(eta))eta=eta1

 end subroutine get_pT_from_de_ryo

 subroutine get_eT_from_dp_ryo(d,p,mue,mui,e,T,eta)
  real(8),intent(in):: d,p,mue,mui
  real(8),intent(inout):: e,T
  real(8),intent(inout),optional:: eta
  real(8):: eta1,u

  call etaT_from_dp(d,p,mue,mui,eta1,T)
  call get_u_from_dT(d,T,mue,mui,u,eta1)
  e = u*d
  if(present(eta))eta=eta1

 end subroutine get_eT_from_dp_ryo

 subroutine get_u_from_dT(d,T,mue,mui,u,eta,dudt)
  real(8),intent(in):: d,T,mue,mui
  real(8),intent(inout):: u
  real(8),intent(inout),optional:: eta
  real(8),intent(out),optional:: dudt
  real(8):: theta,eta1,eele,deeledt,ne

  ne = d/(mue*amu)
  if(T<=Tfloor*Tfermi(ne))then
   theta = fac_theta*Tfloor*Tfermi(ne)
   eta1 = 1d0/Tfloor
  else   
   theta = fac_theta*T
   eta1 = get_eta_fully_ionized(d,mue,theta)
  end if

  if(present(dudt))then
   call get_eele(eta1,theta,eele,nele=ne,deeledt=deeledt)
  else
   call get_eele(eta1,theta,eele,nele=ne)
  end if
  u = (eele + eion(d,T,mui) + erad(T))/d

  if(present(eta))eta = eta1
  if(present(dudt))dudt = (deeledt+deiondT(d,T,mui)+deraddT(T))/d

 end subroutine get_u_from_dT

 subroutine get_p_from_dT(d,T,mue,mui,p,eta,dpdt,dpdd)
  real(8),intent(in):: d,T,mue,mui
  real(8),intent(inout):: p
  real(8),intent(inout),optional:: eta
  real(8),intent(out),optional:: dpdt,dpdd
  real(8):: theta,eta1,pele,dpeledt,dpeledd,ne

  ne = d/(mue*amu)
  if(T<=Tfloor*Tfermi(ne))then
   theta = fac_theta*Tfloor*Tfermi(ne)
   eta1 = 1d0/Tfloor
  else
   theta = fac_theta*T
   eta1 = get_eta_fully_ionized(d,mue,theta)
  end if

  if(present(dpdt).and.present(dpdd))then
   call get_pele(eta1,theta,pele,nele=ne,dpeledt=dpeledt,dpeledne=dpeledd)
  elseif(present(dpdt))then
   call get_pele(eta1,theta,pele,nele=ne,dpeledt=dpeledt)
  elseif(present(dpdd))then
   call get_pele(eta1,theta,pele,nele=ne,dpeledne=dpeledd)
  else
   call get_pele(eta1,theta,pele,nele=ne)
  end if
  if(present(dpdd))dpeledd = dpeledd/(mue*amu)

  p = pele + pion(d,T,mui) + prad(T)

  if(present(eta))eta = eta1
  if(present(dpdt))dpdt = dpeledt+dpiondT(d,T,mui)+dpraddT(T)
  if(present(dpdd))dpdd = dpeledd+pion(d,T,mui)/d

 end subroutine get_p_from_dT

 subroutine get_cs_from_dT_ryo(d,T,mue,mui,cs,gamma1,eta)
  real(8),intent(in):: d,T,mue,mui
  real(8),intent(out):: cs,gamma1
  real(8),intent(inout),optional:: eta
  real(8):: eta1,u,p,dpdt,dpdd,dudt

  if(present(eta))then
   eta1 = eta
  else
   eta1 = 0d0
  end if

  call get_p_from_dT(d,T,mue,mui,p,eta=eta1,dpdt=dpdt,dpdd=dpdd)
  call get_u_from_dT(d,T,mue,mui,u,eta=eta1,dudt=dudt)

  gamma1 = d/p*dpdd + T/(d*p*dudt)*dpdt**2
  cs = sqrt(gamma1*p/d)
  if(present(eta))eta=eta1
  
 end subroutine get_cs_from_dT_ryo

 subroutine get_nele(eta,theta,nele,dnelede)
!   use fermi_dirac_approx
  implicit none
  real(8),intent(in):: eta,theta
  real(8),intent(out):: nele
  real(8),intent(out),optional:: dnelede
  real(8):: theta12,theta32,theta52
  real(8):: f12,f32,f52,df12de,df32de,df52de
  real(8):: deeledt_eta,deeledeta_t,dnedt_eta,dnedeta_t

  theta12 = sqrt(theta)
  theta32 = theta12*theta
  theta52 = theta32*theta

  if(present(dnelede))then
   call dfermi(0.5d0,eta,theta,f12,fdeta=df12de)
   call dfermi(1.5d0,eta,theta,f32,fdeta=df32de)
  else
   call dfermi(1.5d0,eta,theta,f12)
   call dfermi(2.5d0,eta,theta,f32)
  end if

  nele = fac_ne*theta32*(f12+theta*f32)

  if(present(dnelede))then
   dnelede = fac_ne*theta32*(df12de+theta*df32de)
  end if

 end subroutine get_nele

 subroutine get_pele(eta,theta,pele,nele,dpeledt,dpeledne)
!   use fermi_dirac_approx
  implicit none
  real(8),intent(in):: eta,theta
  real(8),intent(out):: pele
  real(8),intent(in),optional:: nele
  real(8),intent(out),optional:: dpeledt,dpeledne
  real(8):: theta12,theta32,theta52
  real(8):: f12,f32,f52,df12de,df32de,df52de,df12dt,df32dt,df52dt
  real(8):: eta1,theta1,ne,pcold, dpeledt_eta,dpeledeta_t,dnedt_eta,dnedeta_t
  real(8):: dpeledne_cold

! Return zero temperature value if theta is small
  if(present(nele))then ! avoid expensive nele calculation if it is provided
   ne=nele
  else ! if not provided, calculate electron number density
   call get_nele(eta,theta,ne)
  end if

!  if(theta<=Tfloor*Tfermi(ne)*fac_theta)then
  if(eta>etamax)then
   pele = Pdeg(ne)
   if(present(dpeledt ))dpeledt  = 0d0
   if(present(dpeledne))then
    dpeledne_cold = (pfermi(ne)*clight)**2/(3d0*(Efermi(ne)+m_e*clight*clight))
    dpeledne = dpeledne_cold
   end if

   return
  end if

! Limit eta value to <etalim and use cold formula instead
  if(eta>etalim)then
   eta1 = etalim
   theta1 = fac_theta*Tfermi(ne)/eta1
  else
   eta1 = eta
   theta1 = theta
  end if

  theta12 = sqrt(theta1)
  theta32 = theta12*theta1
  theta52 = theta32*theta1

  if(present(dpeledt))then
   call dfermi(0.5d0,eta1,theta1,f12,fdeta=df12de,fdtheta=df12dt)
   call dfermi(1.5d0,eta1,theta1,f32,fdeta=df32de,fdtheta=df32dt)
   call dfermi(2.5d0,eta1,theta1,f52,fdeta=df52de,fdtheta=df52dt)
  else
   call dfermi(1.5d0,eta1,theta1,f32)
   call dfermi(2.5d0,eta1,theta1,f52)
  end if

  pele = fac_pe*theta52*(f32+0.5d0*theta1*f52)

  if(eta>etalim)then
   pcold = pdeg(ne)
   pele = (pele-pcold)*etalim/eta + pcold
  end if

  if(present(dpeledt).or.present(dpeledne))then
   dpeledeta_t = fac_pe*theta52*(df32de+0.5d0*theta1*df52de)
   dnedeta_t   = fac_ne*theta32*(df12de+theta1*df32de)

   if(present(dpeledt))then
    dpeledt_eta = fac_pe &
                 *( 2.5d0*theta32*(f32+0.5d0*theta1*f52) &
                 + theta52*(df32dt+0.5d0*f52+0.5d0*theta1*df52dt))
    dnedt_eta = fac_ne*( 1.5d0*theta12*(f12+theta1*f32) &
               + theta32*(df12dt+f32+theta1*df32dt) )
    dpeledt = dpeledt_eta - dpeledeta_t * dnedt_eta/dnedeta_t
    if(eta>etalim)dpeledt = dpeledt*etalim/eta
    dpeledt = dpeledt*fac_theta
   end if

   if(present(dpeledne))then
    dpeledne = dpeledeta_t / dnedeta_t
!    if(eta>etalim)dpeledne = (dpeledne-dpeledne_cold)*etalim/eta + dpeledne_cold
   end if
  end if

 end subroutine get_pele

 subroutine get_eele(eta,theta,eele,nele,deeledt)
!   use fermi_dirac_approx
  real(8),intent(in):: eta,theta
  real(8),intent(out):: eele
  real(8),intent(in),optional:: nele
  real(8),intent(out),optional:: deeledt
  real(8):: theta12,theta32,theta52
  real(8):: f12,f32,f52,df12de,df32de,df52de,df12dt,df32dt,df52dt
  real(8):: eta1,theta1,ecold,ne, deeledt_eta,deeledeta_t,dnedt_eta,dnedeta_t

! Return zero temperature value if theta is small
  if(present(nele))then
   ne = nele
  else
   call get_nele(eta,theta,ne)
  end if

! Limit eta value to <etalim and use cold formula instead
  if(eta>etalim)then
   eta1 = etalim
   theta1 = fac_theta*Tfermi(ne)/eta1
  else
   eta1 = eta
   theta1 = theta
  end if

  theta12 = sqrt(theta1)
  theta32 = theta12*theta1
  theta52 = theta32*theta1

  if(present(deeledt))then
   call dfermi(0.5d0,eta1,theta1,f12,fdeta=df12de,fdtheta=df12dt)
   call dfermi(1.5d0,eta1,theta1,f32,fdeta=df32de,fdtheta=df32dt)
   call dfermi(2.5d0,eta1,theta1,f52,fdeta=df52de,fdtheta=df52dt)
  else
   call dfermi(1.5d0,eta1,theta1,f32)
   call dfermi(2.5d0,eta1,theta1,f52)
  end if

  eele = fac_ue*theta52*(f32+theta1*f52)

  if(eta>etalim)then
   ecold = Edeg(ne)
   eele = (eele-ecold)*etalim/eta + ecold
  end if

  if(present(deeledt))then
   dnedt_eta = ( 1.5d0*theta12*(f12+theta1*f32) &
              + theta32*(df12dt+f32+theta1*df32dt) )
   dnedeta_t = theta32*( df12de+theta1*df32de )

   deeledt_eta = fac_ue &
               *( 2.5d0*theta32*(f32+theta1*f52) &
                + theta52*(df32dt+f52+theta1*df52dt))
   deeledeta_t = fac_ue*theta52*(df32de+theta1*df52de)

   deeledt = deeledt_eta - deeledeta_t * dnedt_eta/dnedeta_t
   if(eta>etalim)deeledt = deeledt*etalim/eta
   deeledt = deeledt*fac_theta
  end if

 end subroutine get_eele

 function pion(d,T,mui)
  real(8),intent(in):: d,T,mui
  real(8):: pion
  pion = d*kbol*T/(amu*mui)
 end function pion

 function dpiondT(d,T,mui)
  real(8),intent(in):: d,T,mui
  real(8):: dpiondT
  dpiondT = d*kbol/(amu*mui)
 end function dpiondT

 function eion(d,T,mui)
  real(8),intent(in):: d,T,mui
  real(8):: eion
  eion = 1.5d0*d*kbol*T/(amu*mui)
 end function eion

 function deiondT(d,T,mui)
  real(8),intent(in):: d,T,mui
  real(8):: deiondT
  deiondT = 1.5d0*d*kbol/(amu*mui)
 end function deiondT

 function prad(T)
  real(8),intent(in):: T
  real(8):: prad
  prad = arad*T**4/3d0
 end function prad

 function dpraddT(T)
  real(8),intent(in):: T
  real(8):: dpraddT
  dpraddT = 4d0*arad*T**3/3d0
 end function dpraddT

 function erad(T)
  real(8),intent(in):: T
  real(8):: erad
  erad = arad*T**4
 end function erad

 function deraddT(T)
  real(8),intent(in):: T
  real(8):: deraddT
  deraddT = 4d0*arad*T**3
 end function deraddT

 function Tfermi(ne) result(T)
  ! Fermi temperature
  real(8),intent(in):: ne
  real(8):: T
  T = Efermi(ne)/kbol
 end function Tfermi

 function Efermi(ne) result(Ef)
  ! Fermi energy without rest mass
  real(8),intent(in):: ne
  real(8):: Ef, f, dfde, corr, m_e_c2, pf_c2
!  Ef = hbar*hbar/(2d0*m_e)*(3d0*pi*pi*ne)**(2d0/3d0)
  m_e_c2 = m_e*clight*clight
  pf_c2 = (pfermi(ne)*clight)**2
  Ef = max(sqrt(pf_c2+m_e_c2**2)-m_e_c2,1d-50)
  do
   f = Ef*Ef+2d0*Ef*m_e_c2-pf_c2
   dfde = 2d0*(Ef+m_e_c2)
   corr = -f/dfde
   Ef = Ef + corr
   if(abs(f/pf_c2)<err)exit
  end do
 end function Efermi

 function pfermi(ne) result(pf)
  ! Fermi momentum
  real(8),intent(in):: ne
  real(8):: pf
  pf = (3d0*pi*pi*ne)**(1d0/3d0)*hbar
 end function pfermi

 function Pdeg(ne) result(P)
  real(8),intent(in):: ne
  real(8):: p, x, sqrtx2p1
! Avoid tedious calculation for extremes
  if(ne>1d34)then
   P = Pdeg_er(ne)
   return
  elseif(ne<1d24)then
   P = Pdeg_nr(ne)
   return
  end if
  x = (fac_pf*ne)**(1d0/3d0)
  sqrtx2p1 = sqrt(1d0+x*x)
  P = pi*(m_e*clight)**4*clight/(3d0*hplanck*hplanck*hplanck)*(x*(2d0*x*x-3d0)*sqrtx2p1+3d0*log(x+sqrtx2p1))
 end function Pdeg

 function Edeg(ne) result(E)
  real(8),intent(in):: ne
  real(8):: e, x, sqrtx2p1
  ! Avoid tedious calculation for extremes
  if(ne>1d34)then
   E = 3d0*Pdeg_er(ne)
   return
  elseif(ne<1d24)then
   E = 1.5d0*Pdeg_nr(ne)
   return
  end if
  x = (fac_pf*ne)**(1d0/3d0)
  sqrtx2p1 = sqrt(1d0+x*x)
  E = pi*(m_e*clight)**4*clight/(3d0*hplanck*hplanck*hplanck)*(8d0*x*x*x*(sqrtx2p1-1d0)-&
  (x*(2d0*x*x-3d0)*sqrtx2p1+3d0*log(x+sqrtx2p1)))
 end function Edeg

 function rhodeg(p,mue) result(d)
!   use constants
  real(8),intent(in):: p,mue
  real(8):: d, t, corr, pnorm
! Avoid tedious calculation for extremes
  if(p>1d28)then
   d = rhodeg_er(p,mue)
   return
  elseif(p<1d17)then
   d = rhodeg_nr(p,mue)
   return
  end if
  pnorm = p*12d0*hplanck**3/(pi*m_e**4*clight**5)
  corr = 1d99;t=100d0
  do while (abs(corr)>t*1d-10)
   corr = (sinh(t)-8d0*sinh(0.5d0*t)+3d0*t-pnorm)/(cosh(t)-4d0*cosh(0.5d0*t)+3d0)
   t = t - corr
   t=abs(t)
  end do
  d = mue*amu*8d0*pi*m_e**3*clight**3/(3d0*hplanck**3)*sinh(t*0.25d0)**3
 end function rhodeg 

 function Pdeg_nr(ne) result(P)
! non-relativistic electron degeneracy pressure
  real(8),intent(in):: ne
  real(8):: P
  P = Knr* ne**(5d0/3d0)
 end function Pdeg_nr

 function rhodeg_nr(P,mue) result(rho)
! non-relativistic electron degeneracy density
  real(8),intent(in):: P,mue
  real(8):: rho
  rho = (P/Knr)**0.6d0 * mue
 end function rhodeg_nr

 function Pdeg_er(ne) result(P)
! extra-relativistic electron degeneracy pressure
  real(8),intent(in):: ne
  real(8):: P
  P = Ker* ne**(4d0/3d0)
 end function Pdeg_er

 function rhodeg_er(P,mue) result(rho)
! extra-relativistic electron degeneracy density
  real(8),intent(in):: P,mue
  real(8):: rho
  rho = (P/Ker)**0.75d0 * mue
 end function rhodeg_er

 function entropy_full(d,T,eta,x,A) result(S)
  real(8),intent(in):: d,T,eta
  real(8),allocatable,intent(in):: x(:),A(:)
  real(8):: S
  S = entropy_ele(d,T,eta) + entropy_ion(d,T,x,A) !+ entropy_rad(d,T)
 end function entropy_full


 function entropy_ele(d,T,eta) result(S)
  real(8),intent(in):: d,T,eta
  real(8):: S, ue,pe,ne,theta
  theta = T*fac_theta
  call get_nele(eta,theta,ne)
  call get_pele(eta,theta,pe)
  call get_eele(eta,theta,ue)
  S = ( (ue+pe)/T - eta*ne*kbol )/d
 end function entropy_ele

 function entropy_ion(d,T,x,A) result(S)
  real(8),intent(in):: d,T
  real(8),allocatable,intent(in):: x(:),A(:)
  real(8):: S,mui, ni, W, eta_ion
  integer:: j,jmax

!  mui = get_mui(x,A)
!  ni = d / (amu*mui)

  jmax = size(x)
  W = 0d0; S = 0d0
  do j = 1, jmax
   if(x(j)<1d-3)cycle
   ni = d*x(j)/(A(j)*amu)
   W = (sqrt(2d0*pi*A(j)*amu*kbol*T)/hplanck)**3/ni
   S = S + x(j)/A(j)*(2.5d0+log(W))
  end do
  S = kbol/amu*S
!!$
!!$  eta_ion = log((sqrt(2d0*pi*mui*amu*kbol*T)/hplanck)**3/ni)
!!$  S = (pion(d,T,mui)+eion(d,T,mui))/(d*T)-eta_ion*kbol*ni/d

 end function entropy_ion

 function entropy_rad(d,T) result(S)
  real(8),intent(in):: d,T
  real(8)::S
  S = 4d0*arad/(3d0*d)*T**3
 end function entropy_rad


 function get_mu(mue,mui) result(mu)
  real(8),intent(in):: mue,mui
  real(8):: mu
  mu = mue*mui/(mue+mui)
 end function get_mu

 function get_mui(x,A) result(mui)
  real(8),allocatable,intent(in):: x(:),A(:)
  real(8):: mui,imui
  integer:: j, jmax
! A couple of checks
  if(abs(sum(x)-1d0)> 10d0*epsilon(1d0))then
   print*,'In get_mui: Mass fractions do not add up to 1!'
   print*,'sum(x_i)-1=',sum(x)-1d0
   stop
  end if
  if(size(x)/=size(A))then
   print*,'In get_mui: Make sure mass fraction array and atomic mass arrays are the same size'
   print*,'size(x_i),size(A_i)=',size(x),size(A)
   stop
  end if

  jmax = size(x)
  imui = 0d0
  do j = 1, jmax
   imui = imui + x(j)/A(j)
  end do
  mui = 1d0/imui
 end function get_mui

 function get_mue(x,A,Z) result(mue)
  real(8),allocatable,intent(in):: x(:),A(:),Z(:)
  real(8):: mue,imue
  integer:: j, jmax
! A couple of checks
  if(abs(sum(x)-1d0)> 10d0*epsilon(1d0))then
   print*,'In get_mue: Mass fractions do not add up to 1!'
   print*,'sum(x_i)-1=',sum(x)-1d0,epsilon(1d0)
   stop
  end if
  if(size(x)/=size(A))then
   print*,'In get_mui: Make sure mass fraction array and atomic mass arrays are the same size'
   print*,'size(x_i),size(A_i)=',size(x),size(A)
   stop
  end if

  jmax = size(x)
  imue = 0d0
  do j = 1, jmax
   imue = imue + x(j)*Z(j)/A(j)
  end do
  mue = 1d0/imue
 end function get_mue

 subroutine getT_from_dp(d,p,mu,T)
! PURPOSE: To calculate temperature from density and pressure for gasrad EoS
  implicit none
  real(8),intent(in):: d,p,mu
  real(8),intent(inout):: T
  real(8):: corr, R_on_mu
  real(8),parameter:: eoserr=1d-15
  integer:: n

  R_on_mu = kbol/amu/mu
  if(T<=0d0) T=1d3
 
  corr = 1d99
  do n = 1, 500
   corr = (p - (arad*T**3/3d0+d*R_on_mu)*T) &
        / (-4d0*arad*T**3/3d0-d*R_on_mu)
   T = T - corr
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>500)then
   write(*,*) 'Error: getT_from_dP did not converge after 500 iterations.'
   print*,'d=',d,'p=',p,'mu=',mu
   stop
  end if

 end subroutine getT_from_dp

 subroutine getT_from_de(d,e,mu,T)
! PURPOSE: To calculate temperature from density and internal energy for gasrad EoS
 implicit none
 real(8),intent(in):: d,e,mu
 real(8),intent(inout):: T
 real(8):: corr, R_on_mu
 real(8),parameter:: eoserr=1d-15
 integer:: n

 R_on_mu = kbol/amu/mu
 if(T<=0d0) T=1d3
 
  corr = 1d99
  do n = 1, 500
   corr = (e - (arad*T**3+1.5d0*d*R_on_mu)*T) &
        / (-4d0*arad*T**3-1.5d0*d*R_on_mu)
   T = T - corr
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>500)then
   write(*,*) 'Error: getT_from_de did not converge after 500 iterations.'
   print*,'d=',d,'e=',e,'mu=',mu
   stop
  end if

 end subroutine getT_from_de


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! fermi_dirac_approx.f90 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine dfermi(k, eta, theta, f, fdeta, fdtheta)
  real(8), intent(in) :: k, eta, theta
  real(8), intent(out) :: f
  real(8), intent(out), optional :: fdeta, fdtheta
  integer:: n
  real(8):: f12, fk, kratio, thratio, detaratio, dthetaratio, tail, bulk
  real(8):: alphak, thtail, thtail_num, thtail_den, thbulk, thsmooth
  real(8):: detabulk, xi0, dthetatail, dthetabulk, weight, dthsmooth
  real(8):: w1, w2, u_prime, term_bridge, term_physics

  f12 = fermi_12(eta)

  select case(nint(k+0.5d0))
  case(1)! k=1/2
   fk = f12
  case(2)! k=3/2
   tail = k
   bulk = k/(k+1d0)*xi(eta,12d0*k)
   kratio = blend(tail,bulk,smooth32)
   fk = f12*kratio
  case(3)! k=5/2
   tail = k
   bulk = k/(k+1d0)*xi(eta,12d0*k)
   kratio = blend(tail,bulk,smooth52)
   tail = k-1d0
   bulk = tail/k*xi(eta,12d0*tail)
   kratio = kratio*blend(tail,bulk,smooth32)
   fk = f12*kratio
  case default
   print*, 'This k value for the Fermi-Dirac integral is not supported'
   print'(a,f4.1)', 'k = ',k
  end select

! Apply thermal correction
  thtail_num = 1d0+(k+1.5d0)*theta+(k+1d0)*(4d0*k+7d0)/16d0*theta**2
  thtail_den = 1d0+(0.5d0*k+1d0)*theta
  thtail = sqrt( thtail_num / thtail_den )
  alphak = 0.5d0*((k+1d0)/(k+1.5d0))**2
  thbulk = sqrt(1d0+alphak*theta*xi(eta,20d0*k))

  thsmooth = 10d0/max(theta,0.02d0)
  thratio = blend(thtail,thbulk,thsmooth,exp(thsmooth))
  f = fk * thratio

  if(present(fdeta))then
   xi0 = xi(eta,10d0*k)
!   detabulk = (1d0/(k+1d0)+alphak/(2d0*(1d0+alphak*xi0)))*xi0
   if(abs(xi0) < epsilon(1d0))then
    detabulk = 0d0
   else
    detabulk = 1d0/((k+1d0)/xi0+alphak*theta/(2d0*(1d0+alphak*theta*xi0)))
   end if
   detaratio = 1d0/blend(1d0,detabulk,2d0,1d0)

   fdeta = f*detaratio
  end if

  if(present(fdtheta))then
    if (theta < 0.02d0) then
       dthsmooth = 0d0
    else
       dthsmooth = -10d0/theta**2
    end if

    dthetatail = 0.5d0*thtail &
               * ( (k+1.5d0+(k+1d0)*(4d0*k+7d0)/8d0*theta)/thtail_num &
                 - (0.5d0*k+1d0)/thtail_den )
    dthetabulk = alphak*xi(eta, 20d0*k)/(2d0*thbulk)

    w1 = exp(thsmooth * (thtail-1d0))
    w2 = exp(thsmooth * (thbulk-1d0))
    u_prime = w1 + w2 - 1.0d0

    term_bridge = (dthsmooth / thsmooth) * ((thtail*w1 + thbulk*w2 - 1.0d0)/u_prime - thratio)
    term_physics = (dthetatail*w1 + dthetabulk*w2) / u_prime

    dthetaratio = (term_bridge + term_physics) / thratio
    fdtheta = f * dthetaratio
  end if
  
 end subroutine dfermi

! Use F_{1/2}(eta,0) as base value
 function fermi_12(eta) result(f)
  real(8),intent(in):: eta
  real(8):: f, expeta

  if(eta>50d0)then
   f = sqrt(eta)**3/1.5d0
  else
   expeta = exp(max(eta,-100d0))
   f = 1d0/(1d0/(gamma32*expeta)+1.5d0/sqrt(log(25d0+expeta))**3)
  end if

 end function fermi_12

 pure function xi(x,offset)
  real(8),intent(in):: x,offset
  real(8):: xi
  xi = 0.5d0*(x+sqrt(x**2+offset))
 end function xi

 function blend(f1,f2,smooth,damp)
  real(8),intent(in):: f1,f2,smooth
  real(8),intent(in),optional:: damp
  real(8):: blend, m

  if((max(f1,f2)-min(f1,f2))*smooth>5d0)then
   blend = max(f1,f2)
   return
  end if

  m = max(f1, f2)
  if(present(damp)) then
   blend = m + log(exp(smooth*(f1-m))+exp(smooth*(f2-m))-damp/exp(smooth*m)) &
              /smooth
  else
   blend = m + log(exp(smooth*(f1-m))+exp(smooth*(f2-m)))/smooth
  end if

 end function blend

end module eos_gasradrec
