!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cons2primsolver
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, metric_tools, utils_gr
!
 use eos, only:ieos,polyk
 implicit none

 public :: conservative2primitive,primitive2conservative

 private :: get_u

 integer, public, parameter :: &
      ien_etotal  = 1, &
      ien_entropy = 2


!!!!!!====================================================
!
!
! NOTE: cons2prim has been written for only adiabatic eos.
!
!
!!!!!!====================================================
 private

contains

!
! A few subroutines to compute stuff to do with the equation of state.
! They assume an adiabatic eos (ideal gas).
! These subroutines will need to be different for a different eos.
! (Should really exist in the eos module)
!
!=========================

pure subroutine get_u(u,P,dens,gamma)
 real, intent(in)  :: dens,P,gamma
 real, intent(out) :: u
 real :: uthermconst

 ! Needed in dust case when dens = 0 causes P/dens = NaN and therefore enth = NaN
 ! or gamma=1 gives divide-by-zero

 if (P==0. .or. abs(p)<tiny(p)) then
    u = 0.
 else
    u = (P/dens)/(gamma-1.)
 endif

 uthermconst = polyk
 if (ieos==4) u=uthermconst

end subroutine get_u

!----------------------------------------------------------------
!+
!  Construct conserved variables from the primitive variables
!  primitive variables are (v^i,d,u,P); v i
!  conserved variables are (rho,pmom_i,en)
!+
!----------------------------------------------------------------
subroutine primitive2conservative(x,metrici,v,dens,u,P,rho,pmom,en,ien_type)
 use utils_gr,     only:get_u0
 use metric_tools, only:unpack_metric
 use io,           only:error
 real, intent(in)  :: x(1:3),metrici(:,:,:)
 real, intent(in)  :: dens,v(1:3),u,P
 real, intent(out) :: rho,pmom(1:3),en
 integer, intent(in) :: ien_type
 real, dimension(0:3,0:3) :: gcov
 real :: sqrtg, enth, gvv, U0, v4U(0:3)
 real :: gam1
 integer :: i, mu, ierror

 v4U(0) = 1.
 v4U(1:3) = v(:)

 enth = 1. + u + P/dens

 ! Hard coded sqrtg=1 since phantom is always in cartesian coordinates
 sqrtg = 1.
 call unpack_metric(metrici,gcov=gcov)

 call get_u0(gcov,v,U0,ierror)
 if (ierror > 0) call error('get_u0 in prim2cons','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')
 rho = sqrtg*dens*U0
 do i=1,3
    pmom(i) = U0*enth*dot_product(gcov(i,:),v4U(:))
 enddo

 gvv = 0.
 do mu=0,3
    do i=1,3
       gvv = gvv + gcov(i,mu)*v4U(mu)*v4U(i)
    enddo
 enddo

 if (ien_type == ien_entropy) then
    if (u > 0) then
       gam1 = 1. + P/(dens*u)
       en = P/(dens**gam1)
    else
       ! handle the case for u = 0
       en = P/dens
    endif
 else
    en = U0*enth*gvv + (1.+u)/U0
 endif

end subroutine primitive2conservative

subroutine conservative2primitive(x,metrici,v,dens,u,P,rho,pmom,en,ierr,ien_type)
 use metric_tools, only:unpack_metric
 use units,        only:unit_ergg,unit_density,unit_pressure
 use eos,          only:calc_temp_and_ene,ieos
 real, intent(in)    :: x(1:3),metrici(:,:,:)
 real, intent(inout) :: dens,P,u
 real, intent(out)   :: v(1:3)
 real, intent(in)    :: rho,pmom(1:3),en
 integer, intent(out) :: ierr
 integer, intent(in)  :: ien_type
 real, dimension(1:3,1:3) :: gammaijUP
 real :: sqrtg,sqrtg_inv,enth,lorentz_LEO,pmom2,alpha,betadown(1:3),betaUP(1:3),enth_old,v3d(1:3)
 real :: f,df,term,lorentz_LEO2,gamfac,pm_dot_b,gamma,temp,sqrt_gamma_inv,gamma_old
 real :: ucgs,Pcgs,denscgs
 integer :: niter, i
 real, parameter :: tol = 1.e-12
 integer, parameter :: nitermax = 100
 logical :: converged
 ierr = 0

 !print*, 'in', dens, u, P, rho, en, pmom ! ppp
 ! Hard coding sqrgt=1 since phantom is always in cartesian coordinates
 sqrtg = 1.
 sqrtg_inv = 1./sqrtg

 ! Get metric components from metric array
 call unpack_metric(metrici,gammaijUP=gammaijUP,alpha=alpha,betadown=betadown,betaUP=betaUP)

 pmom2 = 0.
 do i=1,3
    pmom2 = pmom2 + pmom(i)*dot_product(gammaijUP(:,i),pmom(:))
 enddo

 ! Guess enthalpy (using previous values of dens and pressure)
 enth = 1. + u + P/dens
 if (u > tiny(0.)) then
    gamma = 1. + P/(dens*u)
 else
    gamma = 5./3. ! use gamma for ideal gas
 endif

 niter = 0
 converged = .false.
 sqrt_gamma_inv = alpha*sqrtg_inv ! get determinant of 3 spatial metric
 term = rho*sqrt_gamma_inv
 gamfac = gamma/(gamma-1.)
 pm_dot_b = dot_product(pmom,betaUP)

 do while (.not. converged .and. niter < nitermax)
    enth_old = enth
    !gamma_old = gamma
    lorentz_LEO2 = 1.+pmom2/enth_old**2
    lorentz_LEO = sqrt(lorentz_LEO2)
    dens = term/lorentz_LEO
    !print*, dens, term, lorentz_LEO

    if (ien_type == ien_entropy) then
       !select case(ieos)
       !case(12)
          !print*, 'before', gamma
          !call update_gamma(ien_type,dens,u,en,gamma,ierr)
          !print*, 'update', gamma
          !read*
       !end select
       !print*, 'PT', P
       p = en*dens**gamma
       !print*, 'P', P
    elseif (ieos==4) then
       p = (gamma-1.)*dens*polyk
    else
       p = max(rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-pm_dot_b),0.)
    endif
    !print*, 'iter', niter, 'P', p ! ppp

    enth = 0.
    if (p > 0.) then
      ucgs = u*unit_ergg
      Pcgs = P*unit_pressure
      denscgs = dens*unit_density
      call calc_temp_and_ene(denscgs,Pcgs,ucgs,temp,ierr,guesseint=ucgs)
      u = ucgs/unit_ergg
      enth = 1.+u+P/dens
    endif
    !call get_enthalpy(enth,dens,p,gamma)

    f = enth-enth_old

    !This line is unique to the equation of state - implemented for adiabatic at the moment
    if (ien_type == ien_entropy) then
       df = -1. + (gamma*pmom2*P)/(lorentz_LEO2 * enth_old**3 * dens)
    elseif (ieos==4) then
       df = -1. ! Isothermal, I think...
    else
       df= -1.+gamfac*(1.-pmom2*p/(enth_old**3*lorentz_LEO2*dens))
    endif

    enth = enth_old - f/df

    ! Needed in dust case when f/df = NaN casuses enth = NaN
    if (abs(enth_old-1.)<tiny(enth_old)) enth=1.

    niter = niter + 1

    if (abs(enth-enth_old)/enth < tol) converged = .true.
 enddo

 if (.not.converged) ierr = 1

 lorentz_LEO = sqrt(1.+pmom2/enth**2)
 dens = term/lorentz_LEO

 if (ien_type == ien_entropy) then
    p = en*dens**gamma
 else
    p = max(rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-dot_product(pmom,betaUP)),0.)
 endif

 v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)

! Raise index from down to up
 do i=1,3
    v(i) = dot_product(gammaijUP(:,i),v3d(:))
 enddo

 call get_u(u,P,dens,gamma)
 !print*, 'u', u, P, dens, gamma ! ppp
end subroutine conservative2primitive

subroutine update_gamma(ien_type,dens,u,en,gamma,ierr)
 use physcon,          only:kb_on_mh,radconst
 use eos,              only:gmw,ieos
 use eos_idealplusrad, only:get_idealplusrad_pres,get_idealplusrad_enfromtemp
 use units,            only:unit_pressure,unit_density,unit_ergg
 real, intent(in)     :: dens,en,u
 integer, intent(in)  :: ien_type
 real, intent(inout)  :: gamma
 integer, intent(out) :: ierr
 real, parameter      :: tol = 1.e-12
 real                 :: T,gamma1,gamma2,rhogamma,encgs,Pcgs,denscgs!,ucgs
 real                 :: dPdT,dudT,PT,uT,dgammadT
 real                 :: f,df,T_old
 integer              :: niter
 integer, parameter   :: nitermax = 100
 logical              :: converged
 real::pc,uc

 ! unit conversion
 Pcgs = en*dens**gamma*unit_pressure
 denscgs = dens*unit_density
 !ucgs = u*unit_ergg
 encgs = en*unit_ergg

 ! use pure rad temp as the initial guess
 T = (3.*Pcgs / radconst)**(1./4.)
 !print*, 'guess', T, pcgs, denscgs, gamma
 !read*
 converged = .false.
 niter = 0
 ierr = 0
 !print*, '                n      gamma                t_old,            t'

 do while (.not. converged .and. niter < nitermax)
    T_old = T

    select case(ieos)
    case(12)
      call get_idealplusrad_pres(denscgs,T,gmw,PT)
      !Pc = denscgs*kb_on_mh*T/gmw + 1./3.*radconst*T**4.
      dPdT = denscgs*kb_on_mh/gmw + 4./3.*radconst*T**3.
      call get_idealplusrad_enfromtemp(denscgs,T,gmw,uT)
      !uc = 3./2.*kb_on_mh*T/gmw + radconst*T**4./denscgs
      dudT = 3./2.*kb_on_mh/gmw + 4.*radconst*T**3./denscgs
      !print*, 'Pu', PT, Pc, ut, uc
    case(2)
      PT = denscgs*kb_on_mh*T/gmw
      dPdT = denscgs*kb_on_mh/gmw
      uT = 3./2.*kb_on_mh*T/gmw
      dudT = 3./2.*kb_on_mh/gmw
    end select

    gamma = 1. + PT/(denscgs*uT)
    gamma1 = gamma - 1
    dgammadT = 1/denscgs*(dPdT/uT - PT/uT**2*dudT)

    f = encgs/gamma1*denscgs**gamma1 - uT
    df = encgs/gamma1**2 * dgammadT * denscgs**gamma1 * (gamma1*log(denscgs) - 1) - dudT
    T = T - f/df

    print*, 'val',niter, gamma, T_old, T

    niter = niter + 1
    if (abs(T_old-T)/T < tol) converged = .true.
 enddo
 read*
 if (.not. converged) ierr = 1

end subroutine update_gamma

end module cons2primsolver
