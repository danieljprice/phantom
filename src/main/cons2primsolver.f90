!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cons2primsolver
!
! None
!
! :References: None
!
! :Owner: Fitz) Hu
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, metric_tools, units, utils_gr
!
 use eos, only:ieos,polyk
 use options, only:ien_etotal,ien_entropy
 implicit none

 public :: conservative2primitive,primitive2conservative

 private :: get_u,conservative2primitive_con_gamma,conservative2primitive_var_gamma


!!!!!!====================================================
!
!
! NOTE: cons2prim has been written for only adiabatic & idealplusrad eos.
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

 if (ien_type == ien_etotal) then
    en = U0*enth*gvv + (1.+u)/U0
 else
    if (u > 0) then
       gam1 = 1. + P/(dens*u)
       en = P/(dens**gam1)
    else
       ! handle the case for u = 0
       en = P/dens
    endif
 endif

end subroutine primitive2conservative

!----------------------------------------------------------------
!+
!  solve for primitive variables from the conseerved variables
!  primitive variables are (v^i,d,u,P); v i
!  conserved variables are (rho,pmom_i,en)
!+
!----------------------------------------------------------------
subroutine conservative2primitive(x,metrici,v,dens,u,P,rho,pmom,en,ierr,ien_type)
 use eos,     only:ieos,gamma
 use io,      only:fatal
 real, intent(in)     :: x(1:3),metrici(:,:,:)
 real, intent(inout)  :: dens,P,u
 real, intent(out)    :: v(1:3)
 real, intent(in)     :: rho,pmom(1:3),en
 integer, intent(out) :: ierr
 integer, intent(in)  :: ien_type
 real                 :: enth

 select case (ieos)
 case (12)
    if (ien_type == ien_entropy) then
       call fatal('cons2primsolver','gasplusrad (ieos=12) only works with ien_type=ien_etotal for the moment')
    endif
    call conservative2primitive_var_gamma(x,metrici,v,dens,u,P,rho,pmom,en,ierr,ien_type)
 case default
    call conservative2primitive_con_gamma(x,metrici,v,dens,u,P,gamma,enth,rho,pmom,en,ierr,ien_type)
 end select

end subroutine conservative2primitive

!----------------------------------------------------------------
!+
!  solve for primitive variables from the conserved variables
!  for equations of state where gamma is an OUTPUT
!+
!----------------------------------------------------------------
subroutine conservative2primitive_var_gamma(x,metrici,v,dens,u,P,rho,pmom,en,ierr,ien_type)
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
 real :: f,term,lorentz_LEO2,gamfac,pm_dot_b,gamma,gamma_old,temp,sqrt_gamma_inv
 real :: u_in,P_in,dens_in,ucgs,Pcgs,denscgs,enth0,gamma0,enth_min,enth_max
 real :: enth_rad,enth_gas,gamma_rad,gamma_gas
 integer :: niter,i,ierr1,ierr2
 real, parameter :: tol = 1.e-12
 integer, parameter :: nitermax = 500
 logical :: converged
 ierr = 0

 ! Hard coding sqrgt=1 since phantom is always in cartesian coordinates
 sqrtg = 1.
 sqrtg_inv = 1./sqrtg

 ! Get metric components from metric array
 call unpack_metric(metrici,gammaijUP=gammaijUP,alpha=alpha,betadown=betadown,betaUP=betaUP)

 pmom2 = 0.
 do i=1,3
    pmom2 = pmom2 + pmom(i)*dot_product(gammaijUP(:,i),pmom(:))
 enddo

 niter = 0
 converged = .false.
 sqrt_gamma_inv = alpha*sqrtg_inv ! get determinant of 3 spatial metric
 term = rho*sqrt_gamma_inv
 pm_dot_b = dot_product(pmom,betaUP)

 ! Guess gamma (using previous values of dens and pressure)
 if (u > tiny(0.)) then
    gamma = 1. + P/(dens*u)
 else
    gamma = 5./3. ! use gamma for ideal gas
 endif

 ! NR with const gamma to get boundary value
 P_in = P
 dens_in = dens
 u_in = u

 gamma_gas = 5./3.
 call conservative2primitive_con_gamma(x,metrici,v,dens_in,u_in,P_in,gamma_gas,enth_gas,rho,pmom,en,ierr1,ien_type)

 P_in = P
 dens_in = dens
 u_in = u

 gamma_rad = 4./3.
 call conservative2primitive_con_gamma(x,metrici,v,dens_in,u_in,P_in,gamma_rad,enth_rad,rho,pmom,en,ierr2,ien_type)

 enth_min = enth_rad
 enth_max = enth_gas
 if (enth_min < 1.) enth_min = (enth_max-1.) / 1.6 + 1.

 enth = 0.5*(enth_max+enth_min)

 gamfac = gamma/(gamma-1.)
 enth0 = enth
 gamma0 = gamma

 do while (.not. converged .and. niter < nitermax)
    enth_old = enth
    gamma_old = gamma
    lorentz_LEO2 = 1.+pmom2/enth_old**2
    lorentz_LEO = sqrt(lorentz_LEO2)
    dens = term/lorentz_LEO

    p = max(rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-pm_dot_b),0.)

    if (P > 0.) then
       ucgs = u*unit_ergg
       Pcgs = P*unit_pressure
       denscgs = dens*unit_density

       call calc_temp_and_ene(ieos,denscgs,Pcgs,ucgs,temp,ierr,guesseint=ucgs)
       if (ierr /= 0) stop 'Did not converge'
       u = ucgs/unit_ergg

       gamma = 1. + P/(u*dens)
       gamfac = gamma/(gamma-1.)
    endif

    f = 1. + u + P/dens - enth_old

    if (f < 0) then
       enth_min = enth_old
    else
       enth_max = enth_old
    endif

    enth = 0.5*(enth_min + enth_max)

    ! Needed in dust case when f/df = NaN casuses enth = NaN
    if (abs(enth_old-1.)<tiny(enth_old)) enth=1.

    niter = niter + 1

    if (abs(enth-enth_old)/enth0 < tol) converged = .true.
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

end subroutine conservative2primitive_var_gamma

!----------------------------------------------------------------
!+
!  solve for primitive variables from the conserved variables
!  for equations of state where gamma is constant
!+
!----------------------------------------------------------------
subroutine conservative2primitive_con_gamma(x,metrici,v,dens,u,P,gamma,enth,rho,pmom,en,ierr,ien_type)
 use metric_tools, only:unpack_metric
 use eos,          only:calc_temp_and_ene,ieos
 real, intent(in)    :: x(1:3),metrici(:,:,:),gamma
 real, intent(inout) :: dens,P,u
 real, intent(out)   :: enth,v(1:3)
 real, intent(in)    :: rho,pmom(1:3),en
 integer, intent(out) :: ierr
 integer, intent(in)  :: ien_type
 real, dimension(1:3,1:3) :: gammaijUP
 real :: sqrtg,sqrtg_inv,lorentz_LEO,pmom2,alpha,betadown(1:3),betaUP(1:3),enth_old,v3d(1:3)
 real :: f,df,term,lorentz_LEO2,gamfac,pm_dot_b,sqrt_gamma_inv
 integer :: niter, i
 real, parameter :: tol = 1.e-12
 integer, parameter :: nitermax = 100
 logical :: converged
 ierr = 0

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
 enth = 1 + gamma/(gamma-1.)*P/dens

 niter = 0
 converged = .false.
 sqrt_gamma_inv = alpha*sqrtg_inv ! get determinant of 3 spatial metric
 term = rho*sqrt_gamma_inv
 gamfac = gamma/(gamma-1.)
 pm_dot_b = dot_product(pmom,betaUP)

 do while (.not. converged .and. niter < nitermax)
    enth_old = enth
    lorentz_LEO2 = 1.+pmom2/enth_old**2
    lorentz_LEO = sqrt(lorentz_LEO2)
    dens = term/lorentz_LEO

    if (ien_type == ien_etotal) then
       p = max(rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-pm_dot_b),0.)
    elseif (ieos==4) then
       p = (gamma-1.)*dens*polyk
    else
       p = en*dens**gamma
    endif

    f = 1. + gamfac*P/dens - enth_old

    !This line is unique to the equation of state - implemented for adiabatic at the moment
    if (ien_type == ien_etotal) then
       df= -1.+gamfac*(1.-pmom2*p/(enth_old**3*lorentz_LEO2*dens))
    elseif (ieos==4) then
       df = -1. ! Isothermal, I think...
    else
       df = -1. + (gamma*pmom2*P)/(lorentz_LEO2 * enth_old**3 * dens)
    endif

    enth = enth_old - f/df

    ! Needed in dust case when f/df = NaN casuses enth = NaN
    if (enth-1. < tiny(enth)) enth = 1. + 1.5e-6

    niter = niter + 1

    if (abs(enth-enth_old)/enth < tol) converged = .true.
 enddo

 if (.not.converged) ierr = 1

 lorentz_LEO = sqrt(1.+pmom2/enth**2)
 dens = term/lorentz_LEO

 if (ien_type == ien_etotal) then
    p = max(rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-pm_dot_b),0.)
 elseif (ieos==4) then
    p = (gamma-1.)*dens*polyk
 else
    p = en*dens**gamma
 endif

 v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)

! Raise index from down to up
 do i=1,3
    v(i) = dot_product(gammaijUP(:,i),v3d(:))
 enddo

 call get_u(u,P,dens,gamma)

end subroutine conservative2primitive_con_gamma

end module cons2primsolver
