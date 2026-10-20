!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cons2primsolver
!
! Internal routines containing the GR conservative to
! primitive variable solver, as described in section 7
! of Liptai & Price (2019). For GRMHD the solver follows
! the 1D_W scheme of Noble et al. (2006)
!
! :References:
!   Liptai & Price (2019), MNRAS 485, 819
!   Noble et al. (2006), ApJ 641, 626
!   Tejeda (2012), PhD thesis, IAS Trieste
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, io, metric_tools, part, physcon, units, utils_gr
!
 use eos, only:ieos,polyk
 use part, only:ien_etotal,ien_entropy,ien_entropy_s
 use dim,  only:mhd
 implicit none

 public :: conservative2primitive,primitive2conservative

 private :: get_u

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
!
!  For MHD, bxyz holds the comoving magnetic field (b^1,b^2,b^3,b^0):
!  b^i is input and b^0 is returned, along with the conserved field
!  bevol = B^i/rho* = (b^i - b^0 v^i)/dens
!+
!----------------------------------------------------------------
subroutine primitive2conservative(x,metrici,v,dens,u,P,rho,pmom,en,ien_type,bxyz,bevol)
 use utils_gr,     only:get_u0,get_sqrtg,get_b0,dot_product_gr
 use metric_tools, only:unpack_metric
 use io,           only:error
 use eos,          only:gmw,get_entropy
 real,    intent(in)  :: x(1:3),metrici(:,:,:)
 real,    intent(in)  :: dens,v(1:3),u,P
 real,    intent(out) :: rho,pmom(1:3),en
 integer, intent(in)  :: ien_type
 real,    intent(inout), optional :: bxyz(4)
 real,    intent(out),   optional :: bevol(3)
 real :: gcov(0:3,0:3)
 real :: sqrtg, enth, gvv, U0, v4U(0:3)
 real :: gam1,b4(0:3),b0,b2dens,b0dens
 integer :: i, mu, ierror

 v4U(0) = 1.
 v4U(1:3) = v(:)

 enth = 1. + u + P/dens

 ! get determinant of metric
 call unpack_metric(metrici,gcov=gcov)
 call get_sqrtg(gcov,sqrtg)

 call get_u0(gcov,v,U0,ierror)
 if (ierror > 0) call error('get_u0 in prim2cons','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')
 rho = sqrtg*dens*U0
 do i=1,3
    pmom(i) = U0*enth*dot_product(gcov(i,:),v4U(:))
 enddo

 b2dens = 0.
 if (present(bevol)) bevol = 0.
 if (mhd .and. present(bxyz)) then
    ! p_i = (w + b^2/dens) U_i - b^0 b_i/(dens U^0)
    call get_b0(gcov,v,bxyz(1:3),b0)
    bxyz(4) = b0
    b4(0)   = b0
    b4(1:3) = bxyz(1:3)
    b2dens  = dot_product_gr(b4,b4,gcov)/dens
    b0dens  = b0/(dens*U0)
    do i=1,3
       pmom(i) = pmom(i) + b2dens*U0*dot_product(gcov(i,:),v4U(:)) - b0dens*dot_product(gcov(i,:),b4)
    enddo
    if (present(bevol)) bevol = (bxyz(1:3) - b0*v)/dens
 endif

 gvv = 0.
 do mu=0,3
    do i=1,3
       gvv = gvv + gcov(i,mu)*v4U(mu)*v4U(i)
    enddo
 enddo

 if (ien_type == ien_etotal) then
    en = U0*enth*gvv + (1.+u)/U0
    if (b2dens > 0.) en = dot_product(pmom,v) + (1. + u + 0.5*b2dens)/U0
 elseif (ien_type == ien_entropy_s) then
    en = get_entropy(dens,P,gmw,ieos)
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
!  solve for primitive variables from the conserved variables
!  for equations of state where gamma is constant
!
!  For MHD, pass the conserved field bevol = B^i/rho*; the comoving
!  field (b^1,b^2,b^3,b^0) is then returned in bxyz
!+
!----------------------------------------------------------------
subroutine conservative2primitive(x,metrici,v,dens,u,P,temp,gamma,rho,pmom,en,ierr,ien_type,bevol,bxyz)
 use utils_gr,     only:get_sqrtg,get_sqrt_gamma
 use metric_tools, only:unpack_metric
 use eos,          only:ieos,gmw,get_entropy,get_p_from_rho_s,get_u_from_rho_s,gamma_global=>gamma
 use io,           only:fatal
 use physcon,      only:radconst,Rg
 use units,        only:unit_density,unit_ergg
 real,    intent(in)    :: x(1:3),metrici(:,:,:)
 real,    intent(inout) :: dens,P,u,temp,gamma
 real,    intent(out)   :: v(1:3)
 real,    intent(in)    :: rho,pmom(1:3),en
 integer, intent(out)   :: ierr
 integer, intent(in)    :: ien_type
 real,    intent(in),  optional :: bevol(3)
 real,    intent(out), optional :: bxyz(4)
 real :: gammaijUP(1:3,1:3)
 real :: sqrtg,sqrtg_inv,lorentz_LEO,pmom2,alpha,betadown(1:3),betaUP(1:3),enth_old,v3d(1:3)
 real :: f,df,term,lorentz_LEO2,gamfac,pm_dot_b,sqrt_gamma_inv,enth,gamma1,sqrt_gamma
 real :: dens_eos,p_eos,temp_eos
 real(kind=8) :: cgsdens,cgsu
 integer :: niter, i
 real, parameter :: tol = 1.e-12
 real, parameter :: dens_eps = 1.e-12
 integer, parameter :: nitermax = 100
 logical :: converged,have_eos_cache
 real    :: gcov(0:3,0:3)

 ! GRMHD solver if there is a magnetic field, otherwise the hydro solver below
 if (mhd .and. present(bevol) .and. present(bxyz)) then
    if (any(abs(bevol) > 0.)) then
       call conservative2primitive_mhd(metrici,v,dens,u,P,gamma,rho,pmom,en,bevol,bxyz,ierr,ien_type)
       return
    endif
 endif
 if (present(bxyz)) bxyz = 0.

 ierr = 0

 ! Get metric components from metric array
 call unpack_metric(metrici,gcov=gcov,gammaijUP=gammaijUP,alpha=alpha,betadown=betadown,betaUP=betaUP)

 ! Retrieve sqrt(g)
 call get_sqrtg(gcov,sqrtg)
 sqrtg_inv = 1./sqrtg
 pmom2 = 0.
 do i=1,3
    pmom2 = pmom2 + pmom(i)*dot_product(gammaijUP(:,i),pmom(:))
 enddo

 ! Guess enthalpy (using previous values of dens and pressure)
 enth = 1 + gamma/(gamma-1.)*P/dens

 niter = 0
 converged = .false.
 have_eos_cache = .false.
 call get_sqrt_gamma(gcov,sqrt_gamma)
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
    elseif (ien_type == ien_entropy_s) then
       call get_p_from_rho_s(ieos,en,dens,gmw,P,temp)
       dens_eos = dens
       p_eos    = P
       temp_eos = temp
       have_eos_cache = .true.
       select case(ieos)
       case (10)
          ! inputs and outputs are all in code units
          call get_u_from_rho_s(ieos,en,dens,u)

       case (12)
          cgsdens = dens * unit_density
          cgsu = 1.5*rg*temp/gmw + radconst*temp**4/cgsdens
          u = real(cgsu / unit_ergg)
          if (u > 0.) then
             gamma1 = P/(u*dens)
             gamma = 1. + gamma1
             gamfac = gamma/gamma1
          else
             gamma = gamma_global
             gamfac = gamma/(gamma-1.)
          endif
       case (2)
       case default
          call fatal('cons2primsolver','only implemented for eos 2,10 and 12')
       end select
    else
       p = en*dens**gamma
    endif

    if (ien_type /= ien_entropy_s) f = 1. + gamfac*P/dens - enth_old

    !This line is unique to the equation of state - implemented for adiabatic at the moment
    if (ien_type == ien_etotal) then
       df= -1.+gamfac*(1.-pmom2*p/(enth_old**3*lorentz_LEO2*dens))
    elseif (ieos==4) then
       df = -1. ! Isothermal, I think...
    elseif (ien_type == ien_entropy_s) then
    else
       df = -1. + (gamma*pmom2*P)/(lorentz_LEO2 * enth_old**3 * dens)
    endif

    if (ien_type /= ien_entropy_s) then ! .or. ieos /= 12) then
       enth = enth_old - f/df
    else
       enth = 1. + gamfac*P/dens ! update enth with temp instead of NR
    endif

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
 elseif (ien_type == ien_entropy_s) then
    if (have_eos_cache .and. abs(dens-dens_eos) <= dens_eps*dens) then
       P    = p_eos
       temp = temp_eos
    else
       call get_p_from_rho_s(ieos,en,dens,gmw,P,temp)
    endif
    select case(ieos)
    case (10)
       ! inputs and outputs are all in code units
       call get_u_from_rho_s(ieos,en,dens,u)

    case (12)
       cgsdens = dens * unit_density
       cgsu = 1.5*rg*temp/gmw + radconst*temp**4/cgsdens
       u = real(cgsu / unit_ergg)
       if (u > 0.) then
          gamma = 1. + P/(u*dens)
       else
          gamma = gamma_global
       endif
    case (2)
       call get_u(u,P,dens,gamma)
    end select
 else
    p = en*dens**gamma
 endif

 v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)

! Raise index from down to up
 do i=1,3
    v(i) = dot_product(gammaijUP(:,i),v3d(:))
 enddo

 if (ien_type /= ien_entropy_s) call get_u(u,P,dens,gamma)

end subroutine conservative2primitive

!----------------------------------------------------------------
!+
!  solve for primitive variables from the conserved variables
!  in GRMHD, using the 1D_W scheme of Noble et al. (2006).
!  We solve f(W) = w(W) U^0(W) rho* sqrt(-g) - W = 0 for
!  W = w rho* sqrt(-g) U^0 by Newton-Raphson, where
!  v_i(W) = rho* sqrt(-g) (W p_i + (p_j B^j) B_i)/(W (W + B^2))
!  Only implemented for the entropy variable with an adiabatic eos
!+
!----------------------------------------------------------------
subroutine conservative2primitive_mhd(metrici,v,dens,u,P,gamma,rho,pmom,en,bevol,bxyz,ierr,ien_type)
 use utils_gr,     only:get_sqrtg
 use metric_tools, only:unpack_metric
 use io,           only:fatal
 real,    intent(in)    :: metrici(:,:,:)
 real,    intent(inout) :: dens,P,u,gamma
 real,    intent(out)   :: v(1:3),bxyz(4)
 real,    intent(in)    :: rho,pmom(1:3),en,bevol(3)
 integer, intent(out)   :: ierr
 integer, intent(in)    :: ien_type
 real, parameter :: tol = 1.e-12
 integer, parameter :: nitermax = 100
 real    :: gcov(0:3,0:3),gcon(0:3,0:3)
 real    :: sqrtg,rhog,gam1,gamfac,bcons(3),bdown(3),bcons2,pdotb
 real    :: W,W_old,U0,vdown(3),Kdens,f,df,dvdW(3),vterm(3),enth
 integer :: niter,i,ierror
 logical :: converged

 if (ien_type /= ien_entropy) call fatal('cons2primsolver','GRMHD only implemented for entropy variable')
 ierr = 0
 call unpack_metric(metrici,gcov=gcov,gcon=gcon)
 call get_sqrtg(gcov,sqrtg)
 rhog   = rho*sqrtg
 gam1   = gamma - 1.
 gamfac = gamma/gam1

 ! conserved field B^i, its lowered form B_i, B^2 and p_j B^j
 bcons = bevol*rho
 do i=1,3
    bdown(i) = dot_product(gcov(i,1:3),bcons)
 enddo
 bcons2 = dot_product(bcons,bdown)
 pdotb  = dot_product(pmom,bcons)

 ! initial guess for W from the previous density, U^0 = rho*/(sqrt(-g) dens)
 U0 = rho/(sqrtg*dens)
 W  = rhog*U0*(1. + gamfac*en*dens**gam1)

 niter = 0
 converged = .false.
 do while (.not. converged .and. niter < nitermax)
    W_old = W
    niter = niter + 1
    call get_v_from_W(W,rhog,pmom,bdown,pdotb,bcons2,gcon,vdown,v,U0,ierror)
    if (ierror /= 0) then
       ! W too small gives superluminal v; increase W and try again
       W = 2.*W_old
       cycle
    endif
    dens  = rho/(sqrtg*U0)
    Kdens = en*dens**gam1   ! P/dens

    ! f(W) and f'(W) = rhog dU^0/dW d(U^0 w)/dU^0 - 1
    f = rhog*U0*(1. + gamfac*Kdens) - W
    dvdW  = -rhog*(W**2*pmom + (2.*W + bcons2)*pdotb*bdown)/(W**2*(W + bcons2)**2)
    vterm = v - gcon(0,1:3)/gcon(0,0)
    df = rhog*U0**3*(1. + (gamfac - gamma)*Kdens)*dot_product(vterm,dvdW) - 1.

    W = W_old - f/df
    ! limit the step size, and require w >= 1
    W = min(max(W,0.5*W_old),2.*W_old)
    W = max(W,rhog*U0)
    if (abs(W - W_old) < tol*W) converged = .true.
 enddo
 if (.not.converged) ierr = 1

 ! primitive variables from the converged W
 call get_v_from_W(W,rhog,pmom,bdown,pdotb,bcons2,gcon,vdown,v,U0,ierror)
 if (ierror /= 0) ierr = 2
 dens = rho/(sqrtg*U0)
 P    = en*dens**gamma
 call get_u(u,P,dens,gamma)
 enth = 1. + u + P/dens

 ! comoving field: b^0 = p_j B^j/(sqrt(-g) w) and b^i = B^i/(sqrt(-g) U^0) + b^0 v^i
 bxyz(4)   = pdotb/(sqrtg*enth)
 bxyz(1:3) = bcons/(sqrtg*U0) + bxyz(4)*v

end subroutine conservative2primitive_mhd

!----------------------------------------------------------------
!+
!  coordinate velocity and U^0 as a function of W in the GRMHD
!  primitive solver: v_i from Noble et al. (2006), v_0 from
!  v^0 = g^{0 mu} v_mu = 1, and U^0 = (-g^{mu nu} v_mu v_nu)^(-1/2)
!+
!----------------------------------------------------------------
pure subroutine get_v_from_W(W,rhog,pmom,bdown,pdotb,bcons2,gcon,vdown,v,U0,ierr)
 real,    intent(in)  :: W,rhog,pmom(3),bdown(3),pdotb,bcons2,gcon(0:3,0:3)
 real,    intent(out) :: vdown(3),v(3),U0
 integer, intent(out) :: ierr
 real    :: v0d,vv
 integer :: i

 ierr  = 0
 vdown = rhog*(W*pmom + pdotb*bdown)/(W*(W + bcons2))
 v0d   = (1. - dot_product(gcon(0,1:3),vdown))/gcon(0,0)
 do i=1,3
    v(i) = gcon(i,0)*v0d + dot_product(gcon(i,1:3),vdown)
 enddo
 ! v^mu v_mu with v^0 = 1
 vv = v0d + dot_product(v,vdown)
 if (vv < 0.) then
    U0 = 1./sqrt(-vv)
 else
    U0 = huge(U0)
    ierr = 1
 endif

end subroutine get_v_from_W

end module cons2primsolver
