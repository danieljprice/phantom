!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cons2primsolver
!
! Internal routines containing the GR conservative to
! primitive variable solver, as described in section 7
! of Liptai & Price (2019)
!
! :References:
!   Liptai & Price (2019), MNRAS 485, 819
!   Tejeda (2012), PhD thesis, IAS Trieste
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, metric_tools, part, physcon, units, utils_gr
!
 use eos, only:ieos,polyk
 use part, only:ien_etotal,ien_entropy,ien_entropy_s
 use dim, only:mhd
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
!+
!----------------------------------------------------------------
subroutine primitive2conservative(x,metrici,v,dens,u,P,rho,pmom,en,bprim,bcon,ien_type)
 use utils_gr,     only:get_u0,get_sqrtg,get_sqrt_gamma,four_mag
 use metric_tools, only:unpack_metric
 use io,           only:error
 use eos,          only:gmw,get_entropy
 use part,         only:mhd
 real, intent(in)  :: x(1:3),metrici(:,:,:)
 real, intent(in)  :: dens,v(1:3),u,P,bprim(1:3)
 real, intent(out) :: rho,pmom(1:3),en,bcon(1:3)

 integer, intent(in) :: ien_type
 real, dimension(0:3,0:3) :: gcov
 real :: sqrtg,enth,gvv,U0,v4U(0:3),alpha,sqrt_gamma
 real :: b2,b4up(0:3),b4down(0:3)
 real :: gam1
 integer :: i, mu, ierror

 v4U(0) = 1.
 v4U(1:3) = v(:)

 enth = 1. + u + P/dens

 ! get determinant of metric
 call unpack_metric(metrici,gcov=gcov,alpha=alpha)
 call get_sqrtg(gcov,sqrtg)

 call get_u0(gcov,v,U0,ierror)
 if (ierror > 0) call error('get_u0 in prim2cons','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')
 rho = sqrtg*dens*U0
   
 if (mhd) then
    sqrt_gamma = sqrtg/alpha
    call four_mag(bprim,U0,v,gcov,b4up,b4down,b2)
    enth = enth + 0.5*b2/dens
    bcon = sqrt_gamma*bprim/rho
 endif

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
!+
!----------------------------------------------------------------
subroutine conservative2primitive(x,metrici,v,dens,u,P,temp,gamma,rho,pmom,en,bprim,Bcon,ierr,ien_type)
 use utils_gr,     only:get_sqrtg,get_sqrt_gamma,four_mag
 use metric_tools, only:unpack_metric
 use eos,          only:ieos,gmw,get_entropy,get_p_from_rho_s,gamma_global=>gamma
 use io,           only:fatal
 use physcon,      only:radconst,Rg
 use units,        only:unit_density,unit_ergg
 use inverse4x4,   only:inv4x4
 use vectorutils,  only:matrixinvert3D
 real, intent(in)    :: x(1:3),metrici(:,:,:)
 real, intent(inout) :: dens,P,u,temp,gamma
 real, intent(out)   :: v(1:3),bprim(1:3)
 real, intent(in)    :: rho,pmom(1:3),en,bcon(1:3)
 integer, intent(out) :: ierr
 integer, intent(in)  :: ien_type
 real, dimension(1:3,1:3) :: gammaijUP,gammaijdown
 real :: sqrtg,sqrtg_inv,lorentz_LEO,pmom2,alpha,betadown(1:3),betaUP(1:3),enth_old,v3d(1:3)
 real :: f,df,term,lorentz_LEO2,gamfac,pm_dot_b,sqrt_gamma_inv,enth,gamma1,sqrt_gamma
 real :: vlocalup(3),vlocaldown(3),b2,b4up(0:3),b4down(0:3),bterm(1:3),pterm,U0
 real(kind=8) :: cgsdens,cgsu
 integer :: niter,i
 real, parameter :: tol = 1.e-11
 integer, parameter :: nitermax = 200
 logical :: converged
 real    :: gcov(0:3,0:3)
 ierr = 0

 ! Get metric components from metric array
 call unpack_metric(metrici,gcov=gcov,gammaijUP=gammaijUP,gammaijdown=gammaijdown,alpha=alpha,betadown=betadown,betaUP=betaUP)

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
 call get_sqrt_gamma(gcov,sqrt_gamma)
 sqrt_gamma_inv = alpha*sqrtg_inv ! get determinant of 3 spatial metric
 term = rho*sqrt_gamma_inv
 gamfac = gamma/(gamma-1.)
 pm_dot_b = dot_product(pmom,betaUP)

 if (mhd) then
    bprim = bcon*sqrt_gamma_inv*rho
    vlocalup = (v + betaup)/alpha
    do i = 1,3
       vlocaldown(i) = dot_product(gammaijdown(:,i),vlocalup(:))
    enddo

    lorentz_LEO = 1./sqrt(1. - dot_product(vlocaldown,vlocalup))
    U0 = lorentz_LEO/alpha
    dens = term/lorentz_LEO

    call four_mag(bprim,U0,v,gcov,b4up,b4down,b2)
    enth = enth + 0.5*b2/dens
 endif

 do while (.not. converged .and. niter < nitermax)
    enth_old = enth
    lorentz_LEO2 = 1.+pmom2/enth_old**2
    lorentz_LEO = sqrt(lorentz_LEO2)
    dens = term/lorentz_LEO

    if (ien_type == ien_etotal) then
       p = rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-pm_dot_b)
    elseif (ieos==4) then
       p = (gamma-1.)*dens*polyk
    elseif (ien_type == ien_entropy_s) then
       call get_p_from_rho_s(ieos,en,dens,gmw,P,temp)
       select case(ieos)
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
          call fatal('cons2primsolver','only implemented for eos 2 and 12')
       end select
    else
       p = en*dens**gamma
    endif

    if (ien_type /= ien_entropy_s) then
       f = 1. + gamfac*P/dens - enth_old
       if (mhd) then
          v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)
          do i=1,3
             v(i) = dot_product(gammaijUP(:,i),v3d(:))
          enddo

          U0 = lorentz_LEO/alpha
          call four_mag(bprim,U0,v,gcov,b4up,b4down,b2,pmom=pmom,enth=enth_old)
          f = f + 0.5*b2/dens
        endif
    endif   

    !This line is unique to the equation of state - implemented for adiabatic at the moment
    if (ien_type == ien_etotal) then
       df= -1.+gamfac*(1.-pmom2*p/(enth_old**3*lorentz_LEO2*dens))
    elseif (ieos==4) then
       df = -1. ! Isothermal, I think...
    elseif (ien_type == ien_entropy_s) then
       continue
    else
       pterm = gamma*P
       if (mhd) pterm = pterm - 0.5*b2
       df = -1. + (pmom2*pterm)/(lorentz_LEO2 * enth_old**3 * dens)
       if (mhd) then
          bterm = (b4up(1:3) + b4up(0)*betaUP(:) - 2*alpha*lorentz_LEO*bprim)/(enth_old*lorentz_LEO2)
          df = df - (dot_product(bterm,b4down(1:3)) + b2/enth_old)/dens
       endif
    endif

    if (ien_type /= ien_entropy_s) then ! .or. ieos /= 12) then
       enth = enth_old - f/df
    else
       enth = 1. + gamfac*P/dens ! update enth with temp instead of NR
    endif

    ! Needed in dust case when f/df = NaN casuses enth = NaN
    if (enth-1. < tiny(enth)) then
       enth = 1. + 1.5e-6
       if (mhd) enth = enth + 0.5*b2/dens
    endif

    niter = niter + 1

    if (abs(enth-enth_old)/enth < tol) converged = .true.
 enddo

 if (.not.converged) ierr = 1
 lorentz_LEO = sqrt(1.+pmom2/enth**2)
 dens = term/lorentz_LEO

 if (ien_type == ien_etotal) then
    p = rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-pm_dot_b)
 elseif (ieos==4) then
    p = (gamma-1.)*dens*polyk
 elseif (ien_type == ien_entropy_s) then
    call get_p_from_rho_s(ieos,en,dens,gmw,P,temp)
    select case(ieos)
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

end module cons2primsolver
