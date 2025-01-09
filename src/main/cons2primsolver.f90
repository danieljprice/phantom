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
subroutine primitive2conservative(x,metrici,v,dens,u,P,rho,pmom,en, &
                                  b0,bxyz,bevol,ien_type)
 use utils_gr,     only:get_u0,get_sqrtg,get_sqrt_gamma,dot_product_gr
 use metric_tools, only:unpack_metric
 use io,           only:error
 use eos,          only:gmw,get_entropy
 use part,         only:mhd
 real, intent(in)  :: x(1:3),metrici(:,:,:)
 real, intent(in)  :: dens,v(1:3),u,P
 real, intent(out) :: rho,pmom(1:3),en
 real, intent(in)  :: bxyz(1:3)
 real, intent(out) :: b0,bevol(1:3)
 integer, intent(in) :: ien_type

 real, dimension(0:3,0:3) :: gcov
 real :: sqrtg,enth,gvv,U0,v4U(0:3),alpha,gammaijdown(1:3,1:3),gammaijup(1:3,1:3)
 real :: b2,b4(0:3)
 real :: gam1
 integer :: i, mu, ierror


 v4U(0) = 1.
 v4U(1:3) = v(:)
 !print*, '>>>>>>>>>>>>>>>>>>>'
 !print*, 'dens', dens
    !print*, 'u', u
!print*, 'p', p
!print*, 'v', v4u
!print*, 'b', bxyz
 !print*, '>>>>>>>>>>>>>>>>>>>'

 enth = 1. + u + P/dens

 ! get determinant of metric
 call unpack_metric(metrici,gcov=gcov,gammaijdown=gammaijdown,gammaijup=gammaijup)
 call get_sqrtg(gcov,sqrtg)
 call get_u0(gcov,v,U0,ierror)
 if (ierror > 0) call error('get_u0 in prim2cons','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')
 rho = sqrtg*dens*U0
   
!print*, "U", v4U*U0
!print*, 'UmuUmu', dot_product_gr(v4U*U0,v4U*U0,gcov)
do i=0,3
   !print*, gcov(:,i)
enddo
 if (mhd) then
    b0 = dot_product_gr(bxyz,v,gammaijdown)/(dot_product_gr(v,v,gammaijdown) + 1./U0**2)
   !print*, 'b0',b0, -dot_product_gr(bxyz,v,gammaijdown)/dot_product(v4U,gcov(:,0))
    b4 = (/b0, bxyz/)
    b2 = dot_product_gr(b4,b4,gcov)
    !call four_mag(bxyz,U0,v,gcov,b4,b4d,b2)
    enth = enth + b2/dens
    !bevol = sqrt_gamma*b/rho
    bevol = (bxyz - b0*v)/dens
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
    if (mhd) en = en + 0.5*b2/(dens*U0)
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
 !print*, '<<<<<<<<<<<<<<<<<<<'
 !print*, 'U0', u0
 !print*, 'gamma', gam1
 !print*, 'rho', rho
    !print*, 'pmom', pmom
!print*, 'en', en
!print*, 'enth', enth
!print*, 'Gamma', sqrt(1.+dot_product_gr(pmom,pmom,gammaijup)/enth**2)
!print*, 'b0', b0
!print*, 'bevol', bevol
!print*, 'bmuvmu', dot_product_gr(b4,v4u,gcov)
 !print*, '<<<<<<<<<<<<<<<<<<<'
end subroutine primitive2conservative

!----------------------------------------------------------------
!+
!  solve for primitive variables from the conserved variables
!  for equations of state where gamma is constant
!+
!----------------------------------------------------------------
subroutine conservative2primitive(x,metrici,v,dens,u,P,temp,gamma,rho,pmom,en, &
                                  Bxyz,Bevol,ierr,ien_type)
 use utils_gr,     only:get_sqrtg,get_sqrt_gamma,dot_product_gr
 use metric_tools, only:unpack_metric
 use eos,          only:ieos,gmw,get_entropy,get_p_from_rho_s,gamma_global=>gamma
 use io,           only:fatal
 use physcon,      only:radconst,Rg
 use units,        only:unit_density,unit_ergg
 use inverse4x4,   only:inv4x4
 use vectorutils,  only:matrixinvert3D
 real, intent(in)    :: x(1:3),metrici(:,:,:)
 real, intent(inout) :: dens,P,u,temp,gamma
 real, intent(out)   :: v(1:3)
 real, intent(in)    :: rho,pmom(1:3),en
 real, intent(in)    :: bevol(1:3)
 real, intent(inout) :: bxyz(0:3)
 integer, intent(out) :: ierr
 integer, intent(in)  :: ien_type
 real, dimension(1:3,1:3) :: gammaijUP,gammaijdown
 real :: sqrtg,sqrtg_inv,lorentz_LEO,pmom2,alpha,betadown(1:3),betaUP(1:3),enth_old,v3d(1:3)
 real :: f,df,term,lorentz_LEO2,gamfac,pm_dot_b,sqrt_gamma_inv,enth,gamma1,sqrt_gamma
 real :: vlocalup(3),vlocaldown(3),bterm(1:3),b2,bcons(0:3),v4U(0:3),pterm,U0
 real(kind=8) :: cgsdens,cgsu
 integer :: niter,i
!#ifdef MHD
 !integer, parameter :: nitermax = 300
 !real, parameter :: tol = 1.e-11
!#else
 integer :: nitermax = 100
 real :: tol = 1.e-12
!#endif
 logical :: converged
 real    :: gcov(0:3,0:3)
 ierr = 0

 if (mhd) then
    nitermax = 300
    tol = 5.e-9
 endif

 !print*, '>>>>>>>>>>>>>>>>>>>'
 !print*, 'dens', dens
    !print*, 'u', u
!print*, 'p', p
!print*, 'v', v4u
!print*, 'b', bxyz
 !print*, '>>>>>>>>>>>>>>>>>>>'
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
 if (mhd) then
        !print*, 'bb', bxyz
    enth = enth + dot_product_gr(bxyz,bxyz,gcov)/dens
    bcons(0) = 0.
    bcons(1:3) = bevol(1:3) * rho
    v4U(0) = 1.
    v4U(1:3) = v
 endif

 niter = 0
 converged = .false.
 !call get_sqrt_gamma(gcov,sqrt_gamma)
 sqrt_gamma_inv = alpha*sqrtg_inv ! get determinant of 3 spatial metric
 term = rho*sqrt_gamma_inv
 gamfac = gamma/(gamma-1.)
 pm_dot_b = dot_product(pmom,betaUP)

 !if (mhd) then
 !   bxyz(0) = dot_product_gr(bevol*rho,
 !   bxyz(1:3) = bevol*sqrt_gamma_inv*rho
 !   vlocalup = (v + betaup)/alpha
 !   do i = 1,3
 !      vlocaldown(i) = dot_product(gammaijdown(:,i),vlocalup(:))
 !   enddo

 !   lorentz_LEO = 1./sqrt(1. - dot_product(vlocaldown,vlocalup))
 !   U0 = lorentz_LEO/alpha
 !   dens = term/lorentz_LEO

 !   call four_mag(b,U0,v,gcov,b4,b4d,b2)
 !   enth = enth + 0.5*b2/dens
 !endif
!print*, bcons

 do while (.not. converged .and. niter < nitermax)
    enth_old = enth
    lorentz_LEO2 = 1.+pmom2/enth_old**2
    lorentz_LEO = sqrt(lorentz_LEO2)
!print*, betaup
!print*, 'Gamma',1./sqrt(1.-dot_product_gr((v + betaup)/alpha,(v + betaup)/alpha,gammaijdown))
    dens = term/lorentz_LEO

    if (mhd) then
       v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)

       ! Raise index from down to up
       do i=1,3
          v4u(i) = dot_product(gammaijUP(:,i),v3d(:))
        enddo
    endif

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
          !v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)
          !do i=1,3
          !   v(i) = dot_product(gammaijUP(:,i),v3d(:))
          !enddo

          !U0 = lorentz_LEO/alpha
          !call four_mag(b,U0,v,gcov,b4,b4d,b2,pmom=pmom,enth=enth_old)
        !print*, lorentz_LEO, alpha, sqrt_gamma
        
         bxyz = (bcons/lorentz_LEO + lorentz_LEO/alpha**2 * dot_product(bcons(1:3),v3d) * v4u)*sqrt_gamma_inv
!print*, niter, bxyz
          f = f + dot_product_gr(bxyz,bxyz,gcov)/dens
!print*, 'f', gamfac*P/dens, dot_product_gr(bxyz,bxyz,gcov)/dens
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
       if (mhd) pterm = pterm + dot_product_gr(bxyz,bxyz,gcov)
        !print*, pterm, dot_product_gr(bxyz,bxyz,gcov)
       df = -1. + (pmom2*pterm)/(lorentz_LEO2 * enth_old**3 * dens)
       if (mhd) then
          !bterm = (b4(1:3) + b4(0)*betaUP(:) - 2*alpha*lorentz_LEO*b)/(enth_old*lorentz_LEO2)
          !df = df - (dot_product(bterm,b4d(1:3)) + b2/enth_old)/dens
          df = df + 2.*alpha*bxyz(0)*dot_product_gr(bxyz(1:3),pmom,gammaijup)* &
                    (pmom2/(enth_old**2*lorentz_LEO2) -1.)/(dens*enth_old**2*lorentz_LEO)
       endif
    endif
 !print*, '-------------------'
!print*, 'niter',niter
 !print*, 'Gamma', lorentz_LEO
 !print*, 'enth', enth
 !print*, 'dens', dens
    !print*, 'u', u
!print*, 'p', p
!print*, 'v', v4u
!print*, 'bxyz', bxyz
!print*, 'bcons', bcons
!print*, 'f', f
!print*, 'df', df
 !print*, '-------------------'
!read*

    if (ien_type /= ien_entropy_s) then ! .or. ieos /= 12) then
       enth = enth_old - f/df
    else
       enth = 1. + gamfac*P/dens ! update enth with temp instead of NR
    endif

    ! Needed in dust case when f/df = NaN casuses enth = NaN
    if (enth-1. < tiny(enth)) then
       enth = 1. + 1.5e-6
       if (mhd) enth = enth + dot_product_gr(bxyz,bxyz,gcov)/dens
    endif

    niter = niter + 1

    if (abs(enth-enth_old)/enth < tol) converged = .true.
 enddo

 if (.not.converged) ierr = 1
 lorentz_LEO = sqrt(1.+pmom2/enth**2)
 dens = term/lorentz_LEO

 v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)

! Raise index from down to up
 do i=1,3
    v(i) = dot_product(gammaijUP(:,i),v3d(:))
 enddo

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

 if (mhd) then
    v4u(1:3) = v
         !bxyz = (bcons*alpha/lorentz_LEO + lorentz_LEO/alpha**2 * dot_product(bcons(1:3),v3d) * v4u)/sqrt_gamma
         bxyz = (bcons/lorentz_LEO + lorentz_LEO/alpha**2 * dot_product(bcons(1:3),v3d) * v4u)*sqrt_gamma_inv
    !bxyz = (bcons/lorentz_LEO + lorentz_LEO/alpha**2 * dot_product_gr(bcons,v4u,gcov) * v4u)/sqrt_gamma
 endif

 if (ien_type /= ien_entropy_s) call get_u(u,P,dens,gamma)

 !print*, '-------------------'
!print*, 'niter',niter
 !print*, 'Gamma', lorentz_LEO
 !print*, 'dens', dens
    !print*, 'u', u
!print*, 'p', p
!print*, 'v', v4u
!print*, 'b', bxyz
!print*
!print*, 'f', f
!print*, 'df', df
 !print*, '-------------------'
!read*
end subroutine conservative2primitive

end module cons2primsolver
