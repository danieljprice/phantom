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
 
subroutine levi_civita(sqrtgamma,vec1,vec2,res)
 real, intent(in) :: sqrtgamma,vec1(1:3),vec2(1:3)
 real, intent(out) :: res(1:3)

 res(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
 res(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
 res(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
 res = res*sqrtgamma 
 
end subroutine

subroutine efield(sqrt_gamma,gammaijup,v,b,eup,edown,vup_in)
   real, intent(in) :: v(1:3),b(1:3),sqrt_gamma,gammaijup(1:3,1:3)
   real, intent(in), optional :: vup_in(1:3)
   real, intent(out) :: eup(1:3),edown(1:3)
   integer :: i 
   real :: vup(3)

   if (present(vup_in)) then
      vup = vup_in
   else
      do i = 1,3
         vup(i) = dot_product(gammaijUP(:,i),v(:))
      enddo
   endif

   call levi_civita(sqrt_gamma,vup,b,edown)
   edown = - edown
   do i = 1,3
      eup(i) = dot_product(gammaijUP(:,i),edown(:))
   enddo
end subroutine

subroutine pmom_comp(sqrt_gamma,eup,b,rho,pmom,pmag,phydro)
   real, intent(in) :: sqrt_gamma,eup(1:3),b(1:3),rho,pmom(1:3)
   real, intent(out) :: pmag(1:3)
   real, intent(out), optional :: phydro(1:3)
 
   call levi_civita(sqrt_gamma,eup,b,pmag)
   pmag = pmag/rho
   if (present(phydro)) phydro = pmom - pmag
   
end subroutine

!----------------------------------------------------------------
!+
!  Construct conserved variables from the primitive variables
!  primitive variables are (v^i,d,u,P); v i
!  conserved variables are (rho,pmom_i,en)
!+
!----------------------------------------------------------------
subroutine primitive2conservative(x,metrici,v,dens,u,P,rho,pmom,en,ien_type,B,Bcon)
 use utils_gr,     only:get_u0,get_sqrtg,get_sqrt_gamma
 use metric_tools, only:unpack_metric
 use io,           only:error
 use eos,          only:gmw,get_entropy
 use part,         only:mhd
 real, intent(in)  :: x(1:3),metrici(:,:,:)
 real, intent(in)  :: dens,v(1:3),u,P
 real, intent(out) :: rho,pmom(1:3),en
 real, intent(in), optional :: B(1:3)
 real, intent(out), optional :: Bcon(1:3)

 integer, intent(in) :: ien_type
 real, dimension(0:3,0:3) :: gcov
 real :: sqrtg, enth, gvv, U0, v4U(0:3)
 real :: alpha,betaup(1:3),gammaijup(1:3,1:3),gammaijdown(1:3,1:3),sqrt_gamma
 real :: b2,eup(3),edown(3),pmag(3),vlocal(3)
 real :: gam1
 integer :: i, mu, ierror

 v4U(0) = 1.
 v4U(1:3) = v(:)

 enth = 1. + u + P/dens

 ! get determinant of metric
 call unpack_metric(metrici,gcov=gcov)
 call get_sqrtg(gcov,sqrtg)
 if (mhd) then
    call unpack_metric(metrici,alpha=alpha,betaUP=betaUP,gammaijup=gammaijup,gammaijdown=gammaijdown)
    vlocal = (v+betaUP)/alpha
    sqrt_gamma = sqrtg/alpha
    call get_sqrt_gamma(gcov,sqrt_gamma)
 endif

 call get_u0(gcov,v,U0,ierror)
 if (ierror > 0) call error('get_u0 in prim2cons','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')
 rho = sqrtg*dens*U0
   
 do i=1,3
    pmom(i) = U0*enth*dot_product(gcov(i,:),v4U(:))
 enddo

 if (mhd) then 
    call efield(sqrt_gamma,gammaijup,vlocal,b,eup,edown,vup_in=vlocal)

    b2 = 0.
    do i=1,3
       b2 = b2 + b(i)*dot_product(gammaijdown(i,:),b(:))
    enddo

    call pmom_comp(sqrt_gamma,eup,b,rho,pmom,pmag)
    pmom = pmom + pmag
    bcon = sqrt_gamma*b
 endif

 gvv = 0.
 do mu=0,3
    do i=1,3
       gvv = gvv + gcov(i,mu)*v4U(mu)*v4U(i)
    enddo
 enddo

 if (ien_type == ien_etotal) then
    en = U0*enth*gvv + (1.+u)/U0
    if (mhd) en = en + 0.5*(dot_product(eup,edown)+b2)/rho
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
subroutine conservative2primitive(x,metrici,v,dens,u,P,temp,gamma,rho,pmom,en,ierr,ien_type,B,Bcon)
 use utils_gr,     only:get_sqrtg,get_sqrt_gamma
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
 real, intent(inout), optional :: B(1:3)
 real, intent(in), optional :: Bcon(1:3)
 integer, intent(out) :: ierr
 integer, intent(in)  :: ien_type
 real, dimension(1:3,1:3) :: gammaijUP,gammaijdown
 real :: sqrtg,sqrtg_inv,lorentz_LEO,pmom2,alpha,betadown(1:3),betaUP(1:3),enth_old,v3d(1:3)
 real :: f,df,term,lorentz_LEO2,gamfac,pm_dot_b,sqrt_gamma_inv,enth,gamma1,sqrt_gamma
 real :: vlocal(3),vlocal_old(3),vlocalup(3),eup(3),edown(3),phydro(3),phydroup(3),pmag(3)
 real :: b2,fv(3),dfv(3,3),invdfv(3,3),lev(3),levup(3)
 real(kind=8) :: cgsdens,cgsu
 integer :: niter,i,j
 real, parameter :: tol = 1.e-11
 integer, parameter :: nitermax = 100
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
    b = bcon*sqrt_gamma_inv
    vlocal = (v + betaup)/alpha

    b2 = 0.
    do i = 1,3
       b2 = b2 + 0.5*b(i)*dot_product(gammaijdown(:,i),b(:))
    enddo

    call efield(sqrt_gamma,gammaijup,vlocal,b,eup,edown)
 endif

 do while (.not. converged .and. niter < nitermax)
    enth_old = enth

    if (mhd) then
       vlocal_old = vlocal 
       do i = 1,3
          vlocalup(i) = dot_product(gammaijup(:,i),vlocal(:))
       enddo

       call efield(sqrt_gamma,gammaijup,vlocal,b,eup,edown)
       call pmom_comp(sqrt_gamma,eup,b,rho,pmom,pmag,phydro)

       do i=1,3
          phydroup(i) = dot_product(gammaijUP(:,i),phydro(:))
       enddo
       pmom2 = dot_product(phydro,phydroup)
       pm_dot_b = dot_product(phydro,betaUP)
    endif
    lorentz_LEO2 = 1.+pmom2/enth_old**2
    lorentz_LEO = sqrt(lorentz_LEO2)
    dens = term/lorentz_LEO

    if (ien_type == ien_etotal) then
       p = rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-pm_dot_b)
       if (mhd) p = p + 0.5*dot_product(edown,eup) + b2
       p = max(p,0.)
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

    if (ien_type /= ien_entropy_s) f = 1. + gamfac*P/dens - enth_old
    if (mhd) fv = phydro/(enth*lorentz_LEO) - vlocal_old

    !This line is unique to the equation of state - implemented for adiabatic at the moment
    if (ien_type == ien_etotal) then
       df= -1.+gamfac*(1.-pmom2*p/(enth_old**3*lorentz_LEO2*dens))
    elseif (ieos==4) then
       df = -1. ! Isothermal, I think...
    elseif (ien_type == ien_entropy_s) then
       continue
    else
       df = -1. + (gamma*pmom2*P)/(lorentz_LEO2 * enth_old**3 * dens)
    endif

    ! also root finding on velocity if mhd
    ! same for energy and entropy K
    if (mhd) then
       do i=1,3
          call levi_civita(sqrt_gamma,gammaijup(i,:),b,lev)
          do j=1,3
             levup(j) = dot_product(gammaijup(j,:),lev(:))
          enddo
          call levi_civita(sqrt_gamma,levup,b,lev)
        
          dfv(:,i) = -phydro*lorentz_LEO*vlocalup/enth_old + sqrt_gamma*lev/rho
          dfv(i,i) = dfv(i,i) - 1.
       enddo

       call matrixinvert3D(dfv,invdfv,ierr)
       if (ierr /= 0) call fatal('cons2primsolver','dfv matrix det=0')
    endif

    ! find next step
    if (ien_type /= ien_entropy_s) then ! .or. ieos /= 12) then
       enth = enth_old - f/df
       if (mhd) then
          do i = 1,3
             df = dot_product(invdfv(i,:),fv)
             vlocal(i) = vlocal_old(i) - df
          enddo
       endif
    else
       enth = 1. + gamfac*P/dens ! update enth with temp instead of NR
    endif

    ! Needed in dust case when f/df = NaN casuses enth = NaN
    if (enth-1. < tiny(enth)) enth = 1. + 1.5e-6

    niter = niter + 1

    if (abs(enth-enth_old)/enth < tol) converged = .true.
    if (mhd .and. converged) then
       do i=1,3
          if (abs(vlocal(i)) > tol) then
             if (abs(vlocal(i)-vlocal_old(i))/abs(vlocal(i)) > tol) converged = .false.
          else
             if (abs(vlocal(i)-vlocal_old(i)) > tol) converged = .false.
          endif
       enddo
    endif
 enddo

 if (.not.converged) ierr = 1

 if (mhd) then
    call efield(sqrt_gamma,gammaijup,vlocal,b,eup,edown)
    call pmom_comp(sqrt_gamma,eup,b,rho,pmom,pmag,phydro)

    pmom2 = 0.
    do i=1,3
       pmom2 = pmom2 + phydro(i)*dot_product(gammaijUP(:,i),phydro(:))
    enddo
 endif
 lorentz_LEO = sqrt(1.+pmom2/enth**2)
 dens = term/lorentz_LEO

 if (ien_type == ien_etotal) then
    p = rho*sqrtg_inv*(enth*lorentz_LEO*alpha-en-pm_dot_b)
    if (mhd) p = p + 0.5*dot_product(edown,eup) + b2 + rho*dot_product(pmag,betaup)*sqrtg_inv
    p = max(p,0.)
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
    v3d(:) = alpha*vlocal-betadown(:)
 else
    v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)
 endif

! Raise index from down to up
 do i=1,3
    v(i) = dot_product(gammaijUP(:,i),v3d(:))
 enddo

 if (ien_type /= ien_entropy_s) call get_u(u,P,dens,gamma)

end subroutine conservative2primitive

end module cons2primsolver
