module cons2primsolver
 use eos, only: gamma,ieos,polyk
 implicit none

 public :: conservative2primitive,primitive2conservative

 private :: get_u,get_enthalpy

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

pure subroutine get_u(u,P,dens)
 real, intent(in)  :: dens,P
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

end subroutine

pure subroutine get_enthalpy(enth,dens,P)
 real, intent(in)  :: dens,P
 real, intent(out) :: enth

 ! Needed in dust case when dens = 0 causes P/dens = NaN and therefore enth = NaN
 ! or gamma=1 gives divide-by-zero
 if(P==0. .or. abs(p)<tiny(p)) then
    enth = 1.
 else
    enth = 1.+p/dens*(gamma/(gamma-1.))
 endif

end subroutine
!=========================

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
 integer :: i, mu, ierror

 v4U(0) = 1.
 v4U(1:3) = v(:)

 call get_enthalpy(enth,dens,p) !enth = 1.+ u + P/dens

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
    en = P/(dens**gamma)
 else
    en = U0*enth*gvv + (1.+u)/U0
 endif

end subroutine primitive2conservative

subroutine conservative2primitive(x,metrici,v,dens,u,P,rho,pmom,en,ierr,ien_type)
 use metric_tools, only: unpack_metric
 use io,           only: warning
 real, intent(in)    :: x(1:3),metrici(:,:,:)
 real, intent(inout) :: dens,P
 real, intent(out)   :: v(1:3),u
 real, intent(in)    :: rho,pmom(1:3),en
 integer, intent(out) :: ierr
 integer, intent(in)  :: ien_type
 real, dimension(1:3,1:3) :: gammaijUP
 real :: sqrtg,enth,lorentz_LEO,pmom2,alpha,betadown(1:3),betaUP(1:3),enth_old,v3d(1:3)
 real :: f,df
 integer :: niter, i
 real, parameter :: tol = 1.e-12
 integer, parameter :: nitermax = 100
 logical :: converged
 ierr = 0

 ! Hard coding sqrgt=1 since phantom is always in cartesian coordinates
 sqrtg = 1.

 ! Get metric components from metric array
 call unpack_metric(metrici,gammaijUP=gammaijUP,alpha=alpha,betadown=betadown,betaUP=betaUP)

 pmom2 = 0.
 do i=1,3
    pmom2 = pmom2 + pmom(i)*dot_product(gammaijUP(:,i),pmom(:))
 enddo

 ! Guess enthalpy (using previous values of dens and pressure)
 call get_enthalpy(enth,dens,p)

 niter = 0
 converged = .false.
 do while (.not. converged .and. niter < nitermax)
    enth_old = enth
    lorentz_LEO = sqrt(1.+pmom2/enth_old**2)
    dens = rho*alpha/(sqrtg*lorentz_LEO)

    if (ien_type == ien_entropy) then
       p = en*dens**gamma
    else if (ieos==4) then
       p = (gamma-1.)*dens*polyk
    else
       p = max(rho/sqrtg*(enth*lorentz_LEO*alpha-en-dot_product(pmom,betaUP)),0.)
    endif

    call get_enthalpy(enth,dens,p)

    f = enth-enth_old

    !This line is unique to the equation of state - implemented for adiabatic at the moment
    if (ien_type == ien_entropy) then
       df = -1. + (gamma*pmom2*P)/(lorentz_LEO**2 * enth_old**3 * dens)
    else if (ieos==4) then
       df = -1. ! Isothermal, I think...
    else
       df= -1.+(gamma/(gamma-1.))*(1.-pmom2*p/(enth_old**3*lorentz_LEO**2*dens))
    endif

    enth = enth_old - f/df

    ! Needed in dust case when f/df = NaN casuses enth = NaN
    if (abs(enth_old-1.)<tiny(enth_old)) enth=1.

    niter = niter + 1

    if (abs(enth-enth_old)/enth < tol) converged = .true.
 enddo

 if (.not.converged) then
    call warning('cons2primsolver','enthalpy did not converge. delta enth / enth = ',val=abs(enth-enth_old)/enth)
    ierr = 1
    return
 endif

 lorentz_LEO = sqrt(1.+pmom2/enth**2)
 dens = rho*alpha/(sqrtg*lorentz_LEO)

 if (ien_type == ien_entropy) then
    p = en*dens**gamma
 else
    p = max(rho/sqrtg*(enth*lorentz_LEO*alpha-en-dot_product(pmom,betaUP)),0.)
 endif

 v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-betadown(:)

! Raise index from down to up
 do i=1,3
    v(i) = dot_product(gammaijUP(:,i),v3d(:))
 enddo

 call get_u(u,P,dens)

end subroutine conservative2primitive

end module cons2primsolver
