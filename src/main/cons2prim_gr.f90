module cons2prim_gr
 use eos, only: gamma
implicit none

public :: conservative2primitive,primitive2conservative,get_u,get_pressure
integer, parameter :: ierr_notconverged = 1

private :: get_enthalpy

integer, parameter, private :: eos_type = 2 ! cons2prim has been written for only adiabatic eos.

private

contains

!=========================
subroutine get_pressure(P,dens,u)
use eos, only: equationofstate
   real, intent(out) :: P
   real, intent(in)  :: dens, u
   real :: ponrho, spsound

   call equationofstate(eos_type,ponrho,spsound,dens,0.,0.,0.,u)
   P = ponrho*dens
   ! P = (gamma-1.0)*dens*u

end subroutine

subroutine get_u(u,P,dens)
use eos, only: equationofstate
   real, intent(in)  :: dens,P
   real, intent(out) :: u
   real :: ponrho, spsound

   call equationofstate(eos_type,ponrho,spsound,dens,0.,0.,0.)
   u = ponrho/(gamma-1.)

end subroutine

subroutine get_enthalpy(enth,dens,P)
use eos, only: equationofstate
   real, intent(in)  :: dens,P
   real, intent(out) :: enth
   real :: ponrho, spsound

   call equationofstate(eos_type,ponrho,spsound,dens,0.,0.,0.)
   enth = 1.+ponrho*(gamma/(gamma-1.))

   ! Needed in dust case when dens = 0 causes P/dens = NaN and therefore enth = NaN
   if(abs(p)<tiny(p)) enth=1.

end subroutine
!=========================

!----------------------------------------------------------------
!+
!  Construct conserved variables from the primitive variables
!  primitive variables are (v^i,d,u,P); v i
!  conserved variables are (rho,pmom_i,en)
!+
!----------------------------------------------------------------
subroutine primitive2conservative(x,v,dens,u,P,rho,pmom,en,en_type)
   use utils_gr,     only: get_u0
   use metric_tools, only: get_metric
   real, intent(in)  :: x(1:3)
   real, intent(in)  :: dens,v(1:3),u,P
   real, intent(out) :: rho,pmom(1:3),en
   character(len=*), intent(in) :: en_type
   real, dimension(0:3,0:3) :: gcov, gcon
   real :: sqrtg, enth, gvv, U0, v4U(0:3)
   integer :: i, mu

   v4U(0) = 1.
   v4U(1:3) = v(:)

   call get_enthalpy(enth,dens,p) !enth = 1.+ u + P/dens

   call get_metric(x,gcov,gcon,sqrtg)
   call get_u0(x,v,U0)
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
   en = U0*enth*gvv + (1.+u)/U0

   if (en_type == "entropy") en = P/(dens**gamma)

end subroutine primitive2conservative

subroutine conservative2primitive(x,v,dens,u,P,rho,pmom,en,ierr,en_type)
   use utils_gr, only: dot_product_gr, get_metric3plus1
   use metric_tools, only: get_metric
   real, intent(in)    :: x(1:3)
   real, intent(inout) :: dens,P
   real, intent(out)   :: v(1:3),u
   real, intent(in)    :: rho,pmom(1:3),en
   integer, intent(out) :: ierr
   character(len=*), intent(in) :: en_type
   real, dimension(0:3,0:3) :: gcov,gcon
   real, dimension(1:3,1:3) :: gammaijdown, gammaijUP
   real :: sqrtg,enth,lorentz_LEO,pmom2,alpha,beta(1:3),enth_old,v3d(1:3)
   real :: f,df
   integer :: niter, i,j
   real, parameter :: tol = 1.e-12
   integer, parameter :: nitermax = 100
   logical :: converged
   ierr = 0

   call get_metric3plus1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
   pmom2 = dot_product_gr(pmom,pmom,gammaijUP)

   ! Guess enthalpy (using previous values of dens and pressure)
   call get_enthalpy(enth,dens,p)

   niter = 0
   converged = .false.
   do while (.not. converged .and. niter < nitermax)
      enth_old = enth
      lorentz_LEO = sqrt(1.+pmom2/enth_old**2)
      dens = rho*alpha/(sqrtg*lorentz_LEO)

      p = max(rho/sqrtg*(enth*lorentz_LEO*alpha-en-dot_product_gr(pmom,beta,gammaijUP)),0.)
      if (en_type == 'entropy') p = en*(rho*alpha/(sqrtg*lorentz_LEO))**gamma

      call get_enthalpy(enth,dens,p)

      f = enth-enth_old

      !This line is unique to the equation of state - implemented for adiabatic at the moment
      df= -1.+(gamma/(gamma-1.))*(1.-pmom2*p/(enth_old**3*lorentz_LEO**2*dens))
      if (en_type == 'entropy') df = -1. + (gamma*pmom2*P)/(lorentz_LEO**2 * enth_old**3 * dens)

      enth = enth_old - f/df

      ! Needed in dust case when f/df = NaN casuses enth = NaN
      if (abs(enth_old-1.)<tiny(enth_old)) enth=1.

      niter = niter + 1

      if (abs(enth-enth_old)/enth < tol) converged = .true.
   enddo

   if (.not.converged) then
      print*,' Warning: cons2prim did not converge, reached max number of iterations',niter
      print*,' enthold  =',enth_old
      print*,' enthnew  =',enth
      print*,' error    =',abs(enth-enth_old)/enth
      ierr = ierr_notconverged
      return
   endif

   lorentz_LEO = sqrt(1.+pmom2/enth**2)
   dens = rho*alpha/(sqrtg*lorentz_LEO)

   p = max(rho/sqrtg*(enth*lorentz_LEO*alpha-en-dot_product_gr(pmom,beta,gammaijUP)),0.)
   if (en_type == 'entropy') p = en*(rho*alpha/(sqrtg*lorentz_LEO))**gamma

   v3d(:) = alpha*pmom(:)/(enth*lorentz_LEO)-beta(:)

   v = 0.
   do i=1,3
      do j=1,3
         v(i) = v(i) + gammaijUP(i,j)*v3d(j) ! Raise index from down to up
      enddo
   enddo

   call get_u(u,P,dens)

end subroutine conservative2primitive

subroutine get_v_from_p(pmom,v,x)
   ! Conservative to primitive solver for the dust case (Pressure=0)
   use metric_tools, only: get_metric
   use utils_gr, only: dot_product_gr
   real, intent(in) :: pmom(1:3), x(1:3)
   real, intent(out) :: v(1:3)
   real :: en, rho, P, u, dens
   integer :: ierr

   en  = 1.
   rho = 0.
   P   = 0.
   u   = 0.
   dens= 0.
   call conservative2primitive(x,v,dens,u,P,rho,pmom,en,ierr,'energy')


end subroutine get_v_from_p

subroutine get_p_from_v(pmom,v,x)
   use metric_tools, only: get_metric
   real, intent(in) :: v(1:3), x(1:3)
   real, intent(out) :: pmom(1:3)
   real :: rho, en, dens, u, P

   dens = 0.
   u    = 0.
   P    = 0.
   call primitive2conservative(x,v,dens,u,P,rho,pmom,en,'energy')


end subroutine get_p_from_v

end module cons2prim_gr
