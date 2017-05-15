module force_gr
   implicit none
contains

!----------------------------------------------------------------
!+
!  Compute the source terms required on the right hand side of
!  the relativistic momentum equation. These are of the form:
!   T^\mu\nu dg_\mu\nu/dx^i
!+
!----------------------------------------------------------------
subroutine get_forcegr(x,v,dens,u,p,fterm)
 use metric_tools, only: get_metric, get_metric_derivs
 use utils_gr, only: get_u0
 use eos, only: get_enthalpy
 real,    intent(in)  :: x(3),v(3),dens,u,p
 real,    intent(out) :: fterm(3)
 real    :: gcov(0:3,0:3), gcon(0:3,0:3)
 real    :: sqrtg
 real    :: dgcovdx1(0:3,0:3), dgcovdx2(0:3,0:3), dgcovdx3(0:3,0:3)
 real    :: v4(0:3), term(0:3,0:3)
 real    :: enth, uzero
 integer :: i,j
 ! Note to self: try with potential from Tejeda, Rosswog  2013
 call get_metric(x,gcov,gcon,sqrtg)
 call get_metric_derivs(x,dgcovdx1, dgcovdx2, dgcovdx3)

 call get_enthalpy(enth,dens,p)

 ! lower-case 4-velocity
 v4(0) = 1.
 v4(1:3) = v(:)

 ! first component of the upper-case 4-velocity
 call get_u0(x,v,uzero)

 ! energy-momentum tensor times sqrtg on 2rho*
 do j=0,3
  do i=0,3
   term(i,j) = 0.5*(enth*uzero*v4(i)*v4(j) + P*gcon(i,j)/(dens*uzero))
  enddo
 enddo

 ! source term
 fterm(1) = 0.
 fterm(2) = 0.
 fterm(3) = 0.
 do j=0,3
  do i=0,3
   fterm(1) = fterm(1) + term(i,j)*dgcovdx1(i,j)
   fterm(2) = fterm(2) + term(i,j)*dgcovdx2(i,j)
   fterm(3) = fterm(3) + term(i,j)*dgcovdx3(i,j)
  enddo
 enddo

end subroutine get_forcegr

! Wrapper routine to call get_forcegr for a test particle
subroutine get_sourceterms(x,v,fterm)
 real, intent(in)  :: x(3),v(3)
 real, intent(out) :: fterm(3)
 real :: dens,u,p

 P = 0.
 u = 0.
 dens = 1. ! this value does not matter (will cancel in the momentum equation)
 call get_forcegr(x,v,dens,u,p,fterm)
end subroutine get_sourceterms

end module force_gr
