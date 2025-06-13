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
subroutine primitive2conservative(ipart,x,metrici,v,dens,u,P,rho,pmom,en,ien_type)
 use utils_gr,     only:get_u0,dot_product_gr,get_b0
 use metric_tools, only:unpack_metric
 use io,           only:error
 use eos,          only:gmw
 use part,         only:bxyz,bevol,bxyzd
 real, intent(in)  :: x(1:3),metrici(:,:,:)
 real, intent(in)  :: dens,v(1:3),u,P
 real, intent(out) :: rho,pmom(1:3),en
 integer, intent(in) :: ipart,ien_type

 real, dimension(0:3,0:3) :: gcov
 real :: sqrtg,enth,U0,v4U(0:3)
 real :: b2dens,b4(0:3),gam1,bxyzi(1:3),b0
 integer :: i,ierror

 v4U(0) = 1.
 v4U(1:3) = v(:)

 enth = 1. + u + P/dens

 ! get determinant of metric
 call unpack_metric(metrici,sqrtg=sqrtg,gcov=gcov)
 call get_u0(gcov,v,U0,ierror)
 if (ierror > 0) call error('get_u0 in prim2cons','1/sqrt(-v_mu v^mu) ---> non-negative: v_mu v^mu')
 rho = sqrtg*dens*U0
 bxyzi = bxyz(2:4,ipart)

 call get_b0(gcov,v,bxyzi,b0)
 b4 = (/b0,bxyzi/)

 b2dens = dot_product_gr(b4,b4,gcov)/dens
 bevol(1:3,ipart) = (bxyzi - b0*v)/dens
 bxyz(:,ipart) = b4
 do i=0,3
    bxyzd(i+1,ipart) = dot_product(gcov(i,:),b4)
 enddo
 
 do i=1,3
    pmom(i) = U0*(enth+b2dens)*dot_product(gcov(i,:),v4U(:)) - b0*dot_product(gcov(i,:),b4)/(dens*U0)
 enddo

 if (ien_type == ien_etotal) then
    en = dot_product(pmom,v) + (1.+u+0.5*b2dens)/U0
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
subroutine conservative2primitive(ipart,x,metrici,v,dens,u,P,temp,gamma,rho,pmom,en,ierr,ien_type)
 use utils_gr,     only:get_u0,dot_product_gr,get_b0
 use metric_tools, only:unpack_metric
 use io,           only:fatal,warning
 use part,         only:Bxyz,Bevol,bxyzd
 real, intent(in)    :: x(1:3),metrici(:,:,:)
 real, intent(inout) :: dens,P,u,temp,gamma,v(1:3)
 real, intent(in)    :: rho,pmom(1:3),en
 integer, intent(out) :: ierr
 integer, intent(in)  :: ipart,ien_type
 real :: sqrtg,sqrtg_inv,enth_old,v3d(1:3),gamfac,enth
 real :: bcons(1:3),v4U(0:3),U0,bevoli(1:3),bcons_down(1:3),bcons2
 real(kind=16) :: bconspmom,bconspmom2
 real :: enth2,rhoinvg,gam1,v0d,rhog,reduce_fac,vv,vvold
 real :: f,df,const,vterm(1:3),wterm(1:3)
 integer :: niter,i,ierror
 integer :: nitermax = 200
 real :: tol = 1.e-9
 logical :: converged
 real    :: gcov(0:3,0:3),gcon(0:3,0:3)
 ierr = 0

 ! Get metric components from metric array
 call unpack_metric(metrici,sqrtg=sqrtg,gcov=gcov,gcon=gcon)
 sqrtg_inv = 1./sqrtg

 bevoli = bevol(1:3,ipart)
 bcons = bevoli * rho 
 v4U(0) = 1.
 v4U(1:3) = v

 ! constants
 gam1 = gamma-1.
 gamfac = gamma/(gam1)
 rhog = rho*sqrtg
 rhoinvg = rho*sqrtg_inv

 ! Guess W (using previous values of dens and vel)
 call get_u0(gcov,v,U0,ierror)
 enth = rhog*U0*(1. + gamfac*en*dens**(gamma-1.))

 do i=1,3
    bcons_down(i) = dot_product(gcov(1:3,i),bcons(:))
    v3d(i) = dot_product(gcov(:,i),v4u(:))
 enddo
 bcons2 = dot_product(bcons,bcons_down)
 bconspmom = dot_product(bcons,pmom)
 bconspmom2 = bconspmom**2

 niter = 0
 vvold = -dot_product_gr(v4u,v4u,gcov)
 converged = .false.

 do while (.not. converged .and. niter < nitermax)
    enth_old = enth
    enth2 = enth**2

    v3d = rhog*(enth*pmom + bcons_down*bconspmom)/ (enth*(enth+bcons2))
    v0d = (1. - dot_product(gcon(0,1:3),v3d))/gcon(0,0)
    call get_u0(gcon,v3d,U0,ierror,v0d)
    if (ierror /= 0) then
       nitermax = 5000 ! allow more iteration to converge, give warning
       call warning('conservative2primitive','Close to null like. Will be slow to converge',ipart)
       vv = dot_product_gr((/v0d,v3d/),(/v0d,v3d/),gcon)
       reduce_fac = 0.
       do i=1,3
          reduce_fac = reduce_fac + dot_product(gcon(i,1:3)-gcon(0,i)*gcon(0,1:3)/gcon(0,0),v3d(i)*v3d(:))
       enddo
       reduce_fac = sqrt(1.-(vv+vvold)/reduce_fac)
       v3d = v3d*reduce_fac
       v0d = (1. - dot_product(gcon(0,1:3),v3d))/gcon(0,0)
       call get_u0(gcon,v3d,U0,ierror,v0d)
       ierror = 0
    endif

    do i = 1,3
       v(i) = dot_product(gcon(i,:),(/v0d,v3d/))
    enddo
    dens = rhoinvg/U0

    if (ien_type == ien_entropy) then
       f = rhog*U0*(1. + gamfac*en*dens**gam1) - enth_old

       const = -rhog**2*u0**3*(1. + (gamfac-gamma)*en*dens**gam1)
       vterm = v - gcon(0,1:3)/gcon(0,0)
       wterm = (enth2*pmom + (2.*enth+bcons2)*bconspmom*bcons_down)/(enth2*(enth+bcons2)**2)
       df = const*dot_product(vterm,wterm) - 1.
    else
       call fatal('cons2primsolver','only implemented for entropy')
    endif

    enth = enth_old - f/df ! update enth with temp instead of NR
    if (enth > 1.2*enth_old) then
       enth = 1.2*enth_old
    elseif (enth < 0.8*enth_old) then
       enth = 0.8*enth_old
    endif

    ! Needed in dust case when f/df = NaN casuses enth = NaN
    if (enth/rho-1. < tiny(enth)) enth = (1. + 1.5e-6)*rho

    niter = niter + 1

    if (abs(enth-enth_old)/enth < tol) converged = .true.
 enddo

 if (.not.converged) ierr = 1

 v3d = rhog*(enth*pmom + bcons_down*bconspmom)/ (enth*(enth+bcons2))
 v0d = (1. - dot_product(gcon(0,1:3),v3d))/gcon(0,0)
 call get_u0(gcon,v3d,U0,ierror,v0d)
 if (isnan(u0)) call fatal('conservative2primitive','vv > 0, not converged')

 do i = 1,3
    v(i) = dot_product(gcon(i,:),(/v0d,v3d/))
 enddo
 dens = rhoinvg/U0
 p = en*dens**gamma
 call get_u(u,P,dens,gamma)
 Bxyz(2:4,ipart) = bcons/(sqrtg*U0) + bconspmom*rho*U0*v/enth
 call get_b0(gcov,v,bxyz(2:4,ipart),bxyz(1,ipart))
 do i=0,3
    bxyzd(i+1,ipart) = dot_product(gcov(i,:),Bxyz(:,ipart))
 enddo

end subroutine conservative2primitive

end module cons2primsolver
