!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: directsum
!
!  DESCRIPTION:
!  This module computes self-gravity by direct summation
!  over the particles - used to test the treecode.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, kernel, part
!+
!--------------------------------------------------------------------------
module directsum
 implicit none
 public :: directsum_grav

 private

contains

!----------------------------------------------------------------------------
! Calculates the solution to the gravitational Poisson equation
!
! \nabla^2 \phi = 4 *pi*G \rho
!
! With force softening from the kernel
! by a direct summation over the particles.
!
! Use this to check the accuracy of the tree code
!
! Input:
!
!   xyzh(ndim,ntot)  : coordinates and smoothing length of the particles
!   gradh(2,ntot)    : gradh and gradsoft terms (computed alongside density)
!
! Output:
!
!   phitot          : total potential phi
!   fgrav(ndim,ntot) : gravitational force
!----------------------------------------------------------------------------

subroutine directsum_grav(xyzh,gradh,fgrav,phitot,ntot)
 use kernel,    only:grkern,kernel_softening,radkern2,cnormk
 use part,      only:igas,iamtype,maxphase,maxp,iphase, &
                     iactive,isdead_or_accreted,massoftype,maxgradh
 use dim,       only:maxvxyzu,maxp
 use io,        only:error
 integer,      intent(in)    :: ntot
 real,         intent(in)    :: xyzh(4,ntot)
 real(kind=4), intent(in)    :: gradh(:,:)
 real,         intent(inout) :: fgrav(maxvxyzu,ntot)
 real,         intent(out)   :: phitot
 integer :: i,j,iamtypei,iamtypej
 real :: dx(3),dr(3),fgravi(3),fgravj(3),xi(3)
 real :: rij,rij1,rij2,rij21,pmassi,pmassj
 real :: gradhi,gradsofti,grkerni,grkernj,dsofti,dsoftj
 real :: phii,phij,phiterm,fm,fmi,fmj,phitemp,potensoft0,qi,qj
 real :: hi,hj,hi1,hj1,hi21,hj21,hi41,hj41,q2i,q2j
 logical :: iactivei,iactivej
!
!--reset potential (but not force) initially
!
! phi = 0.
 phitot = 0.

 iactivei = .true.
 iactivej = .true.
 iamtypei = igas
 iamtypej = igas
 pmassi = massoftype(iamtypei)
 pmassj = massoftype(iamtypej)

 call kernel_softening(0.,0.,potensoft0,fmi)
 if (size(gradh(:,1)) < 2) then
    call error('directsum','cannot do direct sum with ngradh < 2')
    return
 endif
 if (maxgradh /= maxp) then
    call error('directsum','gradh not stored (maxgradh /= maxp in part.F90)')
    return
 endif
!
!--calculate gravitational force by direct summation on all particles
!
 overi: do i=1,ntot
    xi(1:3) = xyzh(1:3,i)
    hi      = xyzh(4,i)
    if (isdead_or_accreted(hi)) cycle overi

    if (maxphase==maxp) then
       iamtypei = iamtype(iphase(i))
       iactivei = iactive(iphase(i))
       pmassi = massoftype(iamtypei)
    endif

    hi1  = 1./hi
    hi21 = hi1*hi1
    hi41 = hi21*hi21
    gradhi    = gradh(1,i)
    gradsofti = gradh(2,i)
    fgravi(:) = 0.
    phitemp   = 0.

    overj: do j=i+1,ntot
       dx(1) = xi(1) - xyzh(1,j)
       dx(2) = xi(2) - xyzh(2,j)
       dx(3) = xi(3) - xyzh(3,j)
       hj    = xyzh(4,j)
       if (isdead_or_accreted(hj)) cycle overj
       hj1   = 1./hj
       rij2  = dot_product(dx,dx)
       rij   = sqrt(rij2)
       rij1  = 1./rij
       rij21 = rij1*rij1
       dr(:) = dx(:)*rij1
       hj21  = hj1*hj1
       hj41  = hj21*hj21
       if (maxphase==maxp) then
          iamtypej = iamtype(iphase(j))
          iactivej = iactive(iphase(j))
          pmassj = massoftype(iamtypej)
       endif
       fgravj(:) = 0.
       q2i = rij2*hi21
       q2j = rij2*hj21
       if (q2i < radkern2) then
          qi = sqrt(q2i)
          grkerni = grkern(q2i,qi)
          call kernel_softening(q2i,qi,phii,fmi)
          phii       = phii*hi1
          fmi        = fmi*hi21
          grkerni    = cnormk*grkerni*hi41*gradhi
          dsofti     = 0.5*grkerni*gradsofti
          fgravi(:)  = fgravi(:) - pmassj*dsofti*dr(:)
          fgravj(:)  = fgravj(:) + pmassi*dsofti*dr(:)
       else
          phii = -rij1
          fmi  = rij21
       endif
       if (q2j < radkern2) then
          qj = sqrt(q2j)
          grkernj = grkern(q2j,qj)
          call kernel_softening(q2j,qj,phij,fmj)
          phij       = phij*hj1
          fmj        = fmj*hj21
          grkernj    = cnormk*grkernj*hj41*gradh(1,j)
          dsoftj     = 0.5*grkernj*gradh(2,j)
          fgravi(:)  = fgravi(:) - pmassj*dsoftj*dr(:)
          fgravj(:)  = fgravj(:) + pmassi*dsoftj*dr(:)
       else
          phij = -rij1
          fmj  = rij21
       endif

       phiterm = 0.5*(phii + phij)
       phitemp = phitemp + pmassj*phiterm
       !phi(j) = phi(j) + pmassi*phiterm
       phitot = phitot + pmassj*pmassi*phiterm

       fm = 0.5*(fmi + fmj)
       fgravi(1:3) = fgravi(1:3) - pmassj*dr(1:3)*fm
       if (iactivej) then
          fgrav(1:3,j) = fgrav(1:3,j) + pmassi*dr(1:3)*fm + fgravj(1:3)
       endif
    enddo overj
!
!--add self contribution to potential
!
    if (iactivei) then
       fgrav(1:3,i) = fgrav(1:3,i) + fgravi(1:3)
    endif
    !phi(i) = phitemp + pmassi*potensoft0*hi1
    phitot = phitot + pmassi*phitemp + pmassi*pmassi*potensoft0*hi1
 enddo overi

! phitot = 0.
! do i=1,ntot
!    phitot = phitot + 0.5*pmassi*phi(i)
! enddo
 phitot = 0.5*phitot

 return
end subroutine directsum_grav

end module directsum
