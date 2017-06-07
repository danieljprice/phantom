module cons2prim
 implicit none

 interface conservative_to_primitive
  module procedure conservative_to_primitive_split, conservative_to_primitive_combined
 end interface conservative_to_primitive

 interface primitive_to_conservative
  module procedure primitive_to_conservative_split, primitive_to_conservative_combined
 end interface primitive_to_conservative

 public :: primitive_to_conservative, conservative_to_primitive
 public :: conservative2primitive_combined

 private

contains


subroutine primitive_to_conservative_split(npart,xyzh,dens,v,u,P,rho,pmom,en)
 use part,         only:isdead_or_accreted
 use cons2prim_gr, only:primitive2conservative
 integer, intent(in) :: npart
 real, intent(in) :: xyzh(:,:), v(:,:)
 real, intent(in) :: dens(:),u(:)
 real, intent(out) :: pmom(:,:)
 real, intent(out) :: rho(:), en(:), P(:)
 integer :: i

!$omp parallel do default (none) &
!$omp shared(xyzh,v,dens,u,p,rho,pmom,en,npart) &
!$omp private(i)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call primitive2conservative(xyzh(1:3,i),v(1:3,i),dens(i),u(i),P(i),rho(i),pmom(1:3,i),en(i),'entropy')
    endif
 enddo
!$omp end parallel do

end subroutine primitive_to_conservative_split

subroutine conservative_to_primitive_split(npart,xyzh,rho,pmom,en,dens,v,u,P)
 use part,         only:isdead_or_accreted
 use io,           only:fatal
 use cons2prim_gr, only:conservative2primitive
 integer, intent(in) :: npart
 real, intent(in) :: pmom(:,:),xyzh(:,:)
 real, intent(in) :: rho(:),en(:)
 real, intent(inout) :: v(:,:)
 real, intent(inout) :: dens(:), u(:), P(:)
 integer :: i, ierr

!$omp parallel do default (none) &
!$omp shared(xyzh,v,dens,u,p,rho,pmom,en,npart) &
!$omp private(i,ierr)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call conservative2primitive(xyzh(1:3,i),v(1:3,i),dens(i),u(i),P(i),rho(i),pmom(1:3,i),en(i),ierr,'entropy')
       if (ierr > 0) then
          print*,' pmom =',pmom(1:3,i)
          print*,' rho* =',rho(i)
          print*,' en   =',en(i)
          call fatal('cons2prim','could not solve rootfinding',i)
       endif
    endif
 end do
!$omp end parallel do

end subroutine conservative_to_primitive_split

subroutine primitive_to_conservative_combined(npart,xyzh,vxyzu,dens,pxyzu)
 use part,         only:isdead_or_accreted
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:),vxyzu(:,:)
 real,    intent(inout) :: dens(:)
 real,    intent(out) :: pxyzu(:,:)
 integer :: i

!$omp parallel do default (none) &
!$omp shared(xyzh,vxyzu,dens,pxyzu,npart) &
!$omp private(i)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
      call primitive2conservative_combined(xyzh(:,i),vxyzu(:,i),dens(i),pxyzu(:,i))
    endif
 enddo
!$omp end parallel do

end subroutine primitive_to_conservative_combined

subroutine primitive2conservative_combined(xyzhi,vxyzui,densi,pxyzui)
 use utils_gr,     only:h2dens
 use cons2prim_gr, only:primitive2conservative,get_pressure
 real, dimension(4), intent(in)  :: xyzhi, vxyzui
 real, intent(inout)             :: densi
 real, dimension(4), intent(out) :: pxyzui
 real :: rhoi,Pi,ui,xyzi(1:3),vi(1:3)

 xyzi = xyzhi(1:3)
 vi   = vxyzui(1:3)
 ui   = vxyzui(4)
 call h2dens(densi,xyzhi,vi)
 call get_pressure(pi,densi,ui)
 call primitive2conservative(xyzi,vi,densi,ui,Pi,rhoi,pxyzui(1:3),pxyzui(4),'entropy')

end subroutine primitive2conservative_combined

subroutine conservative_to_primitive_combined(npart,xyzh,pxyzu,vxyzu,dens)
 use part, only:isdead_or_accreted, massoftype, igas, rhoh
 use io,   only:fatal
 integer, intent(in)    :: npart
 real,    intent(in)    :: pxyzu(:,:),xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:),dens(:)
 integer :: i, ierr

!$omp parallel do default (none) &
!$omp shared(xyzh,vxyzu,dens,pxyzu,npart,massoftype) &
!$omp private(i,ierr)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call conservative2primitive_combined(xyzh(:,i),pxyzu(:,i),vxyzu(:,i),dens(i),ierr)
       if (ierr > 0) then
          print*,' pmom =',pxyzu(1:3,i)
          print*,' rho* =',rhoh(xyzh(4,i),massoftype(igas))
          print*,' en   =',pxyzu(4,i)
          call fatal('cons2prim','could not solve rootfinding',i)
       endif
    endif
 end do
!$omp end parallel do

end subroutine conservative_to_primitive_combined

subroutine conservative2primitive_combined(xyzhi,pxyzui,vxyzui,densi,ierr)
 use part,         only:massoftype, igas, rhoh
 use cons2prim_gr, only:conservative2primitive,get_pressure
 use utils_gr,     only:rho2dens
 real,    dimension(4), intent(in)    :: xyzhi,pxyzui
 real,    dimension(4), intent(inout) :: vxyzui
 real, intent(inout)                  :: densi
 integer, intent(out),  optional      :: ierr
 real    :: rhoi, dens_guess, p_guess, xyzi(1:3), v_guess(1:3), u_guess

 rhoi    = rhoh(xyzhi(4),massoftype(igas))
 xyzi    = xyzhi(1:3)
 v_guess = vxyzui(1:3)
 u_guess = vxyzui(4)
 ! call rho2dens(dens_guess,rhoi,xyzi,v_guess)
 call get_pressure(p_guess,dens_guess,u_guess)
 call conservative2primitive(xyzi,vxyzui(1:3),dens_guess,vxyzui(4),p_guess,rhoi,pxyzui(1:3),pxyzui(4),ierr,'entropy')

end subroutine conservative2primitive_combined


end module cons2prim
