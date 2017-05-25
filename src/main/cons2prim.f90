module cons2prim
implicit none

interface conservative_to_primitive
   module procedure conservative_to_primitive_split, conservative_to_primitive_combined
end interface conservative_to_primitive

interface primitive_to_conservative
   module procedure primitive_to_conservative_split, primitive_to_conservative_combined
end interface primitive_to_conservative

public :: primitive_to_conservative, conservative_to_primitive

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

subroutine primitive_to_conservative_combined(npart,xyzh,vxyzu,pxyzu)
   use part,         only:isdead_or_accreted
   use utils_gr,     only:h2dens
   use cons2prim_gr, only:primitive2conservative,get_pressure
   integer, intent(in)  :: npart
   real,    intent(in)  :: xyzh(:,:),vxyzu(:,:)
   real,    intent(out) :: pxyzu(:,:)
   integer :: i
   real    :: densi,rhoi,Pi,ui,xyzi(1:3),vi(1:3)

!$omp parallel do default (none) &
!$omp shared(xyzh,vxyzu,pxyzu,npart) &
!$omp private(i,densi,ui,pi,vi,rhoi,xyzi)
   do i=1,npart
     if (.not.isdead_or_accreted(xyzh(4,i))) then
       xyzi = xyzh(1:3,i)
       vi   = vxyzu(1:3,i)
       ui   = vxyzu(4,i)
       call h2dens(densi,xyzh(:,i),vi)
       call get_pressure(pi,densi,ui)
       call primitive2conservative(xyzi,vi,densi,ui,Pi,rhoi,pxyzu(1:3,i),pxyzu(4,i),'entropy')
     endif
   enddo
!$omp end parallel do
end subroutine primitive_to_conservative_combined

subroutine conservative_to_primitive_combined(npart,xyzh,pxyzu,vxyzu)
   use part,         only:isdead_or_accreted, massoftype, igas, rhoh
   use io,           only:fatal
   use cons2prim_gr, only:conservative2primitive,get_pressure
   use utils_gr,     only:rho2dens
   integer, intent(in)    :: npart
   real,    intent(in)    :: pxyzu(:,:),xyzh(:,:)
   real,    intent(inout) :: vxyzu(:,:)
   integer :: i, ierr
   real    :: rhoi, dens_guess, p_guess, xyzi(1:3), v_guess(1:3), u_guess

!$omp parallel do default (none) &
!$omp shared(xyzh,vxyzu,pxyzu,npart,massoftype) &
!$omp private(i,ierr,rhoi,dens_guess,p_guess,u_guess,v_guess,xyzi)
   do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
         rhoi    = rhoh(xyzh(4,i),massoftype(igas))
         xyzi    = xyzh(1:3,i)
         v_guess = vxyzu(1:3,i)
         u_guess = vxyzu(4,i)
         call rho2dens(dens_guess,rhoi,xyzi,v_guess)
         call get_pressure(p_guess,dens_guess,u_guess)
         call conservative2primitive(xyzi,v_guess,dens_guess,u_guess,p_guess,rhoi,pxyzu(1:3,i),pxyzu(4,i),ierr,'entropy')
         if (ierr > 0) then
            print*,' pmom =',pxyzu(1:3,i)
            print*,' rho* =',rhoi
            print*,' en   =',pxyzu(4,i)
            call fatal('cons2prim','could not solve rootfinding',i)
         endif
      endif
   end do
!$omp end parallel do
end subroutine conservative_to_primitive_combined

end module cons2prim
