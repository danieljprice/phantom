module cons2prim
 use cons2primsolver, only:ien_entropy
 implicit none

 public :: cons2primall,cons2primi
 public :: prim2consall,prim2consi

 private

contains

!-------------------------------------
!
!  Primitive to conservative routines
!
!-------------------------------------

subroutine prim2consall(npart,xyzh,metrics,vxyzu,dens,pxyzu,use_dens)
 use part,         only:isdead_or_accreted
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:),metrics(:,:,:,:),vxyzu(:,:)
 real,    intent(inout) :: dens(:)
 real,    intent(out) :: pxyzu(:,:)
 logical, intent(in), optional :: use_dens
 logical :: usedens
 integer :: i

!  By default, use the smoothing length to compute primitive density, and then compute the conserved variables.
!  (Alternatively, use the provided primitive density to compute conserved variables.
!   Depends whether you have prim dens prior or not.)
 if (present(use_dens)) then
    usedens = use_dens
 else
    usedens = .false.
 endif

!$omp parallel do default (none) &
!$omp shared(xyzh,metrics,vxyzu,dens,pxyzu,npart,usedens) &
!$omp private(i)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call prim2consi(xyzh(:,i),metrics(:,:,:,i),vxyzu(:,i),dens(i),pxyzu(:,i),usedens)
    endif
 enddo
!$omp end parallel do

end subroutine prim2consall

subroutine prim2consi(xyzhi,metrici,vxyzui,dens_i,pxyzui,use_dens)
 use cons2primsolver, only:primitive2conservative
 use utils_gr,        only:h2dens
 use eos,             only:equationofstate,ieos,gamma
 real, dimension(4), intent(in)  :: xyzhi, vxyzui
 real,               intent(in)  :: metrici(:,:,:)
 real, intent(inout)             :: dens_i
 real, dimension(4), intent(out) :: pxyzui
 logical, intent(in), optional   :: use_dens
 logical :: usedens
 real :: rhoi,Pi,ui,xyzi(1:3),vi(1:3),pondensi,spsoundi,densi

 !  By default, use the smoothing length to compute primitive density, and then compute the conserved variables.
 !  (Alternatively, use the provided primitive density to compute conserved variables.
 !   Depends whether you have prim dens prior or not.)
 if (present(use_dens)) then
    usedens = use_dens
 else
    usedens = .false.
 endif

 xyzi = xyzhi(1:3)
 vi   = vxyzui(1:3)
 ui   = vxyzui(4)
 if (usedens) then
    densi = dens_i
 else
    call h2dens(densi,xyzhi,metrici,vi) ! Compute dens from h
    dens_i = densi                      ! Feed the newly computed dens back out of the routine
 endif
 call equationofstate(ieos,pondensi,spsoundi,densi,xyzi(1),xyzi(2),xyzi(3),ui)
 pi = pondensi*densi
 call primitive2conservative(xyzi,metrici,vi,densi,ui,Pi,rhoi,pxyzui(1:3),pxyzui(4),ien_entropy,gamma)

end subroutine prim2consi

!-------------------------------------
!
!  Conservative to primitive routines
!
!-------------------------------------

subroutine cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens)
 use part, only:isdead_or_accreted, massoftype, igas, rhoh
 use io,   only:fatal
 integer, intent(in)    :: npart
 real,    intent(in)    :: pxyzu(:,:),xyzh(:,:),metrics(:,:,:,:)
 real,    intent(inout) :: vxyzu(:,:),dens(:)
 integer :: i, ierr

!$omp parallel do default (none) &
!$omp shared(xyzh,metrics,vxyzu,dens,pxyzu,npart,massoftype) &
!$omp private(i,ierr)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       call cons2primi(xyzh(:,i),metrics(:,:,:,i),pxyzu(:,i),vxyzu(:,i),dens(i),ierr=ierr)
       if (ierr > 0) then
          print*,' pmom =',pxyzu(1:3,i)
          print*,' rho* =',rhoh(xyzh(4,i),massoftype(igas))
          print*,' en   =',pxyzu(4,i)
          call fatal('cons2prim','could not solve rootfinding',i)
       endif
    endif
 end do
!$omp end parallel do

end subroutine cons2primall

! Note: this subroutine needs to be able to return pressure when called before
!       call to getting gr forces, since that requires pressure. Could maybe
!       get around this by calling eos somewhere along the way instead.
subroutine cons2primi(xyzhi,metrici,pxyzui,vxyzui,densi,ierr,pressure)
 use part,            only:massoftype, igas, rhoh
 use cons2primsolver, only:conservative2primitive
 use utils_gr,        only:rho2dens
 use eos,             only:equationofstate,ieos,gamma
 real, dimension(4),         intent(in)    :: xyzhi,pxyzui
 real, dimension(0:3,0:3,2), intent(in)    :: metrici
 real, dimension(4),         intent(inout) :: vxyzui
 real,    intent(inout)                    :: densi
 integer, intent(out)                      :: ierr
 real,    intent(out),  optional      :: pressure
 real    :: rhoi, p_guess, xyzi(1:3), v_guess(1:3), u_guess, pondens, spsound

 rhoi    = rhoh(xyzhi(4),massoftype(igas))
 xyzi    = xyzhi(1:3)
 v_guess = vxyzui(1:3)
 u_guess = vxyzui(4)
 call equationofstate(ieos,pondens,spsound,densi,xyzi(1),xyzi(2),xyzi(3),u_guess)
 p_guess = pondens*densi
 call conservative2primitive(xyzi,metrici,vxyzui(1:3),densi,vxyzui(4),p_guess,rhoi,pxyzui(1:3),pxyzui(4),ierr,ien_entropy,gamma)
 if (present(pressure)) pressure = p_guess

end subroutine cons2primi

end module cons2prim
