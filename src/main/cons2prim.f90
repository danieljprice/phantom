module cons2prim
 use cons2primsolver, only:ien_entropy
 implicit none

 public :: cons2primall
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
 use cons2primsolver, only:conservative2primitive
 use part,            only:isdead_or_accreted,massoftype,igas,rhoh
 use io,              only:fatal
 use eos,             only:equationofstate,ieos,gamma
 integer, intent(in)    :: npart
 real,    intent(in)    :: pxyzu(:,:),xyzh(:,:),metrics(:,:,:,:)
 real,    intent(inout) :: vxyzu(:,:),dens(:)
 integer :: i, ierr
 real    :: p_guess,rhoi,pondens,spsound

!$omp parallel do default (none) &
!$omp shared(xyzh,metrics,vxyzu,dens,pxyzu,npart,massoftype) &
!$omp shared(ieos,gamma) &
!$omp private(i,ierr,spsound,pondens,p_guess,rhoi)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
      ! Construct a guess for pressure (dens is already passed in and is also a guess coming in, but correct value gets passed out)
      call equationofstate(ieos,pondens,spsound,dens(i),xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))
      p_guess = pondens*dens(i)
      rhoi    = rhoh(xyzh(4,i),massoftype(igas))
      call conservative2primitive(xyzh(1:3,i),metrics(:,:,:,i),vxyzu(1:3,i),dens(i),vxyzu(4,i), &
                                  p_guess,rhoi,pxyzu(1:3,i),pxyzu(4,i),ierr,ien_entropy,gamma)
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

end module cons2prim
