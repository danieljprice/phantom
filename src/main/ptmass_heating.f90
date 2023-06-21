!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module ptmass_heating
!
! Heating of particles around softening radius of sink particle
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: kernel, part
!

 implicit none
 public :: energ_sinkheat
 real, public :: Lnuc
 private

contains
!-----------------------------------------------------------------------
!+
!  Heat from point mass
!+
!-----------------------------------------------------------------------
subroutine energ_sinkheat(nptmass,xyzmh_ptmass,xi,yi,zi,dudtheati)
 use part,   only:ihsoft,imassenc,iLum
 use kernel, only:radkern2
 integer, intent(in) :: nptmass
 real, intent(in)    :: xi,yi,zi,xyzmh_ptmass(:,:)
 real, intent(out)   :: dudtheati
 integer             :: i
 real                :: dri2

 dudtheati = 0.
 do i = 1,nptmass
    dri2 = (xi-xyzmh_ptmass(1,i))**2 + (yi-xyzmh_ptmass(2,i))**2 + (zi-xyzmh_ptmass(3,i))**2
    if (dri2 < radkern2*xyzmh_ptmass(ihsoft,i)**2) then
       dudtheati = xyzmh_ptmass(iLum,i) / xyzmh_ptmass(imassenc,i)
    endif
 enddo

end subroutine energ_sinkheat


end module ptmass_heating
