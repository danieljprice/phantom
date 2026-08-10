!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! None
!
! :References: None
!
! :Owner: Antoine Alaguero
!
! :Runtime parameters: None
!
! :Dependencies: io, part, prompting, units
!

 use part,         only:rhoh,igas,kill_particle,shuffle_part
 use prompting,    only:prompt
 use units,        only:umass,udist,utime
 use io,           only:fatal

 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 implicit none
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real   :: pmassi,rhoi,hi
 real   :: rho_threshold
 integer :: i, compt

 rho_threshold    = 5e-8     !in code units

 !rho_threshold = rho_threshold * 1000 / (100*100*100)  !cgs
 !rho_threshold = rho_threshold * udist*udist*udist / umass   !code units
 compt = 0
 do i=1,npart
    hi = xyzh(4,i)
    pmassi = massoftype(igas)
    rhoi = rhoh(hi,pmassi)
    ! write(*,*) rhoi      ! uncomment to have an idea of the density
    if (rhoi < rho_threshold) then
       call kill_particle(i,npartoftype)
       compt = compt+1
       !xyzh(4,i) = -abs(hi)    !call kill_particle(i,npoftype)
    endif
 enddo
 write(*,*) 'Particles deleted :', compt
 call shuffle_part(npart)
 if (npart /= sum(npartoftype)) call fatal('del_dead_part_dens','particles not conserved')

end subroutine modify_dump

end module moddump

