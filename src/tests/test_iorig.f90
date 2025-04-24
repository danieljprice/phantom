!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testiorig
!
! Unit tests particle id
!
! :References:
!
! :Owner: Mats Esseldeurs
!
! :Runtime parameters: None
!
! :Dependencies: io, part, partinject, random, testutils
!
 implicit none

 public :: test_iorig
 private

contains
!-------------------------------------------------
!+
!  unit tests of particle id bookkeeping
!+
!-------------------------------------------------
subroutine test_iorig(ntests,npass)
 use part,      only:iorig,npart,npartoftype,xyzh,vxyzu,kill_particle,shuffle_part
 use io,        only:id,master
 use testutils, only:checkval,update_test_scores,checkvalbuf,checkvalbuf_end
 use random,    only:ran2
 use partinject,        only:add_or_update_particle, update_injected_particles
 integer, intent(inout) :: ntests,npass
 integer :: i, j, iseed, ncheck, ierrmax
 integer :: nfailed(1)
 character(len=10) :: stringi, stringj

 if (id==master) write(*,"(/,a,/)") '--> TESTING PARTICLE ID'

 nfailed(1) = 0
 npart = 0
 iseed = -666

 do i = 1,100
    call add_or_update_particle(1,(/ran2(iseed), ran2(iseed), ran2(iseed)/),(/ran2(iseed), ran2(iseed), ran2(iseed)/), &
    ran2(iseed),ran2(iseed),i,npart,npartoftype,xyzh,vxyzu)
 enddo

 call checkval(npart,100,0,nfailed(1),'Check npart at start')
 call update_test_scores(ntests,nfailed,npass)

 ncheck=0
 nfailed(1)=0
 ierrmax=0
 do i = 1, 10
    call kill_particle(12, npartoftype)
    call kill_particle(3, npartoftype)
    call kill_particle(9, npartoftype)
    call kill_particle(8, npartoftype)

    call shuffle_part(npart)

    write(stringi, "(I2)") i
    call checkvalbuf(npart,100-4*i,0,'Check npart while deleting '//trim(stringi),nfailed(1),ncheck,ierrmax)
 enddo

 call checkvalbuf_end('check npart while deleting', ncheck, nfailed(1), ierrmax, 0)

 do i = npart,npart+100
    call add_or_update_particle(1,(/ran2(iseed), ran2(iseed), ran2(iseed)/),(/ran2(iseed), ran2(iseed), ran2(iseed)/), &
    ran2(iseed),ran2(iseed),i,npart,npartoftype,xyzh,vxyzu)
 enddo

 call checkval(npart,160,0,nfailed(1),'Check npart at end')
 call update_test_scores(ntests,nfailed,npass)

 ncheck=0
 nfailed(1)=0
 do i = 1, npart
    do j = i+1, npart
       write(stringi, "(I2)") i
       write(stringj, "(I2)") j
       call checkvalbuf(iorig(i)==iorig(j),.false.,&
      'Check iorig('//trim(stringi)//' != iorig('//trim(stringj)//')',nfailed(1),ncheck)
    enddo
 enddo
 call checkvalbuf_end('Check iorig',ncheck,nfailed(1))

 call update_test_scores(ntests,nfailed,npass)

 if (id==master) write(*,"(/,a)") '<-- PARTICLE ID TEST COMPLETE'

end subroutine test_iorig


end module testiorig
