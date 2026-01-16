!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Interactively change sink particle properties
!
! :References: None
!
! :Owner: joshcalcino
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part, prompting, ptmass_heating, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,ihsoft,ilum
 use prompting,      only:prompt
 use centreofmass,   only:reset_centreofmass
 use ptmass_heating, only:Lnuc
 use units,          only:unit_energ,utime
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),massoftype(:)

 print*,'Sink particles in dump:'
 do i=1,nptmass
    print "(a,1x,i4,a)",'Sink',i,':'
    print "(7(a5,1x,a,1x,f13.7,/))",&
             'x','=',xyzmh_ptmass(1,i),&
             'y','=',xyzmh_ptmass(2,i),&
             'z','=',xyzmh_ptmass(3,i),&
             'mass','=',xyzmh_ptmass(4,i),&
             'h','=',xyzmh_ptmass(ihsoft,i),&
             'hacc','=',xyzmh_ptmass(ihacc,i),&
             'Lnuc','=',xyzmh_ptmass(ilum,i),&
             '???','=',xyzmh_ptmass(12,i)
    if (i > 10) then
       print*, "The rest of the sink particles are not displayed"
       exit
    endif
 enddo


 ! Aa T = 3900K, L=0.44 Lsun, R=1.5 Rsun
 ! Ab1 T = 3400K, L=0.153 Lsun, R=1.13 Rsun
 ! Ab2 T = 3200K, L=0.077 Lsun, R=0.9 Rsun
 ! sink 1 is GG Tau B
 ! sink 2 is GG Tau Ab
 ! sink 3 is GG Tau Aa
 ! We will move GG Tau Aa to sink 4, and replace sink 2 and 3 with the right binary
 ! Ab1 and Ab2 have a projected separation of 4 au, just pick that
 nptmass = nptmass + 1
 xyzmh_ptmass(:,4) = xyzmh_ptmass(:,3)

 xyzmh_ptmass(:,3) = xyzmh_ptmass(:,2) ! Ab1 = Ab2
 xyzmh_ptmass(1,3) = xyzmh_ptmass(1,3) - 2.0
 xyzmh_ptmass(1,2) = xyzmh_ptmass(1,2) + 2.0 ! sep by 4au
 xyzmh_ptmass(4,2) = 0.3
 xyzmh_ptmass(4,3) = 0.2


 return
end subroutine modify_dump

end module moddump
