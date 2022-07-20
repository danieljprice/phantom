!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! default moddump routine: ports an sphNG dump with sinks to Phantom
!
! :References: None
!
! :Owner: Alison Young
!
! :Runtime parameters: None
!
! :Dependencies: boundary, eos, part, units
!
 implicit none

contains

subroutine modify_dump()
 use boundary, only:set_boundary
 use eos,      only:polyk
 use units,    only:udist,unit_velocity,print_units
 use part,     only:ihsoft,ihacc,nptmass,xyzmh_ptmass
 use prompting, only:prompt
 use physcon,  only:au
 !integer, intent(inout) :: nptmass
! real,    intent(inout) :: xyzmh_ptmass(:,:)
 integer :: i,isinkpart

 print*,' *** Importing sphNG dump file ***'
 print*,' the sound speed in code units is ',sqrt(polyk)
 print*,' the sound speed in cm/s       is ',sqrt(polyk)*unit_velocity

 call print_units
 print*,'Sink particles in dump:'
 do i=1,nptmass
    print*,'Sink ',i,' : ','pos = (',xyzmh_ptmass(1:3,i),') ',&
           'mass = ',xyzmh_ptmass(4,i),' h = ',xyzmh_ptmass(ihsoft,i),&
           'hacc = ',xyzmh_ptmass(ihacc,i)
 enddo
 isinkpart = 2
 !call prompt('Enter the sink particle number to modify:',isinkpart,1,nptmass)
 !Change hacc
 do i=1,nptmass
    print *, "sink no.", i, "old hacc=", xyzmh_ptmass(ihacc,i)
    xyzmh_ptmass(ihacc,i) = 1.5 * au/udist
    print *, "sink no.", i, "new hacc=", xyzmh_ptmass(ihacc,i)
 end do
 return
end subroutine modify_dump

end module moddump

