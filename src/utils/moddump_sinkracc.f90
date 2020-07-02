!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION: Change accretion radius of sink particle
!
!  REFERENCES: None
!
!  OWNER: Mike Lau
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, prompting
!+
!--------------------------------------------------------------------------
module moddump
 use part,      only:xyzmh_ptmass,nptmass,ihacc,ihsoft
 use prompting, only:prompt
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 implicit none
 integer, intent(inout)               :: npart
 integer, dimension(:), intent(inout) :: npartoftype
 real, dimension(:), intent(inout)    :: massoftype
 real, dimension(:,:), intent(inout)  :: xyzh,vxyzu
 integer                              :: i,isinkpart
 real                                 :: racc,hsoft

 print*,'Sink particles in dump:'
 do i=1,nptmass
    print*,'Sink ',i,' : ','pos = (',xyzmh_ptmass(1:3,i),') ',&
           'mass = ',xyzmh_ptmass(4,i),' h = ',xyzmh_ptmass(ihsoft,i),&
           'hacc = ',xyzmh_ptmass(ihacc,i)
 enddo
 isinkpart = 2
 call prompt('Enter the sink particle number to modify:',isinkpart,1,nptmass)

 racc = xyzmh_ptmass(ihacc,isinkpart)
 call prompt('Enter new accretion radius for the sink:',racc,0.)
 print*,'Accretion radius changed to ',racc
 xyzmh_ptmass(ihacc,isinkpart) = racc

 hsoft = xyzmh_ptmass(ihsoft,isinkpart)
 call prompt('Enter new softening length for the sink:',hsoft,0.)
 print*,'Softening length changed to ',hsoft
 xyzmh_ptmass(ihsoft,isinkpart) = hsoft
 return
end subroutine modify_dump

end module moddump

