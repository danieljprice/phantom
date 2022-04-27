!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Give velocity perturbation to gas particles
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: part
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use partinject, only:add_or_update_particle
 use part,              only:igas
 use spherical,         only:set_sphere
 use kernel,       only:hfact_default
 use prompting,    only:prompt
 use units,        only:umass,udist
 use mpidomain,    only:i_belong
 use io,         only:id,master
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, dimension(:,:), allocatable :: xyzh_add,vxyzu_add
 real                   :: totmass,psep,hfact,totmass_sphere,r_sphere,pmass
 real    :: xp(3)
 integer                :: i,np_add,np,ipart,in_disk,r_disk
 integer(kind=8) :: nptot
 hfact = hfact_default
 totmass_sphere = 1.0
 pmass = massoftype(igas)
 in_disk = 0
 r_disk = 100.

 call prompt('Enter total mass in sphere ',totmass_sphere,0.)
 r_sphere = 4.
 call prompt('Enter radius of sphere',r_sphere,0.)
 call prompt('Is the sphere in a disk?',in_disk,0,1)
 xp = (/0.,0.,0./)
 if (in_disk==1) then
    call prompt('Enter radius on x axis for sphere centre',x_sphere,0.)
    xp = (/x_sphere, 0., 0./)
 endif

 np_add = int(totmass_sphere/pmass)
 nptot = np_add + npartoftype(igas)
 np = 0

 allocate(xyzh_add(4,np_add),vxyzu_add(4,np_add))

 call set_sphere('random',id,master,0.,r_sphere,psep,&
                 hfact,np,xyzh_add,xyz_origin=xp,&
                 exactN=.true.,np_requested=np_add,nptot=nptot,mask=i_belong)

                 ipart = npart ! The initial particle number (post shuffle)
write(*,*), "The sphere has been succesfully initialised."

if (in_disk==1) then
   vphi = sqrt((star_m)/x_sphere)
   vxyzu_add(1:3,:) = (/0.,vphi,0./)
else
   vxyzu_add(1:3,:) = (/0.,0.,0./)
end if

 do i = 1,np_add
    ! Add the particle
    ipart = ipart + 1
    call  add_or_update_particle(igas, xyzh_add(1:3,i), vxyzu_add(1:3,i), xyzh_add(4,i), &
                         vxyzu_add(4,i), ipart, npart, npartoftype, xyzh, vxyzu)
 enddo


 deallocate(xyzh_add,vxyzu_add)

 return
end subroutine modify_dump

end module moddump
