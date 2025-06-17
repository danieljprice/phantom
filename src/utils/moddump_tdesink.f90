!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Add a sink particle to the dump
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, io, part, physcon, prompting, sortutils,
!   units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part, only:xyzmh_ptmass,vxyz_ptmass,nptmass,igas,ihacc,ihsoft,rhoh
 use part, only:delete_particles_inside_radius
 use io,             only:fatal
 use prompting,      only:prompt
 use centreofmass,   only:reset_centreofmass
 use units,          only:umass,udist
 use sortutils,      only:set_r2func_origin,indexxfunc,r2func_origin
 use physcon,        only:solarr,solarm

 integer,  intent(inout) :: npart
 integer,  intent(inout) :: npartoftype(:)
 real,     intent(inout) :: massoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)

 real :: mcore,rcore,xpos(3),vpos(3)
 real :: den_all(npart),pmass,r
 integer :: j,n,location,iorder(npart),i
 !
 ! mass of gas particle
 !
 pmass = massoftype(1)
 mcore = 0.0
 !
 ! check gas particles exist
 !
 if (npart <= 0) then
    call fatal('moddump','no gas particles present in file')
 endif
 do j = 1,npart
    den_all(j) = rhoh(xyzh(4,j),pmass)
 enddo

 location = maxloc(den_all,dim=1)
 !
 ! prompt user for rcore/hsoft and check if these are reasonable values
 !
 call prompt('Enter sink radius in Rsun', rcore, 0.)
 if (rcore <= 0.0) then
    call fatal('moddump','Invalid sink radius entered')
 endif
 rcore = rcore * solarr / udist

 !
 ! sort particles by radius from dense core
 !
 xpos(:) = xyzh(1:3,location)
 vpos(:) = vxyzu(1:3,location)
 !
 ! sorting particles by radius
 !
 call set_r2func_origin(xpos(1),xpos(2),xpos(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
 !
 ! find particles within the rcore and determine mcore
 !
 do i = 1, npart
    r = sqrt(dot_product(xyzh(1:3,i)-xpos,xyzh(1:3,i)-xpos))
    if (r < rcore) then
       mcore = mcore + pmass
    endif
 enddo
 print*,'Mass of sink: ', mcore*umass/solarm, 'Msun , Radius of sink: ', rcore*udist/solarr, 'Rsun'
 !
 ! set sink particle
 !
 nptmass = 1
 n = nptmass
 xyzmh_ptmass(:,n)   = 0.
 xyzmh_ptmass(1:3,n) = xpos(:)
 xyzmh_ptmass(4,n)   = mcore
 xyzmh_ptmass(ihsoft,n) = rcore
 vxyz_ptmass(1:3, n) = vpos(:)
 !
 ! delete the gas particles from the region <= rcore
 !
 call delete_particles_inside_radius(xpos,rcore,npart,npartoftype)
 print*, 'Number of gas particles:',npart, ', Number of sink:',nptmass

end subroutine modify_dump

end module moddump
