!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  Input is a relaxed star, output is two relaxed stars in binary orbit
!
!  REFERENCES: None
!
!  OWNER: Terrence Tricco
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, initial_params, part, prompting, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas
 use units,          only: set_units,udist,unit_velocity
 use prompting,      only: prompt
 use centreofmass,   only: reset_centreofmass
 use initial_params, only: get_conserv
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 integer :: opt
 real :: sep,mtot,velocity

 print *, 'Running moddump_binarystar: set up binary star systems in close contact'
 print *, ''
 print *, 'Options:'
 print *, '   1) Duplicate a relaxed star'
 print *, '   2) Adjust separation of existing binary'

 opt = 1
 call prompt('Choice',opt, 1, 2)

 if (opt  /=  1 .and. opt  /=  2) then
    print *, 'Incorrect option selected. Doing nothing.'
    return
 endif

 sep = 10.0
 print *, ''
 print *, 'Distance unit is: ', udist
 call prompt('Enter radial separation between stars (in code unit)', sep, 0.)
 print *, ''

 if (opt == 1) then
    call duplicate_star(npart, npartoftype, massoftype, xyzh, vxyzu)
 endif

 mtot = npart*massoftype(igas)
 velocity = 0.5 * sqrt(1.0 * mtot) / sqrt(sep) ! in code units

 call adjust_sep(npart, npartoftype, massoftype, xyzh, vxyzu, sep)
 call set_velocity(npart, npartoftype, massoftype, xyzh, vxyzu, velocity)

 ! reset centre of mass of the binary system
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 get_conserv = 1.

end subroutine modify_dump

subroutine duplicate_star(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas,temperature
 use units,        only: set_units,udist,unit_velocity
 use prompting,    only: prompt
 use centreofmass, only: reset_centreofmass
 use dim,          only: store_temperature
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 real :: sep,mtot,velocity

 npart = npartoftype(igas)

 sep = 10.0

 ! duplicate relaxed star
 do i = npart+1, 2*npart
    ! place star a distance rad away
    xyzh(1,i) = xyzh(1,i-npart) + sep
    xyzh(2,i) = xyzh(2,i-npart)
    xyzh(3,i) = xyzh(3,i-npart)
    xyzh(4,i) = xyzh(4,i-npart)
    vxyzu(1,i) = vxyzu(1,i-npart)
    vxyzu(2,i) = vxyzu(2,i-npart)
    vxyzu(3,i) = vxyzu(3,i-npart)
    vxyzu(4,i) = vxyzu(4,i-npart)
    if (store_temperature) then
       temperature(i) = temperature(i-npart)
    endif
    call set_particle_type(i,igas)
 enddo

 npart = 2 * npart
 npartoftype(igas) = npart

end subroutine duplicate_star

subroutine adjust_sep(npart,npartoftype,massoftype,xyzh,vxyzu,sep)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas,iamtype,iphase,maxphase,maxp
 use units,        only: set_units,udist,unit_velocity
 use prompting,    only: prompt
 use centreofmass, only: reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: sep
 integer :: i, itype
 real    :: xi, yi, zi, vxi, vyi, vzi
 real    :: x1com, y1com, z1com, x2com, y2com, z2com
 real    :: vx1com, vy1com, vz1com, vx2com, vy2com, vz2com
 real    :: totmass, pmassi, dm
 ! calc centre of mass of each star to form the reference points to adjust the position of the second star


 ! first star
 x1com = 0.
 y1com = 0.
 z1com = 0.
 vx1com = 0.
 vy1com = 0.
 vz1com = 0.
 totmass = 0.
 do i = 1, npart/2
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    vxi = vxyzu(1,i)
    vyi = vxyzu(2,i)
    vzi = vxyzu(3,i)
    if (maxphase == maxp) then
       itype = iamtype(iphase(i))
       if (itype > 0) then
          pmassi = massoftype(itype)
       else
          pmassi = massoftype(igas)
       endif
    else
       pmassi = massoftype(igas)
    endif

    totmass = totmass + pmassi
    x1com = x1com + pmassi * xi
    y1com = y1com + pmassi * yi
    z1com = z1com + pmassi * zi
    vx1com = vx1com + pmassi * vxi
    vy1com = vy1com + pmassi * vyi
    vz1com = vz1com + pmassi * vzi
 enddo

 if (totmass > tiny(totmass)) then
    dm = 1.d0/totmass
 else
    dm = 0.d0
 endif
 x1com = dm * x1com
 y1com = dm * y1com
 z1com = dm * z1com
 vx1com = dm * vx1com
 vy1com = dm * vy1com
 vz1com = dm * vz1com

 ! second star
 x2com = 0.
 y2com = 0.
 z2com = 0.
 vx2com = 0.
 vy2com = 0.
 vz2com = 0.
 totmass = 0.
 do i = npart/2 + 1, npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    vxi = vxyzu(1,i)
    vyi = vxyzu(2,i)
    vzi = vxyzu(3,i)
    if (maxphase == maxp) then
       itype = iamtype(iphase(i))
       if (itype > 0) then
          pmassi = massoftype(itype)
       else
          pmassi = massoftype(igas)
       endif
    else
       pmassi = massoftype(igas)
    endif

    totmass = totmass + pmassi
    x2com = x2com + pmassi * xi
    y2com = y2com + pmassi * yi
    z2com = z2com + pmassi * zi
    vx2com = vx2com + pmassi * vxi
    vy2com = vy2com + pmassi * vyi
    vz2com = vz2com + pmassi * vzi
 enddo

 if (totmass > tiny(totmass)) then
    dm = 1.d0/totmass
 else
    dm = 0.d0
 endif
 x2com = dm * x2com
 y2com = dm * y2com
 z2com = dm * z2com
 vx2com = dm * vx2com
 vy2com = dm * vy2com
 vz2com = dm * vz2com


 ! now we now the centre point of each star, we can set star 1 to origin, star 2 sep away on x axis, then reset com
 do i = 1, npart/2
    xyzh(1,i) = xyzh(1,i) - x1com
    xyzh(2,i) = xyzh(2,i) - y1com
    xyzh(3,i) = xyzh(3,i) - z1com
    vxyzu(1,i) = vxyzu(1,i) - vx1com
    vxyzu(2,i) = vxyzu(2,i) - vy1com
    vxyzu(3,i) = vxyzu(3,i) - vz1com
 enddo

 do i = npart/2 + 1, npart
    xyzh(1,i) = xyzh(1,i) - x2com + sep
    xyzh(2,i) = xyzh(2,i) - y2com
    xyzh(3,i) = xyzh(3,i) - z2com
    vxyzu(1,i) = vxyzu(1,i) - vx2com
    vxyzu(2,i) = vxyzu(2,i) - vy2com
    vxyzu(3,i) = vxyzu(3,i) - vz2com
 enddo

end subroutine adjust_sep


subroutine set_velocity(npart,npartoftype,massoftype,xyzh,vxyzu,velocity)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas
 use units,        only: set_units,udist,unit_velocity
 use prompting,    only: prompt
 use centreofmass, only: reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: velocity
 integer :: i
 real :: mtot

 print *, "Adding a bulk velocity |V| = ", velocity, "( = ", (velocity*unit_velocity), &
                  " physical units) to set stars in mutual orbit"
 print *, ''
 ! Adjust bulk velocity of relaxed star towards second star
 do i = 1, npart/2
    vxyzu(2,i) = vxyzu(2,i) + velocity
 enddo

 do i = npart/2 + 1, npart
    vxyzu(2,i) = vxyzu(2,i) - velocity
 enddo

end subroutine set_velocity


end module moddump

