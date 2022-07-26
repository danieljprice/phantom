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

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use boundary, only:set_boundary
 use eos,      only:polyk
 use units,    only:udist,unit_velocity,print_units,set_units,utime,umass,&
      unit_energ,set_units_extra,unit_ergg
 
 use part,     only:ihsoft,ihacc,nptmass,xyzmh_ptmass,vxyz_ptmass
 use prompting, only:prompt
 use physcon,  only:au,gg
 use readwrite_dumps_fortran, only:dt_read_in_fortran
 use timestep, only:time,dt
 use centreofmass, only: reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real :: massoftype(:)
 real :: xyzh(:,:), vxyzu(:,:)

 !integer, intent(inout) :: nptmass
! real,    intent(inout) :: xyzmh_ptmass(:,:)
 integer :: i,isinkpart
 real    :: newutime,newuvel,temperature1,temperature2
 
 print*,' *** Importing sphNG dump file ***'
 print*,' the sound speed in code units is ',sqrt(polyk)
 print*,' the sound speed in cm/s       is ',sqrt(polyk)*unit_velocity

 call print_units
 print *, "max/min u", maxval(vxyzu(4,:)), minval(vxyzu(4,:))
 temperature1 =  calc_temp(maxval(vxyzu(4,:)))
 print *, "max/min T",temperature1, &
      calc_temp(minval(vxyzu(4,:)))
 print*,'Sink particles in dump:'
 do i=1,nptmass
    print*,'Sink ',i,' : ','pos = (',xyzmh_ptmass(1:3,i),') ',&
           'mass = ',xyzmh_ptmass(4,i),' h = ',xyzmh_ptmass(ihsoft,i),&
           'hacc = ',xyzmh_ptmass(ihacc,i)
 enddo

 print *, "resetting centre of mass"
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 
 !convert units to AU
 ! dist : * udist / au
 newutime = sqrt(au**3/(gg*umass))
 ! time : *utime / newutime
 time = time * utime / newutime
 dt = dt * utime / newutime
 print *, "Converting units to au"
 xyzh(:,:) = xyzh(:,:) * udist / au
 newuvel = (au/newutime) 
 vxyzu(1:3,:) = vxyzu(1:3,:) * unit_velocity / newuvel 
 ! energy
 vxyzu(4,:) = vxyzu(4 ,:) * unit_energ / (umass * newuvel**2) 
 xyzmh_ptmass(1:3,:) =  xyzmh_ptmass(1:3,:) * udist / au
 xyzmh_ptmass(5:6,:) =  xyzmh_ptmass(5:6,:) * udist / au
 !spin angular momentum M L**2 T-1
 xyzmh_ptmass(8:10,:) =  xyzmh_ptmass(8:10,:) * (umass * udist**2/ utime) / (umass * au**2/newutime)  
 vxyz_ptmass(:,:) = vxyz_ptmass(:,:) * unit_velocity / newuvel

 udist = au
 utime = newutime
 call set_units(udist,umass,utime)
 call set_units_extra()
 call print_units
 
 print *, "Converted to au units"
 print *, "max/min x=", maxval(xyzh(1,:)), minval(xyzh(1,:))
 print *, "max/min vel", maxval(vxyzu(1,:)), minval(vxyzu(1,:))
 print *, "max/min u", maxval(vxyzu(4,:)), minval(vxyzu(4,:))
 temperature2 = calc_temp(maxval(vxyzu(4,:)))
 print *, "max/min T",temperature2, &
      calc_temp(minval(vxyzu(4,:)))

 if (temperature1 .ne. temperature2) then
    print *, "Error energy has been changed!"
    stop
 endif
 isinkpart = 2
 !call prompt('Enter the sink particle number to modify:',isinkpart,1,nptmass)
 !Change hacc
 do i=1,nptmass
    print *, "sink no.", i, "old hacc=", xyzmh_ptmass(ihacc,i)
    xyzmh_ptmass(ihacc,i) = 1.5
    print *, "sink no.", i, "new hacc=", xyzmh_ptmass(ihacc,i)
 end do
  
 if (dt_read_in_fortran) then
    print *, "****dt read in: deal with it!****"
    return
 end if

 
 return
end subroutine modify_dump

real function calc_temp(u)
  use eos, only: gmw
  use physcon, only: atomic_mass_unit,kboltz
  use units, only: unit_ergg
  real, intent(in) :: u
  calc_temp = atomic_mass_unit * gmw * u * unit_ergg / kboltz
  return 
end function calc_temp

end module moddump

