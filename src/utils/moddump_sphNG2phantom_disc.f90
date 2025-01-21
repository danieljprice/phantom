!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
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
! :Dependencies: boundary, centreofmass, dim, eos, part, physcon,
!   prompting, readwrite_dumps_fortran, timestep, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use boundary, only:set_boundary
 use eos,      only:gamma
 use dim,      only:maxtypes
 use units,    only:udist,unit_velocity,print_units,set_units,utime,umass,&
      unit_energ,set_units_extra,unit_ergg
 use part,     only:ihsoft,ihacc,nptmass,xyzmh_ptmass,vxyz_ptmass,iphase,&
      igas,istar,iamtype,delete_particles_outside_sphere
 use prompting, only:prompt
 use physcon,  only:au,gg
 use readwrite_dumps_fortran, only:dt_read_in_fortran
 use timestep, only:time,dt,dtmax_max,dtmax_min
 use centreofmass, only: reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real :: massoftype(:),rmax
 real :: xyzh(:,:), vxyzu(:,:)
 integer :: iunit=26,j,npt
 integer :: i,gascount=0,sinkcount=0,othercount=0
 real    :: newutime,newuvel,temperature1,temperature2
 character(len=1) :: trim,addsink
 character(len=25) :: junk

 print*,' *** Importing sphNG dump file ***'

 call print_units
 print *, 'setting gamma=5/3...'
 gamma = 1.6667
 print *, 'max/min u', maxval(vxyzu(4,:)), minval(vxyzu(4,:))
 temperature1 =  calc_temp(maxval(vxyzu(4,:)))
 print *, "unitener", unit_ergg
 print *, 'max/min T (K)',temperature1, calc_temp(minval(vxyzu(4,:)))
 print *,'Sink particles in dump:'
 do i=1,nptmass
    print *, 'Sink ',i,' : ','positon = (',xyzmh_ptmass(1:3,i),') ',&
           'mass = ',xyzmh_ptmass(4,i),' h = ',xyzmh_ptmass(ihsoft,i),&
           'hacc = ',xyzmh_ptmass(ihacc,i)
 enddo

! write sink particle info to file
 if (nptmass > 0) then
    open(unit=iunit,file='sink_properties.dat',status='replace')
    write (iunit,'(I2)') nptmass
    do i=1,nptmass
       write (iunit,'(A)') 'xyzmh_ptmass(1:10,i)'
       write (iunit,'(10E15.6)') (xyzmh_ptmass(j,i),j=1,10)
       write (iunit,'(A)') 'vxyz_ptmass'
       write (iunit,'(3E15.6)') (vxyz_ptmass(j,i),j=1,3)
    enddo
    close(iunit)
 endif

! read sink particle from file
 call prompt('Do you want to add a sink - in original cluster units? (y/n)',addsink)
 if (index(addsink,"y") > 0) then
    print *, 'reading sink_properties.dat'
    open(unit=iunit,file='sink_properties.dat',status='old')
    read (iunit,*) npt
    nptmass = nptmass + npt
    do i=1,npt
       read (iunit,*) junk
       read (iunit,'(10E15.6)') (xyzmh_ptmass(j,i),j=1,10)
       read (iunit,*) junk
       read (iunit,'(3E15.6)') (vxyz_ptmass(j,i),j=1,3)
    enddo
    close(iunit)
 endif

 print *, 'resetting centre of mass'
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 open(unit=iunit,file = 'Temps_in.dat',status='replace')
 do i = 1,npart
    write(iunit,"(E14.7,E14.7)") xyzh(4,i)*udist, calc_temp(vxyzu(4,i))
 enddo
 close(iunit)

 print *, 'dtmax_max,dtmax_min',dtmax_max,dtmax_min
 newutime = sqrt(au**3/(gg*umass))
 print *, "newutime/old", newutime/utime
 time = time * utime / newutime
 dt = dt * utime / newutime
 print *, "Converting units to au"
 xyzh(:,:) = xyzh(:,:) * udist / au
 newuvel = (au/newutime)
 vxyzu(1:3,:) = vxyzu(1:3,:) * unit_velocity / newuvel
 ! energy
 vxyzu(4,:) = vxyzu(4 ,:) * unit_energ / (umass * newuvel**2)
 if (nptmass > 0) then
    xyzmh_ptmass(1:3,:) =  xyzmh_ptmass(1:3,:) * udist / au
    xyzmh_ptmass(5:6,:) =  xyzmh_ptmass(5:6,:) * udist / au
    !spin angular momentum M L**2 T-1
    xyzmh_ptmass(8:10,:) =  xyzmh_ptmass(8:10,:) * (umass * udist**2/ utime) /&
         (umass * au**2/newutime)
    vxyz_ptmass(:,:) = vxyz_ptmass(:,:) * unit_velocity / newuvel
 endif

 udist = au
 utime = newutime
 call set_units(udist,umass,utime)
 call set_units_extra()
 call print_units

 print *, "Converted to au units"
 print *, "max/min x=", maxval(xyzh(1,:)), minval(xyzh(1,:))
 print *, "max/min vel", maxval(vxyzu(1,:)), minval(vxyzu(1,:))
 print *, "max/min u", maxval(vxyzu(4,:)), minval(vxyzu(4,:))
 print *, "max/min h", maxval(xyzh(4,:)), minval(xyzh(4,:))
 temperature2 = calc_temp(maxval(vxyzu(4,:)))
 print *, "unitener", unit_ergg
 print *, "max/min T (K)",temperature2, &
      calc_temp(minval(vxyzu(4,:)))

 if ((temperature1-temperature2)/temperature2 > 0.001) then
    print *, "Error energy has been changed!"
    print *, temperature1/temperature2
    stop
 endif

 ! Trim off stray particles
 call prompt('Do you want to trim? (y/n)',trim)
 if (index(trim,"y") > 0) then
    call prompt('Enter outer radius in au',rmax)
    call delete_particles_outside_sphere((/ 0d0,0d0,0d0/),rmax,npart)
    print *, 'Particles r> ',rmax,' deleted'
 endif

!Change hacc
 do i=1,nptmass
    print *, "sink no.", i, "old hacc=", xyzmh_ptmass(ihacc,i)
    xyzmh_ptmass(ihacc,i) = 1.5
    print *, "sink no.", i, "new hacc=", xyzmh_ptmass(ihacc,i)
 enddo

 if (dt_read_in_fortran) then
    print *, "****dt read in: deal with it!****"
    return
 endif

 print *, "Checking particle types..."

 do i=1, npart
    if (iamtype(iphase(i)) == igas) then
       gascount = gascount + 1
    elseif (iamtype(iphase(i)) == istar) then
       sinkcount = sinkcount + 1
    else
       othercount = othercount + 1
    endif
 enddo

 print *,'Found GAS:', gascount, 'sinks:', sinkcount, &
      'Other:', othercount, 'Total=', gascount+sinkcount+othercount
 print *, 'maxtypes:', maxtypes, 'npartoftype:', npartoftype,&
      'nptmass:', nptmass
 print *, 'gamma=', gamma
 print *, 'Timestep info:'
 print *, 'dtmax_max,dtmax_min', dtmax_max,dtmax_min
 print *, 'utime=', utime

 return
end subroutine modify_dump

real function calc_temp(u)
 use eos, only: gmw,gamma
 use physcon, only: atomic_mass_unit,kboltz
 use units, only: unit_ergg
 real, intent(in) :: u
 ! (gmw = mean molecular weight)
 calc_temp = atomic_mass_unit * gmw * u * unit_ergg / ( kboltz * gamma )
 return
end function calc_temp

end module moddump

