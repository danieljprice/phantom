!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! ports an sphNG dump with sinks to Phantom
!
! :References: None
!
! :Owner: Alison Young
!
! :Runtime parameters:
!   - addsink : *add a sink from sink_properties.dat (in original cluster units)*
!   - do_trim : *trim off stray particles outside a radius*
!   - rmax    : *outer radius to trim to [au] (used if do_trim=T)*
!
! :Dependencies: centreofmass, dim, dynamic_dtmax, eos, infile_utils,
!   moddump_utils, part, physcon, prompting, readwrite_dumps, timestep,
!   units
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 logical :: addsink = .false.   ! add a sink from sink_properties.dat
 logical :: do_trim  = .false.  ! trim off stray particles
 real    :: rmax     = 0.       ! outer radius to trim to in au

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)

 use eos,       only:gamma
 use dim,       only:maxtypes
 use units,     only:udist,unit_velocity,print_units,set_units,utime,umass, &
                     unit_energ,set_units_extra,unit_ergg
 use part,      only:ihsoft,ihacc,nptmass,xyzmh_ptmass,vxyz_ptmass,iphase, &
                     igas,istar,iamtype,delete_particles_outside_sphere

 use physcon,         only:au,gg
 use readwrite_dumps, only:dt_read_in
 use timestep,        only:time,dt
 use dynamic_dtmax,   only:dtmax_max,dtmax_min
 use centreofmass,    only:reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real :: massoftype(:)
 real :: xyzh(:,:), vxyzu(:,:)
 integer :: iunit=26,j,npt
 integer :: i,gascount=0,sinkcount=0,othercount=0
 real    :: newutime,newuvel,temperature1,temperature2

 print*,' *** Importing sphNG dump file ***'
 !
 ! prompt only when the driver did not find a parameter file
 !
 if (prompt_for_params) call read_interactive_moddumpfile()

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
    write (iunit,'(i2)') nptmass
    do i=1,nptmass
       write (iunit,'(a)') 'xyzmh_ptmass(1:10,i)'
       write (iunit,'(10e15.6)') (xyzmh_ptmass(j,i),j=1,10)
       write (iunit,'(a)') 'vxyz_ptmass'
       write (iunit,'(3e15.6)') (vxyz_ptmass(j,i),j=1,3)
    enddo
    close(iunit)
 endif

! read sink particle from file
 if (addsink) then
    print *, 'reading sink_properties.dat'
    open(unit=iunit,file='sink_properties.dat',status='old')
    read (iunit,*) npt
    nptmass = nptmass + npt
    do i=1,npt
       read (iunit,*) ! skip line
       read (iunit,'(10e15.6)') (xyzmh_ptmass(j,i),j=1,10)
       read (iunit,*) ! skip line
       read (iunit,'(3e15.6)') (vxyz_ptmass(j,i),j=1,3)
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
 if (do_trim) then
    call delete_particles_outside_sphere((/ 0d0,0d0,0d0/),rmax,npart)
    print *, 'Particles r> ',rmax,' deleted'
 endif

!Change hacc
 do i=1,nptmass
    print *, "sink no.", i, "old hacc=", xyzmh_ptmass(ihacc,i)
    xyzmh_ptmass(ihacc,i) = 1.5
    print *, "sink no.", i, "new hacc=", xyzmh_ptmass(ihacc,i)
 enddo

 if (dt_read_in) then
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

end subroutine modify_dump

!
!---Interactively set the moddump parameters--------------------------------
!
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Do you want to add a sink - in original cluster units? (y/n)',addsink)
 call prompt('Do you want to trim? (y/n)',do_trim)
 if (do_trim) call prompt('Enter outer radius in au',rmax)

end subroutine read_interactive_moddumpfile

!
!---Write the moddump parameter file----------------------------------------
!
subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 call write_inopt(addsink,'addsink','add a sink from sink_properties.dat (in original cluster units)',iunit)
 call write_inopt(do_trim,'do_trim','trim off stray particles outside a radius',iunit)
 call write_inopt(rmax,'rmax','outer radius to trim to [au] (used if do_trim=T)',iunit)
 close(iunit)

end subroutine write_moddump

!
!---Read the moddump parameter file-----------------------------------------
!
subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 27
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(addsink,'addsink',db,errcount=nerr)
 call read_inopt(do_trim,'do_trim',db,errcount=nerr)
 call read_inopt(rmax,'rmax',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

real function calc_temp(u)
 use eos, only:gmw,gamma
 use physcon, only:atomic_mass_unit,kboltz
 use units, only:unit_ergg
 real, intent(in) :: u
 ! (gmw = mean molecular weight)
 calc_temp = atomic_mass_unit * gmw * u * unit_ergg / ( kboltz * gamma )

end function calc_temp

end module moddump

