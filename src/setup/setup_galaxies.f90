!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! This module initialises two galaxies.  It reads in
!  data from a run made using Hydra (Couchman, Thomas & Pearce 1995;
!  Thacker & Couchman 2006).  Units being read in are 1e5M_sun, kpc, 1e5yr
!  and converted to more reasonable units.
!
! :References:
!    Kuijken K., Dubinski J., 1995, MNRAS, 277, 1341
!    Widrow L. M., Dubinski J., 2005, ApJ, 631, 838
!    Widrow L. M., Pym B., Dubinski J., 2008, ApJ, 679, 1239
!    Wurster J. & Thacker R., 2013, MNRAS, 431, 2513
!    Wurster J. & Thacker R., 2013, MNRAS, 431, 539
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - lowres : *Resolution: T = low res (N~3.4e5), F = fiducial res (N~2.5e6)*
!
! :Dependencies: datafiles, dim, infile_utils, io, part, physcon, timestep,
!   units
!
 implicit none
 public :: setpart

 private
 logical :: lowres

contains
!----------------------------------------------------------------
!+
!  setup for galaxy merger from Hydra data
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxp,maxvxyzu
 use io,           only:master,fatal
 use timestep,     only:tmax,dtmax
 use physcon,      only:solarm,years,kpc
 use units,        only:set_units,udist,utime,umass
 use part,         only:igas,istar,idarkmatter,iamtype,iphase
 use datafiles,    only:find_phantom_datafile
 use infile_utils, only:get_options,infile_exists
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(in)    :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=120)               :: filename
 integer                          :: i,ndark,nstar,ngas,ierr
 real                             :: time_in,dist_in,mass_in
 real                             :: massdark,massstar,massgas
 real                             :: polykset
 real, allocatable                :: utmp(:)

 lowres = .true.

 ! get setup parameters from file
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 ! determine filename based on resolution
 if (lowres) then
    filename = find_phantom_datafile('galaxiesP.dat','galaxy_merger')
 else
    filename = find_phantom_datafile('galaxiesP25e5.dat','galaxy_merger')
 endif

 ! allocate temporary array for internal energy
 allocate(utmp(maxp),stat=ierr)
 if (ierr /= 0) call fatal('setup','not enough memory for utmp array')

 ! read galaxy data from file
 call read_galaxy_data(filename,npart,ndark,nstar,ngas,massdark,massstar,massgas,time,&
                       xyzh,vxyzu,utmp)

 npartoftype              = 0
 npartoftype(idarkmatter) = ndark
 npartoftype(istar)       = nstar
 npartoftype(igas)        = ngas
 massoftype              = 0.0
 massoftype(idarkmatter) = massdark
 massoftype(istar)       = massstar
 massoftype(igas)        = massgas

 ! units
 mass_in = 1.0d5*solarm
 dist_in = kpc
 time_in = 1.0d5*years
 call set_units(dist=10.*kpc,mass=1.0d12*solarm,G = 1.0d0)

 ! convert positions, velocities and masses to code units
 xyzh         = xyzh*dist_in/udist
 vxyzu(1:3,:) = vxyzu(1:3,:)*(utime/udist)/(time_in/dist_in)
 massoftype   = massoftype*mass_in/umass

 ! set energies (if not isothermal)
 polykset = 3.0d5*utime/udist
 polyk = polykset**2
 gamma = 5./3.
 if (maxvxyzu >= 4) then
    do i = 1,npart
       if (iamtype(iphase(i))==igas) then
          if (utmp(i)==0) then
             vxyzu(4,i) = polyk/(gamma * (gamma-1.0))
          else
             vxyzu(4,i) = utmp(i)*(utime/time_in)**2
          endif
       else
          vxyzu(4,i) = 0.0
       endif
    enddo
 endif
 print*,' polyk = ',polyk
 deallocate(utmp)

 ! set general parameters (only if not already done so)
 time = time*(time_in)/utime
 if (.not. infile_exists(fileprefix)) then
    tmax  = (0.1500d10*years)/utime
    dtmax = (0.0005d10*years)/utime
 endif

 write(*,'(1x,a,i10)')     'n_total:             ',npart
 write(*,'(1x,a,3i10)')    'n_dark,n_star,n_gas: ',ndark,nstar,ngas
 write(*,'(1x,a,3es10.2)') 'm_dark,m_star,m_gas: ',massoftype(idarkmatter),massoftype(istar),massoftype(igas)

end subroutine setpart

!----------------------------------------------------------------
!+
!  Read galaxy data from file
!+
!----------------------------------------------------------------
subroutine read_galaxy_data(filename,npart,ndark,nstar,ngas,massdark,massstar,massgas,time,&
                           xyzh,vxyzu,utmp)
 use dim,  only:maxp
 use io,   only:fatal
 use part, only:set_particle_type,idarkmatter,istar,igas
 character(len=*), intent(in)    :: filename
 integer,          intent(out)   :: npart,ndark,nstar,ngas
 real,             intent(out)   :: massdark,massstar,massgas,time
 real,             intent(out)   :: xyzh(:,:),vxyzu(:,:),utmp(:)
 integer :: i,ro,itype,ctrd,ctrs,ctrg,lu,ierr

 ! open file and read header
 open(newunit=lu,file=filename,status='old',action='read',iostat=ierr)
 if (ierr /= 0) call fatal('setup','unable to open '//trim(filename))
 read(lu,*) npart,ndark,nstar,ngas,massdark,massstar,massgas,time

 ! check array bounds
 if (npart > maxp) then
    close(lu)
    call fatal('setup','maxp too small.  Use ./phantomsetup --maxp=',ival=npart)
 endif
 !
 ! read particle data
 !
 ctrd = 0
 ctrs = 0
 ctrg = 0
 i  = 1
 ro = 0
 do while (ro==0)
    read(lu,*,iostat=ro) itype,xyzh(1:3,i),vxyzu(1:3,i),utmp(i),xyzh(4,i)
    if (ro==0) then
       if (itype== 0) then
          call set_particle_type(i,idarkmatter)
          ctrd = ctrd + 1
       elseif (itype==-1) then
          call set_particle_type(i,istar)
          ctrs = ctrs + 1
       elseif (itype== 1) then
          call set_particle_type(i,igas)
          ctrg = ctrg + 1
       else
          i = i - 1  ! read in an unknown particle type, so undo the next addition
       endif
       i = i + 1
    endif
 enddo
 close(lu)

 !
 ! set final particle count and verify counts
 !
 npart = i - 1
 print "(1x,a,/)",'setup: done reading file'
 if (ctrd/=ndark) call fatal('setup','read in incorrect number of dark matter particles')
 if (ctrs/=nstar) call fatal('setup','read in incorrect number of star particles')
 if (ctrg/=ngas ) call fatal('setup','read in incorrect number of gas particles')

end subroutine read_galaxy_data

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for galaxy merger setup routines'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(lowres,'lowres','Resolution: T = low res (N~3.4e5), F = fiducial res (N~2.5e6)',iunit)
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(lowres,'lowres',db,ierr)
 call close_db(db)
 if (ierr > 0) then
    print "(1x,a,i2,a)",'setup_galaxies: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

end module setup
