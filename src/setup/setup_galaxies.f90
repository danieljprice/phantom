!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  This module initialises two galaxies.  It reads in
!  data from a run made using Hydra (Couchman, Thomas & Pearce 1995;
!  Thacker & Couchman 2006).  Units being read in are 1e5M_sun, kpc, 1e5yr
!  and converted to more reasonable units.
!
!  REFERENCES:
!    Kuijken K., Dubinski J., 1995, MNRAS, 277, 1341
!    Widrow L. M., Dubinski J., 2005, ApJ, 631, 838
!    Widrow L. M., Pym B., Dubinski J., 2008, ApJ, 679, 1239
!    Wurster J. & Thacker R., 2013, MNRAS, 431, 2513
!    Wurster J. & Thacker R., 2013, MNRAS, 431, 539
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    lowres -- Resolution: T = low res (N~3.4e5), F = fiducial res (N~2.5e6)
!
!  DEPENDENCIES: boundary, datafiles, dim, infile_utils, io, mpiutils,
!    part, physcon, prompting, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private
 logical :: lowres
contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxp,maxvxyzu
 use io,           only:master,fatal,warning
 use mpiutils,     only:bcast_mpi
 use timestep,     only:tmax,dtmax
 use physcon,      only:solarm,years,kpc
 use units,        only:set_units,udist,utime,umass
 use part,         only:set_particle_type,igas,istar,idarkmatter,iamtype,iphase
 use boundary,     only:set_boundary
 use prompting,    only:prompt
 use datafiles,    only:find_phantom_datafile
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
 integer                          :: i,ro,ndark,nstar,ngas,itype,ctrd,ctrs,ctrg,ierr,lu
 real                             :: time_in,dist_in,mass_in
 real                             :: massdark,massstar,massgas
 real                             :: polykset
 real, allocatable                :: utmp(:)
 logical                          :: iexist
 !
 ! Open setup file (if it exists) to determine resolution
 !
 filename=trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       call fatal('setup','failed to read in all the data from .setup.  Aborting')
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    lowres = .true.
    call prompt('Use the low resolution model (N~3.4e5; yes) or fiducial resolution model (N~2.5e6; no)?',lowres)
    call write_setupfile(filename)
 else
    stop
 endif

 if (lowres) then
    filename = find_phantom_datafile('galaxiesP.dat','galaxy_merger')
 else
    filename = find_phantom_datafile('galaxiesP25e5.dat','galaxy_merger')
 endif
 allocate(utmp(maxp),stat=ierr)
 if (ierr /= 0) call fatal('setup','not enough memory for utmp array')
 !
 !--Open file and read data
 !
 open(newunit=lu,file=filename,status='old',action='read',iostat=ierr)
 if (ierr /= 0) call fatal('setup','unable to open '//trim(filename))
 read(lu,*) npart,ndark,nstar,ngas,massdark,massstar,massgas,time
 ctrd = 0
 ctrs = 0
 ctrg = 0
 if (npart > maxp) then
    close(lu)
    call fatal('setup','maxp too small.  Make bigger than ',ival=npart)
 endif
 i  = 1
 ro = 0
 do while (ro==0)
    read(lu,*,iostat=ro) itype,xyzh(1:3,i),vxyzu(1:3,i),utmp(i),xyzh(4,i)
    if (ro==0) then
       if (itype== 0) then
          call set_particle_type(i,idarkmatter)
          ctrd = ctrd + 1
       else if (itype==-1) then
          call set_particle_type(i,istar)
          ctrs = ctrs + 1
       else if (itype== 1) then
          call set_particle_type(i,igas)
          ctrg = ctrg + 1
       else
          i = i - 1  ! read in an unknown particle type, so undo the next addition
       endif
       i = i + 1
    endif
 enddo
 close(lu)
 npart = i - 1
 print "(1x,a,/)",'setup: done reading file'
 if (ctrd/=ndark) call fatal('setup','read in incorrect number of dark matter particles')
 if (ctrs/=nstar) call fatal('setup','read in incorrect number of star particles')
 if (ctrg/=ngas ) call fatal('setup','read in incorrect number of gas particles')
 npartoftype              = 0
 npartoftype(idarkmatter) = ndark
 npartoftype(istar)       = nstar
 npartoftype(igas)        = ngas
 massoftype              = 0.0
 massoftype(idarkmatter) = massdark
 massoftype(istar)       = massstar
 massoftype(igas)        = massgas
 !
 ! Units
 !
 mass_in = 1.0d5*solarm
 dist_in = kpc
 time_in = 1.0d5*years
 call set_units(dist=10.*kpc,mass=1.0d12*solarm,G = 1.0d0)
 xyzh         = xyzh*dist_in/udist
 vxyzu(1:3,:) = vxyzu(1:3,:)*(utime/udist)/(time_in/dist_in)
 massoftype   = massoftype*mass_in/umass
 !
 ! set energies (if not isothermal)
 !
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
 call bcast_mpi(polykset)
 deallocate(utmp)
 !
 ! set general parameters (only if not already done so)
 !
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 time = time*(time_in)/utime
 if (.not. iexist) then
    tmax  = (0.1500d10*years)/utime
    dtmax = (0.0005d10*years)/utime
 endif

 write(*,'(1x,a,I10)')     'n_total:             ',npart
 write(*,'(1x,a,3I10)')    'n_dark,n_star,n_gas: ',ndark,nstar,ngas
 write(*,'(1x,a,3Es10.2)') 'm_dark,m_star,m_gas: ',massoftype(idarkmatter),massoftype(istar),massoftype(igas)

end subroutine setpart

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
    print "(1x,a,i2,a)",'Setup_galaxies: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile
!----------------------------------------------------------------
end module setup
