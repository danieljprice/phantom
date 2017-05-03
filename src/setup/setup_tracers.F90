!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  The following is an internal module to this file providing
!  interfaces to the c routines which perform the HDF5 read
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, flash_hdf5read, io, part
!+
!--------------------------------------------------------------------------
module flash_hdf5read

 !interface to the c versions
 interface
  subroutine read_flash_hdf5_header(filename,time,npart,ncol,ierr) bind(c)
   use, intrinsic :: iso_c_binding
   character(len=1),    intent(in)  :: filename
   real(kind=c_float),  intent(out) :: time
   integer(kind=c_int), intent(out) :: npart,ncol,ierr
  end subroutine read_flash_hdf5_header

  subroutine read_flash_hdf5_data(filename,npart,ncol,isrequired,ierr) bind(c)
   use, intrinsic :: iso_c_binding
   character(len=1),    intent(in)  :: filename
   integer(kind=c_int), intent(in)  :: isrequired(ncol)
   integer(kind=c_int), intent(in)  :: npart,ncol
   integer(kind=c_int), intent(out) :: ierr
  end subroutine read_flash_hdf5_data

  subroutine write_tracer_particle_density(filename,npart,rho,ierr) bind(c)
   use, intrinsic :: iso_c_binding
   character(len=1),    intent(in)  :: filename
   integer(kind=c_int), intent(in)  :: npart
   real(kind=c_double), intent(in)  :: rho(npart)
   integer(kind=c_int), intent(out) :: ierr
  end subroutine write_tracer_particle_density

 end interface

end module flash_hdf5read

!----------------------------------------------------------------
!+
!  This module sets up the particles for a Phantom run
!+
!----------------------------------------------------------------
module setup
 implicit none
 public :: setpart
!
!--dumpfile is a public variable
! that should be set by the calling routine PRIOR to each call
!
 character(len=120), public :: dumpfile

 private

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use io, only:fatal
 use dim, only:maxp
 use boundary, only:set_boundary
 use flash_hdf5read
 integer,           intent(in)  :: id
 integer,           intent(out) :: npart
 integer,           intent(out) :: npartoftype(:)
 real,              intent(out) :: xyzh(:,:)
 real,              intent(out) :: massoftype(:)
 real,              intent(out) :: vxyzu(:,:)
 real,              intent(out) :: polyk,gamma,hfact,time
 character(len=20), intent(in)  :: fileprefix
 integer :: ncolstep,ierr
 integer :: isrequired(64)
 logical :: iexist
 real(kind=4) :: tread
 real :: totmass
 !
 !--boundaries
 !
 call set_boundary(0.,1.,0.,1.,0.,1.)
 !
 !--check if first data file exists
 !
 inquire(file=dumpfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
    return
 endif

 print "(/,a)",' reading FLASH tracer particles (HDF5) data format '
 write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)

 !
 !--read the file header to get number of particles
 !
 call read_flash_hdf5_header(trim(dumpfile)//achar(0),tread,npart,ncolstep,ierr)
 print "(a,i10,a,es10.3,a,i2)",' npart = ',npart,' time = ',tread,' ncol = ',ncolstep
 if (npart > maxp) call fatal('setpart',' npart > array dimensions, edit dim_tracers.f90 and recompile')

 npartoftype(:) = 0
 npartoftype(1) = npart
 time = tread

 print "(a)",' assuming total mass of all particles is 1.0'
 totmass = 1.0
 massoftype(1) = totmass/real(npart)
 gamma = 1.0
 polyk = 1.0
!
!--now read the timestep data in the dumpfile
!  (to avoid Fortran calling C with the array, we don't actually
!   pass the dat array here - instead we get c to
!   "call back" to fill the dat array, below)
!
 isrequired(:) = 1  ! this is from SPLASH
 call read_flash_hdf5_data(trim(dumpfile)//achar(0),npart,ncolstep,isrequired(1:ncolstep),ierr)

 hfact = 1.2
 print "(a,f5.2,a)",' creating smoothing length using h =',hfact,'(m/rho)^(1/3)'

 return
end subroutine setpart
end module setup

subroutine receive_data_fromc(icol,npart,temparr,id) bind(c)
 use part, only:xyzh,vxyzu,hrho
 use io, only:real4
 integer, intent(in) :: icol,npart
 double precision, dimension(npart), intent(in) :: temparr
 integer, intent(in) :: id(npart)  ! for compatibility with c file only
 integer :: i

 select case(icol)
 case(1)
    print 10,icol,'rho'
    do i=1,npart
       xyzh(4,i) = hrho(real4(temparr(i)))
    enddo
 case(2)
    print 10,icol,'x'
    xyzh(1,1:npart) = temparr(1:npart)
 case(3)
    print 10,icol,'y'
    xyzh(2,1:npart) = temparr(1:npart)
 case(4)
    print 10,icol,'z'
    xyzh(3,1:npart) = temparr(1:npart)
 case(6)
    print 10,icol,'vx'
    vxyzu(1,1:npart) = temparr(1:npart)
 case(7)
    print 10,icol,'vy'
    vxyzu(2,1:npart) = temparr(1:npart)
 case(8)
    print 10,icol,'vz'
    vxyzu(3,1:npart) = temparr(1:npart)
 case(9)
#ifdef NO_ITERATIONS
    print 10,icol,'rho (ignored because NO_ITERATIONS set)'
#else
    print 10,icol,'rho (overwritten)'
    do i=1,npart
       xyzh(4,i) = hrho(real(temparr(i)))
    enddo
#endif
 case default
    print "(a,i2)",' WARNING: ignoring unknown/irrelevant column ',icol
 end select
10 format(' got column ',i2,'->',a)

 return
end subroutine receive_data_fromc
