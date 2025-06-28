!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for the MHD current loop advection problem
!
! :References:
!    Gardiner and Stone (2005)
!    Rosswog and Price (2007)
!    Stone et al. (2008)
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Azero  : *amplitude of vector potential*
!   - nx     : *resolution (number of particles in x)*
!   - przero : *initial pressure*
!   - rloop  : *radius of current loop*
!   - vzero  : *initial velocity*
!
! :Dependencies: boundary, dim, infile_utils, io, part, physcon, prompting,
!   setup_params, slab
!
 implicit none
 public :: setpart

 private

 ! Module variables for setup parameters
 real :: vzero,przero,Azero,rloop
 integer :: nx

contains

!----------------------------------------------------------------
!+
!  setup for MHD current loop advection problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu
 use setup_params, only:rhozero,ihavesetupB,npart_total
 use slab,         only:set_slab
 use boundary,     only:dxbound,dybound
 use part,         only:Bxyz,mhd,igas
 use io,           only:master
 use physcon,      only:pi
 use infile_utils, only:get_options
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(in)    :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: i,ierr
 real    :: rcyl,gam1,uuzero,costheta,sintheta

 ! set default parameters
 vzero = sqrt(5.)
 przero = 1.0
 Azero = 1.e-3
 rloop = 0.3
 nx = 128
 rhozero = 1.0

 if (id==master) print "(/,a)",' MHD current loop advection problem'

 ! get setup parameters from file or interactive setup
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 ! general parameters
 time = 0.
 polyk = 0.
 gamma = 5./3.
 gam1 = gamma - 1.
 uuzero = przero/(gam1*rhozero)
 if (maxvxyzu < 4) polyk = przero/rhozero**gamma

 ! print key parameters
 if (id==master) then
    print "(/,' Azero   = ',f6.3,', radius of loop = ',f6.3,/,"// &
          "' density = ',f6.3,',       pressure = ',f6.3,/)",&
          Azero,rloop,rhozero,przero
 endif

 ! setup particles and boundaries for slab geometry
 call set_slab(id,master,nx,-1.,1.,-0.5,0.5,hfact,npart,npart_total,xyzh,&
               npartoftype,rhozero,massoftype,igas)

 costheta = dxbound/sqrt(dxbound**2 + dybound**2)
 sintheta = dybound/sqrt(dxbound**2 + dybound**2)

 do i=1,npart
    vxyzu(1,i) = vzero*costheta
    vxyzu(2,i) = vzero*sintheta
    vxyzu(3,i) = 0.1*vzero
    if (maxvxyzu >= 4) vxyzu(4,i) = uuzero
    if (mhd) then
       Bxyz(:,i) = 0.
       rcyl = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
       if (rcyl < rloop) then
          Bxyz(1,i) = Azero*(-xyzh(2,i)/rcyl)
          Bxyz(2,i) = Azero*(xyzh(1,i)/rcyl)
       endif
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for MHD current loop advection setup'
 call write_inopt(vzero,'vzero','initial velocity',iunit)
 call write_inopt(przero,'przero','initial pressure',iunit)
 call write_inopt(Azero,'Azero','amplitude of vector potential',iunit)
 call write_inopt(rloop,'rloop','radius of current loop',iunit)
 call write_inopt(nx,'nx','resolution (number of particles in x)',iunit)
 close(iunit)

end subroutine write_setupfile

!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(vzero,'vzero',db,errcount=nerr)
 call read_inopt(przero,'przero',db,min=0.,errcount=nerr)
 call read_inopt(Azero,'Azero',db,errcount=nerr)
 call read_inopt(rloop,'rloop',db,min=0.,errcount=nerr)
 call read_inopt(nx,'nx',db,min=8,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Interactive setup
!+
!-----------------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt

 call prompt('Enter resolution (number of particles in x)',nx,8)

end subroutine setup_interactive

end module setup
