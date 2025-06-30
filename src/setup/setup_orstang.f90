!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for the Orszag-Tang Vortex problem in 3D
!
! :References:
!    Orszag S. A., Tang C.-M., 1979, J. Fluid Mech., 90, 129
!    Dahlburg R. B., Picone J. M., 1989, Physics of Fluids B, 1, 2153
!    Picone J. M., Dahlburg R. B., 1991, Physics of Fluids B, 3, 29
!    Price D. J., Monaghan J. J., 2005, MNRAS, 364, 384
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - betazero : *plasma beta*
!   - bzero    : *magnetic field amplitude*
!   - machzero : *Mach number*
!   - nx       : *number of particles in the x-direction*
!   - vzero    : *velocity amplitude*
!   - xymin    : *xmin ~ ymin*
!
! :Dependencies: boundary, infile_utils, io, part, physcon, prompting,
!   setup_params, slab, timestep, units
!
 use physcon, only:pi
 implicit none
 public :: setpart

 private
 !--private module variables
 integer :: nx = 128
 real    :: xymin = -0.5
 real    :: betazero = 10./3.
 real    :: machzero = 1.0
 real    :: vzero = 1.0
 real    :: bzero = 1.0/sqrt(4.0*pi)

contains

!----------------------------------------------------------------
!+
!  setup for 3D Orszag-Tang vortex problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB,npart_total
 use slab,         only:set_slab
 use boundary,     only:xmin,ymin
 use part,         only:Bxyz,mhd,maxvxyzu,igas
 use io,           only:master
 use physcon,      only:pi,fourpi
 use timestep,     only:dtmax,tmax
 use infile_utils, only:get_options,infile_exists
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer                          :: i,ierr
 real                             :: xymax,przero,uuzero,gam1
!
!--general parameters
!
 time    = 0.
 gamma   = 5./3.
 if (.not. infile_exists(fileprefix)) then
    tmax  = 1.00
    dtmax = 0.05
 endif
!
!--set default parameters
!
 polyk = 0.0
 if (maxvxyzu < 4) stop 'need maxvxyzu=4 for Orszag-Tang Vortex setup'

 print "(/,a)",' Setup for 3D Orszag-Tang vortex problem...'

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

!
!--set calculated parameters
!
 przero   = 0.5*bzero**2*betazero
 rhozero  = gamma*przero*machzero
 gam1     = gamma - 1.
 uuzero   = przero/(gam1*rhozero)
 xymax    = -xymin

 print 10,betazero,machzero,bzero,rhozero,przero
10 format(/,' beta        = ',f6.3,', mach number = ',f6.3,/, &
            ' initial B   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--set particles
!
 call set_slab(id,master,nx,xymin,xymax,xymin,xymax,hfact,npart,npart_total,xyzh,&
               npartoftype,rhozero,massoftype,igas)
!
!--set particle properties
!
 do i=1,npart
    vxyzu(1,i) = -vzero*sin(2.*pi*(xyzh(2,i)-ymin))
    vxyzu(2,i) = vzero*sin(2.*pi*(xyzh(1,i)-xmin))
    vxyzu(3,i) = 0.
    vxyzu(4,i) = uuzero
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(1,i) = -Bzero*sin(2.*pi*(xyzh(2,i)-ymin))
       Bxyz(2,i) = Bzero*sin(4.*pi*(xyzh(1,i)-xmin))

       !--non-zero divB test
       !Bxyz(1,i) = 0.5/pi*sin(2.*pi*(xyzh(1,i)-xmin))
    endif
 enddo

 if (mhd) ihavesetupB = .true.

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
 write(iunit,"(a)") '# input file for the 3D Orszag-Tang vortex problem'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(nx,'nx','number of particles in the x-direction',iunit)
 write(iunit,"(/,a)") '# box size'
 call write_inopt(xymin,'xymin','xmin ~ ymin',iunit)
 write(iunit,"(/,a)") '# initial parameters'
 call write_inopt(betazero,'betazero','plasma beta',iunit)
 call write_inopt(machzero,'machzero','Mach number',iunit)
 call write_inopt(vzero,'vzero','velocity amplitude',iunit)
 call write_inopt(bzero,'bzero','magnetic field amplitude',iunit)
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

 nerr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(nx,      'nx',      db,errcount=nerr)
 call read_inopt(xymin,   'xymin',   db,errcount=nerr)
 call read_inopt(betazero,'betazero',db,errcount=nerr)
 call read_inopt(machzero,'machzero',db,errcount=nerr)
 call read_inopt(vzero,   'vzero',   db,errcount=nerr)
 call read_inopt(bzero,   'bzero',   db,errcount=nerr)
 if (nerr > 0) then
    ierr = nerr
    print "(1x,a,i2,a)",'setup_orstang: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

!----------------------------------------------------------------
!+
!  Interactive setup routine
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt

 call prompt('Enter resolution (number of particles in x)',nx,8)
 call prompt('Enter xmin ~ ymin',xymin)
 call prompt('Enter initial plasma beta',betazero,0.)
 call prompt('Enter initial Mach number',machzero,0.)
 call prompt('Enter initial velocity amplitude',vzero,0.)
 call prompt('Enter initial magnetic field amplitude [default=1/sqrt(4pi)]',bzero,0.)

end subroutine setup_interactive

end module setup
