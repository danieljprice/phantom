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
!   Setup for the Orszag-Tang Vortex problem in 3D
!
!  REFERENCES:
!    Orszag S. A., Tang C.-M., 1979, J. Fluid Mech., 90, 129
!    Dahlburg R. B., Picone J. M., 1989, Physics of Fluids B, 1, 2153
!    Picone J. M., Dahlburg R. B., 1991, Physics of Fluids B, 3, 29
!    Price D. J., Monaghan J. J., 2005, MNRAS, 364, 384
!
!  OWNER: James Wurster
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    betazero -- plasma beta
!    bzero    -- magnetic field amplitude
!    machzero -- Mach number
!    nx       -- number of particles in the x-direction
!    vzero    -- velocity amplitude
!    xymin    -- xmin ~ ymin
!
!  DEPENDENCIES: boundary, infile_utils, io, mpiutils, part, physcon,
!    prompting, setup_params, timestep, unifdis, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private
 !--private module variables
 integer :: nx
 real    :: xymin,machzero,betazero,vzero,bzero

contains

!----------------------------------------------------------------
!+
!  setup for 3D Orszag-Tang vortex problem
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd
 use io,           only:master
 use prompting,    only:prompt
 use mpiutils,     only:bcast_mpi
 use physcon,      only:pi,fourpi
 use timestep,     only:dtmax,tmax
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer                          :: i,ierr,maxvxyzu,maxp
 real                             :: deltax,totmass,dz
 real                             :: xymax,przero,uuzero,gam1
 character(len=26)                :: filename
 logical                          :: use_closepacked,iexist
!
!--general parameters
!
 time    = 0.
 gamma   = 5./3.
 filename= trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then
    tmax  = 1.00
    dtmax = 0.05
 endif
!
!--set particles
!
 maxp     = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 if (maxvxyzu < 4) stop 'need maxvxyzu=4 for orszag-tang setup'
!
!--set default parameters
!
 nx       = 128
 xymin    = -0.5
 betazero = 10./3.
 machzero = 1.0
 vzero    = 1.0
 bzero    = 1.0/sqrt(fourpi)
 polyk    = 0.0
 use_closepacked = .true.
!
!--read parameters from .setup or prompt from user
!
 print "(/,a)",' Setup for 3D Orszag-Tang vortex problem...'
 filename=trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    call prompt('Enter resolution (number of particles in x)',nx,8)
    call prompt('Enter xmin ~ ymin',xymin)
    call prompt('Enter initial plasma beta',betazero,0.)
    call prompt('Enter initial Mach number',machzero,0.)
    call prompt('Enter initial velocity amplitude',vzero,0.)
    call prompt('Enter initial magnetic field amplitude [default=1/sqrt(4pi)]',bzero,0.)
    !
    !--write default input file
    !
    call write_setupfile(filename)
 else
    stop
 endif
!
!--set calculated parameters
!
 przero   = 0.5*bzero**2*betazero
 rhozero  = gamma*przero*machzero
 gam1     = gamma - 1.
 uuzero   = przero/(gam1*rhozero)
 xymax    = -xymin
 deltax   = (xymax - xymin)/nx
 call bcast_mpi(nx)
 !
 print 10,betazero,machzero,bzero,rhozero,przero
10 format(/,' beta        = ',f6.3,', mach number = ',f6.3,/, &
            ' initial B   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--set boundaries
!
 if (use_closepacked) then
    dz = 2.*sqrt(6.)/nx
 else
    dz = 6.*deltax
 endif
 call set_boundary(xymin,xymax,xymin,xymax,-dz,dz)
!
!--set particle lattice
!
 if (use_closepacked) then
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)
 else
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)
 endif
 npartoftype(:) = 0
 npartoftype(1) = npart
!
!--set particle properties
!
 totmass    = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(1)

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
 write(iunit,"(/,a)") '# Box size'
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

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(nx,      'nx',      db,ierr)
 call read_inopt(xymin,   'xymin',   db,ierr)
 call read_inopt(betazero,'betazero',db,ierr)
 call read_inopt(machzero,'machzero',db,ierr)
 call read_inopt(vzero,   'vzero',   db,ierr)
 call read_inopt(bzero,   'bzero',   db,ierr)
 !
 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_sphereinbox: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile
!----------------------------------------------------------------
end module setup

