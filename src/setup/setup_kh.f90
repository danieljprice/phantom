!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for Kelvin-Helmholtz instability from Robertson et al. (2010)
!
! :References:
!   Robertson et al. (2010), MNRAS 401, 2463-2476
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - dy2    : *width of medium 2 (central medium)*
!   - nx     : *number of particles in x direction*
!   - przero : *initial constant pressure*
!   - rho1   : *density of medium 1*
!   - rho2   : *density of medium 2*
!   - v1     : *velocity of medium 1*
!   - v2     : *velocity of medium 2*
!   - xsize  : *size of the box in x-direction*
!   - ysize  : *size of the box in y-direction*
!
! :Dependencies: boundary, infile_utils, io, options, part, physcon,
!   setup_params, slab, timestep
!
 implicit none
 public :: setpart

 private
 ! input parameters
 real :: rho1   =  1.   ! density  of medium 1
 real :: rho2   =  2.   ! density  of medium 2
 real :: v1     = -0.5  ! velocity of medium 1
 real :: v2     =  0.5  ! velocity of medium 2
 real :: przero =  2.5  ! initial constant pressure
 real :: xsize  =  1.0  ! size of the box in the x-direction
 real :: ysize  =  1.0  ! size of the box in the y-direction
 real :: dy2    =  0.5  ! Width of medium 2 (i.e. central medium)
 integer :: nx = 64 ! number of particles in x direction

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:npart_total
 use io,           only:master
 use options,      only:nfulldump
 use slab,         only:set_slab,rho_func
 use boundary,     only:dxbound,dybound,dzbound
 use part,         only:igas,maxvxyzu
 use physcon,      only:pi
 use timestep,     only:dtmax,tmax
 use infile_utils, only:infile_exists,get_options
 use kernel,       only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=26)                :: filename
 integer :: i,ierr
 real    :: totmass
 procedure(rho_func), pointer :: density_func
!
!--general parameters
!
 time  = 0.
 gamma = 5./3
 polyk = 0.
 hfact = hfact_default
 if (.not. infile_exists(filename)) then
    tmax      = 2.00
    dtmax     = 0.1
    nfulldump = 1
 endif

 ! get setup parameters from file or interactive setup
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 npart = 0; npart_total = 0; npartoftype(:) = 0
 density_func => rhofunc

 call set_slab(id,master,nx,0.,xsize,0.,ysize,hfact,npart,npart_total,xyzh,npartoftype,&
               itype=igas,density_func=density_func,dir=2)

 totmass = dy2*dxbound*dzbound*rho2 + (dybound-dy2)*dxbound*dzbound*rho1
 massoftype(igas) = totmass/npart_total
 print*,' particle mass = ',massoftype(igas)

 do i=1,npart
    vxyzu(1,i) = v1 + Rfunc(xyzh(2,i))*(v2 - v1)
    vxyzu(2,i) = 0.1*sin(2.*pi*xyzh(1,i))
    vxyzu(3,i) = 0.
    if (maxvxyzu > 3) vxyzu(4,i) = przero/((gamma - 1.)*rhofunc(xyzh(2,i)))
 enddo

end subroutine setpart

!------------------------------------------
!+
!  desired density profile in y direction
!+
!------------------------------------------
real function rhofunc(y)
 real, intent(in) :: y

 rhofunc = rho1 + Rfunc(y)*(rho2 - rho1)

end function rhofunc

!---------------------------------------------------
!+
!  smoothing function, as per Robertson et al paper
!+
!---------------------------------------------------
real function Rfunc(y)
 real, parameter  :: delta = 0.05
 real, intent(in) :: y
 real :: fac1,fac2,yedgel,yedger

 yedgel = 0.5*ysize - 0.5*dy2
 yedger = 0.5*ysize + 0.5*dy2
 fac1 = (1. - 1./(1. + exp(2.*(y-yedgel)/delta)))
 fac2 = (1. - 1./(1. + exp(2.*(yedger-y)/delta)))
 Rfunc = fac1*fac2

end function Rfunc

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for Kelvin-Helmholtz instability setup'

 write(iunit,"(/,a)") '# resolution'
 call write_inopt(nx,'nx','number of particles in x direction',iunit)

 write(iunit,"(/,a)") '# physical parameters'
 call write_inopt(rho1,'rho1','density of medium 1',iunit)
 call write_inopt(rho2,'rho2','density of medium 2',iunit)
 call write_inopt(v1,'v1','velocity of medium 1',iunit)
 call write_inopt(v2,'v2','velocity of medium 2',iunit)
 call write_inopt(przero,'przero','initial constant pressure',iunit)

 write(iunit,"(/,a)") '# geometry parameters'
 call write_inopt(xsize,'xsize','size of the box in x-direction',iunit)
 call write_inopt(ysize,'ysize','size of the box in y-direction',iunit)
 call write_inopt(dy2,'dy2','width of medium 2 (central medium)',iunit)
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

 print "(a)",'reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 nerr = 0
 call read_inopt(nx,'nx',db,min=8,errcount=nerr)
 call read_inopt(rho1,'rho1',db,min=0.,errcount=nerr)
 call read_inopt(rho2,'rho2',db,min=0.,errcount=nerr)
 call read_inopt(v1,'v1',db,errcount=nerr)
 call read_inopt(v2,'v2',db,errcount=nerr)
 call read_inopt(przero,'przero',db,min=0.,errcount=nerr)
 call read_inopt(xsize,'xsize',db,min=0.,errcount=nerr)
 call read_inopt(ysize,'ysize',db,min=0.,errcount=nerr)
 call read_inopt(dy2,'dy2',db,min=0.,max=ysize,errcount=nerr)
 ierr = nerr
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
 endif
 call close_db(db)
 close(iunit)

end subroutine read_setupfile

end module setup
