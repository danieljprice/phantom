!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for simple MHD sine wave decay test
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - npartx  : *number of particles in x-direction*
!
! :Dependencies: boundary, io, mpidomain, part, physcon, prompting,
!   setup_params, unifdis, infile_utils
!
 implicit none
 public :: setpart

 private
 !--private module variables
 integer :: npartx = 64

contains

!----------------------------------------------------------------
!+
!  Sets up Bx = sin(2*pi*x), zero magnetic and velocity field otherwise
!
!  The amplitude should decay as exp(-eta * 4 * pi^2 * t)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd,periodic,igas
 use io,           only:master
 use physcon,      only:pi
 use infile_utils, only:get_options
 use mpidomain,    only:i_belong
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: i,maxvxyzu,maxp,ierr
 real    :: bzero,przero,uuzero,gam1
 real    :: deltax,totmass
!
!--boundaries
!
 call set_boundary(-0.5,0.5,-0.5,0.5,-0.0625,0.0625)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
!
!--setup parameters
!
 bzero = 1.0
 rhozero = 2.0
 przero = 10.0
 if (maxvxyzu >= 4) then
    gamma = 5./3.
    gam1 = gamma - 1.
    uuzero = przero/(gam1*rhozero)
    polyk = 0.
 else
    gamma = 1.
    polyk = przero/rhozero
    uuzero = 1.5*polyk
 endif

 print "(/,a)",' Setup for MHD sine problem...'

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 deltax = dxbound/npartx

 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                  hfact,npart,xyzh,periodic,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(igas)

 do i=1,npart
    vxyzu(1,i) = 0.
    vxyzu(2,i) = 0.
    vxyzu(3,i) = 0.
    if (maxvxyzu >= 4) vxyzu(4,i) = uuzero
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(2,i) = Bzero*sin(2.*pi*(xyzh(1,i)-xmin))
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
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for MHD sine wave setup routine'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(npartx,'npartx','number of particles in x-direction',iunit)
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartx,'npartx',db,ierr)
 call close_db(db)

end subroutine read_setupfile

end module setup
