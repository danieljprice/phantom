!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for simple MHD wave propagation test
! as per section 5.1 of Iwasaki (2015)
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters:
!   - npartx     : *number of particles in x-direction*
!   - plasmaB    : *plasma beta in the initial blast*
!
! :Dependencies: boundary, infile_utils, io, kernel, mpidomain, mpiutils,
!   options, part, physcon, prompting, setup_params, timestep, unifdis,
!   units
!
 implicit none
 public :: setpart

 !--private module variables
 integer :: npartx = 64
 real    :: plasmabzero = 3.

 private

contains

!----------------------------------------------------------------
!+
!  setup particles in uniform box with velocity perturbation
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,ihavesetupB,npart_total
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd,periodic,igas,maxvxyzu
 use io,           only:master,fatal
 use timestep,     only:dtmax,tmax
 use options,      only:nfulldump
 use physcon,      only:pi
 use mpiutils,     only:bcast_mpi,reduceall_mpi
 use kernel,       only:hfact_default
 use mpidomain,    only:i_belong
 use infile_utils, only:get_options
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer                          :: i,ierr
 real                             :: bzero,przero,uuzero,gam1,hzero
 real                             :: deltax,totmass
!
!--boundaries
!
 call set_boundary(-2.,2.0,-1.0,1.0,-1.0,1.0)
!
!--general parameters
!
 time      = 0.
 tmax      = 0.6
 dtmax     = 0.06
 nfulldump = 1
 hfact     = hfact_default
!
!--setup parameters
!
 rhozero     = 2.0
 przero      = 1.0
 bzero       = sqrt(2.0*przero/plasmabzero)
 if (maxvxyzu >= 4) then
    gamma = 5./3.
    gam1 = gamma - 1.
    uuzero = przero/(gam1*rhozero)
    polyk = przero/rhozero**gamma
 else
    gamma = 1.
    polyk = przero/rhozero
    uuzero = 1.5*polyk
 endif

 print "(/,a)",' Setup for MHD wave problem...'

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                 read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 call bcast_mpi(npartx)
 deltax = dxbound/npartx
 bzero  = sqrt(2.0*przero/plasmabzero)

 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                  hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype(igas) = totmass/npart_total
 print*,'npart = ',npart,' particle mass = ',massoftype(igas)

 Bxyz = 0.0
 hzero = hfact*(massoftype(igas)/rhozero)**(1./3.)
 do i=1,npart
    vxyzu(1,i) = 0.01*exp(-(xyzh(1,i)/(3.0*hzero))**2)
    vxyzu(2,i) = 0.
    vxyzu(3,i) = 0.
    if (maxvxyzu >= 4) vxyzu(4,i) = uuzero
    if (mhd) then
       Bxyz(1,i) = bzero
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

!----------------------------------------------------------------
!+
!  Interactive setup routine
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 
 call prompt(' Enter number of particles in x ',npartx,8,nint((1000000)**(1/3.)))
 call prompt(' Enter initial plasma beta (this will adjust the magnetic field strength) ',plasmabzero)

end subroutine setup_interactive

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
 write(iunit,"(a)") '# input file for MHD wave setup routine'
 write(iunit,"(/,a)") '# dimensions'
 call write_inopt(npartx,'npartx','number of particles in x-direction',iunit)
 write(iunit,"(/,a)") '# magnetic field strength'
 call write_inopt(plasmabzero,'plasmaB','initial plasma beta',iunit)
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
 call read_inopt(plasmabzero,'plasmaB',db,ierr)
 call close_db(db)

end subroutine read_setupfile

end module setup
