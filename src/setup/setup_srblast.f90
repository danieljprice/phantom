!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for the SR blast wave problem
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - boxsize   : *size of the box*
!   - npartx    : *number of particles in x-direction*
!   - pblast    : *pressure in blast*
!   - pmed      : *pressure in medium*
!   - rblast    : *radius of blast*
!   - smoothfac : *IC smoothing factor (in terms of particle spacing)*
!
! :Dependencies: boundary, dim, infile_utils, io, kernel, mpidomain,
!   mpiutils, options, part, physcon, setup_params, timestep, unifdis,
!   units
!
 implicit none
 public :: setpart

 private
 !--private module variables
 integer :: npartx
 real    :: pblast,pmed,boxsize,rblast
 real    :: smoothfac

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu,gr
 use setup_params, only:rhozero
 use unifdis,      only:set_unifdis
 use io,           only:master,fatal
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use physcon,      only:pi
 use timestep,     only:tmax,dtmax
 use options,      only:nfulldump
 use kernel,       only:hfact_default
 use part,         only:igas,periodic
 use mpiutils,     only:reduceall_mpi
 use units,        only:set_units
 use mpidomain,    only:i_belong
 use infile_utils, only:get_options
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real                             :: deltax,totmass,toten
 real                             :: r,del,umed,ublast
 integer                          :: i,ierr
 character(len=100)               :: filename
 logical                          :: iexist

 !
 ! Quit if not properly compiled
 !
 if (.not.periodic) call fatal('setup','require PERIODIC=yes')
 if (.not.gr)       call fatal('setup','require GR=yes')
 if (maxvxyzu < 4)  call fatal('setup','require ISOTHERMAL=no')
 !
 ! Must have G=c=1 in relativity
 !
 call set_units(G=1.d0,c=1.d0)
 !
 ! General parameters
 !
 time        = 0.
 hfact       = hfact_default
 rhozero     = 1.0
 gamma       = 5./3.
 polyk       = 0.
 !
 ! Default setup parameters
 !
 pblast      = 100.0
 pmed        = 0.0
 rblast      = 0.125
 boxsize     = 1.
 npartx      = 40
 smoothfac   = 0.1
 !
 ! Infile
 !
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then
    tmax      = 0.2
    dtmax     = 0.005
    nfulldump = 1
 endif
 !
 ! Read setup parameters from setup file
 !
 print "(/,1x,63('-'),1(/,1x,a),/,1x,63('-'),/)", 'SR Blast Wave.'

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
 !
 ! Set boundaries
 !
 call set_boundary(-boxsize/2.,boxsize/2.,-boxsize/2.,boxsize/2.,-boxsize/2.,boxsize/2.)
 deltax = dxbound/npartx
 !
 ! Put particles on grid
 !
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                  hfact,npart,xyzh,periodic,mask=i_belong)

 del = smoothfac*deltax

 !
 ! Finalise particle properties
 !
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 totmass           = rhozero*dxbound*dybound*dzbound
 massoftype        = totmass/reduceall_mpi('+',npart)
 if (id==master) print*,' particle mass = ',massoftype(igas)

 toten = 0.
 do i=1,npart
    vxyzu(:,i) = 0.
    r          = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
    ublast     = pblast/(rhozero*(gamma - 1.0))
    umed       = pmed/(rhozero*(gamma - 1.0))
    vxyzu(4,i) = (ublast-umed)/(1. + exp((r-rblast)/del)) + umed
 enddo

 write(*,'(2x,a,2(1pg11.4))') 'Pressure in blast, medium: ',pblast,pmed
 write(*,'(2x,a,1pg11.4)')    'Initial blast radius: ',rblast
 write(*,'(2x,a,1pg11.4,/)')  'Initial blast energy: ',toten

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
 write(iunit,"(a)") '# input file for SR Blast Wave setup routine'
 call write_inopt(npartx, 'npartx' ,'number of particles in x-direction',iunit)
 call write_inopt(pblast, 'pblast' ,'pressure in blast' ,iunit)
 call write_inopt(pmed,   'pmed'   ,'pressure in medium',iunit)
 call write_inopt(boxsize,'boxsize','size of the box'   ,iunit)
 call write_inopt(rblast, 'rblast' ,'radius of blast'   ,iunit)
 call write_inopt(smoothfac, 'smoothfac' ,'IC smoothing factor (in terms of particle spacing)'   ,iunit)
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
 integer :: nerr

 nerr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartx ,'npartx' ,db,min=8,errcount=nerr)
 call read_inopt(pblast ,'pblast' ,db,min=0.,errcount=nerr)
 call read_inopt(pmed   ,'pmed'   ,db,min=0.,errcount=nerr)
 call read_inopt(boxsize,'boxsize',db,min=0.,errcount=nerr)
 call read_inopt(rblast ,'rblast' ,db,min=0.,errcount=nerr)
 call read_inopt(smoothfac ,'smoothfac' ,db,min=0.,errcount=nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
