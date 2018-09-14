!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION: Setup for the SR blast wave problem
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    npartx  -- number of particles in x-direction
!
!  DEPENDENCIES: boundary, dim, eos, infile_utils, io, kernel, mpiutils,
!    options, part, physcon, prompting, setup_params, timestep, unifdis,
!    units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private
 !--private module variables
 integer                      :: npartx
 real                         :: Pblast,Pmed,boxsize,Rblast
 real                         :: smoothfac

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxp,maxvxyzu,gr
 use setup_params, only:rhozero
 use unifdis,      only:set_unifdis
 use io,           only:master,fatal
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use physcon,      only:pi
 use timestep,     only:tmax,dtmax
 use options,      only:nfulldump
 use prompting,    only:prompt
 use kernel,       only:hfact_default
 use part,         only:igas,periodic
 use mpiutils,     only:bcast_mpi,reduceall_mpi
 use units,        only:set_units
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
 call set_units(G=1.,c=1.)

 !
 ! General parameters
 !
 time        = 0.
 hfact       = hfact_default
 rhozero     = 1.0
 gamma       = 5./3.

 !
 ! Default setup parameters
 !
 Pblast      = 100.0
 Pmed        = 0.0
 Rblast      = 0.125
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
 ! Read setup parameters from setup file if it exists
 !
 filename=trim(fileprefix)//'.setup'
 print "(/,1x,63('-'),1(/,1x,a),/,1x,63('-'),/)", 'SR Blast Wave.'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       call fatal('setup','failed to read in all the data from .setup.  Aborting')
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    call prompt(' Enter number of particles in x ',npartx,8,nint((maxp)**(1/3.)))
    call prompt(' Enter pressure in blast ',Pblast,0.)
    call prompt(' Enter pressure in medium ',Pmed,0.)
    call prompt(' Enter size of box ',boxsize,0.)
    call prompt(' Enter radius of blast ',Rblast,0.,boxsize/2.)
    call prompt(' IC smoothing factor (in terms of particle spacing) ',smoothfac)
    call write_setupfile(filename)
 endif
 call bcast_mpi(npartx)

 !
 ! Set boundaries
 !
 call set_boundary(-boxsize/2.,boxsize/2.,-boxsize/2.,boxsize/2.,-boxsize/2.,boxsize/2.)
 deltax = dxbound/npartx
 !
 ! Put particles on grid
 !
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)

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
    ublast     = Pblast/(rhozero*(gamma - 1.0))
    umed       = Pmed/(rhozero*(gamma - 1.0))
    vxyzu(4,i) = (ublast-umed)/(1. + exp((r-Rblast)/del)) + umed
 enddo

 write(*,'(2x,a,2Es11.4)')'Pressure in blast, medium: ',Pblast,Pmed
 write(*,'(2x,a, Es11.4)')'Initial blast radius: ',Rblast
 write(*,'(2x,a, Es11.4)')'Initial blast energy: ',toten

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
 write(iunit,"(/,a)") '# dimensions'
 call write_inopt(npartx, 'npartx' ,'number of particles in x-direction',iunit)
 call write_inopt(Pblast, 'Pblast' ,'pressure in blast' ,iunit)
 call write_inopt(Pmed,   'Pmed'   ,'pressure in medium',iunit)
 call write_inopt(boxsize,'boxsize','size of the box'   ,iunit)
 call write_inopt(Rblast, 'Rblast' ,'radius of blast'   ,iunit)
 call write_inopt(smoothfac, 'smoothfac' ,'IC smoothing factor (in terms of particle spacing)'   ,iunit)
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
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartx ,'npartx' ,db,ierr)
 call read_inopt(Pblast ,'Pblast' ,db,ierr)
 call read_inopt(Pmed   ,'Pmed'   ,db,ierr)
 call read_inopt(boxsize,'boxsize',db,ierr)
 call read_inopt(Rblast ,'Rblast' ,db,ierr)
 call read_inopt(smoothfac ,'smoothfac' ,db,ierr)
 call close_db(db)

end subroutine read_setupfile
!----------------------------------------------------------------
end module setup
