!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for the Sedov blast wave problem
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - lattice       : *particle lattice (1=cubic,2=closepacked,3=hcp,4=random)*
!   - npartx        : *number of particles in x-direction*
!   - shuffle_parts : *relax particles by shuffling*
!
! :Dependencies: boundary, infile_utils, io, kernel, mpidomain, mpiutils,
!   options, part, physcon, prompting, setup_params, timestep, unifdis,
!   utils_shuffleparticles
!
 implicit none
 public :: setpart

 private
 !--private module variables
 integer                      :: npartx,ilattice
 logical                      :: shuffle_parts

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact_out,time,fileprefix)
 use setup_params, only:rhozero
 use unifdis,      only:set_unifdis,latticetype,i_random,i_closepacked,get_xyzmin_xyzmax_exact
 use io,           only:iprint,master,fatal
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use physcon,      only:pi
 use timestep,     only:tmax,dtmax
 use options,      only:alphau
 use kernel,       only:wkern,cnormk,radkern2,hfact_default
 use part,         only:hfact,igas,periodic,set_particle_type
 use mpiutils,     only:bcast_mpi,reduceall_mpi
 use mpidomain,    only:i_belong
 use prompting,    only:prompt
 use utils_shuffleparticles, only:shuffleparticles
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact_out
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real                             :: deltax,totmass,toten
 real                             :: enblast,gam1,uui,hsmooth,q2,r2
 integer                          :: i,maxp,maxvxyzu,ierr
 character(len=100)               :: filename
 logical                          :: iexist
!
!--general parameters
!
 time      = 0.
 hfact     = hfact_default
 hfact_out = hfact_default
!
!--setup particles
!
 maxp     = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 npartx   = 50
 ilattice = i_closepacked

 ! Read values from file if it exists, else prompt user for answers
 filename=trim(fileprefix)//'.setup'
 print "(/,1x,63('-'),1(/,1x,a),/,1x,63('-'),/)", 'Sedov Blast Wave.'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       call fatal('setup','failed to read in all the data from .setup.  Aborting')
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    call prompt('Enter number of particles in x ',npartx,8,nint((maxp)**(1/3.)))
    call prompt('Enter the type of particle lattice (1=cubic,2=closepacked,3=hcp,4=random)',ilattice,0,4)
    shuffle_parts = .false.
    if (ilattice==i_random) shuffle_parts = .true.
    call prompt('Relax particles by shuffling?',shuffle_parts)

    call write_setupfile(filename)
    stop 'rerun phantomsetup after editing .setup file'
 else
    stop
 endif
 call bcast_mpi(npartx)
 deltax = dxbound/npartx

 rhozero = 1.0
 polyk   = 0.
 if (maxvxyzu < 4) call fatal('setup','need to compile with ISOTHERMAL=no for sedov problem')
 enblast = 1.0
! rblast  = 2.*hfact*deltax


! prblast = gam1*enblast/(4./3.*pi*rblast**3)
 hsmooth = 2.*hfact*deltax
 gamma   = 5./3.
 gam1    = gamma - 1.

 call get_xyzmin_xyzmax_exact(latticetype(ilattice),xmin,xmax,ymin,ymax,zmin,zmax,ierr,deltax)
 dxbound = xmax-xmin
 dybound = ymax-ymin
 dzbound = zmax-zmin
 call set_unifdis(latticetype(ilattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,periodic,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/reduceall_mpi('+',npart)
 if (id==master) print*,' particle mass = ',massoftype(igas)

 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo
 if (shuffle_parts) then
    call shuffleparticles(iprint,npart,xyzh,massoftype(igas),duniform=rhozero,&
         is_setup=.true.,prefix=trim(fileprefix))
 endif

 toten = 0.
 do i=1,npart
    vxyzu(:,i) = 0.
    r2  = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
    q2  = r2/hsmooth**2
    uui = enblast*cnormk*wkern(q2,sqrt(q2))/hsmooth**3
    if (q2 < radkern2) then
       vxyzu(4,i) = uui
    else
       vxyzu(4,i) = 0.
    endif
    toten = toten + massoftype(igas)*uui
 enddo
!
!--normalise so energy = enblast exactly
!
 vxyzu(4,1:npart) = vxyzu(4,1:npart)*(enblast/toten)
!
!--set default runtime options for this setup (if .in file does not already exist)
!
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then
    tmax   = 0.1
    dtmax  = 0.005
    alphau = 1.
 endif

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
 write(iunit,"(a)") '# input file for Sedov Blast Wave setup routine'
 write(iunit,"(/,a)") '# particle resolution & placement'
 call write_inopt(npartx,'npartx','number of particles in x-direction',iunit)
 call write_inopt(ilattice,'lattice','particle lattice (1=cubic,2=closepacked,3=hcp,4=random)',iunit)
 call write_inopt(shuffle_parts,'shuffle_parts','relax particles by shuffling',iunit)

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
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)
 integer :: nerr

 nerr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartx,'npartx',db,errcount=nerr)
 call read_inopt(ilattice,'lattice',db,ierr,min=1,max=4,errcount=nerr)
 call read_inopt(shuffle_parts,'shuffle_parts',db,errcount=nerr)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif
 call close_db(db)

end subroutine read_setupfile

end module setup
