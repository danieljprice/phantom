!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for the dustybox problem in dust-gas mixtures
!
! :References: Laibe & Price (2011), MNRAS 418, 1491
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - ilattice : *lattice type (1=cubic, 2=closepacked)*
!   - npartx   : *number of particles in x direction*
!   - polykset : *sound speed in code units (sets polyk)*
!   - rhozero  : *initial density (gives particle mass)*
!
! :Dependencies: boundary, infile_utils, io, mpidomain, part, physcon,
!   prompting, setup_params, unifdis, units
!
 use setup_params, only:rhozero
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: deltax,polykset
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:labeltype,set_particle_type,igas,idust,periodic
 use physcon,      only:pi,solarm,au
 use units,        only:set_units
 use mpidomain,    only:i_belong
 use infile_utils, only:get_options
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: totmass
 integer :: i,j,maxp,maxvxyzu,ierr
 integer :: itype,ntypes
 integer :: npart_previous
 logical, parameter :: ishift_box =.true.
!
! units (needed if physical drag is used)
!
 call set_units(mass=solarm,dist=au,G=1.d0)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 1.
 rhozero = 1.
 npartx = 64
 ilattice = 1
 polykset = 1.

 ! read setup parameters from file
 if (id==master) print "(/,a,/)",'  >>> Setting up dustybox problem <<<'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 ! setup particles
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 npart = 0
 npart_total = 0
 npartoftype(:) = 0

 ! set polytropic constant
 if (maxvxyzu < 4) then
    polyk = polykset**2
    if (id==master) print*,' polyk = ',polyk
 else
    polyk = 0.
 endif

 ntypes = 2
 overtypes: do j=1,ntypes
    if (j==1) then
       itype = igas
    else
       itype = idust + (j-2)
    endif
    if (ntypes > 1) then
       print "(/,a)",'  >>> Setting up '//trim(labeltype(itype))//' particles <<<'
    endif

    deltax = dxbound/npartx
    npart_previous = npart

    select case(ilattice)
    case(2)
       call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                        hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
    case default
       if (ishift_box .eqv. .false.) then
          call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                            hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
       else
          if (itype == igas) then
             call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                               hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
          else
             call set_unifdis('cubic',id,master,xmin+0.01*deltax,xmax+0.01*deltax,ymin+0.05*deltax, &
                              ymax+0.05*deltax,zmin+0.05*deltax,zmax+0.05*deltax,deltax, &
                              hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
             !call set_unifdis('cubic',id,master,xmin+0.5*deltax,xmax+0.5*deltax,ymin+0.5*deltax, &
             !                  ymax+0.5*deltax,zmin+0.5*deltax,zmax+0.5*deltax,deltax, &
             !                   hfact,npart,xyzh,nptot=npart_total)
             !--Use this previous setup to how bad the spline is
          endif
       endif
    end select

    !--set which type of particle it is
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)
       if (itype==igas) then
          vxyzu(1,i)   = 0.
          vxyzu(2:3,i) = 0.
       else
          vxyzu(1,i)   = 1.
          vxyzu(2:3,i) = 0.
       endif
    enddo

    npartoftype(itype) = npart - npart_previous
    if (id==master) print*,' npart = ',npart,npart_total

    totmass = rhozero*dxbound*dybound*dzbound
    massoftype(itype) = totmass/npartoftype(itype)
    if (id==master) print*,' particle mass = ',massoftype(itype)

 enddo overtypes

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
 write(iunit,"(a)") '# input file for dustybox setup'
 call write_inopt(npartx,'npartx','number of particles in x direction',iunit)
 call write_inopt(rhozero,'rhozero','initial density (gives particle mass)',iunit)
 call write_inopt(polykset,'polykset','sound speed in code units (sets polyk)',iunit)
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
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
 call read_inopt(npartx,'npartx',db,min=1,errcount=nerr)
 call read_inopt(rhozero,'rhozero',db,min=0.,errcount=nerr)
 call read_inopt(polykset,'polykset',db,min=0.,errcount=nerr)
 call read_inopt(ilattice,'ilattice',db,min=1,max=2,errcount=nerr)
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

 call prompt('enter number of particles in x direction ',npartx,1)
 call prompt('enter density (gives particle mass)',rhozero,0.)
 call prompt('enter sound speed in code units (sets polyk)',polykset,0.)
 call prompt('select lattice type (1=cubic, 2=closepacked)',ilattice,1,2)

end subroutine setup_interactive

end module setup
