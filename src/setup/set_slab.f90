!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module slab
!
! This module sets up particles in a thin slab, i.e. 3D box but with
!   small aspect ratio in the z direction. Useful for performing 2D
!   test problems in 3D
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, mpidomain, unifdis
!
 implicit none
 type slab_t
    integer :: nx
    real :: rhozero
 end type slab_t

 type(slab_t) :: slab_opts
 
 public :: slab_t
 public :: set_slab
 public :: get_options_slab

 private

contains
!----------------------------------------------------------------
!+
!  setup particles in thin slab geometry, mainly useful
!  for performing 2D test problems in 3D
!+
!----------------------------------------------------------------
subroutine set_slab(id,master,nx,xmini,xmaxi,ymini,ymaxi,deltax,hfact,np,nptot,xyzh,lattice)
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound
 use unifdis,      only:set_unifdis
 use mpidomain,    only:i_belong
 integer,          intent(in)    :: id,master,nx
 integer,          intent(inout) :: np
 integer(kind=8),  intent(inout) :: nptot
 real,             intent(in)    :: xmini,xmaxi,ymini,ymaxi,hfact
 real,             intent(out)   :: xyzh(:,:),deltax

 character(len=*), intent(in), optional    :: lattice
 real :: dz
 character(len=20) :: mylattice
!
! use close packed lattice by default
!
 if (present(lattice)) then
    mylattice = lattice
 else
    mylattice = 'closepacked'
 endif
!
! set z direction boundary to achieve nx x nx x 12 particles by default
! the position of the z boundary depends on the lattice choice
!
 select case(mylattice)
 case('closepacked')
    dz = 4.*sqrt(6.)/nx
 case default
    deltax = (xmaxi-xmini)/nx
    dz = 6.*deltax
 end select
 call set_boundary(xmini,xmaxi,ymini,ymaxi,-dz,dz)
 deltax = dxbound/nx
!
! set particle lattice
!
 call set_unifdis(mylattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                  hfact,np,xyzh,.true.,nptot=nptot,mask=i_belong)

end subroutine set_slab

!----------------------------------------------------------------
!+
!  helper routine to get required options for slab setup
!+
!----------------------------------------------------------------
subroutine get_options_slab(fileprefix,id,master,nx,rhozero,ierr)
 use infile_utils, only:get_options
 integer, intent(in) :: id,master
 character(len=*), intent(in) :: fileprefix
 integer, intent(out) :: ierr
 integer, intent(inout) :: nx
 real, intent(inout) :: rhozero

 slab_opts%nx = nx
 slab_opts%rhozero = rhozero

 ! get setup parameters from file or interactive setup
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_options_slab,write_options_slab,setup_interactive_slab)

 nx = slab_opts%nx
 rhozero = slab_opts%rhozero

end subroutine get_options_slab

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_options_slab(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
  
 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for slab setup'
 call write_inopt(slab_opts%nx,'nx','number of particles in x direction',iunit)
 call write_inopt(slab_opts%rhozero,'rhozero','initial density (gives particle mass)',iunit)
 !call write_inopt(polykset,'polykset','sound speed in code units (sets polyk)',iunit)
 !call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
 close(iunit)
  
end subroutine write_options_slab
  
!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
subroutine read_options_slab(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)
  
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(slab_opts%nx,'nx',db,min=1,errcount=nerr)
 call read_inopt(slab_opts%rhozero,'rhozero',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif
  
end subroutine read_options_slab
  
!-----------------------------------------------------------------------
!+
!  Interactive setup
!+
!-----------------------------------------------------------------------
subroutine setup_interactive_slab()
 use prompting, only:prompt
  
 call prompt('enter number of particles in x direction ',slab_opts%nx,8)
 call prompt('enter density (gives particle mass)',slab_opts%rhozero,0.)
 !call prompt('enter sound speed in code units (sets polyk)',polykset,0.)
 !call prompt('select lattice type (1=cubic, 2=closepacked)',ilattice,1,2)
  
end subroutine setup_interactive_slab

end module slab
