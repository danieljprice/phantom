!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for realistic cosmological initial conditions based
! on the Zeldovich approximation.
! Requries velocity files generated from a powerspectrum.
!
! :References: None
!
! :Owner: Spencer Magnall
!
! :Runtime parameters:
!   - cs0                 : *initial sound speed in code units*
!   - ilattice            : *lattice type (1=cubic, 2=closepacked)*
!   - nx                  : *number of particles in x direction*
!   - radiation_dominated : *Radiation dominated universe (yes/no)*
!   - rhozero             : *initial density in code units*
!
! :Dependencies: boundary, dim, eos_shen, infile_utils, io, mpidomain,
!   part, physcon, setunits, setup_params, stretchmap, unifdis, units,
!   utils_gr
!
 use dim,          only:use_dust,gr
 use setup_params, only:rhozero
 use physcon,      only:radconst
 implicit none
 public :: setpart

 integer           :: npartx,ilattice
 real              :: cs0,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,ampl,phaseoffset
 character(len=20) :: dist_unit,mass_unit,perturb_direction,perturb,radiation_dominated
 real              :: perturb_wavelength

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu,gr
 use setup_params, only:npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis,rho_func!,mass_func
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary,cross_boundary
 use part,         only:periodic
 use physcon,      only:years,pc,solarm
 use units,        only:set_units
 use mpidomain,    only:i_belong
 use stretchmap,   only:set_density_profile
 use utils_gr,     only:perturb_metric,get_u0,get_sqrtg
 use infile_utils, only:get_options
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(inout) :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=40) :: lattice,pspec_filename1,pspec_filename2,pspec_filename3
 real    :: totmass,deltax,pi
 integer :: i,ierr,ncross
 logical :: isperiodic(3)
 real    :: length, c1,c3
 real    :: hub
 real    :: last_scattering_temp
 real    :: scale_factor,gradphi(3),vxyz(3),dxgrid,gridorigin
 integer :: nghost, gridres, gridsize
 real, allocatable :: vxgrid(:,:,:),vygrid(:,:,:),vzgrid(:,:,:)
 !
 !--general parameters
 !
 perturb_wavelength = 0.
 time = 0.
 if (maxvxyzu < 4) then
    gamma = 1.
 else
    ! 4/3 for radiation dominated case
    ! irrelevant for
    gamma = 4./3.
 endif
 ! Redefinition of pi to fix numerical error
 pi = 4.*atan(1.)
 !
 ! default units
 !
 mass_unit = 'solarm'
 dist_unit = 'mpc'
 !
 ! set boundaries to default values
 !
 xmini = xmin; xmaxi = xmax
 ymini = ymin; ymaxi = ymax
 zmini = zmin; zmaxi = zmax
 !
 ! set default values for input parameters
 !
 npartx   = 64
 length = 0.
 ilattice = 1
 perturb  = '"no"'
 perturb_direction = '"none"'
 radiation_dominated = '"no"'
 ampl = 0.

 ! Ideally this should read the values of the box length
 ! and initial Hubble parameter from the par file.
 ! Then it should be set using the Friedmann equation:
 !!!!!! rhozero = (3H^2)/(8*pi*a*a)

 hub = 10.553495658357338
 rhozero = 3. * hub**2 / (8. * pi)
 phaseoffset = 0.

 ! Set some default values for the grid
 nghost = 6
 gridres = 64

 gridsize = nghost + gridres
 gridorigin = 0.
 xmax = 1.
 dxgrid = xmax/gridres
 gridorigin = gridorigin-3*dxgrid

 isperiodic = .true.
 ncross = 0

 ! Approx Temp of the CMB in Kelvins
 last_scattering_temp = 3000
 last_scattering_temp = (rhozero/radconst)**(1./4.)*0.99999

 ! Define some parameters for Linear pertubations
 ! We assume ainit = 1, but this may not always be the case
 c1 = 1./(4.*PI*rhozero)
 !c2 = We set g(x^i) = 0 as we only want to extract the growing mode
 c3 = - sqrt(1./(6.*PI*rhozero))

 if (gr) then
    ! 0 Because dust?
    cs0 = 0.
 else
    cs0 = 1.
 endif
 !
 ! get setup parameters from file or interactive setup
 !
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)

 allocate(vxgrid(gridsize,gridsize,gridsize))
 allocate(vygrid(gridsize,gridsize,gridsize))
 allocate(vzgrid(gridsize,gridsize,gridsize))
 !
 ! setup particles
 !
 npart = 0
 npart_total = 0
 length = xmaxi - xmini
 deltax = length/npartx
!
! general parameters
!
! time should be read in from the par file
 time   = 0.18951066686763596 ! z~1000
 rhozero = 3. * hub**2 / (8. * pi)
 lattice = 'cubic'

 call set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,&
      npart,xyzh,periodic,nptot=npart_total,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype(1) = totmass/npart_total
 if (id==master) print*,' particle mass = ',massoftype(1)
 if (id==master) print*,' initial sound speed = ',cs0,' pressure = ',cs0**2/gamma

 if (maxvxyzu < 4 .or. gamma <= 1.) then
    polyk = cs0**2
 else
    polyk = 0.
 endif

 pspec_filename1 = 'init_vel1_64.dat'
 pspec_filename2 = 'init_vel2_64.dat'
 pspec_filename3 = 'init_vel3_64.dat'

 ! Check if files exist otherwise skip and return flat space
 if (.not. check_files(pspec_filename1,pspec_filename2,pspec_filename3)) then
    print*, "Velocity files not found..."
    print*, "Setting up flat space!"
    return
 endif

 call read_veldata(vxgrid,pspec_filename1,gridsize)
 call read_veldata(vygrid,pspec_filename2,gridsize)
 call read_veldata(vzgrid,pspec_filename3,gridsize)

 do i=1,npart
    ! Assign new particle possition + particle velocities here using the Zeldovich approximation:
    ! Valid for Omega = 1
    ! x = q - a grad phi (1), where q is the non perturbed lattice point position
    ! v = -aH grad phi (2)
    ! Get x velocity at particle position
    call interpolate_val(xyzh(1:3,i),vxgrid,gridsize,gridorigin,dxgrid,vxyz(1))
    ! Get y velocity at particle position
    call interpolate_val(xyzh(1:3,i),vygrid,gridsize,gridorigin,dxgrid,vxyz(2))
    ! Get z velocity at particle position
    call interpolate_val(xyzh(1:3,i),vzgrid,gridsize,gridorigin,dxgrid,vxyz(3))

    vxyzu(1:3,i)  = vxyz
    ! solve eqn (2) for grad phi
    ! This is probally not constant??
    scale_factor = 1.
    gradphi = -vxyz/(scale_factor*hub)
    ! Set particle pos
    xyzh(1:3,i) = xyzh(1:3,i) - scale_factor*gradphi
    ! Apply periodic boundary conditions to particle position
    call cross_boundary(isperiodic,xyzh(1:3,i),ncross)
 enddo

end subroutine setpart

!------------------------------------------------------------------------
!
! write setup file
!
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setunits,     only:write_options_units
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(/,a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for uniform setup routine'
 !
 ! units
 !
 call write_options_units(iunit,gr)
 !
 ! boundaries
 !
 write(iunit,"(/,a)") '# boundaries'
 call write_inopt(xmini,'CoordBase::xmin','xmin boundary',iunit)
 call write_inopt(xmaxi,'CoordBase::xmax','xmax boundary',iunit)
 call write_inopt(ymini,'CoordBase::ymin','ymin boundary',iunit)
 call write_inopt(ymaxi,'CoordBase::ymax','ymax boundary',iunit)
 call write_inopt(zmini,'CoordBase::zmin','zmin boundary',iunit)
 call write_inopt(zmaxi,'CoordBase::zmax','zmax boundary',iunit)
 !
 ! other parameters
 !
 write(iunit,"(/,a)") '# setup'
 call write_inopt(npartx,'nx','number of particles in x direction',iunit)
 call write_inopt(rhozero,'rhozero','initial density in code units',iunit)
 call write_inopt(cs0,'cs0','initial sound speed in code units',iunit)
 call write_inopt(perturb,'FLRWSolver::FLRW_perturb','Pertrubations of FLRW?',iunit)
 call write_inopt(ampl,'FLRWSolver::phi_amplitude','Pertubation amplitude',iunit)
 call write_inopt(phaseoffset,'FLRWSolver::phi_phase_offset','Pertubation phase offset',iunit)
 call write_inopt(perturb_direction, 'FLRWSolver::FLRW_perturb_direction','Pertubation direction',iunit)
 call write_inopt(radiation_dominated, 'radiation_dominated','Radiation dominated universe (yes/no)',iunit)
 call write_inopt(perturb_wavelength,'FLRWSolver::single_perturb_wavelength','Perturbation wavelength',iunit)
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!
! read setup file
!
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use setunits,     only:read_options_and_set_units
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 !
 ! units
 !
 call read_options_and_set_units(db,nerr,gr)
 !
 ! boundaries
 !
 call read_inopt(xmini,'CoordBase::xmin',db,errcount=nerr)
 call read_inopt(xmaxi,'CoordBase::xmax',db,min=xmini,errcount=nerr)
 call read_inopt(ymini,'CoordBase::ymin',db,errcount=nerr)
 call read_inopt(ymaxi,'CoordBase::ymax',db,min=ymini,errcount=nerr)
 call read_inopt(zmini,'CoordBase::zmin',db,errcount=nerr)
 call read_inopt(zmaxi,'CoordBase::zmax',db,min=zmini,errcount=nerr)
 !
 ! other parameters
 !
 call read_inopt(npartx,'nx',db,min=8,errcount=nerr)
 call read_inopt(rhozero,'rhozero',db,min=0.,errcount=nerr)
 call read_inopt(cs0,'cs0',db,min=0.,errcount=nerr)

 call read_inopt(perturb_direction,'FLRWSolver::FLRW_perturb_direction',db,errcount=nerr)
 call read_inopt(ampl, 'FLRWSolver::phi_amplitude',db,errcount=nerr)
 call read_inopt(phaseoffset,'FLRWSolver::phi_phase_offset',db,errcount=nerr)
 call read_inopt(ilattice,'ilattice',db,min=1,max=2,errcount=nerr)
 ! TODO Work out why this doesn't read in correctly
 call read_inopt(perturb,'FLRWSolver::FLRW_perturb',db,errcount=nerr)
 call read_inopt(radiation_dominated,'radiation_dominated',db,errcount=nerr)
 call read_inopt(perturb_wavelength,'FLRWSolver::single_perturb_wavelength',db,errcount=nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

subroutine read_veldata(velarray,vfile,gridsize)
 integer, intent(in) :: gridsize
 character(len=20),intent(in) :: vfile
 real,intent(out) :: velarray(:,:,:)
 integer :: i,j,k,iu

 open(newunit=iu,file=vfile,status='old')
 do k=1,gridsize
    do j=1,gridsize
       read(iu,*) (velarray(i,j,k), i=1, gridsize)
    enddo
 enddo
 close(iu)
 print*, "Finished reading ", vfile

end subroutine read_veldata

subroutine interpolate_val(position,valgrid,gridsize,gridorigin,dxgrid,val)
 ! Subroutine to interpolate quanities to particle positions given a cube
 ! Note we have assumed that the grid will always be cubic!!!!
 use eos_shen, only:linear_interpolator_one_d
 real, intent(in)    :: valgrid(:,:,:)
 real, intent(inout)    :: position(3)
 real, intent(inout) :: dxgrid,gridorigin
 integer, intent(in) :: gridsize
 real, intent(out)   :: val
 integer :: xupper,yupper,zupper,xlower,ylower,zlower
 real    :: xlowerpos,ylowerpos,zlowerpos!,xupperpos,yupperpos,zupperpos
 real    :: interptmp(7)
 real    :: xd,yd,zd

 call get_grid_neighbours(position,gridorigin,dxgrid,xlower,ylower,zlower)

 print*,"Neighbours: ", xlower,ylower,zlower
 print*,"Position: ", position
 ! This is not true as upper neighbours on the boundary will be on the side
 ! take a mod of grid size
 xupper = mod(xlower + 1, gridsize)
 yupper = mod(ylower + 1, gridsize)
 zupper = mod(zlower + 1, gridsize)
 ! xupper - xlower should always just be dx provided we are using a uniform grid
 ! xd = (position(1) - xlower)/(xupper - xlower)
 ! yd = (position(2) - ylower)/(yupper - ylower)
 ! zd = (position(3) - zlower)/(zupper - zlower)
 xlowerpos = gridorigin + (xlower-1)*dxgrid
 ylowerpos = gridorigin + (ylower-1)*dxgrid
 zlowerpos = gridorigin + (zlower-1)*dxgrid

 xd = (position(1) - xlowerpos)/(dxgrid)
 yd = (position(2) - ylowerpos)/(dxgrid)
 zd = (position(3) - zlowerpos)/(dxgrid)

 interptmp = 0.

 call linear_interpolator_one_d(valgrid(xlower,ylower,zlower), &
                valgrid(xlower+1,ylower,zlower),xd,interptmp(1))
 call linear_interpolator_one_d(valgrid(xlower,ylower,zlower+1), &
                valgrid(xlower+1,ylower,zlower+1),xd,interptmp(2))
 call linear_interpolator_one_d(valgrid(xlower,ylower+1,zlower), &
                valgrid(xlower+1,ylower+1,zlower),xd,interptmp(3))
 call linear_interpolator_one_d(valgrid(xlower,ylower+1,zlower+1), &
                valgrid(xlower+1,ylower+1,zlower+1),xd,interptmp(4))
 ! Interpolate along y
 call linear_interpolator_one_d(interptmp(1),interptmp(3),yd,interptmp(5))
 call linear_interpolator_one_d(interptmp(2),interptmp(4),yd,interptmp(6))
 ! Interpolate along z
 call linear_interpolator_one_d(interptmp(5),interptmp(6),zd,interptmp(7))

 val = interptmp(7)

end subroutine interpolate_val

subroutine get_grid_neighbours(position,gridorigin,dx,xlower,ylower,zlower)
 ! TODO IDEALLY THIS SHOULDN'T BE HERE AND SHOULD BE IN A UTILS MODULE
 ! WITH THE VERSION USED IN METRIC_ET
 real, intent(in) :: position(3), gridorigin
 real, intent(in) :: dx
 integer, intent(out) :: xlower,ylower,zlower

 ! Get the lower grid neighbours of the position
 ! If this is broken change from floor to int
 ! How are we handling the edge case of a particle being
 ! in exactly the same position as the grid?
 ! Hopefully having different grid sizes in each direction
 ! Doesn't break the lininterp
 xlower = floor((position(1)-gridorigin)/dx)
 ylower = floor((position(2)-gridorigin)/dx)
 zlower = floor((position(3)-gridorigin)/dx)

 ! +1 because fortran
 xlower = xlower + 1
 ylower = ylower + 1
 zlower = zlower + 1

end subroutine get_grid_neighbours

logical function check_files(file1,file2,file3)
 character(len=*), intent(in) :: file1,file2,file3
 logical :: file1_exist, file2_exist, file3_exist

 inquire(file=file1,exist=file1_exist)
 inquire(file=file2,exist=file2_exist)
 inquire(file=file3,exist=file3_exist)

 check_files = file1_exist .and. file2_exist .and. file3_exist

end function check_files

end module setup
