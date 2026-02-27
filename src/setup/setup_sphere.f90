!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for collapsing uniform sphere
!   This is open boundary conditions, so not valid for MHD
!
! :References: None
!
! :Owner: Alison Young
!
! :Runtime parameters:
!   - M_cloud     : *mass of cloud in solar masses*
!   - R_cloud     : *radius of cloud in au*
!   - Temperature : *Temperature*
!   - mu          : *mean molecular mass*
!   - n_particles : *number of particles in sphere*
!
! :Dependencies: centreofmass, cooling, datafiles, dim, dynamic_dtmax, eos,
!   eos_stamatellos, infile_utils, io, io_control, kernel, mpidomain,
!   options, part, physcon, prompting, setup_params, spherical, timestep,
!   units
!
 use dim, only:maxvxyzu,mhd
 implicit none
 public :: setpart

 private
 integer           :: np,ieos_in
 real              :: Rcloud,Mcloud_msun,Temperature,mu
 character(len=32) :: default_cluster

contains

!----------------------------------------------------------------
!
!  Sets up a collapsing sphere calculation
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,years,mass_proton_cgs,au,myr,kb_on_mh
 use setup_params, only:rmax,rhozero,npart_total
 use spherical,    only:set_sphere
 use part,         only:igas,set_particle_type
 use io,           only:fatal,master
 use units,        only:umass,udist,utime,set_units,unit_ergg,unit_density
 use timestep,     only:dtmax,tmax
 use dynamic_dtmax,only:dtmax_dratio,dtmax_min
 use centreofmass, only:reset_centreofmass
 use datafiles,    only:find_phantom_datafile
 use eos,          only:ieos,gmw
 use kernel,       only:hfact_default
 use mpidomain,    only:i_belong
 use cooling,      only:Tfloor
 use options,      only:icooling
 use io_control,   only:rhofinal_cgs
 use infile_utils, only:get_options,infile_exists
 use eos_stamatellos, only:getintenerg_opdep,read_optab,eos_file
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer                      :: i,ierr
 real                         :: r2,totmass,epotgrav,t_ff,psep,uinit
 character(len=16)            :: lattice

 !--ensure this is pure hydro
 if (mhd) call fatal('setup_cluster','This setup is not consistent with MHD.')

 !--Set default values
 np          = 10000
 gamma       = 1.4           ! irrelevant for ieos = 1,8
 Temperature = 10.0          ! Temperature in Kelvin (required for polyK only)
 mu          = 2.46          ! Mean molecular weight (required for polyK only)
 polyk       = 0.

 !single protostellar core with ieos==24
 default_cluster = "Single protostellar core"
 Rcloud   = 10000. ! Input radius [au]
 Mcloud_msun = 1.   ! Input mass [Msun]
 ieos_in     = 24       ! Isothermal equation of state
 icooling    = 9
 Tfloor        = 5.

 !--Read values from .setup
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,get_input_from_prompts)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 !--Set units
 call set_units(dist=au,mass=solarm,G=1.)

 !--Define remaining variables using the inputs
 rmax          = Rcloud
 r2            = rmax*rmax
 totmass       = Mcloud_msun*(solarm/umass)
 rhozero       = totmass/(4./3.*pi*rmax*r2)
 hfact         = hfact_default
 t_ff          = sqrt(3.*pi/(32.*rhozero))  ! free-fall time (the characteristic timescale)
 epotgrav      = 3./5.*totmass**2/rmax      ! Gravitational potential energy
 lattice       = 'random'

 !--Set positions
 call set_sphere(trim(lattice),id,master,0.,rmax,psep,hfact,npart,xyzh,nptot=npart_total, &
                 exactN=.true.,np_requested=np,mask=i_belong)
 npartoftype(:) = 0
 npartoftype(1) = npart
 massoftype(1)  = totmass/real(npart)
 !--Set particle properties
 do i = 1,npart
    call set_particle_type(i,igas)
 enddo

 !--Set velocities (to zero)
 write(*,"(1x,a)") 'Setting up velocity field on the particles...'
 vxyzu(:,:) = 0.

 if (maxvxyzu >= 4) then
    if (ieos_in ==24) then
       call read_optab(eos_file,ierr)
       if (ierr/=0) stop
       call getintenerg_opdep(Temperature,rhozero*unit_density,uinit)
       vxyzu(4,:) = uinit/unit_ergg
       print *, "Temperature:", Temperature, "K, u=",uinit,"erg/g"
    endif
 endif

 !--Setting the centre of mass of the cloud to be zero
 call reset_centreofmass(npart,xyzh,vxyzu)

 !--set options for input file, if .in file does not exist
 if (.not. infile_exists(fileprefix)) then
    tmax          = 2.*t_ff
    dtmax         = 0.01*t_ff
    dtmax_min     = 0.0
    dtmax_dratio  = 1.258
    rhofinal_cgs  = 0.1
    ieos          = ieos_in
    gmw           = mu       ! for consistency; gmw will never actually be used
 endif

 !--Print summary
 if (id==master) then
    write(*,"(1x,a)") 'Collapsing sphere setup: '
    write(*,"(1x,a,es10.3,a)") '       Rcloud = ',Rcloud,' au'
    write(*,"(1x,a,es10.3,a)") '       Mcloud = ',Mcloud_msun,' Msun'
    write(*,"(1x,a,es10.3,a)") ' Mean density = ',rhozero*umass/udist**3,' g/cm^3'
    write(*,"(1x,a,es10.3,a,es10.3,a)") 'Particle mass = ',massoftype(1)*(umass/solarm),' Msun'
    write(*,"(1x,a,es10.3,a,e10.3,a)") 'Freefall time = ',t_ff*(utime/years),' years (',t_ff,' in code units)'
    write(*,"(1x,a,es10.3,a)") '  Sound speed = ', sqrt(gamma*kb_on_mh*temperature/mu), ' cm/s'
 endif

end subroutine setpart
!----------------------------------------------------------------
!
!  Prompt user for inputs
!
!----------------------------------------------------------------
subroutine get_input_from_prompts()
 use prompting, only:prompt

 write(*,'(2a)') 'Default settings: ',trim(default_cluster)
 call prompt('Enter the number of particles in the sphere',np,0,np)
 call prompt('Enter the mass of the cloud (in Msun)',Mcloud_msun)
 call prompt('Enter the radius of the cloud (in au)',Rcloud)
 call prompt('Enter the Temperature of the cloud (used for initial sound speed)',Temperature)
 call prompt('Enter the mean molecular mass (used for initial sound speed)',mu)
 if (maxvxyzu < 4) call prompt('Enter the EOS id (1: isothermal, 8: barotropic, 21: HII region expansion)',ieos_in)

end subroutine get_input_from_prompts
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
 write(iunit,"(a)") '# input file for sphere setup routines'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'n_particles','number of particles in sphere',iunit)
 write(iunit,"(/,a)") '# options for sphere'
 call write_inopt(Mcloud_msun,'M_cloud','mass of cloud in solar masses',iunit)
 call write_inopt(Rcloud,'R_cloud','radius of cloud in au',iunit)
 write(iunit,"(/,a)") '# options required for initial sound speed'
 call write_inopt(Temperature,'Temperature','Temperature',iunit)
 call write_inopt(mu,'mu','mean molecular mass',iunit)
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(np,'n_particles',db,ierr)
 call read_inopt(Mcloud_msun,'M_cloud',db,ierr)
 call read_inopt(Rcloud,'R_cloud',db,ierr)
 call read_inopt(Temperature,'Temperature',db,ierr)
 call read_inopt(mu,'mu',db,ierr)
 call close_db(db)

 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_cluster: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile
!----------------------------------------------------------------
end module setup
