!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup for star cluster formation calculations following
!   Bate, Bonnell & Bromm (2003). Requires pre-calculated velocity cubes.
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters:
!   - M_cloud     : *mass of cloud in solar masses*
!   - R_cloud     : *radius of cloud in pc*
!   - Temperature : *Temperature*
!   - dist_fac    : *distance unit in pc*
!   - mass_fac    : *mass unit in Msun*
!   - mu          : *mean molecular mass*
!   - n_particles : *number of particles in sphere*
!
! :Dependencies: centreofmass, datafiles, dim, domain, eos, infile_utils,
!   io, kernel, part, physcon, prompting, ptmass, setup_params, setvfield,
!   spherical, timestep, units, velfield
!
 use dim, only: maxvxyzu
 implicit none
 public :: setpart

 private
 integer           :: np,ieos_in
 real              :: Rsink_au,Rcloud_pc,Mcloud_msun,Temperature,mu
 real(kind=8)      :: mass_fac,dist_fac
 character(len=32) :: default_cluster

contains

!----------------------------------------------------------------
!
!  Sets up a star cluster formation calculation
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,pc,years,kboltz,mass_proton_cgs,au
 use velfield,     only:set_velfield_from_cubes
 use setup_params, only:rmax,rhozero,npart_total
 use spherical,    only:set_sphere
 use part,         only:igas,set_particle_type
 use io,           only:fatal,master
 use units,        only:umass,udist,utime,set_units
 use setvfield,    only:normalise_vfield
 use timestep,     only:dtmax,tmax
 use centreofmass, only:reset_centreofmass
 use ptmass,       only:h_acc,r_crit,rho_crit_cgs,icreate_sinks
 use datafiles,    only:find_phantom_datafile
 use eos,          only:ieos,gmw
 use kernel,       only:hfact_default
 use domain,       only:i_belong
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
 real                         :: r2,totmass,epotgrav,t_ff,psep
 character(len=20), parameter :: filevx = 'cube_v1.dat'
 character(len=20), parameter :: filevy = 'cube_v2.dat'
 character(len=20), parameter :: filevz = 'cube_v3.dat'
 character(len=16)            :: lattice
 character(len=120)           :: filex,filey,filez,filein,fileset
 logical                      :: inexists,setexists
 logical                      :: BBB03 = .true. ! use the BB03 defaults, else that of a YMC (S. Jaffa)

 !--Check for existence of the .in and .setup files
 filein=trim(fileprefix)//'.in'
 inquire(file=filein,exist=inexists)
 fileset=trim(fileprefix)//'.setup'
 inquire(file=fileset,exist=setexists)

 !--Set default values
 np          = size(xyzh(1,:))
 gamma       = 1.0           ! irrelevant for ieos = 1,8
 Temperature = 10.0          ! Temperature in Kelvin (required for polyK only)
 Rsink_au    = 5.            ! Sink radius [au]
 mu          = 2.46          ! Mean molecular weight (required for polyK only)
 if (BBB03) then
    ! from Bate, Bonnell & Bromm (2003)
    default_cluster = "Bate, Bonnell & Bromm (2003)"
    Rcloud_pc   = 0.1875  ! Input radius [pc]
    Mcloud_msun = 50.     ! Input mass [Msun]
    ieos_in     = 8       ! Barotropic equation of state
    mass_fac    = 1.0     ! mass code unit: mass_fac * solarm
    dist_fac    = 0.1     ! distance code unit: dist_fac * pc
 else
    ! Young Massive Cluster (S. Jaffa, University of Hertfordshire)
    default_cluster = "Young Massive Cluster"
    Rcloud_pc   = 5.0     ! Input radius [pc]
    Mcloud_msun = 1.0d5   ! Input mass [Msun]
    ieos_in     = 1       ! Isothermal equation of state
    mass_fac    = 1.0d5   ! mass code unit: mass_fac * solarm
    dist_fac    = 1.0     ! distance code unit: dist_fac * pc
 endif

 !--Set units
 call set_units(dist=dist_fac*pc,mass=mass_fac*solarm,G=1.)

 if (maxvxyzu >= 4) ieos_in = 2 ! Adiabatic equation of state

 !--Read values from .setup
 if (setexists) then
    call read_setupfile(fileset,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(fileset)
       stop
    endif
    !--Prompt to get inputs and write to file
 elseif (id==master) then
    print "(a,/)",trim(fileset)//' not found: using interactive setup'
    call get_input_from_prompts()
    call write_setupfile(fileset)
 endif

 !--Define remaining variables using the inputs
 polyk         = kboltz*Temperature/(mu*mass_proton_cgs)*(utime/udist)**2
 rmax          = Rcloud_pc*(pc/udist)
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

 !--Set velocities (from pre-made velocity cubes)
 write(*,"(1x,a)") 'Setting up velocity field on the particles...'
 vxyzu(:,:) = 0.
 filex = find_phantom_datafile(filevx,'velfield')
 filey = find_phantom_datafile(filevy,'velfield')
 filez = find_phantom_datafile(filevz,'velfield')

 call set_velfield_from_cubes(xyzh,vxyzu,npartoftype(igas),filex,filey,filez,1.,rmax,.false.,ierr)
 if (ierr /= 0) call fatal('setup','error setting up velocity field')

 !--Normalise the energy
 call normalise_vfield(npart,vxyzu,ierr,ke=epotgrav)
 if (ierr /= 0) call fatal('setup','error normalising velocity field')

 !--Setting the centre of mass of the cloud to be zero
 call reset_centreofmass(npart,xyzh,vxyzu)

 !--set options for input file, if .in file does not exist
 if (.not. inexists) then
    tmax          = 2.*t_ff
    dtmax         = 0.002*t_ff
    h_acc         = Rsink_au*au/udist
    r_crit        = 2.*h_acc
    icreate_sinks = 1
    rho_crit_cgs  = 1.d-10
    ieos          = ieos_in
    gmw           = mu       ! for consistency; gmw will never actually be used
 endif

 !--Print summary
 if (id==master) then
    write(*,"(1x,a)") 'Cluster formation setup: '
    write(*,"(1x,a,es10.3,a)") '       Rcloud = ',Rcloud_pc,' pc'
    write(*,"(1x,a,es10.3,a)") '       Mcloud = ',Mcloud_msun,' Msun'
    write(*,"(1x,a,es10.3,a)") ' Mean density = ',rhozero*umass/udist**3,' g/cm^3'
    write(*,"(1x,a,es10.3,a,es10.3,a)") 'Particle mass = ',massoftype(1)*(umass/solarm),' Msun'
    write(*,"(1x,a,es10.3,a,e10.3,a)") 'Freefall time = ',t_ff*(utime/years),' years (',t_ff,' in code units)'
    write(*,"(1x,a,es10.3,a)") '  Sound speed = ',sqrt(polyk)*(udist/utime), ' cm/s'
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
 call prompt('Enter the distance unit (in pc)',dist_fac)
 call prompt('Enter the mass unit (in Msun)',mass_fac)
 call prompt('Enter the mass of the cloud (in Msun)',Mcloud_msun)
 call prompt('Enter the radius of the cloud (in pc)',Rcloud_pc)
 call prompt('Enter the radius of the sink particles (in au)',Rsink_au)
 call prompt('Enter the Temperature of the cloud (used for initial sound speed)',Temperature)
 call prompt('Enter the mean molecular mass (used for initial sound speed)',mu)
 if (maxvxyzu < 4) call prompt('Enter the EOS id (1: isothermal, 8: barotropic)',ieos_in)

end subroutine get_input_from_prompts
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
 write(iunit,"(a)") '# input file for cluster setup routines'
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'n_particles','number of particles in sphere',iunit)
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_fac,'dist_fac','distance unit in pc',iunit)
 call write_inopt(mass_fac,'mass_fac','mass unit in Msun',iunit)
 write(iunit,"(/,a)") '# options for sphere'
 call write_inopt(Mcloud_msun,'M_cloud','mass of cloud in solar masses',iunit)
 call write_inopt(Rcloud_pc,'R_cloud','radius of cloud in pc',iunit)
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
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
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
 call read_inopt(dist_fac,'dist_fac',db,ierr)
 call read_inopt(mass_fac,'mass_fac',db,ierr)
 call read_inopt(Mcloud_msun,'M_cloud',db,ierr)
 call read_inopt(Rcloud_pc,'R_cloud',db,ierr)
 call read_inopt(Temperature,'Temperature',db,ierr)
 call read_inopt(mu,'mu',db,ierr)
 call close_db(db)

 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_cluster: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile
!----------------------------------------------------------------
end module setup
