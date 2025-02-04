!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for star cluster formation calculations following
!   Bate, Bonnell & Bromm (2003). Requires pre-calculated velocity cubes.
!   This is open boundary conditions, so not valid for MHD
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
!   - relax       : *relax the cloud ?*
!
! :Dependencies: HIIRegion, centreofmass, cooling, datafiles, dim, eos,
!   infile_utils, io, kernel, mpidomain, options, part, physcon, prompting,
!   ptmass, setup_params, setvfield, spherical, subgroup, timestep, units,
!   utils_shuffleparticles, velfield
!
 use dim, only: maxvxyzu,mhd
 implicit none
 public :: setpart

 private
 integer           :: np,ieos_in
 real              :: Rsink_au,Rcloud_pc,Mcloud_msun,Temperature,mu
 real(kind=8)      :: mass_fac,dist_fac
 character(len=32) :: default_cluster
 logical           :: relax

contains

!----------------------------------------------------------------
!
!  Sets up a star cluster formation calculation
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,pc,years,kboltz,mass_proton_cgs,au,myr
 use velfield,     only:set_velfield_from_cubes
 use setup_params, only:rmax,rhozero,npart_total
 use spherical,    only:set_sphere
 use part,         only:igas,set_particle_type
 use io,           only:fatal,master,iprint
 use units,        only:umass,udist,utime,set_units
 use setvfield,    only:normalise_vfield
 use timestep,     only:dtmax,tmax
 use centreofmass, only:reset_centreofmass
 use ptmass,       only:h_acc,r_crit,rho_crit_cgs,icreate_sinks,tmax_acc,h_soft_sinkgas, &
                        r_merge_uncond,use_regnbody,f_crit_override,tseeds
 use datafiles,    only:find_phantom_datafile
 use eos,          only:ieos,gmw
 use kernel,       only:hfact_default
 use mpidomain,    only:i_belong
 use HIIRegion,    only:iH2R
 use subgroup,     only:r_neigh
 use utils_shuffleparticles, only:shuffleparticles
 use cooling,      only:Tfloor
 use options,      only:icooling

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
 integer                      :: icluster = 3 ! BBBO3 = 1, (S. Jaffa) = 2, Embedded = 3

 !--Ensure this is pure hydro
 if (mhd) call fatal('setup_cluster','This setup is not consistent with MHD.')

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
 relax       = .false.

 select case (icluster)
 case (1)
    ! from Bate, Bonnell & Bromm (2003)
    default_cluster = "Bate, Bonnell & Bromm (2003)"
    Rcloud_pc   = 0.1875  ! Input radius [pc]
    Mcloud_msun = 50.     ! Input mass [Msun]
    ieos_in     = 8       ! Barotropic equation of state
    mass_fac    = 1.0     ! mass code unit: mass_fac * solarm
    dist_fac    = 0.1     ! distance code unit: dist_fac * pc
    if (maxvxyzu >= 4) ieos_in = 2 ! Adiabatic equation of state
 case(2)
    ! Young Massive Cluster (S. Jaffa, University of Hertfordshire)
    default_cluster = "Young Massive Cluster"
    Rcloud_pc   = 5.0     ! Input radius [pc]
    Mcloud_msun = 1.0d5   ! Input mass [Msun]
    ieos_in     = 1       ! Isothermal equation of state
    mass_fac    = 1.0d5   ! mass code unit: mass_fac * solarm
    dist_fac    = 1.0     ! distance code unit: dist_fac * pc
    if (maxvxyzu >= 4) ieos_in = 2 ! Adiabatic equation of state

 case(3)
    ! Young Massive Cluster (Yann Bernard, IPAG)
    default_cluster = "Embedded cluster"
    Rcloud_pc   = 10.0    ! Input radius [pc]
    Mcloud_msun = 1.0d4   ! Input mass [Msun]
    ieos_in     = 21      ! Isothermal equation of state + HII
    mass_fac    = 1.0d4   ! mass code unit: mass_fac * solarm
    dist_fac    = 1.0     ! distance code unit: dist_fac * pc
    iH2R        = 1       ! switch HII regions
    Rsink_au    = 4000.   ! Sink radius [au]
    mu          = 2.35    ! mean molecular weight
    if (maxvxyzu >= 4) then
       ieos_in = 22 ! Adiabatic equation of state + HII
       gamma   = 5./3.
       Tfloor  = 6.
       icooling = 6
       Temperature = 40.
    endif


 case default
    ! from Bate, Bonnell & Bromm (2003)
    default_cluster = "Bate, Bonnell & Bromm (2003)"
    Rcloud_pc       = 0.1875  ! Input radius [pc]
    Mcloud_msun     = 50.     ! Input mass [Msun]
    ieos_in         = 8       ! Barotropic equation of state
    mass_fac        = 1.0     ! mass code unit: mass_fac * solarm
    dist_fac        = 0.1     ! distance code unit: dist_fac * pc
    if (maxvxyzu >= 4) ieos_in = 2 ! Adiabatic equation of state
 end select



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

 !--Set units
 call set_units(dist=dist_fac*pc,mass=mass_fac*solarm,G=1.)

 !--Define remaining variables using the inputs
 polyk         = gamma*kboltz*Temperature/(mu*mass_proton_cgs)*(utime/udist)**2
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

 if (relax) then
    call shuffleparticles(iprint,npart,xyzh,massoftype(1),rsphere=rmax,dsphere=rhozero,dmedium=0.,&
                          is_setup=.true.,prefix=trim(fileprefix))
 endif
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

 if (maxvxyzu >= 4) then
    if (gamma > 1.) then
       vxyzu(4,:) = polyk/(gamma*(gamma-1.))
    else
       vxyzu(4,:) = 1.5*polyk
    endif
 endif

 !--Setting the centre of mass of the cloud to be zero
 call reset_centreofmass(npart,xyzh,vxyzu)

 !--set options for input file, if .in file does not exist
 if (.not. inexists) then
    tmax          = 2.*t_ff
    dtmax         = 0.002*t_ff
    h_acc         = Rsink_au*au/udist
    if (icluster == 3) then
       r_crit          = h_acc
       icreate_sinks   = 2
       rho_crit_cgs    = 1.d-18
       h_soft_sinkgas  = h_acc
       tmax_acc        = 0.5*(myr/utime)
       tseeds          = 0.1*(myr/utime)
       r_merge_uncond  = h_acc
       use_regnbody    = .true.
       r_neigh         = 5e-2*h_acc
       f_crit_override = 100.
    else
       r_crit        = 2.*h_acc
       icreate_sinks = 1
       rho_crit_cgs  = 1.d-10
    endif

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
 call prompt('Do you want to relax your cloud',relax)
 if (maxvxyzu < 4) call prompt('Enter the EOS id (1: isothermal, 8: barotropic, 21: HII region expansion)',ieos_in)

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
 call write_inopt(relax, 'relax', 'relax the cloud ?', iunit)
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
 call read_inopt(relax, 'relax',db,ierr)
 call read_inopt(mu,'mu',db,ierr)
 call close_db(db)

 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_cluster: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile
!----------------------------------------------------------------
end module setup
