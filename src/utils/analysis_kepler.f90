!!This code reads the data from the TDE dumpfile and stores the data in the format liked by kepler.

module analysis
  implicit none
  character(len=3), parameter, public :: analysistype = 'tde'
  public :: do_analysis

private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

   use io,              only:warning
   use dump_utils,      only:read_array_from_file
   use units,           only:udist,umass,unit_density,unit_pressure,unit_ergg,unit_velocity !units required to convert to kepler units.
   use prompting,       only:prompt
   use readwrite_dumps, only:opened_full_dump
   use vectorutils,     only:cross_product3D
   use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass,rhoh
   use centreofmass,    only:get_centreofmass
   use sortutils,       only:set_r2func_origin,indexxfunc,r2func_origin
   use eos,             only:equationofstate,entropy
   use physcon,           only:kb_on_mh
   real         :: gmw            = 2.381
   integer :: i
   integer :: no_in_bin                !this stores the number of particles in bin after each loop.
   integer :: ibin
   integer :: number_particle, ieos
   integer, parameter   :: ngrid = 512 !resolution in grid in kepler,
   integer,  intent(in) :: numfile,npart,iunit
   integer :: iorder(npart), j

   real :: pressure(ngrid),ent
   real :: rad_grid(ngrid)         !radius
   real :: mass(ngrid)
   real :: rad_vel(ngrid)          !radial velocity
   real :: density(ngrid)
   real :: temperature(ngrid)
   real :: entropy_array(ngrid)          !entropy
   real :: int_eng(ngrid)          !specific internal energy
   real :: ang_vel(ngrid)          !angular velocity
   real :: bin_mass(ngrid)         !cell mass in kepler
   real :: density_sum,density_i,gamma_inp,eni_input
   real :: u_sum,u_i               !specific internal energy storage
   real :: temperature_i,temperature_sum
   real :: grid
   real :: pressure_i,pressure_sum
   real :: Li(3),pos(3),vel(3),rad !defining angular momentum vector
   real :: xpos(3),vpos(3)         !COM position and velocity
   real :: ponrhoi,spsoundi,vel_i,vel_sum
   real :: ang_vel_sum,ang_i
   real,   intent(in) :: xyzh(:,:),vxyzu(:,:)
   real,   intent(in) :: pmass,time

   character(len=120)             :: output
   character(len=*),   intent(in) :: dumpfile

   !If dumpfile is not a complete dump we don't read it.
   if (.not.opened_full_dump) then
      write(*,'("SKIPPING FILE -- (Not a full dump)")')
      return
   endif

   !get CoM
   call get_centreofmass(xpos,vpos,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

   !use sorting algorithm to sort the particles from the center of mass as a function of radius.
   call set_r2func_origin(xpos(1),xpos(2),xpos(3))
   call indexxfunc(npart,r2func_origin,xyzh,iorder)

    ibin            = 1
    number_particle = npart/ngrid !number of particles per bin
    print*, number_particle, 'particle no'
    no_in_bin       = 0 !this keeps track of the particles added to the bin in the loop implemented.
    density_sum     = 0.
    u_sum           = 0.
    temperature_sum = 0.
    pressure_sum    = 0.
    vel_sum         = 0.
    ang_vel_sum     = 0.
    temperature_i = 0.
    gamma_inp = 5./3.

   !implementing loop for calculating the values we require.
   do j = 1, npart

     i  = iorder(j) !Access the rank of each particle in radius.

     !add 1 to no_in_bin
     no_in_bin = no_in_bin + 1

     !the position of the particle is calculated by subtracting the centre of mass.
     !xyzh is position wrt the black hole present at origin.
     pos(:) = xyzh(1:3,i) - xpos(:)

     !calculate the position which is the location of the particle.
     rad = sqrt(dot_product(pos(:),pos(:)))

     !radial velocity
     vel(:)  = vxyzu(1:3,i) - vpos(:) !relative velocity
     vel_i   = dot_product(vel(:),pos(:))/rad !this calculates (v.r)/|r| which is the dot production of velocity with position.
     vel_sum = vel_sum + vel_i

     !angular velocity
     call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),Li)
     ang_i       = sqrt(dot_product(Li, Li))
     ang_vel_sum = ang_vel_sum + ang_i

     !density
     density_i   = rhoh(xyzh(4,i),pmass)  !density of particle
     density_sum = density_sum + density_i

     !internal energy
     u_i   = vxyzu(4,i)
     u_sum = u_sum + u_i
     !print*, i, j , rad, no_in_bin, ibin, density_sum, density_i

     !using the adiabatic equation of state.
     ieos = 2
     !call eos routine
     eni_input = u_i*unit_ergg
     call equationofstate(ieos,ponrhoi,spsoundi,density_i,xyzh(1,i),xyzh(2,i),xyzh(3,i),eni=eni_input, tempi=temperature_i)

     ponrhoi = ((5./3.)-1)*u_i*unit_ergg
     !temperature_i = (1/kb_on_mh)*gmw*ponrhoi
     !pressure and temperature calculation.
     pressure_i      = ponrhoi*density_i
     pressure_sum    = pressure_sum + pressure_i
     temperature_sum = temperature_sum + temperature_i

     if (no_in_bin >= number_particle) then

       !make the bin properties.
       rad_grid(ibin)    = rad !last particle
       mass(ibin)        = j*pmass !mass of paricles < r. Calculates cell outer total mass required by kepler.
       bin_mass(ibin)    = number_particle*pmass !every bin has same mass.
       density(ibin)     = density_sum / no_in_bin
       temperature(ibin) = temperature_sum / no_in_bin
       pressure(ibin)    = pressure_sum / no_in_bin
       int_eng(ibin)     = u_sum / no_in_bin
       ang_vel(ibin)     = ang_vel_sum / no_in_bin
       rad_vel(ibin)     = vel_sum / no_in_bin

       ent = kb_on_mh * (1/gmw) * log(temperature(ibin)**1.5/(density(ibin)*unit_density))

       entropy_array(ibin)     = ent
       print*, entropy_array(ibin), 'entropy'

       !print*, 'Created bin', ibin, rad_grid(ibin) , mass(ibin), density(ibin)
       no_in_bin       = 0
       ibin            = ibin + 1
       density_sum     = 0.
       u_sum           = 0.
       temperature_sum = 0.
       pressure_sum    = 0.
       vel_sum         = 0.
       ang_vel_sum     = 0.

     end if

   end do
   ! Print the analysis being done
    write(*,'("Performing analysis type ",A)') analysistype
    write(*,'("Input file name is ",A)') dumpfile

    write(output,"(a4,i5.5)") 'ptok',numfile
    write(*,'("Output file name is ",A)') output

    !open the output file and save the data in the format kepler likes. Using same labels as kepler.
    open(iunit,file=output)
    write(iunit,'("# ",es20.12,"   # TIME")') time
    write(iunit,"('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
    !write(iunit,*) &
          1,'grid',                        &  !grid number/ bin number
          2,'mass',                        &  !bin mass
          3,'tot mass',                    &  !total mass < r
          4,'outer rad',                   &  !position
          5,'o. velocity',                 &  !velocity
          6,'density',                     &  !density
          7,'temperature',                 &  !temperature
          8,'pressure',                    &  !pressure
          9,'spec.int.energy' ,            &  !specific internal energy
          10,'spec.entropy',               &  !entropy
          11,'ang. vel.'                      !angular velocity

    do i = 1, ngrid
      grid = i
       write(iunit,'(11(es18.10,1X))') &
              grid,                            &
              bin_mass(i)*umass,               &
              mass(i)*umass,                   &
              rad_grid(i)*udist,               &
              rad_vel(i)*unit_velocity,        &
              density(i)*unit_density,         &
              temperature(i),                  &
              pressure(i)*unit_pressure,       &
              int_eng(i)*unit_ergg,            &
              entropy_array(i)*unit_ergg,            &
              ang_vel(i)
    enddo



 end subroutine do_analysis

end module analysis
