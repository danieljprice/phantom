!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!

module analysis
  !
  ! Module for generating KEPLER file from a TDE dumpfile.
  !
  ! :References: None
  !
  !
  ! :Dependencies:dump_utils,units,io,prompting,readwrite_dumps,vectorutils,
  !               part,centreofmass,sortutils,eos,physcon
  !
  implicit none
  character(len=3), parameter, public :: analysistype = 'tde'
  public :: do_analysis

private

contains
  !----------------------------------------------------------------
  !+
  !  routine to write an input file for KEPLER.
  !  uses phantom_to_kepler_arrays subroutine.
  !+
  !----------------------------------------------------------------
subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

   use io,              only:warning
   use dump_utils,      only:read_array_from_file
   use units,           only:udist,umass,unit_density,unit_pressure,unit_ergg,unit_velocity !units required to convert to kepler units.
   use prompting,       only:prompt
   use readwrite_dumps, only:opened_full_dump
   use part,            only: abundance

   integer, parameter   :: ngrid = 512 !resolution in grid in kepler
   integer,  intent(in) :: numfile,npart,iunit
   real :: grid
   integer :: i

   real,   intent(in) :: xyzh(:,:),vxyzu(:,:)
   real,   intent(in) :: pmass,time
   real , dimension(ngrid)::pressure,rad_grid,mass,rad_vel,density,temperature,entropy_array,&
                            int_eng,ang_vel,bin_mass,y_e,a_bar
  real, dimension(3,ngrid) :: angular_vel_3D

   character(len=120)             :: output
   character(len=*),   intent(in) :: dumpfile

   !If dumpfile is not a complete dump we don't read it.
   if (.not.opened_full_dump) then
      write(*,'("SKIPPING FILE -- (Not a full dump)")')
      return
   endif

   !if dumpfile is a full dump, we call the subroutine for getting the arrays we need
   call phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,ngrid,pressure,rad_grid,mass,rad_vel,&
                                  density,temperature,entropy_array,int_eng,ang_vel,bin_mass,&
                                  y_e,a_bar,angular_vel_3D)

    !Print the analysis being done
    write(*,'("Performing analysis type ",A)') analysistype
    write(*,'("Input file name is ",A)') dumpfile

    write(output,"(a4,i5.5)") 'ptok',numfile
    write(*,'("Output file name is ",A)') output

    !open the output file and save the data in the format kepler likes. Using same labels as kepler.
    open(iunit,file=output)
    write(iunit,'("# ",es20.12,"   # TIME")') time
    write(iunit,"('#',13(1x,'[',i2.2,1x,a28,']',2x))") &
          1,'grid',                                    &  !grid number/ bin number
          2,'cell mass',                               &  !bin mass
          3,'cell outer total mass',                   &  !total mass < r
          4,'cell outer radius',                       &  !position
          5,'cell outer velocity',                     &  !velocity
          6,'cell density',                            &  !density
          7,'cell temperature',                        &  !temperature
          8,'cell pressure',                           &  !pressure
          9,'cell spec. int. energy' ,                 &  !specific internal energy
          10,'cell specific entropy',                  &  !entropy
          11,'cell angular velocity',                  &  !angular velocity
          12,'cell A_bar',                             &  !average molecular mass
          13,'cell Y_e'

    do i = 1, ngrid
      grid = i
       write(iunit,'(13(es18.10,1X))')         &
              grid,                            &
              bin_mass(i)*umass,               &
              mass(i)*umass,                   &
              rad_grid(i)*udist,               &
              rad_vel(i)*unit_velocity,        &
              density(i)*unit_density,         &
              temperature(i),                  &
              pressure(i)*unit_pressure,       &
              int_eng(i)*unit_ergg,            &
              entropy_array(i),                &
              angular_vel_3D(3,i),             &
              a_bar(i),                        &
              y_e(i)
    enddo

    call find_mass_lost(xyzh,vxyzu,pmass,npart)
    !print*,mass(ngrid)*umass,'total particle mass',pmass,'pmass',umass,'umass'
 end subroutine do_analysis

 !----------------------------------------------------------------
 !+
 !  routine for binning the data as a function of radius.
 !  The arrays generated are used by do_analysis subroutine.
 !+
 !----------------------------------------------------------------
 subroutine phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,ngrid,pressure,rad_grid,mass,&
                                    rad_vel,density,temperature,entropy_array,int_eng,ang_vel,bin_mass,&
                                    y_e,a_bar,angular_vel_3D)

   use units,           only:unit_ergg,udist,umass,unit_density,unit_pressure!units required to convert to kepler units.
   use vectorutils,     only:cross_product3D
   use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass,rhoh
   use centreofmass,    only:get_centreofmass
   use sortutils,       only:set_r2func_origin,indexxfunc,r2func_origin
   use eos,             only:equationofstate,entropy,X_in,Z_in,entropy,gmw,init_eos
   use physcon,         only:kb_on_mh,kboltz,atomic_mass_unit,avogadro

   integer :: i
   integer,  intent(in) :: npart,ngrid
   integer :: no_in_bin  !this stores the number of particles in bin after each loop.
   integer :: ibin
   integer :: iorder(npart), j
   integer :: number_particle, ieos,ierr

   real,intent(in) :: xyzh(:,:),vxyzu(:,:)
   real,intent(in) :: pmass,time
   real,intent(out):: rad_grid(ngrid),mass(ngrid),rad_vel(ngrid),density(ngrid)!rad_grid stores radius, rad_vel stores radial velocity
   real,intent(out):: temperature(ngrid),entropy_array(ngrid),int_eng(ngrid),ang_vel(ngrid)
   real,intent(out):: pressure(ngrid),y_e(ngrid),a_bar(ngrid),angular_vel_3D(3,ngrid)

   real :: bin_mass(ngrid)         !cell mass
   real :: density_sum,density_i,eni_input
   real :: u_sum,u_i,omega_sum(3)               !specific internal energy storage
   real :: temperature_i,temperature_sum
   real :: pressure_i,pressure_sum
   real :: Li(3),pos(3),vel(3),rad !defining angular momentum vector
   real :: xpos(3),vpos(3),omega(3) !COM position and velocity
   real :: ponrhoi,spsoundi,vel_i,vel_sum
   real :: ang_vel_sum,ang_i,number_density,moment_of_inertia
   real :: min_pos(3),magn_rad
   real :: ent1, ent2,ent,Y_in

   ! we use the equation number 12 from eos file.
   ieos = 12
   call init_eos(ieos,ierr)
   !The star is not on the origin as BH exists at that point. Hence we need to get COM.
   call get_centreofmass(xpos,vpos,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
   print*, '(x,y,z)','->', '(',xpos(1),xpos(2),xpos(3),')','COM in code units'
   print*, '(x,y,z)','->', '(',xpos(1)*udist,xpos(2)*udist,xpos(3)*udist,')','COM in cm'
   !use sorting algorithm to sort the particles from the center of mass as a function of radius.
   call set_r2func_origin(xpos(1),xpos(2),xpos(3))
   call indexxfunc(npart,r2func_origin,xyzh,iorder)

    ibin            = 1
    number_particle = npart/ngrid !number of particles per bin
    print*, number_particle, 'particle no',npart,'npart',ngrid,'ngrid'
    no_in_bin       = 0 !this keeps track of the particles added to the bin in the loop implemented.
    density_sum     = 0.
    u_sum           = 0.
    temperature_sum = 0.
    pressure_sum    = 0.
    vel_sum         = 0.
    ang_vel_sum     = 0.
    temperature_i   = 0.
    omega_sum(:)        = 0.
    moment_of_inertia = 0.
    X_in = 0.71492308
    Y_in = 0.27112283
    Z_in = 1.-X_in-Y_in
    gmw = 0.61 !mean molecular weight
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
     vel_i   = dot_product(vel(:),pos(:))/rad
     !this calculates (v.r)/|r| which is the dot production of velocity with position.
     vel_sum = vel_sum + vel_i

     !angular velocity
     !it is calculated by taking cross product of position and velocity.
     !Then the magntitude of this vector is taken to get the ang_i.
     call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),Li)
     ang_i       = sqrt(dot_product(Li, Li))
     ang_vel_sum = ang_vel_sum + ang_i
     moment_of_inertia = rad**2 + moment_of_inertia
     omega(:)  = Li(:)
     omega_sum(:) = omega_sum(:) + omega(:)

     !density
     density_i   = rhoh(xyzh(4,i),pmass)  !density of particle
     density_sum = density_sum + density_i

     !internal energy
     u_i   = vxyzu(4,i)
     u_sum = u_sum + u_i

     !using the adiabatic equation of state.

     eni_input = u_i
     !call eos routine
     call equationofstate(ieos,ponrhoi,spsoundi,density_i,xyzh(1,i),xyzh(2,i),xyzh(3,i),eni=eni_input, tempi=temperature_i)

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
       angular_vel_3D(:,ibin)    = omega_sum(:) / moment_of_inertia

       !calculating Y_e = X_e /(A_e*m_u*N_A)
       y_e(ibin)         = (X_in/(1.*avogadro*atomic_mass_unit)) + (Y_in/(4.*avogadro*atomic_mass_unit))
       a_bar(ibin)       = X_in+(4.*Y_in) !average atomic mass in each bin.
       !calculating entropy
       !because we are using adibatic eos, the change in entropy is 0, as, dQ = 0.
       !Using the eos from eos module. I think this is a monoatomic gas, Sackur-Tetrode equation.
       !ent                 = kb_on_mh * (1/gmw) * log(temperature(ibin)**1.5/(density(ibin)))
       !entropy_array(ibin) = ent
       number_density = density(ibin)/(gmw*atomic_mass_unit)
       !implementing entropy from the Sackur-Tetrode equation found on https://scholar.harvard.edu/files/schwartz/files/6-entropy.pdf
       !ent1 = number_density*kboltz*(log((bin_mass(ibin)*umass*pmass)/(density(ibin)*unit_density)))
       !ent2 = number_density*kboltz*(1.5*log(1.5*kboltz*temperature(ibin)))
       !ent = ent1 + ent2


       entropy_array(ibin) = entropy(density(ibin)*unit_density,pressure(ibin)*unit_pressure,2,ierr)
       entropy_array(ibin) = entropy_array(ibin)/(kboltz*avogadro)
       if (ierr/=0) then
         print*, 'Entropy is calculated incorrectly'
       end if
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

 end subroutine phantom_to_kepler_arrays

 !----------------------------------------------------------------
 !+
 !  routine for finding the number of partilces that have been
 !  stripped from the star we started with. This can be used to
 !  calculate the mass lost by the star.
 !  Might be useful in correctly finding the COM when the star has
 !  lost mass.
 !+
 !----------------------------------------------------------------
 subroutine find_mass_lost(xyzh,vxyzu,pmass,npart)
   use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass
   use centreofmass,    only:get_centreofmass
   use sortutils,       only:set_r2func_origin,indexxfunc,r2func_origin
   use units ,          only:umass,unit_velocity,udist
   use physcon,         only:gg,solarm
   integer,  intent(in) :: npart
   real,   intent(in) :: xyzh(:,:),vxyzu(:,:)
   real,   intent(in) :: pmass
   real :: xpos(3),vpos(3),pos(3),rad,percentage,count
   real :: v2,eps(npart)

   integer :: i, iorder(npart),j
   real :: mh

   !sort the particles based on radius.
   !compare the
   !get CoM
   call get_centreofmass(xpos,vpos,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

   !use sorting algorithm to sort the particles from the center of mass as a function of radius.
   call set_r2func_origin(xpos(1),xpos(2),xpos(3))
   call indexxfunc(npart,r2func_origin,xyzh,iorder)

   do j = 1,npart
     !the position of the particle is calculated by subtracting the centre of mass.
     !xyzh is position wrt the black hole present at origin.
     i  = iorder(j) !Access the rank of each particle in radius.
     pos(:) = xyzh(1:3,i) - xpos(:)

    !calculate the position which is the location of the particle.
    rad = sqrt(dot_product(pos(:),pos(:)))
    v2     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))

    !Energy = KE + PE = v^2/2 - G mh/r.
    eps(i) = (v2**2*unit_velocity**2)/2. - (1e6*solarm*gg)/(rad*udist) !in code units

  end do

 count = 0.

   do i = 1,npart
     !if energy is positive, then particle has been removed from star.
     if (eps(i) > 0) then
       count = count + 1.
     end if
   end do

   print*,count,'count',npart,'npart'
   print*, count*umass*pmass, 'mass lost', npart*umass*pmass,umass,'mass unit'
   percentage = (count/npart)*100.
   print*, percentage,'mass lost percentage'
  end subroutine find_mass_lost

end module analysis
