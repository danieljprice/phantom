!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! analysis
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, dump_utils, eos, fileutils, io, part,
!   physcon, prompting, readwrite_dumps, sortutils, units, vectorutils
!

 !
 ! Module for generating KEPLER file from a TDE dumpfile.
 !
 ! :References: None
 !
 !
 ! :Dependencies:dump_utils,units,io,prompting,readwrite_dumps,vectorutils,
 !               part,centreofmass,sortutils,eos,physcon,fileutils
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

 use io,              only : warning
 use dump_utils,      only : read_array_from_file
 use units,           only : udist,umass,unit_density,unit_pressure,unit_ergg,unit_velocity !units required to convert to kepler units.
 use prompting,       only : prompt
 use readwrite_dumps, only : opened_full_dump

 integer,  intent(in) :: numfile,npart,iunit
 integer              :: i,j,n_comp
 integer              :: ngrid = 0

 real                              :: grid
 real,intent(in)                   :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)                   :: pmass,time
 real , allocatable,dimension(:)   :: pressure,rad_grid,mass,density,temperature,entropy_array,&
                                        int_eng,bin_mass,y_e,a_bar,rad_mom
 real, allocatable                 :: velocity_3D(:,:),angular_vel_3D(:,:)
 real, allocatable                 :: composition_kepler(:,:)

 character(len=20),allocatable     :: comp_label(:)
 character(len=120)                :: output
 character(len=*),intent(in)       :: dumpfile
 integer :: max_pos

 !If dumpfile is not a complete dump we don't read it.
 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif

 !if dumpfile is a full dump, we call the subroutine for getting the arrays we need
 call phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,pressure,rad_grid,mass,&
                                  density,temperature,entropy_array,int_eng,velocity_3D,bin_mass,&
                                  y_e,a_bar,composition_kepler,comp_label,n_comp,ngrid,rad_mom,&
                                  angular_vel_3D)

 !allocate for composition_kepler
 !Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'ptok',numfile
 write(*,'("Output file name is ",A)') output

 !open the output file and save the data in the format kepler likes. Using same labels as kepler.
 open(iunit,file=output)
 write(iunit,'("# ",es20.12,"   # Version")') time
 !  write(iunit,"('#',50(1x,'[',1x,a22,']',2x))")    &
 write(iunit,"('#',50(a22,1x))")    &
          'grid',                                    &  !grid number/ bin number
          'cell mass',                               &  !bin mass
          'cell outer tot. mass',                    &  !total mass < r
          'cell outer radius',                       &  !position
          'cell density',                            &  !density
          'cell temperature',                        &  !temperature
          'cell radial momentum',                    &  !radial momentum
          'angular vel (x)',                         &  !velocity x component
          'angular vel (y)',                         &  !velocity y component
          'angular vel (z)',                         &  !velocity z component
          comp_label                                    !chemical composition

 print*, shape(composition_kepler),'kepler compo'
 do i = 1, ngrid
    grid = i
    write(iunit,'(50(es18.10,1X))')                     &
              grid,                                       &
              bin_mass(i)*umass,                          &
              mass(i)*umass,                              &
              rad_grid(i)*udist,                          &
              density(i)*unit_density,                    &
              temperature(i),                             &
              rad_mom(i)*unit_velocity*umass,             &
              (angular_vel_3D(j,i)*unit_velocity, j=1,3), &
              (composition_kepler(j,i), j=1,n_comp)
 enddo

end subroutine do_analysis

 !----------------------------------------------------------------
 !+
 !  routine for binning the data as a function of radius.
 !  The arrays generated are used by do_analysis subroutine.
 !+
 !----------------------------------------------------------------
subroutine phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,pressure,rad_grid,mass,&
                                    density,temperature,entropy_array,int_eng,velocity_3D,bin_mass,&
                                    y_e,a_bar,composition_kepler,comp_label,columns_compo,ngrid,rad_mom,&
                                    angular_vel_3D)

 use units,           only : udist,unit_density,unit_pressure!units required to convert to kepler units.
 use vectorutils,     only : cross_product3D
 use part,            only : nptmass,xyzmh_ptmass,vxyz_ptmass,rhoh
 use centreofmass,    only : get_centreofmass
 use sortutils,       only : set_r2func_origin,indexxfunc,r2func_origin
 use eos,             only : equationofstate,entropy,X_in,Z_in,entropy,gmw,init_eos
 use physcon,         only : kb_on_mh,kboltz,atomic_mass_unit,avogadro

 integer,intent(in)   :: npart
 integer,intent(out)  :: ngrid
 real,intent(in)      :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)      :: pmass,time
 real,intent(out),allocatable    :: rad_grid(:),mass(:),density(:)!rad_grid stores radius
 real,intent(out),allocatable     :: temperature(:),entropy_array(:),int_eng(:),bin_mass(:),rad_mom(:)
 real,intent(out),allocatable    :: pressure(:),y_e(:),a_bar(:),velocity_3D(:,:)
 real,intent(out),allocatable              :: composition_kepler(:,:),angular_vel_3D(:,:)
 character(len=20),allocatable,intent(out) :: comp_label(:)

 integer :: no_in_bin !this stores the number of particles in bin after each loop.
 integer :: ibin
 integer :: iorder(npart),j,i,s,m,index_val
 integer :: number_particle,ieos,ierr
 integer :: columns_compo,location
 integer :: particle_sum
 integer::  number_bins,number_tot, number_per_bin
 real :: p_no
 real :: density_sum,density_i,eni_input
 real :: u_sum,u_i !specific internal energy storage
 real :: omega_sum(3),omega(3)
 real :: moment_of_inertia
 real :: temperature_i,temperature_sum
 real :: rad_velocity,rad_vel_sum,momentum
 real :: pressure_i,pressure_sum
 real :: pos(3),vel(3),rad,rad_next
 real :: xpos(3),vpos(3),star_centre(3) !COM position and velocity
 real :: ponrhoi,spsoundi,vel_sum(3),Li(3)
 real :: Y_in
 real :: mu_i
 real,allocatable    :: interpolate_comp(:,:),composition_i(:),composition_sum(:)
 real, dimension(17) :: z_value,a_value

 !The star is not on the origin as BH exists at that point.
 !minimum h value corresponds to position of maximum density.
 !COM is not a good option as it does not work for severe disruptione events.
 location = minloc(xyzh(4,:),dim=1)
 star_centre(:) = xyzh(1:3,location)
 print*,'density at center',rhoh(xyzh(4,location),pmass),xyzh(4,location)
 print*, xyzh(1:3,location),'center of star in code units'
 !we use the equation number 12 from eos file.
 ieos = 12
 call init_eos(ieos,ierr)

 !use sorting algorithm to sort the particles from the center of star as a function of radius.
 xpos(:) = star_centre(:)
 vpos(:) = vxyzu(1:3,location)
 call set_r2func_origin(xpos(1),xpos(2),xpos(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
 print*, xpos(:),'xpos',vpos(:),'vpos'
 !Call composition_array subroutine to get the composition.
 call composition_array(interpolate_comp,columns_compo,comp_label) !rows correspond to the particle used and column correspond to the elements.

 if (columns_compo /= 0) then
    if (size(interpolate_comp(1,:))/= npart) then
       print*, 'Error: Number of particles does not match with the composition size'
       return
    endif
 endif

 ibin = 1
 !instead of setting the bin number, the goal is to fix particle per bin for each run and calculate
 !the ngrid using it.
 number_particle = 201
 number_per_bin = number_particle
 number_tot = npart
 !number_particle    = (npart/ngrid)!number of particles per bin
 number_bins = number_tot/number_per_bin
 ngrid = npart/number_particle

 print*, mod(number_tot,number_per_bin),'mod(number_tot,number_per_bin)'
 if (mod(number_tot,number_per_bin) > 0 ) then
    ngrid = ngrid + 1
 endif
 print*, ngrid, 'ngrid'

 allocate(rad_grid(ngrid),mass(ngrid),density(ngrid))!rad_grid stores radius, stores radial velocity
 allocate(temperature(ngrid),entropy_array(ngrid),int_eng(ngrid),bin_mass(ngrid),rad_mom(ngrid))
 allocate(pressure(ngrid),y_e(ngrid),a_bar(ngrid),velocity_3D(3,ngrid),angular_vel_3D(3,ngrid))
 no_in_bin          = 0 !this keeps track of the particles added to the bin in the loop implemented.
 density_sum        = 0.
 rad_vel_sum        = 0.
 u_sum              = 0.
 temperature_sum    = 0.
 pressure_sum       = 0.
 vel_sum            = 0.
 temperature_i      = 0.
 moment_of_inertia  = 0.
 omega_sum(:)       = 0.
 X_in               = 0.71492308
 Y_in               = 0.27112283
 Z_in               = 1.-X_in-Y_in
 gmw                = 0.61 !mean molecular weight
 z_value = (/1,2,2,6,7,8,10,12,14,16,18,20,22,24,26,26,28/)
 a_value = (/1,3,4,12,14,16,20,24,28,32,36,40,44,48,52,56,56/)
 !allocating storage for composition of one particle.
 allocate(composition_i(columns_compo))
 allocate(composition_sum(columns_compo))
 allocate(composition_kepler(columns_compo,ngrid))
 composition_sum(:) = 0.
 !implementing loop for calculating the values we require.
 do j = 1, npart

    i  = iorder(j) !Access the rank of each particle in radius.

    !add 1 to no_in_bin
    no_in_bin = no_in_bin + 1

    !the position of the particle is calculated by subtracting the point of highest density.
    !xyzh is position wrt the black hole present at origin.
    pos(:) = xyzh(1:3,i) - xpos(:)
    !calculate the position which is the location of the particle.
    rad    = sqrt(dot_product(pos(:),pos(:)))

    !veloctiy
    vel(:)     = vxyzu(1:3,i) - vpos(:)
    vel_sum(:) = vel_sum(:) + vel(:)
    if (j==1) then
       print*, vel(:),pos(:)
    endif
    !radial momentum
    !we skip the first particle as its the one that exists at the center of star and hence will give infinite rad_vel as rad = 0.
    if (rad > 0.) then
       rad_velocity = dot_product(vel(:),pos(:))/rad
       momentum     = rad_velocity*pmass
       rad_vel_sum  = rad_vel_sum + momentum
    endif

    !density
    density_i   = rhoh(xyzh(4,i),pmass)
    density_sum = density_sum + density_i

    !internal energy
    u_i   = vxyzu(4,i)
    u_sum = u_sum + u_i

    !composition
    if (columns_compo /= 0) then
       composition_i(:)   = interpolate_comp(:,i)
    endif
    composition_sum(:) = composition_sum(:) + composition_i(:)

    !calculate mean molecular weight that is required by the eos module using the
    !mass fractions for each particle.
    !do not consider neutron which is the first element of the composition_i array.
    mu_i = 0
    if (columns_compo /= 0) then
       do index_val = 1, columns_compo-1
          mu_i = mu_i + (composition_i(index_val+1)*(1+z_value(index_val)))/a_value(index_val)
       enddo
       gmw = 1/mu_i
    endif

    !using the adiabatic equation of state.
    eni_input = u_i
    !call eos routine
    call equationofstate(ieos,ponrhoi,spsoundi,density_i,xyzh(1,i),xyzh(2,i),xyzh(3,i),eni=eni_input, tempi=temperature_i)

    pressure_i      = ponrhoi*density_i
    pressure_sum    = pressure_sum + pressure_i
    temperature_sum = temperature_sum + temperature_i

    !angular velocity of the cells. In kepler, we need it for the interface.
    call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),Li)
    moment_of_inertia = moment_of_inertia + rad**2*pmass*(2./3.)  !moment of interia of shell
    omega(:)          = Li(:)*pmass
    omega_sum(:)      = omega_sum(:) + omega(:)

    if (no_in_bin >= number_particle .or. ibin == ngrid)  then
       !make the bin properties.
       !the radius is calculated as the mean of i and i+1 bins,
       if (ibin /= ngrid) then
          s = j + 1 !we consider the first particle in the i+1 bin and then use it.
          m = iorder(s)
          pos(:) = xyzh(1:3,m) - xpos(:)
          !calculate the position which is the location of the particle.
          rad_next       = sqrt(dot_product(pos(:),pos(:)))
          rad_grid(ibin) = (rad + rad_next)/2
       else
          if (j < npart) then
             cycle
          else
             rad_grid(ibin) = rad
          endif
       endif

       bin_mass(ibin)             = no_in_bin*pmass !every bin has same mass.
       mass(ibin)                 = j*pmass !mass of paricles < r. Calculates cell outer total mass required by kepler.
       density(ibin)              = density_sum / no_in_bin
       temperature(ibin)          = temperature_sum / no_in_bin
       pressure(ibin)             = pressure_sum / no_in_bin
       int_eng(ibin)              = u_sum / no_in_bin
       velocity_3D(:,ibin)        = vel_sum(:) / no_in_bin !in cartesian coordinates
       composition_kepler(:,ibin) = composition_sum(:) / no_in_bin
       !radial momentum of the bin is calculated by using -
       rad_mom(ibin)     = rad_vel_sum
       y_e(ibin)         = (X_in/(1.*avogadro*atomic_mass_unit)) + (Y_in/(4.*avogadro*atomic_mass_unit))
       a_bar(ibin)       = X_in + (4.*Y_in) !average atomic mass in each bin.
       angular_vel_3D(:,ibin) = omega_sum(:) / moment_of_inertia !omega_sum is L_sum
       !calculating entropy
       !implementing entropy from the Sackur-Tetrode equation.
       entropy_array(ibin) = entropy(density(ibin)*unit_density,pressure(ibin)*unit_pressure,2,ierr)
       entropy_array(ibin) = entropy_array(ibin)/(kboltz*avogadro)
       if (ierr/=0) then
          print*, 'Entropy is calculated incorrectly'
       endif
       print*, no_in_bin,'no in bin',ibin,'ibin',number_particle,'no particle',no_in_bin,'no in bin'
       no_in_bin          = 0
       ibin               = ibin + 1
       density_sum        = 0.
       u_sum              = 0.
       temperature_sum    = 0.
       pressure_sum       = 0.
       vel_sum(:)         = 0.
       composition_sum(:) = 0.
       rad_vel_sum        = 0.
       omega_sum(:)       = 0.
       moment_of_inertia  = 0.
    endif
 enddo
end subroutine phantom_to_kepler_arrays
 !----------------------------------------------------------------
 !+
 !  This routine reads the output file that contains composition
 !  and saves it as a composition array that can be passed to the
 !  phantom_to_kepler_arrays subroutine.
 !+
 !----------------------------------------------------------------
subroutine composition_array(interpolate_comp,columns_compo,comp_label)
 !first read the file with compositon and save that into an array.
 use fileutils,only : get_nlines,skip_header,get_column_labels

 real, allocatable, intent(out)           :: interpolate_comp(:,:)
 character(len=20),allocatable,intent(out):: comp_label(:)
 integer                                  :: n_cols
 integer                                  :: n_rows,ierr,k,nheader
 integer, intent(out)                     :: columns_compo
 integer                                  :: n_labels
 character(len=10000)                     :: line
 character(len=120)                       :: filename
 logical                                  :: iexist

 columns_compo = 0
 n_rows = 0
 iexist = .false.

 filename = 'kepler.comp'
 !First check if kepler.comp exists.
 !This file will only be generated if KEPLER file had composition stored in it.
 inquire(file=filename,exist=iexist)
 if (.not.iexist) then
    write(*,'("No file with name ",A)') filename
    write(*,'("No composition to save  ",A)')
    allocate(interpolate_comp(columns_compo,n_rows))
    interpolate_comp(:,:) = 0.
    allocate(comp_label(columns_compo))
 else
    write(*,'("Reading composition from ",A)') filename
    n_rows = get_nlines(filename,skip_comments=.true.,n_columns=n_cols,n_headerlines=nheader)
    columns_compo = n_cols

    !Save composition read from file.
    allocate(interpolate_comp(columns_compo,n_rows))
    open(12, file=filename)
    ierr = 0
    !get column labels and send them back.
    read(12, '(a)', iostat=ierr) line
    allocate(comp_label(columns_compo))
    call get_column_labels(line,n_labels,comp_label)
    close(12)

    open(13, file=filename)
    call skip_header(13,nheader,ierr)
    do k = 1, n_rows
       read(13,*,iostat=ierr) interpolate_comp(:,k)
    enddo
    close(13)

    print*, '>>>>>> done'
 endif

end subroutine composition_array

 !-------

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
 use units,           only:umass,unit_velocity,udist
 use physcon,         only:gg,solarm
 integer,intent(in) :: npart
 real,intent(in)    :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)    :: pmass
 real :: xpos(3),vpos(3),pos(3),rad,percentage,count
 real :: v2,eps(npart)

 integer :: i, iorder(npart),j
 real    :: mh

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
    v2  = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))

    !Energy = KE + PE = v^2/2 - G mh/r.
    eps(i) = (v2**2*unit_velocity**2)/2. - (1e6*solarm*gg)/(rad*udist)

 enddo

 count = 0.

 do i = 1,npart
    !if energy is positive, then particle has been removed from star.
    if (eps(i) > 0) then
       count = count + 1.
    endif
 enddo

 ! print*,count,'count',npart,'npart'
 ! print*, count*umass*pmass, 'mass lost', npart*umass*pmass,umass,'mass unit'
 ! percentage = (count/npart)*100.
 ! print*, percentage,'mass lost percentage'
end subroutine find_mass_lost

end module analysis