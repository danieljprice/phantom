!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! analysis of surviving remnants after a partial tidal disruption event
!
! :References: Sharma, Price & Heger (2024), MNRAS 532, 89
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, dump_utils, eos, fileutils, io, linalg,
!   orbits, part, physcon, prompting, readwrite_dumps, sortutils, units,
!   vectorutils
!
 implicit none
 character(len=*), parameter, public :: analysistype = 'partial tde'
 public :: do_analysis

 private

contains

!-------------------------------------------------------------------------
!+
!  analysis of surviving remnants after a partial tidal disruption event
!+
!-------------------------------------------------------------------------
subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,              only:warning
 use dump_utils,      only:read_array_from_file
 use units,           only:udist,umass,unit_density,unit_velocity,utime
 use prompting,       only:prompt
 use readwrite_dumps, only:opened_full_dump
 integer,          intent(in) :: numfile,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: pmass,time
 character(len=*), intent(in) :: dumpfile
 integer :: i,j,columns_compo,ngrid
 real :: grid
 real, allocatable :: density(:),rad_grid(:),mass_enclosed(:),bin_mass(:),temperature(:)
 real, allocatable :: rad_vel(:),angular_vel_3D(:,:),composition_kepler(:,:)
 character(len=20), allocatable :: comp_label(:)
 character(len=120) :: output

 ngrid = 0
 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif

 write(*,'("Performing analysis type ",a)') analysistype
 write(*,'("Input file name is ",a)') dumpfile

 call phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,density,rad_grid,mass_enclosed,bin_mass,&
                              temperature,rad_vel,angular_vel_3D,composition_kepler,comp_label,columns_compo,ngrid,numfile)
 write(output,"(a4,i5.5)") 'ptok',numfile

 write(*,'("Output file name is ",a)') output
 open(iunit,file=output)
 write(iunit,'("# ",i5," # Version")') 10000
 write(iunit,'("# ",es20.12,"   # Dump file number")') time
 write(iunit,"('#',50(a22,1x))")                     &
          'grid',                                    &  ! grid number/bin number
          'cell mass',                               &  ! bin mass
          'cell outer tot. mass',                    &  ! total mass < r
          'radius',                                  &
          'cell density',                            &
          'cell temperature',                        &
          'cell radial vel',                    &
          'angular vel (x)',                         &  ! ang velocity x component
          'angular vel (y)',                         &  ! ang velocity y component
          'angular vel (z)',                         &  ! velocity z component
          comp_label                                    ! chemical composition
 print*,' Composition array shape: ',shape(composition_kepler)
 do i = 1, ngrid
    grid = i
    write(iunit,'(50(es18.10,1x))')                       &
              grid,                                       &
              bin_mass(i)*umass,                          &
              mass_enclosed(i)*umass,                     &
              rad_grid(i)*udist,                          &
              density(i)*unit_density,                    &
              temperature(i),                             &
              rad_vel(i)*unit_velocity,                   &
              (angular_vel_3D(j,i)/utime, j=1,3),         &
              (composition_kepler(j,i), j=1,columns_compo)
 enddo
 close(iunit)

end subroutine do_analysis

!----------------------------------------------------------------
!+
!  This subroutine bins the particles
!  Max density particle is considered as the centre of the remnant
!+
!----------------------------------------------------------------
subroutine phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,density,rad_grid,mass_enclosed,bin_mass,&
                                   temperature,rad_vel,angular_vel_3D,composition_kepler,comp_label,columns_compo,ibin,numfile)
 use units,        only:udist,umass,utime,unit_density
 use vectorutils,  only:cross_product3D
 use part,         only:rhoh
 use centreofmass, only:get_centreofmass
 use sortutils,    only:sort_by_radius
 use eos,          only:equationofstate,gmw,init_eos
 use physcon,      only:kb_on_mh,kboltz,atomic_mass_unit,avogadro,gg,pi,pc,years
 use linalg  ,     only:inverse
 integer, intent(in)               :: npart,numfile
 integer, intent(out)              :: ibin,columns_compo
 real, intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real, intent(in)                  :: pmass,time
 real, allocatable, intent(out)     :: rad_grid(:),density(:),mass_enclosed(:),bin_mass(:),temperature(:)
 real, allocatable, intent(out)     :: composition_kepler(:,:),rad_vel(:),angular_vel_3D(:,:)
 character(len=20), allocatable, intent(out) :: comp_label(:)
 integer :: i,j,location,ieos,ierr
 integer :: last_particle_with_neg_e,energy_verified_no
 integer :: dummy_bins,number_per_bin,count_particles,big_bins_no,tot_binned_particles
 integer :: index_sort,double_count
 integer :: iu_compo,iu_bigloop,iu_vphi,iu_vbra,iu_rot,iu_remnant
 integer :: iorder(npart)
 integer, allocatable :: sorted_index_npart(:),bound_index(:),sorted_index(:)
 real :: den_all(npart),xpos(3),vpos(3)
 real :: pos_mag_next
 real :: density_i,density_sum,rad_inner,radius_star
 real :: omega_particle
 real :: temperature_i,temperature_sum
 real :: rad_vel_i,momentum_i,rad_mom_sum
 real :: bhmass
 real :: i_matrix(3,3),I_sum(3,3),Li(3),L_i(3),L_sum(3),inverse_of_i(3,3)
 real :: L_reshape(3,1),matrix_result(3,1),omega(3)
 real :: pos_i,vel_i,pos_vec_i(3),vel_vec_i(3),ke_i,pe_i,tot_e_sum
 real :: tot_rem_mass,pos_com(3),vel_com(3),pos_com_mag,vel_com_mag
 real :: vphi_i,R_mag_i,vphi_sum,R_vec(2),vphi_avg,omega_vec(3),rad_cyl,breakup
 real :: ke_star,u_star,total_star,distance_from_bh,vel_at_infinity
 real, allocatable :: interpolate_comp(:,:),composition_i(:),composition_sum(:)
 real, allocatable :: temp_npart(:),den_npart(:),pos_npart(:),vel_npart(:)
 real, allocatable :: pos_vec_npart(:,:),vel_vec_npart(:,:)
 real, allocatable :: h_npart(:),tot_eng_npart(:),ke_npart(:),pe_npart(:)
 real, allocatable :: pos_wrt_bh(:,:),vel_wrt_bh(:,:),interp_comp_npart(:,:)
 logical :: double_the_no
 character(len=120) :: output

 dummy_bins = 5000
 ! use adiabatic EOS
 ieos = 2
 call init_eos(ieos,ierr)
 gmw=0.61
 ! Set mass of black hole in code units
 bhmass = 1
 ! performing a loop to determine maximum density particle position
 do j = 1,npart
    den_all(j) = rhoh(xyzh(4,j),pmass)
 enddo

 ! Save the location of max density particle
 location = maxloc(den_all,dim=1)
 print*,' Maximum density is on particle ',location
 ! Determining centre of star as max density particle.
 xpos(:) = xyzh(1:3,location)
 vpos(:) = vxyzu(1:3,location)
 distance_from_bh = sqrt(dot_product(xpos(:),xpos(:)))

 ! sorting particles by radius. Letting the max density particle be the centre of the star.
 ! This here for npart particles
 call sort_by_radius(npart,xyzh,iorder,xpos)

 ! Get the composition array for all the particles
 call composition_array(interpolate_comp,columns_compo,comp_label)

 ! Call the following function to obtain important information about each particle
 call calculate_npart_quantities(npart,iorder,numfile,xyzh,vxyzu,pmass,xpos,vpos,comp_label,&
                                       interpolate_comp,columns_compo,temp_npart,den_npart,pos_npart,vel_npart,&
                                       pos_vec_npart,vel_vec_npart,tot_eng_npart,sorted_index_npart,ke_npart,pe_npart,&
                                       pos_wrt_bh,vel_wrt_bh,h_npart,interp_comp_npart)

 ! This determines the particles bound to the star. Also removes the streams from the data
 call particles_bound_to_star(pos_npart,temp_npart,tot_eng_npart,npart,sorted_index_npart,&
                              bound_index,sorted_index,energy_verified_no,&
                              last_particle_with_neg_e,ke_npart,pe_npart,den_npart)

 ! This determines how many particles would be added to the bins. We set a min no of bins as 500 in the code and using that along with the particles
 ! that we consider as being part of the remnant implies that the model returns what the max number of particles I would have to add to a bin to get 500 bins
 call particles_per_bin(energy_verified_no,number_per_bin)
 big_bins_no          = number_per_bin
 tot_binned_particles = 0
 ibin                 = 1
 double_the_no        = .True.
 print*,' Bound particles: ',size(bound_index),' (energy verified: ',energy_verified_no,')'
 ! Define bins arrays
 allocate(density(dummy_bins),rad_grid(dummy_bins),mass_enclosed(dummy_bins),bin_mass(dummy_bins))
 allocate(temperature(dummy_bins),rad_vel(dummy_bins),angular_vel_3D(3,dummy_bins))
 allocate(composition_i(columns_compo),composition_sum(columns_compo),composition_kepler(columns_compo,dummy_bins))

 density_sum     = 0.
 temperature_sum = 0.
 rad_mom_sum     = 0.
 L_sum(:)        = 0.
 I_sum(:,:)      = 0.
 count_particles = 0
 composition_sum(:) = 0.
 composition_i(:)   = 0.
 pos_com(:) = 0.
 vel_com(:) = 0.
 tot_e_sum = 0.
 vphi_sum = 9.
 ! Write a comp file that includes information about the remnant only
 write(output,"(a4,i5.5)") 'compo',numfile
 open(newunit=iu_compo,file=output)
 write(iu_compo,"(32(a22,1x))") &
          "i",      &
          "ibin",   &
          "radius", &
          "x", &
          "y", &
          "z", &
          "radial_vel",&
          'temp',&
          'density',&
          comp_label,&
          'omega',&
          'breakup',&
          'j',&
          'index_sort'

 ! this will determine when sorted indices are part of the star. We would also need the normal i indices of the sorted particles
 ! Using this we can determine which sorted particles are part of the array and then use the sorted information to calculate all the
 ! quantities we require for the project
 print*,' Processing ',size(bound_index),' bound particles (',energy_verified_no,' energy verified)'
 open(newunit=iu_bigloop,file="big_loop_clean.txt")
 write(output,"(a4,i5.5)") 'vphi',numfile
 open(newunit=iu_vphi,file=output)
 write(iu_vphi,"(2(a22,1x))") &
         "rad",&
         "vphi"

 write(output,"(a4,i5.5)") 'vbra',numfile
 open(newunit=iu_vbra,file=output)
 write(iu_vbra,"(2(a22,1x))") &
         "rad",&
         "vbreak"
 write(output,"(a4,i5.5)") 'rot_i',numfile
 open(newunit=iu_rot,file=output)
 open(newunit=iu_remnant,file="remnant")
 write(iu_remnant,"(6(a22,1x))") &
          "x", &
          "y", &
          "z", &
          "m", &
          "h", &
          'rho'

 do i = 1,energy_verified_no
    ! get the sorted index of bound particles
    j = bound_index(i)
    index_sort = sorted_index(i)
    count_particles = count_particles + 1
    ! Calculate the position and velocity VEC of COM
    pos_com(:) = pos_com(:) + xyzh(1:3,index_sort)*pmass
    vel_com(:) = vel_com(:) + vxyzu(1:3,index_sort)*pmass

    ! Obtain the values for each particle that is bound/ part of remnant
    density_i     = den_npart(j)
    temperature_i = temp_npart(j)
    pos_i         = pos_npart(j)
    vel_i         = vel_npart(j)
    pos_vec_i(:)  = pos_vec_npart(:,j)
    vel_vec_i(:)  = vel_vec_npart(:,j)
    R_vec(:) = pos_vec_i(1:2)
    R_mag_i  = norm2(R_vec)
    ke_i = ke_npart(j)
    pe_i = pe_npart(i)
    write(iu_bigloop,*) i,j,pos_i,vel_i,pos_vec_i(:)*udist,temperature_i,density_i*unit_density,sorted_index(i)

    ! Calculate the angular velocity in cylindrical coordinates
    vphi_i = vel_vec_i(1)*(-pos_vec_i(2)/R_mag_i) + vel_vec_i(2)*(pos_vec_i(1)/R_mag_i)
    vphi_i = vphi_i/R_mag_i

    ! Position magnitude of the next bound particle
    if (i  /=  energy_verified_no) then
       pos_mag_next = pos_npart(j+1)
    endif

    ! composition
    if (columns_compo /= 0) then
       composition_i(:) = interp_comp_npart(:,j)
    endif

    if (index_sort == 13) then
       print*,' Debug: particle 13 composition = ',composition_i(:),' (j=',j,', i=',i,', index_sort=',index_sort,')'
    endif

    ! Calculate extra quantities
    if (pos_i > 0.) then
       ! Radial velocity
       rad_vel_i = dot_product(vel_vec_i(:),pos_vec_i)/pos_i
       momentum_i = rad_vel_i*pmass
    endif
    ! Angular momentum
    call cross_product3D(pos_vec_i(:),vel_vec_i(:),Li(:))
    L_i(:)   = Li(:)*pmass
    ! Moment of Inertia Matrix
    call moment_of_inertia(pos_vec_i,pos_i,pmass,i_matrix)

    if (pos_i == 0.) then
       omega_particle = 0.
    else
       omega_vec(:) = Li(:)/pos_i**2
       omega_particle = norm2(omega_vec)
    endif
    breakup = ((gg*i*pmass*umass)/(pos_i*udist)**3)**(0.5)
    write(iu_rot,*) pos_i,omega_particle/utime,vphi_i/utime,pos_vec_i(:)
    write(iu_compo,'(i9,1x,i5,1x,27(e18.10,1x),1x,i10,1x,i10)') &
              i, ibin, pos_i*udist, pos_vec_i(:)*udist, &
              rad_vel_i, temperature_i, density_i*unit_density, &
              composition_i(:), omega_particle/utime, breakup, j, index_sort

    write(iu_remnant,'(6(e18.10,1x))') pos_vec_i(:), pmass, h_npart(j), density_i

    ! Count particles keeps track of particles in a bin.
    ! Rad_inner is the radius of the first particle that is added to a bin
    if (count_particles == 1) then
       rad_inner = pos_i
    endif
    ! Calculate how many particles will go in a bin
    call no_per_bin(i,count_particles,double_the_no,number_per_bin,big_bins_no,&
      energy_verified_no,pos_mag_next,rad_inner,double_count)

    ! We sum the quantities we want to save for the particles
    density_sum        = density_sum + density_i
    temperature_sum    = temperature_sum + temperature_i
    rad_mom_sum        = rad_mom_sum + momentum_i
    L_sum(:)           = L_sum(:) + L_i(:)
    I_sum(:,:)         = I_sum(:,:) + i_matrix(:,:)
    composition_sum(:) = composition_sum(:) + composition_i(:)
    tot_e_sum          = ke_i + pe_i + tot_e_sum
    vphi_sum           = vphi_sum + vphi_i

    ! We check id the count_particles is the same as the number_per_bin
    ! If true then we save the bin information
    if (count_particles==number_per_bin) then
       ! Total particles binned. Should be the same as energy_verified_no at the end
       tot_binned_particles = tot_binned_particles+count_particles

       ! Calculate the bin quantities
       call radius_of_remnant(bound_index,count_particles,number_per_bin,i,energy_verified_no,&
        pos_npart,radius_star,pos_vec_npart,rad_cyl)

       rad_grid(ibin)      = radius_star
       density(ibin)       = density_sum/count_particles
       mass_enclosed(ibin) = tot_binned_particles*pmass
       bin_mass(ibin)      = count_particles*pmass
       ! Change the temperature of particles if its < 1.e3 to 1.e3
       if (temperature_sum < 1.e3) then
          print*,' WARNING! bin ',ibin,' has temperature < 1000 K (',temperature_sum,' K)'
       endif
       temperature(ibin)   = max(temperature_sum/count_particles,1e3)
       rad_vel(ibin)       = rad_mom_sum/bin_mass(ibin) ! radial vel of each bin is summation(vel_rad_i*m_i)/summation(m_i)
       if (count_particles == 1) then
          if (rad_grid(ibin) < tiny(0.)) then
             angular_vel_3D(:,ibin)  = L_sum(:)
          else
             angular_vel_3D(:,ibin) = L_sum(:) / (pos_i**2*pmass)
          endif
       else
          inverse_of_i  = inverse(I_sum, 3)
          L_reshape     = reshape(L_sum(:),(/3,1/))
          matrix_result = matmul(inverse_of_i,L_reshape)
          omega         = reshape(matrix_result,(/3/))
          angular_vel_3D(:,ibin) = omega
       endif
       composition_kepler(:,ibin) = composition_sum(:)/count_particles
       vphi_avg = vphi_sum/count_particles
       breakup = ((gg*mass_enclosed(ibin)*umass)/(rad_grid(ibin)*udist)**3)**(0.5)
       if (norm2(angular_vel_3D(:,ibin)) > 0) then
          write(iu_vphi,*) udist*rad_grid(ibin),norm2(angular_vel_3D(:,ibin))/utime
          write(iu_vbra,*) udist*rad_grid(ibin),breakup
       endif

       ! Reset the sum values
       count_particles = 0
       density_sum     = 0.
       temperature_sum = 0.
       rad_mom_sum     = 0.
       L_sum(:)        = 0.
       I_sum(:,:)      = 0.
       composition_sum(:) = 0.
       vphi_sum        = 0.
       ibin            = ibin+1
       number_per_bin =  big_bins_no
    endif
 enddo
 close(iu_vbra)
 close(iu_compo)
 close(iu_bigloop)
 close(iu_vphi)
 close(iu_rot)
 close(iu_remnant)
 ! We want to set the radial and angular velocity of the first bin as the same as the second bin so that they are not zero
 angular_vel_3D(:,1) = angular_vel_3D(:,2)
 rad_vel(1) = rad_vel(2)
 ibin = ibin-1
 tot_rem_mass = mass_enclosed(ibin)

 ! Get the COM pos and vel magnitudes
 call determine_pos_vel_com(vel_com,pos_com,pos_com_mag,vel_com_mag,tot_rem_mass)
 print*,' Radius of last particle: ',pos_i,' (code units)'
 print*,' Centre of mass position: ',pos_com
 print*,' Centre of mass velocity: ',vel_com
 print*,' Max density position: ',xpos
 print*,' Max density velocity: ',vpos
 print*,' Remnant mass: ',mass_enclosed(ibin)*umass
 ! Next we calculate the energy for the COM and determine if its bound or unbound
 call determine_bound_unbound(vel_com,pos_com,pos_com_mag,vel_com_mag,bhmass,tot_rem_mass,pmass,&
                              total_star,ke_star,u_star,vel_at_infinity)
 print*,' Total energy: ',tot_e_sum
 call write_dump_info(numfile,density(1),temperature(1),mass_enclosed(ibin),xpos,rad_grid(ibin),distance_from_bh,&
                         pos_com_mag,vel_com_mag,total_star,ke_star,u_star,time,vel_at_infinity)

end subroutine phantom_to_kepler_arrays

!----------------------------------------------------------------
!+
!  This subroutine returns the magnitude of the COM pos and vel
!+
!----------------------------------------------------------------
subroutine determine_pos_vel_com(vel_com,pos_com,pos_com_mag,vel_com_mag,tot_rem_mass)
 real, intent(inout),dimension(3) :: vel_com,pos_com
 real, intent(in)  :: tot_rem_mass
 real, intent(out) :: vel_com_mag,pos_com_mag

 ! Divide the pos_com and vel_com with the total mass enclosed
 pos_com(:) = pos_com(:)/tot_rem_mass
 vel_com(:) = vel_com(:)/tot_rem_mass

 pos_com_mag = norm2(pos_com)
 vel_com_mag = norm2(vel_com)

end subroutine determine_pos_vel_com

!----------------------------------------------------------------
!+
!  This subroutine returns if remnant is bound or unbound
!+
!----------------------------------------------------------------
subroutine determine_bound_unbound(vel_com,pos_com,pos_com_mag,vel_com_mag,bhmass,tot_rem_mass,pmass,&
                                   tot_energy_remnant_com,ke_star,pe_star,vel_at_infinity)
 use units,   only:udist,umass,unit_velocity
 use physcon, only:gg,au
 real, intent(in) :: vel_com_mag,pos_com_mag,bhmass,tot_rem_mass,pmass
 real, intent(in) :: pos_com(3),vel_com(3)
 real, intent(out) :: ke_star,pe_star,tot_energy_remnant_com,vel_at_infinity
 real :: bhmass_cgs,rem_mass,period_val
 real :: vel_com_cgs(3),pos_com_cgs(3)
 real :: er, ar

 bhmass_cgs = bhmass*umass
 rem_mass   = tot_rem_mass*umass
 vel_com_cgs(:) = vel_com(:)*unit_velocity
 pos_com_cgs(:) = pos_com(:)*udist

 ! Check if Total specific Energy of COM is < 0 or not (in cgs units)
 ke_star = 0.5*(vel_com_mag*unit_velocity)**2
 pe_star = -gg*bhmass_cgs/(pos_com_mag*udist)
 tot_energy_remnant_com = ke_star + pe_star
 print*,' Centre of mass velocity (cgs): ',vel_com_cgs
 print*,' Centre of mass position (cgs): ',pos_com_cgs

 if (tot_energy_remnant_com < 0.) then
    print*,' Remnant is BOUND to the black hole (energy = ',tot_energy_remnant_com,')'
    call determine_orbital_params(rem_mass,bhmass_cgs,pos_com_cgs,vel_com_cgs,period_val)
    ar = -gg*0.5*(bhmass_cgs + rem_mass)/tot_energy_remnant_com
    er = 1 - (56.77892268*udist)/ar
    print*,' Semi-major axis: ',ar/au,' AU, eccentricity: ',er
 elseif (tot_energy_remnant_com == 0.) then
    print*,' Parabolic orbit!'
 else
    print*,' Remnant is UNBOUND'
    call determine_inf_vel(tot_energy_remnant_com,vel_at_infinity)
    print*,' Velocity at infinity: ',vel_at_infinity*1e-5,' km/s'
    ar = gg*0.5*(bhmass_cgs + rem_mass)/tot_energy_remnant_com
    er = 1 + (56.77892268*udist)/ar
    print*,' Semi-major axis: ',ar/au,' AU, eccentricity: ',er
 endif
 print*,' Energy of centre of mass: ',pmass*(0.5*vel_com_mag**2 - (1/pos_com_mag))

end subroutine determine_bound_unbound

!----------------------------------------------------------------
!+
!  This subroutine returns the vel infinity for the remnant
!  if its unbound
!+
!----------------------------------------------------------------
subroutine determine_orbital_params(rem_mass,bhmass_cgs,pos_com,vel_com,period_val)
 use physcon, only:gg
 use orbits,  only:get_orbital_period,get_eccentricity
 real, intent(in)  :: rem_mass,bhmass_cgs,pos_com(3),vel_com(3)
 real, intent(out) :: period_val
 real :: ecc_val,mu

 mu = gg*(rem_mass+bhmass_cgs)  ! standard gravitational parameter G*(m1+m2), here in cgs units
 ecc_val = get_eccentricity(mu,pos_com,vel_com)
 print*,' Eccentricity: ',ecc_val
 print*,' Remnant mass: ',rem_mass,', BH mass: ',bhmass_cgs
 period_val = get_orbital_period(mu,pos_com,vel_com)
 print*,' Orbital period: ',period_val

end subroutine determine_orbital_params

!----------------------------------------------------------------
!+
!  This subroutine returns the orbital properties
!+
!----------------------------------------------------------------
subroutine determine_inf_vel(tot_energy_remnant_com,vel_at_infinity)
 real, intent(in)  :: tot_energy_remnant_com
 real, intent(out) :: vel_at_infinity

 vel_at_infinity = sqrt(2.*tot_energy_remnant_com)

end subroutine determine_inf_vel

!----------------------------------------------------------------
!+
!  This subroutine returns the position and velocity of a
!  particle wrt to the centre of star/max density point
!+
!----------------------------------------------------------------
subroutine particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos,vel,i,pos_mag,vel_mag)
 real, intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real, intent(in)                  :: xpos(3),vpos(3)
 integer, intent(in)               :: i
 real, intent(out)                 :: pos(3),vel(3),pos_mag,vel_mag

 pos(:) = xyzh(1:3,i) - xpos(:)
 vel(:) = vxyzu(1:3,i) - vpos(:)
 pos_mag = sqrt(dot_product(pos(:),pos(:)))
 vel_mag = sqrt(dot_product(vel(:),vel(:)))

end subroutine particle_pos_and_vel_wrt_centre

!----------------------------------------------------------------
!+
!  This subroutine returns which particles are bound to the star
!+
!----------------------------------------------------------------
subroutine particles_bound_to_star(pos_npart,temp_npart,tot_eng_npart,npart,sorted_index_npart,&
                                   bound_index,sorted_index,bound_particles_no,&
                                   last_particle_with_neg_e,ke_npart,pe_npart,den_npart)
 real, intent(in)    :: temp_npart(:),tot_eng_npart(:),ke_npart(:),pe_npart(:),pos_npart(:),den_npart(:)
 integer, intent(in) :: sorted_index_npart(:)
 integer, intent(in) :: npart
 integer, allocatable, intent(out) :: bound_index(:),sorted_index(:)
 integer, intent(out) ::  bound_particles_no,last_particle_with_neg_e
 integer :: energy_verified_no,i,count_loops_temp,iu
 integer, allocatable :: index_particle_star(:),index_bound(:),index_bound_sorted(:),index_bound_new(:)
 real :: max_temp,temp_cut
 real, allocatable :: temp_bound(:),temp_particles(:)
 logical :: temp_found,implement_temp_cut

 ! Implement temp cut would try to remove the streams. But if you only want
 ! to consider what is bound based on energy condition set this parameter to False
 max_temp = 8000.
 implement_temp_cut = .true.
 count_loops_temp = 0
 bound_particles_no = 0
 temp_found = .false.
 energy_verified_no = 0
 allocate(index_particle_star(npart),index_bound(npart),temp_particles(npart))
 open(newunit=iu,file="particle_index_clean")
 ! Use the sorted array information and check the energy condition first
 do i=1,npart
    ! if energy is less than 0, we have bound system. We can accept these particles.
    if (tot_eng_npart(i) < 0. .and. ke_npart(i) < 0.5*abs(pe_npart(i))) then
       write(iu,*) i,temp_npart(i),pos_npart(i),sorted_index_npart(i)
       energy_verified_no = energy_verified_no + 1
       ! Save the index of these particles
       ! this is because sometimes even if a particle is farther it could be could but the one before could be unbound
       last_particle_with_neg_e = i
       index_particle_star(energy_verified_no) = sorted_index_npart(i)
       index_bound(energy_verified_no) = i
       temp_particles(energy_verified_no) = temp_npart(i)
       print*,' Particle ',i,' is bound'
    endif
 enddo
 close(iu)
 allocate(temp_bound(energy_verified_no), index_bound_sorted(energy_verified_no),index_bound_new(energy_verified_no))
 do i = 1,energy_verified_no
    temp_bound(i) = temp_particles(i)
    ! This is the sorted index
    index_bound_sorted(i) = index_particle_star(i)
    index_bound_new(i) = index_bound(i)
 enddo
 if (implement_temp_cut) then
    ! next we loop over the bound particles based on energy condition to find the temp_cut
    ! As the models would need ages to evolve and I can not do that due to how slow some models run, we have streams around the remnants
    ! Hence, we bin the temperature particles and try to find the cut in temperature
    ! But using a temperature cut of 8000 K implies that if I use a model that has streams at high temperature (>1e4 K) because the remnant has just formed
    ! then I would not get rid of the correct particles
    ! Hence, we keep looping until the temperature being returned is the same as the max_temp
    count_loops_temp = count_loops_temp + 1
    call calculate_temp_cut(temp_bound,energy_verified_no,temp_cut,max_temp,temp_found,count_loops_temp,den_npart)
    max_temp = max_temp + 1000

    allocate(bound_index(energy_verified_no),sorted_index(energy_verified_no))
    ! use temp_cut to ignore the streams
    do i = 1,energy_verified_no
       if (temp_bound(i) > temp_cut) then
          bound_particles_no = bound_particles_no + 1
          ! Save the sorted array indices only
          bound_index(bound_particles_no) = index_bound_new(i)
          sorted_index(bound_particles_no) = index_bound_sorted(i)
          if (sorted_index(bound_particles_no) == 13) then
             print*,' Debug: particle 13 bound_index = ',bound_index(bound_particles_no)
          endif
       endif
    enddo
 else
    bound_particles_no = energy_verified_no
    allocate(bound_index(energy_verified_no),sorted_index(energy_verified_no))
    bound_index(:) = index_bound_new(:)
    sorted_index(:) = index_bound_sorted(:)
 endif

end subroutine particles_bound_to_star

!----------------------------------------------------------------
!+
!  This subroutine returns number of particles that can be put into
!  big bins
!+
!----------------------------------------------------------------
subroutine particles_per_bin(energy_verified_no,number_per_bin)
 integer, intent(in) :: energy_verified_no
 integer, intent(out) :: number_per_bin
 integer, parameter :: number_bins = 500

 ! Calculate the number of particles per bin
 number_per_bin = (energy_verified_no/number_bins)
 if (mod(energy_verified_no,number_bins) /=  0) then
    number_per_bin = 1+ number_per_bin
 endif

end subroutine particles_per_bin

!----------------------------------------------------------------
!+
!  This subroutine returns number of particles for each bin based
!  on some conditions
!+
!----------------------------------------------------------------
subroutine no_per_bin(j,count_particles,double_the_no,number_per_bin,big_bins_no,&
    energy_verified_no,pos_mag_next,rad_inner,double_count)
 integer, intent(inout) :: number_per_bin,double_count
 logical, intent(inout) :: double_the_no
 integer, intent(in)    :: count_particles,big_bins_no,j,energy_verified_no
 real, intent(in)       :: pos_mag_next,rad_inner
 real, parameter :: min_no=5
 integer :: i,iu
 real :: avg_val,diff_val

 avg_val = (pos_mag_next+rad_inner)/2
 diff_val = (pos_mag_next-rad_inner)

 open(newunit=iu,file="rad_to_bin",status='old',action='write',iostat=i)
 if (i /= 0) then
    ! File does not exist, create it
    open(newunit=iu,file="rad_to_bin",status='new',action='write',iostat=i)
 endif

 if (j==1) then
    number_per_bin = 1
    double_count = 1
 elseif (double_the_no .and. count_particles==1) then
    double_count = double_count*2
    number_per_bin = double_count
    if (number_per_bin >= big_bins_no) then
       number_per_bin = big_bins_no
       double_the_no = .False.
    endif
 else
    if (.not.double_the_no .and. j /= count_particles) then
       if (100*(pos_mag_next-rad_inner)/rad_inner > 30) then
          write(iu,*) pos_mag_next,rad_inner,j,number_per_bin
          number_per_bin=count_particles
       endif
    endif
 endif
 if (j==energy_verified_no) then
    number_per_bin = count_particles
 endif

end subroutine no_per_bin

!----------------------------------------------------------------
!+
!  This subroutine returns radius of the remnant
!+
!----------------------------------------------------------------
subroutine radius_of_remnant(bound_index,count_particles,number_per_bin,i,energy_verified_no,&
                             pos_npart,radius_star,pos_vec_npart,rad_cyl)
 integer, intent(in)    :: count_particles,number_per_bin,i,energy_verified_no,bound_index(:)
 real, intent(in)       :: pos_npart(:),pos_vec_npart(:,:)
 real, intent(out)      :: radius_star,rad_cyl
 integer :: index_val_next,index_val
 real    :: pos_mag_next,pos_mag
 real    :: pos_cyl,pos_cyl_next
 real    :: pos_cyl_vec(3),pos_cyl_vec_next(3)

 index_val = bound_index(i)
 index_val_next = bound_index(i+1)
 pos_mag = pos_npart(index_val)
 pos_cyl_vec(:) = pos_vec_npart(:,index_val)
 pos_cyl = sqrt(pos_cyl_vec(1)**2 + pos_cyl_vec(2)**2)
 if (count_particles == number_per_bin .and. i  /=  energy_verified_no) then
    pos_mag_next = pos_npart(index_val_next)
    pos_cyl_vec_next(:) = pos_vec_npart(:,index_val_next)

    radius_star = (pos_mag+pos_mag_next)/2
    pos_cyl_next = sqrt(pos_cyl_vec_next(1)**2 + pos_cyl_vec_next(2)**2)
    rad_cyl = (pos_cyl + pos_cyl_next)/2
 else
    radius_star = pos_mag
    rad_cyl = pos_cyl
 endif

end subroutine radius_of_remnant
!----------------------------------------------------------------
!+
!  This subroutine calculates the moment of inertia of each particle
!+
!----------------------------------------------------------------
subroutine moment_of_inertia(pos,pos_mag,pmass,i_matrix)
 real, intent(in)  :: pos(3),pos_mag,pmass
 real, intent(out) :: i_matrix(3,3)
 real :: delta(3,3),matrix1(3,1),matrix2(1,3),result_matrix(3,3)

 delta = reshape((/1,0,0,0,1,0,0,0,1/),shape(delta))
 i_matrix(:,:) = 0.
 matrix1=reshape(pos,shape(matrix1))
 matrix2=reshape(pos,shape(matrix2))
 result_matrix = matmul(matrix1,matrix2)
 i_matrix =  pmass*(pos_mag**2*delta - result_matrix)

end subroutine moment_of_inertia

!----------------------------------------------------------------
!+
! This subroutine calculates the temperature of all particles
! by sorting them out with radius
! Density is also sorted and saved. Along with the radius
!+
!----------------------------------------------------------------
subroutine calculate_npart_quantities(npart,iorder,numfile,xyzh,vxyzu,pmass,xpos,vpos,comp_label,&
                                       interpolate_comp,columns_compo,temp_npart,den_npart,pos_npart,vel_npart,&
                                       pos_vec_npart,vel_vec_npart,tot_eng_npart,sorted_index_npart,ke_npart,pe_npart,&
                                       pos_wrt_bh,vel_wrt_bh,h_npart,interp_comp_npart)
 use vectorutils,     only : cross_product3D
 use part,            only : rhoh,poten
 use eos,             only : equationofstate,gmw,init_eos
 use physcon,         only : gg
 integer, intent(in)               :: npart,numfile
 integer, intent(in)               :: iorder(:)
 real, intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real, intent(in)                  :: pmass
 real, intent(inout)               :: xpos(:),vpos(:)
 character(len=20), intent(in)     :: comp_label(:)
 real, intent(in)                  :: interpolate_comp(:,:)
 integer, intent(in)               :: columns_compo
 real, allocatable, intent(out)     :: temp_npart(:),den_npart(:),pos_npart(:),vel_npart(:)
 real, allocatable, intent(out)     :: pos_wrt_bh(:,:),vel_wrt_bh(:,:),h_npart(:)
 real, allocatable, intent(out)     :: pos_vec_npart(:,:),vel_vec_npart(:,:),tot_eng_npart(:)
 real, allocatable, intent(out)     :: ke_npart(:),pe_npart(:),interp_comp_npart(:,:)
 integer, allocatable, intent(out)  :: sorted_index_npart(:)
 integer             :: i,j,ierr,ieos
 real                :: pos(3),vel(3)
 real                :: potential_i, kinetic_i,energy_i,pos_mag,vel_mag
 real                :: density_i,temperature_i,eni_input,u_i
 real                :: ponrhoi,spsoundi,mu
 real, allocatable    :: composition_i(:)
 real, allocatable    :: A_array(:), Z_array(:)

 ieos = 2
 gmw = 0.61
 call init_eos(ieos,ierr)
 allocate(composition_i(columns_compo))
 call assign_atomic_mass_and_number(comp_label,A_array,Z_array)
 ! Allocate arrays to save the sorted index,density,temperature,radius,total energy of particle wrt centre, velocity_npart
 allocate(temp_npart(npart),den_npart(npart),pos_npart(npart),vel_npart(npart))
 allocate(sorted_index_npart(npart),tot_eng_npart(npart),ke_npart(npart),pe_npart(npart))
 allocate(pos_vec_npart(3,npart),vel_vec_npart(3,npart),pos_wrt_bh(3,npart),vel_wrt_bh(3,npart))
 allocate(h_npart(npart),interp_comp_npart(columns_compo,npart))

 do j = 1, npart
    ! access the rank of each particle in radius and save the sorted index
    i  = iorder(j)
    sorted_index_npart(j) = i

    ! the position of the particle is calculated by subtracting the point of
    ! highest density.
    ! xyzh is position wrt the black hole present at origin.
    call particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos,vel,i,pos_mag,vel_mag)
    potential_i  = poten(i)
    kinetic_i    = 0.5*pmass*vel_mag**2
    density_i    = rhoh(xyzh(4,i),pmass)
    energy_i     = potential_i + kinetic_i + vxyzu(4,i)*pmass
    print*,' Potential: ',potential_i,', Kinetic: ',kinetic_i

    ! composition
    if (columns_compo /= 0) then
       composition_i(:) = interpolate_comp(:,i)
    endif
    if (i == 13) then
       print*,' Debug: particle 13 composition = ',composition_i(:),' (i=',i,', j=',j,')'
    endif
    ! calculate mean molecular weight that is required by the eos module using
    ! the mass fractions for each particle.
    ! do not consider neutron which is the first element of the composition_i array.
    call calculate_mu(A_array,Z_array,composition_i,columns_compo,mu)
    gmw = 1./mu
    u_i       = vxyzu(4,i)
    eni_input = u_i
    call equationofstate(ieos,ponrhoi,spsoundi,density_i,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi=temperature_i,eni=eni_input)
    ! Save the information for each particle that we need
    den_npart(j) = density_i
    temp_npart(j) = temperature_i
    pos_npart(j) = pos_mag
    vel_npart(j) = vel_mag
    vel_vec_npart(:,j) = vel(:)
    pos_vec_npart(:,j) = pos(:)
    tot_eng_npart(j) = energy_i
    ke_npart(j) = kinetic_i
    pe_npart(j) = potential_i
    pos_wrt_bh(:,j) = xyzh(1:3,i)
    vel_wrt_bh(:,j) = vxyzu(1:3,i)
    h_npart(j) = xyzh(4,i)

    interp_comp_npart(:,j) = interpolate_comp(:,i)
 enddo

end subroutine calculate_npart_quantities

!----------------------------------------------------------------
!+
!  This routine reads the output file that contains composition
!  and saves it as a composition array that can be passed to the
!  phantom_to_kepler_arrays subroutine.
!+
!----------------------------------------------------------------
subroutine composition_array(interpolate_comp,columns_compo,comp_label)
 ! First read the file with composition and save that into an array.
 use fileutils, only:get_nlines,skip_header,get_column_labels
 real, allocatable, intent(out)           :: interpolate_comp(:,:)
 character(len=20), allocatable, intent(out) :: comp_label(:)
 integer, intent(out)                     :: columns_compo
 integer                                  :: n_cols,n_rows,ierr,k,nheader,n_labels,iu
 character(len=10000)                     :: line
 character(len=120)                       :: filename
 logical                                  :: iexist

 columns_compo = 0
 n_rows = 0
 iexist = .false.

 filename = 'tde.comp'
 ! first check if kepler.comp exists.
 ! this file will only be generated if KEPLER file had composition stored in it.
 inquire(file=filename,exist=iexist)
 if (.not.iexist) then
    write(*,"(a)") ' No file with name '//trim(filename)
    write(*,"(a)") ' -> no composition to save'
    allocate(interpolate_comp(columns_compo,n_rows))
    interpolate_comp(:,:) = 0.
    allocate(comp_label(columns_compo))
 else
    write(*,"(a)") ' Reading composition from '//trim(filename)
    n_rows = get_nlines(filename,skip_comments=.true.,n_columns=n_cols,n_headerlines=nheader)
    columns_compo = n_cols

    ! save composition read from file.
    allocate(interpolate_comp(columns_compo,n_rows))
    open(newunit=iu,file=filename)
    ierr = 0
    ! get column labels and send them back.
    read(iu, '(a)',iostat=ierr) line
    allocate(comp_label(columns_compo))
    call get_column_labels(line,n_labels,comp_label)
    close(iu)
    print*,' Composition labels: ',comp_label

    open(newunit=iu,file=filename)
    call skip_header(iu,nheader,ierr)
    do k = 1, n_rows
       read(iu,*,iostat=ierr) interpolate_comp(:,k)
    enddo
    close(iu)
    print*, '>>>>>> done reading composition file'
 endif

end subroutine composition_array

!----------------------------------------------------------------
!+
!  This routine is for assigning A and Z values based on the
!  element found in kepler.comp file
!+
!----------------------------------------------------------------
subroutine assign_atomic_mass_and_number(comp_label,A_array,Z_array)
 character(len=20), intent(in)   :: comp_label(:)
 real, allocatable, intent(out)  :: A_array(:), Z_array(:)
 integer                        :: size_to_allocate, i
 character(len=20), allocatable :: new_comp_label(:)

 if ( any( comp_label=="nt1" ) ) then
    size_to_allocate = size(comp_label(:))-1
 else
    size_to_allocate = size(comp_label(:))
 endif

 allocate(A_array(size_to_allocate), Z_array(size_to_allocate),new_comp_label(size_to_allocate))
 ! we skip nt1 as its neutron.
 new_comp_label = pack(comp_label,comp_label/="nt1")

 do i = 1, size_to_allocate
    select case (trim(new_comp_label(i)))
    case("h1")
       A_array(i) = 1
       Z_array(i) = 1
    case("he3")
       A_array(i) = 3
       Z_array(i) = 2
    case("he4")
       A_array(i) = 4
       Z_array(i) = 2
    case("c12")
       A_array(i) = 12
       Z_array(i) = 6
    case("n14")
       A_array(i) = 14
       Z_array(i) = 7
    case("o16")
       A_array(i) = 16
       Z_array(i) = 8
    case("ne20")
       A_array(i) = 20
       Z_array(i) = 10
    case("mg24")
       A_array(i) = 24
       Z_array(i) = 12
    case("si28")
       A_array(i) = 28
       Z_array(i) = 14
    case("s32")
       A_array(i) = 32
       Z_array(i) = 16
    case("ar36")
       A_array(i) = 36
       Z_array(i) = 18
    case("ca40")
       A_array(i) = 40
       Z_array(i) = 20
    case("ti44")
       A_array(i) = 44
       Z_array(i) = 22
    case("cr48")
       A_array(i) = 48
       Z_array(i) = 24
    case("fe52")
       A_array(i) = 52
       Z_array(i) = 26
    case("fe54")
       A_array(i) = 54
       Z_array(i) = 26
    case("ni56")
       A_array(i) = 56
       Z_array(i) = 28
    end select
 enddo
 print*,' A and Z arrays assigned'

end subroutine assign_atomic_mass_and_number
!----------------------------------------------------------------
!+
!  This routine is for calculating the gmw value by using
!  1/mu_i = Summation_a ((mass fraction)_a/A_a )*(1+Z_a)
!+
!----------------------------------------------------------------
subroutine calculate_mu(A_array,Z_array,composition_i,columns_compo,mu)
 real, allocatable, intent(in) :: A_array(:), Z_array(:), composition_i(:)
 integer, intent(in)         :: columns_compo
 real,    intent(out)        :: mu
 integer                     :: index_val

 mu = 0.
 if (columns_compo /= 0) then
    do index_val = 1, columns_compo-1
       mu = mu + (composition_i(index_val+1)*(1+Z_array(index_val)))/A_array(index_val)
    enddo
 endif

end subroutine calculate_mu

!----------------------------------------------------------------
!+
!  This routine updates the dump_info file with the information
!  for full dumps files
!+
!----------------------------------------------------------------
subroutine write_dump_info(fileno,density,temperature,mass,xpos,rad,distance,pos_mag_star,vel_mag_star,&
                tot_energy,kinetic_energy,potential_energy,time,vel_at_infinity)
 use units,   only:udist,umass,unit_velocity,utime,unit_density
 use physcon, only:years,km
 integer, intent(in) :: fileno
 real, intent(in) :: vel_at_infinity,density,time,temperature,mass
 real, intent(in) :: xpos(3),rad,distance,pos_mag_star,vel_mag_star
 real, intent(in) :: tot_energy,kinetic_energy,potential_energy
 integer :: status,file_id
 character(len=10) :: filename
 logical :: file_exists

 ! set a file name
 filename = 'dump_info'

 ! check if the file exists
 inquire(file=filename, exist=file_exists)

 ! open the file for appending or creating
 if (file_exists) then
    open(unit=file_id,file=filename,status='old', position="append",action="write",iostat=status)
    if (status /= 0) then
       write(*,*) 'Error opening file: ', filename
       stop
    endif

 else
    open(unit=file_id,file=filename,status='new',action='write',iostat=status)
    if (status /= 0) then
       write(*,*) 'Error creating file: ', filename
       stop
    endif
    ! Write headers to file
    write(file_id,'(17(a22,1x))') &
              "FileNo", &
              "Density",&
              "Temperature",&
              "Mass",&
              "x",&
               "y",&
               "z",&
               "radius",&
               "DistanceFromBH",&
               "Pos Mag",&
               "Vel Star",&
               "specTotEn",&
               "specKE",&
               "specPE",&
               "time",&
               "Escape_in",&
               "Accretion_r"
 endif
 write(file_id,'(i5,1x,16(e18.10,1x))') fileno,density*unit_density,temperature,mass*umass,&
       xpos(:)*udist,rad*udist,distance*udist,pos_mag_star*udist,&
       vel_mag_star*unit_velocity,tot_energy,kinetic_energy,potential_energy,time*utime,&
       vel_at_infinity/km,(mass*umass)/(time/years*utime)
 close(file_id)

end subroutine write_dump_info

!----------------------------------------------------------------
!+
!  This subroutine is to get the temperature cut
!
!+
!----------------------------------------------------------------
subroutine calculate_temp_cut(temperature_array,count_bound,temp_cut,max_temp,temp_found,&
                              count_loops_temp,density_array)
 real, intent(in) :: temperature_array(:),max_temp,density_array(:)
 integer, intent(in) :: count_bound,count_loops_temp
 real, intent(out)   :: temp_cut
 logical, intent(inout) :: temp_found
 integer, parameter :: nbins = 20000
 integer :: i,count_possible_temp,m,count_cut_index
 real :: temp_start,count_temp_particles,dtemp
 real :: mean,variance,std,cut_off
 real :: count_cut,lower_limit,upper_limit
 real, dimension(nbins) ::temp_array_test
 real, allocatable :: avg_density(:)
 real, allocatable :: temp_array_new(:),count_particles_temp(:)

 ! First we create an array of possible temperature from max_temp to 0 with a step size of 100.
 temp_start = 0.
 dtemp = 100.

 count_temp_particles = 0
 count_cut_index = 0
 count_cut = 0.
 count_possible_temp=1+nint(max_temp/dtemp)

 ! Create array with the temperatures ranging from 0 to max_temp
 do m=1,nbins
    if (temp_start <= max_temp) then
       temp_array_test(m) = temp_start
       temp_start = temp_start + dtemp
    endif
 enddo

 ! Allocate arrays to save the number of particles per bin
 allocate(temp_array_new(count_possible_temp),count_particles_temp(count_possible_temp))
 allocate(avg_density(count_possible_temp))

 count_particles_temp(:) = 0

 ! Next we create the same size array as count_possible_temp
 do m=1,count_possible_temp
    temp_array_new(m) = temp_array_test(m)
 enddo

 ! this will count the particles for each temperature and then save them into a new array
 do i =1,count_bound
    do m=1,size(temp_array_new)-1
       if (temperature_array(i) >= temp_array_new(m) .and. temperature_array(i) < temp_array_new(m+1) )  then
          count_temp_particles = count_particles_temp(m) + 1
          count_particles_temp(m) =  count_temp_particles
          avg_density(m) = density_array(i)
       endif
    enddo
 enddo

 print*,' Temperature cut analysis:'
 print*,'     Temperature array size: ',size(temp_array_new)
 print*,'  Particle count array size: ',size(count_particles_temp)
 ! Calculate the mean, std of the data
 call statistics(count_particles_temp,mean,variance,std)

 ! Using 2 sigma as the data sample is small to determine the outlier
 cut_off = std*2
 lower_limit = mean - cut_off
 upper_limit = mean + cut_off

 ! This loops and find the last element which is outside the limits based on 2 sigma
 do i=1,size(count_particles_temp)
    if (count_particles_temp(i) > upper_limit .or. count_particles_temp(i) < lower_limit) then
       count_cut = count_particles_temp(i)
       count_cut_index = i
    endif
 enddo
 print*,' Initial temperature cut: ',count_cut,' at index ',count_cut_index
 ! this starts from the cound_cut_index found earlier but then tries to make sure that the cut is done when the gaussian bins
 ! have less than 5% particles compared to the max_temp_cut found above
 do i=count_cut_index,size(count_particles_temp)
    if ((count_particles_temp(i)/count_cut)*100 < 1.) then
       count_cut = count_particles_temp(i)
       print*,' Temperature cut updated: ',count_cut,' (',(count_particles_temp(i)/count_cut)*100,'% of max)'
       count_cut_index = i
       exit
    endif
 enddo


 ! Define the temperature to cut the model at
 temp_cut = temp_array_new(count_cut_index)

 if (temp_cut  /=   max_temp) then
    temp_found = .true.
 endif

 ! If we get the temp_cut as 0. K and the count_loops_temp is 1, then we accept that as a true value
 if (temp_cut == 0.0 .and. count_loops_temp /= 1) then
    temp_found = .false.
 endif
 print*,' Final temperature cut: ',temp_cut,' K'
end subroutine calculate_temp_cut

!--------------------------------------------------------------------
!+
!  computes the mean, variance and standard deviation
!+
!--------------------------------------------------------------------
subroutine statistics(array_data,mean,variance,std)
 real, allocatable, intent(in) :: array_data(:)
 real, intent(out) :: mean,variance,std
 integer :: size_array,i
 real :: var,sum_val

 sum_val = 0.
 var = 0.
 size_array = size(array_data)
 do i=1,size_array
    sum_val = sum_val + array_data(i)
 enddo
 mean = sum_val/size_array

 do i=1,size_array
    var = var + (array_data(i) - mean)**2
 enddo

 variance = var/(size_array-1)
 std = sqrt(variance)

end subroutine statistics

end module analysis
