!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
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
! :Dependencies: centreofmass, dump_utils, eos, fileutils, io, linalg,
!   orbits_data, part, physcon, prompting, readwrite_dumps, sortutils,
!   units, vectorutils
!
 implicit none
 character(len=3), parameter, public :: analysistype = 'tde'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

 use io,              only : warning
 use dump_utils,      only : read_array_from_file
 use units,           only : udist,umass,unit_density,unit_ergg,unit_velocity,utime !units required to convert to kepler units.
 use prompting,       only : prompt
 use readwrite_dumps, only : opened_full_dump

 integer,  intent(in) :: numfile,npart,iunit
 integer              :: i,j,columns_compo
 integer              :: ngrid = 0

 real                              :: grid
 real,intent(in)                   :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)                   :: pmass,time
 real , allocatable                :: density(:),rad_grid(:),mass_enclosed(:),bin_mass(:),temperature(:),rad_vel(:),angular_vel_3D(:,:)
 real, allocatable                 :: composition_kepler(:,:)
 character(len=20),allocatable     :: comp_label(:)
 character(len=120)                :: output
 character(len=*),intent(in)       :: dumpfile

 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif

 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 call phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,density,rad_grid,mass_enclosed,bin_mass,&
                              temperature,rad_vel,angular_vel_3D,composition_kepler,comp_label,columns_compo,ngrid,numfile)
 write(output,"(a4,i5.5)") 'ptok',numfile

 write(*,'("Output file name is ",A)') output
 open(iunit,file=output)
 write(iunit,'("# ",i5," # Version")') 10000
 write(iunit,'("# ",es20.12,"   # Dump file number")') time
 write(iunit,"('#',50(a22,1x))")                     &
          'grid',                                    &  !grid number/bin number
          'cell mass',                               &  !bin mass
          'cell outer tot. mass',                    &  !total mass < r
          'radius',                                  &
          'cell density',                            &
          'cell temperature',                        &
          'cell radial vel',                    &
          'angular vel (x)',                         &  !ang velocity x component
          'angular vel (y)',                         &  !ang velocity y component
          'angular vel (z)',                         &  !velocity z component
          comp_label                                    !chemical composition
 print*, shape(composition_kepler),'kepler compo'
 do i = 1, ngrid
    grid = i
    write(iunit,'(50(es18.10,1X))')                       &
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
 !  Max density particle is considered as the centre of the remanant
 !
 !+
 !----------------------------------------------------------------
subroutine phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,density,rad_grid,mass_enclosed,bin_mass,&
                                   temperature,rad_vel,angular_vel_3D,composition_kepler,comp_label,columns_compo,ibin,numfile)
 use units , only: udist,umass,unit_velocity,utime,unit_energ,unit_density
 use vectorutils,     only : cross_product3D
 use part,            only : rhoh,poten
 use centreofmass,    only : get_centreofmass
 use sortutils,       only : set_r2func_origin,indexxfunc,r2func_origin
 use eos,             only : equationofstate,entropy,X_in,Z_in,gmw,init_eos
 use physcon,         only : kb_on_mh,kboltz,atomic_mass_unit,avogadro,gg,pi,pc,years
 use orbits_data,     only : escape,semimajor_axis,period_star
 use linalg  ,        only : inverse
 integer,intent(in)               :: npart,numfile
 integer,intent(out)              :: ibin,columns_compo
 real,intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)                  :: pmass,time
 real,intent(out),allocatable     :: rad_grid(:),density(:),mass_enclosed(:),bin_mass(:),temperature(:),rad_vel(:),angular_vel_3D(:,:)
 real,allocatable,intent(out)     :: composition_kepler(:,:)
 character(len=20),allocatable,intent(out) :: comp_label(:)
 real :: den_all(npart),xpos(3),vpos(3)
 integer :: i,j,location,iorder(npart),next_particle,ieos,ierr
 character(len=120)                :: output
 real :: potential_wrt_bh,kinetic_wrt_bh,tot_wrt_bh
 real :: pos(3),vel(3),kinetic_i,potential_i,energy_i,vel_mag,pos_mag,pos_next(3),vel_next(3),vel_mag_next,pos_mag_next
 integer ::particle_bound_bh,last_particle_with_neg_e,energy_verified_no,index_val,i_next,iu1,iu2,iu3,i_prev
 integer,allocatable :: index_particle_star(:),array_particle_j(:),array_bh_j(:)
 integer :: dummy_size,dummy_bins=5000,number_per_bin,count_particles,number_bins,no_particles,big_bins_no,tot_binned_particles
 real :: density_i,density_sum,rad_inner,rad_outer,radius_star
 logical :: double_the_no,escape_star
 real :: omega_particle,omega_bin,pos_mag_star,vel_mag_star
 real :: eni_input,u_i,temperature_i,temperature_sum,mu
 real :: ponrhoi,spsoundi,rad_vel_i,momentum_i,rad_mom_sum
 real :: bhmass,pos_prev(3),vel_prev(3),pos_mag_prev,vel_mag_prev
 real :: i_matrix(3,3),I_sum(3,3),Li(3),L_i(3),L_sum(3),inverse_of_i(3,3),L_reshape(3,1),matrix_result(3,1),omega(3)
 real,allocatable    :: A_array(:), Z_array(:)
 real,allocatable    :: interpolate_comp(:,:),composition_i(:),composition_sum(:)
 real :: ke_star,u_star,total_star,distance_from_bh,vel_from_bh,vel_at_infinity
 real :: period_val
 real,allocatable :: count_particles_temp(:),temp_array_new(:),temp_array_diff(:),temp_all_particles(:)
 real :: max_temp=8000.,temp_cut=0.
 real, dimension(200) :: temp_array_test
 logical :: temp_found=.false.
 integer :: count_loops_temp = 0
 real,allocatable :: temp_npart(:),den_npart(:),pos_npart(:),vel_npart(:),pos_vec_npart(:,:),vel_vec_npart(:,:)
 real,allocatable :: h_npart(:),tot_eng_npart(:),ke_npart(:),pe_npart(:)
 integer,allocatable :: sorted_index_npart(:),bound_index(:),sorted_index(:)
 real :: pos_i,vel_i,pos_vec_i(3),vel_vec_i(3),ke_i,pe_i,tot_e_sum
 real :: tot_rem_mass,pos_com(3),vel_com(3),pos_com_mag,vel_com_mag
 integer :: index_sort,double_count
 real,allocatable :: pos_wrt_bh(:,:),vel_wrt_bh(:,:),interp_comp_npart(:,:)
 real :: vphi_i,R_mag_i,vphi_sum,R_vec(2),vphi_avg,omega_vec(3),rad_cyl,breakup

 ! use adiabatic EOS
 ieos = 2
 call init_eos(ieos,ierr)
 gmw=0.61
 ! Set mass of black hole in code units
 bhmass = 1
 ! Set initial cut based on temperature as zero K
 temp_cut = 0.
 allocate(temp_all_particles(npart))

 ! performing a loop to determine maximum density particle position
 do j = 1,npart
    den_all(j) = rhoh(xyzh(4,j),pmass)
 enddo

 ! Save the location of max density particle
 location = maxloc(den_all,dim=1)
 print*,location,"LOCATION OF MAX"
 ! Determining centre of star as max density particle.
 xpos(:) = xyzh(1:3,location)
 vpos(:) = vxyzu(1:3,location)
 distance_from_bh = sqrt(dot_product(xpos(:),xpos(:)))
 vel_from_bh = sqrt(dot_product(vpos(:),vpos(:)))

 ! sorting particles by radius. Letting the max density particle be the centre of the star.
 ! This here for npart particles
 call set_r2func_origin(xpos(1),xpos(2),xpos(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)

 ! Get the composition array for all the particles
 call composition_array(interpolate_comp,columns_compo,comp_label)

 ! Call the following function to obatin important information about each particle
 call calculate_npart_quantities(npart,iorder,numfile,xyzh,vxyzu,pmass,xpos,vpos,comp_label,&
                                       interpolate_comp,columns_compo,temp_npart,den_npart,pos_npart,vel_npart,&
                                       pos_vec_npart,vel_vec_npart,tot_eng_npart,sorted_index_npart,ke_npart,pe_npart,&
                                       pos_wrt_bh,vel_wrt_bh,h_npart,interp_comp_npart)

 ! This determines the particles bound to the star. Also removes the streams from the data
 call particles_bound_to_star(pos_npart,temp_npart,tot_eng_npart,npart,sorted_index_npart,bound_index,sorted_index,energy_verified_no,&
                                    last_particle_with_neg_e,ke_npart,pe_npart,den_npart)

 ! This determins how many particles would be added to the bins. We set a min no of bins as 500 in the code and using that along with the particles
 ! that we consider as being part of the remnant implies that the model returns what the max number of particles I would have to add to a bin to get 500 bins
 call particles_per_bin(energy_verified_no,number_per_bin)
 big_bins_no          = number_per_bin
 tot_binned_particles = 0
 no_particles         = 1
 ibin                 = 1
 double_the_no        = .True.
 print*,size(bound_index),"BOUND INDEX ARRAY",energy_verified_no,"energy verified no"
 ! Define bins arrays
 allocate(density(dummy_bins),rad_grid(dummy_bins),mass_enclosed(dummy_bins),bin_mass(dummy_bins),temperature(dummy_bins),rad_vel(dummy_bins),angular_vel_3D(3,dummy_bins))
 allocate(composition_i(columns_compo),composition_sum(columns_compo),composition_kepler(columns_compo,dummy_bins))

 print*,"WORKED1"
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
 ! Write a comp file that includes informtation about the remnant only
 write(output,"(a4,i5.5)") 'compo',numfile
 open(4,file=output)
 write(4,"(32(a22,1x))") &
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



 ! this will determine when sorted indices are part of the star. We would also need the normal i indicies of the sorted particles
 ! Using this we can determine which sorted particles are part of the array and then use the sorted information to calculate all the
 ! quantities we require for the project
 print*,size(bound_index),"bound index size",energy_verified_no,"energy_verified_no"
 open(14,file="big_loop_clean.txt")
 write(output,"(a4,i5.5)") 'vphi',numfile
 open(10,file=output)
 write(10,"(2(a22,1x))") &
         "rad",&
         "vphi"

 write(output,"(a4,i5.5)") 'vbra',numfile
 open(111,file=output)
 write(111,"(2(a22,1x))") &
         "rad",&
         "vbreak"
 write(output,"(a4,i5.5)") 'rot_i',numfile
 open(11,file=output)
 open(41,file="remnant")
 write(41,"(6(a22,1x))") &
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
    write(14,*) i,j,pos_i,vel_i,pos_vec_i(1)*udist,pos_vec_i(2)*udist,pos_vec_i(3)*udist,temperature_i,density_i*unit_density,sorted_index(i)

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
       print*,composition_i(:),"compo in big look",j,"j",i,"i",index_sort,"index_sort2"
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
    write(11,*) pos_i,omega_particle/utime,vphi_i/utime,pos_vec_i(1),pos_vec_i(2),pos_vec_i(3)
    write(4,'(i9,1x,i5,1x,27(e18.10,1x),1x,i10,1x,i10)') &
              i, &
              ibin, &
              pos_i*udist, &
              pos_vec_i(1)*udist, &
              pos_vec_i(2)*udist, &
              pos_vec_i(3)*udist, &
              rad_vel_i,&
              temperature_i,&
              density_i*unit_density,&
              composition_i(:),&
              omega_particle/utime,&
              breakup,&
              j,&
              index_sort

    write(41,'(6(e18.10,1x))') &
           pos_vec_i(1), &
           pos_vec_i(2), &
           pos_vec_i(3), &
           pmass,  &
           h_npart(j), &
           density_i

    ! Count particles keeps track of particles in a bin.
    ! Rad_inner is the radius of the first particle that is added to a bin
    if (count_particles == 1) then
       rad_inner = pos_i
    endif
    ! Calculate how many particles will go in a bin
    call no_per_bin(i,count_particles,double_the_no,number_per_bin,big_bins_no,energy_verified_no,pos_mag_next,rad_inner,double_count)

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
    !if (count_particles==number_per_bin .or. i==energy_verified_no) then
    if (count_particles==number_per_bin) then
       ! Total particles binned. Should be the same as energy_verified_no at the end
       tot_binned_particles = tot_binned_particles+count_particles

       ! Calculate the bin quantities
       call radius_of_remnant(bound_index,count_particles,number_per_bin,i,energy_verified_no,pos_npart,radius_star,pos_vec_npart,rad_cyl)

       rad_grid(ibin)      = radius_star
       density(ibin)       = density_sum/count_particles
       mass_enclosed(ibin) = tot_binned_particles*pmass
       bin_mass(ibin)      = count_particles*pmass
       ! Change the temperature of particles if its < 1.e3 to 1.e3
       if (temperature_sum < 1.e3) then
          print*,"THIS BIN HAS TEMP LESS THAN 1000 K",temperature_sum
       endif
       temperature(ibin)   = max(temperature_sum/count_particles,1e3)
       rad_vel(ibin)       = rad_mom_sum/bin_mass(ibin) !Radial vel of each bin is summation(vel_rad_i*m_i)/summation(m_i)
       if (count_particles == 1) then
          if (rad_grid(ibin)==0.) then

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
          write(10,*) udist*rad_grid(ibin),norm2(angular_vel_3D(:,ibin))/utime
          write(111,*) udist*rad_grid(ibin),breakup
       endif

       !print*,count_particles,"count particles",ibin,"ibin",rad_grid(ibin),"rad",number_per_bin,"number per bin"
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
 close(111)
 close(4)
 close(14)
 close(10)
 close(11)
 close(41)
 ! We want to set the radial and angular velocity of the first bin as the same as the second bin so that they are not zero
 angular_vel_3D(:,1) = angular_vel_3D(:,2)
 rad_vel(1) = rad_vel(2)
 ibin = ibin-1
 tot_rem_mass = mass_enclosed(ibin)

 ! Get the COM pos and vel magnitudes
 call determine_pos_vel_com(vel_com,pos_com,pos_com_mag,vel_com_mag,tot_rem_mass)
 print*,pos_i,"Radius of last particle in code units"
 print*,pos_com,"POS COM",vel_com,"VEL COM"
 print*,xpos,"XPOS",vpos,"VPOS",mass_enclosed(ibin)*umass,"mass"
 print*,norm2(xpos),norm2(pos_com),"pos mag",norm2(pos_com)/170.33676805,"how far from rt?"
 print*,norm2(vpos),norm2(vel_com),"vel mag"
 ! Next we calculate the energy for the COM and determine if its bound or unbound
 call determine_bound_unbound(vel_com,pos_com,pos_com_mag,vel_com_mag,bhmass,tot_rem_mass,pmass,&
                              total_star,ke_star,u_star,vel_at_infinity)
 print*,tot_e_sum,"TOT E SUM"
 call write_dump_info(numfile,density(1),temperature(1),mass_enclosed(ibin),xpos,rad_grid(ibin),distance_from_bh,&
                         pos_com_mag,vel_com_mag,total_star,ke_star,u_star,time,vel_at_infinity)

end subroutine phantom_to_kepler_arrays
 !----------------------------------------------------------------
 !+
 !  This subroutine returns the magntitude of the COM pos and vel
 !+
 !----------------------------------------------------------------
subroutine determine_pos_vel_com(vel_com,pos_com,pos_com_mag,vel_com_mag,tot_rem_mass)
 real,intent(inout),dimension(3) :: vel_com,pos_com
 real,intent(in)  :: tot_rem_mass
 real,intent(out) :: vel_com_mag,pos_com_mag

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
 use units , only : udist,umass,unit_velocity
 use physcon,only : gg

 real,intent(in) :: vel_com_mag,pos_com_mag,bhmass,tot_rem_mass,pmass
 real,intent(in) :: pos_com(3),vel_com(3)
 real,intent(out):: ke_star,pe_star,tot_energy_remnant_com,vel_at_infinity
 real :: bhmass_cgs,rem_mass
 real :: period_val,vel_com_cgs(3),pos_com_cgs(3)
 real :: er, ar

 bhmass_cgs = bhmass*umass
 rem_mass   = tot_rem_mass*umass
 vel_com_cgs(:) = vel_com(:)*unit_velocity
 pos_com_cgs(:) = pos_com(:)*udist
 ! Check if Total specific Energy of COM is < 0 or not (in cgs units)
 ke_star = 0.5*(vel_com_mag*unit_velocity)**2
 pe_star = -gg*bhmass_cgs/(pos_com_mag*udist)
 tot_energy_remnant_com = ke_star + pe_star
 print*,vel_com_cgs,"CGS vel com",pos_com_cgs,"CGS pos com"

 if (tot_energy_remnant_com < 0.) then
    print*, "REMNANT IS BOUND TO THE BLACKHOLE",tot_energy_remnant_com,"energy val"
    call determine_orbital_params(rem_mass,bhmass_cgs,pos_com_cgs,vel_com_cgs,period_val)
    ar = -gg*0.5*(bhmass_cgs + rem_mass)/tot_energy_remnant_com
    er = 1 - (56.77892268*udist)/ar
    print*,"******************"
    print*,ar/1.496e13,"ar",er,"er"
 elseif (tot_energy_remnant_com == 0.) then
    print*, "Parabolic orbit!"
 else
    print*, "REMNANT IS UNBOUND"
    call determine_inf_vel(tot_energy_remnant_com,vel_at_infinity)
    print*,"VELOCITY OF REMNANT IN kms/s :",vel_at_infinity*1e-5
    ar = gg*0.5*(bhmass_cgs + rem_mass)/tot_energy_remnant_com
    er = 1 + (56.77892268*udist)/ar
    print*,"******************"
    print*,ar/1.496e13,"ar",er,"er"
 endif

 print*,pmass*(0.5*vel_com_mag**2 - (1/pos_com_mag)),"ENERGY OF COM"
end subroutine determine_bound_unbound
 !----------------------------------------------------------------
 !+
 !  This subroutine returns the vel infinity for the remnant
 !  if its unbound
 !+
 !----------------------------------------------------------------
subroutine determine_orbital_params(rem_mass,bhmass_cgs,pos_com,vel_com,period_val)
 use orbits_data,     only : escape,semimajor_axis,period_star,eccentricity_star
 real,intent(in) :: rem_mass,bhmass_cgs,pos_com(3),vel_com(3)
 real,intent(out):: period_val
 real :: ecc_val

 ecc_val = eccentricity_star(rem_mass,bhmass_cgs,pos_com,vel_com)
 print*,ecc_val,"ECCENTRICITY VALUE!!!!",rem_mass,"rem mass", bhmass_cgs,"bhmass cgs",pos_com,"com pos",vel_com,"com vel"
 period_val = period_star(rem_mass,bhmass_cgs,pos_com,vel_com)
 print*,period_val,"PERIOD OF STAR"

end subroutine determine_orbital_params
 !----------------------------------------------------------------
 !+
 !  This subroutine returns the oribital properties
 !+
 !----------------------------------------------------------------
subroutine determine_inf_vel(tot_energy_remnant_com,vel_at_infinity)
 real,intent(in) :: tot_energy_remnant_com
 real,intent(out):: vel_at_infinity

 vel_at_infinity = sqrt(2.*tot_energy_remnant_com)

end subroutine determine_inf_vel
 !----------------------------------------------------------------
 !+
 !  This subroutine returns the position and velocity of a
 !  particle wrt to the centre of star/max density point
 !+
 !----------------------------------------------------------------
subroutine particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos,vel,i,pos_mag,vel_mag)
 real,intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)                  :: xpos(3),vpos(3)
 integer,intent(in)               :: i
 real,intent(out)                 :: pos(3),vel(3),pos_mag,vel_mag

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
subroutine particles_bound_to_star(pos_npart,temp_npart,tot_eng_npart,npart,sorted_index_npart,bound_index,sorted_index,bound_particles_no,&
                                    last_particle_with_neg_e,ke_npart,pe_npart,den_npart)

 real,intent(in)    :: temp_npart(:),tot_eng_npart(:),ke_npart(:),pe_npart(:),pos_npart(:),den_npart(:)
 integer,intent(in) :: sorted_index_npart(:)
 integer,intent(in) :: npart

 integer,allocatable,intent(out) :: bound_index(:),sorted_index(:)
 integer,intent(out) ::  bound_particles_no,last_particle_with_neg_e
 integer :: energy_verified_no,i
 real,allocatable :: index_particle_star(:),temp_bound(:),temp_particles(:)
 integer,allocatable :: index_bound(:),index_bound_sorted(:),index_bound_new(:)
 real :: max_temp=8000.,index_val
 integer :: count_loops_temp=0
 logical :: temp_found,implement_temp_cut
 real :: temp_cut

 ! Implement temp cut would try to remove the strems. But if you only want
 ! to consider what is bound based on energy condition set this parameter to False
 implement_temp_cut = .true.
 bound_particles_no = 0
 temp_found = .false.
 energy_verified_no = 0
 allocate(index_particle_star(npart),index_bound(npart),temp_particles(npart))
 open(unit=10,file="particle_index_clean")
 ! Use the sorted array information and check the energy condition first
 do i=1,npart
    !if energy is less than 0, we have bound system. We can accept these particles.
    if (tot_eng_npart(i) < 0. .and. ke_npart(i) < 0.5*abs(pe_npart(i))) then
       write(10,*) i,temp_npart(i),pos_npart(i),sorted_index_npart(i)
       energy_verified_no = energy_verified_no + 1
       ! Save the index of these particles
       ! this is because sometimes even if a particle is farther it could be could but the one before could be unbound
       last_particle_with_neg_e = i
       index_particle_star(energy_verified_no) = sorted_index_npart(i)
       index_bound(energy_verified_no) = i
       temp_particles(energy_verified_no) = temp_npart(i)
       print*,"YES BOUND",i,"i"
    endif
 enddo
 close(10)
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
             print*, bound_index(bound_particles_no),"bound_index(bound_particles_no)"
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
 integer,intent(in) :: energy_verified_no
 integer,intent(out):: number_per_bin
 integer :: number_bins

 !calculate the number of particles per bin
 number_bins = 500
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
subroutine no_per_bin(j,count_particles,double_the_no,number_per_bin,big_bins_no,energy_verified_no,pos_mag_next,rad_inner,double_count)
 integer,intent(inout) :: number_per_bin,double_count
 logical,intent(inout) :: double_the_no
 integer,intent(in)    :: count_particles,big_bins_no,j,energy_verified_no
 real,intent(in)       :: pos_mag_next,rad_inner
 real,parameter :: min_no=5
 integer :: i
 real :: avg_val,diff_val

 avg_val = (pos_mag_next+rad_inner)/2
 diff_val = (pos_mag_next-rad_inner)

 open(15,file="rad_to_bin",status='old',action='write',iostat=i)
 if (i /= 0) then
    ! File does not exist, create it
    open(unit=15,file="rad_to_bin",status='new',action='write',iostat=i)
 endif

 if (j==1) then
    number_per_bin = 1
    double_count = 1
 elseif (double_the_no==.True. .and. count_particles==1) then
    double_count = double_count*2
    number_per_bin = double_count
    if (number_per_bin >= big_bins_no) then
       number_per_bin = big_bins_no
       double_the_no = .False.
    endif
 else
    if  (double_the_no == .False. .and. j  /=  count_particles) then
       if (100*(pos_mag_next-rad_inner)/rad_inner > 30) then
          !print*,(((pos_mag_next-rad_inner)/rad_inner)*100),"per inc",j,"j",count_particles,"count_particles"
          write(15,*) pos_mag_next,rad_inner,j,number_per_bin
          number_per_bin=count_particles
          !if (number_per_bin < min_no) then
          !   number_per_bin = min_no
          ! endif
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
subroutine radius_of_remnant(bound_index,count_particles,number_per_bin,i,energy_verified_no,pos_npart,radius_star,pos_vec_npart,rad_cyl)
 integer,intent(in)    :: count_particles,number_per_bin,i,energy_verified_no,bound_index(:)
 real,intent(in)       :: pos_npart(:),pos_vec_npart(:,:)
 real,intent(out)      :: radius_star,rad_cyl

 real    :: pos_mag_next,pos_mag
 integer :: index_val_next,index_val
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
! print*,norm2(pos_cyl_vec),"mag of pos_cyl vector",pos_mag,"pos_mag"

end subroutine radius_of_remnant
!----------------------------------------------------------------
!+
!  This subroutine calculates the moment of inertia of each particle
!+
!----------------------------------------------------------------
subroutine moment_of_inertia(pos,pos_mag,pmass,i_matrix)
 real,intent(in)  :: pos(3),pos_mag,pmass
 real,intent(out) :: i_matrix(3,3)

 real ::delta(3,3),matrix1(3,1),matrix2(1,3),result_matrix(3,3)

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

 use units ,          only : udist,umass,unit_velocity,unit_energ
 use vectorutils,     only : cross_product3D
 use part,            only : rhoh,poten
 use sortutils,       only : set_r2func_origin,indexxfunc,r2func_origin
 use eos,             only : equationofstate,entropy,X_in,Z_in,gmw,init_eos
 use physcon,         only : gg

 integer,intent(in)               :: npart,iorder(:),numfile
 real,intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)                  :: pmass
 real,intent(inout)               :: xpos(:),vpos(:)
 character(len=20),intent(in)     :: comp_label(:)
 real,intent(in)                  :: interpolate_comp(:,:)
 integer,intent(in)               :: columns_compo
 real,allocatable,intent(out)     :: temp_npart(:),den_npart(:),pos_npart(:),vel_npart(:),pos_wrt_bh(:,:),vel_wrt_bh(:,:),h_npart(:)
 real,allocatable,intent(out)     :: pos_vec_npart(:,:),vel_vec_npart(:,:),tot_eng_npart(:)
 real,allocatable,intent(out)     :: ke_npart(:),pe_npart(:),interp_comp_npart(:,:)
 integer,allocatable,intent(out)  :: sorted_index_npart(:)

 integer             :: i,j,ierr,ieos
 real                :: pos(3),vel(3)
 real                :: potential_i, kinetic_i,energy_i,pos_mag,vel_mag
 real                :: density_i,temperature_i,eni_input,u_i
 real                :: ponrhoi,spsoundi,mu
 real,allocatable    :: composition_i(:)
 real,allocatable    :: A_array(:), Z_array(:)

 ieos = 2
 gmw = 0.61
 call init_eos(ieos,ierr)
 allocate(composition_i(columns_compo))
 call assign_atomic_mass_and_number(comp_label,A_array,Z_array)
 ! Allocate arrays to save the sorted index,density,temperature,radius,total energy of particle wrt centre, velocity_npart
 allocate(temp_npart(npart),den_npart(npart),pos_npart(npart),vel_npart(npart),sorted_index_npart(npart),tot_eng_npart(npart),ke_npart(npart),pe_npart(npart))
 allocate(pos_vec_npart(3,npart),vel_vec_npart(3,npart),pos_wrt_bh(3,npart),vel_wrt_bh(3,npart),h_npart(npart),interp_comp_npart(columns_compo,npart))

 do j = 1, npart
    !Access the rank of each particle in radius and save the sorted index
    i  = iorder(j)
    sorted_index_npart(j) = i

    !if (columns_compo /= 0) then
    !  composition_i(:) = interpolate_comp(:,i)
    !endif

    !the position of the particle is calculated by subtracting the point of
    !highest density.
    !xyzh is position wrt the black hole present at origin.
    call particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos,vel,i,pos_mag,vel_mag)
    !calculate the position which is the location of the particle.
    potential_i  = poten(i)
    kinetic_i    = 0.5*pmass*vel_mag**2
    density_i    = rhoh(xyzh(4,i),pmass)
    energy_i     = potential_i + kinetic_i + vxyzu(4,i)*pmass
    print*,potential_i,"POTENTIAL I",kinetic_i,"Kinetic I"

    ! composition
    if (columns_compo /= 0) then
       composition_i(:) = interpolate_comp(:,i)
    endif
    if (i == 13) then
       print*,composition_i(:),"compo",i,"i before",j,"j"
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

 filename = 'tde.comp'
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
    open(12,file=filename)
    ierr = 0
    !get column labels and send them back.
    read(12, '(a)',iostat=ierr) line
    allocate(comp_label(columns_compo))
    call get_column_labels(line,n_labels,comp_label)
    close(12)
    print*,"comp_label ",comp_label

    open(13,file=filename)
    call skip_header(13,nheader,ierr)
    do k = 1, n_rows
       read(13,*,iostat=ierr) interpolate_comp(:,k)
    enddo
    close(13)
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

 character(len=20),intent(in)   :: comp_label(:)
 character(len=20), allocatable :: new_comp_label(:)
 real,allocatable               :: A_array(:), Z_array(:)
 integer                        :: size_to_allocate, i

 if ( any( comp_label=="nt1" ) ) then
    size_to_allocate = size(comp_label(:))-1

 else
    size_to_allocate = size(comp_label(:))
 endif

 allocate(A_array(size_to_allocate), Z_array(size_to_allocate),new_comp_label(size_to_allocate))
 !we skip nt1 as its neutron.
 new_comp_label = pack(comp_label,comp_label/="nt1")

 do i = 1, size_to_allocate
    if (new_comp_label(i)=="h1") then
       A_array(i) = 1
       Z_array(i) = 1
    endif

    if (new_comp_label(i)=="he3") then
       A_array(i) = 3
       Z_array(i) = 2
    endif

    if (new_comp_label(i)=="he4") then
       A_array(i) = 4
       Z_array(i) = 2
    endif

    if (new_comp_label(i)=="c12") then
       A_array(i) = 12
       Z_array(i) = 6
    endif

    if (new_comp_label(i)=="n14") then
       A_array(i) = 14
       Z_array(i) = 7
    endif

    if (new_comp_label(i)=="o16") then
       A_array(i) = 16
       Z_array(i) = 8
    endif

    if (new_comp_label(i)=="ne20") then
       A_array(i) = 20
       Z_array(i) = 10
    endif

    if (new_comp_label(i)=="mg24") then
       A_array(i) = 24
       Z_array(i) = 12
    endif
    if (new_comp_label(i)=="si28") then
       A_array(i) = 28
       Z_array(i) = 14
    endif

    if (new_comp_label(i)=="s32") then
       A_array(i) = 32
       Z_array(i) = 16
    endif

    if (new_comp_label(i)=="ar36") then
       A_array(i) = 36
       Z_array(i) = 18
    endif

    if (new_comp_label(i)=="ca40") then
       A_array(i) = 40
       Z_array(i) = 20
    endif

    if (new_comp_label(i)=="ti44") then
       A_array(i) = 44
       Z_array(i) = 22
    endif

    if (new_comp_label(i)=="cr48") then
       A_array(i) = 48
       Z_array(i) = 24
    endif

    if (new_comp_label(i)=="fe52") then
       A_array(i) = 52
       Z_array(i) = 26
    endif

    if (new_comp_label(i)=="fe54") then
       A_array(i) = 54
       Z_array(i) = 26
    endif

    if (new_comp_label(i)=="ni56") then
       A_array(i) = 56
       Z_array(i) = 28
    endif

 enddo
 print*, "A and Z arrays assigned"

end subroutine assign_atomic_mass_and_number
!----------------------------------------------------------------
!+
!  This routine is for calculating the gmw value by using
!  1/mu_i = Summation_a ((mass fraction)_a/A_a )*(1+Z_a)
!+
!----------------------------------------------------------------
subroutine calculate_mu(A_array,Z_array,composition_i,columns_compo,mu)

 real,allocatable,intent(in) :: A_array(:), Z_array(:), composition_i(:)
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

 use units , only: udist,umass,unit_velocity,utime,unit_energ,unit_density
 real, intent(in) :: vel_at_infinity,density,time,temperature,mass,xpos(3),rad,distance,pos_mag_star,vel_mag_star,tot_energy,kinetic_energy,potential_energy
 integer, intent(in) :: fileno
 integer :: status, file_id,iostat
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
 write(file_id,'(i5,1x,16(e18.10,1x))')fileno,density*unit_density,temperature,mass*umass,xpos(1)*udist,xpos(2)*udist,xpos(3)*udist,rad*udist,distance*udist,pos_mag_star*udist,&
                      vel_mag_star*unit_velocity,tot_energy,kinetic_energy,potential_energy,time*utime,vel_at_infinity*1e-5,(mass*umass)/(time/(365*24*3600)*utime)
 close(file_id)

end subroutine write_dump_info


 !----------------------------------------------------------------
 !+
 !  This subroutine can write a file with composition of particles wrt black hole
 !
 !+
 !----------------------------------------------------------------
subroutine write_compo_wrt_bh(xyzh,vxyzu,xpos,vpos,pmass,npart,iorder,array_bh_j,interpolate_comp,columns_compo,comp_label,energy_verified_no,last_particle_with_neg_e)
 use units , only: udist

 real,intent(in)    :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)    :: xpos(3),vpos(3),pmass
 integer,intent(in) :: npart,iorder(:),columns_compo
 integer,allocatable,intent(in) :: array_bh_j(:)
 integer,intent(in) :: energy_verified_no,last_particle_with_neg_e
 character(len=20),intent(in) :: comp_label(:)
 real,intent(in)    :: interpolate_comp(:,:)

 integer,allocatable :: array_particle_j(:)
 real,allocatable    :: composition_i(:)
 integer             :: i,j
 real                :: pos_to_bh
 character(len=120)  :: output

 !call particles_bound_to_star(xpos,vpos,xyzh,vxyzu,pmass,npart,iorder,energy_verified_no,last_particle_with_neg_e,array_particle_j,array_bh_j)
 !call composition_array(interpolate_comp,columns_compo,comp_label)
 write(output,"(a8)") 'compo_bh'
 open(4,file=output)
 write(4,"(19(a22,1x))") &
          "posToBH",      &
          comp_label

 allocate(composition_i(columns_compo))
 do j = 1, size(array_bh_j)
    i  = iorder(j) !Access the rank of each particle in radius.
    pos_to_bh = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    if (columns_compo /= 0) then
       composition_i(:) = interpolate_comp(:,i)
    endif
    write(4,'(19(e18.10,1x))') &
              pos_to_bh*udist,&
              composition_i(:)
 enddo
 close(4)

end subroutine write_compo_wrt_bh

 !----------------------------------------------------------------
 !+
 !  This subroutine is to get the temperature cut
 !
 !+
 !----------------------------------------------------------------
subroutine calculate_temp_cut(temperature_array,count_bound,temp_cut,max_temp,temp_found,count_loops_temp,density_array)
 real,intent(in) :: temperature_array(:),max_temp,density_array(:)
 integer,intent(in) :: count_bound,count_loops_temp
 real,intent(out)   :: temp_cut
 integer :: i,count_possible_temp,m
 integer,parameter :: nbins=20000
 real, dimension(nbins)::temp_array_test
 real,allocatable :: avg_density(:)
 real,allocatable :: temp_array_new(:),count_particles_temp(:),diff_count_particles(:),diff2_count_particles(:),diff3_count_particles(:),array_input(:)
 real :: temp_start,count_temp_particles=0,dtemp
 integer :: index_val,avg_inde
 real :: mean,variance,std,cut_off
 real :: count_cut,count_cut_index,lower_limit,upper_limit
 logical, intent(inout) :: temp_found


 ! First we create an array of possible temperature from max_temp to 0 with a step size of 100.
 temp_start = 0.
 dtemp = 100.

 count_cut_index = 0
 count_cut = 0.
 count_possible_temp=1+(max_temp/dtemp)

 ! Create array with the temperatures ranging from 0 to max_temp
 do m=1,nbins
    if (temp_start <= max_temp) then
       temp_array_test(m) = temp_start
       temp_start = temp_start + dtemp
    endif
 enddo

 ! Allocate arrays to save the number of particles per bin
 allocate(temp_array_new(count_possible_temp),count_particles_temp(count_possible_temp), array_input(count_possible_temp),avg_density(count_possible_temp))

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

 print*,"***-------------------------------------"
 print*,temp_array_new,"TEMP ARRAY",size(temp_array_new)
 print*,count_particles_temp,"COUNT PARTICLES TEMP",size(count_particles_temp)
 print*,avg_density,"AVG DENSITY FOR EACH BIN"
 print*,"***-------------------------------------"
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
 print*,count_cut,"count cut first",count_cut_index,"count_cut_index"
 ! this starts from the cound_cut_index found earlier but then tries to make sure that the cut is done when the gaussian bins
 ! have less than 5% particles compared to the max_temp_cut found above
 do i=count_cut_index,size(count_particles_temp)
    if ((count_particles_temp(i)/count_cut)*100 < 1.) then
       count_cut = count_particles_temp(i)
       print*,count_cut,"count_cut",(count_particles_temp(i)/count_cut)*100,"(count_particles_temp(i)/count_cut)*100"
       count_cut_index = i
       exit
    endif
 enddo

 !print*,count_cut_index,"final cut index"

 ! Define the temperature to cut the model at
 temp_cut = temp_array_new(count_cut_index)

 if (temp_cut  /=   max_temp) then
    temp_found = .true.
 endif

 ! If we get the temp_cut as 0. K and the count_loops_temp is 1, then we accept that as a true value
 if (temp_cut == 0.0 .and. count_loops_temp /= 1) then
    temp_found = .false.
 endif
 print*,temp_cut,"TEMP CUT"
end subroutine calculate_temp_cut

! --------------------------------------------------------------------
!  This subroutine calculates the mean, variance and standard deviation
!
! --------------------------------------------------------------------
subroutine statistics(array_data,mean,variance,std)
 real,allocatable,intent(in) :: array_data(:)
 real,intent(out) :: mean,variance
 integer :: size_array,i
 real :: var,sum_val,std

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
