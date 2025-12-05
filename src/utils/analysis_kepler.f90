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
! :Owner: Daniel Price
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
 integer :: i,j,ncomp,ngrid
 real, allocatable :: density(:),rad_grid(:),mass_enclosed(:),bin_mass(:),temperature(:)
 real, allocatable :: rad_vel(:),angular_vel_3D(:,:),comp_kepler(:,:)
 character(len=20), allocatable :: comp_label(:)
 character(len=120) :: output

 ngrid = 0
 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif

 call phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,density,rad_grid,mass_enclosed,bin_mass,&
                              temperature,rad_vel,angular_vel_3D,comp_kepler,comp_label,ncomp,ngrid,numfile)
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
          'cell radial vel',                         &
          'angular vel (x)',                         &  ! ang velocity x component
          'angular vel (y)',                         &  ! ang velocity y component
          'angular vel (z)',                         &  ! velocity z component
          comp_label                                    ! chemical composition
 print*,' Composition array shape: ',shape(comp_kepler)
 do i = 1, ngrid
    write(iunit,'(50(es18.10,1x))')                       &
             real(i),                                    &
             bin_mass(i)*umass,                          &
             mass_enclosed(i)*umass,                     &
             rad_grid(i)*udist,                          &
             density(i)*unit_density,                    &
             temperature(i),                             &
             rad_vel(i)*unit_velocity,                   &
             (angular_vel_3D(j,i)/utime, j=1,3),         &
             (comp_kepler(j,i), j=1,ncomp)
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
                                   temperature,rad_vel,angular_vel_3D,comp_kepler,comp_label,ncomp,ibin,numfile)
 use units,        only:udist,umass,utime,unit_density
 use vectorutils,  only:cross_product3D
 use part,         only:rhoh
 use centreofmass, only:get_centreofmass
 use sortutils,    only:sort_by_radius
 use eos,          only:equationofstate,gmw,init_eos
 use physcon,      only:kb_on_mh,kboltz,atomic_mass_unit,avogadro,gg,pi,pc,years
 use linalg  ,     only:inverse
 integer, intent(in)               :: npart,numfile
 integer, intent(out)              :: ibin,ncomp
 real, intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real, intent(in)                  :: pmass,time
 real, allocatable, intent(out)     :: rad_grid(:),density(:),mass_enclosed(:),bin_mass(:),temperature(:)
 real, allocatable, intent(out)     :: comp_kepler(:,:),rad_vel(:),angular_vel_3D(:,:)
 character(len=20), allocatable, intent(out) :: comp_label(:)
 integer :: i,j,location,ieos,ierr,last_bound,nbound,dummy_bins,nper_bin,n_in_bin,n_big,ntot_bin
 integer :: index_sort,bin_mult,iu_comp,iu_loop,iu_vphi,iu_vbra,iu_rot,iu_remnant
 integer :: iorder(npart)
 integer, allocatable :: isort(:),bound_index(:),sorted_index(:)
 real :: den_all(npart),xpos(3),vpos(3),pos_mag_next
 real :: density_i,density_sum,rad_inner,radius_star,omega_particle
 real :: temperature_i,temperature_sum,rad_vel_i,momentum_i,rad_mom_sum,bhmass
 real :: i_matrix(3,3),I_sum(3,3),Li(3),L_i(3),L_sum(3),inverse_of_i(3,3)
 real :: L_reshape(3,1),matrix_result(3,1),omega(3),pos_i,vel_i,pos_vec_i(3),vel_vec_i(3)
 real :: ke_i,pe_i,tot_e_sum,tot_rem_mass,pos_com(3),vel_com(3),pos_com_mag,vel_com_mag
 real :: vphi_i,R_mag_i,vphi_sum,R_vec(2),omega_vec(3),rad_cyl,breakup
 real :: ke_star,u_star,total_star,distance_from_bh,vel_at_infinity
 real, allocatable :: comp_interp(:,:),comp_i(:),comp_sum(:)
 real, allocatable :: temp(:),den(:),r(:),v(:)
 real, allocatable :: rvec(:,:),vvec(:,:)
 real, allocatable :: h(:),etot(:),ke(:),pe(:)
 real, allocatable :: pos_wrt_bh(:,:),vel_wrt_bh(:,:),comp(:,:)
 logical :: double_bin
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
 call composition_array(comp_interp,ncomp,comp_label)

 ! Call the following function to obtain important information about each particle
 call calc_particles(npart,iorder,numfile,xyzh,vxyzu,pmass,xpos,vpos,comp_label,&
                     comp_interp,ncomp,temp,den,r,v,rvec,vvec,etot,isort,ke,pe,&
                     pos_wrt_bh,vel_wrt_bh,h,comp)

 ! This determines the particles bound to the star. Also removes the streams from the data
 call find_bound(r,temp,etot,npart,isort,bound_index,sorted_index,nbound,&
                 last_bound,ke,pe,den)

 ! This determines how many particles would be added to the bins. We set a min no of bins as 500 in the code and using that along with the particles
 ! that we consider as being part of the remnant implies that the model returns what the max number of particles I would have to add to a bin to get 500 bins
 call calc_nbin(nbound,nper_bin)
 n_big = nper_bin
 bin_mult = 1
 ntot_bin = 0
 ibin = 1
 double_bin = .True.
 print*,' Bound particles: ',size(bound_index),' (energy verified: ',nbound,')'
 ! Define bins arrays
 allocate(density(dummy_bins),rad_grid(dummy_bins),mass_enclosed(dummy_bins),bin_mass(dummy_bins))
 allocate(temperature(dummy_bins),rad_vel(dummy_bins),angular_vel_3D(3,dummy_bins))
 allocate(comp_i(ncomp),comp_sum(ncomp),comp_kepler(ncomp,dummy_bins))

 density_sum     = 0.
 temperature_sum = 0.
 rad_mom_sum     = 0.
 L_sum(:)        = 0.
 I_sum(:,:)      = 0.
 n_in_bin = 0
 comp_sum(:) = 0.
 comp_i(:)   = 0.
 pos_com(:) = 0.
 vel_com(:) = 0.
 tot_e_sum = 0.
 vphi_sum = 9.

 ! write a comp file that includes information about the remnant only
 write(output,"(a4,i5.5)") 'compo',numfile
 open(newunit=iu_comp,file=output)
 write(iu_comp,"(32(a22,1x))") 'i','ibin','radius','x','y','z','radial_vel','temp','density',&
          comp_label,'omega','breakup','j','index_sort'

 ! this will determine when sorted indices are part of the star. We would also need the normal i indices of the sorted particles
 ! Using this we can determine which sorted particles are part of the array and then use the sorted information to calculate all the
 ! quantities we require for the project
 print*,' Processing ',size(bound_index),' bound particles (',nbound,' energy verified)'
 open(newunit=iu_loop,file='big_loop_clean.txt')

 write(output,"(a4,i5.5)") 'vphi',numfile
 open(newunit=iu_vphi,file=output)
 write(iu_vphi,"(2(a22,1x))") 'rad','vphi'

 write(output,"(a4,i5.5)") 'vbra',numfile
 open(newunit=iu_vbra,file=output)
 write(iu_vbra,"(2(a22,1x))") 'rad','vbreak'

 write(output,"(a4,i5.5)") 'rot_i',numfile
 open(newunit=iu_rot,file=output)

 open(newunit=iu_remnant,file='remnant')
 write(iu_remnant,"(6(a22,1x))") 'x','y','z','m','h','rho'

 do i = 1,nbound
    ! get the sorted index of bound particles
    j = bound_index(i)
    index_sort = sorted_index(i)
    n_in_bin = n_in_bin + 1
    ! Calculate the position and velocity VEC of COM
    pos_com(:) = pos_com(:) + xyzh(1:3,index_sort)*pmass
    vel_com(:) = vel_com(:) + vxyzu(1:3,index_sort)*pmass

    ! Obtain the values for each particle that is bound/ part of remnant
    density_i     = den(j)
    temperature_i = temp(j)
    pos_i         = r(j)
    vel_i         = v(j)
    pos_vec_i(:)  = rvec(:,j)
    vel_vec_i(:)  = vvec(:,j)
    R_vec(:)      = pos_vec_i(1:2)
    R_mag_i       = norm2(R_vec)
    ke_i          = ke(j)
    pe_i          = pe(j)
    write(iu_loop,*) i,j,pos_i,vel_i,pos_vec_i(:)*udist,temperature_i,density_i*unit_density,sorted_index(i)

    ! Calculate the angular velocity in cylindrical coordinates
    vphi_i = vel_vec_i(1)*(-pos_vec_i(2)/R_mag_i) + vel_vec_i(2)*(pos_vec_i(1)/R_mag_i)
    vphi_i = vphi_i/R_mag_i

    ! Position magnitude of the next bound particle
    if (i  /=  nbound) pos_mag_next = r(j+1)

    ! composition
    if (ncomp /= 0) comp_i(:) = comp(:,j)

    if (index_sort == 13) then
       print*,' Debug: particle 13 composition = ',comp_i(:),' (j=',j,', i=',i,', index_sort=',index_sort,')'
    endif

    ! Calculate extra quantities
    momentum_i = 0.0
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
    write(iu_comp,'(i9,1x,i5,1x,27(e18.10,1x),1x,i10,1x,i10)') &
              i, ibin, pos_i*udist, pos_vec_i(:)*udist, &
              rad_vel_i, temperature_i, density_i*unit_density, &
              comp_i(:), omega_particle/utime, breakup, j, index_sort

    write(iu_remnant,'(6(e18.10,1x))') pos_vec_i(:), pmass, h(j), density_i

    ! Count particles keeps track of particles in a bin.
    ! Rad_inner is the radius of the first particle that is added to a bin
    if (n_in_bin == 1) rad_inner = pos_i

    ! calculate how many particles will go in a bin
    call no_per_bin(i,n_in_bin,double_bin,nper_bin,n_big,&
                    nbound,pos_mag_next,rad_inner,bin_mult)

    ! we sum the quantities we want to save for the particles
    density_sum        = density_sum + density_i
    temperature_sum    = temperature_sum + temperature_i
    rad_mom_sum        = rad_mom_sum + momentum_i
    L_sum(:)           = L_sum(:) + L_i(:)
    I_sum(:,:)         = I_sum(:,:) + i_matrix(:,:)
    comp_sum(:) = comp_sum(:) + comp_i(:)
    tot_e_sum          = ke_i + pe_i + tot_e_sum
    vphi_sum           = vphi_sum + vphi_i

    ! we check id the n_in_bin is the same as the nper_bin
    ! If true then we save the bin information
    if (n_in_bin==nper_bin) then
       ! total particles binned. Should be the same as nbound at the end
       ntot_bin = ntot_bin+n_in_bin

       ! calculate the bin quantities
       call calc_rbin(bound_index,n_in_bin,nper_bin,i,nbound,&
        r,radius_star,rvec,rad_cyl)

       rad_grid(ibin)      = radius_star
       density(ibin)       = density_sum/n_in_bin
       mass_enclosed(ibin) = ntot_bin*pmass
       bin_mass(ibin)      = n_in_bin*pmass
       ! Change the temperature of particles if its < 1.e3 to 1.e3
       if (temperature_sum < 1.e3) then
          print*,' WARNING! bin ',ibin,' has temperature < 1000 K (',temperature_sum,' K)'
       endif
       temperature(ibin)   = max(temperature_sum/n_in_bin,1e3)
       rad_vel(ibin)       = rad_mom_sum/bin_mass(ibin) ! radial vel of each bin is summation(vel_rad_i*m_i)/summation(m_i)
       if (n_in_bin == 1) then
          if (rad_grid(ibin) < tiny(0.)) then
             angular_vel_3D(:,ibin) = L_sum(:)
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
       comp_kepler(:,ibin) = comp_sum(:)/n_in_bin
       breakup = ((gg*mass_enclosed(ibin)*umass)/(rad_grid(ibin)*udist)**3)**(0.5)
       if (norm2(angular_vel_3D(:,ibin)) > 0) then
          write(iu_vphi,*) udist*rad_grid(ibin),norm2(angular_vel_3D(:,ibin))/utime
          write(iu_vbra,*) udist*rad_grid(ibin),breakup
       endif

       ! Reset the sum values
       n_in_bin        = 0
       density_sum     = 0.
       temperature_sum = 0.
       rad_mom_sum     = 0.
       L_sum(:)        = 0.
       I_sum(:,:)      = 0.
       comp_sum(:)     = 0.
       vphi_sum        = 0.
       ibin            = ibin+1
       nper_bin        = n_big
    endif
 enddo
 close(iu_vbra)
 close(iu_comp)
 close(iu_loop)
 close(iu_vphi)
 close(iu_rot)
 close(iu_remnant)

 ! We want to set the radial and angular velocity of the first bin as the same as the second bin
 angular_vel_3D(:,1) = angular_vel_3D(:,2)
 rad_vel(1) = rad_vel(2)
 ibin = ibin-1
 tot_rem_mass = mass_enclosed(ibin)

 ! Get the centre of mass position and velocity and their magnitudes
 call calc_com(vel_com,pos_com,pos_com_mag,vel_com_mag,tot_rem_mass)
 print*,' Radius of last particle: ',pos_i,' (code units)'
 print*,' Centre of mass position: ',pos_com
 print*,' Centre of mass velocity: ',vel_com
 print*,' Max density position: ',xpos
 print*,' Max density velocity: ',vpos
 print*,' Remnant mass: ',mass_enclosed(ibin)*umass

 ! calculate the energy for the COM and determine if its bound or unbound
 call check_bound(vel_com,pos_com,pos_com_mag,vel_com_mag,bhmass,tot_rem_mass,pmass,&
                  total_star,ke_star,u_star,vel_at_infinity)

 print*,' total energy: ',tot_e_sum

 ! write the dump info file
 call write_dump_info(numfile,density(1),temperature(1),mass_enclosed(ibin),xpos,rad_grid(ibin),distance_from_bh,&
                      pos_com_mag,vel_com_mag,total_star,ke_star,u_star,time,vel_at_infinity)

end subroutine phantom_to_kepler_arrays

!----------------------------------------------------------------
!+
!  This subroutine returns the magnitude of the COM pos and vel
!+
!----------------------------------------------------------------
subroutine calc_com(vel_com,pos_com,pos_com_mag,vel_com_mag,tot_rem_mass)
 real, intent(inout),dimension(3) :: vel_com,pos_com
 real, intent(in)  :: tot_rem_mass
 real, intent(out) :: vel_com_mag,pos_com_mag

 ! Divide the pos_com and vel_com with the total mass enclosed
 pos_com(:) = pos_com(:)/tot_rem_mass
 vel_com(:) = vel_com(:)/tot_rem_mass

 pos_com_mag = norm2(pos_com)
 vel_com_mag = norm2(vel_com)

end subroutine calc_com

!----------------------------------------------------------------
!+
!  This subroutine returns if remnant is bound or unbound
!+
!----------------------------------------------------------------
subroutine check_bound(vel_com,pos_com,pos_com_mag,vel_com_mag,bhmass,tot_rem_mass,pmass,&
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
    call calc_orbit(rem_mass,bhmass_cgs,pos_com_cgs,vel_com_cgs,period_val)
    ar = -gg*0.5*(bhmass_cgs + rem_mass)/tot_energy_remnant_com
    er = 1 - (56.77892268*udist)/ar
    print*,' Semi-major axis: ',ar/au,' AU, eccentricity: ',er
 elseif (tot_energy_remnant_com == 0.) then
    print*,' Parabolic orbit!'
 else
    print*,' Remnant is UNBOUND'
    call calc_vinf(tot_energy_remnant_com,vel_at_infinity)
    print*,' Velocity at infinity: ',vel_at_infinity*1e-5,' km/s'
    ar = gg*0.5*(bhmass_cgs + rem_mass)/tot_energy_remnant_com
    er = 1 + (56.77892268*udist)/ar
    print*,' Semi-major axis: ',ar/au,' AU, eccentricity: ',er
 endif
 print*,' Energy of centre of mass: ',pmass*(0.5*vel_com_mag**2 - (1/pos_com_mag))

end subroutine check_bound

!----------------------------------------------------------------
!+
!  This subroutine returns the vel infinity for the remnant
!  if its unbound
!+
!----------------------------------------------------------------
subroutine calc_orbit(rem_mass,bhmass_cgs,pos_com,vel_com,period_val)
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

end subroutine calc_orbit

!----------------------------------------------------------------
!+
!  This subroutine returns the orbital properties
!+
!----------------------------------------------------------------
subroutine calc_vinf(tot_energy_remnant_com,vel_at_infinity)
 real, intent(in)  :: tot_energy_remnant_com
 real, intent(out) :: vel_at_infinity

 vel_at_infinity = sqrt(2.*tot_energy_remnant_com)

end subroutine calc_vinf

!----------------------------------------------------------------
!+
!  This subroutine returns the position and velocity of a
!  particle wrt to the centre of star/max density point
!+
!----------------------------------------------------------------
subroutine particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos,vel,i,pos_mag,vel_mag)
 real,    intent(in)  :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)  :: xpos(3),vpos(3)
 integer, intent(in)  :: i
 real,    intent(out) :: pos(3),vel(3),pos_mag,vel_mag

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
subroutine find_bound(r,temp,etot,npart,isort,bound_index,sorted_index,nbound,&
                      last_bound,ke,pe,den)
 real, intent(in)    :: temp(:),etot(:),ke(:),pe(:),r(:),den(:)
 integer, intent(in) :: isort(:)
 integer, intent(in) :: npart
 integer, allocatable, intent(out) :: bound_index(:),sorted_index(:)
 integer, intent(out) ::  nbound,last_bound
 integer :: nbound_e,i,nloops_temp,iu
 integer, allocatable :: index_particle_star(:),index_bound(:),index_bound_sorted(:),index_bound_new(:)
 real :: max_temp,temp_cut
 real, allocatable :: temp_bound(:),temp_part(:)
 logical :: temp_found,use_tempcut

 ! Implement temp cut would try to remove the streams. But if you only want
 ! to consider what is bound based on energy condition set this parameter to False
 max_temp = 8000.
 use_tempcut = .true.
 nloops_temp = 0
 nbound = 0
 temp_found = .false.
 nbound_e = 0

 allocate(index_particle_star(npart),index_bound(npart),temp_part(npart))

 open(newunit=iu,file="particle_index_clean")
 ! Use the sorted array information and check the energy condition first
 do i=1,npart
    ! if energy is less than 0, we have bound system. We can accept these particles.
    if (etot(i) < 0. .and. ke(i) < 0.5*abs(pe(i))) then
       write(iu,*) i,temp(i),r(i),isort(i)
       nbound_e = nbound_e + 1
       ! Save the index of these particles
       ! this is because sometimes even if a particle is farther it could be could but the one before could be unbound
       last_bound = i
       index_particle_star(nbound_e) = isort(i)
       index_bound(nbound_e) = i
       temp_part(nbound_e) = temp(i)
       print*,' Particle ',i,' is bound'
    endif
 enddo
 close(iu)

 allocate(temp_bound(nbound_e), index_bound_sorted(nbound_e),index_bound_new(nbound_e))

 do i = 1,nbound_e
    temp_bound(i) = temp_part(i)
    ! This is the sorted index
    index_bound_sorted(i) = index_particle_star(i)
    index_bound_new(i) = index_bound(i)
 enddo

 if (use_tempcut) then
    ! next we loop over the bound particles based on energy condition to find the temp_cut
    ! As the models would need ages to evolve and I can not do that due to how slow some models run, we have streams around the remnants
    ! Hence, we bin the temperature particles and try to find the cut in temperature
    ! But using a temperature cut of 8000 K implies that if I use a model that has streams at high temperature (>1e4 K) because the remnant has just formed
    ! then I would not get rid of the correct particles
    ! Hence, we keep looping until the temperature being returned is the same as the max_temp
    nloops_temp = nloops_temp + 1
    call find_tempcut(temp_bound,nbound_e,temp_cut,max_temp,temp_found,nloops_temp,den)
    max_temp = max_temp + 1000

    allocate(bound_index(nbound_e),sorted_index(nbound_e))
    ! use temp_cut to ignore the streams
    do i = 1,nbound_e
       if (temp_bound(i) > temp_cut) then
          nbound = nbound + 1
          ! Save the sorted array indices only
          bound_index(nbound) = index_bound_new(i)
          sorted_index(nbound) = index_bound_sorted(i)
          if (sorted_index(nbound) == 13) print*,' Debug: particle 13 bound_index = ',bound_index(nbound)
       endif
    enddo
 else
    nbound = nbound_e
    allocate(bound_index(nbound_e),sorted_index(nbound_e))
    bound_index(:) = index_bound_new(:)
    sorted_index(:) = index_bound_sorted(:)
 endif

end subroutine find_bound

!----------------------------------------------------------------
!+
!  This subroutine returns number of particles that can be put into
!  big bins
!+
!----------------------------------------------------------------
subroutine calc_nbin(nbound,nper_bin)
 integer, intent(in) :: nbound
 integer, intent(out) :: nper_bin
 integer, parameter :: number_bins = 500

 ! Calculate the number of particles per bin
 nper_bin = (nbound/number_bins)
 if (mod(nbound,number_bins) /= 0) nper_bin = nper_bin + 1

end subroutine calc_nbin

!----------------------------------------------------------------
!+
!  This subroutine returns number of particles for each bin based
!  on some conditions
!+
!----------------------------------------------------------------
subroutine no_per_bin(j,n_in_bin,double_bin,nper_bin,n_big,&
                      nbound,pos_mag_next,rad_inner,bin_mult)
 integer, intent(inout) :: nper_bin,bin_mult
 logical, intent(inout) :: double_bin
 integer, intent(in)    :: n_in_bin,n_big,j,nbound
 real, intent(in)       :: pos_mag_next,rad_inner
 integer :: i,iu

 open(newunit=iu,file="rad_to_bin",status='old',action='write',iostat=i)
 if (i /= 0) then
    ! File does not exist, create it
    open(newunit=iu,file="rad_to_bin",status='new',action='write',iostat=i)
 endif

 if (j==1) then
    nper_bin = 1
    bin_mult = 1
 elseif (double_bin .and. n_in_bin==1) then
    bin_mult = bin_mult*2
    nper_bin = bin_mult
    if (nper_bin >= n_big) then
       nper_bin = n_big
       double_bin = .False.
    endif
 else
    if (.not.double_bin .and. j /= n_in_bin) then
       if (100*(pos_mag_next-rad_inner)/rad_inner > 30) then
          write(iu,*) pos_mag_next,rad_inner,j,nper_bin
          nper_bin=n_in_bin
       endif
    endif
 endif
 if (j==nbound) nper_bin = n_in_bin

end subroutine no_per_bin

!----------------------------------------------------------------
!+
!  This subroutine returns radius of the remnant
!+
!----------------------------------------------------------------
subroutine calc_rbin(bound_index,n_in_bin,nper_bin,i,nbound,&
                     r,radius_star,rvec,rad_cyl)
 integer, intent(in)    :: n_in_bin,nper_bin,i,nbound,bound_index(:)
 real, intent(in)       :: r(:),rvec(:,:)
 real, intent(out)      :: radius_star,rad_cyl
 integer :: index_val_next,index_val
 real :: pos_mag_next,pos_mag,pos_cyl,pos_cyl_next,pos_cyl_vec(3),pos_cyl_vec_next(3)

 index_val = bound_index(i)
 index_val_next = bound_index(i+1)
 pos_mag = r(index_val)
 pos_cyl_vec(:) = rvec(:,index_val)
 pos_cyl = sqrt(pos_cyl_vec(1)**2 + pos_cyl_vec(2)**2)
 if (n_in_bin == nper_bin .and. i  /=  nbound) then
    pos_mag_next = r(index_val_next)
    pos_cyl_vec_next(:) = rvec(:,index_val_next)

    radius_star = (pos_mag+pos_mag_next)/2
    pos_cyl_next = sqrt(pos_cyl_vec_next(1)**2 + pos_cyl_vec_next(2)**2)
    rad_cyl = (pos_cyl + pos_cyl_next)/2
 else
    radius_star = pos_mag
    rad_cyl = pos_cyl
 endif

end subroutine calc_rbin
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
subroutine calc_particles(npart,iorder,numfile,xyzh,vxyzu,pmass,xpos,vpos,comp_label,&
                          comp_interp,ncomp,temp,den,r,v,rvec,vvec,etot,isort,ke,pe,&
                          pos_wrt_bh,vel_wrt_bh,h,comp)
 use part,            only:rhoh,poten
 use eos,             only:equationofstate,gmw,init_eos
 use physcon,         only:gg
 integer, intent(in)               :: npart,numfile
 integer, intent(in)               :: iorder(:)
 real, intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real, intent(in)                  :: pmass
 real, intent(inout)               :: xpos(:),vpos(:)
 character(len=20), intent(in)     :: comp_label(:)
 real, intent(in)                  :: comp_interp(:,:)
 integer, intent(in)               :: ncomp
 real, allocatable, intent(out)     :: temp(:),den(:),r(:),v(:)
 real, allocatable, intent(out)     :: pos_wrt_bh(:,:),vel_wrt_bh(:,:),h(:)
 real, allocatable, intent(out)     :: rvec(:,:),vvec(:,:),etot(:)
 real, allocatable, intent(out)     :: ke(:),pe(:),comp(:,:)
 integer, allocatable, intent(out)  :: isort(:)
 integer             :: i,j,ierr,ieos
 real                :: pos(3),vel(3)
 real                :: potential_i, kinetic_i,energy_i,pos_mag,vel_mag
 real                :: density_i,temperature_i,eni_input,u_i
 real                :: ponrhoi,spsoundi,mu
 real, allocatable   :: comp_i(:)
 real, allocatable   :: A_array(:), Z_array(:)

 ieos = 2
 gmw = 0.61
 call init_eos(ieos,ierr)
 allocate(comp_i(ncomp))
 call get_atomic_data(comp_label,A_array,Z_array)
 ! Allocate arrays to save the sorted index,density,temperature,radius,total energy of particle wrt centre, velocity
 allocate(temp(npart),den(npart),r(npart),v(npart))
 allocate(isort(npart),etot(npart),ke(npart),pe(npart))
 allocate(rvec(3,npart),vvec(3,npart),pos_wrt_bh(3,npart),vel_wrt_bh(3,npart))
 allocate(h(npart),comp(ncomp,npart))

 do j = 1, npart
    ! access the rank of each particle in radius and save the sorted index
    i  = iorder(j)
    isort(j) = i

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
    if (ncomp /= 0) comp_i(:) = comp_interp(:,i)
    if (i == 13) print*,' Debug: particle 13 composition = ',comp_i(:),' (i=',i,', j=',j,')'
    ! calculate mean molecular weight that is required by the eos module using
    ! the mass fractions for each particle.
    ! do not consider neutron which is the first element of the comp_i array.
    call calculate_mu(A_array,Z_array,comp_i,ncomp,mu)
    gmw = 1./mu
    u_i = vxyzu(4,i)
    eni_input = u_i
    call equationofstate(ieos,ponrhoi,spsoundi,density_i,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi=temperature_i,eni=eni_input)
    ! Save the information for each particle that we need
    den(j) = density_i
    temp(j) = temperature_i
    r(j) = pos_mag
    v(j) = vel_mag
    vvec(:,j) = vel(:)
    rvec(:,j) = pos(:)
    etot(j) = energy_i
    ke(j) = kinetic_i
    pe(j) = potential_i
    pos_wrt_bh(:,j) = xyzh(1:3,i)
    vel_wrt_bh(:,j) = vxyzu(1:3,i)
    h(j) = xyzh(4,i)

    comp(:,j) = comp_interp(:,i)
 enddo

end subroutine calc_particles

!----------------------------------------------------------------
!+
!  This routine reads the output file that contains composition
!  and saves it as a composition array that can be passed to the
!  phantom_to_kepler_arrays subroutine.
!+
!----------------------------------------------------------------
subroutine composition_array(comp_interp,ncomp,comp_label)
 ! First read the file with composition and save that into an array.
 use fileutils, only:get_nlines,skip_header,get_column_labels
 real, allocatable, intent(out)           :: comp_interp(:,:)
 character(len=20), allocatable, intent(out) :: comp_label(:)
 integer, intent(out)                     :: ncomp
 integer                                  :: n_cols,n_rows,ierr,k,nheader,n_labels,iu
 character(len=10000)                     :: line
 character(len=120)                       :: filename
 logical                                  :: iexist

 ncomp = 0
 n_rows = 0
 iexist = .false.

 filename = 'tde.comp'
 ! first check if kepler.comp exists.
 ! this file will only be generated if KEPLER file had composition stored in it.
 inquire(file=filename,exist=iexist)
 if (.not.iexist) then
    write(*,"(a)") ' No file with name '//trim(filename)
    write(*,"(a)") ' -> no composition to save'
    allocate(comp_interp(ncomp,n_rows))
    comp_interp(:,:) = 0.
    allocate(comp_label(ncomp))
 else
    write(*,"(a)") ' Reading composition from '//trim(filename)
    n_rows = get_nlines(filename,skip_comments=.true.,n_columns=n_cols,n_headerlines=nheader)
    ncomp = n_cols

    ! save composition read from file.
    allocate(comp_interp(ncomp,n_rows))
    open(newunit=iu,file=filename)
    ierr = 0
    ! get column labels and send them back.
    read(iu, '(a)',iostat=ierr) line
    allocate(comp_label(ncomp))
    call get_column_labels(line,n_labels,comp_label)
    close(iu)
    print*,' Composition labels: ',comp_label

    open(newunit=iu,file=filename)
    call skip_header(iu,nheader,ierr)
    do k = 1, n_rows
       read(iu,*,iostat=ierr) comp_interp(:,k)
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
subroutine get_atomic_data(comp_label,A_array,Z_array)
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

end subroutine get_atomic_data
!----------------------------------------------------------------
!+
!  This routine is for calculating the gmw value by using
!  1/mu_i = Summation_a ((mass fraction)_a/A_a )*(1+Z_a)
!+
!----------------------------------------------------------------
subroutine calculate_mu(A_array,Z_array,comp_i,ncomp,mu)
 real, allocatable, intent(in) :: A_array(:), Z_array(:), comp_i(:)
 integer, intent(in)         :: ncomp
 real,    intent(out)        :: mu
 integer                     :: index_val

 mu = 0.
 if (ncomp /= 0) then
    do index_val = 1, ncomp-1
       mu = mu + (comp_i(index_val+1)*(1+Z_array(index_val)))/A_array(index_val)
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
    open(newunit=file_id,file=filename,status='old',position='append',action='write',iostat=status)
    if (status /= 0) then
       write(*,*) 'Error opening file: ', filename
       stop
    endif
 else
    open(newunit=file_id,file=filename,status='new',action='write',iostat=status)
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
!+
!----------------------------------------------------------------
subroutine find_tempcut(temp_arr,nbound,temp_cut,max_temp,temp_found,&
                        nloops_temp,density_array)
 real, intent(in) :: temp_arr(:),max_temp,density_array(:)
 integer, intent(in) :: nbound,nloops_temp
 real, intent(out)   :: temp_cut
 logical, intent(inout) :: temp_found
 integer, parameter :: nbins = 20000
 integer :: i,npossible_temp,m,icut
 real :: temp_start,ntemp_part,dtemp
 real :: mean,variance,std,cut_off
 real :: ncut,lower_limit,upper_limit
 real, dimension(nbins) ::temp_test
 real, allocatable :: avg_density(:)
 real, allocatable :: temp_bins(:),npart_temp(:)

 ! First we create an array of possible temperature from max_temp to 0 with a step size of 100.
 temp_start = 0.
 dtemp = 100.

 ntemp_part = 0
 icut = 0
 ncut = 0.
 npossible_temp=1+nint(max_temp/dtemp)

 ! Create array with the temperatures ranging from 0 to max_temp
 do m=1,nbins
    if (temp_start <= max_temp) then
       temp_test(m) = temp_start
       temp_start = temp_start + dtemp
    endif
 enddo

 ! Allocate arrays to save the number of particles per bin
 allocate(temp_bins(npossible_temp),npart_temp(npossible_temp))
 allocate(avg_density(npossible_temp))

 npart_temp(:) = 0

 ! create an array of the same size as npossible_temp
 do m=1,npossible_temp
    temp_bins(m) = temp_test(m)
 enddo

 ! this will count the particles for each temperature and then save them into a new array
 do i=1,nbound
    do m=1,npossible_temp-1
       if (temp_arr(i) >= temp_bins(m) .and. temp_arr(i) < temp_bins(m+1) )  then
          ntemp_part = npart_temp(m) + 1
          npart_temp(m) =  ntemp_part
          avg_density(m) = density_array(i)
       endif
    enddo
 enddo

 print*,' Temperature cut analysis:'
 print*,'     Temperature array size: ',size(temp_bins)
 print*,'  Particle count array size: ',size(npart_temp)

 ! calculate the mean, std of the data
 call statistics(npart_temp,mean,variance,std)

 ! using 2 sigma as the data sample is small to determine the outlier
 cut_off = std*2
 lower_limit = mean - cut_off
 upper_limit = mean + cut_off

 ! This loops and find the last element which is outside the limits based on 2 sigma
 do i=1,size(npart_temp)
    if (npart_temp(i) > upper_limit .or. npart_temp(i) < lower_limit) then
       ncut = npart_temp(i)
       icut = i
    endif
 enddo
 print*,' Initial temperature cut: ',ncut,' at index ',icut

 ! this starts from the icut found earlier but then tries to make sure that the cut is done when the gaussian bins
 ! have less than 5% particles compared to the max_temp_cut found above
 do i=icut,size(npart_temp)
    if ((npart_temp(i)/ncut)*100 < 1.) then
       ncut = npart_temp(i)
       print*,' Temperature cut updated: ',ncut,' (',(npart_temp(i)/ncut)*100,'% of max)'
       icut = i
       exit
    endif
 enddo

 ! Define the temperature to cut the model at
 temp_cut = temp_bins(icut)

 if (temp_cut < max_temp .or. temp_cut > max_temp) temp_found = .true.

 ! If we get the temp_cut as 0. K and the nloops_temp is 1, then we accept that as a true value
 if (temp_cut < tiny(0.0) .and. nloops_temp /= 1) temp_found = .false.
 print*,' Final temperature cut: ',temp_cut,' K'

end subroutine find_tempcut

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
