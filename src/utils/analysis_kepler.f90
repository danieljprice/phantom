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
          'cell radial momentum',                    &
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
 !  This subroutine returns the position and velocity of a
 !  particle wrt to the centre of star/max density point
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
 use orbits_data,     only : escape
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
 integer :: dummy_size,dummy_bins,number_per_bin,count_particles,number_bins,no_particles,big_bins_no,tot_binned_particles
 real :: density_i,density_sum,rad_inner,rad_outer,radius_star
 logical :: double_the_no,escape_star
 real :: omega_particle,omega_bin,pos_com(3),vel_com(3),pos_mag_star,vel_mag_star
 real :: eni_input,u_i,temperature_i,temperature_sum,mu
 real :: ponrhoi,spsoundi,rad_vel_i,momentum_i,rad_mom_sum
 real :: bhmass,pos_prev(3),vel_prev(3),pos_mag_prev,vel_mag_prev
 real :: i_matrix(3,3),I_sum(3,3),Li(3),L_i(3),L_sum(3),inverse_of_i(3,3),L_reshape(3,1),matrix_result(3,1),omega(3)
 real,allocatable    :: A_array(:), Z_array(:)
 real,allocatable    :: interpolate_comp(:,:),composition_i(:),composition_sum(:)
 real :: ke_star,u_star,total_star,distance_from_bh,vel_from_bh,vel_at_infinity
 ieos = 2
 call init_eos(ieos,ierr)
 gmw=0.61
 bhmass=1.

 ! performing a loop to determine maximum density particle position
 do j = 1,npart
    den_all(j) = rhoh(xyzh(4,j),pmass)
 enddo
 location = maxloc(den_all,dim=1)
 print*,location,"location of max density"
 ! Determining centre of star as max density particle
 xpos(:) = xyzh(1:3,location)
 vpos(:) = vxyzu(1:3,location)
 distance_from_bh = sqrt(dot_product(xpos(:),xpos(:)))
 vel_from_bh = sqrt(dot_product(vpos(:),vpos(:)))
 print*,"******************"
 print*,distance_from_bh*udist,"distance from bh",vel_from_bh*unit_velocity,"velfrom bh"
 print*,"*******************"
 ! sorting particles
 call set_r2func_origin(xpos(1),xpos(2),xpos(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
 call composition_array(interpolate_comp,columns_compo,comp_label)
 call particles_bound_to_star(xpos,vpos,xyzh,vxyzu,pmass,npart,iorder,energy_verified_no,last_particle_with_neg_e,array_particle_j,array_bh_j,interpolate_comp,columns_compo,comp_label,numfile)
 call assign_atomic_mass_and_number(comp_label,A_array,Z_array)
 print*,array_particle_j(energy_verified_no),"Last particle indes",last_particle_with_neg_e
 print*,energy_verified_no,"energy_verified_no",size(array_particle_j)
 call particles_per_bin(energy_verified_no,number_per_bin)
 tot_binned_particles = 0
 big_bins_no          = number_per_bin
 no_particles         = 1
 dummy_bins           = 5000
 ibin                 = 1
 double_the_no        = .True.

 allocate(density(dummy_bins),rad_grid(dummy_bins),mass_enclosed(dummy_bins),bin_mass(dummy_bins),temperature(dummy_bins),rad_vel(dummy_bins),angular_vel_3D(3,dummy_bins))
 density_sum     = 0.
 temperature_sum = 0.
 rad_mom_sum     = 0.
 L_sum(:)        = 0.
 I_sum(:,:)      = 0.
 count_particles = 0
 allocate(composition_i(columns_compo))
 allocate(composition_sum(columns_compo))
 allocate(composition_kepler(columns_compo,dummy_bins))
 composition_sum(:) = 0.
 composition_i(:)   = 0.

 ! writing files with angular velocity info
 open(1,file="particleOmega.info")
 write(1,*) "[pos]"," ","[omega]"
 open(2,file="binOmega.info")
 write(2,*) "[pos]"," ","[omega]"
 open(3,file="radius_of_bins.info")
 write(3,*) "[ibin]"," ","[rad_inner]"," ","[rad_outer]"," ","[Position rad next]"," ","[particles in bin]"
 write(output,"(a4,i5.5)") 'compo',numfile
 open(4,file=output)
 write(4,"(25(a22,1x))") &
          "i",      &
          "ibin",   &
          "radius", &
          "x", &
          "y", &
          "x", &
          comp_label

 pos_com(:) = 0.
 vel_com(:) = 0.
 ! Now we calculate the different quantities of the particles and bin them
 do j=1,energy_verified_no
    i      = iorder(array_particle_j(j))
    if (j /= energy_verified_no) then
       i_next = iorder(array_particle_j(j+1))
    else
       i_next = iorder(array_particle_j(j))
    endif

    call particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos,vel,i,pos_mag,vel_mag)

    ! calculating centre of mass position and velocity wrt black hole
    pos_com(:) = pos_com(:) + xyzh(1:3,i)*pmass
    vel_com(:) = vel_com(:) + vxyzu(1:3,i)*pmass

    if (j  /=  energy_verified_no) then
       call particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos_next,vel_next,i_next,pos_mag_next,vel_mag_next)
    endif

    count_particles = count_particles + 1
    if (count_particles == 1) then
       rad_inner = pos_mag
       !print*,j,"j","first",rad_inner,"rad_inner",count_particles
    endif

    !print*,pos_mag_next-pos_mag,"difference in pos mag of next and current particle",i,"i",i_next,"i_next",j,"j",j+1,"jnext"
    call  no_per_bin(j,count_particles,double_the_no,number_per_bin,big_bins_no,energy_verified_no,pos_mag_next,rad_inner)
    if (number_per_bin == count_particles) then
       rad_outer = pos_mag
    endif
    ! composition
    if (columns_compo /= 0) then
       composition_i(:) = interpolate_comp(:,i)
    endif
    composition_sum(:) = composition_sum(:) + composition_i(:)

    write(4,'(i9,1x,i5,1x,23(e18.10,1x))') &
               i, &
               ibin, &
               pos_mag*udist, &
               pos(1)*udist, &
               pos(2)*udist, &
               pos(3)*udist, &
               composition_i(:)
    ! calculate mean molecular weight that is required by the eos module using
    ! the mass fractions for each particle.
    ! do not consider neutron which is the first element of the composition_i array.

    call calculate_mu(A_array,Z_array,composition_i,columns_compo,mu)

    gmw = 1./mu
    ! Density
    density_i   = rhoh(xyzh(4,i),pmass)
    density_sum = density_sum + density_i

    ! Temperature
    u_i       = vxyzu(4,i)
    eni_input = u_i
    call equationofstate(ieos,ponrhoi,spsoundi,density_i,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi=temperature_i,eni=eni_input)
    temperature_sum = temperature_sum + temperature_i

    ! Radial momentum
    ! we skip the first particle as its the one that exists at the center of
    ! star and hence will give infinite rad_vel as rad = 0.
    if (pos_mag > 0.) then
       rad_vel_i    = dot_product(vel(:),pos(:))/pos_mag
       momentum_i   = rad_vel_i*pmass
       rad_mom_sum  = rad_mom_sum + momentum_i
    endif

    ! Angular momentum
    call cross_product3D(pos(:),vel(:),Li(:))
    L_i(:)   = Li(:)*pmass
    L_sum(:) = L_sum(:) + L_i(:)
    if (pos_mag == 0.) then
       omega_particle = 0.
    else
       omega_particle = sqrt(dot_product(Li(:)/(pos_mag**2),Li(:)/(pos_mag**2)))
    endif

    write(1,*)pos_mag,omega_particle
    ! Moment of inertia
    call moment_of_inertia(pos,pos_mag,pmass,i_matrix)
    I_sum(:,:) = I_sum(:,:) + i_matrix(:,:)

    if (count_particles==number_per_bin .or. j==energy_verified_no) then
       tot_binned_particles = tot_binned_particles+count_particles
       call radius_of_remnant(array_particle_j,count_particles,number_per_bin,j,energy_verified_no,xpos,vpos,xyzh,vxyzu,iorder,pos_mag,radius_star)
       rad_grid(ibin)      = radius_star
       density(ibin)       = density_sum/count_particles
       mass_enclosed(ibin) = tot_binned_particles*pmass
       bin_mass(ibin)      = count_particles*pmass
       temperature(ibin)   = max(temperature_sum/count_particles,1e3)
       rad_vel(ibin)       = rad_mom_sum/bin_mass(ibin) !Radial vel of each bin is summation(vel_rad_i*m_i)/summation(m_i)
       if (count_particles == 1) then
          if (pos_mag==0.) then
             angular_vel_3D(:,ibin)  = L_sum(:)
          else
             angular_vel_3D(:,ibin) = L_sum(:) / (pos_mag**2*pmass)
          endif
       else
          inverse_of_i  = inverse(I_sum, 3)
          L_reshape     = reshape(L_sum(:),(/3,1/))
          matrix_result = matmul(inverse_of_i,L_reshape)
          omega         = reshape(matrix_result,(/3/))
          angular_vel_3D(:,ibin) = omega
       endif
       omega_bin = sqrt(dot_product(angular_vel_3D(:,ibin),angular_vel_3D(:,ibin)))
       write(2,*)pos_mag,omega_bin
       composition_kepler(:,ibin) = composition_sum(:)/count_particles
       write(3,*) ibin,rad_inner,rad_outer,pos_mag_next,count_particles
       count_particles = 0
       density_sum     = 0.
       temperature_sum = 0.
       rad_mom_sum     = 0.
       L_sum(:)        = 0.
       I_sum(:,:)      = 0.
       composition_sum(:) = 0.
       ibin            = ibin+1
    endif
 enddo
 print*,iu1,"iu",iu2,"iu2"
 close(1)
 close(2)
 close(3)
 close(4)
 ibin = ibin-1
 print*,mass_enclosed(ibin)*umass,"enclodsed mass",pos_com,"pos com"
 print*,rad_grid(ibin),"Radius MAX"
 pos_com(:) = pos_com(:)/mass_enclosed(ibin)
 print*,pos_com,"pos of com"
 vel_com(:) = vel_com(:)/mass_enclosed(ibin)
 print*,pos_com,"pos_com",xpos,"xpos",vel_com,"vel com",vpos,"vpos"
 print*,ibin,"ibin",tot_binned_particles
 pos_mag_star = sqrt(dot_product(pos_com,pos_com))
 vel_mag_star = sqrt(dot_product(vel_com,vel_com))
 print*,"bhmass*umass",bhmass*umass
 print*,"-------------------------------------------------------------------------"
 print*,escape(vel_mag_star*unit_velocity,bhmass*umass,pos_mag_star*udist)
 print*,"-------------------------------------------------------------------------"
 ke_star = 0.5*(vel_mag_star*unit_velocity)**2
 u_star = -gg*(bhmass*umass)/(pos_mag_star*udist)
 print*,"--------------"
 print*,bhmass*umass,"BH mass", pos_mag_star*udist, "Pos Mag star",vel_mag_star*unit_velocity,"Vel Mag Star"
 print*,"--------------"
 total_star = u_star+ke_star
 print*,umass,"umass",gg,"gg"
 print*,total_star,"total_star",u_star,"ustar",ke_star,"ke_star"
 print*,mass_enclosed(ibin),"/mass_enclosed(ibin)",mass_enclosed(ibin)*umass
 vel_at_infinity = sqrt(2.*total_star)
 if (isnan(vel_at_infinity)) then
    vel_at_infinity = 0.
 endif
 print*,vel_at_infinity*1.e-5,"vel at infinity in Km/s"
 print*,umass,"umass",udist,"udist",unit_density,"unit_density",unit_velocity,"unit_velocity",utime,"utime"
 ! write information to the dump_info file
 call write_compo_wrt_bh(xyzh,vxyzu,xpos,vpos,pmass,npart,iorder,array_bh_j,interpolate_comp,columns_compo,comp_label,energy_verified_no,last_particle_with_neg_e)
 call write_dump_info(numfile,density(1),temperature(1),mass_enclosed(ibin),xpos,rad_grid(ibin),distance_from_bh,&
                         pos_mag_star,vel_mag_star,total_star,ke_star,u_star,time,vel_at_infinity)
end subroutine phantom_to_kepler_arrays

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
subroutine particles_bound_to_star(xpos,vpos,xyzh,vxyzu,pmass,npart,iorder,energy_verified_no,last_particle_with_neg_e,array_particle_j,array_bh_j,interpolate_comp,columns_compo,comp_label,numfile)
 use units ,          only : udist,umass,unit_velocity,unit_energ
 use vectorutils,     only : cross_product3D
 use part,            only : rhoh,poten
 use centreofmass,    only : get_centreofmass
 use sortutils,       only : set_r2func_origin,indexxfunc,r2func_origin
 use eos,             only : equationofstate,entropy,X_in,Z_in,gmw,init_eos
 use physcon,         only : gg

 integer,intent(in)               :: npart,iorder(:),numfile
 real,intent(in)                  :: xyzh(:,:),vxyzu(:,:),xpos(:),vpos(:)
 real,intent(in)                  :: pmass
 character(len=20),intent(in)     :: comp_label(:)
 real,intent(in)                  :: interpolate_comp(:,:)
 integer,intent(in)               :: columns_compo
 integer,intent(out)              :: energy_verified_no,last_particle_with_neg_e
 integer,allocatable,intent(out)  :: array_particle_j(:),array_bh_j(:)

 character(len=120)  :: output
 integer,allocatable :: index_particle_star(:),index_particle_bh(:)
 integer :: i,j,dummy_size,index_val,particle_bound_bh,index_val_bh,count_val,count_val_unbound,count_bound_both
 real :: potential_wrt_bh,kinetic_wrt_bh,tot_wrt_bh,pos(3),vel(3)
 real :: potential_i, kinetic_i,energy_i,pos_mag,vel_mag
 logical :: bound_to_bh,bound_to_star
 real,allocatable    :: composition_i(:)

 bound_to_bh = .false.
 bound_to_star = .false.
 particle_bound_bh = 0
 energy_verified_no = 0
 index_val = 1
 index_val_bh = 1
 dummy_size = 1e8
 count_val = 0
 count_val_unbound = 0
 count_bound_both = 0

 write(output,"(a8,i5.5)") 'compfull',numfile
 open(5,file=output)
 write(5,"(18(a22,1x))") &
          comp_label

 allocate(index_particle_star(dummy_size),index_particle_bh(dummy_size))
 allocate(composition_i(columns_compo))
 do j = 1, npart

    i  = iorder(j) !Access the rank of each particle in radius.
    potential_wrt_bh = -(gg*umass*pmass*umass)/(sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))*udist)
    kinetic_wrt_bh = 0.5*pmass*umass*dot_product(vxyzu(1:3,i),vxyzu(1:3,i))*unit_velocity**2
    tot_wrt_bh = potential_wrt_bh+ kinetic_wrt_bh+vxyzu(4,i)*pmass*unit_energ
    if (tot_wrt_bh < 0.) then
       bound_to_bh = .True.
    endif
    if (columns_compo /= 0) then
       composition_i(:) = interpolate_comp(:,i)
    endif
    write(5,'(18(e18.10,1x))') &
       composition_i(:)
    !the position of the particle is calculated by subtracting the point of
    !highest density.
    !xyzh is position wrt the black hole present at origin.
    call particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos,vel,i,pos_mag,vel_mag)

    !calculate the position which is the location of the particle.
    potential_i = poten(i)
    kinetic_i     = 0.5*pmass*vel_mag**2

    energy_i = potential_i + kinetic_i + vxyzu(4,i)*pmass

    !if energy is less than 0, we have bound system. We can accept these particles.
    if (energy_i < 0. .and. kinetic_i < 0.5*abs(potential_i)) then
       bound_to_star = .True.
       energy_verified_no = energy_verified_no + 1
       last_particle_with_neg_e = j
       index_particle_star(index_val) = j
       index_val = index_val+1
    endif

    if (bound_to_bh == .True. .and. bound_to_star == .false.) then
       count_val = count_val + 1
       index_particle_bh(index_val_bh) = j
       particle_bound_bh = particle_bound_bh +1
       index_val_bh = index_val_bh+1
    endif
    if (bound_to_bh == .True. .and. bound_to_star == .True.) then
       count_bound_both = count_bound_both + 1
    endif

    if (bound_to_bh == .false. .and. bound_to_star == .false.) then
       count_val_unbound = count_val_unbound + 1
    endif
    bound_to_bh = .false.
    bound_to_star = .false.
 enddo
 close(5)
 print*,"==================================="
 print*,count_val,"count val",count_val_unbound,"unbound count",count_bound_both,"count_bound_both"
 print*,"===================================="
 !next we save the index of particles which are part of star into a new array
 allocate(array_particle_j(energy_verified_no))
 do i=1,energy_verified_no
    array_particle_j(i) = index_particle_star(i)
 enddo
 print*,"-------"
 print*,npart,"npart in determing particles"
 print*,"-------"
 ! we save the index of the particles bound to the SMBH
 allocate(array_bh_j(particle_bound_bh))
 do i=1,particle_bound_bh
    array_bh_j(i) = index_particle_bh(i)
 enddo

 print*,"--------"
 print*,particle_bound_bh,"particle bound to the bh",size(array_bh_j),"array bh j"
 print*,"Size of array with particles",size(array_particle_j)
 print*,"--------"
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
    number_per_bin = 1 + number_per_bin
 endif

end subroutine particles_per_bin
!----------------------------------------------------------------
!+
!  This subroutine returns number of particles for each bin based
!  on some conditions
!+
!----------------------------------------------------------------
subroutine no_per_bin(j,count_particles,double_the_no,number_per_bin,big_bins_no,energy_verified_no,pos_mag_next,rad_inner)
 integer,intent(inout) :: number_per_bin
 logical,intent(inout) :: double_the_no
 integer,intent(in)    :: count_particles,big_bins_no,j,energy_verified_no
 real,intent(in)       :: pos_mag_next,rad_inner


 if (j==1) then
    number_per_bin = 1
 elseif (double_the_no==.True. .and. count_particles==1) then
    number_per_bin = number_per_bin*2
    if (number_per_bin >= big_bins_no) then
       number_per_bin = big_bins_no
       double_the_no = .False.
    endif
 else
    if (pos_mag_next - rad_inner > 0.1) then
       number_per_bin=count_particles
       if (number_per_bin < 10) then
          number_per_bin = 10
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
subroutine radius_of_remnant(array_particle_j,count_particles,number_per_bin,j,energy_verified_no,xpos,vpos,xyzh,vxyzu,iorder,pos_mag,radius_star)
 integer,intent(in)    :: count_particles,number_per_bin,j,energy_verified_no,iorder(:),array_particle_j(:)
 real,intent(in)       :: xyzh(:,:),vxyzu(:,:),xpos(:),vpos(:),pos_mag
 real,intent(out)      :: radius_star

 real :: pos_mag_next,vel_mag_next,pos_next(3),vel_next(3)
 integer :: i_next

 if (count_particles==number_per_bin .and. j  /=  energy_verified_no) then
    i_next = iorder(array_particle_j(j+1))
    call particle_pos_and_vel_wrt_centre(xpos,vpos,xyzh,vxyzu,pos_next,vel_next,i_next,pos_mag_next,vel_mag_next)
    radius_star = (pos_mag+pos_mag_next)/2
 else
    radius_star = pos_mag
 endif

end subroutine radius_of_remnant
!----------------------------------------------------------------
!+
!  This subroutine calculates the moment of inertia
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
    print*,"comp_label ",comp_label

    open(13, file=filename)
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

 if ( ANY( comp_label=="nt1" ) ) then
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
    open(unit=file_id, file=filename, status='old', position="append", action="write", iostat=status)
    if (status /= 0) then
       write(*,*) 'Error opening file: ', filename
       stop
    endif

 else
    open(unit=file_id, file=filename, status='new', action='write', iostat=status)
    if (status /= 0) then
       write(*,*) 'Error creating file: ', filename
       stop
    endif
    ! Write headers to file
    write(file_id,'(16(a22,1x))') &
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
               "Escape_in"
 endif
 write(file_id,'(i5,1x,15(e18.10,1x))')fileno,density*unit_density,temperature,mass*umass,xpos(1)*udist,xpos(2)*udist,xpos(3)*udist,rad*udist,distance*udist,pos_mag_star*udist,&
                      vel_mag_star*unit_velocity,tot_energy,kinetic_energy,potential_energy,time*utime,vel_at_infinity*1e-5
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

end module analysis
