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
 use units,           only : udist,umass,unit_density,unit_ergg,unit_velocity,utime !units required to convert to kepler units.
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


 !If dumpfile is not a complete dump we don't read it.
 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif


 !allocate for composition_kepler
 !Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile


 !if dumpfile is a full dump, we call the subroutine for getting the arrays we need
 call phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,pressure,rad_grid,mass,&
                                  density,temperature,entropy_array,int_eng,velocity_3D,bin_mass,&
                                  y_e,a_bar,composition_kepler,comp_label,n_comp,ngrid,rad_mom,&
                                  angular_vel_3D,numfile)

 write(output,"(a4,i5.5)") 'ptok',numfile
 write(*,'("Output file name is ",A)') output
print*,ngrid,"NGRID*****"
 !open the output file and save the data in the format kepler likes. Using same labels as kepler.
 open(iunit,file=output)
 write(iunit,'("# ",i5,"   # Version")') 10000
 write(iunit,'("# ",es20.12,"   # Dump file number")') time
 write(iunit,"('#',50(a22,1x))")                     &
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

 !as the units units in GR setup and newtonian setup are not same,
 !we need to consider both units while writing the file.
 print*,udist,'udist',umass,'umass'

 do i = 1, ngrid
    grid = i
    write(iunit,'(50(es18.10,1X))')                       &
              grid,                                       &
              bin_mass(i)*umass,                          &
              mass(i)*umass,                              &
              rad_grid(i)*udist,                          &
              density(i)*unit_density,                    &
              temperature(i),                             &
              rad_mom(i)*unit_velocity*umass,             &
              (angular_vel_3D(j,i)/utime, j=1,3),         &
              (composition_kepler(j,i), j=1,n_comp)
 enddo
 close(iunit)
 print*,'----------------------------------------------------------------------'
 print*,"Important summary"
 print*,rad_grid(ngrid)*udist, 'radius of star',rad_grid(ngrid)
 print*,mass(ngrid)*umass, 'total mass of star',mass(ngrid)
 print*,density(1)*unit_density,'maximum density of star'
end subroutine do_analysis

 !----------------------------------------------------------------
 !+
 !  routine for binning the data as a function of radius.
 !  The arrays generated are used by do_analysis subroutine.
 !+
 !----------------------------------------------------------------
subroutine phantom_to_kepler_arrays(xyzh,vxyzu,pmass,npart,time,pressure,rad_grid,mass,&
                                    density,temperature,entropy_array,int_eng,velocity_3D,bin_mass,&
                                    y_e,a_bar,composition_kepler,comp_label,columns_compo,correct_ngrid,rad_mom,&
                                    angular_vel_3D,numfile)
 use units , only: udist,umass,unit_velocity,utime,unit_energ,unit_density
 use vectorutils,     only : cross_product3D
 use part,            only : rhoh,poten
 use centreofmass,    only : get_centreofmass
 use sortutils,       only : set_r2func_origin,indexxfunc,r2func_origin
 use eos,             only : equationofstate,entropy,X_in,Z_in,gmw,init_eos
 use physcon,         only : kb_on_mh,kboltz,atomic_mass_unit,avogadro,gg,pi,pc,years
 use orbits_data,     only : escape, orbital_parameters
! use linalg,          only : inverse
use linearalgebra , only : inverse
 integer,intent(in)               :: npart,numfile
 integer,intent(out)              :: correct_ngrid
 real,intent(in)                  :: xyzh(:,:),vxyzu(:,:)
 real,intent(in)                  :: pmass,time
 real,intent(out),allocatable     :: rad_grid(:),mass(:),density(:)!rad_grid stores radius
 real,intent(out),allocatable     :: temperature(:),entropy_array(:),int_eng(:),bin_mass(:),rad_mom(:)
 real,intent(out),allocatable     :: pressure(:),y_e(:),a_bar(:),velocity_3D(:,:)
 real,intent(out),allocatable              :: composition_kepler(:,:),angular_vel_3D(:,:)
 character(len=20),allocatable,intent(out) :: comp_label(:)

 integer :: no_in_bin !this stores the number of particles in bin after each loop.
 integer :: ibin
 integer :: iorder(npart),j,i,s,m,last_particle_with_neg_e
 integer :: number_particle,ieos,ierr
 integer :: columns_compo,location
 integer :: ngrid
 integer::  number_bins,number_tot, number_per_bin
 integer :: count=0
 integer :: c_particle,index_i,energy_verified_no
 real :: density_sum,density_i,eni_input
 real :: u_sum,u_i !specific internal energy storage
 real :: omega_sum(3),omega(3)
 real :: moment_of_inertia
 real :: temperature_i,temperature_sum
 real :: rad_velocity,rad_vel_sum,momentum
 real :: pressure_i,pressure_sum
 real :: pos(3),vel(3),rad,rad_next
 real :: xpos(3),vpos(3),star_centre(3) !COM position and velocity
  !COM position and velocity
 real :: ponrhoi,spsoundi,vel_sum(3),Li(3)
 real :: velocity_norm,escape_vel,kinetic_add
 real :: Y_in,mu
 real :: bh_mass,L(3)
 real :: tot_energy
 real :: potential_i,kinetic_i,energy_i,energy_total,angular_momentum_h(3)
 real :: velocity_wrt_bh(3),rad_test,velocity_bh,com_star(3),position_bh
 real,allocatable    :: interpolate_comp(:,:),composition_i(:),composition_sum(:)
 real,allocatable    :: energy_tot(:)
 real,allocatable    :: A_array(:), Z_array(:)
 real :: mass_star
 real :: omega_val(3),val_omega, den_all(npart)
 integer :: count_new,skip_breakup
 real :: tot_thermal_energy, tot_internal_energy
 real ::delta(3,3),matrix1(3,1),matrix2(1,3),result_matrix(3,3),final_val_omega(3)
 real ::i_matrix(3,3),inverse_of_i(3,3),omega_reshape(3,1),matrix_result(3,1)
 real :: radius_last 
 logical                 :: iexist

print*,utime,'time!!','this is analysis_test file!'
 !The star is not on the origin as BH exists at that point.
 !minimum h value corresponds to position of maximum density.
 !COM is not a good option as it does not work for severe disruptione events.
print*,"Before location is determined"
 location = minloc(xyzh(4,:),dim=1)
 print*,"This is where issue happend",location,"location",size(xyzh)
 star_centre(:) = xyzh(1:3,location)
print*,"issue happened"
print*, location, "INITIAL LOCATION BASED ON h"
do j = 1,npart 
  den_all(j) = rhoh(xyzh(4,j),pmass)
enddo
location = maxloc(den_all,dim=1)
star_centre(:) = xyzh(1:3,location)
print*,"**********"
print*,location
print*,maxval(den_all),"MAX DEN VAL IN ARRAY",location,star_centre,"centre"
 !we use the equation number 12 for Newtonian and 2 for GR analysis.
print*,rhoh(xyzh(4,location),pmass),"MAX DENSITY",location,"location"
print*,"__________________________________________"
inquire(file="Distance_from_BH",exist=iexist)
 if (.not. iexist) then

    open(121,file="Distance_from_BH",status='new',action='write',form='formatted')
    write(121,"(5(a))") "[x]"," ","[y]"," ","[z]"
    write(121,"(3(1x,es12.5))") xyzh(1,location)*udist,xyzh(2,location)*udist,xyzh(3,location)*udist 
    close(121)

  else

    open(121,file="Distance_from_BH",status='old',action='write',form='formatted',position="append")
    write(121,"(3(1x,es12.5))") xyzh(1,location)*udist,xyzh(2,location)*udist,xyzh(3,location)*udist  

  close(121)

  endif
 
 ieos = 2
!print*,maxval(rhoh(xyzh(4,:),pmass)),"MAX DENSITU FROM MAX RHO"
 call init_eos(ieos,ierr)
 print*,ieos,"ieos"
 !use sorting algorithm to sort the particles from the center of star as a function of radius.
 xpos(:) = star_centre(:)
 vpos(:) = vxyzu(1:3,location)
 velocity_norm = sqrt(dot_product(vpos(:),vpos(:)))

 print*,'com velocity',vpos(:)
 call set_r2func_origin(xpos(1),xpos(2),xpos(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
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
 !number_particle = 201
 !number_bins = number_tot/number_per_bin
 !ngrid = npart/number_particle
 number_particle = 0
 energy_verified_no = 0
 !if (mod(npart,number_particle) > 0) then
  !  ngrid = ngrid + 1
 !endif
print*,"**********************************"
 print"(a,i5)", 'number of bins = ',ngrid

 do j = 1, npart

    i  = iorder(j) !Access the rank of each particle in radius.

    !the position of the particle is calculated by subtracting the point of highest density.
    !xyzh is position wrt the black hole present at origin.
    pos(:) = xyzh(1:3,i) - xpos(:)

    !calculate the position which is the location of the particle.
    rad_test    = sqrt(dot_product(pos(:),pos(:)))
    potential_i = poten(i)

    !velocity
    vel(:)     = vxyzu(1:3,i) - vpos(:)

    velocity_norm = dot_product(vel(:),vel(:))
    kinetic_i     = 0.5*pmass*(velocity_norm)

    !total energy of single particle
    energy_i      = potential_i + kinetic_i

    !if energy is less than 0, we have bound system. We can accept these particles.
    if (energy_i < 0. .and. kinetic_i < 0.5*abs(potential_i)) then
      energy_verified_no = energy_verified_no + 1
      last_particle_with_neg_e = j
    endif

 enddo
 
 !based on the number of particles that satisfy the energy distribution formula, we will chose the
 !number of particles per bin
 !mod(npart,number_particle) > 0
print*,mod(energy_verified_no,500),"mod(energy_verified_no,500)"
 if (mod(energy_verified_no,500)>0) then
    number_particle = energy_verified_no/500 + 1
 print*,number_particle,"indie"
    print*,energy_verified_no/500,"energy verified no/500"
    if ((energy_verified_no/number_particle>=500.) .or. (energy_verified_no < 500) ) then 
        number_particle=number_particle
    else 
        number_particle = energy_verified_no/500
    endif
    print*,"We use 2 particles?)",number_particle
 elseif (mod(energy_verified_no,500)==0) then
    print*,"HERE WE ARE" 
    number_particle = energy_verified_no/500
 else  
    number_particle = 1
    print*,"HEREMwemare"
 endif
 
 print*,number_particle,"NUmber Particle", "Number of particles veritifed",energy_verified_no
 !we allocate arrays that can save 2000 data points.
 ngrid = 5000
 print*,"YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"
 print*,number_particle,"number_particle"
 allocate(rad_grid(ngrid),mass(ngrid),density(ngrid),energy_tot(npart))!rad_grid stores radius, stores radial velocity
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
 velocity_wrt_bh(:) = 0.
 X_in               = 0.71492308
 Y_in               = 0.27112283
 Z_in               = 1.-X_in-Y_in
 gmw                = 0.61 !mean molecular weight. Setting default value as gmw of sun.
 com_star(:)        = 0.
 bh_mass            = 1e6 !change the mass of the black hole here.
 count_new          = 1
 delta = reshape((/1,0,0,0,1,0,0,0,1/),shape(delta)) 
 i_matrix(:,:) = 0. 
L(:) = 0.
print*,"columns_compo",columns_compo,ngrid,"ngrid"
 skip_breakup = 0
 call assign_atomic_mass_and_number(comp_label,A_array,Z_array)
 print*,A_array,"A_array",Z_array,"Z_array"

 kinetic_add = 0.
 !allocating storage for composition of one particle.
 allocate(composition_i(columns_compo))
 allocate(composition_sum(columns_compo))
 allocate(composition_kepler(columns_compo,ngrid))
 composition_sum(:) = 0.
 composition_i(:) = 0.
 c_particle = 0 !no of particles that are bound the star.
 tot_thermal_energy = 0.
 tot_internal_energy = 0. 
 open(21,file="saving_pos",form='formatted')
 write(21,*) "[x]"," ","[y]", " ","[z]"
open(121,file="bin_angular_data") 
write(121,*) "[angular]"," ","[radius]"
open(10,file="Angular_data")
write(10,*) "[ANgular vel]"," ","[Radius]"
open(101,file="angular_analysis") 
write(101,"(15(a))") "[vx]"," ","[vy]"," ","[vz]", " ","[x]"," ","[y]"," ","[z]"," ","[pmass]"," ","[density]"

do j = 1, npart
   ! print*,rhoh(xyzh(4,j),pmass),"rhoh(xyzh(4,j),pmass",j,"j"
    !print*,"*******************************"
    i  = iorder(j) !Access the rank of each particle in radius.

    !the position of the particle is calculated by subtracting the point of highest density.
    !xyzh is position wrt the black hole present at origin.
    pos(:) = xyzh(1:3,i) - xpos(:)

    !calculate the position which is the location of the particle.
    rad_test    = sqrt(dot_product(pos(:),pos(:)))
    potential_i = poten(i)
    !velocity
    vel(:)     = vxyzu(1:3,i) - vpos(:)
    vel_sum(:) = vel_sum(:) + vel(:)
    velocity_norm = dot_product(vel(:),vel(:))
    kinetic_i     = 0.5*pmass*(velocity_norm)

    !total energy of single particle
    energy_i      = potential_i + kinetic_i
    
   call cross_product3D(pos(:),vel(:),Li)
    omega_val(:) =(Li(:)*unit_velocity*udist)/((rad_test*udist)**2)
    val_omega = dot_product(omega_val(:),omega_val(:))
    

    !if (j .ne. 1) then 
         if (val_omega < (gg*pmass*j*umass)/((rad_test*udist)**3)) then 
           count_new = count_new+1
          endif 
    !endif 
    
    if (j==1) then 
       val_omega = 0.
    endif
    !if energy is less than 0, we have bound system. We can accept these particles.
   ! print*, energy_i,"energy_i",kinetic_i,"kinetic_i",potential_i,"potential_i"
    if (energy_i < 0. .and. kinetic_i < 0.5*abs(potential_i) ) then
    tot_thermal_energy = tot_thermal_energy + kinetic_i
    write(21,*) xyzh(1:3,i)

       rad  = rad_test
       c_particle = 1+c_particle
       !print*,"c_particle",j,"j"
       !add 1 to no_in_bin if energy is < 0.
       no_in_bin  = no_in_bin + 1
       !Calculating COM of star and velocity of the COM. This is only done for particles with
       !negative energy wrt black hole.
       com_star(:) = com_star(:) + pmass*xyzh(1:3,i)
       !momentum of star wrt black hole. need this to find com velocity.
       velocity_wrt_bh(:)   = velocity_wrt_bh(:) + pmass*vxyzu(1:3,i)

       !radial momentum
       !we skip the first particle as its the one that exists at the center of star and hence will give infinite rad_vel as rad = 0.
       if (rad > 0.) then
          rad_velocity = dot_product(vel(:),pos(:))/rad
          momentum     = rad_velocity*pmass
          rad_vel_sum  = rad_vel_sum + momentum
       endif

       !angular velocity of the shells of the star.
       call cross_product3D(pos(:),vel(:),Li)
       moment_of_inertia = moment_of_inertia + rad**2*pmass  !moment of intertia
       L(:)          = Li(:)*pmass
       omega_sum(:)      = omega_sum(:) + L(:)
      
       !print*,"++++++++++++++++++++++++++++++++++"
      
        val_omega = sqrt(dot_product(L(:),L(:)))

       !moment of inertia matrix
       matrix1=reshape(pos,shape(matrix1))
       matrix2=reshape(pos,shape(matrix2))
       result_matrix = matmul(matrix1,matrix2)
       !print*,"i_matrix",i_matrix,"j",j
       i_matrix = i_matrix + pmass*(rad**2*delta - result_matrix) 
       
!print*,"**********************************"
 !print*,i_matrix*umass*udist**2,"i_matrix"
!print*,inverse(i_matrix*umass*udist**2, 3),"inverse(i_matrix, 3)"
 ! print*,"**********",j,"j" 

      !density
      
       density_i   = rhoh(xyzh(4,i),pmass)
       print*,density_i,"density_i",j,"j"
       density_sum = density_sum + density_i
       if (j==1) then 
          val_omega = 0. 
       endif
       write(10,*) (val_omega/(rad**2*pmass))*(1/utime),rad*udist
!print*,"PPPPPPPPPPPPPPPPP"
      ! print*,val_omega/(rad**2*pmass),"val_omega/(rad**2*pmass)",j,"j"
       !internal energy
       u_i   = vxyzu(4,i)
       u_sum = u_sum + u_i
       
       tot_internal_energy = tot_internal_energy + u_i
       
       !composition
       if (columns_compo /= 0) then
          composition_i(:) = interpolate_comp(:,i)
       endif
       composition_sum(:) = composition_sum(:) + composition_i(:)

       !calculate mean molecular weight that is required by the eos module using the
       !mass fractions for each particle.
       !do not consider neutron which is the first element of the composition_i array.

       call calculate_mu(A_array,Z_array,composition_i,columns_compo,mu)

       gmw = 1./mu
        write(101,"(8(1x,es12.5))") vel(:)*unit_velocity,pos(:)*udist,pmass*umass,density_i*unit_density

       eni_input = u_i
       !call eos routine
       !call equationofstate(ieos,ponrhoi,spsoundi,density_i,xyzh(1,i),xyzh(2,i),xyzh(3,i),eni=eni_input, tempi=temperature_i,mu_local=1./mu_i)
       call equationofstate(ieos,ponrhoi,spsoundi,density_i,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi=temperature_i, eni=eni_input)

       pressure_i      = ponrhoi*density_i
       pressure_sum    = pressure_sum + pressure_i
       temperature_sum = temperature_sum + temperature_i
       !print*,j,'j'
    endif
    !now, if the bin has the same number of particles as we specified initially, we save the values
    if ((no_in_bin == number_particle) .and. (j .ne. npart) .and. (j .ne. last_particle_with_neg_e)) then

       !calculate the position which is the location of the particle.
       do index_i = j,npart-1
          !we check that the first particle in the next bin is bound and if it is we calculate the radius as average
          s              = index_i + 1 !we consider the first particle in the i+1 bin and then use it.
          m              = iorder(s)
          pos(:)         = xyzh(1:3,m) - xpos(:)
          potential_i    = poten(m)
          vel(:)         = vxyzu(1:3,m) - vpos(:)
          velocity_norm  = dot_product(vel(:),vel(:))
          kinetic_i      = 0.5*pmass*(velocity_norm)
          energy_i       = potential_i + kinetic_i
          if (energy_i < 0. .and. kinetic_i < 0.5*abs(potential_i)) then
             rad_next     = sqrt(dot_product(pos(:),pos(:)))
             exit
          endif
       enddo
       print*,"_______________________________________"
       print*,j,"j",rad*udist,"rad",rad_next*udist,"rad_next",s,"s"
       print*,"______________________________________"
       rad_grid(ibin)   = (rad + rad_next)/2
   elseif (j == last_particle_with_neg_e)  then
       rad_grid(ibin) = rad
       print*,"__________________________"
       print*,rad*udist,"last rad"
       print*,"__________________________"
       !if number of bin does not have the same number of particles we loop again. so cycle for j<npart
    elseif (j<npart) then
       cycle
endif 
       !elseif j is last then save the last bin if the number of particle /= no_in_bin
    !elseif (j == last_particle_with_neg_e)  then
       !rad_grid(ibin) = rad
!print*,"HERE"
 !      print*,"------------------------"
  !     print*,rad,"rad",j,"j",ibin,"ibin"
   ! endif

    if (no_in_bin == number_particle .or. j == last_particle_with_neg_e) then
       bin_mass(ibin)             = no_in_bin*pmass !bin mass
       mass(ibin)                 = (c_particle)*pmass !mass of paricles < r.
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
       inverse_of_i      = inverse(i_matrix, 3)
       !print*,ibin,"ibin","///////////////////"
       !print*,"-------------------------------------------"
     !  print*,j,"j"
     !print*,"i_matrix"
      ! print*,i_matrix*umass*udist**2
      
       omega_reshape = reshape(omega_sum(:),(/3,1/))
       !print*,"----------------------------------------"
       !print*,"omega reshape"
       !print*,omega_reshape
       !print*,"_---------------------------------------"
       !print*,"inverse of i"
       !print*,inverse_of_i
        
       matrix_result = matmul(inverse_of_i,omega_reshape)  !omega_sum is L_sum
       final_val_omega = reshape(matrix_result,(/3/))
       !print*,"MATRIX INVERSE"
    
     ! print*,inverse(i_matrix*umass*udist**2,3)
    !  print*,matmul(inverse_of_i,i_matrix)
       if (no_in_bin == 1) then 
          if (rad==0.) then 
             angular_vel_3D(:,ibin) = omega_sum(:)
          else
              angular_vel_3D(:,ibin) = omega_sum(:) / (rad**2*pmass*no_in_bin)
          endif
       else
           angular_vel_3D(:,ibin) = final_val_omega
       endif
       val_omega = sqrt(dot_product(angular_vel_3D(:,ibin),angular_vel_3D(:,ibin)))
       write(121,*) val_omega/utime,rad_grid(ibin)*udist
        
     !  print*,"ANgular vel 3D"
     !  print*,angular_vel_3D(:,ibin)/utime
       !calculating entropy
       !implementing entropy from the Sackur-Tetrode equation.
       entropy_array(ibin) = entropy(density(ibin),pressure(ibin),mu,ientropy=2,ierr=ierr)
       entropy_array(ibin) = entropy_array(ibin)/(kboltz*avogadro)
       if (ierr/=0) then
          print*, 'Entropy is calculated incorrectly'
       endif
       !print*,"YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"
       !print*,ibin,"IBIN",j,"j"
       !print*,"YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"
       tot_energy         = 0
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
       kinetic_add        = 0.
       i_matrix(:,:)      = 0.
    endif
 enddo
 !close(1)
close(21)
close(121)
close(10)
close(101)
 correct_ngrid = ibin-1
 print*,c_particle,'c_particle',npart-c_particle
 print*,'----------------------------------------------------------------------'
 print*,'orbit parameters'
print*,";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
print*,"PMASS",pmass
 !this is centre of mass of the star
 print*,star_centre*udist,'center of star'
 com_star(:) = (com_star(:)/(c_particle*pmass))*udist
 position_bh = sqrt(dot_product(com_star(:),com_star(:)))
 !velocity of centre of mass of star
 velocity_wrt_bh(:) = (velocity_wrt_bh(:) /(c_particle*pmass))*unit_velocity
 velocity_bh = sqrt(dot_product(velocity_wrt_bh(:),velocity_wrt_bh(:)))
 mass_star = c_particle*pmass

 !we calculate angular momentum vector on the orbit around black hole.
 call cross_product3D(com_star,velocity_wrt_bh,angular_momentum_h)
print*,mass(ibin-1),"MASS OFBIN"
print*,mass_star,"mass_star",mass_star*umass,ngrid,"ngrid"
print*,"------"
 print*,com_star(:),'position of star wrt bh'
 print*,velocity_wrt_bh,'velocity of star wrt bh'
 print*,angular_momentum_h(:),'angular momentum h'

 !energy total is total energy of star. If negative, star is bound to BH but if positive then its unbound.
 !Calculate orbital paramters
 energy_total = 0.5*velocity_bh**2 - (gg*1e6*umass)/position_bh
 print*,"here"
 if (energy_total < 0.) then
    print*, escape(velocity_bh,bh_mass,position_bh)
    print*,"negative energy"
    call orbital_parameters(angular_momentum_h,bh_mass,mass_star,com_star,position_bh,velocity_wrt_bh)
print*,"here2"
 else
    print*, escape(velocity_bh,bh_mass,position_bh)
    print*,'positive energy of the star on orbit'
 endif
 print*,npart-c_particle,'unbound number of particles',correct_ngrid
 print*,npart-count_new,"count_new"
 print*,'----------------------------------------------------------------------'
print*,"here3",mass_star,"mass_star"
print*,"*******************"
print*,tot_thermal_energy,"tot_thermal"
 call write_mass_of_star_and_no(numfile,mass_star,star_centre(:),maxval(den_all(:)),tot_thermal_energy, tot_internal_energy,rad_grid(correct_ngrid),energy_verified_no)
print*,"here4"
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
    print*,"comp_label",comp_label

    open(13, file=filename)
    call skip_header(13,nheader,ierr)
    do k = 1, n_rows
       read(13,*,iostat=ierr) interpolate_comp(:,k)
    enddo
    close(13)
    print*, '>>>>>> done'
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

 print*,size_to_allocate,"size_to_allocate"
 allocate(A_array(size_to_allocate), Z_array(size_to_allocate), new_comp_label(size_to_allocate))
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
!  This routine for writing the mass and time of each particle.
!  This will help me plot them to see if mass becomes constant with
!  evolution time.
!+
!----------------------------------------------------------------
subroutine write_mass_of_star_and_no(numfile,mass_star,central_position,central_den,tot_thermal_energy, tot_internal_energy,radius,energy_verified_no)
 use units, only:umass,udist,unit_energ
 integer, intent(in)  :: numfile,energy_verified_no
 real,    intent(in)  ::mass_star,tot_thermal_energy,tot_internal_energy,radius
 character(len=120)   :: filename
logical                 :: iexist
real :: central_position(:),central_den
 filename = 'all_analysis.dat'
print*,udist,"UDIST"
print*,mass_star,"mass_star"
print*,"________"
print*,energy_verified_no,"energy verified no"
!  if (numfile==00000) then
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then 

    open(21,file=filename,status='new',action='write',form='formatted')
    write(21,"(13(a))") "[Mass of remnant]"," ", "[File numbe]"," ","[Central Density]"," ","[Tot thermal Energy]"," ","[Tot Internal Energy]"," ","[Radius in cm]"," ","[No of particles in star]"
    write(21,"(es12.5,1x,i4,4(1x,es12.5),1x,i4)") mass_star*umass, numfile,central_den,tot_thermal_energy*unit_energ, tot_internal_energy*unit_energ,radius*udist,energy_verified_no
    close(21)
    
  else

    open(21,file=filename,status='old',action='write',form='formatted',position="append")
     write(21,"(es12.5,1x,i4,4(1x,es12.5),1x,i4)") mass_star*umass, numfile,central_den,tot_thermal_energy*unit_energ, tot_internal_energy*unit_energ,radius*udist,energy_verified_no
    close(21)

  endif
print*,mass_star*umass, "star"
end subroutine write_mass_of_star_and_no

end module analysis
