!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for common envelope simulations
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, dust_formation, energies, eos,
!   eos_gasradrec, eos_mesa, extern_corotate, io, ionization_mod, kernel,
!   mesa_microphysics, part, physcon, prompting, ptmass, setbinary,
!   sortutils, table_utils, units, vectorutils
!

 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass,poten,ihsoft,ihacc,&
                        rhoh,nsinkproperties,maxvxyzu,maxptmass,isdead_or_accreted
 use units,        only:print_units,umass,utime,udist,unit_ergg,unit_density,&
                        unit_pressure,unit_velocity,unit_Bfield,unit_energ
 use physcon,      only:gg,pi,c,Rg
 use io,           only:fatal
 use prompting,    only:prompt
 use centreofmass, only:get_centreofmass, reset_centreofmass
 use energies,     only:compute_energies,ekin,etherm,epot,etot
 use ptmass,       only:get_accel_sink_gas,get_accel_sink_sink
 use kernel,       only:kernel_softening,radkern,wkern,cnormk
 use eos,          only:equationofstate,ieos,init_eos,X_in,Z_in,gmw,get_spsound,done_init_eos
 use eos_gasradrec,only:irecomb
 use eos_mesa,     only:get_eos_kappa_mesa,get_eos_pressure_temp_mesa,&
                        get_eos_various_mesa,get_eos_pressure_temp_gamma1_mesa
 use setbinary,    only:Rochelobe_estimate,L1_point
 use sortutils,    only:set_r2func_origin,r2func_origin,indexxfunc
 use table_utils,  only:logspace
 implicit none
 character(len=20), parameter, public :: analysistype = 'common_envelope'
 integer                              :: analysis_to_perform
 integer                              :: dump_number = 0
 real                                 :: omega_corotate=0,init_radius,rho_surface,gamma
 logical, dimension(5)                :: switch = .false.
 public                               :: do_analysis
 public                               :: tconv_profile,get_interior_mass ! public = no unused fn warning
 public                               :: planet_destruction,total_dust_mass ! make public to avoid compiler warning
 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in)    :: particlemass,time
 integer                         :: unitnum,i,ncols
 logical                         :: requires_eos_opts

 !case 5 variables
 real                         :: rhopart

 !case 7 variables
 character(len=17), allocatable :: columns(:)

 !case 12 variables
 real                         :: etoti, ekini, einti, epoti, phii

 real, dimension(3)           :: com_xyz, com_vxyz
 real, dimension(3)           :: xyz_a, vxyz_a
 real, allocatable            :: histogram_data(:,:)
 real                         :: ang_vel

 real :: pres_1i, proint_1i, peint_1i, temp_1i
 real :: troint_1i, teint_1i, entrop_1i, abad_1i, gamma1_1i, gam_1i

 !case 16 variables
 real, allocatable :: thermodynamic_quantities(:,:)
 real, allocatable :: radius_1i, dens_1i


 !chose analysis type
 if (dump_number==0) then
    print "(41(a,/))", &
            ' 1) Sink separation', &
            ' 2) Bound and unbound quantities', &
            ' 3) Energies', &
            ' 4) Profile from centre of mass', &
            ' 5) Roche-lobe utils', &
            ' 6) Star stabilisation suite', &
            ' 7) Simulation units and particle properties', &
            ' 8) Output .divv', &
            ' 9) EoS testing', &
            '11) Profile of newly unbound particles', &
            '12) Sink properties', &
            '13) MESA EoS compute total entropy and other average td quantities', &
            '14) MESA EoS save on file thermodynamical quantities for all particles', &
            '15) Gravitational drag on sinks', &
            '16) CoM of gas around primary core', &
            '17) Miscellaneous', &
            '18) J-E plane', &
            '19) Rotation profile', &
            '20) Energy profile', &
            '21) Recombination statistics', &
            '22) Optical depth profile', &
            '23) Particle tracker', &
            '24) Unbound ion fraction', &
            '25) Optical depth at recombination', &
            '26) Envelope binding energy', &
            '27) Print dumps number matching separation', &
            '28) Companion mass coordinate vs. time', &
            '29) Energy histogram',&
            '30) Analyse disk',&
            '31) Recombination energy vs time',&
            '32) Binding energy profile',&
            '33) planet_rvm',&
            '34) Velocity histogram',&
            '35) Unbound temperature',&
            '36) Planet mass distribution',&
            '37) Planet profile',&
            '38) Velocity profile',&
            '39) Angular momentum profile',&
            '40) Keplerian velocity profile',&
            '41) Total dust mass'
    analysis_to_perform = 1
    call prompt('Choose analysis type ',analysis_to_perform,1,41)
 endif

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 call adjust_corotating_velocities(npart,particlemass,xyzh,vxyzu,&
                                   xyzmh_ptmass,vxyz_ptmass,omega_corotate,dump_number)

 ! List of analysis options that require specifying EOS options
 requires_eos_opts = any((/2,3,4,6,8,9,11,13,14,15,20,21,22,23,24,25,26,29,30,31,32,33,35,41/) == analysis_to_perform)
 if (dump_number == 0 .and. requires_eos_opts) call set_eos_options(analysis_to_perform)

 select case(analysis_to_perform)
 case(1) !sink separation
    call separation_vs_time(time)
 case(2) !bound and unbound quantities
    call bound_mass(time,npart,particlemass,xyzh,vxyzu)
 case(3) !Energies and bound mass
    call calculate_energies(time,npart,particlemass,xyzh,vxyzu)
 case(4) !Profile from COM (can be used for stellar profile)
    call create_profile(time, num, npart, particlemass, xyzh, vxyzu)
 case(5) !Mass within roche lobes
    call roche_lobe_values(time,npart,particlemass,xyzh,vxyzu)
 case(6) !Star stabilisation suite
    call star_stabilisation_suite(time,npart,particlemass,xyzh,vxyzu)
 case(7) !Units
    call print_simulation_parameters(npart,particlemass)
 case(8) !Output .divv
    call output_divv_files(time,dumpfile,npart,particlemass,xyzh,vxyzu)
 case(9) !EoS testing
    call eos_surfaces
 case(11) !New unbound particle profiles in time
    call unbound_profiles(time,num,npart,particlemass,xyzh,vxyzu)
 case(19) ! Rotation profile
    call rotation_profile(time,num,npart,xyzh,vxyzu)
 case(20) ! Energy profile
    call energy_profile(time,npart,particlemass,xyzh,vxyzu)
 case(21) ! Recombination statistics
    call recombination_stats(time,num,npart,particlemass,xyzh,vxyzu)
 case(22) ! Optical depth profile
    call tau_profile(time,num,npart,particlemass,xyzh)
 case(23) ! Particle tracker
    call track_particle(time,particlemass,xyzh,vxyzu)
 case(24) ! Unbound ion fractions
    call unbound_ionfrac(time,npart,particlemass,xyzh,vxyzu)
 case(25) ! Optical depth at recombination
    call recombination_tau(time,npart,particlemass,xyzh,vxyzu)
 case(26) ! Calculate binding energy outside core
    call env_binding_ene(npart,particlemass,xyzh,vxyzu)
 case(27) ! Print dump number corresponding to given set of sink-sink separations
    call print_dump_numbers(dumpfile)
 case(28) ! Companion mass coordinate (spherical mass shells) vs. time
    call m_vs_t(time,npart,particlemass,xyzh)
 case(29) ! Energy histogram
    call energy_hist(time,npart,particlemass,xyzh,vxyzu)
 case(30) ! Analyse disk around companion
    call analyse_disk(num,npart,particlemass,xyzh,vxyzu)
 case(31) ! Recombination energy vs. time
    call erec_vs_t(time,npart,particlemass,xyzh)
 case(32) ! Binding energy profile
    call create_bindingEnergy_profile(time,num,npart,particlemass,xyzh,vxyzu)
 case(33) ! Planet coordinates and mass
    call planet_rvm(time,particlemass,xyzh,vxyzu)
 case(34) ! Velocity histogram
    call velocity_histogram(time,num,npart,particlemass,xyzh,vxyzu)
 case(35) ! Unbound temperatures
    call unbound_temp(time,npart,particlemass,xyzh,vxyzu)
 case(36) ! Planet mass distribution
    call planet_mass_distribution(time,num,npart,xyzh)
 case(37) ! Calculate planet profile
    call planet_profile(num,dumpfile,particlemass,xyzh,vxyzu)
 case(38) ! Velocity profile
    call velocity_profile(time,num,npart,particlemass,xyzh,vxyzu)
 case(39) ! Angular momentum profile
    call angular_momentum_profile(time,num,npart,particlemass,xyzh,vxyzu)
 case(40) ! Keplerian velocity profile
    call vkep_profile(time,num,npart,particlemass,xyzh,vxyzu)
 case(41) !Total dust mass
    call total_dust_mass(time,npart,particlemass,xyzh)
 case(12) !sink properties
    call sink_properties(time,npart,particlemass,xyzh,vxyzu)
 case(13) !MESA EoS compute total entropy and other average thermodynamical quantities
    call bound_unbound_thermo(time,npart,particlemass,xyzh,vxyzu)
 case(14) !MESA EoS save on file thermodynamical quantities for all particles
    allocate(thermodynamic_quantities(5,npart))
    do i=1,npart

       !particle radius
       radius_1i = distance(xyzh(1:3,i)) * udist

       !particles density in code units
       rhopart = rhoh(xyzh(4,i), particlemass)
       dens_1i = rhopart * unit_density

       !gets entropy for the current particle
       call get_eos_various_mesa(rhopart*unit_density,vxyzu(4,i) * unit_ergg, &
                                 pres_1i,proint_1i,peint_1i,temp_1i,troint_1i, &
                                 teint_1i,entrop_1i,abad_1i,gamma1_1i,gam_1i)

       !stores everything in an array
       thermodynamic_quantities(1,i) = radius_1i
       thermodynamic_quantities(2,i) = dens_1i
       thermodynamic_quantities(3,i) = pres_1i
       thermodynamic_quantities(4,i) = temp_1i
       thermodynamic_quantities(5,i) = entrop_1i

    enddo
    ncols = 5
    allocate(columns(ncols))
    columns = (/'      radius', &
                '     density', &
                '    pressure', &
                ' temperature', &
                '     entropy'/)
    call write_file('td_quantities', 'thermodynamics', columns, thermodynamic_quantities, npart, ncols, num)

    unitnum = unitnum + 1
    deallocate(thermodynamic_quantities)

 case(15) !Gravitational drag on sinks
    call gravitational_drag(time,npart,particlemass,xyzh,vxyzu)

 case(16)
    call get_core_gas_com(time,npart,xyzh,vxyzu)

 case(17)
    ncols = 6
    allocate(columns(ncols))
    columns = (/'           x', &
                '           y', &
                '           z', &
                '           r', &
                'spec. energy', &
                ' omega ratio'/)

    call orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)

    ang_vel = 0.

    do i=1,nptmass
       if (xyzmh_ptmass(4,i) > 0.) then
          xyz_a(1:3) = xyzmh_ptmass(1:3,i) - com_xyz(1:3)
          vxyz_a(1:3) = vxyz_ptmass(1:3,i) - com_vxyz(1:3)
          ang_vel = ang_vel + (-xyz_a(2) * vxyz_a(1) + xyz_a(1) * vxyz_a(2)) / dot_product(xyz_a(1:2), xyz_a(1:2))
       endif
    enddo

    ang_vel = ang_vel / 2.

    allocate(histogram_data(6,npart))

    do i=1,npart
       xyz_a(1:3) = xyzh(1:3,i) - com_xyz(1:3)
       vxyz_a(1:3) = vxyzu(1:3,i) - com_vxyz(1:3)

       call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)
       histogram_data(1:3,i) = xyzh(1:3,i)
       histogram_data(4,i) = distance(xyz_a(1:3))
       histogram_data(5,i) = epoti + ekini
       histogram_data(6,i) = (-xyz_a(2) * vxyz_a(1) + xyz_a(1) * vxyz_a(2)) / dot_product(xyz_a(1:2), xyz_a(1:2))
       histogram_data(6,i) = (histogram_data(6,i) - ang_vel) / ang_vel
    enddo

    call write_file('specific_energy_particles', 'histogram', columns, histogram_data, size(histogram_data(1,:)), ncols, num)

    deallocate(histogram_data)

 case(18)
    call J_E_plane(num,npart,particlemass,xyzh,vxyzu)
 end select
 !increase dump number counter
 dump_number = dump_number + 1

end subroutine do_analysis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                Analysis  routines                !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine total_dust_mass(time,npart,particlemass,xyzh)
 use part,           only:nucleation,idK3,idK0,idK1, idJstar
 use dust_formation, only:set_abundances, mass_per_H
 use physcon, only:atomic_mass_unit
 real, intent(in)               :: time,particlemass,xyzh(:,:)
 integer, intent(in)            :: npart
 integer                        :: i,ncols,j
 real, dimension(2)             :: dust_mass
 character(len=17), allocatable :: columns(:)
 real, allocatable              :: temp(:) !npart
 real                           :: median,mass_factor,grain_size
 real, parameter :: a0 = 1.28e-4 !radius of a carbon atom in micron

 call set_abundances !initialize mass_per_H
 dust_mass = 0.
 ncols = 2
 print *,'size(nucleation,1) = ',size(nucleation,1)
 print *,'size(nucleation,2) = ',size(nucleation,2)
 allocate(columns(ncols),temp(npart))
 columns = (/'Dust mass [Msun]', &
             'median size [um]'/)
 j=0
 mass_factor = 12.*atomic_mass_unit*particlemass/mass_per_H
 do i = 1,npart
    if (.not. isdead_or_accreted(xyzh(4,i))) then
       dust_mass(1) = dust_mass(1) + nucleation(idK3,i) *mass_factor
       grain_size = a0*nucleation(idK1,i)/(nucleation(idK0,i)+1.0E-99) !in micron
       if (grain_size > a0) then
          j = j+1
          temp(j) = grain_size
       endif
    endif
 enddo

 call sort(temp,j)
 if (mod(j,2)==0) then !npart
    median = (temp(j/2)+temp(j/2+1))/2.0 !(temp(npart/2)+temp(npart/2+1))/2.0
 else
    median = (temp(j/2)+temp(j/2+1))/2.0 !temp(npart/2+1)
 endif

 dust_mass(2) = median

 call write_time_file('total_dust_mass_vs_time', columns, time, dust_mass, ncols, dump_number)
 !after execution of the analysis routine, a file named "total_dust_mass_vs_time.ev" appears
 deallocate(columns,temp)

end subroutine total_dust_mass

! --------------------------------------------------------------------
! integer function  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

integer function  FindMinimum(x, Start, Fin)
 implicit  none
 integer, intent(in)                   :: start, fin
 real, dimension(Fin), intent(in) :: x
 real                            :: minimum
 integer                            :: location
 integer                            :: i

 minimum  = x(start)          ! assume the first is the min
 location = start             ! record its position
 do i = start+1, fin          ! start with next elements
    if (x(i) < minimum) then  !   if x(i) less than the min?
       minimum  = x(i)        !      yes, a new minimum found
       location = i                !      record its position
    endif
 enddo
 findminimum = location            ! return the position
end function FindMinimum

! --------------------------------------------------------------------
! subroutine  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

subroutine  Sort(x, longitud)
 implicit  none
 integer, intent(in)                   :: longitud
 real, dimension(longitud), intent(inout) :: x
 integer                               :: i
 integer                               :: location

 do i = 1, longitud-1             ! except for the last
    location = findminimum(x, i, longitud)  ! find min from this to last
    call swap(x(i), x(location))  ! swap this and the minimum
 enddo
end subroutine Sort


!----------------------------------------------------------------
!+
!  Separation vs. time
!+
!----------------------------------------------------------------
subroutine separation_vs_time(time)
 real, intent(in)               :: time
 character(len=17), allocatable :: columns(:)
 real                           :: sink_separation(4,nptmass-1)
 integer                        :: i,ncols
 ncols = 4*(nptmass-1)
 allocate(columns(ncols))

 do i=1,(nptmass-1)
    call separation_vector(xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,i+1),sink_separation(1:4,i))

    write(columns((i*4)-3), '(A11,I1)') '    x sep. ', i
    write(columns((i*4)-2), '(A11,I1)') '    y sep. ', i
    write(columns((i*4)-1), '(A11,I1)') '    z sep. ', i
    write(columns((i*4)),   '(A11,I1)') '      sep. ', i
 enddo

 call write_time_file('separation_vs_time', columns, time, sink_separation, ncols, dump_number)
 deallocate(columns)
end subroutine separation_vs_time


!----------------------------------------------------------------
!+
!  Output planet position (x,y,z,r) and velocity (vx,vy,vz,|v|)
!  relative to core, instantaneous mass according to different
!  criteria (m1,m2,m3,m4,m5), max. density, and min. entropy
!
!  For small dumps, only (x,y,z,r) and rhomax may be determined.
!  All other quantities will be outputted as zero.
!+
!----------------------------------------------------------------
subroutine planet_rvm(time,particlemass,xyzh,vxyzu)
 use eos, only:entropy
 real, intent(in)               :: time,xyzh(:,:),vxyzu(:,:),particlemass
 character(len=17), allocatable :: columns(:)
 real, dimension(3)             :: planet_com,planet_vel,sep,vel
 real                           :: rhoi,rhoprev,sepi,si,smin,presi,Rthreshold
 real, allocatable              :: data_cols(:),mass(:),vthreshold(:)
 integer                        :: i,j,ncols,maxrho_ID,ientropy,Nmasks
 integer, save                  :: nplanet
 integer, allocatable, save     :: planetIDs(:)
 logical                        :: isfulldump

 if (.not. done_init_eos) call fatal("planet_rvm","EOS has not been initialised.")

 ncols = 15
 allocate(data_cols(ncols),columns(ncols))
 columns = (/'       x sep', &
             '       y sep', &
             '       z sep', &
             '         sep', &
             '          vx', &
             '          vy', &
             '          vz', &
             '           v', &
             '          m1', &
             '          m2', &
             '          m3', &
             '          m4', &
             '          m5', &
             '      rhomax', &
             '        smin'/)

 if (dump_number == 0) call get_planetIDs(nplanet,planetIDs)
 isfulldump = (vxyzu(4,1) > 0.)

 ! Find highest density and lowest entropy in planet
 rhoprev = 0.
 maxrho_ID = 1
 smin = huge(0.)
 ientropy = 1
 ieos = 2
 gamma = 5./3.
 do i = 1,nplanet
    rhoi = rhoh(xyzh(4,planetIDs(i)), particlemass)
    if (rhoi > rhoprev) then
       maxrho_ID = planetIDs(i)
       rhoprev = rhoi
    endif

    if (isfulldump) then
       presi = (gamma-1.)*vxyzu(4,i)
       si = entropy(rhoi*unit_density,presi*unit_pressure,gmw,ientropy)
       smin = min(smin,si)
    endif
 enddo

 planet_com = xyzh(1:3,maxrho_ID)
 sep = planet_com - xyzmh_ptmass(1:3,1)

 if (isfulldump) then
    planet_vel = vxyzu(1:3,maxrho_ID)
    vel = planet_vel - vxyz_ptmass(1:3,1)
 else
    vel = 0.
    smin = 0.
 endif

 ! Sum planet mass according to criterion
 Nmasks = 5  ! Number of velocity thresholds for calculating planet mass
 allocate(mass(Nmasks),vthreshold(Nmasks))
 mass = 0.
 if (isfulldump) then
    Rthreshold = 0.21  ! Radius criterion to be considered part of planet
    vthreshold = (/0.1,0.3,0.5,0.7,0.9/) ! Allowed fractional deviation in particle velocity from velocity of densest planet particle
    do i = 1,nplanet
       sepi = separation(xyzh(1:3,planetIDs(i)), planet_com)
       do j = 1,Nmasks
          if ( (sepi < Rthreshold) .and. (abs(1. - dot_product(vxyzu(1:3,planetIDs(i)),planet_vel)/&
                dot_product(planet_vel,planet_vel)) < vthreshold(j)) ) then ! vi dot vp / vp^2 > threshold
             mass(j:Nmasks) = mass(j:Nmasks) + 1.
             exit
          endif
       enddo
    enddo
    mass = mass * particlemass
 endif

 data_cols = (/ sep(1), sep(2), sep(3), distance(planet_com),&
                vel(1), vel(2), vel(3), distance(vel),&
                mass(1), mass(2), mass(3), mass(4), mass(5), rhoprev, smin /)
 call write_time_file('planet_rvm', columns, time, data_cols, ncols, dump_number)

 deallocate(data_cols,columns,mass,vthreshold)

end subroutine planet_rvm


!----------------------------------------------------------------
!+
!  Output radial distribution of planetary material
!+
!----------------------------------------------------------------
subroutine planet_mass_distribution(time,num,npart,xyzh)
 integer, intent(in)          :: npart,num
 real, intent(in)             :: time
 real, intent(inout)          :: xyzh(:,:)
 real, allocatable            :: rad_part(:),dist_part(:),hist_var(:)
 real                         :: mina,maxa,xyz_origin(3)
 character(len=17)            :: filename
 character(len=100)           :: data_formatter,headerline
 integer                      :: i,iu,nbins
 integer, save                :: nplanet
 integer, allocatable, save   :: planetIDs(:)

 if (dump_number == 0) call get_planetIDs(nplanet,planetIDs)

 nbins = 1000 ! Radial bins
 mina = 0.
 maxa = 4.2

 allocate(rad_part(nplanet),dist_part(nplanet),hist_var(nbins))
 filename = ' planet_m_dist.ev'
 xyz_origin = xyzmh_ptmass(1:3,1)

 dist_part = 0.
 rad_part = 0.
 do i = 1,nplanet
    rad_part(i) = separation(xyzh(1:3,planetIDs(i)),xyz_origin)
    dist_part(i) = 1.
 enddo

 call histogram_setup(rad_part,dist_part,hist_var,nplanet,maxa,mina,nbins,.false.,.false.)

 write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
 if (num == 0) then
    open(newunit=iu, file=trim(adjustl(filename)), status='replace')
    write(headerline, "(a,i5,a,f5.2,a,f5.2)") "# Planet mass distribution, nbins = ", nbins,", min a = ", mina, ", max a = ", maxa
    write(iu, "(a)") headerline
    close(unit=iu)
 endif
 open(newunit=iu, file=trim(adjustl(filename)), position='append')
 write(iu,data_formatter) time,hist_var(:)
 close(unit=iu)

 deallocate(rad_part,dist_part,hist_var)

end subroutine planet_mass_distribution


!----------------------------------------------------------------
!+
!  Companion mass coordinate (spherical mass shells) vs. time
!+
!----------------------------------------------------------------
subroutine m_vs_t(time,npart,particlemass,xyzh)
 integer, intent(in) :: npart
 real, intent(in)    :: time,particlemass,xyzh(:,:)
 character(len=17)   :: colname
 real                :: sinksinksep,mass(1)
 integer             :: i,k
 integer, allocatable :: iorder(:)

 allocate(iorder(npart))

 call set_r2func_origin(xyzmh_ptmass(1,1),xyzmh_ptmass(2,1),xyzmh_ptmass(3,1)) ! Order particles by distance from core
 call indexxfunc(npart,r2func_origin,xyzh,iorder)

 sinksinksep = separation(xyzmh_ptmass(1:3,1), xyzmh_ptmass(1:3,2))
 do i=1,npart
    k = iorder(i)
    if (separation(xyzh(1:3,k), xyzmh_ptmass(1:3,1)) > sinksinksep) exit
 enddo

 mass = i*particlemass + xyzmh_ptmass(4,1)
 write(colname, '(A11)') ' mass coord'
 call write_time_file('            m_vs_t',colname,time,mass,1,dump_number)

 deallocate(iorder)

end subroutine m_vs_t


!----------------------------------------------------------------
!+
!  Bound mass
!+
!----------------------------------------------------------------
subroutine bound_mass(time,npart,particlemass,xyzh,vxyzu)
 use part,           only:eos_vars,itemp
 use ptmass,         only:get_accel_sink_gas
 use ionization_mod, only:calc_thermal_energy
 use vectorutils,    only:cross_product3D
 integer, intent(in)            :: npart
 real, intent(in)               :: time,particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 real                           :: etoti,ekini,epoti,phii,einti,ethi
 real                           :: E_H2,E_HI,E_HeI,E_HeII
 real, save                     :: Xfrac,Yfrac,Zfrac
 real                           :: rhopart,ponrhoi,spsoundi,tempi,dum1,dum2,dum3
 real, dimension(3)             :: rcrossmv
 real, dimension(28)            :: bound
 integer                        :: i,bound_i,ncols
 integer, parameter             :: ib=1,ibt=9,ibe=17
 character(len=17), allocatable :: columns(:)

 if (.not. done_init_eos) call fatal("bound_mass","EOS has not been initialised.")

 ncols = 28
 bound = 0.
 allocate(columns(ncols))
 columns = (/'  b num part', & ! Total bound number of particles
             '      b mass', & ! Total bound gas mass
             '   b ang mom', & ! Total bound gas angular momentum wrt CoM of entire system
             '    b tot en', & ! Total bound energy of gas
             ' ub num part', &
             '     ub mass', &
             '  ub ang mom', &
             '   ub tot en', &
             ' bt num part', & ! As in comments above, but including thermal energy in criterion
             '     bt mass', &
             '  bt ang mom', &
             '   bt tot en', &
             'ubt num part', &
             '    ubt mass', &
             ' ubt ang mom', &
             '  ubt tot en', &
             ' be num part', &
             '     be mass', &
             '  be ang mom', &
             '   be tot en', &
             'ube num part', &
             '    ube mass', &
             ' ube ang mom', &
             '  ube tot en', &
             '     HeII bm', & ! Bound mass including recombination energy of HeII
             ' HeII+HeI bm', & ! Bound mass including recombination energy of HeII, HeI
             '    He+HI bm', & ! Bound mass including recombination energy of HeII, HeI, HI
             ' He+HI+H2 bm'/)  ! Bound mass including recombination energy of HeII, HeI, HI, H2

 Zfrac = 0.
 if (dump_number == 0) then
    if (ieos /= 10 .and. ieos /= 20) then ! For MESA EoS, just use X_in and Z_in from eos module
       Xfrac = 0.69843
       Zfrac = 0.01426
       call prompt('Enter hydrogen mass fraction to assume for recombination:',Xfrac,0.,1.)
       call prompt('Enter metallicity to assume for recombination:',Zfrac,0.,1.)
    else
       Xfrac = X_in
       Zfrac = Z_in
    endif
    Yfrac = 1. - Xfrac - Zfrac
 endif

 ! Ionisation energies per particle (in code units)
 E_H2   = 0.5*Xfrac*0.0022866 * particlemass
 E_HI   = Xfrac*0.0068808 * particlemass
 E_HeI  = 0.25*Yfrac*0.012442 * particlemass
 E_HeII = 0.25*Yfrac*0.027536 * particlemass

 do i = 1,npart
    if (.not. isdead_or_accreted(xyzh(4,i))) then
       call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)
       call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,dum1,dum2,dum3,phii)
       rhopart = rhoh(xyzh(4,i), particlemass)
       tempi = eos_vars(itemp,i)
       call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
       call cross_product3D(xyzh(1:3,i), particlemass * vxyzu(1:3,i), rcrossmv)  ! Angular momentum w.r.t. CoM
       call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,tempi,gamma,ethi)
       etoti = ekini + epoti + ethi ! Overwrite etoti outputted by calc_gas_energies to use ethi instead of einti
    else
       ! Output 0 for quantities pertaining to accreted particles
       etoti   = 0.
       epoti   = 0.
       ekini   = 0.
       einti   = 0.
       ethi    = 0.
       phii    = 0.
       ponrhoi = 0.
       rcrossmv = (/ 0., 0., 0. /)
    endif

    ! Bound criterion
    if ((epoti + ekini < 0.) .or. isdead_or_accreted(xyzh(4,i))) then
       bound_i = ib
    else
       bound_i = ib + 4 ! Unbound
    endif

    bound(bound_i)     = bound(bound_i)     + 1
    bound(bound_i + 1) = bound(bound_i + 1) + particlemass
    bound(bound_i + 2) = bound(bound_i + 2) + distance(rcrossmv)
    bound(bound_i + 3) = bound(bound_i + 3) + etoti

    ! Bound criterion INCLUDING thermal energy
    if ((epoti + ekini + ethi < 0.) .or. isdead_or_accreted(xyzh(4,i))) then
       bound_i = ibt
    else
       bound_i = ibt + 4
    endif

    bound(bound_i)     = bound(bound_i)     + 1
    bound(bound_i + 1) = bound(bound_i + 1) + particlemass
    bound(bound_i + 2) = bound(bound_i + 2) + distance(rcrossmv)
    bound(bound_i + 3) = bound(bound_i + 3) + etoti

    ! Bound criterion using enthalpy
    if ((epoti + ekini + ethi + ponrhoi*particlemass < 0.)  .or. isdead_or_accreted(xyzh(4,i))) then
       bound_i = ibe
    else
       bound_i = ibe + 4
    endif

    bound(bound_i)     = bound(bound_i)     + 1
    bound(bound_i + 1) = bound(bound_i + 1) + particlemass
    bound(bound_i + 2) = bound(bound_i + 2) + distance(rcrossmv)
    bound(bound_i + 3) = bound(bound_i + 3) + etoti

    ! Bound criterion including HeI + HeII ionisation energy
    if ((epoti + ekini + ethi + E_HeII < 0.)  .or. isdead_or_accreted(xyzh(4,i))) then
       bound(25) = bound(25) + particlemass
    endif

    ! Bound criterion including HeI + HeII ionisation energy
    if ((epoti + ekini + ethi + E_HeII + E_HeI < 0.)  .or. isdead_or_accreted(xyzh(4,i))) then
       bound(26) = bound(26) + particlemass
    endif

    ! Bound criterion including HeI + HeII + HI ionisation energy
    if ((epoti + ekini + ethi + E_HeII + E_HeI + E_HI < 0.)  .or. isdead_or_accreted(xyzh(4,i))) then
       bound(27) = bound(27) + particlemass
    endif

    ! Bound criterion including HeI + HeII + HI + H2 ionisation energy
    if ((epoti + ekini + ethi + E_HeII + E_HeI + E_HI + E_H2 < 0.)  .or. isdead_or_accreted(xyzh(4,i))) then
       bound(28) = bound(28) + particlemass
    endif
 enddo

 call write_time_file('boundunbound_vs_time', columns, time, bound, ncols, dump_number)
 deallocate(columns)

end subroutine bound_mass


!----------------------------------------------------------------
!+
!  Calculate energies
!+
!----------------------------------------------------------------
subroutine calculate_energies(time,npart,particlemass,xyzh,vxyzu)
 use vectorutils, only:cross_product3D
 integer, intent(in)            :: npart
 real, intent(in)               :: time,particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 real                           :: etoti,ekini,einti,epoti,phii,phii1,jz,fxi,fyi,fzi
 real                           :: rhopart,ponrhoi,spsoundi,tempi,r_ij,radvel
 real, dimension(3)             :: rcrossmv
 character(len=17), allocatable :: columns(:)
 integer                        :: i, j, ncols
 logical                        :: inearsink
 integer, parameter             :: ie_tot        = 1
 integer, parameter             :: ie_pot        = ie_tot + 1
 integer, parameter             :: ie_kin        = ie_pot + 1
 integer, parameter             :: ie_therm      = ie_kin + 1
 integer, parameter             :: ipot_sink     = ie_therm + 1
 integer, parameter             :: ikin_sink     = ipot_sink + 1
 integer, parameter             :: iorb_sink     = ikin_sink + 1
 integer, parameter             :: iorb_comp     = iorb_sink + 1
 integer, parameter             :: ipot_env      = iorb_comp + 1
 integer, parameter             :: ie_env        = ipot_env + 1
 integer, parameter             :: ikin_bound    = ie_env + 1
 integer, parameter             :: ikin_unbound  = ikin_bound + 1
 integer, parameter             :: imass_bound   = ikin_unbound + 1
 integer, parameter             :: imass_unbound = imass_bound + 1
 integer, parameter             :: ipot_pp       = imass_unbound + 1
 integer, parameter             :: ipot_ps       = ipot_pp + 1
 integer, parameter             :: ijz_tot       = ipot_ps + 1
 integer, parameter             :: ijz_bound     = ijz_tot + 1
 integer, parameter             :: ijz_unbound   = ijz_bound + 1
 integer, parameter             :: ijz_orb       = ijz_unbound + 1
 integer, parameter             :: ie_gas        = ijz_orb + 1
 integer, parameter             :: fallbackmass  = ie_gas + 1
 integer, parameter             :: fallbackmom   = fallbackmass + 1
 real, dimension(fallbackmom)   :: encomp

 ncols = 23
 allocate(columns(ncols))
 columns = (/'total energy',&
             '  pot energy',&
             '  kin energy',&
             'therm energy',&
             '    sink pot',&
             '    sink kin',&
             '    sink orb',&
             '    comp orb',&
             '     env pot',&
             '  env energy',&
             '   bound kin',&
             ' unbound kin',&
             '  bound mass',&
             'unbound mass',&
             '     p-p pot',&
             '     p-s pot',&
             ' tot ang mom',&
             '   b ang mom',&
             '  ub ang mom',&
             ' orb ang mom',&
             '  gas energy',&
             '    fallback',&
             'fallback mom'/)

 encomp(5:) = 0.
 call compute_energies(time)
 ekin = 0.

 do i=1,npart
    encomp(ipot_pp)  = encomp(ipot_pp) + poten(i) ! poten already includes factor of 1/2 to correct for double counting
    encomp(ipot_env) = encomp(ipot_env) + poten(i)

    call cross_product3D(xyzh(1:3,i), particlemass * vxyzu(1:3,i), rcrossmv)
    jz = rcrossmv(3)
    encomp(ijz_tot) = encomp(ijz_tot) + jz

    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

    encomp(ipot_ps) = encomp(ipot_ps) + particlemass * phii

    phii1 = 0.
    call get_accel_sink_gas(1,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,fxi,fyi,fzi,phii1)
    encomp(ipot_env) = encomp(ipot_env) + phii1 * particlemass

    do j=1,nptmass
       if (xyzmh_ptmass(4,j) > 0.) then
          r_ij = separation(xyzmh_ptmass(1:3,j),xyzh(1:3,i))
          if (r_ij < 80.) then
             inearsink = .true.
          endif
       endif
    enddo

    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))

    if (etoti < 0) then
       encomp(ikin_bound) = encomp(ikin_bound) + ekini
       encomp(imass_bound) = encomp(imass_bound) + particlemass
       encomp(ijz_bound) = encomp(ijz_bound) + jz
       radvel = dot_product(vxyzu(1:3,i),xyzh(1:3,i)) / distance(xyzh(1:3,i))

       if (inearsink .eqv. .false.) then
          if (radvel < 0.) then
             encomp(fallbackmass) = encomp(fallbackmass) + particlemass
             encomp(fallbackmom) = encomp(fallbackmom) + particlemass * radvel
          endif
       endif

    else
       encomp(ikin_unbound) = encomp(ikin_unbound) + ekini
       encomp(imass_unbound) = encomp(imass_unbound) + particlemass
       encomp(ijz_unbound) = encomp(ijz_unbound) + jz
    endif
 enddo

 do i=1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       call cross_product3D(xyzmh_ptmass(1:3,i), xyzmh_ptmass(4,i)*vxyz_ptmass(1:3,i), rcrossmv)
       jz = rcrossmv(3)
       encomp(ijz_tot) = jz + encomp(ijz_tot)
       encomp(ijz_orb) = jz + encomp(ijz_orb)
       encomp(ikin_sink) = encomp(ikin_sink) + 0.5 * xyzmh_ptmass(4,i) * distance(vxyz_ptmass(1:3,i))**2
       if (i==2) encomp(iorb_comp) = encomp(iorb_comp) + 0.5 * xyzmh_ptmass(4,i) * distance(vxyz_ptmass(1:3,i))**2
    endif
 enddo

 do i=1,nptmass-1
    if (xyzmh_ptmass(4,i) > 0.) then
       do j=i+1,nptmass
          if (xyzmh_ptmass(4,j) > 0.) then
             r_ij = separation(xyzmh_ptmass(1:3,i),xyzmh_ptmass(1:3,j))
             encomp(ipot_sink) = encomp(ipot_sink) - xyzmh_ptmass(4,i) * xyzmh_ptmass(4,j) / r_ij
             if (i==1 .and. j==2) encomp(iorb_comp) = encomp(iorb_comp) - xyzmh_ptmass(4,i) * xyzmh_ptmass(4,j) / r_ij
          endif
       enddo
    endif
 enddo

 ekin = encomp(ikin_bound) + encomp(ikin_unbound) + encomp(ikin_sink)
 encomp(iorb_sink) = encomp(ipot_sink) + encomp(ikin_sink)
 encomp(ie_env) = encomp(ipot_env) + etherm + encomp(ikin_bound)
 epot = encomp(ipot_pp) + encomp(ipot_ps) + encomp(ipot_sink)
 etot = epot + ekin + etherm
 encomp(ie_gas) = encomp(ikin_bound) + encomp(ikin_unbound) + encomp(ipot_ps)

 encomp(ie_tot) = etot
 encomp(ie_pot) = epot
 encomp(ie_kin) = ekin
 encomp(ie_therm) = etherm

 call write_time_file('energy', columns, time, encomp, ncols, dump_number)
 deallocate(columns)

end subroutine calculate_energies


!!!!! Create profile !!!!!
subroutine create_profile(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)            :: npart, num
 real, intent(in)               :: time, particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 real, save                     :: profile_vector(3)
 integer                        :: ncols
 character(len=15)              :: name_in
 real, allocatable              :: profile(:,:)

 if (dump_number == 0) then
    profile_vector=(/1.,0.,0./)
    call prompt('Would you like simple profiles?', switch(1), .true.)
    call prompt('Choose profile vector x-component ',profile_vector(1))
    call prompt('Choose profile vector y-component ',profile_vector(2))
    call prompt('Choose profile vector z-component ',profile_vector(3))
 endif

 if (switch(1)) then
    ncols = 8
 else
    ncols = 18
 endif

 if (all(profile_vector <= tiny(profile_vector))) then
    write(*,*)'Using all particles!'
    call stellar_profile(time,ncols,particlemass,npart,xyzh,vxyzu,profile,switch(1))
    write(name_in, "(a)") 'part_profile'
 else
    write(*,*)'Profile_vector is:',profile_vector
    call stellar_profile(time,ncols,particlemass,npart,xyzh,vxyzu,profile,switch(1),profile_vector)
    write(name_in, "(a,i1,i1,i1)") 'ray_profile_',int(profile_vector(1:3))
 endif

 allocate(columns(18))
 columns = (/'      radius',&
             '  mass coord',&
             '     azimuth',&
             '     density',&
             '    velocity',&
             '   rad. vel.',&
             '    vxy tan.',&
             '       omega',& !Simple creates up to here
             ' int. energy',&
             '    pressure',&
             ' sound speed',&
             '        temp',&
             '       kappa',&
             '         mfp',&
             '      energy',&
             '    HII frac',&
             '   HeII frac',&
             '  HeIII frac'/)

 call write_file(name_in, 'profile', columns, profile, size(profile(1,:)), ncols, num)

 deallocate(profile,columns)
end subroutine create_profile


!!!!! Roche lobe values !!!!!
subroutine roche_lobe_values(time,npart,particlemass,xyzh,vxyzu)
 use vectorutils, only:cross_product3D
 integer, intent(in)            :: npart
 real, intent(in)               :: time, particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 integer                        :: i, j, nFB, nR1T, ncols
 integer, parameter             :: iRL1   = 1
 integer, parameter             :: iMRL1  = 2
 integer, parameter             :: iBMRL1 = 3
 integer, parameter             :: ijzRL1 = 4
 integer, parameter             :: iRL2   = 5
 integer, parameter             :: iMRL2  = 6
 integer, parameter             :: iBMRL2 = 7
 integer, parameter             :: ijzRL2 = 8
 integer, parameter             :: iR1    = 9
 integer, parameter             :: iR1T   = 10
 integer, parameter             :: iRej   = 11
 integer, parameter             :: iMej   = 12
 integer, parameter             :: iBMej  = 13
 integer, parameter             :: ijzej  = 14
 integer, parameter             :: iBjzej = 15
 integer, parameter             :: iMF    = 16
 integer, parameter             :: ijzMF  = 17
 integer, parameter             :: iDR    = 18
 integer, parameter             :: iFB    = 19
 integer, parameter             :: iFBV   = 20
 integer, parameter             :: iFBJz  = 21
 real, dimension(iFBJz)         :: MRL
 real                           :: etoti, ekini, einti, epoti, phii, jz
 logical, dimension(:), allocatable, save:: transferred
 real, save                     :: m1, m2
 real                           :: sep, sep1, sep2
 real                           :: rhovol, rhomass, rhopart, R1, rad_vel, sepCoO
 real                           :: temp_const, ponrhoi, spsoundi, tempi
 real, dimension(3)             :: rcrossmv, CoO, com_xyz, com_vxyz
 real, allocatable              :: xyz_a(:,:)
 integer                        :: npart_a, mean_rad_num
 integer, allocatable           :: iorder(:)

 allocate(iorder(npart),xyz_a(3,npart))

 MRL = 0.
 rhovol = 0.
 rhomass = 0.
 nFB = 0
 nR1T = 0
 temp_const = (unit_pressure / unit_density) * 1.34 / Rg

 if (dump_number == 0) then
    m1 = npart * particlemass + xyzmh_ptmass(4,1)
    m2 = xyzmh_ptmass(4,2)
    allocate(transferred(npart))
    transferred(1:npart) = .false.

    rho_surface = rhoh(xyzh(4,1), particlemass)
    do i=1,npart
       rhopart = rhoh(xyzh(4,i), particlemass)
       if (rhopart < rho_surface) then
          rho_surface = rhopart
       endif
    enddo
 endif

 mean_rad_num = npart / 200
 npart_a = 0

 do i=1,npart
    rhopart = rhoh(xyzh(4,i), particlemass)
    if (rhopart > rho_surface) then
       if (separation(xyzh(1:3,i), xyzmh_ptmass(1:3,1)) < &
              separation(xyzh(1:3,i), xyzmh_ptmass(1:3,2))) then
          rhomass = rhomass + particlemass
          rhovol = rhovol + particlemass / rhopart
          npart_a = npart_a + 1
          xyz_a(1:3,npart_a) = xyzh(1:3,i)
       endif
    endif
 enddo

 call set_r2func_origin(xyzmh_ptmass(1,1),xyzmh_ptmass(2,1),xyzmh_ptmass(3,1))
 call indexxfunc(npart_a,r2func_origin,xyz_a,iorder)

 R1 = 0
 do i=npart_a-mean_rad_num,npart_a
    j = iorder(i)
    R1 = R1 + separation(xyz_a(1:3,j),xyzmh_ptmass(1:3,1))
 enddo

 R1 = R1 / real(mean_rad_num)

 sep = separation(xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,2))
 MRL(iRL1) = Rochelobe_estimate(m2,m1,sep)
 MRL(iRL2) = Rochelobe_estimate(m1,m2,sep)

 !R1 = (3. * rhovol/(4. * pi))**(1./3.)
 CoO(1:3) = (xyzmh_ptmass(1:3,1) + xyzmh_ptmass(1:3,2)) / 2.
 MRL(iR1) = R1
 MRL(iRej) = separation(CoO(1:3),xyzmh_ptmass(1:3,1)) + R1

 call orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)

 do i=1,npart
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

    sep1 = separation(xyzmh_ptmass(1:3,1),xyzh(1:3,i))
    sep2 = separation(xyzmh_ptmass(1:3,2),xyzh(1:3,i))
    sepCoO = separation(CoO(1:3),xyzh(1:3,i))

    call cross_product3D(xyzh(1:3,i)-com_xyz(1:3), particlemass * vxyzu(1:3,i), rcrossmv)
    jz = rcrossmv(3)

    if (sep1 < MRL(iRL1)) then
       MRL(iMRL1) = MRL(iMRL1) + particlemass
       MRL(ijzRL1) = MRL(ijzRL1) + jz
       if (etoti < 0) then
          MRL(iBMRL1) = MRL(iBMRL1) + particlemass
       endif
    endif

    if (sep2 < MRL(iRL2)) then
       MRL(iMRL2) = MRL(iMRL2) + particlemass
       MRL(ijzRL2) = MRL(ijzRL2) + jz

       if (transferred(i) .eqv. .false.) then
          MRL(iMF) = MRL(iMF) + particlemass
          MRL(ijzMF) = MRL(ijzMF) + jz
          transferred(i) = .true.
       endif

       if (etoti < 0) then
          MRL(iBMRL2) = MRL(iBMRL2) + particlemass
       endif
    endif

    if ((sep1 - xyzh(4,i) < R1) .and. (sep1 + xyzh(4,i) > R1)) then !!!!FIX THIS
       call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
       MRL(iR1T) = MRL(iR1T) + ponrhoi * temp_const
       nR1T = nR1T + 1
    endif

    if (sepCoO > MRL(iRej)) then
       rad_vel = dot_product(vxyzu(1:3,i),xyzh(1:3,i)) / distance(xyzh(1:3,i))

       MRL(iMej) = MRL(iMej) + particlemass
       MRL(ijzej) = MRL(ijzej) + jz

       if (etoti < 0) then
          MRL(iBMej) = MRL(iBMej) + particlemass
          MRL(iBjzej) = MRL(iBjzej) + jz
       endif

       if (rad_vel < 0) then
          MRL(iFB) = MRL(iFB) + particlemass
          MRL(iFBV) = MRL(iFBV) + rad_vel
          MRL(iFBJz) = MRL(iFBJz) + jz
          nFB = nFB + 1
       endif
    endif
 enddo

 MRL(iR1T) = MRL(iR1T) / real(nR1T)
 MRL(iFBV) = MRL(iFBV) / real(nFB)

 MRL(iMRL1) = MRL(iMRL1) + xyzmh_ptmass(4,1)
 MRL(iMRL2) = MRL(iMRL2) + xyzmh_ptmass(4,2)

 MRL(iDR) = (R1 - MRL(iRL1)) / R1

 call cross_product3D(xyzmh_ptmass(1:3,1) - com_xyz(1:3),xyzmh_ptmass(4,1) * vxyz_ptmass(1:3,1),rcrossmv)
 MRL(ijzRL1) = MRL(ijzRL1) + rcrossmv(3)

 call cross_product3D(xyzmh_ptmass(1:3,2) - com_xyz(1:3),xyzmh_ptmass(4,2) * vxyz_ptmass(1:3,2),rcrossmv)
 MRL(ijzRL2) = MRL(ijzRL2) + rcrossmv(3)

 m1 = rhomass + xyzmh_ptmass(4,1)
 m2 = MRL(iMRL2)

 ncols = 21
 allocate(columns(ncols))
 columns = (/'         RL1',&
             ' Mass in RL1',&
             '  B Mass RL1',&
             '   jz in RL1',&
             '         RL2',&
             ' Mass in RL2',&
             '  B Mass RL2',&
             '   jz in RL2',&
             '          R1',&
             '     R1 temp',&
             '    R_ejecta',&
             'Mass ejected',&
             'B Mass eject',&
             '  jz ejected',&
             '  B jz eject',&
             '   Mass flow',&
             'Mass flow jz',&
             '   R1-RL1/R1',&
             '    Fallback',&
             'Fallback vel',&
             ' Fallback Jz'/)

 call write_time_file('roche_lobes', columns, time, MRL, ncols, dump_number)
 deallocate(columns,iorder)

end subroutine roche_lobe_values

!----------------------------------------------------------------
!+
!  Star stabilisation
!+
!----------------------------------------------------------------
subroutine star_stabilisation_suite(time,npart,particlemass,xyzh,vxyzu)
 use part,   only:fxyzu
 use eos,    only:equationofstate
 integer, intent(in)            :: npart
 real, intent(in)               :: time, particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 integer                        :: i,j,k,ncols,mean_rad_num,npart_a
 integer, allocatable           :: iorder(:),iorder_a(:)
 real, allocatable              :: star_stability(:)
 real                           :: total_mass,rhovol,totvol,rhopart,virialpart,virialfluid
 real                           :: phii,ponrhoi,spsoundi,tempi,epoti,ekini,einti,etoti,totekin,totepot,virialintegral,gamma
 integer, parameter             :: ivoleqrad    = 1
 integer, parameter             :: idensrad     = 2
 integer, parameter             :: imassout     = 3
 integer, parameter             :: imassfracout = 4
 integer, parameter             :: ipartrad     = 5
 integer, parameter             :: ipart2hrad   = 6
 integer, parameter             :: ipdensrad    = 7
 integer, parameter             :: ip2hdensrad  = 8
 integer, parameter             :: ivirialpart  = 9
 integer, parameter             :: ivirialfluid = 10

 ncols = 10
 allocate(columns(ncols),star_stability(ncols),iorder(npart),iorder_a(npart))
 columns = (/'vol. eq. rad',&
             ' density rad',&
             'mass outside',&
             'frac outside',&
             '    part rad',&
             ' part 2h rad',&
             '  p dens rad',&
             'p2h dens rad',&
             'part. virial',& ! Residual of virial theorem for self-gravitating particles
             'fluid virial'/) ! Residual of virial theorem for fluid

 ! Get order of particles by distance from sink particle core
 call set_r2func_origin(xyzmh_ptmass(1,1),xyzmh_ptmass(2,1),xyzmh_ptmass(3,1))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)

 ! Get density of outermost particle in initial star dump
 if (dump_number == 0) then
    rho_surface = rhoh(xyzh(4,iorder(npart)), particlemass)
 endif

 npart_a = 0
 totvol = 0.
 rhovol = 0.
 virialpart = 0.
 totekin = 0.
 totepot = 0.
 virialintegral= 0.
 do i = 1,npart
    rhopart = rhoh(xyzh(4,i), particlemass)
    totvol = totvol + particlemass / rhopart ! Sum "volume" of all particles
    virialpart = virialpart + particlemass * ( dot_product(fxyzu(1:3,i),xyzh(1:3,i)) + dot_product(vxyzu(1:3,i),vxyzu(1:3,i)) )
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)
    totekin = totekin + ekini
    totepot = totepot + 0.5*epoti ! Factor of 1/2 to correct for double counting
    if (rhopart > rho_surface) then
       ! Sum "volume" of particles within "surface" of initial star dump
       rhovol = rhovol + particlemass / rhopart
       npart_a = npart_a + 1 ! Count number of particles within "surface" of initial star dump
    endif
    ! Calculate residual of Virial theorem for fluid
    if (ieos == 2) then
       call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,gamma_local=gamma)
    else
       call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    endif
    virialintegral = virialintegral + 3. * ponrhoi * particlemass
 enddo
 virialpart = virialpart / (abs(totepot) + 2.*abs(totekin)) ! Normalisation for the virial
 virialfluid = (virialintegral + totepot) / (abs(virialintegral) + abs(totepot))

 ! Sort particles within "surface" by radius
 call indexxfunc(npart_a,r2func_origin,xyzh,iorder_a)

 mean_rad_num = npart / 200 ! 0.5 percent of particles
 star_stability = 0.
 ! Loop over the outermost npart/200 particles that are within the "surface"
 do i = npart_a - mean_rad_num,npart_a
    j = iorder(i)
    k = iorder_a(i)
    star_stability(ipartrad)    = star_stability(ipartrad)    + separation(xyzh(1:3,j),xyzmh_ptmass(1:3,1))
    star_stability(ipart2hrad)  = star_stability(ipart2hrad)  + separation(xyzh(1:3,j),xyzmh_ptmass(1:3,1)) + xyzh(4,j)
    star_stability(ipdensrad)   = star_stability(ipdensrad)   + separation(xyzh(1:3,k),xyzmh_ptmass(1:3,1))
    star_stability(ip2hdensrad) = star_stability(ip2hdensrad) + separation(xyzh(1:3,k),xyzmh_ptmass(1:3,1)) + xyzh(4,j)
 enddo

 star_stability(ipartrad)    = star_stability(ipartrad)    / real(mean_rad_num)
 star_stability(ipart2hrad)  = star_stability(ipart2hrad)  / real(mean_rad_num)
 star_stability(ipdensrad)   = star_stability(ipdensrad)   / real(mean_rad_num)
 star_stability(ip2hdensrad) = star_stability(ip2hdensrad) / real(mean_rad_num)
 star_stability(ivoleqrad)   = (3. * totvol/(4. * pi))**(1./3.)
 star_stability(idensrad)    = (3. * rhovol/(4. * pi))**(1./3.)
 star_stability(ivirialpart) = virialpart
 star_stability(ivirialfluid)= virialfluid

 if (dump_number == 0) then
    init_radius = star_stability(ivoleqrad)
 endif

 star_stability(imassout) = 0.
 total_mass = xyzmh_ptmass(4,1)
 do i = 1,npart
    if (separation(xyzmh_ptmass(1:3,1),xyzh(1:3,i)) > init_radius) then
       star_stability(imassout) = star_stability(imassout) + particlemass
    endif
    total_mass = total_mass + particlemass
 enddo

 star_stability(imassfracout) = star_stability(imassout) / total_mass
 call write_time_file('star_stability', columns, time, star_stability, ncols, dump_number)
 deallocate(columns,star_stability,iorder,iorder_a)

end subroutine star_stabilisation_suite


!----------------------------------------------------------------
!+
!  Print simulation parameters
!+
!----------------------------------------------------------------
subroutine print_simulation_parameters(npart,particlemass)
 integer, intent(in)            :: npart
 real, intent(in)               :: particlemass
 integer                        :: i

 write(*,"(/,3(a,es10.3,1x),a)") '     Mass: ',umass,    'g       Length: ',udist,  'cm     Time: ',utime,'s'
 write(*,"(3(a,es10.3,1x),a)") '  Density: ',unit_density, 'g/cm^3  Energy: ',unit_energ,'erg    En/m: ',unit_ergg,'erg/g'
 write(*,"(3(a,es10.3,1x),a)") ' Velocity: ',unit_velocity,'cm/s    Bfield: ',unit_Bfield,'G  Pressure: ',&
                                     unit_pressure,'g/cm s^2'
 write(*,"(2(a,es10.3,1x),/)")   '        G: ', gg*umass*utime**2/udist**3,'             c: ',c*utime/udist

 do i=1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       write(*,'(A,I2,A,ES10.3,A,ES10.3)') 'Point mass ',i,': M = ',xyzmh_ptmass(4,i),' and h_soft = ',xyzmh_ptmass(ihsoft,i)
    endif
 enddo
 write(*,"(A,ES10.3)")  'Sink-sink separation: ', separation(xyzmh_ptmass(1:3,1), xyzmh_ptmass(1:3,2))

 write(*,'(A,I7,A,ES10.3)') 'Gas particles : ',npart,' particles, each of mass ',particlemass

end subroutine print_simulation_parameters


!----------------------------------------------------------------
!+
!  Write quantities (up to four) to divv file
!+
!----------------------------------------------------------------
subroutine output_divv_files(time,dumpfile,npart,particlemass,xyzh,vxyzu)
 use part,              only:eos_vars,itemp,nucleation,idK0,idK1,idK2,idK3,idJstar,idmu,idgamma
 use eos,               only:entropy
 use eos_mesa,          only:get_eos_kappa_mesa
 use mesa_microphysics, only:getvalue_mesa
 use sortutils,         only:set_r2func_origin,r2func_origin,indexxfunc
 use ionization_mod,    only:calc_thermal_energy,ionisation_fraction
 use dust_formation,    only:psat_C,eps,set_abundances,mass_per_H, chemical_equilibrium_light, calc_nucleation!, Scrit
 !use dim,     only:nElements
 integer, intent(in)          :: npart
 character(len=*), intent(in) :: dumpfile
 real, intent(in)             :: time,particlemass
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 integer                      :: i,k,Nquantities,ierr,iu
 integer, save                :: quantities_to_calculate(4)
 integer, allocatable         :: iorder(:)
 real                         :: ekini,einti,epoti,ethi,phii,rho_cgs,ponrhoi,spsoundi,tempi,&
                                 omega_orb,kappai,kappat,kappar,pgas,mu,entropyi,rhopart,&
                                 dum1,dum2,dum3,dum4,dum5
 real, allocatable, save      :: init_entropy(:)
 real, allocatable            :: quant(:,:)
 real, dimension(3)           :: com_xyz,com_vxyz,xyz_a,vxyz_a
 real                         :: pC, pC2, pC2H, pC2H2, nH_tot, epsC, S
 real                         :: taustar, taugr, JstarS
 real, parameter :: Scrit = 2. ! Critical saturation ratio
 logical :: verbose = .false.

 allocate(quant(4,npart))
 Nquantities = 13
 if (dump_number == 0) then
    print "(13(a,/))",&
           '1) Total energy (kin + pot + therm)', &
           '2) Mach number', &
           '3) Opacity from MESA tables', &
           '4) Gas omega w.r.t. effective CoM', &
           '5) Fractional difference between gas and orbital omega', &
           '6) MESA EoS specific entropy', &
           '7) Fractional entropy gain', &
           '8) Specific recombination energy', &
           '9) Total energy (kin + pot)', &
           '10) Mass coordinate', &
           '11) Gas omega w.r.t. CoM', &
           '12) Gas omega w.r.t. sink 1',&
           '13) JstarS' !option to calculate JstarS

    quantities_to_calculate = (/1,2,4,5/)
    call prompt('Choose first quantity to compute ',quantities_to_calculate(1),0,Nquantities)
    call prompt('Choose second quantity to compute ',quantities_to_calculate(2),0,Nquantities)
    call prompt('Choose third quantity to compute ',quantities_to_calculate(3),0,Nquantities)
    call prompt('Choose fourth quantity to compute ',quantities_to_calculate(4),0,Nquantities)
 endif

 ! Calculations performed outside loop over particles
 call compute_energies(time)
 omega_orb = 0.
 com_xyz = 0.
 com_vxyz = 0.
 do k=1,4
    select case (quantities_to_calculate(k))
    case(0,1,2,3,6,8,9,13) ! Nothing to do
    case(4,5,11,12) ! Fractional difference between gas and orbital omega
       if (quantities_to_calculate(k) == 4 .or. quantities_to_calculate(k) == 5) then
          com_xyz  = (xyzmh_ptmass(1:3,1)*xyzmh_ptmass(4,1) + xyzmh_ptmass(1:3,2)*xyzmh_ptmass(4,2)) &
                  / (xyzmh_ptmass(4,1) + xyzmh_ptmass(4,2))
          com_vxyz = (vxyz_ptmass(1:3,1)*xyzmh_ptmass(4,1)  + vxyz_ptmass(1:3,2)*xyzmh_ptmass(4,2))  &
                  / (xyzmh_ptmass(4,1) + xyzmh_ptmass(4,2))
       elseif (quantities_to_calculate(k) == 11 .or. quantities_to_calculate(k) == 12) then
          com_xyz = xyzmh_ptmass(1:3,1)
          com_vxyz = vxyz_ptmass(1:3,1)
       endif
       do i=1,nptmass
          xyz_a(1:3) = xyzmh_ptmass(1:3,i) - com_xyz(1:3)
          vxyz_a(1:3) = vxyz_ptmass(1:3,i) - com_vxyz(1:3)
          omega_orb = omega_orb + 0.5 * (-xyz_a(2) * vxyz_a(1) + xyz_a(1) * vxyz_a(2)) / dot_product(xyz_a(1:2), xyz_a(1:2))
       enddo
    case(7)
       if (dump_number==0) allocate(init_entropy(npart))
    case(10)
       call set_r2func_origin(0.,0.,0.)
       allocate(iorder(npart))
       call indexxfunc(npart,r2func_origin,xyzh,iorder)
       deallocate(iorder)
    case default
       print*,"Error: Requested quantity is invalid."
       stop
    end select
 enddo

 !set initial abundances to get mass_per_H
 call set_abundances
 ! Calculations performed in loop over particles
 do i=1,npart
    do k=1,4
       select case (quantities_to_calculate(k))
       case(13) !to calculate JstarS
          rhopart = rhoh(xyzh(4,i), particlemass)
          rho_cgs = rhopart*unit_density
          !call equationofstate to obtain temperature and store it in tempi
          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
          JstarS = 0.
          !nH_tot is needed to normalize JstarS
          nH_tot = rho_cgs/mass_per_H
          epsC   = eps(3) - nucleation(idK3,i)
          if (epsC < 0.) then
             print *,'eps(C) =',eps(3),', K3=',nucleation(idK3,i),', epsC=',epsC,', T=',tempi,' rho=',rho_cgs
             print *,'JKmuS=',nucleation(:,i)
             stop '[S-dust_formation] epsC < 0!'
          endif
          if (tempi > 450.) then
             !call chemical_equilibrium_light to obtain pC, and pC2H2
             call chemical_equilibrium_light(rho_cgs, tempi, epsC, pC, pC2, pC2H, pC2H2, nucleation(idmu,i), nucleation(idgamma,i))
             S = pC/psat_C(tempi)
             if (S > Scrit) then
                !call nucleation_function to obtain JstarS
                call calc_nucleation(tempi, pC, 0., 0., 0., pC2H2, S, JstarS, taustar, taugr)
                JstarS = JstarS/ nH_tot
             endif
          endif
          !Check if the variables have meaningful values close to condensation temperatures
          if (tempi >= 1400. .and. tempi <= 1500. .and. verbose ) then
             print *,'size(nucleation,1) = ',size(nucleation,1)
             print *,'size(nucleation,2) = ',size(nucleation,2)
             print *,'nucleation(idK3,i) = ',nucleation(idK3,i)
             print *,'epsC = ',epsC
             print *,'tempi = ',tempi
             print *,'S = ',S
             print *,'pC =',pC
             print *,'psat_C(tempi) = ',psat_C(tempi)
             print *,'nucleation(idmu,i) = ',nucleation(idmu,i)
             print *,'nucleation(idgamma,i) = ',nucleation(idgamma,i)
             print *,'taustar = ',taustar
             print *,'eps = ',eps
             print *,'JstarS = ',JstarS
          endif
          quant(k,i) = JstarS

       case(0) ! Skip
          quant(k,i) = 0.

       case(1,9) ! Total energy (kin + pot + therm)
          rhopart = rhoh(xyzh(4,i), particlemass)
          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
          call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum1)
          if (quantities_to_calculate(k)==1) then
             call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)
             quant(k,i) = (ekini + epoti + ethi) / particlemass ! Specific energy
          elseif (quantities_to_calculate(k)==9) then
             quant(k,i) = (ekini + epoti) / particlemass ! Specific energy
          endif

       case(2) ! Mach number
          rhopart = rhoh(xyzh(4,i), particlemass)
          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
          quant(k,i) = distance(vxyzu(1:3,i)) / spsoundi

       case(3) ! Opacity from MESA tables
          rhopart = rhoh(xyzh(4,i), particlemass)
          call ionisation_fraction(rhopart*unit_density,eos_vars(itemp,i),X_in,1.-X_in-Z_in,dum1,dum2,dum3,dum4,dum5)
          if (ieos == 10) then
             call get_eos_kappa_mesa(rhopart*unit_density,eos_vars(itemp,i),kappai,kappat,kappar)
             quant(k,i) = kappai
          else
             quant(k,i) = 0.
          endif

       case(4,11,12) ! Gas omega w.r.t. effective CoM
          xyz_a  = xyzh(1:3,i)  - com_xyz(1:3)
          vxyz_a = vxyzu(1:3,i) - com_vxyz(1:3)
          quant(k,i) = (-xyz_a(2) * vxyz_a(1) + xyz_a(1) * vxyz_a(2)) / dot_product(xyz_a(1:2), xyz_a(1:2))

       case(5) ! Fractional difference between gas and orbital omega
          xyz_a  = xyzh(1:3,i)  - com_xyz(1:3)
          vxyz_a = vxyzu(1:3,i) - com_vxyz(1:3)
          quant(k,i) = (-xyz_a(2) * vxyz_a(1) + xyz_a(1) * vxyz_a(2)) / dot_product(xyz_a(1:2), xyz_a(1:2))
          quant(k,i) = (quant(k,i) - omega_orb) / omega_orb

       case(6,7) ! Calculate MESA EoS entropy
          entropyi = 0.
          rhopart = rhoh(xyzh(4,i), particlemass)
          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
          if (ieos==10) then
             call getvalue_mesa(rhopart*unit_density,vxyzu(4,i)*unit_ergg,3,pgas,ierr) ! Get gas pressure
             mu = rhopart*unit_density * Rg * eos_vars(itemp,i) / pgas
             entropyi = entropy(rhopart*unit_density,ponrhoi*rhopart*unit_pressure,mu,3,vxyzu(4,i)*unit_ergg,ierr)
          elseif (ieos==2) then
             entropyi = entropy(rhopart*unit_density,ponrhoi*rhopart*unit_pressure,gmw,1)
          endif

          if (quantities_to_calculate(k) == 7) then
             if (dump_number == 0) then
                init_entropy(i) = entropyi ! Store initial entropy on each particle
             endif
             quant(k,i) = entropyi/init_entropy(i) - 1.
          elseif (quantities_to_calculate(k) == 6) then
             quant(k,i) = entropyi
          endif

       case(8) ! Specific recombination energy
          rhopart = rhoh(xyzh(4,i), particlemass)
          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
          call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)
          quant(k,i) = vxyzu(4,i) - ethi / particlemass ! Specific energy

       case(10) ! Mass coordinate
          quant(k,iorder(i)) = real(i,kind=kind(time)) * particlemass

       case default
          print*,"Error: Requested quantity is invalid."
          stop
       end select
    enddo
 enddo

 open(newunit=iu,file=trim(dumpfile)//".divv",status='replace',form='unformatted')
 do k=1,4
    write(iu) (quant(k,i),i=1,npart)
 enddo
 close(iu)
 deallocate(quant)

end subroutine output_divv_files



!!!!! EoS surfaces !!!!!
subroutine eos_surfaces
 integer :: i, j, ierr
 real    :: rho_array(1000) = (/ (10**(i/10000.), i=-180000,-30150,150) /)
 real    :: eni_array(1000) = (/ (10**(i/10000.), i=120000,149970,30) /)
 real    :: temp_array(400) = (/ (10**(i/1000.), i=3000,6990,10) /)
 real    :: kappa_array(1000,400)
 real    :: gam1_array(1000,1000)
 real    :: pres_array(1000,1000)
 real    :: dum(1000,1000)
 real    :: kappat, kappar


 do i=1,size(rho_array)
    do j=1,size(eni_array)
       if (j < size(temp_array) + 1) then
          call get_eos_kappa_mesa(rho_array(i),temp_array(j),kappa_array(i,j),kappat,kappar)
       endif
       call get_eos_pressure_temp_gamma1_mesa(rho_array(i),eni_array(j),pres_array(i,j),dum(i,j),gam1_array(i,j),ierr)
       !call get_eos_pressure_temp_mesa(rho_array(i),eni_array(j),pres_array(i,j),temp)
       !pres_array(i,j) = eni_array(j)*rho_array(i)*0.66667 / pres_array(i,j)
    enddo
 enddo

 open(unit=1000, file='mesa_eos_pressure.out', status='replace')

 !Write data to file
 do i=1,1000
    write(1000,"(1000(3x,es18.11e2,1x))") pres_array(i,:)
 enddo

 close(unit=1000)

 open(unit=1002, file='mesa_eos_gamma.out', status='replace')

 !Write data to file
 do i=1,1000
    write(1002,"(1000(3x,es18.11e2,1x))") gam1_array(i,:)
 enddo

 close(unit=1002)

 open(unit=1001, file='mesa_eos_kappa.out', status='replace')

 !Write data to file
 do i=1,1000
    write(1001,"(400(3x,es18.11e2,1x))") kappa_array(i,:)
 enddo

 close(unit=1001)

end subroutine eos_surfaces


!----------------------------------------------------------------
!+
!  Particle tracker: Paint the life of a particle
!+
!----------------------------------------------------------------
subroutine track_particle(time,particlemass,xyzh,vxyzu)
 use part, only:eos_vars,itemp
 use eos,  only:entropy
 use mesa_microphysics, only:getvalue_mesa
 use ionization_mod, only:calc_thermal_energy,ionisation_fraction
 real, intent(in)        :: time,particlemass
 real, intent(inout)     :: xyzh(:,:),vxyzu(:,:)
 integer, parameter      :: nparttotrack=10,ncols=17
 real                    :: r,v,rhopart,ponrhoi,Si,spsoundi,tempi,machi,xh0,xh1,xhe0,xhe1,xhe2,&
                            ekini,einti,epoti,ethi,etoti,dum,phii,pgas,mu
 real, dimension(ncols)  :: datatable
 character(len=17)       :: filenames(nparttotrack),columns(ncols)
 integer                 :: i,k,partID(nparttotrack),ientropy,ierr

 partID = (/ 1,2,3,4,5,6,7,8,9,10 /)
 columns = (/ '      r',&
              '      v',&
              '    rho',&
              '   temp',&
              'entropy',&
              'spsound',&
              '   mach',&
              '   ekin',&
              '   epot',&
              '    eth',&
              '   eint',&
              '   etot',&
              '    xHI',&
              '   xHII',&
              '   xHeI',&
              '  xHeII',&
              ' xHeIII' /)

 call compute_energies(time)

 do i=1,nparttotrack
    write (filenames(i),"(a1,i7.7)") "p", partID(i)
 enddo

 do k=1,nparttotrack
    i = partID(k)
    r = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
    v = separation(vxyzu(1:3,i),vxyz_ptmass(1:3,1))
    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    machi = v / spsoundi
    select case(ieos)
    case(2)
       ientropy = 1
    case(10,12)
       ientropy = 2
    case default
       ientropy = -1
    end select
    if (ieos==10) then
       call getvalue_mesa(rhopart*unit_density,vxyzu(4,i)*unit_ergg,3,pgas,ierr) ! Get gas pressure
       mu = rhopart*unit_density * Rg * eos_vars(itemp,i) / pgas
    else
       mu = gmw
    endif
    ! MESA ENTROPY
    Si = 0.
    if (ieos==10) then
       Si = entropy(rhopart*unit_density,ponrhoi*rhopart*unit_pressure,mu,3,vxyzu(4,i)*unit_ergg,ierr)
    endif
    ! MESA ENTROPY
    !  Si = entropy(rhopart*unit_density,ponrhoi*rhopart*unit_pressure,mu,ientropy,vxyzu(4,i)*unit_ergg,ierr)
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
    call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)
    etoti = ekini + epoti + ethi
    call ionisation_fraction(rhopart*unit_density,eos_vars(itemp,i),X_in,1.-X_in-Z_in,xh0,xh1,xhe0,xhe1,xhe2)

    ! Write file
    datatable = (/ r,v,rhopart,eos_vars(itemp,i),Si,spsoundi,machi,ekini,epoti,ethi,einti,etoti,xh0,xh1,xhe0,xhe1,xhe2 /)
    call write_time_file(trim(adjustl(filenames(k))),columns,time,datatable,ncols,dump_number)
 enddo

end subroutine track_particle


!----------------------------------------------------------------
!+
!  Optical depth profile
!+
!----------------------------------------------------------------
subroutine tau_profile(time,num,npart,particlemass,xyzh)
 use part, only:eos_vars,itemp
 integer, intent(in)    :: npart,num
 real, intent(in)       :: time,particlemass
 real, intent(inout)    :: xyzh(:,:)
 integer                :: nbins
 real, allocatable      :: rad_part(:),kappa_part(:),rho_part(:)
 real, allocatable      :: kappa_hist(:),rho_hist(:),tau_r(:),sepbins(:)
 real                   :: maxloga,minloga,kappa,kappat,kappar
 character(len=17)      :: filename
 character(len=40)      :: data_formatter
 integer                :: i,unitnum

 call compute_energies(time)
 nbins      = 500

 allocate(rad_part(npart),kappa_part(npart),rho_part(npart))
 rad_part   = 0.
 kappa_part = 0.
 rho_part   = 0.
 minloga    = 0.5
 maxloga    = 4.3

 allocate(rho_hist(nbins),kappa_hist(nbins),sepbins(nbins),tau_r(nbins))
 filename = '      grid_tau.ev'

 do i=1,npart
    rho_part(i) = rhoh(xyzh(4,i), particlemass)
    rad_part(i) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
    call get_eos_kappa_mesa(rho_part(i)*unit_density,eos_vars(itemp,i),kappa,kappat,kappar)
    kappa_part(i) = kappa ! In cgs units?
 enddo

 call histogram_setup(rad_part(1:npart),kappa_part,kappa_hist,npart,maxloga,minloga,nbins,.true.,.true.)
 call histogram_setup(rad_part(1:npart),rho_part,rho_hist,npart,maxloga,minloga,nbins,.true.,.true.)


 ! Integrate optical depth inwards
 sepbins = (/ (10**(minloga + (i-1) * (maxloga-minloga)/real(nbins)), i=1,nbins) /) ! Create log-uniform bins
 ! Convert to cgs units (kappa has already been outputted in cgs)
 rho_hist = rho_hist * unit_density
 sepbins = sepbins * udist ! udist should be Rsun in cm

 tau_r(nbins) = 0.
 do i=nbins,2,-1
    tau_r(i-1) = tau_r(i) + kappa_hist(i) * rho_hist(i) * (sepbins(i+1) - sepbins(i))
 enddo

 ! Write data row
 write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
 if (num == 0) then
    unitnum = 1000
    open(unit=unitnum, file=trim(adjustl(filename)), status='replace')
    write(unitnum, "(a)") '# Optical depth profile'
    close(unit=unitnum)
 endif
 unitnum=1002
 open(unit=unitnum, file=trim(adjustl(filename)), position='append')
 write(unitnum,data_formatter) time,tau_r
 close(unit=unitnum)
 deallocate(rad_part,kappa_part,rho_part)
 deallocate(rho_hist,kappa_hist,sepbins,tau_r)

end subroutine tau_profile

!----------------------------------------------------------------
!+
!  Sound crossing time profile
!+
!----------------------------------------------------------------
subroutine tconv_profile(time,num,npart,particlemass,xyzh,vxyzu)
 use part,  only:itemp
 use eos,   only:get_spsound
 use units, only:unit_velocity
 integer, intent(in)    :: npart,num
 real, intent(in)       :: time,particlemass
 real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                :: nbins
 real, allocatable      :: rad_part(:),cs_part(:)
 real, allocatable      :: cs_hist(:),tconv(:),sepbins(:)
 real                   :: maxloga,minloga,rhoi
 character(len=17)      :: filename
 character(len=40)      :: data_formatter
 integer                :: i,unitnum

 call compute_energies(time)
 nbins      = 500
 allocate(rad_part(npart),cs_part(npart))
 rad_part   = 0.
 cs_part   = 0.
 minloga    = 0.5
 maxloga    = 4.3

 allocate(cs_hist(nbins),sepbins(nbins),tconv(nbins))
 filename = '    grid_tconv.ev'

 do i=1,npart
    rhoi = rhoh(xyzh(4,i), particlemass)
    rad_part(i) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
    cs_part(i) = get_spsound(eos_type=ieos,xyzi=xyzh(:,i),rhoi=rhoi,vxyzui=vxyzu(:,i),gammai=gamma,mui=gmw,Xi=X_in,Zi=Z_in)
 enddo

 call histogram_setup(rad_part(1:npart),cs_part,cs_hist,npart,maxloga,minloga,nbins,.true.,.true.)

 ! Integrate sound-crossing time from surface inwards
 sepbins = (/ (10**(minloga + (i-1) * (maxloga-minloga)/real(nbins)), i=1,nbins) /) ! Create log-uniform bins
 ! Convert to cgs units
 cs_hist = cs_hist * unit_velocity
 sepbins = sepbins * udist ! udist should be Rsun in cm

 tconv(nbins) = 0.
 do i=nbins,2,-1
    if (cs_hist(i) < tiny(1.)) then
       tconv(i-1) = tconv(i)
    else
       tconv(i-1) = tconv(i) + (sepbins(i+1) - sepbins(i)) / cs_hist(i)
    endif
 enddo

 ! Write data row
 write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
 if (num == 0) then
    unitnum = 1000
    open(unit=unitnum, file=trim(adjustl(filename)), status='replace')
    write(unitnum, "(a)") '# Sound crossing time profile'
    close(unit=unitnum)
 endif
 unitnum=1002
 open(unit=unitnum, file=trim(adjustl(filename)), position='append')
 write(unitnum,data_formatter) time,tconv
 close(unit=unitnum)

 deallocate(rad_part,cs_part)

end subroutine tconv_profile


!----------------------------------------------------------------
!+
!  Histogram of optical depth at hydrogen recombination
!+
!----------------------------------------------------------------
subroutine recombination_tau(time,npart,particlemass,xyzh,vxyzu)
 use part, only:eos_vars,itemp
 use ionization_mod, only:calc_thermal_energy,ionisation_fraction
 integer, intent(in)    :: npart
 real,    intent(in)    :: time,particlemass
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: nbins
 integer, allocatable   :: recombined_pid(:)
 real, allocatable      :: rad_part(:),kappa_part(:),rho_part(:)
 real, allocatable, save:: tau_recombined(:)
 real, allocatable      :: kappa_hist(:),rho_hist(:),tau_r(:),sepbins(:),sepbins_cm(:)
 logical, allocatable, save :: prev_recombined(:)
 real                   :: maxloga,minloga,kappa,kappat,kappar,xh0,xh1,xhe0,xhe1,xhe2,&
                           ponrhoi,spsoundi,tempi,etoti,ekini,einti,epoti,ethi,phii,dum
 real, parameter        :: recomb_th=0.9
 integer                :: i,j,nrecombined,bin_ind

 call compute_energies(time)
 allocate(rad_part(npart),kappa_part(npart),rho_part(npart),recombined_pid(npart))
 rad_part   = 0.
 kappa_part = 0.
 rho_part   = 0.
 nbins      = 300 ! Number of radial bins
 minloga    = 0.5
 maxloga    = 4.3
 allocate(rho_hist(nbins),kappa_hist(nbins),sepbins(nbins),sepbins_cm(nbins),tau_r(nbins))
 if (dump_number == 0) then
    allocate(tau_recombined(npart),prev_recombined(npart))
    tau_recombined = -1. ! Store tau of newly-recombined particles. -ve data means particle never recombined]
    prev_recombined = .false. ! All hydrogen is ionised at the start
 endif

 j=0
 do i=1,npart
    rho_part(i) = rhoh(xyzh(4,i), particlemass)
    rad_part(i) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
    call equationofstate(ieos,ponrhoi,spsoundi,rho_part(i),xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    call get_eos_kappa_mesa(rho_part(i)*unit_density,eos_vars(itemp,i),kappa,kappat,kappar)
    kappa_part(i) = kappa ! In cgs units
    call ionisation_fraction(rho_part(i)*unit_density,eos_vars(itemp,i),X_in,1.-X_in-Z_in,xh0,xh1,xhe0,xhe1,xhe2)
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum) ! Calculate total energy
    call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rho_part(i),eos_vars(itemp,i),gamma,ethi)
    etoti = ekini + epoti + ethi
    if ((xh0 > recomb_th) .and. (.not. prev_recombined(i)) .and. (etoti < 0.)) then ! Recombination event and particle is still bound
       j=j+1
       recombined_pid(j) = i
       prev_recombined(i) = .true.
    else
       prev_recombined(i) = .false.
    endif
 enddo
 nrecombined = j

 call histogram_setup(rad_part(1:npart),kappa_part,kappa_hist,npart,maxloga,minloga,nbins,.true.,.true.)
 call histogram_setup(rad_part(1:npart),rho_part,rho_hist,npart,maxloga,minloga,nbins,.true.,.true.)

 ! Integrate optical depth inwards
 sepbins = (/ (10.**(minloga + (i-1) * (maxloga-minloga)/real(nbins)), i=1,nbins) /) ! Create log-uniform bins

 ! Convert to cgs units (kappa has already been outputted in cgs)
 rho_hist = rho_hist * unit_density
 sepbins_cm = sepbins * udist ! udist should be Rsun in g

 ! Integrate bins in tau(r)
 tau_r(nbins) = 0.
 do i=nbins,2,-1
    tau_r(i-1) = tau_r(i) + kappa_hist(i) * rho_hist(i) * (sepbins_cm(i+1) - sepbins_cm(i))
 enddo

 ! Integrate optical depth for each newly recombined particle
 do j=1,nrecombined
    i = recombined_pid(j)
    bin_ind = 1 + nint( nbins * ( log10(rad_part(i))-minloga ) / (maxloga-minloga) )   ! Find radial bin of recombined particle
    tau_recombined(i) = tau_r(bin_ind)
 enddo
 ! Trick write_time_file into writing my data table
 if (dump_number == 320) then
    do i=1,npart
       call write_time_file("recombination_tau",(/'          tau'/),-1.,tau_recombined(i),1,i-1) ! Set num = i-1 so that header will be written for particle 1 and particle 1 only
    enddo
 endif
 deallocate(recombined_pid,rad_part,kappa_part,rho_part)

end subroutine recombination_tau


!----------------------------------------------------------------
!+
!  Energy histogram
!+
!----------------------------------------------------------------
subroutine energy_hist(time,npart,particlemass,xyzh,vxyzu)
 use part, only:eos_vars,itemp
 use ionization_mod, only:calc_thermal_energy
 integer, intent(in)            :: npart
 real, intent(in)               :: time,particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: filename(:)
 character(len=40)              :: data_formatter
 integer                        :: nbins,nhists,i,unitnum
 real, allocatable              :: hist(:),coord(:,:),Emin(:),Emax(:)
 real                           :: rhopart,ponrhoi,spsoundi,tempi,phii,epoti,ekini,einti,ethi,dum
 real, allocatable              :: quant(:)
 logical                        :: ilogbins

 nhists = 3
 nbins = 500
 allocate(filename(nhists),coord(npart,nhists),hist(nbins),Emin(nhists),Emax(nhists))
 Emin = (/ -0.0446, 0., 0. /)
 Emax = (/ 0.0315, 0.0105, 0.0105 /)
 ilogbins = .false.
 filename = (/ '       hist_kp.ev', &
               '     hist_erec.ev', &
               '      hist_eth.ev' /)

 allocate(quant(npart))
 quant = (/ (1., i=1,npart) /)
 do i=1,npart
    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
    if (ieos==10 .or. ieos==20) then
       call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)
    else
       ethi = einti
    endif
    coord(i,1) = (ekini + epoti)/particlemass
    coord(i,2) = vxyzu(4,i) - ethi/particlemass
    coord(i,3) = ethi/particlemass
 enddo

 write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
 do i=1,nhists
    call histogram_setup(coord(:,i),quant,hist,npart,Emax(i),Emin(i),nbins,.false.,ilogbins)
    if (dump_number == 0) then
       unitnum = 1000
       open(unit=unitnum, file=trim(adjustl(filename(i))), status='replace')
       close(unit=unitnum)
    endif
    unitnum=1001+i
    open(unit=unitnum, file=trim(adjustl(filename(i))), status='old', position='append')
    write(unitnum,data_formatter) time,hist
    close(unit=unitnum)
 enddo
 deallocate(filename,coord,hist,Emin,Emax,quant)

end subroutine energy_hist


!----------------------------------------------------------------
!+
!  Energy profile
!+
!----------------------------------------------------------------
subroutine energy_profile(time,npart,particlemass,xyzh,vxyzu)
 use part,              only:eos_vars,itemp
 use eos,               only:entropy
 use mesa_microphysics, only:getvalue_mesa
 use ionization_mod,    only:calc_thermal_energy,ionisation_fraction
 integer, intent(in)    :: npart
 real, intent(in)       :: time,particlemass
 real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                :: nbins
 real, allocatable      :: coord(:)
 real, allocatable      :: hist(:),quant(:,:)
 real                   :: ekini,einti,epoti,ethi,phii,pgas,mu,dum,rhopart,ponrhoi,spsoundi,tempi,&
                           maxcoord,mincoord,xh0,xh1,xhe0,xhe1,xhe2
 character(len=17), allocatable :: filename(:),headerline(:)
 character(len=40)      :: data_formatter
 integer                :: i,k,unitnum,ierr,ientropy,nvars
 integer, allocatable   :: iorder(:)
 integer, save          :: iquantity
 logical                :: ilogbins
 logical, save          :: use_mass_coord

 if (dump_number==0) then
    iquantity = 1
    use_mass_coord = .false.
    print "(4(/,a))",'1. Energy',&
                     '2. Entropy',&
                     '3. Bernoulli energy',&
                     '4. Ion fractions'
    call prompt("Select quantity to calculate",iquantity,1,4)
    call prompt("Bin in mass coordinates instead of radius?",use_mass_coord)
 endif

 nbins = 500
 allocate(hist(nbins))
 if (use_mass_coord) then
    mincoord  = 3.8405  ! Min. mass coordinate
    maxcoord  = 12.0 ! Max. mass coordinate
    ilogbins = .false.
 else
    mincoord  = 0.5  ! Min. log(r)
    maxcoord  = 4.3  ! Max. log(r)
    ilogbins = .true.
 endif

 call compute_energies(time)

 ! Allocate arrays for single variable outputs
 if ( (iquantity==1) .or. (iquantity==2) .or. (iquantity==3) ) then
    nvars = 1
 else
    nvars = 5
 endif
 allocate(filename(nvars),headerline(nvars),quant(npart,nvars),coord(npart))

 coord = 0.
 quant = 0.
 select case (iquantity)
 case(1) ! Energy
    filename = '     grid_Etot.ev'
    headerline = '# Energy profile '
 case(2) ! Entropy
    filename = '  grid_entropy.ev'
    headerline = '# Entropy profile'
    select case(ieos)
    case(2)
       ientropy = 1
    case(12)
       ientropy = 2
    case(10,20)
       ientropy = 3
    case default
       ientropy = -1
    end select
 case(3) ! Bernoulli energy (per unit mass)
    filename = 'grid_bernoulli.ev'
    headerline = '# Bernoulli prof.'
 case(4) ! Ion fraction profiles
    filename = (/ '       grid_HI.ev', &
                  '      grid_HII.ev', &
                  '      grid_HeI.ev', &
                  '     grid_HeII.ev', &
                  '    grid_HeIII.ev' /)
    headerline = (/ '             # HI', &
                    '            # HII', &
                    '            # HeI', &
                    '           # HeII', &
                    '          # HeIII' /)
 end select

 allocate(iorder(npart))
 if (use_mass_coord) then
    call set_r2func_origin(xyzmh_ptmass(1,1),xyzmh_ptmass(2,1),xyzmh_ptmass(3,1)) ! Order particles by distance from core
    call indexxfunc(npart,r2func_origin,xyzh,iorder)
 else
    iorder = (/(i, i=1,npart, 1)/) ! Have iorder(k) be same as k
 endif

 do k=1,npart
    i = iorder(k) ! Loop from innermost to outermost particle
    if (use_mass_coord) then
       coord(i) = real(k-1) ! Number of particles interior to particle k
    else
       coord(i) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
    endif

    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    select case (iquantity)
    case(1) ! Energy
       call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
       call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)
       quant(i,1) = ekini + epoti + ethi
    case(2) ! Entropy
       if ((ieos==10) .and. (ientropy==2)) then
          call getvalue_mesa(rhopart*unit_density,vxyzu(4,i)*unit_ergg,3,pgas,ierr) ! Get gas pressure
          mu = rhopart*unit_density * Rg * eos_vars(itemp,i) / pgas
       else
          mu = gmw
       endif
       if ((ieos==10) .and. (ientropy==3)) then
          quant(i,1) = entropy(rhopart*unit_density,ponrhoi*rhopart*unit_pressure,mu,ientropy,vxyzu(4,i)*unit_ergg,ierr)
       else
          quant(i,1) = entropy(rhopart*unit_density,ponrhoi*rhopart*unit_pressure,mu,ientropy,ierr=ierr)
       endif
    case(3) ! Bernoulli energy (per unit mass)
       call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
       quant(i,1) = 0.5*dot_product(vxyzu(1:3,i),vxyzu(1:3,i)) + ponrhoi + vxyzu(4,i) + epoti/particlemass ! 1/2 v^2 + P/rho + phi
    case(4) ! Ion fraction
       call ionisation_fraction(rhopart*unit_density,eos_vars(itemp,i),X_in,1.-X_in-Z_in,xh0,xh1,xhe0,xhe1,xhe2)
       quant(i,1) = xh0
       quant(i,2) = xh1
       quant(i,3) = xhe0
       quant(i,4) = xhe1
       quant(i,5) = xhe2
    end select
 enddo

 if (use_mass_coord) coord = coord * particlemass + xyzmh_ptmass(4,1)

 write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
 do i=1,nvars
    call histogram_setup(coord,quant(:,i),hist,npart,maxcoord,mincoord,nbins,.true.,ilogbins)
    if (dump_number == 0) then
       unitnum = 1000
       open(unit=unitnum, file=trim(adjustl(filename(i))), status='replace')
       write(unitnum, "(a)") trim(headerline(i))
       close(unit=unitnum)
    endif
    unitnum=1001+i
    open(unit=unitnum, file=trim(adjustl(filename(i))), status='old', position='append')
    write(unitnum,data_formatter) time,hist
    close(unit=unitnum)
 enddo
 deallocate(iorder,coord,headerline,filename,quant,hist)

end subroutine energy_profile


!----------------------------------------------------------------
!+
!  Rotation profiles
!+
!----------------------------------------------------------------
subroutine rotation_profile(time,num,npart,xyzh,vxyzu)
 use vectorutils, only:cross_product3D
 integer, intent(in)          :: npart,num
 real, intent(in)             :: time
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 integer                      :: nbins
 real, allocatable            :: rad_part(:)
 real, allocatable            :: hist_var(:),dist_part(:,:)
 real                         :: minloga,maxloga,sep_vector(3),vel_vector(3),J_vector(3),xyz_origin(3),vxyz_origin(3),omega,vphi
 character(len=17), allocatable :: grid_file(:)
 character(len=40)            :: data_formatter
 integer                      :: i,unitnum,nfiles,iradius

 nbins = 500
 minloga = 0.5
 maxloga = 4.3
 iradius = 1 ! 1: Bin by cylindrical radius; 2: Bin by spherical radius; 3: Bin by cylindrical radius from CM

 nfiles = 2
 allocate(hist_var(nbins),grid_file(nfiles),dist_part(nfiles,npart),rad_part(npart))
 rad_part = 0.
 dist_part = 0.
 grid_file = (/ '    grid_omega.ev', &
                '       grid_Jz.ev' /)

 select case (iradius)
 case(1,2) ! Take donor core as origin
    xyz_origin = xyzmh_ptmass(1:3,1)
    vxyz_origin = vxyz_ptmass(1:3,1)
 case(3) ! Take sink CM as origin
    xyz_origin = (xyzmh_ptmass(1:3,1)*xyzmh_ptmass(4,1) + xyzmh_ptmass(1:3,2)*xyzmh_ptmass(4,2)) / (xyzmh_ptmass(4,1) + &
                  xyzmh_ptmass(4,2))
    vxyz_origin = (vxyz_ptmass(1:3,1)*xyzmh_ptmass(4,1) + vxyz_ptmass(1:3,2)*xyzmh_ptmass(4,2)) / (xyzmh_ptmass(4,1) + &
                  xyzmh_ptmass(4,2))
 end select

 do i=1,npart
    select case (iradius)
    case(1,3) ! Bin by cylindrical radius
       rad_part(i) = sqrt( dot_product(xyzh(1:2,i) - xyz_origin(1:2), xyzh(1:2,i) - xyz_origin(1:2)) )
    case(2) ! Bin by spherical radius
       rad_part(i) = separation(xyzh(1:3,i),xyz_origin)
    end select

    call get_gas_omega(xyz_origin,vxyz_origin,xyzh(1:3,i),vxyzu(1:3,i),vphi,omega)
    dist_part(1,i) = vphi

    sep_vector = xyzh(1:3,i) - xyz_origin(1:3)
    vel_vector = vxyzu(1:3,i) - vxyz_origin(1:3)
    call cross_product3D(vel_vector, sep_vector, J_vector)
    dist_part(2,i) = dot_product(J_vector, (/0.,0.,1./))
 enddo

 do i=1,nfiles
    call histogram_setup(rad_part(1:npart),dist_part(i,1:npart),hist_var,npart,maxloga,minloga,nbins,.true.,.true.)
    write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
    if (num == 0) then
       unitnum = 1000
       open(unit=unitnum, file=trim(adjustl(grid_file(i))), status='replace')
       write(unitnum, "(a)") '# z-component of angular velocity'
       close(unit=unitnum)
    endif
    unitnum=1001+i
    open(unit=unitnum, file=trim(adjustl(grid_file(i))), position='append')
    write(unitnum,data_formatter) time,hist_var(:)
    close(unit=unitnum)
 enddo
 deallocate(hist_var,grid_file,dist_part,rad_part)

end subroutine rotation_profile


!----------------------------------------------------------------
!+
!  Velocity distribution
!+
!----------------------------------------------------------------
subroutine velocity_histogram(time,num,npart,particlemass,xyzh,vxyzu)
 use part,           only:eos_vars,itemp
 use ionization_mod, only:calc_thermal_energy
 real, intent(in)    :: time,particlemass
 integer, intent(in) :: npart,num
 real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=40)   :: data_formatter
 character(len=40)   :: file_name1,file_name2
 integer             :: i,iu1,iu2,ncols
 real                :: ponrhoi,rhopart,spsoundi,phii,epoti,ekini,einti,tempi,ethi,dum
 real, allocatable   :: vbound(:),vunbound(:),vr(:)

 allocate(vbound(npart),vunbound(npart),vr(npart))
 do i = 1,npart
    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
    call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)
    vr(i) = dot_product(xyzh(1:3,i),vxyzu(1:3,i)) / sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))

    if (ekini+epoti > 0.) then
       vbound(i) = -1.e15
       vunbound(i) = vr(i)
    else
       vbound(i) = vr(i)
       vunbound(i) = -1.e15
    endif
 enddo

 ncols = npart
 write(data_formatter, "(a,I6.6,a)") "(", ncols+1, "(2x,es18.11e2))"
 file_name1 = "vel_bound.ev"
 file_name2 = "vel_unbound.ev"

 if (dump_number == 0) then
    open(newunit=iu1, file=file_name1, status='replace')
    open(newunit=iu2, file=file_name2, status='replace')
 else
    open(newunit=iu1, file=file_name1, position='append')
    open(newunit=iu2, file=file_name2, position='append')
 endif

 write(iu1,data_formatter) time,vbound
 write(iu2,data_formatter) time,vunbound
 close(unit=iu1)
 close(unit=iu2)

 deallocate(vbound,vunbound,vr)

end subroutine velocity_histogram


!----------------------------------------------------------------
!+
!  Velocity profile
!+
!----------------------------------------------------------------
subroutine velocity_profile(time,num,npart,particlemass,xyzh,vxyzu)
 real, intent(in)    :: time,particlemass
 integer, intent(in) :: npart,num
 real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=40)   :: data_formatter
 character(len=40)   :: file_name
 integer             :: i,nbins,iu,count
 real                :: rmin,rmax,xyz_origin(3),vxyz_origin(3),vphi,omega,&
                        theta1,theta2,tantheta1,tantheta2,tantheta
 real, allocatable, dimension(:) :: rad_part,dist_part,hist

 nbins = 500
 rmin = 0.
 rmax = 5.

 allocate(hist(nbins),dist_part(npart),rad_part(npart))
 dist_part = 0.
 file_name = '  vphi_profile.ev'

 ! Select origin
 xyz_origin = xyzmh_ptmass(1:3,1)
 vxyz_origin = vxyz_ptmass(1:3,1)

 ! Masking in polar angle
 theta1 = 75.   ! Polar angle in deg
 theta2 = 105.
 tantheta1 = tan(theta1*3.14159/180.)
 tantheta2 = tan(theta2*3.14159/180.)

 count = 0
 do i = 1,npart
    rad_part(i) = sqrt( dot_product(xyzh(1:2,i) - xyz_origin(1:2), xyzh(1:2,i) - xyz_origin(1:2)) )  ! Cylindrical radius

    ! Masking in polar angle
    tantheta = rad_part(i)/(xyzh(3,i) - xyzmh_ptmass(3,1))
    if ( (tantheta>0. .and. tantheta<tantheta1) .or. (tantheta<0. .and. tantheta>tantheta2) ) cycle

    call get_gas_omega(xyz_origin,vxyz_origin,xyzh(1:3,i),vxyzu(1:3,i),vphi,omega)
    dist_part(i) = vphi
    count = count + 1
 enddo

 call histogram_setup(rad_part,dist_part,hist,count,rmax,rmin,nbins,.true.,.false.)
 write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
 if (num == 0) then
    open(newunit=iu, file=trim(adjustl(file_name)), status='replace')
    write(iu, "(a)") '# Azimuthal velocity profile'
    close(unit=iu)
 endif
 open(newunit=iu, file=trim(adjustl(file_name)), position='append')
 write(iu,data_formatter) time,hist
 close(unit=iu)
 deallocate(hist,dist_part,rad_part)

end subroutine velocity_profile


!----------------------------------------------------------------
!+
!  Specific z-angular momentum profile
!+
!----------------------------------------------------------------
subroutine angular_momentum_profile(time,num,npart,particlemass,xyzh,vxyzu)
 real, intent(in)    :: time,particlemass
 integer, intent(in) :: npart,num
 real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=40)   :: data_formatter
 character(len=40)   :: file_name
 integer             :: i,nbins,iu,count
 real                :: rmin,rmax,xyz_origin(3),vxyz_origin(3),&
                        theta1,theta2,tantheta1,tantheta2,tantheta
 real, allocatable, dimension(:) :: rad_part,dist_part,hist

 nbins = 500
 rmin = 0.
 rmax = 5.

 allocate(hist(nbins),dist_part(npart),rad_part(npart))
 dist_part = 0.
 file_name = '    jz_profile.ev'

 ! Select origin
 xyz_origin = xyzmh_ptmass(1:3,1)
 vxyz_origin = vxyz_ptmass(1:3,1)

 ! Masking in polar angle
 theta1 = 75.   ! Polar angle in deg
 theta2 = 105.
 tantheta1 = tan(theta1*3.14159/180.)
 tantheta2 = tan(theta2*3.14159/180.)

 count = 0
 do i = 1,npart
    rad_part(i) = sqrt( dot_product(xyzh(1:2,i) - xyz_origin(1:2), xyzh(1:2,i) - xyz_origin(1:2)) )  ! Cylindrical radius

    ! Masking in polar angle
    tantheta = rad_part(i)/(xyzh(3,i) - xyzmh_ptmass(3,1))
    if ( (tantheta>0. .and. tantheta<tantheta1) .or. (tantheta<0. .and. tantheta>tantheta2) ) cycle

    dist_part(i) = ( (xyzh(1,i)-xyz_origin(1))*(vxyzu(2,i)-vxyz_origin(2)) - &
                   (xyzh(2,i)-xyz_origin(2))*(vxyzu(1,i)-vxyz_origin(1)) )
    count = count + 1
 enddo

 call histogram_setup(rad_part,dist_part,hist,count,rmax,rmin,nbins,.true.,.false.)
 write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
 if (num == 0) then
    open(newunit=iu, file=trim(adjustl(file_name)), status='replace')
    write(iu, "(a)") '# z-angular momentum profile'
    close(unit=iu)
 endif
 open(newunit=iu, file=trim(adjustl(file_name)), position='append')
 write(iu,data_formatter) time,hist
 close(unit=iu)

end subroutine angular_momentum_profile


!----------------------------------------------------------------
!+
!  Keplerian velocity profile
!+
!----------------------------------------------------------------
subroutine vkep_profile(time,num,npart,particlemass,xyzh,vxyzu)
 use sortutils, only:set_r2func_origin,r2func_origin,find_rank
 use part,      only:iorder=>ll
 real, intent(in)    :: time,particlemass
 integer, intent(in) :: npart,num
 real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=40)   :: data_formatter,file_name
 integer             :: i,nbins,iu
 real                :: rmin,rmax,massi,Mtot
 real, allocatable   :: hist(:),rad_part(:),dist_part(:)

 nbins = 500
 rmin = 0.
 rmax = 5.

 allocate(hist(nbins),rad_part(npart),dist_part(npart))
 dist_part = 0.
 file_name = '  vkep_profile.ev'

 call set_r2func_origin(xyzmh_ptmass(1,1),xyzmh_ptmass(2,1),xyzmh_ptmass(3,1))
 call find_rank(npart,r2func_origin,xyzh(1:3,:),iorder)

 Mtot = npart*particlemass
 do i = 1,npart
    massi = Mtot * real(iorder(i)-1) / real(npart) + xyzmh_ptmass(4,1)
    rad_part(i) = separation( xyzh(1:3,i), xyzmh_ptmass(1:3,1) )
    dist_part(i) = sqrt(massi/rad_part(i))
 enddo

 call histogram_setup(rad_part,dist_part,hist,npart,rmax,rmin,nbins,.true.,.false.)
 write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))"
 if (num == 0) then
    open(newunit=iu, file=trim(adjustl(file_name)), status='replace')
    write(iu, "(a)") '# Keplerian velocity profile'
    close(unit=iu)
 endif
 open(newunit=iu, file=trim(adjustl(file_name)), position='append')
 write(iu,data_formatter) time,hist
 close(unit=iu)
 deallocate(hist,dist_part,rad_part)

end subroutine vkep_profile


!----------------------------------------------------------------
!+
!  Planet profile
!+
!----------------------------------------------------------------
subroutine planet_profile(num,dumpfile,particlemass,xyzh,vxyzu)
 character(len=*), intent(in) :: dumpfile
 integer, intent(in)        :: num
 real, intent(in)           :: particlemass
 real, intent(inout)        :: xyzh(:,:),vxyzu(:,:)
 character(len=40)          :: file_name
 integer                    :: i,maxrho_ID,iu
 integer, save              :: nplanet
 integer, allocatable, save :: planetIDs(:)
 real                       :: rhoprev
 real, dimension(3)         :: planet_com,planet_vcom,vnorm,ri,Rvec
 real, allocatable          :: R(:),z(:),rho(:)

 if (dump_number ==0 ) call get_planetIDs(nplanet,planetIDs)
 allocate(R(nplanet),z(nplanet),rho(nplanet))

 ! Find highest density in planet
 rhoprev = 0.
 maxrho_ID = planetIDs(1)
 do i = 1,nplanet
    rho(i) = rhoh(xyzh(4,planetIDs(i)), particlemass)
    if (rho(i) > rhoprev) then
       maxrho_ID = planetIDs(i)
       rhoprev = rho(i)
    endif
 enddo
 planet_com = xyzh(1:3,maxrho_ID)
 planet_vcom = vxyzu(1:3,maxrho_ID)
 vnorm = planet_vcom / sqrt(dot_product(planet_vcom,planet_vcom))

 ! Write to file
 file_name = trim(dumpfile)//".planetpart"
 open(newunit=iu, file=file_name, status='replace')

 ! Record R and z cylindrical coordinates w.r.t. planet_com
 do i = 1,nplanet
    ri = xyzh(1:3,planetIDs(i)) - planet_com ! Particle position w.r.t. planet_com
    z(i) = dot_product(ri, vnorm)
    Rvec = ri - z(i)*vnorm
    R(i) = sqrt(dot_product(Rvec,Rvec))
   !  write(iu,"(es13.6,2x,es13.6,2x,es13.6)") R(i),z(i),rho(i)
    write(iu,"(es13.6,2x,es13.6,2x,es13.6,2x,es13.6,2x,es13.6)") xyzh(1,i),xyzh(2,i),xyzh(3,i),rho(i),vxyzu(4,i)
 enddo

 close(unit=iu)
 deallocate(R,z,rho)

end subroutine planet_profile


!----------------------------------------------------------------
!+
!  Unbound profiles
!+
!----------------------------------------------------------------
subroutine unbound_profiles(time,num,npart,particlemass,xyzh,vxyzu)
 use ionization_mod, only:calc_thermal_energy
 integer, intent(in)                          :: npart,num
 real,    intent(in)                          :: time,particlemass
 real,    intent(inout)                       :: xyzh(:,:),vxyzu(:,:)
 integer, dimension(4)                        :: npart_hist
 real,    dimension(5,npart)                  :: dist_part,rad_part
 real,    dimension(:), allocatable           :: hist_var
 real                                         :: etoti,ekini,einti,epoti,ethi,phii,dum,rhopart,ponrhoi,spsoundi,tempi
 real                                         :: maxloga,minloga
 character(len=18), dimension(4)              :: grid_file
 character(len=40)                            :: data_formatter
 logical, allocatable, save                   :: prev_unbound(:,:),prev_bound(:,:)
 integer                                      :: i,unitnum,nbins

 call compute_energies(time)
 npart_hist = 0     ! Stores number of particles fulfilling each of the four bound/unbound criterion
 nbins      = 500
 rad_part   = 0.    ! (4,npart_hist)-array storing separations of particles
 dist_part  = 0.
 minloga    = 0.5
 maxloga    = 4.3

 allocate(hist_var(nbins))
 grid_file = (/ 'grid_unbound_th.ev', &
                'grid_unbound_kp.ev', &
                ' grid_bound_kpt.ev', &
                '  grid_bound_kp.ev' /)

 if (dump_number == 0) then
    allocate(prev_bound(2,npart))
    allocate(prev_unbound(2,npart))
    prev_bound   = .false.
    prev_unbound = .false.
 endif


 do i=1,npart
    if (.not. isdead_or_accreted(xyzh(4,i))) then
       rhopart = rhoh(xyzh(4,i), particlemass)
       call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
       call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
       call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,tempi,gamma,ethi)
       etoti = ekini + epoti + ethi

       ! Ekin + Epot + Eth > 0
       if ((etoti > 0.) .and. (.not. prev_unbound(1,i))) then
          npart_hist(1) = npart_hist(1) + 1 ! Keeps track of number of particles that have become newly unbound in this dump
          rad_part(1,npart_hist(1)) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          dist_part(1,npart_hist(1)) = 1. ! Array of ones with size of npart_hist(1)?
          prev_unbound(1,i) = .true.
       elseif (etoti < 0.) then
          prev_unbound(1,i) = .false.
       endif

       ! Ekin + Epot > 0
       if ((ekini + epoti > 0.) .and. (.not. prev_unbound(2,i))) then
          npart_hist(2) = npart_hist(2) + 1
          rad_part(2,npart_hist(2)) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          dist_part(2,npart_hist(2)) = 1.
          prev_unbound(2,i) = .true.
       elseif (ekini + epoti < 0.) then
          prev_unbound(2,i) = .false.
       endif

       ! Ekin + Epot + Eth < 0
       if ((etoti < 0.) .and. (.not. prev_bound(1,i))) then
          npart_hist(3) = npart_hist(3) + 1
          rad_part(3,npart_hist(3)) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          dist_part(3,npart_hist(3)) = 1.
          prev_bound(1,i) = .true.
       elseif (etoti > 0.) then
          prev_bound(1,i) = .false.
       endif

       ! Ekin + Epot < 0
       if ((ekini + epoti < 0.) .and. (.not. prev_bound(2,i))) then
          npart_hist(4) = npart_hist(4) + 1
          rad_part(4,npart_hist(4)) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          dist_part(4,npart_hist(4)) = 1.
          prev_bound(2,i) = .true.
       elseif (ekini + epoti > 0.) then
          prev_bound(2,i) = .false.
       endif
    endif
 enddo

 do i=1,4
    call histogram_setup(rad_part(i,1:npart_hist(i)),dist_part(i,1:npart_hist(i)),hist_var,npart_hist(i),maxloga,minloga,nbins,&
                        .false.,.true.)

    write(data_formatter, "(a,I5,a)") "(", nbins+1, "(3x,es18.10e3,1x))" ! Time column plus nbins columns

    if (num == 0) then ! Write header line
       unitnum = 1000
       open(unit=unitnum, file=trim(adjustl(grid_file(i))), status='replace')
       write(unitnum, "(a)") '# Newly bound/unbound particles'
       close(unit=unitnum)
    endif

    unitnum=1001+i

    open(unit=unitnum, file=trim(adjustl(grid_file(i))), position='append')

    write(unitnum,"()")
    write(unitnum,data_formatter) time,hist_var(:)

    close(unit=unitnum)
 enddo
 deallocate(hist_var)

end subroutine unbound_profiles


!----------------------------------------------------------------
!+
!  Unbound ion fractions: Look at distribution of ion fraction when given particle is unbound
!+
!----------------------------------------------------------------
subroutine unbound_ionfrac(time,npart,particlemass,xyzh,vxyzu)
 use ionization_mod, only:calc_thermal_energy,get_xion,ionisation_fraction
 integer, intent(in)       :: npart
 real,    intent(in)       :: time,particlemass
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 character(len=17)         :: columns(5)
 integer                   :: i
 real                      :: etoti,ekini,einti,epoti,ethi,phii,dum,rhopart,xion(1:4),&
                              ponrhoi,spsoundi,tempi,xh0,xh1,xhe0,xhe1,xhe2
 logical, allocatable, save :: prev_unbound(:),prev_bound(:)
 real, allocatable, save :: ionfrac(:,:)

 columns = (/'        xion1', &
             '        xion2', &
             '        xion3', &
             '        xion4', &
             '      1-xion3' /)

 if (dump_number == 0) then
    allocate(prev_unbound(npart),prev_bound(npart))
    prev_bound   = .false.
    prev_unbound = .false.
    allocate(ionfrac(npart,5))
    ionfrac = -1. ! Initialise ion states to -1
 endif

 call compute_energies(time)
 do i=1,npart
    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
    call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,tempi,gamma,ethi)
    etoti = ekini + epoti + ethi

    if ((etoti > 0.) .and. (.not. prev_unbound(i))) then
       if (ieos == 10) then  ! MESA EoS
          call ionisation_fraction(rhopart*unit_density,tempi,X_in,1.-X_in-Z_in,xh0,xh1,xhe0,xhe1,xhe2)
       elseif (ieos == 20) then  ! Gas + radiation + recombination EoS
          call get_xion(log10(rhopart*unit_density),tempi,X_in,1.-X_in-Z_in,xion)
          xh0 = xion(1)  ! H2 ionisation fraction
          xh1 = xion(2)  ! H ionisation fraction
          xhe1 = xion(3) ! He ionisation to He+ fraction
          xhe2 = xion(4) ! He+ ionisation to He++ fraction
          xhe0 = 1.-xion(3)
       else  ! Not supported
          print*,"Error, not sensible to use unbound_ionfrac when not using MESA EoS (ieos=10) or gas+rad+rec EoS (ieos=20)"
          stop
       endif
       ionfrac(i,1) = xh0
       ionfrac(i,2) = xh1
       ionfrac(i,3) = xhe1
       ionfrac(i,4) = xhe2
       ionfrac(i,5) = xhe0
       prev_unbound(i) = .true.
    elseif (etoti < 0.) then
       prev_unbound(i) = .false.
    endif
 enddo

 ! Trick write_time_file into writing my data table
 print*,'Dump number is ',dump_number
 if (dump_number == 258) then
    do i=1,npart
       call write_time_file("unbound_ionfrac",columns,-1.,ionfrac(i,1:5),5,i-1) ! Set num = i-1 so that header will be written for particle 1 and particle 1 only
    enddo
 endif

end subroutine unbound_ionfrac


!----------------------------------------------------------------
!+
!  Unbound temperature
!+
!----------------------------------------------------------------
subroutine unbound_temp(time,npart,particlemass,xyzh,vxyzu)
 use part,           only:eos_vars,itemp
 use ionization_mod, only:calc_thermal_energy,get_xion
 integer, intent(in)        :: npart
 real,    intent(in)        :: time,particlemass
 real,    intent(inout)     :: xyzh(:,:),vxyzu(:,:)
 character(len=17)          :: columns(1)
 integer                    :: i,final_count(7)
 real                       :: etoti,ekini,einti,epoti,ethi,phii,dum,rhopart,&
                               ponrhoi,spsoundi,temp_bins(7)
 logical, allocatable, save :: prev_unbound(:),prev_bound(:)
 real, allocatable, save    :: temp_unbound(:)

 columns = (/'         temp'/)

 if (dump_number == 0) then
    allocate(prev_unbound(npart),prev_bound(npart),temp_unbound(npart))
    prev_bound   = .false.
    prev_unbound = .false.
    temp_unbound = 0. ! Initialise temperatures to 0.
 endif

 do i=1,npart
    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),eos_vars(itemp,i),vxyzu(4,i))
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
    call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)
    etoti = ekini + epoti + ethi

    if ((etoti > 0.) .and. (.not. prev_unbound(i))) then
       temp_unbound(i) = eos_vars(itemp,i)
       prev_unbound(i) = .true.
    elseif (etoti < 0.) then
       prev_unbound(i) = .false.
    endif
 enddo

 print*,'dump_number=',dump_number
 ! Trick write_time_file into writing my data table
 if (dump_number == 167) then
    temp_bins = (/ 2.e3, 5.5e3, 8.e3, 1.5e4, 2.e4, 4.e4, 1.e15 /)
    final_count = 0
    do i=1,npart
       if (temp_unbound(i) > 1.e-15) then
          if (temp_unbound(i) < temp_bins(1)) then
             final_count(1:7) = final_count(1:7) + 1
          elseif (temp_unbound(i) < temp_bins(2)) then
             final_count(2:7) = final_count(2:7) + 1
          elseif (temp_unbound(i) < temp_bins(3)) then
             final_count(3:7) = final_count(3:7) + 1
          elseif (temp_unbound(i) < temp_bins(4)) then
             final_count(4:7) = final_count(4:7) + 1
          elseif (temp_unbound(i) < temp_bins(5)) then
             final_count(5:7) = final_count(5:7) + 1
          elseif (temp_unbound(i) < temp_bins(6)) then
             final_count(6:7) = final_count(6:7) + 1
          elseif (temp_unbound(i) < temp_bins(7)) then
             final_count(7) = final_count(7) + 1
          endif
       endif
       call write_time_file("unbound_temp",columns,-1.,temp_unbound(i),1,i-1) ! Set num = i-1 so that header will be written for particle 1 and particle 1 only
    enddo

    print*,final_count
 endif

end subroutine unbound_temp


!----------------------------------------------------------------
!+
!  Recombination statistics
!+
!----------------------------------------------------------------
subroutine recombination_stats(time,num,npart,particlemass,xyzh,vxyzu)
 use part, only:eos_vars,itemp
 use ionization_mod, only:calc_thermal_energy,ionisation_fraction
 integer, intent(in)       :: npart,num
 real,    intent(in)       :: time,particlemass
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 real                      :: etoti,ekini,einti,epoti,ethi,phii,dum,rhopart,&
                              ponrhoi,spsoundi,tempi,pressure,temperature,xh0,xh1,xhe0,xhe1,xhe2
 character(len=40)         :: data_formatter,logical_format
 logical, allocatable      :: isbound(:)
 integer, allocatable      :: H_state(:),He_state(:)
 integer                   :: i
 real, parameter           :: recomb_th=0.05

 call compute_energies(time)

 allocate(isbound(npart),H_state(npart),He_state(npart))
 do i=1,npart
    ! Calculate total energy
    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,dum)
    call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)
    etoti = ekini + epoti + ethi

    call get_eos_pressure_temp_mesa(rhopart*unit_density,vxyzu(4,i)*unit_ergg,pressure,temperature) ! This should depend on ieos
    call ionisation_fraction(rhopart*unit_density,temperature,X_in,1.-X_in-Z_in,xh0,xh1,xhe0,xhe1,xhe2)

    ! Is unbound?
    if (etoti > 0.) then
       isbound(i) = .false.
    else
       isbound(i) = .true.
    endif

    ! H ionisation state
    if (xh0 > recomb_th) then
       H_state(i) = 1
    elseif (xh1 > recomb_th) then
       H_state(i) = 2
    else
       H_state(i) = 0 ! This should not happen
    endif

    ! H ionisation state
    if (xhe0 > recomb_th) then
       He_state(i) = 1
    elseif (xhe1 > recomb_th) then
       He_state(i) = 2
    elseif (xhe2 > recomb_th) then
       He_state(i) = 3
    else
       He_state(i) = 0 ! This should not happen
    endif
 enddo

 write(data_formatter, "(a,I5,a)") "(es18.10e3,", npart, "(1x,i1))" ! Time column plus npart columns
 write(logical_format, "(a,I5,a)") "(es18.10e3,", npart, "(1x,L))" ! Time column plus npart columns

 if (num == 0) then ! Write header line
    open(unit=1000, file="H_state.ev", status='replace')
    write(1000, "(a)") '# Ion fraction statistics'
    close(unit=1000)
    open(unit=1001, file="He_state.ev", status='replace')
    write(1001, "(a)") '# Ion fraction statistics'
    close(unit=1001)
    open(unit=1002, file="isbound.ev", status='replace')
    write(1002, "(a)") '# Ion fraction statistics'
    close(unit=1002)
 endif

 open(unit=1000, file="H_state.ev", position='append')
 write(1000,data_formatter) time,H_state(:)
 close(unit=1000)

 open(unit=1000, file="He_state.ev", position='append')
 write(1000,data_formatter) time,He_state(:)
 close(unit=1000)

 open(unit=1000, file="isbound.ev", position='append')
 write(1000,logical_format) time,isbound(:)
 close(unit=1000)

 deallocate(isbound,H_state,He_state)

end subroutine recombination_stats


!----------------------------------------------------------------
!+
!  Sink properties
!+
!----------------------------------------------------------------
subroutine sink_properties(time,npart,particlemass,xyzh,vxyzu)
 use vectorutils, only:cross_product3D
 integer, intent(in)          :: npart
 real, intent(in)             :: time, particlemass
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 character(len=17)            :: filename
 real                         :: sinkcomp(35)
 real                         :: ang_mom(3)
 real                         :: phitot, dtsinksink, fonrmax
 real                         :: fxi, fyi, fzi, phii
 real, dimension(4,maxptmass) :: fssxyz_ptmass
 real, dimension(4,maxptmass) :: fxyz_ptmass
 real, dimension(3)           :: com_xyz,com_vxyz
 integer                      :: i,ncols,merge_n,merge_ij(nptmass)

 ncols = 31
 allocate(columns(ncols))
 columns = (/'            x', &
             '            y', &
             '            z', &
             '            r', &
             '           vx', &
             '           vy', &
             '           vz', &
             '          |v|', &
             '           px', &
             '           py', &
             '           pz', &
             '          |p|', &
             '         fssx', &
             '         fssy', &
             '         fssz', &
             '        |fss|', &
             '          fsx', &
             '          fsy', &
             '          fsz', &
             '         |fs|', &
             '    ang mom x', &
             '    ang mom y', &
             '    ang mom z', &
             '    |ang mom|', &
             '       kin en', &
             '       CoM x ', &
             '       CoM y ', &
             '       CoM z ', &
             '       CoM vx', &
             '       CoM vy', &
             '       CoM vz' /)

 fxyz_ptmass = 0.
 call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,phitot,dtsinksink,0,0.,merge_ij,merge_n)
 fssxyz_ptmass = fxyz_ptmass
 do i=1,npart
    call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                            fxi,fyi,fzi,phii,particlemass,fxyz_ptmass,fonrmax)
 enddo

 ! Determine position and velocity of the CoM
 call orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)

 do i=1,nptmass
    sinkcomp = 0.
    write (filename, "(A16,I0)") "sink_properties_", i

    ! position xyz
    sinkcomp(1:3)   = xyzmh_ptmass(1:3,i)
    ! position modulus
    sinkcomp(4)     = distance(xyzmh_ptmass(1:3,i))
    ! velocity xyz
    sinkcomp(5:7)   = vxyz_ptmass(1:3,i)
    ! velocity modulus
    sinkcomp(8)     = distance(vxyz_ptmass(1:3,i))
    ! momentum xyz
    sinkcomp(9:11)  = xyzmh_ptmass(4,i)*vxyz_ptmass(1:3,i)
    ! momentum modulus
    sinkcomp(12)    = xyzmh_ptmass(4,i)*sinkcomp(8)
    ! force xyz
    sinkcomp(13:15) = fssxyz_ptmass(1:3,i)
    ! force modulus
    sinkcomp(16)    = distance(fssxyz_ptmass(1:3,i))
    ! tot force xyz
    sinkcomp(17:19) = fxyz_ptmass(1:3,i)
    ! tot force modulus
    sinkcomp(20)    = distance(fxyz_ptmass(1:3,i))
    ! angular momentum xyz
    call cross_product3D(xyzmh_ptmass(1:3,i), xyzmh_ptmass(4,i)*vxyz_ptmass(1:3,i), ang_mom)
    sinkcomp(21:23) = ang_mom
    ! angular momentum modulus
    sinkcomp(24)    = distance(ang_mom(1:3))
    ! kinetic energy
    sinkcomp(25)    = 0.5*xyzmh_ptmass(4,i)*sinkcomp(8)**2
    ! CoM position
    sinkcomp(26:28) = com_xyz(1:3)
    ! CoM velocity
    sinkcomp(29:31) = com_vxyz(1:3)

    call write_time_file(filename, columns, time, sinkcomp, ncols, dump_number)
 enddo
 deallocate(columns)

end subroutine sink_properties



subroutine env_binding_ene(npart,particlemass,xyzh,vxyzu)
 use part, only:eos_vars,itemp
 use ionization_mod, only:calc_thermal_energy
 integer, intent(in)    :: npart
 real, intent(in)       :: particlemass
 real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                :: i
 real                   :: ethi,phii,rhoi,ponrhoi,spsoundi,tempi,dum1,dum2,dum3
 real                   :: bind_g,bind_th,bind_int,eth_tot,eint_tot

 bind_g = 0.
 bind_th = 0.
 bind_int = 0.
 eint_tot = 0.
 eth_tot = 0.
 do i=1,npart
    ! Gas-gas potential
    bind_g = bind_g + poten(i) ! Double counting factor of 1/2 already included in poten

    ! Sink-sink potential
    phii = 0.
    call get_accel_sink_gas(1,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass(:,1),dum1,dum2,dum3,phii) ! Include only core particle; no companion
    bind_g = bind_g + particlemass * phii

    rhoi = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhoi,eos_vars(itemp,i),gamma,ethi)

    eth_tot = eth_tot + ethi
    eint_tot = eint_tot + particlemass * vxyzu(4,i)
 enddo
 bind_th  = bind_g + eth_tot
 bind_int = bind_g + eint_tot

 print*,bind_g*unit_energ, bind_th*unit_energ, bind_int*unit_energ

end subroutine env_binding_ene


subroutine bound_unbound_thermo(time,npart,particlemass,xyzh,vxyzu)
 integer, intent(in)          :: npart
 real, intent(in)             :: time, particlemass
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 integer                      :: i, ncols
 real, dimension(8)           :: entropy_array
 real                         :: etoti, ekini, einti, epoti, phii, rhopart
 real                         :: pres_1, proint_1, peint_1, temp_1
 real                         :: troint_1, teint_1, entrop_1, abad_1, gamma1_1, gam_1
 integer, parameter           :: ient_b   = 1
 integer, parameter           :: ient_ub  = 2
 integer, parameter           :: itemp_b  = 3
 integer, parameter           :: itemp_ub = 4
 integer, parameter           :: ipres_b  = 5
 integer, parameter           :: ipres_ub = 6
 integer, parameter           :: idens_b  = 7
 integer, parameter           :: idens_ub = 8

 !zeroes the entropy variable and others
 entropy_array = 0.

 !setup
 if (dump_number == 0) then
    call prompt('Would you like to use thermal energy in the computation of the bound/unbound status?', switch(1),.false.)
 endif

 call compute_energies(time)

 do i=1,npart
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

    rhopart = rhoh(xyzh(4,i), particlemass)

    !gets entropy for the current particle
    call get_eos_various_mesa(rhopart*unit_density,vxyzu(4,i) * unit_ergg, &
                              pres_1,proint_1,peint_1,temp_1,troint_1, &
                              teint_1,entrop_1,abad_1,gamma1_1,gam_1)

    !sums entropy and other quantities for bound particles and unbound particles

    if (.not. switch(1)) then
       etoti = etoti - einti
    endif

    if (etoti < 0.0) then !bound
       entropy_array(ient_b)  = entropy_array(ient_b) + entrop_1
       entropy_array(itemp_b) = entropy_array(itemp_b) + temp_1
       entropy_array(ipres_b) = entropy_array(ipres_b) + pres_1
       entropy_array(idens_b) = entropy_array(idens_b) + rhopart*unit_density

    else !unbound
       entropy_array(ient_ub)  = entropy_array(ient_ub) + entrop_1
       entropy_array(itemp_ub) = entropy_array(itemp_ub) + temp_1
       entropy_array(ipres_ub) = entropy_array(ipres_ub) + pres_1
       entropy_array(idens_ub) = entropy_array(idens_ub) + rhopart*unit_density

    endif

 enddo

 !average
 entropy_array(itemp_b)  = entropy_array(itemp_b) / npart
 entropy_array(itemp_ub) = entropy_array(itemp_ub) / npart
 entropy_array(ipres_b)  = entropy_array(ipres_b) / npart
 entropy_array(ipres_ub) = entropy_array(ipres_ub) / npart
 entropy_array(idens_b)  = entropy_array(idens_b) / npart
 entropy_array(idens_ub) = entropy_array(idens_ub) / npart

 !writes on file
 ncols = 8
 allocate(columns(ncols))
 columns = (/'       b entr',&
                '     unb entr',&
                '   avg b temp',&
                ' avg unb temp',&
                '   avg b pres',&
                ' avg unb pres',&
                '   avg b dens',&
                ' avg unb dens'/)
 call write_time_file('entropy_vs_time', columns, time, entropy_array, ncols, dump_number)
 deallocate(columns)
end subroutine bound_unbound_thermo


!----------------------------------------------------------------
!+
!  Gravitational drag
!+
!----------------------------------------------------------------
subroutine gravitational_drag(time,npart,particlemass,xyzh,vxyzu)
 use prompting, only:prompt
 use vectorutils, only:cross_product3D
 integer, intent(in)                   :: npart
 real,    intent(in)                   :: time,particlemass
 real,    intent(inout)                :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable        :: columns(:)
 character(len=17)                     :: filename
 integer                               :: i,j,k,ncols,sizeRcut,vol_npart,merge_ij(nptmass),merge_n
 integer, allocatable                  :: iorder(:)
 real, dimension(:), allocatable, save :: ang_mom_old,time_old
 real, dimension(:,:), allocatable     :: drag_force
 real, dimension(4,maxptmass)          :: fxyz_ptmass,fxyz_ptmass_sinksink
 real, dimension(3)                    :: avg_vel,avg_vel_par,avg_vel_perp,&
                                          com_xyz,com_vxyz,unit_vel,unit_vel_perp,&
                                          pos_wrt_CM,vel_wrt_CM,ang_mom,com_vec,&
                                          unit_sep,unit_sep_perp,vel_contrast_vec,Fgrav
 real                                  :: drag_perp,drag_par,drag_perp_proj,&
                                          vel_contrast,mdot,sep,Jdot,R2,&
                                          rho_avg,cs,racc,fonrmax,fxi,fyi,fzi,&
                                          phii,phitot,dtsinksink,interior_mass,sinksinksep,&
                                          volume,vol_mass,vKep,omega,maxsep,cos_psi,mass_coregas,&
                                          com_sink_sep,Fgrav_mag
 real, dimension(:), allocatable       :: Rcut
 real, dimension(:,:,:), allocatable   :: force_cut_vec
 logical, save                         :: iacc,icentreonCM
 integer, save                         :: iavgopt

 ! OPTIONS
 if (dump_number == 0) then
    print*,'Options for averaging gas properties:'
    print "(6(/,a))",'1. Average over sphere centred on the companion (not recommended)',&
                     '2. Average over sphere centred on opposite side of orbit',&
                     '3. Average over annulus',&
                     '4. Average over annulus but excluding sphere centred on the companion',&
                     '5. Average over sphere twice as far on the opposite side of the orbit',&
                     '6. Average over sphere half as far on the opposite side of the orbit'
    iavgopt = 2
    call prompt('Select option above : ',iavgopt,1,6)
    icentreonCM = .false.
    select case (iavgopt)
    case(2,5,6)
       call prompt('Centre averaging sphere on the CM (otherwise, centre on primary core)?: ',icentreonCM)
    case(3,4)
       call prompt('Centre annulus on the CM (otherwise, centre on primary core)?: ',icentreonCM)
    end select

    write(*,"(a,i2)") 'Using ieos = ',ieos
    if ( xyzmh_ptmass(ihacc,2) > 0 ) then
       write(*,"(a,f13.7,a)") 'Companion has accretion radius = ', xyzmh_ptmass(ihacc,2), '(code units)'
       write(*,"(a)") 'Will analyse accretion'
       iacc = .true.
    else
       iacc = .false.
    endif
 endif

 ncols = 31
 allocate(columns(ncols),iorder(npart),force_cut_vec(4,maxptmass,5))
 allocate(drag_force(ncols,nptmass))
 columns = (/'   drag_perp', & ! 1  Component of net force (excluding sink-sink) perpendicular to sink separation (projection on (r2-r1) x z)
             '    drag_par', & ! 2  Component of net force (excluding sink-sink) projected along sink separation, -(r2-r1)
             'drag_perp_pr', & ! 3  'drag_perp' projected along the -v direction
             '     F_dot_v', & ! 4  Dot product of 'drag_perp_pr' and sink velocity (<0 means energy dissipation)
             ' drag_torque', & ! 5  torque / r of sink
             '     cos_psi', & ! 6  Cosine of angle between (r2-r1) x z and -v
             '       Fgrav', & ! 7  Magnitude of gravitational force from core and gas around it inferred from net force minus drag
             'mass_coregas', & ! 8  Mass of core+gas inferred from net force minus drag
             '    drag_BHL', & ! 9  Bondi-Hoyle-Lyttleton drag force
             '    mdot_BHL', & ! 10 Bond-Hoyle-Lyttleton mass accretion rate
             '       v_con', & ! 11 Magnitude of average background gas velocity minus sink velocity, positive when vsink dot vgas < 0
             '   v_con_par', & ! 12 Projection of velocity contrast on -vsink
             '       v_Kep', & ! 13 Keplerian velocity of companion, sqrt(M(<r) / |r2-r1|)
             ' mass_inside', & ! 14 Mass interior to current companion radius
             '  donor_spin', & ! 15 Spin angular velocity of gas near donor along z-direction
             ' sound_speed', & ! 16 Averaged sound speed
             ' rho_at_sink', & ! 17 Averaged gas density
             '        Racc', & ! 18 Bondi-Hoyle-Lyttleton accretion radius
             'com_sink_sep', & ! 19 Separation between companion and CoM calculated via gas+core gravity
             'orbitcom_sep', & ! 20 Separation between companion and CoM calculated with orbit_com subroutine
             ' sinksinksep', & ! 21 Sink-sink separation
             '   par_cut_1', & ! 22 Same as 'par. drag', but limited to particles within some radius from the sink
             '  perp_cut_1', & ! 23 Same as 'perp. drag', but limited to particles within some radius from the sink
             '   par_cut_2', & ! 24
             '  perp_cut_2', & ! 25
             '   par_cut_3', & ! 26
             '  perp_cut_3', & ! 27
             '   par_cut_4', & ! 28
             '  perp_cut_4', & ! 29
             '   par_cut_5', & ! 30
             '  perp_cut_5' /) ! 31

 drag_force = 0.

 do i = 1,2
    ! Note: The analysis below is performed for both the companion (i=2) and the donor core (i=1). We
    !       comment on the case of the companion only for clarity.

    ! Initialise output
    rho_avg          = 0.
    mdot             = 0.
    avg_vel          = 0.
    avg_vel_par      = 0.
    avg_vel_perp     = 0.
    vel_contrast     = 0.
    vel_contrast_vec = 0.
    racc             = 0.
    cs               = 0.
    vol_mass         = 0.
    omega            = 0.

    ! Calculate unit vectors
    call unit_vector(vxyz_ptmass(1:3,i), unit_vel)
    call cross_product3D(unit_vel, (/0.,0.,1./), unit_vel_perp)                     ! Direction perpendicular to sink velocity (pointing away from the CoM)
    call unit_vector( xyzmh_ptmass(1:3,i) - xyzmh_ptmass(1:3,3-i) , unit_sep)       ! Unit vector of r2-r1
    call cross_product3D((/0.,0.,1./), unit_sep, unit_sep_perp)                     ! z x (r2-r1)

    ! Calculate orbit CoM by calculating the CoM of the stellar cores plus with the inclusion
    ! of a number of particles around the primary.
    call orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)
    ! Calculate z-angular momentum of the sink about the orbital CoM in first dump that is analysed
    if (dump_number == 0) then
       if (i == 1) then ! Ensures only allocated once throughout the analysis
          allocate(ang_mom_old(nptmass))
          allocate(time_old(nptmass))
       endif
       pos_wrt_CM = xyzmh_ptmass(1:3,i) - com_xyz(1:3)
       vel_wrt_CM = vxyz_ptmass(1:3,i) - com_vxyz(1:3)
       call cross_product3D(pos_wrt_CM, xyzmh_ptmass(4,i) * vel_wrt_CM, ang_mom)
       ang_mom_old(i) = ang_mom(3)
       time_old = -50. ! Denotes time difference between (full) dumps, s.t. time - time_old is time in current dump
       ! This should actually be -dtmax in the infile
    endif


    ! Calculate volume averages
    call average_in_vol(xyzh,vxyzu,npart,particlemass,com_xyz,com_vxyz,i,icentreonCM,iavgopt,avg_vel,cs,omega,volume,vol_mass,&
                     vol_npart)
    if (vol_npart > 0.) then
       rho_avg           = vol_mass / volume
       avg_vel_par(1:3)  = dot_product(avg_vel, unit_vel) * unit_vel
       avg_vel_perp(1:3) = avg_vel(1:3) - avg_vel_par(1:3)
       vel_contrast_vec  = avg_vel - vxyz_ptmass(1:3,i)
       vel_contrast      = sign( distance(vel_contrast_vec), -dot_product(vxyz_ptmass(1:3,i), avg_vel) )
       racc              = 2. * xyzmh_ptmass(4,i) / (vel_contrast**2 + cs**2) ! Accretion radius
       mdot              = 4.*pi * xyzmh_ptmass(4,i)**2 * rho_avg / (cs**2 + vel_contrast**2)**1.5 ! BHL mass accretion rate
    endif


    ! Sum acceleration (fxyz_ptmass) on companion due to gravity of gas particles
    force_cut_vec = 0.
    fxyz_ptmass = 0.
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,phitot,dtsinksink,0,0.,merge_ij,merge_n)

    sizeRcut = 5
    if (i == 1) allocate(Rcut(sizeRcut))
    call logspace(Rcut,0.4,2.5)
    !Rcut = Rcut * racc ! Bin by fraction of accretion radius
    Rcut = Rcut * separation( xyzmh_ptmass(1:3,1), xyzmh_ptmass(1:3,2) ) ! Bin by fraction of sink-sink separation

    do j = 1,npart
       if (.not. isdead_or_accreted(xyzh(4,j))) then
          ! Get total gravitational force from gas
          call get_accel_sink_gas(nptmass,xyzh(1,j),xyzh(2,j),xyzh(3,j),xyzh(4,j),xyzmh_ptmass,&
                                  fxi,fyi,fzi,phii,particlemass,fxyz_ptmass,fonrmax)
          ! Get force from gas within distance cutoff
          do k = 1,sizeRcut
             if ( separation(xyzh(1:3,j), xyzmh_ptmass(1:4,i)) < Rcut(k) ) then
                call get_accel_sink_gas(nptmass,xyzh(1,j),xyzh(2,j),xyzh(3,j),xyzh(4,j),xyzmh_ptmass,&
                                        fxi,fyi,fzi,phii,particlemass,force_cut_vec(1:4,:,k),fonrmax)
             endif
          enddo
       endif
    enddo

    ! Calculate angular momentum of companion wrt orbit CoM
    pos_wrt_CM = xyzmh_ptmass(1:3,i) - com_xyz(1:3)
    vel_wrt_CM = vxyz_ptmass(1:3,i) - com_vxyz(1:3)
    call cross_product3D(pos_wrt_CM, xyzmh_ptmass(4,i) * vel_wrt_CM, ang_mom)
    Jdot = (ang_mom(3) - ang_mom_old(i)) / (time - time_old(i)) ! Average change in angular momentum
    R2 = distance(xyzmh_ptmass(1:3,i) - com_xyz(1:3))
    ang_mom_old(i) = ang_mom(3) ! Set ang_mom_old for next dump
    time_old(i) = time

    ! Calculate mass interior to companion
    call set_r2func_origin(xyzmh_ptmass(1,3-i),xyzmh_ptmass(2,3-i),xyzmh_ptmass(3,3-i)) ! Order particles by distance from donor core
    call indexxfunc(npart,r2func_origin,xyzh,iorder)
    sinksinksep = separation(xyzmh_ptmass(1:3,1), xyzmh_ptmass(1:3,2))
    interior_mass = xyzmh_ptmass(4,3-i) ! Include mass of donor core
    select case(iavgopt)
    case(5)       ! Calculate mass interior to R/2
       maxsep = 2.*sinksinksep
    case(6)       ! Calculate mass interior to 2*R
       maxsep = 0.5*sinksinksep
    case default  ! Calculate mass interior to R
       maxsep = sinksinksep
    end select
    do j = 1,npart
       k = iorder(j)
       sep = separation(xyzmh_ptmass(1:3,3-i), xyzh(1:3,k))
       if (sep > maxsep) exit
       interior_mass = interior_mass + particlemass
    enddo
    vKep = sqrt(interior_mass / sinksinksep)

    ! Calculate perpendicular force projected along -v
    cos_psi        = cos_vector_angle(-unit_sep_perp, -vxyz_ptmass(1:3,i))                ! Theta is angle between (r2-r1) x z and -v
    drag_par       = - dot_product(fxyz_ptmass(1:3,i),unit_sep) * xyzmh_ptmass(4,i)       ! Total force projected along -(r2-r1)
    drag_perp      = dot_product(fxyz_ptmass(1:3,i),-unit_sep_perp) * xyzmh_ptmass(4,i)   ! Total force projected along -(r2-r1) x z
    drag_perp_proj = drag_perp / cos_psi                                                  ! Perpendicular force projected along -v

    ! Calculate core + gas mass based on projected gravitational force
    Fgrav = fxyz_ptmass(1:3,i) * xyzmh_ptmass(4,i) - drag_perp_proj * (-unit_vel)                               ! Ftot,gas + Fsinksink = Fdrag + Fgrav
    call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass_sinksink,phitot,dtsinksink,0,0.,merge_ij,merge_n)
    Fgrav = Fgrav + fxyz_ptmass_sinksink(1:3,i) * xyzmh_ptmass(4,i)
    Fgrav_mag = distance(Fgrav)
    mass_coregas = Fgrav_mag * sinksinksep**2 / xyzmh_ptmass(4,i)

    ! Calculate CoM inferred from core + gas mass
    com_vec = (mass_coregas * xyzmh_ptmass(1:3,3-i) + xyzmh_ptmass(4,i) * xyzmh_ptmass(1:3,i)) / (mass_coregas + xyzmh_ptmass(4,i))
    com_sink_sep = separation(com_vec, xyzmh_ptmass(1:3,i))

    drag_force(1,i)  = drag_perp
    drag_force(2,i)  = drag_par
    drag_force(3,i)  = drag_perp_proj
    drag_force(4,i)  = drag_perp_proj * (-distance(vxyz_ptmass(1:3,i)))
    drag_force(5,i)  = Jdot / R2
    drag_force(6,i)  = cos_psi
    drag_force(7,i)  = Fgrav_mag
    drag_force(8,i)  = mass_coregas
    drag_force(9,i)  = mdot * vel_contrast ! BHL drag force
    drag_force(10,i) = mdot
    drag_force(11,i) = vel_contrast
    drag_force(12,i) = dot_product(vel_contrast_vec, -unit_vel)
    drag_force(13,i) = vKep
    drag_force(14,i) = interior_mass
    drag_force(15,i) = omega
    drag_force(16,i) = cs
    drag_force(17,i) = rho_avg
    drag_force(18,i) = racc
    drag_force(19,i) = com_sink_sep
    drag_force(20,i) = separation(com_xyz(1:3),xyzmh_ptmass(1:3,i))
    drag_force(21,i) = sinksinksep
    drag_force(22,i) = - dot_product(force_cut_vec(1:3,i,1),unit_sep)      * xyzmh_ptmass(4,i)
    drag_force(23,i) = - dot_product(force_cut_vec(1:3,i,1),unit_sep_perp) * xyzmh_ptmass(4,i)
    drag_force(24,i) = - dot_product(force_cut_vec(1:3,i,2),unit_sep)      * xyzmh_ptmass(4,i)
    drag_force(25,i) = - dot_product(force_cut_vec(1:3,i,2),unit_sep_perp) * xyzmh_ptmass(4,i)
    drag_force(26,i) = - dot_product(force_cut_vec(1:3,i,3),unit_sep)      * xyzmh_ptmass(4,i)
    drag_force(27,i) = - dot_product(force_cut_vec(1:3,i,3),unit_sep_perp) * xyzmh_ptmass(4,i)
    drag_force(28,i) = - dot_product(force_cut_vec(1:3,i,4),unit_sep)      * xyzmh_ptmass(4,i)
    drag_force(29,i) = - dot_product(force_cut_vec(1:3,i,4),unit_sep_perp) * xyzmh_ptmass(4,i)
    drag_force(30,i) = - dot_product(force_cut_vec(1:3,i,5),unit_sep)      * xyzmh_ptmass(4,i)
    drag_force(31,i) = - dot_product(force_cut_vec(1:3,i,5),unit_sep_perp) * xyzmh_ptmass(4,i)

    ! Write to output
    write (filename, "(A16,I0)") "sink_drag_", i
    call write_time_file(trim(adjustl(filename)), columns, time, drag_force(:,i), ncols, dump_number)
 enddo
 deallocate(columns,drag_force,force_cut_vec,Rcut)

end subroutine gravitational_drag


subroutine J_E_plane(num,npart,particlemass,xyzh,vxyzu)
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart,num
 real,    intent(in) :: particlemass,xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 integer :: ncols,i
 real :: com_xyz(3),com_vxyz(3),dum1,dum2,dum3,dum4,etoti,angmom_com(3),angmom_core(3)
 real, allocatable :: data(:,:)

 ncols = 7
 allocate(columns(ncols),data(ncols,npart))
 columns = (/'          E',&
             '      Jxcom',&
             '      Jycom',&
             '      Jzcom',&
             '     Jxcore',&
             '     Jycore',&
             '     Jzcore'/)

 call get_centreofmass(com_xyz,com_vxyz,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 do i=1,npart
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,dum1,dum2,dum3,dum4,etoti)
    data(1,i) = etoti
    call cross_product3D(xyzh(1:3,i)-xyzmh_ptmass(1:3,1), vxyzu(1:3,i)-vxyz_ptmass(1:3,1), angmom_core)
    data(5:7,i) = angmom_core
    call cross_product3D(xyzh(1:3,i)-com_xyz(1:3), vxyz_ptmass(1:3,i)-com_vxyz(1:3), angmom_com)
    data(2:4,i) = angmom_com
 enddo

 data(1,:) = data(1,:) / particlemass ! specific energy

 call write_file('JEplane','JEplane',columns,data,size(data(1,:)),ncols,num)
 deallocate(columns,data)

end subroutine J_E_plane

!-------------------------------------------------------------------
!+
! Planet destruction
!+
!-------------------------------------------------------------------
subroutine planet_destruction(time,npart,particlemass,xyzh,vxyzu)
 use kernel, only:wkern
 integer, intent(in)              :: npart
 real, intent(in)                 :: time,particlemass
 real, intent(inout)              :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable   :: columns(:)
 character(len=18)                :: filename
 real, allocatable                :: planetDestruction(:)
 integer                          :: ncols,i,j
 real, save                       :: time_old
 real, allocatable, save          :: particleRho(:)
 character(len=50)                :: planetRadiusPromptString
 real, allocatable, save          :: planetRadii(:) !In units of Rsun

 real, dimension(3)               :: currentGasVel, currentVelContrast
 real                             :: currentRho(1) !Is a one element array because sphInterpolation returns a 1 dimensional array.
 real                             :: currentRhoScaled,currentVelContrastScaled,currentPlanetRhoScaled
 real                             :: currentPlanetMassScaled,currentPlanetRadiusScaled
 real, allocatable, save          :: currentKhAblatedMass(:)

 ncols=5
 allocate(columns(ncols),planetDestruction(ncols))
 columns=(/"      rhoGas", &
           "  kh_rhoCrit", &
           "     kh_lmax", &
           "     kh_mdot", &
           " kh_ablatedM" /)

 !Kelvin-Helmholtz instability planet destruction as described in "On the survival of brown dwarfs
 !and planets by their giant host star" (https://arxiv.org/abs/1210.0879). Description of columns:
 !rhoGas: background gas density at sink. In units of g/cm^3.
 !kh_rhoCrit: paper equation 5. In units of g/cm^3.
 !kh_lmax: paper equation 6. In units of Jupiter radii.
 !kh_mdot: paper equation 7. In units of Jupiter mass/year.
 !kh_ablatedM: kh_mdot integrated over time. In units of Jupiter masses.

 currentRho = 0.
 do i=1,nptmass
    if (i==1) cycle !The first sink is assumed to be the core.

    if ((dump_number==0) .and. (i==2)) then !This is only done once.
       allocate(planetRadii(nptmass))
       planetRadii=0.1
       do j=2,nptmass
          write(planetRadiusPromptString,"(A13,I0,A32)") "Enter planet ",j-1," radius in units of solar radii"
          call prompt(planetRadiusPromptString,planetRadii(i),0.0,1.0)
       enddo

       allocate(particleRho(npart))
       allocate(currentKhAblatedMass(nptmass))

       time_old=0.0
       particleRho=getParticleRho(xyzh(4,:),particlemass)
       currentKhAblatedMass=0.0
    endif


    currentRho=sphInterpolation(npart,particlemass,particleRho,xyzh,xyzmh_ptmass(1:3,i),reshape(particleRho,(/1,npart/)))
    currentGasVel=sphInterpolation(npart,particlemass,particleRho,xyzh,xyzmh_ptmass(1:3,i),vxyzu(1:3,:))
    currentVelContrast=vxyz_ptmass(1:3,i)-currentGasVel

    currentPlanetRadiusScaled=planetRadii(i)/0.1 !In units of 0.1 Rsun.
    currentPlanetMassScaled=xyzmh_ptmass(4,i)*104.74 !In units of 10 jupiter masses.
    currentPlanetRhoScaled=(xyzmh_ptmass(4,i)/((4.0/3.0)*pi*(planetRadii(i)**3.0)))*0.44 !In units of 13.34 g/cm^3
    currentRhoScaled=currentRho(1)*59000.0 !In units of 10^-4 g/cm^3.
    currentVelContrastScaled=distance(currentVelContrast)*4.37 !In units of 100 km/s.

    planetDestruction(1)=currentRho(1)*5.9
    planetDestruction(2)=3.82*(currentPlanetRhoScaled**(4.0/3.0))*(currentPlanetMassScaled**(2.0/3.0))&
                             *(currentVelContrastScaled**(-2.0))
    planetDestruction(3)=0.0000263*(currentVelContrastScaled**2.0)*currentRhoScaled*(currentPlanetRhoScaled**((-5.0)/3.0))&
                                   *(currentPlanetMassScaled**((-1.0)/3.0))
    planetDestruction(4)=11.0*currentVelContrastScaled*currentRhoScaled*(currentPlanetRadiusScaled**2.0)&
                             *(planetDestruction(3)/(currentPlanetRadiusScaled*0.973))

    currentKhAblatedMass(i)=currentKhAblatedMass(i)+((time-time_old)*planetDestruction(4)*0.0000505)
    planetDestruction(5)=currentKhAblatedMass(i)


    write(filename, "(A17,I0)") "sink_destruction_",i
    call write_time_file(filename, columns, time, planetDestruction, ncols, dump_number)
 enddo

 time_old=time

 deallocate(columns,planetDestruction)
end subroutine planet_destruction

!-----------------------------------------------------------------------------------------
!+
!Binding energy profile
!+
!-----------------------------------------------------------------------------------------
subroutine create_bindingEnergy_profile(time,num,npart,particlemass,xyzh,vxyzu)
 real, intent(in)     :: time,particlemass
 integer, intent(in)  :: num,npart
 real, intent(in)     :: xyzh(4,npart),vxyzu(4,npart)

 character(len=17), allocatable :: columns(:)
 real, allocatable              :: profile(:,:)
 integer                        :: ncols,i,j
 integer, allocatable           :: iorder(:)
 real                           :: currentInteriorMass,currentParticleGPE,currentCoreParticleSeparation
 real                           :: previousBindingEnergy,previousBindingEnergyU

 ncols=3
 allocate(columns(ncols),iorder(npart))
 allocate(profile(ncols,npart))
 columns=(/"      radius",&
           "     bEnergy",& !Binding energy without internal energy.
           " bEnergy (u)"/) !Binding energy with internal energy.


 call set_r2func_origin(xyzmh_ptmass(1,1),xyzmh_ptmass(2,1),xyzmh_ptmass(3,1))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
 currentInteriorMass=xyzmh_ptmass(4,1)+(npart*particlemass) !Initally set to the entire mass of the star.

 do i=npart,1,-1 !Loops over all particles from outer to inner.
    j=iorder(i)
    currentInteriorMass=currentInteriorMass-particlemass
    currentCoreParticleSeparation=separation(xyzmh_ptmass(1:3,1),xyzh(1:3,j))
    currentParticleGPE=(currentInteriorMass*particlemass)/currentCoreParticleSeparation

    !The binding energy at a particular radius is the sum of the gravitational potential energies
    !(and internal energies in the case of the third column) of all particles beyond that radius.
    if (i==npart) then
       previousBindingEnergy=0.0
       previousBindingEnergyU=0.0
    else
       previousBindingEnergy=profile(2,i+1)
       previousBindingEnergyU=profile(3,i+1)
    endif

    profile(1,i)=currentCoreParticleSeparation
    profile(2,i)=previousBindingEnergy+currentParticleGPE
    profile(3,i)=previousBindingEnergyU+currentParticleGPE-(vxyzu(4,j)*particlemass)
 enddo

 call write_file('bEnergyProfile','bEnergyProfiles',columns,profile,npart,ncols,num)
 deallocate(columns,iorder,profile)

end subroutine create_bindingEnergy_profile


subroutine get_core_gas_com(time,npart,xyzh,vxyzu)
 use sortutils, only:set_r2func_origin,r2func_origin,indexxfunc
 integer, intent(in)                   :: npart
 real,    intent(in)                   :: time
 real,    intent(inout)                :: xyzh(:,:),vxyzu(:,:)
 real                                  :: sep,maxsep,core_gas_com(3),core_gas_vcom(3),xyz_gas(4,npart),vxyz_gas(3,npart)
 real, allocatable                     :: mytable(:)
 character(len=17), allocatable        :: columns(:)
 character(len=17)                     :: filename
 integer, save                         :: ngas
 integer, allocatable, save            :: iorder(:)
 integer                               :: ncols,j,k

 ncols = 12
 allocate(columns(ncols))
 allocate(mytable(ncols))
 mytable = 0.
 columns = (/'   gas_com_x', &
             '   gas_com_y', &
             '   gas_com_z', &
             '  gas_com_vx', &
             '  gas_com_vy', &
             '  gas_com_vz', &
             '      core_x', &
             '      core_y', &
             '      core_z', &
             '     core_vx', &
             '     core_vy', &
             '     core_vz' /)


 ! Record particles that are closest to primary core
 if (dump_number == 0) then
    allocate(iorder(npart))
    maxsep = 10. ! 10 Rsun
    ngas = 0
    call set_r2func_origin(xyzmh_ptmass(1,1),xyzmh_ptmass(2,1),xyzmh_ptmass(3,1)) ! Order particles by distance from donor core
    call indexxfunc(npart,r2func_origin,xyzh,iorder)

    do j=1,npart
       k = iorder(j)
       if (j < 10) print*,k
       sep = separation(xyzmh_ptmass(1:3,1), xyzh(1:3,k))
       if (sep > maxsep) exit
       ngas = ngas + 1
    enddo
 endif

 print*,'ngas=',ngas

 do j=1,ngas
    k = iorder(j)
    xyz_gas(1:4,j)  = xyzh(1:4,k)
    vxyz_gas(1:3,j) = vxyzu(1:3,k)
 enddo

 call get_centreofmass(core_gas_com,core_gas_vcom,ngas,xyz_gas,vxyz_gas) ! Do not include sinks

 mytable(1:3)   = core_gas_com(1:3)
 mytable(4:6)   = core_gas_vcom(1:3)
 mytable(7:9)   = xyzmh_ptmass(1:3,1)
 mytable(10:12) = vxyz_ptmass(1:3,1)

 write (filename, "(A16,I0)") "core_gas_com"
 call write_time_file(trim(adjustl(filename)),columns,time,mytable,ncols,dump_number)
end subroutine get_core_gas_com


!----------------------------------------------------------------
!+
!  Print dump numbers corresponding to given sink-sink separations
!+
!----------------------------------------------------------------
subroutine print_dump_numbers(dumpfile)
 character(len=*), intent(in) :: dumpfile
 character(len=50), allocatable, save :: dumpfiles(:)
 integer :: nseps
 integer, save :: i
 real, allocatable :: sinksinksep(:)
 real :: sep

 nseps = 2
 allocate(sinksinksep(nseps))
 if (dump_number == 0) then
    allocate(dumpfiles(nseps))
    i=1
 endif
 sinksinksep = (/ 938., 67. /)

 sep = separation(xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,2))
 if ( sep < sinksinksep(i) ) then
    dumpfiles(i) = trim(dumpfile)
    i=i+1
 endif
 if (i==nseps+1) then
    print "(5(a,/))",'../',dumpfiles
    return
 endif

end subroutine print_dump_numbers


!----------------------------------------------------------------
!+
!  Analyse disk
!+
!----------------------------------------------------------------
subroutine analyse_disk(num,npart,particlemass,xyzh,vxyzu)
 use part,            only:eos_vars,itemp
 use extern_corotate, only:get_companion_force
 use ionization_mod,  only:calc_thermal_energy
 use vectorutils,     only:cross_product3D
 integer, intent(in)             :: num,npart
 real, intent(in)                :: particlemass
 real, intent(inout)             :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable  :: columns(:)
 real, allocatable               :: data(:,:)
 real                            :: diskz,diskR2,diskR1,R,omegai,phii,rhopart,ponrhoi,spsoundi,tempi,&
                                    epoti,ekini,ethi,Ji(3),vrel2,fxi,fyi,fzi,vphi
 integer                         :: ncols,i

 ncols = 9
 allocate(columns(ncols),data(ncols,npart))
 data = -1.
 columns = (/'          R',&  ! cylindrical radius w.r.t companion
             '          E',&  ! specific energy (kin+pot only) w.r.t. companion
             '      Omega',&  ! angular momentum w.r.t. companion
             '         Jx',&  ! specific angular momentum components
             '         Jy',&
             '         Jz',&
             '       ekin',&
             '       epot',&  ! gravitational potential energy due to companion only
             '     etherm'/)

 ! Set disk dimensions
 diskz  = 50.  ! disk half-thickness
 diskR1 = 5.   ! disk inner radius
 diskR2 = 150.  ! disk outer radius

 do i=1,npart
    ! Skip if particle is not within the defined disk
    if (abs(xyzh(3,i) - xyzmh_ptmass(3,2)) > diskz) cycle
    R = sqrt( (xyzh(1,i) - xyzmh_ptmass(1,2))**2 + (xyzh(2,i) - xyzmh_ptmass(2,2))**2 )
    if ( (R > diskR2) .or. (R < diskR1) ) cycle

    vrel2 = (vxyzu(1,i) - vxyz_ptmass(1,2))**2 + (vxyzu(2,i) - vxyz_ptmass(2,2))**2 + (vxyzu(3,i) - vxyz_ptmass(3,2))**2
    ekini = 0.5*particlemass*vrel2

    ! Calculate gravitational potential due to companion only
    phii = 0.
    call get_companion_force(xyzh(1:3,i),fxi,fyi,fzi,phii)
    epoti = phii*particlemass

    ! Calculate thermal energy
    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))
    call calc_thermal_energy(particlemass,ieos,xyzh(:,i),vxyzu(:,i),ponrhoi*rhopart,eos_vars(itemp,i),gamma,ethi)

    call get_gas_omega(xyzmh_ptmass(1:3,2),vxyz_ptmass(1:3,2),xyzh(1:3,i),vxyzu(1:3,i),vphi,omegai)
    call cross_product3D(xyzh(1:3,i)-xyzmh_ptmass(1:3,2), vxyzu(1:3,i)-vxyz_ptmass(1:3,2), Ji)

    data(1,i) = R
    data(2,i) = (ekini+epoti) / particlemass
    data(3,i) = omegai
    data(4:6,i) = Ji
    data(7,i) = ekini
    data(8,i) = epoti
    data(9,i) = ethi
 enddo
 call write_file('companion_disk','companion_disk',columns,data,npart,ncols,num)
 deallocate(columns)

end subroutine analyse_disk


!----------------------------------------------------------------
!+
!  Recombination energy vs. time
!+
!----------------------------------------------------------------
subroutine erec_vs_t(time,npart,particlemass,xyzh)
 use ionization_mod, only:get_erec_components
 integer, intent(in) :: npart
 real, intent(in)    :: time,particlemass
 real, intent(inout) :: xyzh(:,:)
 character(len=17)   :: filename,columns(4)
 integer             :: i
 real                :: ereci(4),erec(4),tempi,rhoi

 columns = (/'          H2', &
             '          HI', &
             '         HeI', &
             '        HeII'/)

 erec = 0.
 do i = 1,npart
    rhoi = rhoh(xyzh(4,i), particlemass)
    call get_erec_components( log10(rhoi*unit_density), tempi, X_in, 1.-X_in-Z_in, ereci)
    erec = erec + ereci
 enddo

 write (filename, "(A16,I0)") "erec_vs_t"
 call write_time_file(trim(adjustl(filename)),columns,time,erec,4,dump_number)

end subroutine erec_vs_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!        Routines used in analysis routines        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------
!+
!  Calculate the angular velocity of an envelope gas particle
!  relative to a reference point
!+
!----------------------------------------------------------------
subroutine get_gas_omega(xyz_centre,vxyz_centre,xyzi,vxyzi,vphi,omega)
 use vectorutils, only:cross_product3D
 real, intent(in)  :: xyz_centre(3),vxyz_centre(3),xyzi(3),vxyzi(3)
 real, intent(out) :: vphi,omega
 real              :: Rmag,R(3),phi_unitvec(3),R_unitvec(3)

 ! xyz_centre: Position vector of reference point
 ! vxyz_centre: Velocity vector of reference point
 ! R: Cylindrical radius vector
 R(1:2) = xyzi(1:2) - xyz_centre(1:2) ! Separation in x-y plane
 R(3) = 0.
 Rmag = sqrt(dot_product(R,R))
 R_unitvec = R / Rmag
 call cross_product3D((/0.,0.,1./), R_unitvec, phi_unitvec) ! phi = z x R
 vphi = dot_product(vxyzi - vxyz_centre, phi_unitvec)
 omega = vphi / Rmag
end subroutine get_gas_omega


!----------------------------------------------------------------
!+
!  Calculate kinetic, gravitational potential (gas-gas and sink-gas),
!  and internal energy of a gas particle.
!+
!----------------------------------------------------------------
subroutine calc_gas_energies(particlemass,poten,xyzh,vxyzu,xyzmh_ptmass,phii,epoti,ekini,einti,etoti)
 ! Warning: Do not sum epoti or etoti as it is to obtain a total energy; this would not give the correct
 !          total energy due to complications related to double-counting.
 use ptmass, only:get_accel_sink_gas
 use part,   only:nptmass
 real, intent(in)                       :: particlemass
 real(4), intent(in)                    :: poten
 real, dimension(4), intent(in)         :: xyzh,vxyzu
 real, dimension(5,nptmass), intent(in) :: xyzmh_ptmass
 real, intent(out)                      :: phii,epoti,ekini,einti,etoti
 real                                   :: fxi,fyi,fzi

 phii = 0.0

 call get_accel_sink_gas(nptmass,xyzh(1),xyzh(2),xyzh(3),xyzh(4),xyzmh_ptmass,fxi,fyi,fzi,phii)

 epoti = 2.*poten + particlemass * phii ! For individual particles, need to multiply 2 to poten to get \sum_j G*mi*mj/r
 ekini = particlemass * 0.5 * dot_product(vxyzu(1:3),vxyzu(1:3))
 einti = particlemass * vxyzu(4)
 etoti = epoti + ekini + einti
end subroutine calc_gas_energies


subroutine adjust_corotating_velocities(npart,particlemass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,omega_c,dump_number)
 use vectorutils, only:cross_product3D
 real, dimension(:,:), intent(in)    :: xyzmh_ptmass,xyzh
 real, dimension(:,:), intent(inout) :: vxyzu,vxyz_ptmass
 real, intent(inout) :: omega_c
 real, intent(in)    :: particlemass
 integer, intent(in) :: npart, dump_number

 logical             :: switch
 real                :: sep, mtot
 real, dimension(3)  :: omega_vec, omegacrossr
 integer             :: i

 if (dump_number == 0) then
    call prompt('Was this in a corotating frame?',switch,.false.)

    if (switch) then
       sep = separation(xyzmh_ptmass(1:3,1), xyzmh_ptmass(1:3,2))
       mtot = sum(xyzmh_ptmass(4,:)) + npart*particlemass
       omega_c = sqrt(mtot / sep**3)
    else
       omega_c = -1
    endif
 endif

 if (omega_c > 0.) then
    omega_vec = (/ 0.,0.,omega_c /)

    do i=1,npart
       call cross_product3D(omega_vec,xyzh(1:3,i),omegacrossr)
       vxyzu(1:3,i) = vxyzu(1:3,i) + omegacrossr(1:3)
    enddo

    do i=1,nptmass
       call cross_product3D(omega_vec,xyzmh_ptmass(1:3,i),omegacrossr)
       vxyz_ptmass(1:3,i) = vxyz_ptmass(1:3,i) + omegacrossr(1:3)
    enddo
 endif
end subroutine adjust_corotating_velocities


! returns a profile from the centre of mass
! profile can either use all particles or can find particles within 2h of a given ray
! if simple flag is set to true, it will only produce a limited subset
subroutine stellar_profile(time,ncols,particlemass,npart,xyzh,vxyzu,profile,simple,ray)
 use eos,          only:ieos,equationofstate,X_in, Z_in
 use eos_mesa,     only:get_eos_kappa_mesa,get_eos_pressure_temp_mesa
 use physcon,      only:kboltz,mass_proton_cgs
 use centreofmass, only:get_centreofmass
 use energies,     only:compute_energies
 use part,         only:xyzmh_ptmass,rhoh,ihsoft,poten
 use units,        only:udist,unit_ergg,unit_density,unit_pressure,unit_velocity,unit_energ
 use kernel,       only:kernel_softening,radkern
 use ptmass,       only:get_accel_sink_gas
 use ionization_mod, only:ionisation_fraction

 real,    intent(in)    :: time
 integer, intent(in)    :: ncols
 real,    intent(in)    :: particlemass
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real, intent(out), allocatable :: profile(:,:)
 logical, intent(in)    :: simple
 real, intent(in), optional :: ray(3)
 integer           :: i,iprofile
 real              :: proj(3),orth(3),proj_mag,orth_dist,orth_ratio
 real              :: rhopart,ponrhoi,spsoundi,tempi
 real              :: temp,kappa,kappat,kappar,pres
 real              :: ekini,epoti,einti,etoti,phii
 real              :: xh0, xh1, xhe0, xhe1, xhe2
 real              :: temp_profile(ncols,npart)
 logical           :: criteria

 call compute_energies(time)

 iprofile = 0

 do i=1,npart
    if (xyzh(4,i)  >=  0) then

       if (present(ray)) then
          proj_mag = dot_product(xyzh(1:3,i),ray(1:3))
          proj = proj_mag * ray
          orth(1:3) = xyzh(1:3,i) - proj(1:3)
          orth_dist = separation(orth,(/0.,0.,0./))
          orth_ratio = orth_dist / xyzh(4,i)
          if (orth_ratio < radkern .and. proj_mag > 0.) then
             criteria = .true.
          else
             criteria = .false.
          endif
       else
          criteria = .true.
       endif

       if (criteria) then

          iprofile = iprofile + 1

          rhopart = rhoh(xyzh(4,i), particlemass)

          temp_profile(1,iprofile)  = distance(xyzh(1:3,i)) * udist
          temp_profile(3,iprofile)  = atan2(xyzh(2,i),xyzh(1,i))
          temp_profile(4,iprofile)  = rhopart * unit_density
          temp_profile(5,iprofile)  = distance(vxyzu(1:3,i)) * unit_velocity
          temp_profile(6,iprofile)  = dot_product(vxyzu(1:3,i),xyzh(1:3,i)) / distance(xyzh(1:3,i)) * unit_velocity
          temp_profile(7,iprofile)  = sqrt(distance(vxyzu(1:2,i))**2 - (dot_product(vxyzu(1:2,i),xyzh(1:2,i)) &
                                      / distance(xyzh(1:2,i)))**2) * unit_velocity
          temp_profile(8,iprofile)  = temp_profile(7,iprofile) / (distance(xyzh(1:2,i)) * udist)
          if (simple .eqv. .false.) then
             call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),tempi,vxyzu(4,i))

             if (ieos == 10) then
                call get_eos_pressure_temp_mesa(rhopart*unit_density,vxyzu(4,i) * unit_ergg,pres,temp)
                call get_eos_kappa_mesa(rhopart*unit_density,temp,kappa,kappat,kappar)
             else
                temp = (ponrhoi * (unit_pressure/unit_density) * 2.381 * mass_proton_cgs) / kboltz
                kappa = 1.
             endif

             call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

             call ionisation_fraction(rhopart*unit_density,temp,X_in,1.-X_in-Z_in,xh0,xh1,xhe0,xhe1,xhe2)

             temp_profile(9,iprofile)  = vxyzu(4,i) * unit_ergg
             temp_profile(10,iprofile) = ponrhoi * rhopart * unit_pressure
             temp_profile(11,iprofile) = spsoundi * unit_velocity
             temp_profile(12,iprofile) = temp
             temp_profile(13,iprofile) = kappa
             temp_profile(14,iprofile) = 1. / (kappa * rhopart * unit_density)
             temp_profile(15,iprofile) = etoti * unit_energ
             temp_profile(16,iprofile) = xh1
             temp_profile(17,iprofile) = xhe1
             temp_profile(18,iprofile) = xhe2
          endif
       endif
    endif
 enddo

 allocate(profile(ncols,iprofile))
 profile(1:ncols,1:iprofile) = temp_profile(1:ncols,1:iprofile)

 call quicksort(profile, 1, iprofile, ncols, 1)

 do i=1,iprofile
    if (i==1) profile(2,i) = particlemass
    if (i > 1) profile(2,i) = profile(2,i-1) + particlemass
 enddo

 deallocate(profile)
 print*, "Profile completed"

end subroutine stellar_profile

!----------------------------------------------------------------
!+
!  Calculate mass interior to companion
!+
!----------------------------------------------------------------
subroutine get_interior_mass(xyzh,vxyzu,donor_xyzm,companion_xyzm,particlemass,npart,iavgopt,interior_mass,com_xyz,com_vxyz)
 real, intent(in) :: xyzh(:,:),vxyzu(:,:),donor_xyzm(4),companion_xyzm(4),particlemass
 real, intent(out) :: interior_mass,com_xyz(3),com_vxyz(3)
 integer, intent(in) :: npart,iavgopt
 real :: sinksinksep,maxsep,sep,xyz_int(3,npart),vxyz_int(3,npart)
 integer :: j,k,npart_int
 integer, allocatable :: iorder(:)

 ! Calculate mass interior to companion
 allocate(iorder(npart))
 call set_r2func_origin(donor_xyzm(1),donor_xyzm(2),donor_xyzm(3)) ! Order particles by distance from donor core
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
 sinksinksep = separation(donor_xyzm(1:3), companion_xyzm(1:3))
 interior_mass = donor_xyzm(4) ! Include mass of donor core
 select case(iavgopt)
 case(5)       ! Calculate mass interior to R/2
    maxsep = 2.*sinksinksep
 case(6)       ! Calculate mass interior to 2*R
    maxsep = 0.5*sinksinksep
 case default  ! Calculate mass interior to R
    maxsep = sinksinksep
 end select
 npart_int = 0
 do j = 1,npart
    k = iorder(j)
    sep = separation(donor_xyzm(1:3), xyzh(1:3,k))
    if (sep > maxsep) exit
    npart_int = npart_int + 1
    xyz_int(1:3,npart_int) = xyzh(1:3,k)
    vxyz_int(1:3,npart_int) = vxyzu(1:3,k)
 enddo
 interior_mass = npart_int * particlemass

 call get_centreofmass(com_xyz,com_vxyz,npart_int,xyz_int,vxyz_int,nptmass,xyzmh_ptmass,vxyz_ptmass)
 deallocate(iorder)

end subroutine get_interior_mass

!----------------------------------------------------------------
!+
!  Get CoM position and velocity of the two point masses plus
!  gas particles radius = 2*sep from the donor, where sep is the
!  distance between the donor and the CoM of just the point masses.
!+
!----------------------------------------------------------------
subroutine orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)
 integer, intent(in)             :: npart,nptmass
 real, intent(in)                :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(out), dimension(3) :: com_xyz,com_vxyz
 real, allocatable               :: xyz_a(:,:)
 real, allocatable               :: vxyz_a(:,:)
 integer, allocatable            :: iorder(:)
 integer                         :: npart_a
 real                            :: sep
 integer                         :: i,j,k

 allocate(iorder(npart),xyz_a(4,npart),vxyz_a(3,npart))

 ! Get order of particles by distance from CoM of point masses
 com_xyz(1) = sum(xyzmh_ptmass(1,:)*xyzmh_ptmass(4,:))/nptmass
 com_xyz(2) = sum(xyzmh_ptmass(2,:)*xyzmh_ptmass(4,:))/nptmass
 com_xyz(3) = sum(xyzmh_ptmass(3,:)*xyzmh_ptmass(4,:))/nptmass
 call set_r2func_origin(com_xyz(1),com_xyz(2),com_xyz(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
 ! Displacement of donor core from the CoM of point masses
 sep = separation(xyzmh_ptmass(1:3,1),com_xyz(1:3))

 ! Calculate CoM of orbit, including only gas particles within radius = 2*sep from donor core
 ! The point is that by including some gas particles around the donor core, we get a more accurate
 ! position of the CoM about which the stellar cores orbit
 i = 1
 k = 1
 do while (i < npart+1)
    j = iorder(i) ! Loop from particles closest to farthest from CoM
    if (isdead_or_accreted(xyzh(4,j))) then
       i = i + 1
    else
       if (separation(xyzh(1:3,j),com_xyz(1:3)) > 2.*sep) exit
       xyz_a(1:4,k)  = xyzh(1:4,j)
       vxyz_a(1:3,k) = vxyzu(1:3,j)
       i = i + 1
       k = k + 1
    endif
 enddo
 npart_a = k - 1
 call get_centreofmass(com_xyz,com_vxyz,npart_a,xyz_a,vxyz_a,nptmass,xyzmh_ptmass,vxyz_ptmass)
 deallocate(iorder,xyz_a,vxyz_a)

end subroutine orbit_com

subroutine average_in_vol(xyzh,vxyzu,npart,particlemass,com_xyz,com_vxyz,isink,icentreonCM,iavgopt,vel,cs,omega,volume,vol_mass,&
                          vol_npart)
 real,    intent(in) :: xyzh(:,:),vxyzu(:,:),com_xyz(:),com_vxyz(:),particlemass
 logical, intent(in) :: icentreonCM
 real,    intent(out) :: vel(:),cs,omega,volume,vol_mass
 integer, intent(out) :: vol_npart
 integer, intent(in) :: npart,isink,iavgopt
 real :: orbit_centre(3),orbit_centre_vel(3),sphere_centre(3),Rarray(size(xyzh(1,:))),zarray(size(xyzh(1,:))),vxyzu_copy(4)
 real :: Rsphere,sep,omega_out,Rsinksink,dR,dz,vphi
 integer :: i,j,k,iorder(size(xyzh(1,:)))

 i = isink
 if (icentreonCM) then   ! Centre on orbit CoM
    orbit_centre     = com_xyz
    orbit_centre_vel = com_vxyz
 else                     ! Centre on primary core
    orbit_centre     = xyzmh_ptmass(1:3,3-i)
    orbit_centre_vel = vxyz_ptmass(1:3,3-i)
 endif

 Rsphere = 0.2 * separation(orbit_centre, xyzmh_ptmass(1:3,i))
 Rsinksink = separation(xyzmh_ptmass(1:2,i), xyzmh_ptmass(1:2,3-i))                          ! [(x2-x1)^2 + (y2-y1)^2]^0.5
 dR = 0.2*Rsinksink
 dz = 0.2*Rsinksink
 vol_npart = 0
 vol_mass = 0.
 omega = 0.
 cs = 0.

 ! If averaging over a sphere, get order of particles from closest to farthest from sphere centre
 dr = 0.
 dz = 0.
 Rsinksink = 0.
 vol_npart = 0
 Rsphere = 0.
 select case(iavgopt)
 case(1,2,5,6)
    select case (iavgopt)
    case(1) ! Use companion position
       sphere_centre = xyzmh_ptmass(1:3,i)
    case(2) ! Use companion position on the opposite side of orbit
       sphere_centre = 2.*orbit_centre - xyzmh_ptmass(1:3,i) ! Just r1 - (r2 - r1)
    case(5) ! Averaging twice as far on opposite side of orbit
       sphere_centre = 2.*(orbit_centre - xyzmh_ptmass(1:3,i)) ! Just r1 - 2(r2 - r1)
    case(6) ! Averaging half as far on opposite side of orbit
       sphere_centre = 1.5*orbit_centre - 0.5*xyzmh_ptmass(1:3,i) ! Just r1 - 0.5*(r2 - r1)
    end select
    call set_r2func_origin(sphere_centre(1),sphere_centre(2),sphere_centre(3))
    call indexxfunc(npart,r2func_origin,xyzh,iorder)

    ! Sum velocities, cs, and densities of all particles within averaging sphere
    do j = 1,npart
       k = iorder(j) ! Only use particles within the averaging sphere
       if (.not. isdead_or_accreted(xyzh(4,k))) then
          sep = separation(xyzh(1:3,k), sphere_centre)
          if (sep > Rsphere) exit
          vel(1:3) = vel(1:3) + vxyzu(1:3,k)
          vxyzu_copy = vxyzu(:,k)
          cs       = cs + get_spsound(ieos,xyzh(1:3,k),rhoh(xyzh(4,k),particlemass),vxyzu_copy)
          call get_gas_omega(orbit_centre,orbit_centre_vel,xyzh(1:3,k),vxyzu(1:3,k),vphi,omega_out)
          omega    = omega + omega_out
       endif
    enddo
    vol_npart = j-1 ! Number of (unaccreted) particles in the sphere
    vol_mass  = vol_npart * particlemass
    if ((iavgopt == 2) .or. (iavgopt == 5) .or. (iavgopt == 6)) vel = -vel ! To-do: get rid of this line

    ! Averaging in annulus
 case(3,4)
    Rarray = sqrt( (xyzh(1,:) - xyzmh_ptmass(1,3-i))**2 + (xyzh(2,:) - xyzmh_ptmass(2,3-i))**2) ! [(x-x1)^2 + (y-y1)^2]^0.5
    zarray = xyzh(3,:) - xyzmh_ptmass(3,3-i)
    if (iavgopt == 4) Rsphere = 0.2*separation(xyzmh_ptmass(1:3,3-i),xyzmh_ptmass(1:3,i))
    do k = 1,npart
       if ( (iavgopt == 4) .and. (separation(xyzh(1:3,k), xyzmh_ptmass(1:3,i)) < Rsphere) ) cycle
       if ( (abs(Rarray(k) - Rsinksink) < 0.5*dR) .and.&
              (abs(zarray(k) - xyzmh_ptmass(3,3-i)) < 0.5*dz) ) then
          vel   = vel + vxyzu(1:3,k)
          vxyzu_copy = vxyzu(:,k)
          cs    = cs + get_spsound(ieos,xyzh(1:3,k),rhoh(xyzh(4,k),particlemass),vxyzu_copy)
          call get_gas_omega(orbit_centre,orbit_centre_vel,xyzh(1:3,k),vxyzu(1:3,k),vphi,omega_out)
          omega = omega + omega_out
          vol_npart = vol_npart + 1
       endif
    enddo
    vol_mass = vol_npart * particlemass
 end select

 ! Calculate averaging volume based on averaging option
 select case (iavgopt)
 case (1,2,5,6) ! Spheres
    volume = 4./3.*pi*Rsphere**3
 case(3) ! Annulus
    volume  = 2.*pi * Rsinksink * dR * dz
 case(4) ! Annulus with sphere subtracted
    volume  = 2.*pi * Rsinksink * dR * dz
    volume  = volume - 0.4*dR*dz*Rsinksink
 case default
    volume = 0.
    print*,'Unknown averaging option'
    return
 end select

 ! Calculate volume averages
 if (vol_npart > 0) then
    vel(1:3) = vel(1:3) / float(vol_npart)
    omega  = omega / float(vol_npart)
    cs     = cs / float(vol_npart)
 endif

end subroutine average_in_vol


!----------------------------------------------------------------
!+
!  Returns hist, the radial or mass-coordinate profile of a
!  quantity.
!
!  Inputs:
!   coord: Array of radius or mass-coordinate of each particle
!   quant: Array containing quantity for each particle to be binned
!   bin_min: Lower bin edge for coord
!   bin_max: Upper bin edge for coord
!   nbins: Number of bins for coord
!   logbins: If true, produce log-uniform bins
!   normalise_by_bincount: If true, normalises histogram by bin
!            count, thus averaging the quantity
!+
!----------------------------------------------------------------
subroutine histogram_setup(coord,quant,hist,npart,bin_max,bin_min,nbins,normalise_by_bincount,logbins)
 integer, intent(in)    :: npart,nbins
 real, intent(in)       :: coord(npart),quant(npart),bin_max, bin_min
 logical, intent(in)    :: normalise_by_bincount,logbins
 real, intent(out)      :: hist(nbins)
 integer                :: i,j,bincount(nbins)
 real                   :: bins(nbins)

 if (logbins) then ! Create log-uniform bins
    bins = (/ (10**(bin_min + (i-1) * (bin_max-bin_min)/real(nbins)), i=1,nbins) /)
 else ! Create linear bins
    bins = (/ (bin_min + (i-1) * (bin_max-bin_min)/real(nbins), i=1,nbins) /)
 endif

 hist = 0.
 bincount = 0

 do j=1,npart
    do i=1,nbins-1
       if (coord(j) >= bins(i) .and. coord(j) < bins(i+1)) then
          bincount(i) = bincount(i) + 1
          hist(i) = hist(i) + quant(j)
          exit ! Move onto next particle
       endif
    enddo
 enddo

 if (normalise_by_bincount) then
    do i=1,nbins
       if (bincount(i) > 0) then
          hist(i) = hist(i) / real(bincount(i))
       endif
    enddo
 endif

end subroutine histogram_setup

subroutine write_file(name_in, dir_in, cols, data_in, npart, ncols, num)
 !outputs a file from a single dump
 character(len=*), intent(in) :: name_in, dir_in
 integer, intent(in)          :: npart, ncols, num
 character(len=*), dimension(ncols), intent(in) :: cols
 character(len=20), dimension(ncols) :: columns
 character(len=40)             :: data_formatter, column_formatter
 character(len(name_in)+9)    :: file_name

 real, dimension(ncols,npart), intent(in) :: data_in
 integer                      :: i, unitnum

 unitnum = 1000 + num
 if (dump_number == 0) then
    call system('mkdir ' // dir_in )
 endif

 write(file_name, "(2a,i5.5,a)") trim(name_in), "_", num, ".ev"

 open(unit=unitnum, file='./'//dir_in//'/'//file_name, status='replace')

 write(column_formatter, "(a,I2.2,a)") "('#',2x,", ncols, "('[',a15,']',3x))"
 write(data_formatter, "(a,I2.2,a)") "(", ncols, "(2x,es19.11e3))"

 do i=1,ncols
    write(columns(i), "(I2.2,a)") i, cols(i)
 enddo

 !set column headings
 write(unitnum, column_formatter) columns(:)

 !Write data to file
 do i=1,npart
    write(unitnum,data_formatter) data_in(:ncols,i)
 enddo

 close(unit=unitnum)
end subroutine write_file


subroutine write_time_file(name_in, cols, time, data_in, ncols, num)
 !outputs a file over a series of dumps
 character(len=*), intent(in) :: name_in
 integer, intent(in)          :: ncols, num
 character(len=*), dimension(ncols), intent(in) :: cols
 character(len=20), dimension(ncols) :: columns
 character(len=40)             :: data_formatter, column_formatter
 character(len(name_in)+9)    :: file_name
 real, intent(in)             :: time
 real, dimension(ncols), intent(in) :: data_in
 integer                      :: i, unitnum

 write(column_formatter, "(a,I2.2,a)") "('#',2x,", ncols+1, "('[',a15,']',3x))"
 write(data_formatter, "(a,I2.2,a)") "(", ncols+1, "(2x,es18.11e2))"
 write(file_name,"(2a,i3.3,a)") name_in, '.ev'

 if (num == 0) then
    unitnum = 1000

    open(unit=unitnum, file=file_name, status='replace')
    do i=1,ncols
       write(columns(i), "(I2,a)") i+1, cols(i)
    enddo

    !set column headings
    write(unitnum, column_formatter) '1         time', columns(:)
    close(unit=unitnum)
 endif

 unitnum=1001+num

 open(unit=unitnum, file=file_name, position='append')

 write(unitnum,data_formatter) time, data_in(:ncols)

 close(unit=unitnum)

end subroutine write_time_file

real function distance(a)
 ! Return norm of a vector of arbitrary dimension
 real, intent(in), dimension(:) :: a

 distance = sqrt(dot_product(a,a))
end function distance

subroutine unit_vector(a,b)
 real, intent(in), dimension(3)  :: a
 real, intent(out), dimension(3) :: b

 b(1:3) = a(1:3) / distance(a(1:3))
end subroutine unit_vector

real function cos_vector_angle(a,b)
 real, intent(in), dimension(3) :: a,b
 if (distance(a) == 0 .or. distance(b) == 0) then
    cos_vector_angle = 1.
 else
    cos_vector_angle = dot_product(a,b) / (distance(a) * distance(b))
 endif
end function cos_vector_angle

subroutine separation_vector(a,b,c)
 !return difference between two vectors
 real, intent(in), dimension(3) :: a,b
 real, intent(out), dimension(4) :: c

 c(1) = a(1) - b(1)
 c(2) = a(2) - b(2)
 c(3) = a(3) - b(3)
 c(4) = distance(c(1:3))
end subroutine separation_vector

real function separation(a,b)
 !return the distance between two vectors
 real, intent(in), dimension(:) :: a,b

 separation = distance(a - b)
end function separation

!Creates an array of SPH particle densities for each value of h.
elemental real function getParticleRho(h,particlemass)
 real, intent(in) :: h,particlemass
 getParticleRho=rhoh(h,particlemass)
end function getParticleRho

!Performs SPH interpolation on the SPH particle property toInterpolate at the location interpolateXyz.
!The smoothing length used is the smoothing length of the closest SPH particle to interpolateXyz.
function sphInterpolation(npart,particlemass,particleRho,particleXyzh,interpolateXyz,toInterpolate) result(interpolatedData)
 use kernel, only:wkern
 integer, intent(in) :: npart
 real, intent(in)    :: particlemass
 real, intent(in)    :: particleRho(npart)
 real, intent(in)    :: particleXyzh(4,npart)
 real, intent(in)    :: interpolateXyz(3)
 real, intent(in)    :: toInterpolate(:,:)
 real                :: interpolatedData(size(toInterpolate,1))

 integer              :: i,j
 integer, allocatable :: iorder(:)
 real                :: currentR,currentQ,currentQ2
 real                :: nearestSphH
 real                :: currentParticleRho,currentSphSummandFactor

 interpolatedData=0.0
 allocate(iorder(npart))
 call set_r2func_origin(interpolateXyz(1),interpolateXyz(2),interpolateXyz(3))
 call indexxfunc(npart,r2func_origin,particleXyzh,iorder) !Gets the order of SPH particles from the interpolation point.
 nearestSphH=particleXyzh(4,iorder(1)) !The smoothing length of the nearest SPH particle to the ineterpolation point.

 do i=1,npart
    j=iorder(i)

    currentR=separation(interpolateXyz,particleXyzh(1:3,j))
    currentQ=currentR/nearestSphH !currentR is scaled in units of nearestSphH
    currentQ2=currentQ**2.0

    !All SPH particles beyond 2 smoothing lengths are ignored.
    if (currentQ>2) then
       exit
    endif

    !SPH interpolation is done below.
    currentParticleRho=particleRho(j)
    currentSphSummandFactor=(particlemass/currentParticleRho)*((1.0/((nearestSphH**3.0)*pi))*wkern(currentQ2,currentQ))
    interpolatedData=interpolatedData+(currentSphSummandFactor*toInterpolate(:,j))
 enddo
 deallocate(iorder)

end function sphInterpolation

!Sorting routines
recursive subroutine quicksort(a, first, last, ncols, sortcol)
 integer, intent(in)                                :: first, last, ncols, sortcol
 real, dimension(ncols,last-first+1), intent(inout) :: a
 real                                               :: x
 integer                                            :: i, j, k

 x = a(sortcol, (first+last) / 2 )
 i = first
 j = last
 do
    do while (a(sortcol, i) < x)
       i=i+1
    enddo

    do while (x < a(sortcol, j))
       j=j-1
    enddo

    if (i >= j) exit

    do k=1,ncols
       call swap(a(k,i),a(k,j))
    enddo

    i=i+1
    j=j-1
 enddo
 if (first < i-1) call quicksort(a, first, i-1, ncols, sortcol)
 if (j+1 < last)  call quicksort(a, j+1, last, ncols, sortcol)
end subroutine quicksort

subroutine swap(a,b)
 real, intent(inout) :: a,b
 real                :: c

 c = a
 a = b
 b = c

end subroutine swap


!----------------------------------------------------------------
!+
!  Determine ID of planet particles based on distance from host star core
!+
!----------------------------------------------------------------
subroutine get_planetIDs(nplanet,planetIDs)
 integer, allocatable, intent(out) :: planetIDs(:)
 integer, intent(out)              :: nplanet
 integer                           :: i

 ! Determine planet particle IDs (the nplanet particles initially farthest from the donor star)
 nplanet = 1262
 call prompt('Enter number of planet particles:',nplanet,0)
 allocate(planetIDs(nplanet))
 do i = 1,nplanet
    planetIDs(i) = i
 enddo

end subroutine get_planetIDs


!----------------------------------------------------------------
!+
!  Set EOS options for analysis
!+
!----------------------------------------------------------------
subroutine set_eos_options(analysis_to_perform)
 integer, intent(in) :: analysis_to_perform
 integer             :: ierr

 ieos = 2
 call prompt('Enter ieos:',ieos)
 select case(ieos)
 case(2,12)
    gamma = 5./3.
    call prompt('Enter gamma:',gamma,0.)
    if (ieos==12) then
       gmw = 0.618212823
       call prompt('Enter mean molecular weight for gas+rad EoS:',gmw,0.)
    endif
 case(10,20)
    gamma = 5./3.
    X_in = 0.69843
    Z_in = 0.01426
    call prompt('Enter hydrogen mass fraction:',X_in,0.,1.)
    call prompt('Enter metallicity:',Z_in,0.,1.)
    irecomb = 0
    if (ieos==20) call prompt('Using gas+rad+rec EoS. Enter irecomb:',irecomb,0,2)
 case default
    call fatal('analysis_common_envelope',"EOS type not supported")
 end select
 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('analysis_common_envelope',"Failed to initialise EOS")

end subroutine set_eos_options

end module analysis
