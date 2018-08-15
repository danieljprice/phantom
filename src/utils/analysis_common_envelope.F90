!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!   Analysis routine for common envelope simulations
!
!  REFERENCES: None
!
!  OWNER: Thomas Reichardt
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, energies, eos, eos_mesa, kernel, part,
!    physcon, prompting, ptmass, setbinary, sortutils, units
!+
!--------------------------------------------------------------------------

module analysis

 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass,poten,ihsoft,rhoh,nsinkproperties,maxvxyzu,maxptmass
 use units,        only:print_units,umass,utime,udist,unit_ergg,unit_density,unit_pressure,unit_velocity,unit_Bfield,unit_energ
 use physcon,      only:gg,pi,c,kb_on_mh
 use prompting,    only:prompt
 use centreofmass, only:get_centreofmass, reset_centreofmass
 use energies,     only:compute_energies,ekin,etherm,epot,etot
 use ptmass,       only:get_accel_sink_gas,get_accel_sink_sink
 use kernel,       only:kernel_softening,radkern,wkern,cnormk
 use eos,          only:equationofstate,ieos,init_eos,finish_eos,X_in,Z_in,get_spsound
 use eos_mesa,     only:get_eos_kappa_mesa,get_eos_pressure_temp_mesa, get_eos_various_mesa
 use setbinary,    only:Rochelobe_estimate,L1_point
 use sortutils,    only:set_r2func_origin,r2func_origin,indexxfunc
 implicit none
 character(len=20), parameter, public :: analysistype = 'common_envelope'
 integer                              :: analysis_to_perform
 integer                              :: dump_number = 0

 real                                 :: omega_corotate=0, init_radius, rho_surface
 logical, dimension(5)                :: switch = .false.


 public                               :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 !general variables
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: unitnum,i,ierr,ncols


 !case 5 variables
 real                         :: rhopart

 !case 7 variables
 character(len=17), dimension(:), allocatable :: columns

 !case 12 variables
 real                         :: etoti, ekini, einti, epoti, phii


 real, dimension(3)           :: com_xyz, com_vxyz
 real, dimension(3)           :: xyz_a, vxyz_a
 real, dimension(6,npart)     :: histogram_data
 real                         :: ang_vel


 real, dimension(npart)      :: pres_1, proint_1, peint_1, temp_1, troint_1, teint_1, entrop_1, abad_1, gamma1_1, gam_1

 !case 16 variables
 real                        :: thermodynamic_quantities(5,npart)
 real, dimension(npart)      :: radius_1, dens_1


 !chose analysis type
 if (dump_number==0) then

    print "(16(a,/))", &
            ' 1) Sink separation', &
            ' 2) Bound and unbound quantities', &
            ' 3) Energies', &
            ' 4) Profile from centre of mass', &
            ' 5) Roche-lobe utils', &
            ' 6) Star stabilisation suite', &
            ' 7) Simulation units and particle properties', &
            ' 8) Output .divv', &
            ' 9) EoS testing', &
            '10) Ion fraction profiles in time', &
            '11) New unbound particle profiles in time', &
            '12) Sink properties', &
            '13) MESA EoS compute total entropy and other average td quantities', &
            '14) MESA EoS save on file thermodynamical quantities for all particles', &
            '15) Gravitational drag on sinks', &
            '16) Miscellaneous'

    analysis_to_perform = 1

    call prompt('Choose analysis type ',analysis_to_perform,1,16)

 endif

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 call adjust_corotating_velocities(npart,particlemass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,omega_corotate,dump_number)

 if ( ANY((/ 2, 3, 4, 6, 8, 9, 10, 11, 13, 14 /) == analysis_to_perform) .and. dump_number == 0 ) call init_eos(ieos,ierr)

 !analysis
 select case(analysis_to_perform)

 case(1) !sink separation
    call separation_vs_time(time,num)

 case(2) !bound and unbound quantities
    call bound_mass(time, num, npart, particlemass, xyzh, vxyzu)

 case(3) !Energies and bound mass
    call calculate_energies(time, num, npart, particlemass, xyzh, vxyzu)

 case(4) !Profile from COM (can be used for stellar profile)
    call create_profile(time, num, npart, particlemass, xyzh, vxyzu)

 case(5) !Mass within roche lobes
    call roche_lobe_values(time, num, npart, particlemass, xyzh, vxyzu)

 case(6) !Star stabilisation suite
    call star_stabilisation_suite(time, num, npart, particlemass, xyzh, vxyzu)

 case(7) !Units
    call print_simulation_parameters(num, npart, particlemass)

 case(8) !Output .divv
    call output_divv_files(time, dumpfile, num, npart, particlemass, xyzh, vxyzu)

 case(9) !EoS testing
    call eos_surfaces

 case(10) !Ion fraction profiles in time
    call ion_profiles(time, num, npart, particlemass, xyzh, vxyzu)

 case(11) !New unbound particle profiles in time
    !If you want to use this, remove the "/ real(n)" in the histogram_setup routine
    !This case can be somewhat tidied up
    call unbound_profiles(time, num, npart, particlemass, xyzh, vxyzu)

 case(12) !sink properties
    call sink_properties(time, num, npart, particlemass, xyzh, vxyzu)

 case(13) !MESA EoS compute total entropy and other average thermodynamical quantities
    call bound_unbound_thermo(time, num, npart, particlemass, xyzh, vxyzu)

 case(14) !MESA EoS save on file thermodynamical quantities for all particles

    do i=1,npart

       !particle radius
       radius_1(i) = distance(xyzh(1:3,i)) * udist

       !particles density in code units
       rhopart = rhoh(xyzh(4,i), particlemass)
       dens_1(i) = rhopart * unit_density

       !gets entropy for the current particle
       call get_eos_various_mesa(rhopart*unit_density,vxyzu(4,i) * unit_ergg, &
                                    pres_1(i),proint_1(i),peint_1(i),temp_1(i),troint_1(i), &
                                    teint_1(i),entrop_1(i),abad_1(i),gamma1_1(i),gam_1(i))

       !stores everything in an array
       thermodynamic_quantities(1,i) = radius_1(i)
       thermodynamic_quantities(2,i) = dens_1(i)
       thermodynamic_quantities(3,i) = pres_1(i)
       thermodynamic_quantities(4,i) = temp_1(i)
       thermodynamic_quantities(5,i) = entrop_1(i)

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

 case(15) !Gravitational drag on sinks
    call gravitational_drag(time, num, npart, particlemass, xyzh, vxyzu)

 case(16)
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
       xyz_a(1:3) = xyzmh_ptmass(1:3,i) - com_xyz(1:3)
       vxyz_a(1:3) = vxyz_ptmass(1:3,i) - com_vxyz(1:3)
       ang_vel = ang_vel + (-xyz_a(2) * vxyz_a(1) + xyz_a(1) * vxyz_a(2)) / dot_product(xyz_a(1:2), xyz_a(1:2))
    enddo

    ang_vel = ang_vel / 2.

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
 end select
 !increase dump number counter
 dump_number = dump_number + 1

end subroutine do_analysis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                Analysis  routines                !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! Separation vs time !!!!!
subroutine separation_vs_time(time,num)
 integer, intent(in)            :: num
 real, intent(in)               :: time
 character(len=17), allocatable :: columns(:)
 real                           :: sink_separation(4,nptmass-1)
 integer                        :: i, ncols
 ncols = 4 * (nptmass-1)
 allocate(columns(ncols))

 do i=1,(nptmass-1)
    call separation_vector(xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,i+1),sink_separation(1:4,i))

    write(columns((i*4)-3), '(A11,I1)') '    x sep. ', i
    write(columns((i*4)-2), '(A11,I1)') '    y sep. ', i
    write(columns((i*4)-1), '(A11,I1)') '    z sep. ', i
    write(columns((i*4)), '(A11,I1)')   '      sep. ', i
 enddo

 call write_time_file('separation_vs_time', columns, time, sink_separation, ncols, dump_number)
 deallocate(columns)
end subroutine separation_vs_time


!!!!! Bound mass !!!!!
subroutine bound_mass(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)            :: npart, num
 real, intent(in)               :: time, particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 real                           :: etoti, ekini, einti, epoti, phii
 real                           :: rhopart, ponrhoi, spsoundi
 real, dimension(3)             :: rcrossmv
 real, dimension(24)            :: bound
 integer                        :: i, bound_i, ncols
 character(len=17), allocatable :: columns(:)

 bound = 0.

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 do i=1,npart

    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

    rhopart = rhoh(xyzh(4,i), particlemass)

    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))

    call cross(xyzh(1:3,i), particlemass * vxyzu(1:3,i), rcrossmv)

    !bound criterion
    if (epoti + ekini < 0.0) then
       bound_i = 1
    else
       bound_i = 5
    endif

    bound(bound_i)     = bound(bound_i)     + 1
    bound(bound_i + 1) = bound(bound_i + 1) + particlemass
    bound(bound_i + 2) = bound(bound_i + 2) + distance(rcrossmv)
    bound(bound_i + 3) = bound(bound_i + 3) + etoti

    !bound criterion INCLUDING thermal energy
    if (etoti < 0.0) then
       bound_i = 9
    else
       bound_i = 13
    endif

    bound(bound_i)     = bound(bound_i)     + 1
    bound(bound_i + 1) = bound(bound_i + 1) + particlemass
    bound(bound_i + 2) = bound(bound_i + 2) + distance(rcrossmv)
    bound(bound_i + 3) = bound(bound_i + 3) + etoti

    !bound criterion INCLUDING enthalpy
    if (etoti + ponrhoi*particlemass < 0.0) then
       bound_i = 17
    else
       bound_i = 21
    endif

    bound(bound_i)     = bound(bound_i)     + 1
    bound(bound_i + 1) = bound(bound_i + 1) + particlemass
    bound(bound_i + 2) = bound(bound_i + 2) + distance(rcrossmv)
    bound(bound_i + 3) = bound(bound_i + 3) + etoti

 enddo

 ncols = 24
 allocate(columns(ncols))
 columns = (/'  b num part', &
             '      b mass', &
             '   b ang mom', &
             '    b tot en', &
             ' ub num part', &
             '     ub mass', &
             '  ub ang mom', &
             '   ub tot en', &
             ' bt num part', &
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
             '  ube tot en'/)

 call write_time_file('boundunbound_vs_time', columns, time, bound, ncols, dump_number)
 deallocate(columns)
end subroutine bound_mass


!!!!! Calculate energies !!!!!
subroutine calculate_energies(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)            :: npart, num
 real, intent(in)               :: time, particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 real                           :: etoti, ekini, einti, epoti, phii, phii1, jz, fxi, fyi, fzi
 real                           :: rhopart, ponrhoi, spsoundi, r_ij, radvel
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

 encomp(5:) = 0.

 call compute_energies(time)

 ekin = 0.

 do i=1,npart
    encomp(ipot_pp) = encomp(ipot_pp) + poten(i)
    encomp(ipot_env) = encomp(ipot_env) + poten(i)

    call cross(xyzh(1:3,i),particlemass * vxyzu(1:3,i),rcrossmv)

    jz = rcrossmv(3)
    encomp(ijz_tot) = encomp(ijz_tot) + jz

    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

    encomp(ipot_ps) = encomp(ipot_ps) + particlemass * phii

    phii1 = 0.
    call get_accel_sink_gas(1,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,fxi,fyi,fzi,phii1)
    encomp(ipot_env) = encomp(ipot_env) + phii1 * particlemass

    do j=1,nptmass
       r_ij = separation(xyzmh_ptmass(1:3,j),xyzh(1:3,i))

       if (r_ij < 80.) then
          inearsink = .true.
       endif
    enddo

    rhopart = rhoh(xyzh(4,i), particlemass)
    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))

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

    call cross(xyzmh_ptmass(1:3,i), xyzmh_ptmass(4,i)*vxyz_ptmass(1:3,i), rcrossmv)

    jz = rcrossmv(3)
    encomp(ijz_tot) = jz + encomp(ijz_tot)
    encomp(ijz_orb) = jz + encomp(ijz_orb)
    encomp(ikin_sink) = encomp(ikin_sink) + 0.5 * xyzmh_ptmass(4,i) * distance(vxyz_ptmass(1:3,i))**2
    if (i==2) encomp(iorb_comp) = encomp(iorb_comp) + 0.5 * xyzmh_ptmass(4,i) * distance(vxyz_ptmass(1:3,i))**2
 enddo

 do i=1,nptmass-1
    do j=i+1,nptmass
       r_ij = separation(xyzmh_ptmass(1:3,i),xyzmh_ptmass(1:3,j))

       encomp(ipot_sink) = encomp(ipot_sink) - xyzmh_ptmass(4,i) * xyzmh_ptmass(4,j) / r_ij
       if (i==1 .and. j==2) encomp(iorb_comp) = encomp(iorb_comp) - xyzmh_ptmass(4,i) * xyzmh_ptmass(4,j) / r_ij
    enddo
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

 deallocate(profile)
 deallocate(columns)
end subroutine create_profile


!!!!! Roche lobe values !!!!!
subroutine roche_lobe_values(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)            :: npart, num
 real, intent(in)               :: time, particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 integer                        :: i, nFB, nR1T, ncols
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
 real                           :: temp_const, ponrhoi, spsoundi
 real, dimension(3)             :: rcrossmv, CoO, com_xyz, com_vxyz

 MRL = 0.
 rhovol = 0.
 rhomass = 0.
 nFB = 0
 nR1T = 0
 temp_const = (unit_pressure / unit_density) * 1.34 / kb_on_mh

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

 do i=1,npart
    rhopart = rhoh(xyzh(4,i), particlemass)
    if (rhopart > rho_surface) then
       if (separation(xyzh(1:3,i), xyzmh_ptmass(1:3,1)) / 2. < &
              separation(xyzmh_ptmass(1:3,2), xyzmh_ptmass(1:3,1))) then
          rhomass = rhomass + particlemass
          rhovol = rhovol + particlemass / rhopart
       endif
    endif
 enddo

 sep = separation(xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,2))
 MRL(iRL1) = Rochelobe_estimate(m2,m1,sep)
 MRL(iRL2) = Rochelobe_estimate(m1,m2,sep)

 R1 = (3. * rhovol/(4. * pi))**(1./3.)
 CoO(1:3) = (xyzmh_ptmass(1:3,1) + xyzmh_ptmass(1:3,2)) / 2.
 MRL(iR1) = R1
 MRL(iRej) = separation(CoO(1:3),xyzmh_ptmass(1:3,1)) + R1

 call orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)

 do i=1,npart
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

    sep1 = separation(xyzmh_ptmass(1:3,1),xyzh(1:3,i))
    sep2 = separation(xyzmh_ptmass(1:3,2),xyzh(1:3,i))
    sepCoO = separation(CoO(1:3),xyzh(1:3,i))

    call cross(xyzh(1:3,i) - com_xyz(1:3),particlemass * vxyzu(1:3,i),rcrossmv)
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
       call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))
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

 call cross(xyzmh_ptmass(1:3,1) - com_xyz(1:3),xyzmh_ptmass(4,1) * vxyz_ptmass(1:3,1),rcrossmv)
 MRL(ijzRL1) = MRL(ijzRL1) + rcrossmv(3)

 call cross(xyzmh_ptmass(1:3,2) - com_xyz(1:3),xyzmh_ptmass(4,2) * vxyz_ptmass(1:3,2),rcrossmv)
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
 deallocate(columns)
end subroutine roche_lobe_values

!!!!! Star stabilisation suite !!!!!
subroutine star_stabilisation_suite(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)            :: npart, num
 real, intent(in)               :: time, particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 integer                        :: i, ncols
 real                           :: star_stability(4)
 real                           :: total_mass, rhovol, totvol, rhopart
 integer, parameter             :: ivoleqrad    = 1
 integer, parameter             :: idensrad     = 2
 integer, parameter             :: imassout     = 3
 integer, parameter             :: imassfracout = 4

 totvol = 0
 rhovol = 0

 if (dump_number == 0) then
    rho_surface = rhoh(xyzh(4,1), particlemass)
    do i=1,npart
       rhopart = rhoh(xyzh(4,i), particlemass)
       if (rhopart < rho_surface) then
          rho_surface = rhopart
       endif
    enddo
 endif

 do i=1,npart
    rhopart = rhoh(xyzh(4,i), particlemass)
    totvol = totvol + particlemass / rhopart
    if (rhopart > rho_surface) then
       if (nptmass > 1) then
          if (separation(xyzh(1:3,i), xyzmh_ptmass(1:3,1)) / 2. < &
              separation(xyzmh_ptmass(1:3,2), xyzmh_ptmass(1:3,1))) then
             rhovol = rhovol + particlemass / rhopart
          endif
       else
          rhovol = rhovol + particlemass / rhopart
       endif
    endif
 enddo

 star_stability(ivoleqrad) = (3. * totvol/(4. * pi))**(1./3.)
 star_stability(idensrad)  = (3. * rhovol/(4. * pi))**(1./3.)

 if (dump_number == 0) then
    init_radius = star_stability(ivoleqrad)
 endif

 star_stability(imassout) = 0.
 total_mass = xyzmh_ptmass(4,1)
 do i=1,npart
    if (separation(xyzmh_ptmass(1:3,1),xyzh(1:3,i)) > init_radius) then
       star_stability(imassout) = star_stability(imassout) + particlemass
    endif
    total_mass = total_mass + particlemass
 enddo

 star_stability(imassfracout) = star_stability(imassout) / total_mass

 ncols = 4
 allocate(columns(ncols))
 columns = (/'vol. eq. rad',&
             ' density rad',&
             'mass outside',&
             'frac outside'/)

 call write_time_file('star_stability', columns, time, star_stability, ncols, dump_number)
 deallocate(columns)
end subroutine star_stabilisation_suite

!!!!! Print simulation parameters !!!!!
subroutine print_simulation_parameters(num, npart, particlemass)
 integer, intent(in)            :: npart, num
 real, intent(in)               :: particlemass
 integer                        :: i

 write(*,"(/,3(a,es10.3,1x),a)") '     Mass: ',umass,    'g       Length: ',udist,  'cm     Time: ',utime,'s'
 write(*,"(3(a,es10.3,1x),a)") '  Density: ',unit_density, 'g/cm^3  Energy: ',unit_energ,'erg    En/m: ',unit_ergg,'erg/g'
 write(*,"(3(a,es10.3,1x),a)") ' Velocity: ',unit_velocity,'cm/s    Bfield: ',unit_Bfield,'G  Pressure: ',&
                                     unit_pressure,'g/cm s^2'
 write(*,"(2(a,es10.3,1x),/)")   '        G: ', gg*umass*utime**2/udist**3,'             c: ',c*utime/udist

 do i=1,nptmass
    write(*,'(A,I2,A,ES10.3,A,ES10.3)') 'Point mass ',i,': M = ',xyzmh_ptmass(4,i),' and h_soft = ',xyzmh_ptmass(ihsoft,i)
 enddo

 write(*,'(A,I7,A,ES10.3)') 'Gas particles : ',npart,' particles, each of mass ',particlemass

end subroutine print_simulation_parameters


!!!!! Output .divv files !!!!!
subroutine output_divv_files(time, dumpfile, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)            :: npart, num
 character(len=*), intent(in)   :: dumpfile
 real, intent(in)               :: time, particlemass
 real, intent(inout)            :: xyzh(:,:),vxyzu(:,:)
 integer                        :: i
 real                           :: etoti, ekini, einti, epoti, phii
 real                           :: rhopart, ponrhoi, spsoundi
 real, dimension(npart)         :: bound_energy, mach_number, omega_ratio
 real, dimension(3)             :: xyz_a, vxyz_a, com_xyz, com_vxyz
 real                           :: ang_vel

 call compute_energies(time)

 call orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)

 ang_vel = 0.

 do i=1,nptmass
    xyz_a(1:3) = xyzmh_ptmass(1:3,i) - com_xyz(1:3)
    vxyz_a(1:3) = vxyz_ptmass(1:3,i) - com_vxyz(1:3)
    ang_vel = ang_vel + 0.5 * (-xyz_a(2) * vxyz_a(1) + xyz_a(1) * vxyz_a(2)) / dot_product(xyz_a(1:2), xyz_a(1:2))
 enddo

 do i=1,npart
    call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

    bound_energy(i) = etoti/particlemass
    rhopart = rhoh(xyzh(4,i), particlemass)

    call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))
    mach_number(i) = distance(vxyzu(1:3,i)) / spsoundi

    xyz_a(1:3) = xyzh(1:3,i) - com_xyz(1:3)
    vxyz_a(1:3) = vxyzu(1:3,i) - com_vxyz(1:3)

    omega_ratio(i) = (-xyz_a(2) * vxyz_a(1) + xyz_a(1) * vxyz_a(2)) / dot_product(xyz_a(1:2), xyz_a(1:2))
    omega_ratio(i) = (omega_ratio(i) - ang_vel) / ang_vel

    !pressure(i) = ponrhoi * rhopart * unit_pressure


    !if (ieos==10) then
    !   call get_eos_pressure_temp_mesa(rhopart*unit_density,vxyzu(4,i) * unit_ergg,pressure(i),temp(i))
    !   call get_eos_kappa_mesa(rhopart*unit_density,temp(i),kappa(i),kappat,kappar)
       !call ionisation_fraction(rhopart * unit_density, temp(i), X_in, 1.-X_in-Z_in, &
       !                         ions(1,i),ions(2,i),ions(3,i),ions(4,i),ions(5,i))
    !else
    !   kappa(i) = 1
    !endif
 enddo

 open(26,file=trim(dumpfile)//".divv",status='replace',form='unformatted')
 write(26) (real(bound_energy(i),kind=4),i=1,npart)
 !write(26) (real(pressure(i),kind=4),i=1,npart)
 write(26) (real(omega_ratio(i),kind=4),i=1,npart)
 write(26) (real(mach_number(i),kind=4),i=1,npart)

 close(26)

end subroutine output_divv_files


!!!!! EoS surfaces !!!!!
subroutine eos_surfaces
 integer :: i, j
 real    :: rho_array(1000) = (/ (10**(i/10000.), i=-180000,-30150,150) /)
 real    :: eni_array(1000) = (/ (10**(i/10000.), i=120000,149970,30) /)
 real    :: temp_array(400) = (/ (10**(i/1000.), i=3000,6990,10) /)
 real    :: kappa_array(1000,400)
 real    :: pres_array(1000,1000)
 real    :: kappat, kappar, temp


 do i=1,size(rho_array)
    do j=1,size(eni_array)
       if (j < size(temp_array) + 1) then
          call get_eos_kappa_mesa(rho_array(i),temp_array(j),kappa_array(i,j),kappat,kappar)
       endif
       call get_eos_pressure_temp_mesa(rho_array(i),eni_array(j),pres_array(i,j),temp)
       !pres_array(i,j) = eni_array(j)*rho_array(i)*0.66667 / pres_array(i,j)
    enddo
 enddo


 open(unit=1000, file='mesa_eos_pressure.out', status='replace')

 !Write data to file
 do i=1,1000
    write(1000,"(1000(3x,es18.11e2,1x))") pres_array(i,:)
 enddo

 close(unit=1000)


 open(unit=1001, file='mesa_eos_kappa.out', status='replace')

 !Write data to file
 do i=1,1000
    write(1001,"(400(3x,es18.11e2,1x))") kappa_array(i,:)
 enddo

 close(unit=1001)

end subroutine eos_surfaces


!!!!! Ion profiles !!!!!
subroutine ion_profiles(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)          :: npart, num
 real, intent(in)             :: time, particlemass
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 integer                      :: npart_hist, nbins
 real, dimension(5,npart)     :: dist_part, rad_part
 real, dimension(:), allocatable :: hist_var
 real                         :: xh0, xh1, xhe0, xhe1, xhe2
 real                         :: pressure, temperature, rhopart
 character(len=17), dimension(:), allocatable :: grid_file
 character(len=40)            :: data_formatter
 integer                      :: i, unitnum

 call compute_energies(time)

 npart_hist = 0
 nbins = 300
 rad_part = 0.
 dist_part = 0.

 do i=1,npart
    if (xyzh(4,i)  >=  0) then
       rhopart = rhoh(xyzh(4,i), particlemass)
       npart_hist = npart_hist + 1
       rad_part(1,npart_hist) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
       call get_eos_pressure_temp_mesa(rhopart*unit_density,vxyzu(4,i) * unit_ergg,pressure,temperature)
       call ionisation_fraction(rhopart * unit_density, temperature, X_in, 1.-X_in-Z_in,xh0,xh1,xhe0,xhe1,xhe2)
       dist_part(1,npart_hist) = xh0
       dist_part(2,npart_hist) = xh1
       dist_part(3,npart_hist) = xhe0
       dist_part(4,npart_hist) = xhe1
       dist_part(5,npart_hist) = xhe2
    endif
 enddo
 allocate(hist_var(nbins))
 grid_file = (/ '   grid_HIfrac.ev', &
                '  grid_HIIfrac.ev', &
                '  grid_HeIfrac.ev', &
                ' grid_HeIIfrac.ev', &
                'grid_HeIIIfrac.ev' /)


 do i=1,5
    call histogram_setup(rad_part(1,1:npart_hist),dist_part(i,1:npart_hist),hist_var,npart_hist,3.,0.5,nbins,.true.)

    write(data_formatter, "(a,I5,a)") "(", nbins, "(3x,es18.10e3,1x))"

    if (num == 0) then
       unitnum = 1000

       open(unit=unitnum, file=trim(adjustl(grid_file(i))), status='replace')
       write(unitnum, "(a)") '# Ion fraction - look at the name of the file'
       close(unit=unitnum)
    endif

    unitnum=1001+i

    open(unit=unitnum, file=trim(adjustl(grid_file(i))), position='append')

    write(unitnum,data_formatter) hist_var(:)

    close(unit=unitnum)
 enddo

end subroutine ion_profiles


!!!!! Unbound profiles !!!!!
subroutine unbound_profiles(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)          :: npart, num
 real, intent(in)             :: time, particlemass
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 integer, dimension(4)        :: npart_hist
 real, dimension(5,npart)     :: dist_part, rad_part
 real, dimension(:), allocatable :: hist_var
 real                         :: etoti, ekini, einti, epoti, phii
 character(len=17), dimension(:), allocatable :: grid_file
 logical, allocatable, save   :: prev_unbound(:,:), prev_bound(:,:)
 character(len=40)            :: data_formatter
 integer                      :: i, unitnum, nbins
 call compute_energies(time)

 npart_hist = 0
 nbins = 300
 rad_part = 0.
 dist_part = 0.

 if (dump_number == 0) then
    allocate(prev_bound(2,npart))
    allocate(prev_unbound(2,npart))
    prev_bound = .false.
    prev_unbound = .false.
 endif


 do i=1,npart
    if (xyzh(4,i)  >=  0) then

       call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

       if ((etoti > 0.) .and. (prev_unbound(1,i) .eqv. .false.)) then
          npart_hist(1) = npart_hist(1) + 1
          rad_part(1,npart_hist(1)) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          dist_part(1,npart_hist(1)) = 1.
          prev_unbound(1,i) = .true.
       elseif (etoti < 0.) then
          prev_unbound(1,i) = .false.
       endif


       if ((ekini + epoti > 0.) .and. (prev_unbound(2,i) .eqv. .false.)) then
          npart_hist(2) = npart_hist(2) + 1
          rad_part(2,npart_hist(2)) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          dist_part(2,npart_hist(2)) = 1.
          prev_unbound(2,i) = .true.
       elseif (ekini + epoti < 0.) then
          prev_unbound(2,i) = .false.
       endif

       if ((etoti < 0.) .and. (prev_bound(1,i) .eqv. .false.)) then
          npart_hist(3) = npart_hist(3) + 1
          rad_part(3,npart_hist(3)) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          dist_part(3,npart_hist(3)) = 1.
          prev_bound(1,i) = .true.
       elseif (etoti > 0.) then
          prev_bound(1,i) = .false.
       endif

       if ((ekini + epoti < 0.) .and. (prev_bound(2,i) .eqv. .false.)) then
          npart_hist(4) = npart_hist(4) + 1
          rad_part(4,npart_hist(4)) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          dist_part(4,npart_hist(4)) = 1.
          prev_bound(2,i) = .true.
       elseif (ekini + epoti > 0.) then
          prev_bound(2,i) = .false.
       endif
    endif
 enddo
 allocate(hist_var(nbins))
 grid_file = (/ 'gridpkiunbound.ev', &
                ' gridpkunbound.ev', &
                '  gridpkibound.ev', &
                '   gridpkbound.ev' /)

 do i=1,4
    call histogram_setup(rad_part(i,1:npart_hist(i)),dist_part(i,1:npart_hist(i)),hist_var,npart_hist(i),3.,0.5,nbins,.false.)

    write(data_formatter, "(a,I5,a)") "(", nbins, "(3x,es18.10e3,1x))"

    if (num == 0) then
       unitnum = 1000

       open(unit=unitnum, file=trim(adjustl(grid_file(i))), status='replace')
       write(unitnum, "(a)") '# Newly bound/unbound particles'
       close(unit=unitnum)
    endif

    unitnum=1001+i

    open(unit=unitnum, file=trim(adjustl(grid_file(i))), position='append')

    write(unitnum,data_formatter) hist_var(:)

    close(unit=unitnum)
 enddo
end subroutine unbound_profiles


!!!!! Sink properties !!!!!
subroutine sink_properties(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)          :: npart, num
 real, intent(in)             :: time, particlemass
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 character(len=17)            :: filename
 real                         :: sinkcomp(30)
 real                         :: ang_mom(3)
 real                         :: phitot, dtsinksink, fonrmax
 real                         :: fxi, fyi, fzi, phii
 real, dimension(4,maxptmass) :: fssxyz_ptmass
 real, dimension(4,maxptmass) :: fxyz_ptmass
 integer                      :: i, ncols

 ncols = 25
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
             '    ang mom x', &
             '    ang mom y', &
             '    ang mom z', &
             '    |ang mom|', &
             '       kin en'/)

 fxyz_ptmass = 0.
 call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,phitot,dtsinksink,0,0.)
 fssxyz_ptmass = fxyz_ptmass
 do i=1,npart
    call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                            fxi,fyi,fzi,phii,particlemass,fxyz_ptmass,fonrmax)
 enddo

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
    call cross(xyzmh_ptmass(1:3,i), xyzmh_ptmass(4,i)*vxyz_ptmass(1:3,i), ang_mom)
    sinkcomp(21:23) = ang_mom
    ! angular momentum modulus
    sinkcomp(24)    = distance(ang_mom(1:3))
    ! kinetic energy
    sinkcomp(25)    = 0.5*xyzmh_ptmass(4,i)*sinkcomp(8)**2

    call write_time_file(filename, columns, time, sinkcomp, ncols, dump_number)
 enddo
 deallocate(columns)
end subroutine sink_properties

subroutine bound_unbound_thermo(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)          :: npart, num
 real, intent(in)             :: time, particlemass
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 integer                      :: i, ncols
 real, dimension(8)           :: entropy_array
 real                         :: etoti, ekini, einti, epoti, phii, rhopart
 real, dimension(npart)       :: pres_1, proint_1, peint_1, temp_1, troint_1, teint_1, entrop_1, abad_1, gamma1_1, gam_1
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
                                    pres_1(i),proint_1(i),peint_1(i),temp_1(i),troint_1(i), &
                                    teint_1(i),entrop_1(i),abad_1(i),gamma1_1(i),gam_1(i))

          !sums entropy and other quantities for bound particles and unbound particles

          if (.not. switch(1)) then
             etoti = etoti - einti
          endif

          if (etoti < 0.0) then !bound
             entropy_array(ient_b)  = entropy_array(ient_b) + entrop_1(i)
             entropy_array(itemp_b) = entropy_array(itemp_b) + temp_1(i)
             entropy_array(ipres_b) = entropy_array(ipres_b) + pres_1(i)
             entropy_array(idens_b) = entropy_array(idens_b) + rhopart*unit_density

          else !unbound
             entropy_array(ient_ub)  = entropy_array(ient_ub) + entrop_1(i)
             entropy_array(itemp_ub) = entropy_array(itemp_ub) + temp_1(i)
             entropy_array(ipres_ub) = entropy_array(ipres_ub) + pres_1(i)
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

!!!!! Gravitational drag !!!!!
subroutine gravitational_drag(time, num, npart, particlemass, xyzh, vxyzu)
 integer, intent(in)          :: npart, num
 real, intent(in)             :: time, particlemass
 real, intent(inout)          :: xyzh(:,:),vxyzu(:,:)
 character(len=17), allocatable :: columns(:)
 character(len=17)            :: filename
 integer                      :: nvc, i, j, k, iorder(npart), ncols
 real, dimension(:), allocatable, save :: ang_mom_old
 real, save                   :: time_old
 real, dimension(11,nptmass)   :: drag_force
 real, dimension(3)           :: avg_vel, avg_vel_par, avg_vel_per, com_xyz, com_vxyz, unit_vel, unit_vel_per, ang_mom
 real                         :: vel_contrast, sep, centre_sep
 real                         :: rhopart, cs=0., racc, fonrmax, fxi, fyi, fzi, phii
 real, dimension(4,maxptmass) :: fxyz_ptmass


    if (dump_number == 0) then
       allocate(ang_mom_old(nptmass))
       do i=1,nptmass
          call cross(xyzmh_ptmass(1:3,i), xyzmh_ptmass(4,i)*vxyz_ptmass(1:3,i), ang_mom)
          ang_mom_old(i) = ang_mom(3)
       enddo
       time_old = -50.
    endif

    drag_force = 0.

    ncols = 11
    allocate(columns(ncols))
    columns = (/'   par. num.', &
                '  perp. num.', &
                '    from dJz', &
                '  analytical', &
                'vel contrast', &
                'par. v. con.', &
                'per. v. con.', &
                ' sound speed', &
                ' rho at sink', &
                '        racc', &
                'com-sink sep'/)

    call orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)

    do i=2,2
       fxyz_ptmass = 0.
       avg_vel = 0
       avg_vel_par = 0
       avg_vel_per = 0
       vel_contrast = 0.
       nvc = 0
       racc = 0.
       cs = 0.
       rhopart = 0.

       call unit_vector(vxyz_ptmass(1:3,i), unit_vel(1:3))

       call set_r2func_origin(xyzmh_ptmass(1,i),xyzmh_ptmass(2,i),xyzmh_ptmass(3,i))
       call indexxfunc(npart,r2func_origin,xyzh,iorder)

       !k = iorder(1)
       !cs = get_spsound(ieos,xyzh(1:3,k),rhoh(xyzh(4,k), particlemass),vxyzu(1:3,k))
       !rhopart = rhoh(xyzh(4,k), particlemass)

       centre_sep = separation(com_xyz(1:3),xyzmh_ptmass(1:3,i))

       do j=1,npart
          k = iorder(j)
          sep = separation(xyzh(1:3,k),xyzmh_ptmass(1:3,i))
          if (sep > centre_sep)exit
          avg_vel(1:3) = avg_vel(1:3) + vxyzu(1:3,k)
          !vel_contrast = vel_contrast + distance(vxyz_ptmass(1:3,i)) - dot_product(vxyzu(1:3,k), unit_vel)
          cs = cs + get_spsound(ieos,xyzh(1:3,k),rhoh(xyzh(4,k), particlemass),vxyzu(1:3,k))
          !rhopart = rhopart + rhoh(xyzh(4,k), particlemass)
          rhopart = rhopart + particlemass
       enddo

       if (j > 1) then
          avg_vel(1:3) = avg_vel(1:3) / (j-1)
          avg_vel_par(1:3) = dot_product(avg_vel, unit_vel) * unit_vel
          avg_vel_per(1:3) = avg_vel(1:3) - avg_vel_par(1:3)
          vel_contrast = distance(vxyz_ptmass(1:3,i)) - cos_vector_angle(unit_vel, avg_vel_par) * distance(avg_vel_par(1:3))
          !vel_contrast = separation(vxyz_ptmass(1:3,i),avg_vel_par(1:3))
          cs = cs / float(j-1)
          rhopart = rhopart / float(j-1)
          racc = 2. * xyzmh_ptmass(4,i) / (vel_contrast**2 + cs**2)
       else
          racc = 0.
       endif

       do j=1,npart
          k = iorder(j)
          if (separation(xyzh(1:3,k),xyzmh_ptmass(1:3,i)) > centre_sep) exit
          call get_accel_sink_gas(nptmass,xyzh(1,k),xyzh(2,k),xyzh(3,k),xyzh(4,k),xyzmh_ptmass,&
                                  fxi,fyi,fzi,phii,particlemass,fxyz_ptmass,fonrmax)
       enddo

       call cross(unit_vel, (/ 0., 0., 1. /), unit_vel_per)

       drag_force(1,i) = dot_product(fxyz_ptmass(1:3,i),unit_vel)
       drag_force(2,i) = dot_product(fxyz_ptmass(1:3,i),unit_vel_per)
       drag_force(4,i) = - rhopart * (vel_contrast * abs(vel_contrast)) * pi * racc**2
       drag_force(5,i) = vel_contrast
       drag_force(6,i) = cos_vector_angle(unit_vel, avg_vel_par) * distance(avg_vel_par)
       drag_force(7,i) = distance(avg_vel_per)
       drag_force(8,i) = cs
       drag_force(9,i) = rhopart
       drag_force(10,i) = racc
       drag_force(11,i) = centre_sep
    enddo

    do i=2,2
       call cross(xyzmh_ptmass(1:3,i) - com_xyz(1:3), xyzmh_ptmass(4,i)*vxyz_ptmass(1:3,i), ang_mom)
       drag_force(3,i) = (ang_mom(3) - ang_mom_old(i)) / ((time - time_old) * distance(xyzmh_ptmass(1:3,i) - com_xyz(1:3)))
       ang_mom_old(i) = ang_mom(3)

       write (filename, "(A16,I0)") "sink_drag_", i

       call write_time_file(trim(adjustl(filename)), columns, time, drag_force(:,i), ncols, dump_number)
    enddo
    time_old = time
    deallocate(columns)
end subroutine gravitational_drag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!        Routines used in analysis routines        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_gas_energies(particlemass,poten,xyzh,vxyzu,xyzmh_ptmass,phii,epoti,ekini,einti,etoti)
 !calculates kinetic, potential and internal energy of a gas particle
 use ptmass, only:get_accel_sink_gas
 use part,   only:nptmass
 real, intent(in)                       :: particlemass
 real(4), intent(in)                    :: poten
 real, dimension(4), intent(in)         :: xyzh, vxyzu
 real, dimension(5,nptmass), intent(in) :: xyzmh_ptmass
 real, intent(out)                      :: phii, epoti, ekini, einti, etoti
 real                                   :: fxi, fyi, fzi

 phii = 0.0

 call get_accel_sink_gas(nptmass,xyzh(1),xyzh(2),xyzh(3),xyzh(4),xyzmh_ptmass,fxi,fyi,fzi,phii)

 epoti = poten + particlemass * phii
 ekini = particlemass * 0.5 * distance(vxyzu(1:3))**2
 einti = particlemass * vxyzu(4)
 etoti = epoti + ekini + einti
end subroutine calc_gas_energies


subroutine adjust_corotating_velocities(npart,particlemass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,omega_c,dump_number)
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
          call cross(omega_vec,xyzh(1:3,i),omegacrossr)
          vxyzu(1:3,i) = vxyzu(1:3,i) + omegacrossr(1:3)
       enddo

       do i=1,nptmass
          call cross(omega_vec,xyzmh_ptmass(1:3,i),omegacrossr)
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
 real              :: rhopart,ponrhoi,spsoundi
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
             call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))

             if (ieos == 10) then
                call get_eos_pressure_temp_mesa(rhopart*unit_density,vxyzu(4,i) * unit_ergg,pres,temp)
                call get_eos_kappa_mesa(rhopart*unit_density,temp,kappa,kappat,kappar)
             else
                temp = (ponrhoi * (unit_pressure/unit_density) * 2.381 * mass_proton_cgs) / kboltz
                kappa = 1.
             endif

             call calc_gas_energies(particlemass,poten(i),xyzh(:,i),vxyzu(:,i),xyzmh_ptmass,phii,epoti,ekini,einti,etoti)

             call ionisation_fraction(rhopart * unit_density, temp, X_in, 1.-X_in-Z_in,xh0, xh1, xhe0, xhe1, xhe2)

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

 print*, "Profile completed"

end subroutine stellar_profile

subroutine orbit_com(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,com_xyz,com_vxyz)
 integer, intent(in) :: npart, nptmass
 real, intent(in) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 real, intent(out), dimension(3) :: com_xyz, com_vxyz


 real, dimension(4,npart) :: xyz_a
 real, dimension(3,npart) :: vxyz_a
 integer                  :: iorder(npart)
 real                     :: sep
 integer                  :: i,j


 com_xyz(1) = sum(xyzmh_ptmass(1,:))/nptmass
 com_xyz(2) = sum(xyzmh_ptmass(2,:))/nptmass
 com_xyz(3) = sum(xyzmh_ptmass(3,:))/nptmass

 call set_r2func_origin(com_xyz(1),com_xyz(2),com_xyz(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
    
 sep = separation(xyzmh_ptmass(1:3,1),com_xyz(1:3))

 do i=1,npart
    j = iorder(i)
    if (separation(xyzh(1:3,j),com_xyz(1:3)) > 2.*sep) exit
    xyz_a(1:4,i) = xyzh(1:4,j)
    vxyz_a(1:3,i) = vxyzu(1:3,j)
 enddo

 call get_centreofmass(com_xyz,com_vxyz,i-1,xyz_a,vxyz_a,nptmass,xyzmh_ptmass,vxyz_ptmass)

end subroutine orbit_com

subroutine ionisation_fraction(dens,temp,X,Y,xh0, xh1, xhe0, xhe1, xhe2)
 !solves three Saha equations simultaneously to return ion fractions of hydrogen and helium
 use physcon, only:twopi, kboltz, eV, planckh, mass_electron_cgs, mass_proton_cgs
 real, intent(in) :: dens, temp, X, Y
 real, intent(out):: xh0, xh1, xhe0, xhe1, xhe2
 real             :: n, nh, nhe
 real             :: A, B, C, const
 real             :: xh1g, xhe1g, xhe2g
 real             :: f, g, h
 real, parameter  :: chih0 = 13.6, chihe0 = 24.6, chihe1 = 54.4
 real, dimension(3,3) :: M, M_inv
 real, dimension(3) :: dx
 integer          :: i

 nh = X * dens / mass_proton_cgs
 nhe = Y * dens / (4. * mass_proton_cgs)
 n = nh + nhe

 const = (sqrt(twopi * mass_electron_cgs * kboltz) / planckh)**3 / n

 A = 1. * const * temp**(1.5) * exp(-chih0 * eV / (kboltz * temp))
 B = 4. * const * temp**(1.5) * exp(-chihe0 * eV / (kboltz * temp))
 C = 1. * const * temp**(1.5) * exp(-chihe1 * eV / (kboltz * temp))

 xh1g = 0.4
 xhe1g = 0.3
 xhe2g = 0.2

 do i=1,50
    f = xh1g * (xh1g + xhe1g + 2*xhe2g) - A * ((nh/n) - xh1g)
    g = xhe1g * (xh1g + xhe1g + 2*xhe2g) - B * ((nhe/n) - xhe1g - xhe2g)
    h = xhe2g * (xh1g + xhe1g + 2*xhe2g) - C * xhe1g

    M(1,:) = (/ 2*xh1g + xhe1g + 2*xhe2g + A, xh1g, 2*xh1g /)
    M(2,:) = (/ xhe1g, xh1g + 2*xhe1g + 2*xhe2g + B, 2*xhe1g + B /)
    M(3,:) = (/ xhe2g, xhe2g - C, xh1g + xhe1g + 4*xhe2g /)

    call minv(M, M_inv)

    dx = matmul(M_inv, (/ -f, -g, -h/))

    xh1g = xh1g + dx(1)
    xhe1g = xhe1g + dx(2)
    xhe2g = xhe2g + dx(3)
 enddo

 xh1 = xh1g * n / nh
 xhe1 = xhe1g * n / nhe
 xhe2 = xhe2g * n / nhe
 xh0 = ((nh/n) - xh1g) * n / nh
 xhe0 = ((nhe/n) - xhe1g - xhe2g) * n / nhe
end subroutine ionisation_fraction

subroutine histogram_setup(dist_var,avg_var,hist_var,npart,max_value,min_value,nbins,average)
 !returns a radial histogram of a given distribution variable
 integer, intent(in)                :: npart,nbins
 real, dimension(npart), intent(in) :: dist_var, avg_var
 real, dimension(nbins), intent(out) :: hist_var
 real, intent(in)                   :: max_value, min_value
 logical, intent(in)                :: average
 real                               :: bins(nbins)
 integer                            :: i,j,n

 bins = (/ (10**(min_value + (i-1) * (max_value-min_value)/real(nbins)), i=1,nbins) /)
 !bins = (/ (min_value + (i-1) * (max_value-min_value)/real(nbins), i=1,nbins) /)
 hist_var(:) = 0.

 do i=1,nbins-1
    n = 0
    do j=1,npart
       if (dist_var(j) >= bins(i) .and. dist_var(j) < bins(i+1)) then
          n = n + 1
          hist_var(i) = hist_var(i) + avg_var(j)
       endif
    enddo
    if (n>0) then
       if (average) then
          hist_var(i) = hist_var(i) / real(n)
       endif
    endif
 enddo

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

 write(column_formatter, "(a,I2,a)") "('#',3x,", ncols, "('[',a15,']',5x))"
 write(data_formatter, "(a,I2,a)") "(", ncols, "(3x,es19.11e3,1x))"

 do i=1,ncols
    write(columns(i), "(I2,a)") i, cols(i)
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

 write(column_formatter, "(a,I2,a)") "('#',3x,", ncols+1, "('[',a15,']',5x))"
 write(data_formatter, "(a,I2,a)") "(", ncols+1, "(3x,es18.11e2,1x))"
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

end subroutine

subroutine cross(a,b,c)
 ! Return the vector cross product of two 3d vectors
 real,intent(in),dimension(3)  :: a,b
 real,intent(out),dimension(3) :: c

 c(1) = a(2)*b(3)-b(2)*a(3)
 c(2) = a(3)*b(1)-b(3)*a(1)
 c(3) = a(1)*b(2)-b(1)*a(2)

end subroutine cross

real function distance(a)
 !return modulus of an arbitrary vector
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
    cos_vector_angle = 1
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
 real, intent(in), dimension(3) :: a,b

 separation = distance(a(1:3) - b(1:3))
end function separation

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

subroutine minv (M, M_inv)

 implicit none

 real, dimension(3,3), intent(in)  :: M
 real, dimension(3,3), intent(out) :: M_inv

 real :: det
 real, dimension(3,3) :: cofactor


 det =   M(1,1)*M(2,2)*M(3,3)  &
       - M(1,1)*M(2,3)*M(3,2)  &
       - M(1,2)*M(2,1)*M(3,3)  &
       + M(1,2)*M(2,3)*M(3,1)  &
       + M(1,3)*M(2,1)*M(3,2)  &
       - M(1,3)*M(2,2)*M(3,1)

 cofactor(1,1) = +(M(2,2)*M(3,3)-M(2,3)*M(3,2))
 cofactor(1,2) = -(M(2,1)*M(3,3)-M(2,3)*M(3,1))
 cofactor(1,3) = +(M(2,1)*M(3,2)-M(2,2)*M(3,1))
 cofactor(2,1) = -(M(1,2)*M(3,3)-M(1,3)*M(3,2))
 cofactor(2,2) = +(M(1,1)*M(3,3)-M(1,3)*M(3,1))
 cofactor(2,3) = -(M(1,1)*M(3,2)-M(1,2)*M(3,1))
 cofactor(3,1) = +(M(1,2)*M(2,3)-M(1,3)*M(2,2))
 cofactor(3,2) = -(M(1,1)*M(2,3)-M(1,3)*M(2,1))
 cofactor(3,3) = +(M(1,1)*M(2,2)-M(1,2)*M(2,1))

 M_inv = transpose(cofactor) / det

 return

end subroutine minv


end module
