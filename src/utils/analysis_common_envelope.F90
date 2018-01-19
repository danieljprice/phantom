!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, energies, eos, eos_mesa, kernel, part,
!    physcon, prompting, ptmass, setbinary, units
!+
!--------------------------------------------------------------------------

module analysis


 implicit none
 character(len=20), parameter, public :: analysistype = 'common_envelope'
 integer                              :: analysis_to_perform, histogram_number
 integer                              :: dump_number = 0
 character(len=32)                    :: val1, val2
 real                                 :: constraint_max, constraint_min, omega_corotate=0, init_radius, surfacedensity
 integer                              :: bin_res
 logical, dimension(5)                :: switch = .false.

 real                                 :: hsoftk, hsoft1, hsoft21, q2i, qi, psoft, fsoft
 real, allocatable                    :: prev_bound_energy(:)
 logical, allocatable                 :: prev_unbound(:)


 public                               :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 !phantom modules
 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass,poten,ihsoft,rhoh
 use units,        only:print_units,umass,utime,udist,unit_ergg,unit_density,unit_pressure,unit_velocity,unit_Bfield,unit_energ
 use physcon,      only:gg,pi,c
 use prompting,    only:prompt
 use centreofmass, only:get_centreofmass
 use energies,     only:compute_energies,ekin,etherm,epot,etot
 use ptmass,       only:get_accel_sink_sink, get_accel_sink_gas
 use kernel,       only:kernel_softening,radkern,wkern,cnormk
 use eos,          only:gamma,equationofstate,ieos,init_eos,finish_eos,X_in,Z_in
 use eos_mesa,     only:get_eos_kappa_mesa,get_eos_pressure_temp_mesa, get_eos_various_mesa
 use setbinary,    only:Rochelobe_estimate

 !general variables
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character(len=30)            :: fileout
 integer                      :: unitnum,i,j,k,ierr
 real                         :: hsoft = 3.0

 !case 1 variables
 real                         :: sink_separation(4,nptmass-1)

 !case 2 variables
 real                         :: com_pos(3),com_v(3)
 real                         :: etot_current_part,etherm_current_part
 real                         :: r_sink_part, ppart
 real                         :: ltot_current_part(3)
 real                         :: bound(24)

 integer, parameter           :: in_b1    = 1
 integer, parameter           :: imass_b1 = 2
 integer, parameter           :: iltot_b1 = 3
 integer, parameter           :: ietot_b1 = 4
 integer, parameter           :: in_u1    = 5
 integer, parameter           :: imass_u1 = 6
 integer, parameter           :: iltot_u1 = 7
 integer, parameter           :: ietot_u1 = 8
 integer, parameter           :: in_b2    = 9
 integer, parameter           :: imass_b2 = 10
 integer, parameter           :: iltot_b2 = 11
 integer, parameter           :: ietot_b2 = 12
 integer, parameter           :: in_u2    = 13
 integer, parameter           :: imass_u2 = 14
 integer, parameter           :: iltot_u2 = 15
 integer, parameter           :: ietot_u2 = 16
 integer, parameter           :: in_b3    = 17
 integer, parameter           :: imass_b3 = 18
 integer, parameter           :: iltot_b3 = 19
 integer, parameter           :: ietot_b3 = 20
 integer, parameter           :: in_u3    = 21
 integer, parameter           :: imass_u3 = 22
 integer, parameter           :: iltot_u3 = 23
 integer, parameter           :: ietot_u3 = 24

 !case 3 variables
 real                         :: encomp(23)
 real                         :: sink_vel(4)=0., r_sink_sink, radvel
 real                         :: part_kin, part_ps_pot, part_therm
 real                         :: ponrhoi, spsoundi
 real                         :: e=0, a, mtot, sep, vmag, v_part_com(4)=0, omega_vec(3)=0, omegacrossr(3)
 real                         :: r_part_com(4), rcrossmv(3), jz, boundparts(8,npart+nptmass)
 logical                      :: inearsink

 integer, parameter           :: ie_tot        = 1
 integer, parameter           :: ie_pot        = 2
 integer, parameter           :: ie_kin        = 3
 integer, parameter           :: ie_therm      = 4
 integer, parameter           :: ipot_sink     = 5
 integer, parameter           :: ikin_sink     = 6
 integer, parameter           :: iorb_sink     = 7
 integer, parameter           :: iorb_comp     = 8
 integer, parameter           :: ipot_env      = 9
 integer, parameter           :: ie_env        = 10
 integer, parameter           :: ikin_bound    = 11
 integer, parameter           :: ikin_unbound  = 12
 integer, parameter           :: imass_bound   = 13
 integer, parameter           :: imass_unbound = 14
 integer, parameter           :: ipot_pp       = 15
 integer, parameter           :: ipot_ps       = 16
 integer, parameter           :: ijz_tot       = 17
 integer, parameter           :: ijz_bound     = 18
 integer, parameter           :: ijz_unbound   = 19
 integer, parameter           :: ijz_orb       = 20
 integer, parameter           :: ie_gas        = 21
 integer, parameter           :: fallbackmass  = 22
 integer, parameter           :: fallbackmom   = 23

 !case 4 variables
 real, allocatable            :: profile(:,:)
 real                         :: ray(3), xyz(3,npart), xyz_profile(3), proj(3), proj_orth(3)
 real                         :: proj_mag, proj_dist, proj_orth_dist, far_dist, Wj, kappa(npart), kappat, kappar
 real                         :: pres, temp(npart)
 integer                      :: far_i, n_sample
 real, dimension(:,:), allocatable    :: ray_profile
 real                         :: profile_vector(3)
 character(len=15)            :: name_in

 !case 5 variables
 real                         :: totvol, rhopart, VER(1)

 !case 6 variables
 !real                         :: m1,m2
 real                         :: RL1,RL2
 real, dimension(6)           :: MRL=0

 !case 7 variables
 real, dimension(4,npart)                     :: distance_from_com, velocity_wrt_com
 real                                         :: vel_rad_comp(4,npart), vel_tan_comp(4,npart), projection, ang_vel
 real, dimension(:,:), allocatable            :: histogram_bins
 real                                         :: bin_size
 character(len=17), dimension(:), allocatable :: columns

 !case 8 variables
 real                         :: star_stabilisation_params(4), total_mass=0.0, densvol

 integer, parameter           :: ivoleqrad    = 1
 integer, parameter           :: idensrad     = 2
 integer, parameter           :: imassout     = 3
 integer, parameter           :: imassfracout = 4

 !case 9 variables
 real                         :: C_force, dtsg, f2i, hi, dtsgi, fextx, fexty, fextz, phii, fonrmaxi, dtphi2i
 real                         :: fxyz_ptmass_thread(4,nptmass), minsep1, minsep2, minsep1i, minsep2i

 !case 12 variables
 real, dimension(npart)       :: bound_energy

 !case 13 variables
 real  :: rho_array(10000) = (/ (10**(i/10000.), i=-180000,-30015,15) /)
 real  :: temp_array(4000) = (/ (10**(i/1000.), i=3000,6999) /)
 real  :: kappa_array(10000,4000)

 !case 14 variables
 integer                      :: npart_hist, nbins
 real, dimension(5,npart)     :: dist_part
 real, dimension(npart)       ::rad_part
 real, dimension(:), allocatable :: hist_var
 character(len=40)            :: data_formatter
 real                         :: orth_ratio, temperature
 logical                      :: criteria
 real                         :: xh0, xh1, xhe0, xhe1, xhe2
 character(len=17), dimension(5) :: grid_file

 !case 15 variables
 real                        :: entropy_bound = 0.0, entropy_unbound = 0.0, entropy_array(6)
 real                        :: avgtemp_bound = 0.0, avgtemp_unbound = 0.0
 real                        :: avgpres_bound = 0.0, avgpres_unbound = 0.0
 real                        :: boundparts_1(8,npart)
 real, dimension(npart)      :: pres_1, temp_1, proint_1, peint_1, troint_1, teint_1, entrop_1, abad_1, gamma1_1, gam_1

 !chose analysis type
 if (dump_number==0) then
    print "(15(/,a,/))", ' 1) Sink separation', &
            ' 2) Bound and unbound quantities', &
            ' 3) Energies', &
            ' 4) Profile from centre of mass', &
            ' 5) Volume equivalent radius of star', &
            ' 6) Mass within initial Roche-lobes', &
            ' 7) Histograms ', &
            ' 8) Star stabilisation suite', &
            ' 9) Output simulation units', &
            '10) Something', &
            '11) Interpolated ray profile', &
            '12) Output .divv', &
            '13) Print sink particle m and h', &
            '14) Another something', &
            '15) Entropy MESA EoS'

    analysis_to_perform = 1

    call prompt('Choose analysis type ',analysis_to_perform,1,15)

 endif

 if ( ANY((/ 2, 3, 4, 8, 11, 12, 13, 14, 15 /) == analysis_to_perform) .and. dump_number == 0 ) call init_eos(ieos,ierr)

 !analysis
 select case(analysis_to_perform)

 case(1) !sink separation

    allocate(columns(4 * (nptmass-1)))
    !computes sink separation
    do i=1,(nptmass-1)
       call separation_vector(xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,i+1),sink_separation(1:4,i))

       write(columns((i*4)-3), '(A11,I1)') '    x sep. ', i
       write(columns((i*4)-2), '(A11,I1)') '    y sep. ', i
       write(columns((i*4)-1), '(A11,I1)') '    z sep. ', i
       write(columns((i*4)), '(A11,I1)')   '      sep. ', i
    enddo

    call write_time_file('separation_vs_time', columns, time, sink_separation, 4*(nptmass-1), num)
    deallocate(columns)

 case(2) !bound and unbound quantities
    bound = 0
    !get the position of the COM to compute all the vectorial quantities wrt it
    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    !loop on particles
    do i=1,npart

       !part to sinks pot en
       part_ps_pot = 0.0

       do k=1,nptmass

          r_sink_part = separation(xyzmh_ptmass(1:3,k),xyzh(1:3,i))

          part_ps_pot = part_ps_pot - ((gg*umass*utime**2/udist**3) * xyzmh_ptmass(4,k)) / r_sink_part

       enddo

       part_ps_pot = particlemass * part_ps_pot

       !kin en
       part_kin = particlemass * 0.5 * separation(vxyzu(1:3,i),com_v(1:3))**2

       !therm en
       etherm_current_part = particlemass * vxyzu(4,i)

       !tot en
       etot_current_part = poten(i) + part_ps_pot + part_kin + etherm_current_part

       rhopart = rhoh(xyzh(4,i), particlemass)
       ppart = (gamma - 1.) * rhopart * etherm_current_part

       call separation_vector(xyzh(1:3,i),com_pos(1:3),r_part_com(1:4))
       call separation_vector(vxyzu(1:3,i),com_v(1:3),v_part_com(1:4))

       call cross(r_part_com(1:3),v_part_com(1:3)*particlemass,ltot_current_part)

       !bound criterion NOT INCLUDING thermal energy
       if (poten(i) + part_ps_pot + part_kin < 0.0) then

          !num parts, total mass, ang mom and total energy respectively - bound
          bound(in_b1) = bound(in_b1) + 1
          bound(imass_b1) = bound(imass_b1) + particlemass
          bound(iltot_b1) = bound(iltot_b1) + sqrt(dot_product(ltot_current_part,ltot_current_part))
          bound(ietot_b1) = bound(ietot_b1) + etot_current_part

       else

          !num parts, total mass, ang mom and total energy respectively - unbound
          bound(in_u1) = bound(in_u1) + 1
          bound(imass_u1) = bound(imass_u1) + particlemass
          bound(iltot_u1) = bound(iltot_u1) + sqrt(dot_product(ltot_current_part,ltot_current_part))
          bound(ietot_u1) = bound(ietot_u1) + etot_current_part

       endif

       !bound criterion INCLUDING thermal energy
       if (poten(i) + part_ps_pot + part_kin + etherm_current_part < 0.0) then

          !num parts, total mass, ang mom and total energy respectively - bound
          bound(in_b2) = bound(in_b2) + 1
          bound(imass_b2) = bound(imass_b2) + particlemass
          bound(iltot_b2) = bound(iltot_b2) + sqrt(dot_product(ltot_current_part,ltot_current_part))
          bound(ietot_b2) = bound(ietot_b2) + etot_current_part

       else

          !num parts, total mass, ang mom and total energy respectively - unbound
          bound(in_u2) = bound(in_u2) + 1
          bound(imass_u2) = bound(imass_u2) + particlemass
          bound(iltot_u2) = bound(iltot_u2) + sqrt(dot_product(ltot_current_part,ltot_current_part))
          bound(ietot_u2) = bound(ietot_u2) + etot_current_part

       endif

       !bound criterion INCLUDING enthalpy
       if (poten(i) + part_ps_pot + part_kin &
+ etherm_current_part + ppart/rhopart < 0.0) then

          !num parts, total mass, ang mom and total energy respectively - bound
          bound(in_b3) = bound(in_b3) + 1
          bound(imass_b3) = bound(imass_b3) + particlemass
          bound(iltot_b3) = bound(iltot_b3) + sqrt(dot_product(ltot_current_part,ltot_current_part))
          bound(ietot_b3) = bound(ietot_b3) + etot_current_part

       else

          !num parts, total mass, ang mom and total energy respectively - unbound
          bound(in_u3) = bound(in_u3) + 1
          bound(imass_u3) = bound(imass_u3) + particlemass
          bound(iltot_u3) = bound(iltot_u3) + sqrt(dot_product(ltot_current_part,ltot_current_part))
          bound(ietot_u3) = bound(ietot_u3) + etot_current_part

       endif

    enddo

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

    call write_time_file('boundunbound_vs_time', columns, time, bound, 24, num)

 case(3) !Energies and bound mass

    if (dump_number == 0) then
       print "(4(/,a))", 'Who would want CE energies', &
                            'Must answer me', &
                            'These questions three', &
                            'Ere the new results he see'

       call prompt('1. Was this in a corotating frame?',switch(1),.false.)

       if (switch(1)) then
          call prompt('  1a. Enter value of initial eccentricity : ',e)
          a = xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1)
          sep = a*(1. + e)
          mtot = sum(xyzmh_ptmass(4,:)) + npart*particlemass
          vmag = sqrt(a*(1.-e**2)*mtot)/sep
          omega_corotate = vmag/a
          omega_corotate = Dnint(omega_corotate*10000000.0)/10000000.0
          omega_vec = (/ 0.,0.,omega_corotate /)
       else
          omega_vec = (/ 0.,0.,0. /)
       endif

       call prompt('2. Would you like to output bound particle files?', switch(2))

       if (switch(2)) then
          call prompt('3. Would you like to use thermal energy in the computation of the bound/unbound status?', switch(3),.false.)
       endif
    endif

    encomp(5:) = 0.

    call compute_energies(time)

    ekin = 0.

    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    do i=1,npart
       if (xyzh(4,i)  >=  0) then
          encomp(ipot_pp) = encomp(ipot_pp) + poten(i)
          encomp(ipot_env) = encomp(ipot_env) + poten(i)

          call separation_vector(xyzh(1:3,i),com_pos(1:3),r_part_com(1:4))

          rhopart = rhoh(xyzh(4,i), particlemass)

          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,&
                                  r_part_com(1),r_part_com(2),r_part_com(3),vxyzu(4,i))

          call cross(omega_vec,r_part_com(1:3),omegacrossr)

          call separation_vector(vxyzu(1:3,i) + omegacrossr(1:3),com_v(1:3),v_part_com(1:4))

          part_kin = 0.5 * particlemass * v_part_com(4)**2
          part_ps_pot = 0.

          call cross(r_part_com(1:3),particlemass * v_part_com(1:3),rcrossmv)

          jz = rcrossmv(3)
          encomp(ijz_tot) = encomp(ijz_tot) + jz

          do k=1,nptmass
             r_sink_part = separation(xyzmh_ptmass(1:3,k),xyzh(1:3,i))**2

             hsoftk = max(xyzmh_ptmass(ihsoft,k),xyzh(4,i))
             hsoft1 = 1.0/hsoftk
             hsoft21= hsoft1**2
             q2i    = r_sink_part*hsoft21
             qi     = sqrt(q2i)

             call kernel_softening(q2i,qi,psoft,fsoft)

             part_ps_pot = part_ps_pot + psoft*hsoft1*xyzmh_ptmass(4,k)*particlemass

             if (k == 1) then
                encomp(ipot_env) = encomp(ipot_env) + psoft*hsoft1*xyzmh_ptmass(4,k)*particlemass
             endif

             if (sqrt(r_sink_part) < 80) then
                inearsink = .true.
             endif
          enddo

          encomp(ipot_ps) = encomp(ipot_ps) + part_ps_pot

          part_therm = particlemass * vxyzu(4,i)

          if (switch(2)) then
             boundparts(1:3,i) = xyzh(1:3,i)
             boundparts(4:6,i) = vxyzu(1:3,i)
          endif

          if (switch(3)) then
             boundparts(7,i) = (part_kin + part_ps_pot + poten(i) + part_therm)/particlemass

             if (part_kin + part_ps_pot + poten(i) + part_therm < 0) then
                encomp(ikin_bound) = encomp(ikin_bound) + part_kin
                encomp(imass_bound) = encomp(imass_bound) + particlemass
                encomp(ijz_bound) = encomp(ijz_bound) + jz
                boundparts(8,i) = -1
                radvel = dot_product(v_part_com(1:3),r_part_com(1:3)) / r_part_com(4)

                if (inearsink .eqv. .false.) then

                   if (radvel < 0.) then
                      encomp(fallbackmass) = encomp(fallbackmass) + particlemass
                      encomp(fallbackmom) = encomp(fallbackmom) + particlemass * radvel
                   endif
                endif

             else
                encomp(ikin_unbound) = encomp(ikin_unbound) + part_kin
                encomp(imass_unbound) = encomp(imass_unbound) + particlemass
                encomp(ijz_unbound) = encomp(ijz_unbound) + jz
                boundparts(8,i) = 1
             endif

          else
             boundparts(7,i) = (part_kin + part_ps_pot + poten(i))/particlemass

             if (part_kin + part_ps_pot + poten(i) < 0) then
                encomp(ikin_bound) = encomp(ikin_bound) + part_kin
                encomp(imass_bound) = encomp(imass_bound) + particlemass
                encomp(ijz_bound) = encomp(ijz_bound) + jz
                boundparts(8,i) = -1
                radvel = dot_product(v_part_com(1:3),r_part_com(1:3)) / r_part_com(4)

                if (inearsink .eqv. .false.) then

                   if (radvel < 0.) then
                      encomp(fallbackmass) = encomp(fallbackmass) + particlemass
                      encomp(fallbackmom) = encomp(fallbackmom) + particlemass * radvel
                   endif
                endif

             else
                encomp(ikin_unbound) = encomp(ikin_unbound) + part_kin
                encomp(imass_unbound) = encomp(imass_unbound) + particlemass
                encomp(ijz_unbound) = encomp(ijz_unbound) + jz
                boundparts(8,i) = 1
             endif
          endif
       endif
    enddo

    !open(26,file=trim(dumpfile)//".divv",status='replace',form='unformatted')
    !write(26) (real(boundparts(7,i),kind=4),i=1,npart)
    !close(26)

    do i=1,nptmass
       call cross(omega_vec,xyzmh_ptmass(:3,i),omegacrossr)

       call separation_vector(vxyz_ptmass(1:3,i) + omegacrossr(1:3),com_v(1:3),sink_vel(1:4))
       call separation_vector(xyzmh_ptmass(1:3,i),com_pos(1:3),r_part_com(1:4))

       call cross(r_part_com(:3),xyzmh_ptmass(4,i) * sink_vel(:3),rcrossmv)

       jz = rcrossmv(3)
       encomp(ijz_tot) = jz + encomp(ijz_tot)
       encomp(ijz_orb) = jz + encomp(ijz_orb)
       encomp(ikin_sink) = encomp(ikin_sink) + 0.5 * xyzmh_ptmass(4,i) * sink_vel(4)**2
       if (i==2) encomp(iorb_comp) = encomp(iorb_comp) + 0.5 * xyzmh_ptmass(4,i) * sink_vel(4)**2
    enddo

    do i=1,nptmass-1
       do j=i+1,nptmass
          r_sink_sink = separation(xyzmh_ptmass(1:3,i),xyzmh_ptmass(1:3,j))

          encomp(ipot_sink) = encomp(ipot_sink) - xyzmh_ptmass(4,i) * xyzmh_ptmass(4,j) / r_sink_sink
          if (i==1 .and. j==2) encomp(iorb_comp) = encomp(iorb_comp) - xyzmh_ptmass(4,i) * xyzmh_ptmass(4,j) / r_sink_sink
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

    call write_time_file('energy', columns, time, encomp, 23, num)

    if (switch(2)) then
       columns = (/'     x-coord', '     y-coord', '     z-coord', &
                   '       x-vel', '       y-vel', '       z-vel', &
                   ' spec energy', '  bound bool'/)
       call write_file('boundparts', 'bound', columns, boundparts, npart, 8, num)
    endif

    unitnum = unitnum + 1

 case(4) !Profile from COM (can be used for stellar profile)
    profile_vector=(/1.,0.,0./)
    call prompt('Choose profile vector x-component ',profile_vector(1))
    call prompt('Choose profile vector y-component ',profile_vector(2))
    call prompt('Choose profile vector z-component ',profile_vector(3))

    call stellar_profile(time,particlemass,npart,xyzh,vxyzu,profile,profile_vector)

    columns = (/'      radius',&
                '  mass coord',&
                '     azimuth',&
                ' int. energy',&
                '     density',&
                '    pressure',&
                '    velocity',&
                '   rad. vel.',&
                ' sound speed',&
                '        temp',&
                '       kappa',&
                '         mfp',&
                '      energy',&
                '    HII frac',&
                '   HeII frac',&
                '  HeIII frac'/)

    print*, 'Writing file...'
    write(name_in, "(a,i1,i1,i1)") 'ray_profile_',int(profile_vector(1:3))
    call write_file(name_in, 'profile', columns, profile, size(profile(1,:)), 16, num)

    print*, 'File written'

    deallocate(profile)

    unitnum = unitnum + 1

 case(5) !Finding volume equivalent radius of particles (only useful for finding radius of star)

    totvol = 0

    do i=1,npart
       rhopart = rhoh(xyzh(4,i), particlemass)
       totvol = totvol + particlemass / rhopart
    enddo

    VER(1) = (3. * totvol/(4. * pi))**(1./3.)

    print*, 'Total volume of star is ', totvol
    print*, 'Volume equivalent radius is ', VER(1)

    columns = (/'vol. eq. rad.'/)

    call write_time_file('VERStar', columns, time, VER, 1, num)

 case(6) !Mass within roche lobes

    MRL(1:6) = 0

    !do i=1,npart
    !   total_mass = total_mass + particlemass
    !enddo
    !m2 = xyzmh_ptmass(4,2)
    !m1 = total_mass + xyzmh_ptmass(4,1)

    !separation(1:3) = xyzmh_ptmass(1:3,1) - xyzmh_ptmass(1:3,2)
    !separation(4) = sqrt(dot_product(separation(1:3),separation(1:3)))

    !RL1 = Rochelobe_estimate(m1,m2,separation(4))
    !RL2 = Rochelobe_estimate(m2,m1,separation(4))


    RL1 = 10
    RL2 = 10

    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    do i=1,npart
       if (xyzh(4,i)  >=  0) then
          v_part_com(1:3) = vxyzu(1:3,i) - com_v(1:3)
          v_part_com(4) = sqrt(dot_product(v_part_com(1:3),v_part_com(1:3)))
          part_kin = 0.5 * particlemass * v_part_com(4)**2
          part_ps_pot = 0.

          do k=1,nptmass
             r_sink_part = (xyzmh_ptmass(1,k) - xyzh(1,i))**2 + &
                              (xyzmh_ptmass(2,k) - xyzh(2,i))**2 + &
                              (xyzmh_ptmass(3,k) - xyzh(3,i))**2

             hsoftk = max(xyzmh_ptmass(ihsoft,k),xyzh(4,i))
             hsoft1 = 1.0/hsoftk
             hsoft21= hsoft1**2
             q2i    = r_sink_part*hsoft21
             qi     = sqrt(q2i)

             call kernel_softening(q2i,qi,psoft,fsoft)

             part_ps_pot = part_ps_pot + psoft*hsoft1*xyzmh_ptmass(4,k)*particlemass

          enddo

          do k=1,nptmass
             r_sink_part = sqrt((xyzmh_ptmass(1,k) - xyzh(1,i))**2 + &
                                   (xyzmh_ptmass(2,k) - xyzh(2,i))**2 + &
                                   (xyzmh_ptmass(3,k) - xyzh(3,i))**2)
             if (r_sink_part < RL1 .and. k==1) then
                MRL(1) = MRL(1) + particlemass
                if (part_kin + part_ps_pot + poten(i) < 0) then
                   MRL(2) = MRL(2) + particlemass
                endif
             elseif (r_sink_part < RL2 .and. k==2) then
                MRL(3) = MRL(3) + particlemass
                if (part_kin + part_ps_pot + poten(i) < 0) then
                   MRL(4) = MRL(4) + particlemass
                endif
             else
                MRL(5) = MRL(5) + particlemass
                if (part_kin + part_ps_pot + poten(i) < 0) then
                   MRL(6) = MRL(6) + particlemass
                endif
             endif
          enddo

       endif
    enddo

    columns = (/' Mass in RL1',&
                   '  B Mass RL1',&
                   ' Mass in RL2',&
                   '  B Mass RL2',&
                   'Mass ejected',&
                   'B Mass eject'/)

    call write_time_file('MassDist', columns, time, MRL, 6, num)

 case(7) !Histograms
    if (dump_number == 0) then
       print "(4(/,a,/))",' 1) Radius (from centre of mass)', &
                ' 2) Velocity (with respect to centre of mass)', &
                ' 3) Radial velocity (with respect to centre of mass)', &
                ' 4) Velocity with respect to sink as function of radius from sink'
       call prompt('Choose histogram type ',histogram_number,1,4)
    endif

    select case(histogram_number)

    case(1) !Number of particles vs radius
       call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

       do i=1,npart
          velocity_wrt_com(1,i) = sqrt((vxyzu(1,i) - com_v(1))**2)
          velocity_wrt_com(2,i) = sqrt((vxyzu(2,i) - com_v(2))**2)
          velocity_wrt_com(3,i) = sqrt((vxyzu(3,i) - com_v(3))**2)
          velocity_wrt_com(4,i) = sqrt(velocity_wrt_com(1,i)**2 + velocity_wrt_com(2,i)**2 + velocity_wrt_com(3,i)**2)
       enddo

       if (dump_number == 0) then
          print*, 'Write the minimum velocity limit (default=minimum):'
          write(*, '(2x, a26, F15.10)') 'Minimum value in data is ', minval(velocity_wrt_com(4,:))
          read(*,'(a)') val1
          print*, 'Write the maximum velocity limit (default=maximum):'
          write(*, '(2x, a26, F15.10)') 'Maximum value in data is ', maxval(velocity_wrt_com(4,:))
          read(*,'(a)') val2
       endif

       if (val1 == '') then
          constraint_min = minval(velocity_wrt_com(4,:))
       else
          read(val1,*) constraint_min
       endif
       if (val2 == '') then
          constraint_max = maxval(velocity_wrt_com(4,:))
       else
          read(val2,*) constraint_max
       endif

       !Set the distances of particles from the centre of mass (x,y,z and total distances)
       do i=1,npart
          if ((velocity_wrt_com(4,i)  >  constraint_min) .and. (velocity_wrt_com(4,i)  <  constraint_max)) then
             distance_from_com(1,i) = sqrt((xyzh(1,i) - com_pos(1))**2)
             distance_from_com(2,i) = sqrt((xyzh(2,i) - com_pos(2))**2)
             distance_from_com(3,i) = sqrt((xyzh(3,i) - com_pos(3))**2)
             distance_from_com(4,i) = sqrt(distance_from_com(1,i)**2 + &
distance_from_com(2,i)**2 + distance_from_com(3,i)**2)
          endif
       enddo

       if (dump_number == 0) then
          call histogram_setup(distance_from_com,npart,bin_res,bin_size,histogram_bins,.true.)
       else
          call histogram_setup(distance_from_com,npart,bin_res,bin_size,histogram_bins,.false.)
       endif

       !Count number of particles within a increasing radial shells from COM
       do i=1,npart
          do j=1,bin_res
             if ((distance_from_com(4,i)  >  (histogram_bins(1,j) - bin_size)) .and. (distance_from_com(4,i)&
  <=  histogram_bins(1,j))) then
                histogram_bins(2,j) = histogram_bins(2,j) + 1
             endif
          enddo
       enddo

       allocate(columns(2))
       columns = (/'radius','number'/)
       call write_file('rad_dist','histograms', columns, histogram_bins, bin_res, 2, num)

    case(2) !Number of particles vs velocity

       call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

       do i=1,npart
          distance_from_com(1,i) = sqrt((xyzh(1,i) - com_pos(1))**2)
          distance_from_com(2,i) = sqrt((xyzh(2,i) - com_pos(2))**2)
          distance_from_com(3,i) = sqrt((xyzh(3,i) - com_pos(3))**2)
          distance_from_com(4,i) = sqrt(distance_from_com(1,i)**2 + distance_from_com(2,i)**2 + distance_from_com(3,i)**2)
       enddo

       if (dump_number == 0) then
          print*, 'Write the minimum radius limit (default=minimum):'
          write(*, '(2x, a26, F15.10)') 'Minimum value in data is ', minval(distance_from_com(4,:))
          read(*,'(a)') val1
          print*, 'Write the maximum radius limit (default=maximum):'
          write(*, '(2x, a26, F15.10)') 'Maximum value in data is ', maxval(distance_from_com(4,:))
          read(*,'(a)') val2
       endif

       if (val1 == '') then
          constraint_min = minval(distance_from_com(4,:))
       else
          read(val1,*) constraint_min
       endif
       if (val2 == '') then
          constraint_max = maxval(distance_from_com(4,:))
       else
          read(val2,*) constraint_max
       endif

       !Set the distances of particles from the centre of mass (x,y,z and total distances)
       do i=1,npart
          if ((distance_from_com(4,i)  >  constraint_min) .and. (distance_from_com(4,i)  <  constraint_max)) then
             velocity_wrt_com(1,i) = sqrt((vxyzu(1,i) - com_v(1))**2)
             velocity_wrt_com(2,i) = sqrt((vxyzu(2,i) - com_v(2))**2)
             velocity_wrt_com(3,i) = sqrt((vxyzu(3,i) - com_v(3))**2)
             velocity_wrt_com(4,i) = sqrt(velocity_wrt_com(1,i)**2 + velocity_wrt_com(2,i)**2 + velocity_wrt_com(3,i)**2)
          endif
       enddo

       if (dump_number == 0) then
          call histogram_setup(velocity_wrt_com,npart,bin_res,bin_size,histogram_bins,.true.)
       else
          call histogram_setup(velocity_wrt_com,npart,bin_res,bin_size,histogram_bins,.false.)
       endif

       do i=1,npart
          do j=1,bin_res
             if ((velocity_wrt_com(4,i)  >  (histogram_bins(1,j) - bin_size)) .and. (velocity_wrt_com(4,i)&
  <=  histogram_bins(1,j))) then
                histogram_bins(2,j) = histogram_bins(2,j) + 1
             endif
          enddo
       enddo

       allocate(columns(2))
       columns = (/'velocity','  number'/)
       call write_file('vel_dist','histograms', columns, histogram_bins, bin_res, 2, num)

    case(3)!Radial velocity

       call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

       do i=1,npart
          distance_from_com(1,i) = sqrt((xyzh(1,i) - com_pos(1))**2)
          distance_from_com(2,i) = sqrt((xyzh(2,i) - com_pos(2))**2)
          distance_from_com(3,i) = sqrt((xyzh(3,i) - com_pos(3))**2)
          distance_from_com(4,i) = sqrt(distance_from_com(1,i)**2 + distance_from_com(2,i)**2 + distance_from_com(3,i)**2)
       enddo

       if (dump_number == 0) then
          print*, 'Write the minimum radius limit (default=minimum):'
          write(*, '(2x, a26, F15.10)') 'Minimum value in data is ', minval(distance_from_com(4,:))
          read(*,'(a)') val1
          print*, 'Write the maximum radius limit (default=maximum):'
          write(*, '(2x, a26, F15.10)') 'Maximum value in data is ', maxval(distance_from_com(4,:))
          read(*,'(a)') val2
       endif

       if (val1 == '') then
          constraint_min = minval(distance_from_com(4,:))
       else
          read(val1,*) constraint_min
       endif
       if (val2 == '') then
          constraint_max = maxval(distance_from_com(4,:))
       else
          read(val2,*) constraint_max
       endif

       !Set the distances of particles from the centre of mass (x,y,z and total distances)
       do i=1,npart
          if ((distance_from_com(4,i)  >  constraint_min) .and. (distance_from_com(4,i)  <  constraint_max)) then
             velocity_wrt_com(1,i) = sqrt((vxyzu(1,i) - com_v(1))**2)
             velocity_wrt_com(2,i) = sqrt((vxyzu(2,i) - com_v(2))**2)
             velocity_wrt_com(3,i) = sqrt((vxyzu(3,i) - com_v(3))**2)
             velocity_wrt_com(4,i) = (vxyzu(1,i) - com_v(1))*(xyzh(1,i) - com_pos(1))/distance_from_com(4,i) +&
                                               (vxyzu(2,i) - com_v(2))*(xyzh(2,i) - com_pos(2))/distance_from_com(4,i) +&
                                               (vxyzu(3,i) - com_v(3))*(xyzh(3,i) - com_pos(3))/distance_from_com(4,i)
          endif
       enddo

       if (dump_number == 0) then
          call histogram_setup(velocity_wrt_com,npart,bin_res,bin_size,histogram_bins,.true.)
       else
          call histogram_setup(velocity_wrt_com,npart,bin_res,bin_size,histogram_bins,.false.)
       endif

       do i=1,npart
          do j=1,bin_res
             if ((velocity_wrt_com(4,i)  >  (histogram_bins(1,j) - bin_size)) .and. (velocity_wrt_com(4,i)&
  <=  histogram_bins(1,j))) then
                histogram_bins(2,j) = histogram_bins(2,j) + 1
             endif
          enddo
       enddo

       allocate(columns(2))
       columns = (/'radial velocity','         number'/)
       call write_file('rad_vel_dist','histograms', columns, histogram_bins, bin_res, 2, num)


    case(4)!Corotation
       call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

       do i=1,nptmass
          distance_from_com(1,i) = sqrt((xyzmh_ptmass(1,i) - com_pos(1))**2)
          distance_from_com(2,i) = sqrt((xyzmh_ptmass(2,i) - com_pos(2))**2)
          distance_from_com(3,i) = sqrt((xyzmh_ptmass(3,i) - com_pos(3))**2)
          distance_from_com(4,i) = sqrt(distance_from_com(1,i)**2 + distance_from_com(2,i)**2 + distance_from_com(3,i)**2)


          velocity_wrt_com(1,i) = sqrt((vxyz_ptmass(1,i) - com_v(1))**2)
          velocity_wrt_com(2,i) = sqrt((vxyz_ptmass(2,i) - com_v(2))**2)
          velocity_wrt_com(3,i) = sqrt((vxyz_ptmass(3,i) - com_v(3))**2)
          velocity_wrt_com(4,i) = sqrt(velocity_wrt_com(1,i)**2 + velocity_wrt_com(2,i)**2 + velocity_wrt_com(3,i)**2)
          projection = ((velocity_wrt_com(1,i))*(distance_from_com(1,i))+&
                                 (velocity_wrt_com(2,i))*(distance_from_com(2,i))+&
                                 (velocity_wrt_com(3,i))*(distance_from_com(3,i)))/distance_from_com(4,i)**2

          vel_rad_comp(1,i) = projection * distance_from_com(1,i)
          vel_rad_comp(2,i) = projection * distance_from_com(2,i)
          vel_rad_comp(3,i) = projection * distance_from_com(3,i)
          vel_rad_comp(4,i) = sqrt(vel_rad_comp(1,i)**2 + vel_rad_comp(2,i)**2 + vel_rad_comp(3,i)**2)

          vel_tan_comp(1,i) = velocity_wrt_com(1,i) - vel_rad_comp(1,i)
          vel_tan_comp(2,i) = velocity_wrt_com(2,i) - vel_rad_comp(2,i)
          vel_tan_comp(3,i) = velocity_wrt_com(3,i) - vel_rad_comp(3,i)
          vel_tan_comp(4,i) = sqrt(vel_tan_comp(1,i)**2 + vel_tan_comp(2,i)**2 + vel_tan_comp(3,i)**2)

       enddo

       ang_vel = ((vel_tan_comp(4,1)/distance_from_com(4,1)) + (vel_tan_comp(4,2)/distance_from_com(4,2)))/2.0
       if (dump_number == 0) then
          print*, 'Radius from COM within which particles will be counted (default=100):'
          write(*, '(2x, a15, F15.10)') 'The default is ', 100.0
          read(*,FMT='(a)') val1
       endif

       if (val1 == '') then
          constraint_max = 100
       else
          read(val1,*) constraint_max
       endif


       do i=1,npart
          distance_from_com(1,i) = sqrt((xyzh(1,i) - com_pos(1))**2)
          distance_from_com(2,i) = sqrt((xyzh(2,i) - com_pos(2))**2)
          distance_from_com(3,i) = sqrt((xyzh(3,i) - com_pos(3))**2)
          distance_from_com(4,i) = sqrt(distance_from_com(1,i)**2 + distance_from_com(2,i)**2 + distance_from_com(3,i)**2)

          if (distance_from_com(4,i)  >  constraint_max) then
             velocity_wrt_com(1,i) = sqrt((vxyzu(1,i) - com_v(1))**2)
             velocity_wrt_com(2,i) = sqrt((vxyzu(2,i) - com_v(2))**2)
             velocity_wrt_com(3,i) = sqrt((vxyzu(3,i) - com_v(3))**2)
             velocity_wrt_com(4,i) = sqrt(velocity_wrt_com(1,i)**2 + velocity_wrt_com(2,i)**2 + velocity_wrt_com(3,i)**2)
             projection = ((velocity_wrt_com(1,i))*(distance_from_com(1,i))+&
                                     (velocity_wrt_com(2,i))*(distance_from_com(2,i))+&
                                     (velocity_wrt_com(3,i))*(distance_from_com(3,i)))/distance_from_com(4,i)**2

             vel_rad_comp(1,i) = projection * distance_from_com(1,i)
             vel_rad_comp(2,i) = projection * distance_from_com(2,i)
             vel_rad_comp(3,i) = projection * distance_from_com(3,i)
             vel_rad_comp(4,i) = sqrt(vel_rad_comp(1,i)**2 + vel_rad_comp(2,i)**2 + vel_rad_comp(3,i)**2)

             vel_tan_comp(1,i) = velocity_wrt_com(1,i) - vel_rad_comp(1,i)
             vel_tan_comp(2,i) = velocity_wrt_com(2,i) - vel_rad_comp(2,i)
             vel_tan_comp(3,i) = velocity_wrt_com(3,i) - vel_rad_comp(3,i)
             vel_tan_comp(4,i) = sqrt(vel_tan_comp(1,i)**2 + vel_tan_comp(2,i)**2 + vel_tan_comp(3,i)**2) &
                                           / (distance_from_com(4,i) * ang_vel)
          endif
       enddo


       if (dump_number == 0) then
          call histogram_setup(vel_tan_comp,npart,bin_res,bin_size,histogram_bins,.true.)
       else
          call histogram_setup(vel_tan_comp,npart,bin_res,bin_size,histogram_bins,.false.)
       endif

       do i=1,npart
          do j=1,bin_res
             if ((vel_tan_comp(4,i)  >  (histogram_bins(1,j) - bin_size)) .and. (velocity_wrt_com(4,i)&
  <=  histogram_bins(1,j))) then
                histogram_bins(2,j) = histogram_bins(2,j) + 1
             endif
          enddo
       enddo

       print*, histogram_bins
       columns = (/'ang vel frac','      number'/)
       call write_file('ang_vel_frac', 'histograms', columns, histogram_bins, bin_res, 2, num)
    end select


 case(8)
    totvol = 0
    densvol = 0
    if (dump_number == 0) then
       surfacedensity = rhoh(xyzh(4,1), particlemass)
       do i=1,npart
          rhopart = rhoh(xyzh(4,i), particlemass)
          if (rhopart < surfacedensity) then
             surfacedensity = rhopart
          endif
       enddo
    endif

    do i=1,npart
       rhopart = rhoh(xyzh(4,i), particlemass)
       totvol = totvol + particlemass / rhopart
       if (rhopart > surfacedensity) then
          densvol = densvol + particlemass / rhopart
       endif
    enddo

    star_stabilisation_params(ivoleqrad) = (3. * totvol/(4. * pi))**(1./3.)
    star_stabilisation_params(idensrad)  = (3. * densvol/(4. * pi))**(1./3.)

    if (dump_number == 0) then
       init_radius = star_stabilisation_params(ivoleqrad)
    endif

    star_stabilisation_params(imassout) = 0.
    total_mass = xyzmh_ptmass(4,1)
    do i=1,npart
       r_sink_part = sqrt((xyzmh_ptmass(1,1) - xyzh(1,i))**2.0 + &
                         (xyzmh_ptmass(2,1) - xyzh(2,i))**2.0 + &
                         (xyzmh_ptmass(3,1) - xyzh(3,i))**2.0)
       if (r_sink_part > init_radius) then
          star_stabilisation_params(imassout) = star_stabilisation_params(imassout) + particlemass
       endif
       total_mass = total_mass + particlemass
    enddo

    star_stabilisation_params(imassfracout) = star_stabilisation_params(imassout) / total_mass

    columns = (/'vol. eq. rad',&
               ' density rad',&
               'mass outside',&
               'frac outside'/)

    call write_time_file('star_stabilisation', columns, time, star_stabilisation_params, 4, dump_number)

    if (dump_number == 0) then
       call prompt('Make profiles too (this will make the analysis take much longer)?',switch(1),.false.)
    endif
    if (switch(1)) then
       call stellar_profile(time,particlemass,npart,xyzh,vxyzu,profile)
       columns = (/'      radius',&
                  '  mass coord',&
                  '     azimuth',&
                  ' int. energy',&
                  '     density',&
                  '    pressure',&
                  '    velocity',&
                  '   rad. vel.',&
                  ' sound speed'/)

       call write_file('profile', 'profile', columns, profile, npart, 9, num)
    endif
    unitnum = unitnum + 1

 case(9) !Units
    write(*,"(/,3(a,es10.3,1x),a)") '     Mass: ',umass,    'g       Length: ',udist,  'cm     Time: ',utime,'s'
    write(*,"(3(a,es10.3,1x),a)") '  Density: ',unit_density, 'g/cm^3  Energy: ',unit_energ,'erg    En/m: ',unit_ergg,'erg/g'
    write(*,"(3(a,es10.3,1x),a)") ' Velocity: ',unit_velocity,'cm/s    Bfield: ',unit_Bfield,'G  Pressure: ',&
                                     unit_pressure,'g/cm s^2'
    write(*,"(2(a,es10.3,1x),/)")   '        G: ', gg*umass*utime**2/udist**3,'             c: ',c*utime/udist

 case(10) !It's a thing

    C_force = 0.25
    dtsg    = huge(dtsg)
    minsep1  = huge(minsep1)
    minsep2  = huge(minsep2)
    fextx = 0.
    fexty = 0.
    fextz = 0.
    do i = 1,npart
       call get_accel_sink_gas(nptmass,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),xyzmh_ptmass,&
                        fextx,fexty,fextz,phii,particlemass,fxyz_ptmass_thread,fonrmaxi,dtphi2i)
       f2i   = fextx**2 + fexty**2 + fextz**2
       hi    = min(hsoft,xyzh(4,i))
       dtsgi = C_force*sqrt(hi/sqrt(f2i))
       dtsg  = min(dtsg,dtsgi)


       minsep1i = sqrt((xyzh(1,i) - xyzmh_ptmass(1,1))**2 &
                      + (xyzh(2,i) - xyzmh_ptmass(2,1))**2 &
                      + (xyzh(3,i) - xyzmh_ptmass(3,1))**2)

       minsep2i = sqrt((xyzh(1,i) - xyzmh_ptmass(1,2))**2 &
                      + (xyzh(2,i) - xyzmh_ptmass(2,2))**2 &
                      + (xyzh(3,i) - xyzmh_ptmass(3,2))**2)

       minsep1 = min(minsep1,minsep1i)
       minsep2 = min(minsep2,minsep2i)
    enddo

    unitnum = 1001
    write(fileout,"(2a,i3.3,a)") 'timestep.ev'

    if (num == 0) then
       open(unit=unitnum, file=fileout, status='replace')

       write(unitnum,"('#',1x,4('[',i2.1,a18,']',2x))") &
              1,'time',            &
              2,'sg timestep',     &
              3,'minimum sep 1',   &
              4,'minimum sep 2'

       close(unit=unitnum)

       unitnum = unitnum + 1
    endif

    open(unit=unitnum, file=fileout, position='append')

    write(unitnum,"(4(4x,es18.11e2,2x))") &
            time,        &
            dtsg,        &
            minsep1,     &
            minsep2

    close(unit=unitnum)

 case(11) !Interpolated Ray Profile

    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    ray = (/ 1, 0, 0 /)

    ray = ray / sqrt(dot_product(ray,ray))
    print*, ray
    far_dist = 0.

    do i=1,npart
       xyz(1:3,i) = xyzh(1:3,i) - com_pos(1:3)
       proj_mag = dot_product(xyz(1:3,i),ray(1:3))
       proj = proj_mag * ray
       proj_orth(1:3) = xyz(1:3,i) - proj(1:3)
       proj_dist = separation(proj,(/0.,0.,0./))
       proj_orth_dist = separation(proj_orth,(/0.,0.,0./))
       qi = proj_orth_dist / xyzh(4,i)
       if (proj_dist > far_dist .and. qi < radkern .and. proj_mag > 0.) then
          far_dist = proj_dist
          far_i = i
       endif
    enddo

    n_sample = int(10**nint(log10(real(npart))))/10

    allocate(ray_profile(9,n_sample))

    ray_profile = 0.
    do i=1,n_sample
       ray_profile(1,i) = i * far_dist / n_sample
       xyz_profile(1:3) = ray_profile(1,i) * ray(1:3)
       do j=1,npart
          qi = separation(xyz_profile(1:3),xyz(1:3,j)) / xyzh(4,j)
          if (qi < radkern) then
             Wj = (cnormk * wkern(qi**2,qi) / (xyzh(4,j)**3))
             ray_profile(2,i) = ray_profile(2,i) + particlemass * Wj
             ray_profile(3,i) = ray_profile(3,i) + particlemass * (vxyzu(4,j) / rhoh(xyzh(4,j), particlemass)) * Wj
             call separation_vector(xyzh(1:3,j),com_pos(1:3),r_part_com)
             call separation_vector(vxyzu(1:3,j),com_v(1:3),v_part_com)
             ray_profile(5,i) = ray_profile(5,i) + particlemass * (v_part_com(4) / rhoh(xyzh(4,j), particlemass)) * Wj
             radvel = dot_product(v_part_com(1:3),r_part_com(1:3)) / r_part_com(4)
             ray_profile(6,i) = ray_profile(6,i) + particlemass * (radvel / rhoh(xyzh(4,j), particlemass)) * Wj

          endif
       enddo

       call get_eos_pressure_temp_mesa(X_in,ray_profile(2,i)*unit_density,&
                                       ray_profile(3,i) * unit_ergg,ray_profile(4,i),ray_profile(7,i))

       call get_eos_kappa_mesa(X_in,ray_profile(2,i)*unit_density,ray_profile(7,i),ray_profile(8,i),kappat,kappar,ierr)

       !ray_profile(4,i) = ponrhoi * ray_profile(2,i)
       ray_profile(9,i) = 1. / (ray_profile(8,i) * ray_profile(2,i) * unit_density)
       if (mod(i*10,n_sample)==0) write(*,"(2A,I3,A)",advance='no') achar(13), "Profile is ", i * 100 / n_sample, "% complete."
    enddo

    ray_profile(1,:) = ray_profile(1,:) * udist
    ray_profile(2,:) = ray_profile(2,:) * unit_density
    ray_profile(3,:) = ray_profile(3,:) * unit_ergg
    ray_profile(5:6,:) = ray_profile(5:6,:) * unit_velocity

    columns = (/'      radius',&
                '     density',&
                ' int. energy',&
                '    pressure',&
                '    velocity',&
                '     rad vel',&
                ' temperature',&
                '       kappa',&
                '         MFP'/)

    call write_file('rayprofile', 'profile', columns, ray_profile, n_sample, 9, num)

    print*, 'File written'

    unitnum = unitnum + 1


 case(12) !Output .divv

    call compute_energies(time)

    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    do i=1,npart
       if (xyzh(4,i)  >=  0) then
          call separation_vector(xyzh(1:3,i),com_pos(1:3),r_part_com(1:4))
          call separation_vector(vxyzu(1:3,i) + omegacrossr(1:3),com_v(1:3),v_part_com(1:4))

          part_kin = 0.5 * particlemass * v_part_com(4)**2
          part_ps_pot = 0.

          do k=1,nptmass
             r_sink_part = separation(xyzmh_ptmass(1:3,k),xyzh(1:3,i))**2

             hsoftk = max(xyzmh_ptmass(ihsoft,k),xyzh(4,i))
             hsoft1 = 1.0/hsoftk
             hsoft21= hsoft1**2
             q2i    = r_sink_part*hsoft21
             qi     = sqrt(q2i)

             call kernel_softening(q2i,qi,psoft,fsoft)

             part_ps_pot = part_ps_pot + psoft*hsoft1*xyzmh_ptmass(4,k)*particlemass
          enddo

          bound_energy(i) = (part_kin + part_ps_pot + poten(i))/particlemass
          rhopart = rhoh(xyzh(4,i), particlemass)

          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,&
                                  r_part_com(1),r_part_com(2),r_part_com(3),vxyzu(4,i))
          if (ieos==10) then
             call get_eos_pressure_temp_mesa(X_in,rhopart*unit_density,vxyzu(4,i) * unit_ergg,pres,temp(i))
             call get_eos_kappa_mesa(X_in,rhopart*unit_density,temp(i),kappa(i),kappat,kappar,ierr)
          else
             kappa(i) = 1
          endif
       endif
    enddo

    open(26,file=trim(dumpfile)//".divv",status='replace',form='unformatted')
    write(26) (real(bound_energy(i),kind=4),i=1,npart)
    write(26) (real(kappa(i),kind=4),i=1,npart)
    write(26) (real(temp(i),kind=4),i=1,npart)
    close(26)

 case(13)
    do i=1,nptmass
       write(*,'(A,I2,A,ES10.3,A,ES10.3)') 'Point mass ',i,': M = ',xyzmh_ptmass(4,i),' and h_soft = ',xyzmh_ptmass(ihsoft,i)
    enddo

    do i=1,10000
       do j=1,4000
          call get_eos_kappa_mesa(X_in,rho_array(i),temp_array(j),kappa_array(i,j),kappat,kappar,ierr)
          !kappa_array(i,j) = kappa
       enddo
    enddo

    unitnum = 1000

    open(unit=unitnum, file='kappa_array.out', status='replace')

    !Write data to file
    do i=1,10000
       write(unitnum,"(4000(3x,es18.11e2,1x))") kappa_array(i,:)
    enddo

    close(unit=unitnum)

 case(14)

    call compute_energies(time)

    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    npart_hist = 0
    nbins = 300
    rad_part = 0.
    dist_part = 0.

    if (dump_number == 0) then
       !allocate(prev_bound_energy(npart))
       !prev_bound_energy = 0.
       allocate(prev_unbound(npart))
       prev_unbound = .false.
    endif
    do i=1,npart
       if (xyzh(4,i)  >=  0) then
          call separation_vector(vxyzu(1:3,i),com_v(1:3),v_part_com)
          call separation_vector(xyzh(1:3,i),com_pos(1:3),r_part_com)
          part_kin = 0.5 * v_part_com(4)**2
          part_ps_pot = 0.

          do k=1,nptmass
             r_sink_part = separation(xyzmh_ptmass(1:3,k),xyzh(1:3,i))**2

             hsoftk = max(xyzmh_ptmass(ihsoft,k),xyzh(4,i))
             hsoft1 = 1.0/hsoftk
             hsoft21= hsoft1**2
             q2i    = r_sink_part*hsoft21
             qi     = sqrt(q2i)

             call kernel_softening(q2i,qi,psoft,fsoft)

             part_ps_pot = part_ps_pot + psoft*hsoft1*xyzmh_ptmass(4,k)
          enddo

          !ray = (/ 1., 0., 0. /)
          !proj_mag = dot_product(r_part_com(1:3),ray)
          !proj = proj_mag * ray
          !proj_orth(1:3) = r_part_com(1:3) - proj(1:3)
          !proj_orth_dist = separation(proj_orth,(/0.,0.,0./))
          !orth_ratio = proj_orth_dist / xyzh(4,i)
          !if (orth_ratio < radkern .and. proj_mag > 0.) then
          !   criteria = .true.
          !else
          !   criteria = .false.
          !endif
          !criteria = .true.

          !bound_energy(i) = ((part_kin + part_ps_pot + vxyzu(4,i)) * particlemass + poten(i)) * unit_energ
          !if ((bound_energy(i) < 0.) .and. (prev_unbound(i) .eqv. .false.)) then
          !if (criteria) then
          rhopart = rhoh(xyzh(4,i), particlemass)
          npart_hist = npart_hist + 1
          rad_part(npart_hist) = separation(xyzh(1:3,i),xyzmh_ptmass(1:3,1))
          !dist_part(npart_hist) = rhopart * unit_density
          call get_eos_pressure_temp_mesa(X_in,rhopart*unit_density,vxyzu(4,i) * unit_ergg,pres,temperature)
          call ionisation_fraction(rhopart * unit_density, temperature, X_in, 1.-X_in-Z_in,xh0,xh1,xhe0,xhe1,xhe2)
          dist_part(1,npart_hist) = xh0
          dist_part(2,npart_hist) = xh1
          dist_part(3,npart_hist) = xhe0
          dist_part(4,npart_hist) = xhe1
          dist_part(5,npart_hist) = xhe2
          !dist_part(npart_hist) = 1.
          !prev_unbound(i) = .true.
          !elseif ((bound_energy(i) < 0.) .and. (prev_unbound(i) .eqv. .true.)) then
          !   prev_unbound(i) = .false.
          !endif
          !prev_bound_energy(i) = bound_energy(i)
       endif
    enddo
    allocate(hist_var(nbins))
    grid_file = (/ '   grid_HIfrac.ev', &
                   '  grid_HIIfrac.ev', &
                   '  grid_HeIfrac.ev', &
                   ' grid_HeIIfrac.ev', &
                   'grid_HeIIIfrac.ev' /)
    do i=1,5
       call avg_histogram_setup(rad_part(1:npart_hist),dist_part(i,1:npart_hist),hist_var,npart_hist,4.,0.5,nbins)

       write(data_formatter, "(a,I5,a)") "(", nbins, "(3x,es18.10e3,1x))"

       if (num == 0) then
          unitnum = 1000

          open(unit=unitnum, file=grid_file(i), status='replace')
          write(unitnum, "(a)") '#why isnt this being displayed?'
          close(unit=unitnum)
       endif

       unitnum=1001+i

       open(unit=unitnum, file=grid_file(i), position='append')

       write(unitnum,data_formatter) hist_var(:)

       close(unit=unitnum)
    enddo

 case(15) !Entropy MESA EoS
    !zeroes the entropy variable and others
    entropy_bound = 0.0
    entropy_unbound = 0.0
    avgtemp_bound = 0.0
    avgtemp_unbound = 0.0
    avgpres_bound = 0.0
    avgpres_unbound = 0.0

    !calling again part of the procedure in case(3) to obtain bound and unbound particles
    !setup
    if (dump_number == 0) then
       print "(4(/,a))", 'Who would want CE energies', &
                            'Must answer me', &
                            'These questions three', &
                            'Ere the new results he see'

       call prompt('1. Was this in a corotating frame?',switch(1),.false.)

       if (switch(1)) then
          call prompt('  1a. Enter value of initial eccentricity : ',e)
          a = xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1)
          sep = a*(1. + e)
          mtot = sum(xyzmh_ptmass(4,:)) + npart*particlemass
          vmag = sqrt(a*(1.-e**2)*mtot)/sep
          omega_corotate = vmag/a
          omega_corotate = Dnint(omega_corotate*10000000.0)/10000000.0
          omega_vec = (/ 0.,0.,omega_corotate /)
       else
          omega_vec = (/ 0.,0.,0. /)
       endif

       call prompt('2. Would you like to use thermal energy in the computation of the bound/unbound status?', switch(2),.false.)

    endif

    encomp(5:) = 0.

    call compute_energies(time)

    ekin = 0.

    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    do i=1,npart
       if (xyzh(4,i)  >=  0) then
          !compute energies
          encomp(ipot_pp) = encomp(ipot_pp) + poten(i)
          encomp(ipot_env) = encomp(ipot_env) + poten(i)

          call separation_vector(xyzh(1:3,i),com_pos(1:3),r_part_com(1:4))

          rhopart = rhoh(xyzh(4,i), particlemass)

          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,&
                                  r_part_com(1),r_part_com(2),r_part_com(3),vxyzu(4,i))

          call cross(omega_vec,r_part_com(1:3),omegacrossr)

          call separation_vector(vxyzu(1:3,i) + omegacrossr(1:3),com_v(1:3),v_part_com(1:4))

          part_kin = 0.5 * particlemass * v_part_com(4)**2
          part_ps_pot = 0.

          call cross(r_part_com(1:3),particlemass * v_part_com(1:3),rcrossmv)

          jz = rcrossmv(3)
          encomp(ijz_tot) = encomp(ijz_tot) + jz

          do k=1,nptmass
             r_sink_part = separation(xyzmh_ptmass(1:3,k),xyzh(1:3,i))**2

             hsoftk = max(xyzmh_ptmass(ihsoft,k),xyzh(4,i))
             hsoft1 = 1.0/hsoftk
             hsoft21= hsoft1**2
             q2i    = r_sink_part*hsoft21
             qi     = sqrt(q2i)

             call kernel_softening(q2i,qi,psoft,fsoft)

             part_ps_pot = part_ps_pot + psoft*hsoft1*xyzmh_ptmass(4,k)*particlemass

             if (k == 1) then
                encomp(ipot_env) = encomp(ipot_env) + psoft*hsoft1*xyzmh_ptmass(4,k)*particlemass
             endif

             if (sqrt(r_sink_part) < 80) then
                inearsink = .true.

             endif

          enddo

          encomp(ipot_ps) = encomp(ipot_ps) + part_ps_pot

          part_therm = particlemass * vxyzu(4,i)

          !compute bound/unbound status
          boundparts_1(1:3,i) = xyzh(1:3,i)
          boundparts_1(4:6,i) = vxyzu(1:3,i)

          if (switch(2)) then
             boundparts_1(7,i) = (part_kin + part_ps_pot + poten(i) + part_therm)/particlemass

             if (part_kin + part_ps_pot + poten(i) + part_therm < 0) then
                boundparts_1(8,i) = -1

             else
                boundparts_1(8,i) = 1

             endif

          else
             boundparts_1(7,i) = (part_kin + part_ps_pot + poten(i))/particlemass

             if (part_kin + part_ps_pot + poten(i) < 0) then
                boundparts_1(8,i) = -1

             else
                boundparts_1(8,i) = 1

             endif

          endif

          !gets entropy for the current particle
          call get_eos_various_mesa(X_in,rhopart*unit_density,vxyzu(4,i) * unit_ergg, &
                                    pres_1(i),temp_1(i),proint_1(i),peint_1(i),troint_1(i), &
                                    teint_1(i),entrop_1(i),abad_1(i),gamma1_1(i),gam_1(i))

          !sums entropy and other quantities for bound particles and unbound particles
          if (boundparts_1(8,i) < 0.0) then !bound
             entropy_bound = entropy_bound + entrop_1(i)
             avgtemp_bound = avgtemp_bound + temp_1(i)
             avgpres_bound = avgpres_bound + pres_1(i)

          else !unbound
             entropy_unbound = entropy_unbound + entrop_1(i)
             avgtemp_unbound = avgtemp_unbound + temp_1(i)
             avgpres_unbound = avgpres_unbound + pres_1(i)

          endif

          entropy_array(1) = entropy_bound
          entropy_array(2) = entropy_unbound
          entropy_array(3) = avgtemp_bound !/ npart !make average
          entropy_array(4) = avgtemp_unbound !/ npart !make average
          entropy_array(5) = avgpres_bound !/ npart !make average
          entropy_array(6) = avgpres_unbound !/ npart !make average

       endif

    enddo

    !writes on file
    columns = (/'       b entr',&
                '     unb entr',&
	        '   avg b temp',&
	        ' avg unb temp',&
	        '   avg b pres',&
	        ' avg unb pres'/)
    call write_time_file('entropy_vs_time', columns, time, entropy_array, 6, dump_number)




    unitnum = unitnum + 1

 end select

 !increase dump number counter
 dump_number = dump_number + 1

 !if ( ANY((/ 2, 3, 4, 8, 11, 12, 13, 14 /) == analysis_to_perform) ) call finish_eos(ieos, ierr)

end subroutine do_analysis

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stellar_profile(time,particlemass,npart,xyzh,vxyzu,profile,ray)
 use eos,          only:ieos,equationofstate,X_in, Z_in
 use eos_mesa,     only:get_eos_kappa_mesa,get_eos_pressure_temp_mesa
 use physcon,      only:kboltz,mass_proton_cgs
 use centreofmass, only:get_centreofmass
 use energies,     only:compute_energies
 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass,rhoh,ihsoft,poten
 use units,        only:udist,unit_ergg,unit_density,unit_pressure,unit_velocity,unit_energ
 use kernel,       only:kernel_softening,radkern

 integer           :: i,j,iprofile,ierr
 integer, parameter:: ncols = 19
 real              :: proj(3),orth(3),proj_mag,orth_dist,orth_ratio
 real              :: rhopart,ponrhoi,spsoundi
 real              :: ideal_temp,temp,kappa,kappat,kappar,pres,ideal_kappa
 real              :: part_kin,part_ps_pot,r_sink_part
 real              :: v_part_com(4),r_part_com(4)
 real              :: com_pos(3),com_v(3)
 real              :: xh0, xh1, xhe0, xhe1, xhe2
 real, intent(in), optional :: ray(3)
 real, intent(in)  :: xyzh(:,:),vxyzu(:,:)
 real, intent(in)  :: particlemass,time
 integer, intent(in) :: npart
 real, intent(out), allocatable :: profile(:,:)
 real              :: temp_profile(ncols,npart)
 logical           :: criteria

 call compute_energies(time)

 call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 iprofile = 0

 do i=1,npart
    if (xyzh(4,i)  >=  0) then

       call separation_vector(xyzh(1:3,i),com_pos(1:3),r_part_com)
       call separation_vector(vxyzu(1:3,i),com_v(1:3),v_part_com)


       if (present(ray)) then
          proj_mag = dot_product(r_part_com(1:3),ray(1:3))
          proj = proj_mag * ray
          orth(1:3) = r_part_com(1:3) - proj(1:3)
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

          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))

          call get_eos_pressure_temp_mesa(X_in,rhopart*unit_density,vxyzu(4,i) * unit_ergg,pres,temp)

          call get_eos_kappa_mesa(X_in,rhopart*unit_density,temp,kappa,kappat,kappar,ierr)

          !kappa = 1.

          part_kin = 0.5 * particlemass * v_part_com(4)**2
          part_ps_pot = 0.
          do j=1,nptmass
             r_sink_part = separation(xyzmh_ptmass(1:3,j),xyzh(1:3,i))**2

             hsoftk = max(xyzmh_ptmass(ihsoft,j),xyzh(4,i))
             hsoft1 = 1.0/hsoftk
             hsoft21= hsoft1**2
             q2i    = r_sink_part*hsoft21
             qi     = sqrt(q2i)

             call kernel_softening(q2i,qi,psoft,fsoft)

             part_ps_pot = part_ps_pot + psoft*hsoft1*xyzmh_ptmass(4,j)*particlemass

          enddo

          call ionisation_fraction(rhopart * unit_density, temp, X_in, 1.-X_in-Z_in,xh0, xh1, xhe0, xhe1, xhe2)

          temp_profile(1,iprofile)  = r_part_com(4) * udist
          temp_profile(3,iprofile)  = atan2(xyzh(2,i),xyzh(1,i))
          temp_profile(4,iprofile)  = vxyzu(4,i) * unit_ergg
          temp_profile(5,iprofile)  = rhopart * unit_density
          temp_profile(6,iprofile)  = ponrhoi * rhopart * unit_pressure
          temp_profile(7,iprofile)  = v_part_com(4) * unit_velocity
          temp_profile(8,iprofile)  = (v_part_com(1) * r_part_com(1) + &
                                  v_part_com(2) * r_part_com(2) + &
                                  v_part_com(3) * r_part_com(3)) / r_part_com(4) * unit_velocity
          temp_profile(9,iprofile)  = spsoundi * unit_velocity
          temp_profile(10,iprofile) = temp
          temp_profile(11,iprofile) = kappa
          temp_profile(12,iprofile) = 1. / (kappa * rhopart * unit_density)
          temp_profile(13,iprofile) = (part_kin + part_ps_pot + poten(i) + (vxyzu(4,i)*particlemass)) * unit_energ
          temp_profile(14,iprofile) = xh1
          temp_profile(15,iprofile) = xhe1
          temp_profile(16,iprofile) = xhe2
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

subroutine ionisation_fraction(dens,temp,X,Y,xh0, xh1, xhe0, xhe1, xhe2)
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

subroutine histogram_setup(distribution_variable,npart,bin_res,bin_size,histogram_bins,allow_bin_res)
 integer, intent(in)          :: npart
 real, dimension(:,:), allocatable, intent(out) :: histogram_bins
 real, dimension(4,npart)     :: distribution_variable
 real                         :: max_value, min_value
 real, intent(out)            :: bin_size
 integer, intent(inout)       :: bin_res
 integer                      :: i
 character(len=30)            :: output_format, val3
 logical                      :: allow_bin_res

 if (allow_bin_res .eqv. .true.) then
    print*, 'Write the number of bins to set resolution (default=1000):'
    read(*,FMT='(a)') val3
    if (val3 == '') then
       bin_res = 1000
    else
       read(val3,*) bin_res
    endif
 elseif ((allow_bin_res .eqv. .false.) .and. (dump_number == 0)) then
    bin_res = 1000
 endif

 allocate(histogram_bins(2,bin_res))

 histogram_bins = 0
 max_value = maxval(distribution_variable(4,:))
 min_value = minval(distribution_variable(4,:))
 bin_size = (max_value - min_value) / bin_res
 output_format = "(5x, a10, a2, F20.10)"
 write(*, output_format) 'Max value',':', max_value
 write(*, output_format) 'Min value',':', min_value
 write(*, output_format) 'Bin size',':', bin_size

 do i=1,bin_res
    histogram_bins(1,i) = min_value + i * bin_size
 enddo

 return
end subroutine histogram_setup

subroutine avg_histogram_setup(dist_var,avg_var,hist_var,npart,max_value,min_value,nbins)
 integer, intent(in)                :: npart,nbins
 real, dimension(npart), intent(in) :: dist_var, avg_var
 real, dimension(nbins), intent(out) :: hist_var
 real, intent(in)                   :: max_value, min_value
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
       hist_var(i) = hist_var(i) / real(n)
    endif
 enddo

end subroutine avg_histogram_setup

subroutine write_file(name_in, dir_in, cols, data_in, npart, ncols, num)
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

 write(file_name, "(2a,i5.5,a)") name_in, "_", num, ".ev"

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

subroutine separation_vector(a,b,c)
 real, intent(in), dimension(3) :: a,b
 real, intent(out), dimension(4) :: c

 c(1) = a(1) - b(1)
 c(2) = a(2) - b(2)
 c(3) = a(3) - b(3)
 c(4) = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
end subroutine separation_vector

real function separation(a,b)
 real, intent(in), dimension(3) :: a,b

 separation = sqrt((a(1)-b(1))**2 +&
                   (a(2)-b(2))**2 +&
                   (a(3)-b(3))**2)
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
