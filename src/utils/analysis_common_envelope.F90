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
!  DEPENDENCIES: centreofmass, energies, eos, kernel, part, physcon,
!    prompting, ptmass, setbinary, units
!+
!--------------------------------------------------------------------------

module analysis

 implicit none
 character(len=20), parameter, public :: analysistype = 'common_envelope'
 integer                              :: analysis_to_perform, histogram_number
 integer                              :: dump_number = 0
 character(len=32)                    :: val1, val2
 real                                 :: constraint_max, constraint_min, omega_corotate=0, init_radius
 integer                              :: bin_res
 logical                              :: boundfile = .false.


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
 use kernel,       only:kernel_softening
 use eos,          only:gamma,equationofstate,ieos,init_eos,finish_eos
 use setbinary,   only:Rochelobe_estimate

 !general variables
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character(len=30)            :: fileout
 integer                      :: unitnum,i,j,k,ierr
 real                         :: hsoft = 3.0

 !case 1 variables
 real                         :: separation(4)

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
 real                         :: encomp(20)
 real                         :: sink_vel(4)=0, r_sink_sink
 real                         :: part_kin, part_ps_pot
 real                         :: hsoftk, hsoft1, hsoft21, q2i, qi, psoft, fsoft, ponrhoi, spsoundi
 real                         :: e=0, a, mtot, sep, vmag, v_part_com(4)=0, omega_vec(3)=0, omegacrossr(3)
 real                         :: r_part_com(4), rcrossmv(3), jz, boundparts(8,npart+nptmass)
 logical                      :: corotating

 integer, parameter           :: ie_tot        = 1
 integer, parameter           :: ie_pot        = 2
 integer, parameter           :: ie_kin        = 3
 integer, parameter           :: ie_therm      = 4
 integer, parameter           :: ipot_sink     = 5
 integer, parameter           :: ikin_sink     = 6
 integer, parameter           :: iorb_sink     = 7
 integer, parameter           :: ipot_env      = 8
 integer, parameter           :: ie_env        = 9
 integer, parameter           :: ikin_bound    = 10
 integer, parameter           :: ikin_unbound  = 11
 integer, parameter           :: imass_bound   = 12
 integer, parameter           :: imass_unbound = 13
 integer, parameter           :: ipot_pp       = 14
 integer, parameter           :: ipot_ps       = 15
 integer, parameter           :: ijz_tot       = 16
 integer, parameter           :: ijz_bound     = 17
 integer, parameter           :: ijz_unbound   = 18
 integer, parameter           :: ijz_orb       = 19
 integer, parameter           :: ie_gas        = 20

 !case 5 variables
 real                         :: totvol, rhopart, VER(1), profile(8,npart)!, avgPressure(npart)
 integer                      :: avgRange, SGNorm
 integer, dimension(:), allocatable :: SGCoeffs

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
 real                         :: mass_outside, mass_frac_out, total_mass=0.0


 !case 9 variables
 real                         :: C_force, dtsg, f2i, hi, dtsgi, fextx, fexty, fextz, phii, fonrmaxi, dtphi2i
 real                         :: fxyz_ptmass_thread(4,nptmass), minsep1, minsep2, minsep1i, minsep2i

 !chose analysis type
 if (dump_number==0) then
    print "(9(/,a,/))", ' 1) Sink separation', &
            ' 2) Bound and unbound quantities', &
            ' 3) Energies', &
            ' 4) Profile from centre of mass', &
            ' 5) Volume equivalent radius of star', &
            ' 6) Mass within initial Roche-lobes', &
            ' 7) Histograms ', &
            ' 8) Mass outside radius', &
            ' 9) Output simulation units'

    analysis_to_perform = 1

    call prompt('Choose analysis type ',analysis_to_perform,1,9)

 endif

 call init_eos(ieos,ierr)

 !analysis
 select case(analysis_to_perform)

 case(1) !sink separation

    !computes sink separation
    separation(1:3) = xyzmh_ptmass(1:3,1) - xyzmh_ptmass(1:3,2)
    separation(4) = sqrt(dot_product(separation(1:3),separation(1:3)))

    columns = (/'x separation',&
                   'y separation',&
                   'z separation',&
                   '  separation'/)

    call write_time_file('separation_vs_time', columns, time, separation, 4, num)

 case(2) !bound and unbound quantities
    bound = 0
    !get the position of the COM to compute all the vectorial quantities wrt it
    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    !loop on particles
    do i=1,npart

       !part to sinks pot en
       part_ps_pot = 0.0

       do k=1,nptmass

          r_sink_part = sqrt((xyzmh_ptmass(1,k) - xyzh(1,i))**2.0 + &
                                (xyzmh_ptmass(2,k) - xyzh(2,i))**2.0 + &
                                (xyzmh_ptmass(3,k) - xyzh(3,i))**2.0)

          part_ps_pot = part_ps_pot - ((gg*umass*utime**2/udist**3) * xyzmh_ptmass(4,k)) / r_sink_part

       enddo

       part_ps_pot = particlemass * part_ps_pot

       !kin en
       part_kin = particlemass * 0.5 * ((vxyzu(1,i) - com_v(1))**2.0 + &
                                                    (vxyzu(2,i) - com_v(2))**2.0 + &
                                                    (vxyzu(3,i) - com_v(3))**2.0)
       !therm en
       etherm_current_part = particlemass * vxyzu(4,i)

       !tot en
       etot_current_part = poten(i) + part_ps_pot + part_kin + etherm_current_part

       rhopart = rhoh(xyzh(4,i), particlemass)
       ppart = (gamma - 1.) * rhopart * etherm_current_part

       r_part_com(1:3) = xyzh(1:3,i) - com_pos(1:3)
       r_part_com(4) = sqrt(dot_product(r_part_com(1:3),r_part_com(1:3)))
       v_part_com(1:3) = vxyzu(1:3,i) - com_v(1:3)
       v_part_com(4) = sqrt(dot_product(v_part_com(1:3),v_part_com(1:3)))

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

       call prompt('1. Was this in a corotating frame?',corotating,.false.)

       if (corotating) then
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

       call prompt('2. Would you like to output bound particle files?', boundfile)
    endif

    encomp(5:) = 0.

    call compute_energies(time)

    ekin = 0.

    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    do i=1,npart
       if (xyzh(4,i)  >=  0) then
          encomp(ipot_pp) = encomp(ipot_pp) + poten(i)
          encomp(ipot_env) = encomp(ipot_env) + poten(i)

          r_part_com(1:3) = xyzh(1:3,i) - com_pos(1:3)
          r_part_com(4) = sqrt(dot_product(r_part_com(1:3),r_part_com(1:3)))

          rhopart = rhoh(xyzh(4,i), particlemass)

          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,&
                                  r_part_com(1),r_part_com(2),r_part_com(3),vxyzu(4,i))

          call cross(omega_vec,r_part_com(1:3),omegacrossr)

          v_part_com(1:3) = vxyzu(1:3,i) + omegacrossr(1:3) - com_v(1:3)
          v_part_com(4) = sqrt(dot_product(v_part_com(1:3),v_part_com(1:3)))

          part_kin = 0.5 * particlemass * v_part_com(4)**2
          part_ps_pot = 0.

          call cross(r_part_com(:3),particlemass * v_part_com(:3),rcrossmv)

          jz = rcrossmv(3)
          encomp(ijz_tot) = encomp(ijz_tot) + jz

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

             if (k == 1) then
                encomp(ipot_env) = encomp(ipot_env) + psoft*hsoft1*xyzmh_ptmass(4,k)*particlemass
             endif
          enddo

          encomp(ipot_ps) = encomp(ipot_ps) + part_ps_pot

          if (boundfile) then
             boundparts(1:3,i) = xyzh(1:3,i)
             boundparts(4:6,i) = vxyzu(1:3,i)
          endif

          boundparts(7,i) = (part_kin + part_ps_pot + poten(i))/particlemass

          if (part_kin + part_ps_pot + poten(i) < 0) then
             encomp(ikin_bound) = encomp(ikin_bound) + part_kin
             encomp(imass_bound) = encomp(imass_bound) + particlemass
             encomp(ijz_bound) = encomp(ijz_bound) + jz
             boundparts(8,i) = -1
          else
             encomp(ikin_unbound) = encomp(ikin_unbound) + part_kin
             encomp(imass_unbound) = encomp(imass_unbound) + particlemass
             encomp(ijz_unbound) = encomp(ijz_unbound) + jz
             boundparts(8,i) = 1
          endif

       endif
    enddo

    open(26,file=trim(dumpfile)//".divv",status='replace',form='unformatted')
    write(26) (real(boundparts(7,i),kind=4),i=1,npart)
    write(26) (real(spsoundi,kind=4),i=1,npart)
    close(26)

    do i=1,nptmass
       call cross(omega_vec,xyzmh_ptmass(:3,i),omegacrossr)

       sink_vel(1:3) = vxyz_ptmass(1:3,i) + omegacrossr(1:3) - com_v(1:3)
       sink_vel(4) = sqrt(dot_product(sink_vel(1:3),sink_vel(1:3)))

       r_part_com(1:3) = xyzmh_ptmass(1:3,i) - com_pos(1:3)
       r_part_com(4) = sqrt(dot_product(r_part_com(1:3),r_part_com(1:3)))

       call cross(r_part_com(:3),xyzmh_ptmass(4,i) * sink_vel(:3),rcrossmv)

       jz = rcrossmv(3)
       encomp(ijz_tot) = jz + encomp(ijz_tot)
       encomp(ijz_orb) = jz + encomp(ijz_orb)
       encomp(ikin_sink) = encomp(ikin_sink) + 0.5 * xyzmh_ptmass(4,i) * sink_vel(4)**2
    enddo

    r_sink_sink = sqrt((xyzmh_ptmass(1,1) - xyzmh_ptmass(1,2))**2 + &
                          (xyzmh_ptmass(2,1) - xyzmh_ptmass(2,2))**2 + &
                          (xyzmh_ptmass(3,1) - xyzmh_ptmass(3,2))**2)

    encomp(ipot_sink) = - xyzmh_ptmass(4,1) * xyzmh_ptmass(4,2) / r_sink_sink

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
                   '  gas energy'/)

    call write_time_file('energy', columns, time, encomp, 20, num)

    if (boundfile) then
       columns = (/'     x-coord', '     y-coord', '     z-coord', &
                       '       x-vel', '       y-vel', '       z-vel', &
                       ' spec energy', '  bound bool'/)
       call write_file('boundparts', 'bound', columns, boundparts, npart, 8, num)
    endif

    unitnum = unitnum + 1

 case(4) !Profile from COM (can be used for stellar profile)

    call compute_energies(time)

    call get_centreofmass(com_pos,com_v,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    !mu = ((2. * X_in) + (0.75 * (1. - X_in - Z_in)) + (0.5 * Z_in))**(-1)

    do i=1,npart
       if (xyzh(4,i)  >=  0) then

          rhopart = rhoh(xyzh(4,i), particlemass)

          call equationofstate(ieos,ponrhoi,spsoundi,rhopart,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))

          r_part_com(1) = xyzh(1,i) - com_pos(1)
          r_part_com(2) = xyzh(2,i) - com_pos(2)
          r_part_com(3) = xyzh(3,i) - com_pos(3)
          r_part_com(4) = sqrt(r_part_com(1)**2 + &
                                  r_part_com(2)**2 + &
                                  r_part_com(3)**2)

          v_part_com(1) = vxyzu(1,i) - com_v(1)
          v_part_com(2) = vxyzu(2,i) - com_v(2)
          v_part_com(3) = vxyzu(3,i) - com_v(3)
          v_part_com(4) = sqrt(v_part_com(1)**2 + &
                                  v_part_com(2)**2 + &
                                  v_part_com(3)**2)

          profile(1,i)  = r_part_com(4) * udist
          profile(3,i)  = atan2(xyzh(2,i),xyzh(1,i))
          profile(4,i)  = vxyzu(4,i) * unit_ergg
          profile(5,i)  = rhopart * unit_density
          profile(6,i)  = ponrhoi * rhopart * unit_pressure
          !profile(9,i)  = profile(6,i) / (profile(1,i) * udist)
          profile(7,i)  = v_part_com(4) * unit_velocity
          profile(8,i) = (v_part_com(1) * r_part_com(1) + &
                              v_part_com(2) * r_part_com(2) + &
                              v_part_com(3) * r_part_com(3)) / r_part_com(4) * unit_velocity

       endif
    enddo

    do i=1,npart
       if (xyzh(4,i)  >=  0) then
          rhopart = rhoh(xyzh(4,i), particlemass)
          profile(2,i) = 0
          !profile(10,i) = 0
          do j=1,npart
             if (profile(1,i) > profile(1,j)) then
                profile(2,i) = profile(2,i) + particlemass
             endif
          enddo
          !profile(10,i) = -(gg * profile(2,i) * umass * rhopart * unit_density) / (profile(1,i) * udist)**2.
       endif
    enddo

    !call Sort(profile,npart,10,1)

    avgRange = 10
    SGCoeffs = (/-21,14,39,54,59,54,39,14,-21/)
    SGNorm   = 231

    !do i=1,npart
    !   if (i < avgRange + 1) then
    !      avgPressure(i) = sum(profile(6,1:i+avgRange))/(i+avgRange)
    !   elseif (i > npart - avgRange) then
    !      avgPressure(i) = sum(profile(6,i-avgRange:npart))/(npart-i+avgRange+1)
    !   else
    !      avgPressure(i) = sum(profile(6,i-avgRange:i+avgRange))/(2*avgRange + 1)
    !   endif
    !enddo

    do i=1,npart-1
       !profile(9,i) = (avgPressure(i+1) - avgPressure(i)) / ((profile(1,i+1) - profile(1,i)) * udist)
    enddo

    columns = (/'      radius',&
                   '  mass coord',&
                   '     azimuth',&
                   ' int. energy',&
                   '     density',&
                   '    pressure',&
                   '    velocity',&
                   '   rad. vel.'/)
    !'       dP/dr',&
    !'  -Gmrho/r^2'/)


    call write_file('profile', 'profile', columns, profile, npart, 8, num)

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
    if (dump_number == 0) then
       call prompt('Input the initial radius of the star', init_radius)
    endif
    mass_outside = 0.
    total_mass = xyzmh_ptmass(4,1)
    do i=1,npart
       r_sink_part = sqrt((xyzmh_ptmass(1,1) - xyzh(1,i))**2.0 + &
                             (xyzmh_ptmass(2,1) - xyzh(2,i))**2.0 + &
                             (xyzmh_ptmass(3,1) - xyzh(3,i))**2.0)
       if (r_sink_part > init_radius) then
          mass_outside = mass_outside + particlemass
       endif
       total_mass = total_mass + particlemass
    enddo

    mass_frac_out = mass_outside / total_mass

    unitnum = 1001
    write(fileout,"(2a,i3.3,a)") 'massoutside.ev'

    if (num == 0) then
       open(unit=unitnum, file=fileout, status='replace')

       write(unitnum,"('#',1x,3('[',i2.1,a18,']',2x))") &
               1,'time',            &
               2,'mass outside',    &
               3,'mass frac out'

       close(unit=unitnum)

       unitnum = unitnum + 1
    endif

    open(unit=unitnum, file=fileout, position='append')

    write(unitnum,"(3(4x,es18.11e2,2x))") &
            time,          &
            mass_outside,  &
            mass_frac_out

    close(unit=unitnum)

    unitnum = unitnum + 1

 case(9) !Units
    write(*,"(/,3(a,es10.3,1x),a)") '     Mass: ',umass,    'g       Length: ',udist,  'cm     Time: ',utime,'s'
    write(*,"(3(a,es10.3,1x),a)") '  Density: ',unit_density, 'g/cm^3  Energy: ',unit_energ,'erg    En/m: ',unit_ergg,'erg/g'
    write(*,"(3(a,es10.3,1x),a)") ' Velocity: ',unit_velocity,'cm/s    Bfield: ',unit_Bfield,'G  Pressure: ',&
                                     unit_pressure,'g/cm s^2'
    write(*,"(2(a,es10.3,1x),/)")   '        G: ', gg*umass*utime**2/udist**3,'             c: ',c*utime/udist

 case(10)

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

 end select

 !increase dump number counter
 dump_number = dump_number + 1

 call finish_eos(ieos, ierr)

end subroutine do_analysis

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine histogram_setup(distribution_variable,npart,bin_res,bin_size,histogram_bins,allow_bin_res)
 implicit none
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
 write(data_formatter, "(a,I2,a)") "(", ncols, "(3x,es18.11e2,1x))"

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
 implicit none
 real,intent(in),dimension(3)  :: a,b
 real,intent(out),dimension(3) :: c

 c(1) = a(2)*b(3)-b(2)*a(3)
 c(2) = a(3)*b(1)-b(3)*a(1)
 c(3) = a(1)*b(2)-b(1)*a(2)

end subroutine cross


!Sorting routines

! --------------------------------------------------------------------
! subroutine  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

subroutine Swap(a, b)
 implicit none
 real, intent(inout) :: a, b
 real                :: Temp

 Temp = a
 a    = b
 b    = Temp
end subroutine Swap

! --------------------------------------------------------------------
! subroutine  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

subroutine Sort(x, npart, nrows, isort)
 implicit none
 integer, intent(in)                   :: npart,nrows,isort
 integer                               :: i,j
 integer                               :: Location
 real, dimension(nrows,npart), intent(inout)    :: x

 do i = 1, npart-1                 ! except for the last
    Location = minloc(x(isort,i:npart),DIM=1) + i - 1  ! find min from this to last
    do j=1,nrows
       call Swap(x(j,i), x(j,Location)) ! swap this and the minimum
    enddo
    if (mod(i,10000) == 0) then
       print*, "Sorting particle ", i
    endif
 enddo
end subroutine Sort

end module
