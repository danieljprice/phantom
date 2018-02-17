!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!                  1D-SPH Example Programme: Main Code                 !
!                                                                      !
!                 Copyright (c) 2015-2017 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!+
! This is an example 1D SPH programme to outline the general procedure
! of how to implement NICIL into a pre-existing SPH code.
!
! This is the main code, which can is to be used as a template to
! implement NICIL into the user's code
!
! WARNING! This example code is intented to be an example and not used
! for actual computations.
! Comments indicated with '!**' indicate code that is required for NICIL
!+
!----------------------------------------------------------------------!
program nicil_ex_sph
 use nicil,  only:nicil_initialise,nicil_get_ion_n,nicil_get_eta
 use nicil,  only:nicil_translate_error
 use nicil,  only:nimhd_get_jcbcb,nimhd_get_dBdt,nimhd_get_dudt,nimhd_get_dt
 use nicil,  only:nicil_get_vion,nicil_get_halldrift
 use sphsup, only:initialise_particles,calculate_rhoh,calculate_jcurrent
 use sphsup, only:dkernel,fatal,iprint,iprintdat,iprintwarn
 implicit none
 integer, parameter :: nmax            = 10000           ! maximum array size (not the actual number of particles)
 real,    parameter :: kboltz          =  1.38066d-16    ! Boltzmann constant  [erg/K]
 real,    parameter :: mass_proton_cgs =  1.67262158d-24 ! Proton mass [g]
 real,    parameter :: fourpi          = 12.5663706144d0 ! 4pi
 real,    parameter :: cgsmu0          = fourpi          ! Vacuum permeability [cm/s]
 real,    parameter :: meanmolmass     = 2.38095236      ! Mean molecular mass; calculated with default values in NICIL
 integer            :: npart,ixmin,ixmax,n,i,k,p,ierr
 integer            :: c0,c1,count_rate,count_max
 real               :: utime,udist,umass
 real               :: unit_density,unit_charge,unit_Bfield,unit_energy,unit_velocity,unit_eta
 real               :: xyzh(4,nmax),vxyz(3,nmax),rho(nmax),omega(nmax),B(3,nmax),u(nmax),jcurrent(3,nmax)
 real               :: dBdt(3,nmax),dudt(nmax),n_R(4,nmax),n_electronT(nmax)
 real               :: dt,mass,cs,temperature
 real               :: Bxi,Byi,Bzi,Bi,Bi1, Bxk,Byk,Bzk,Bk,Bk1, dx,dy,dz,rik2,rik,r1,dxr1,dyr1,dzr1,qi,qk
 real               :: vxi,vyi,vzi
 real               :: n_rand,r_rand
 real               :: eta(3,nmax)
 !**terms required for NICIL
 real               :: dBnonideal(3),dBnonidealk(3),jcbcbi(3),jcbi(3),jcbcbk(3),jcbk(3)
 real               :: eta_ohmi,eta_halli,eta_ambii, eta_ohmk,eta_hallk,eta_ambik
 real               :: dudtnonideal,dtohmi,dthalli,dtambii
 !**optional term for NICIL
 real               :: vhall(3,nmax),vdrift(3,nmax),vion(3,nmax)
 !
 !--Open files
 open(unit=iprintdat ,file="data/sph_density.dat")
 open(unit=iprintwarn,file="data/sph_warning.log")
 write(iprintdat,"('#',10(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'density',    &
       2,'B_z',        &
       3,'J_y',        &
       4,'dB/dt_ni,z', &
       5,'eta_OR',     &
       6,'eta_HE',     &
       7,'eta_AD',     &
       8,'du/dt_ni',   &
       9,'v_iondft,x', &
      10,'v_hall,y'
 write(iprintwarn,"(a)") "NICIL: SPH TEST: WARNINGS LOG"
 !
 !--Initialise parameters
 !  Units
 udist         = 1.000d16                  ! = 1 code length unit
 umass         = 1.989d33                  ! = 1 M_sun = 1 code mass unit
 utime         = 8.681d10                  ! = 1 code time unit (Chosen such that G=1)
 unit_velocity = udist/utime               ! = 1 code velocity unit
 unit_density  = umass/udist**3            ! = 1 code length-density unit
 unit_energy   = unit_velocity**2          ! = 1 code internal energy unit
 unit_charge   = sqrt(umass*udist/cgsmu0)  ! = 1 code charge unit
 unit_Bfield   = umass/(utime*unit_charge) ! = 1 code magentic field unit *may be defined differently in user's code*
 unit_eta      = udist**2/utime            ! = 1 code eta unit
 !  Additional parameters
 cs            = 2.19d4                    ! sound speed (cm/s)
 temperature   = cs**2*(meanmolmass*mass_proton_cgs)/kboltz ! Temperature as calculated from sound speed
 cs            = cs/unit_velocity          ! sound speed (code units)
 !  Zero Arrays
 n_R           = 0.0
 n_electronT   = 0.0
 !
 !--Initialise your code; include initialisation of NICIL here.
 call initialise_particles(nmax,npart,ixmin,ixmax,umass,udist,unit_density,unit_Bfield &
                          ,unit_energy,unit_velocity,xyzh,vxyz,rho,omega,B,u,cs,mass)
 !**Initialise NICIL; abort if there is an error in the setup
 call nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint,iprintwarn)
 if (ierr/=0) call fatal(ierr)
 !
 !--Main Loop (The programme will run for n steps; this will likely be
 !  a time constraint in the user's code)
 do n = 1,3
    !--To time one loop
    call system_clock(c0,count_rate,count_max)
    !
    !--Update density prior to main force loop; since n_R & n_electronT are not
    !  dependent on neighbours, they can be calculated here
!$omp parallel default(none) &
!$omp shared(npart,rho,omega,xyzh,mass,B,temperature,n_R,n_electronT) &
!$omp private(i,Bi,ierr)
!$omp do schedule(runtime)
    do i = 1,npart
       call calculate_rhoh(nmax,npart,xyzh,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i),rho(i),omega(i),mass)
       Bi = sqrt( dot_product(B(1:3,i),B(1:3,i)) )
       call nicil_get_ion_n(rho(i),temperature,n_R(1:4,i),n_electronT(i),ierr)
       if (ierr/=0) then
          call nicil_translate_error(ierr)
          call fatal(ierr,rho(i),temperature)
       endif
    enddo
!$omp enddo
!$omp end parallel
    !call calculate_density(nmax,npart,r,rho,mass,h)
    !--Update J prior to main force loop
    call calculate_jcurrent(nmax,npart,mass,xyzh,rho,B,jcurrent,omega)
    !--Initialise the relevant values
    dBdt = 0.0
    dudt = 0.0
    dt   = 1.0d30
    !
    !--Update forces in the main loop
!$omp parallel default(none) &
!$omp shared(npart,ixmin,ixmax,xyzh,vxyz,B,rho,omega,jcurrent,cs,mass,temperature) &
!$omp shared(dBdt,dudt,udist,utime,unit_energy,unit_density,unit_Bfield,unit_eta,eta) &
!$omp shared(n_R,n_electronT,vion,vdrift,vhall) &                                                    !for NICIL
!$omp private(i,k,Bxi,Byi,Bzi,Bi,Bi1,Bxk,Byk,Bzk,Bk,Bk1,vxi,vyi,vzi) &
!$omp private(dx,dy,dz,rik2,rik,r1,dxr1,dyr1,dzr1,qi,qk) &
!$omp private(eta_ohmi,eta_halli,eta_ambii,eta_ohmk,eta_hallk,eta_ambik,dtohmi,dthalli,dtambii) &    !for NICIL
!$omp private(jcbcbi,jcbi,jcbk,jcbcbk,dBnonideal,dBnonidealk,dudtnonideal,ierr) &                    !for NICIL
!$omp reduction(min: dt)
!$omp do schedule(runtime)
    do i = ixmin,ixmax  ! should be i = 1,npart, but for speed, only update the middle values
       !
       !--Define v-terms for simplicity
       vxi = vxyz(1,i)
       vyi = vxyz(2,i)
       vzi = vxyz(3,i)
       !--Define B-terms for simplicity
       Bxi = B(1,i)
       Byi = B(2,i)
       Bzi = B(3,i)
       Bi  = sqrt(Bxi*Bxi + Byi*Byi + Bzi*Bzi)
       Bi1 = 1.0/Bi
       !**Calculate coefficients for particle i
       call nicil_get_eta(eta_ohmi,eta_halli,eta_ambii,Bi,rho(i),temperature,n_R(1:4,i),n_electronT(i),ierr)
       !**Calculate diagnostic velocities for particle i
       call nicil_get_vion(eta_ambii,vxi,vyi,vzi,Bxi,Byi,Bzi,jcurrent(1:3,i),vion(1:3,i),ierr,vdrift(1:3,i))
       call nicil_get_halldrift(eta_halli,Bxi,Byi,Bzi,jcurrent(1:3,i),vhall(1:3,i))
       if (ierr/=0) then
          call nicil_translate_error(ierr)
          call fatal(ierr,rho(i),temperature)
       endif
       !**Calculate JxB and JxBxB for particle i
       call nimhd_get_jcbcb(jcbcbi,jcbi,jcurrent(1:3,i),Bxi,Byi,Bzi,Bi1)
       !
       !--Calculate forces dependent on neighbours
       do k = 1,npart
          if (i/=k) then
             !--Calculate normalised particle separations
             dx   = xyzh(1,i) - xyzh(1,k)
             dy   = xyzh(2,i) - xyzh(2,k)
             dz   = xyzh(3,i) - xyzh(3,k)
             rik2 = dx*dx + dy*dy + dz*dz
             if (rik2 < 4.0*(max(xyzh(4,i),xyzh(4,k)))**2) then
                rik  = sqrt(rik2)
                r1   = 1.0/rik
                dxr1 = dx*r1
                dyr1 = dy*r1
                dzr1 = dz*r1
                qi   = rik/xyzh(4,i)
                qk   = rik/xyzh(4,k)
                !--Define B-terms for simplicity
                Bxk = B(1,k)
                Byk = B(2,k)
                Bzk = B(3,k)
                Bk  = sqrt(Bxk*Bxk + Byk*Byk + Bzk*Bzk)
                Bk1 = 1.0/Bk
                !--Calculate acceleration here
                !--Calculate ideal MHD terms here
                !
                !**Calculate the i-term of the conjugate B-operator
                call nimhd_get_dBdt(dBnonideal,eta_ohmi,eta_halli,eta_ambii,jcurrent(1:3,i),jcbi,jcbcbi,dxr1,dyr1,dzr1)
                dBnonideal = dBnonideal*mass/(rho(i)**2*omega(i))*dkernel(qi)/xyzh(4,i)**4
                !
                !**Calculate coefficients for particle k
                call nicil_get_eta(eta_ohmk,eta_hallk,eta_ambik,Bk,rho(k),temperature,n_R(1:4,k),n_electronT(k),ierr)
                if (ierr/=0) then
                   call nicil_translate_error(ierr)
                   call fatal(ierr,rho(k),temperature)
                endif
                !**Calculate JxB and JxBxB for particle k
                call nimhd_get_jcbcb(jcbcbk,jcbk,jcurrent(1:3,k),Bxk,Byk,Bzk,Bk1)
                !**Calculate the k-term of the conjugate B-operator
                call nimhd_get_dBdt(dBnonidealk,eta_ohmk,eta_hallk,eta_ambik,jcurrent(1:3,k),jcbk,jcbcbk,dxr1,dyr1,dzr1)
                dBnonidealk = dBnonidealk*mass/(rho(k)**2*omega(k))*dkernel(qk)/xyzh(4,k)**4
                !**Sum the two terms & add to the dBdt
                dBnonideal = dBnonideal + dBnonidealk
                dBdt(1,i)  = dBdt(1,i)  - dBnonideal(1)
                dBdt(2,i)  = dBdt(2,i)  - dBnonideal(2)
                dBdt(3,i)  = dBdt(3,i)  - dBnonideal(3)
                !
             endif
          endif
       enddo
       !
       !--Calculate additional values independent of neighbours
       !**Update internal energy due to non-ideal MHD
       call nimhd_get_dudt(dudtnonideal,eta_ohmi,eta_ambii,rho(i),jcurrent(1:3,i),B(1:3,i))
       dudt(i) = dudt(i) + dudtnonideal
       !
       !**Determine minimum possible timestep
       call nimhd_get_dt(dtohmi,dthalli,dtambii,xyzh(4,i),eta_ohmi,eta_halli,eta_ambii)
       dt = min(dt,dtohmi,dthalli,dtambii)
       !
       !--Save the non-ideal mhd coefficients so that they can later be written to file outside the parallel loop
       !  under normal circumstances, it will not be necessary to save these values
       eta(1,i) = eta_ohmi
       eta(2,i) = eta_halli
       eta(3,i) = eta_ambii
    enddo
!$omp enddo
!$omp end parallel
    !
    !--Update to current values
!$omp parallel default(none) &
!$omp shared(ixmin,ixmax,B,u,dt,dBdt,dudt,n_R,n_electronT) &
!$omp private(i,n_rand,r_rand)
!$omp do schedule(runtime)
    do i = ixmin,ixmax
       ! Do not modify in this test since B and u are only being affected by non-ideal MHD
       !B(1:3,i) = B(1:3,i) + dt*dBdt(1:3,i)
       !u(i)     = u(i)     + dt*dudt(i)
       ! The following lines are not to be used in production codes,
       ! but are included for more realistic approximations on runtime
       call random_number(r_rand)
       n_rand      = 1.0+(1.0-r_rand)*1.0d-5
       n_R         = n_R*n_rand
       n_electronT = n_electronT*n_rand
    enddo
!$omp enddo
!$omp end parallel
    !
    !--To time the step
    call system_clock(c1,count_rate,count_max)
    write(iprint,'(a,I1,a,F6.3,a)') "NICIL: SPH TEST: Runtime on loop ",n,": ",(c1-c0)/float(count_rate)," seconds"
 enddo
 !
 !--Write the important values to file
 do i = ixmin,ixmax
    write(iprintdat,'(10(1pe18.10,1x))') &
         rho(i)*unit_density,B(3,i)*unit_Bfield,jcurrent(2,i)*unit_Bfield/udist**3, &
         dBdt(3,i)*unit_Bfield/utime,eta(:,i)*unit_eta,dudt(i)*unit_energy/utime, &
         vdrift(1,i)*unit_velocity,vhall(2,i)*unit_velocity
 enddo
 write(iprint,'(a)') "NICIL: SPH TEST: properties vs number density written to data/sph_density.dat"
 !
 !--Print Statements & outputs
 write(iprint,'(a)')        "NICIL: SPH TEST: Complete."
 write(iprint,'(a,I8)'    ) "NICIL: SPH TEST: Total number of particles =  ", npart
 write(iprint,'(a,I8)'    ) "NICIL: SPH TEST: Number of active particles = ", ixmax-ixmin + 1
 write(iprint,'(a, F10.3)') "NICIL: SPH TEST: Temperature (K) =       ", temperature
 write(iprint,'(a, F10.3)') "NICIL: SPH TEST: Total distance (AU) =   ", (xyzh(1,ixmax)-xyzh(1,ixmin))*udist/1.49597871d13
 write(iprint,'(a,Es10.3)') "NICIL: SPH TEST: Particle Mass (M_sun) = ", mass
 write(iprint,'(a,Es10.3)') "NICIL: SPH TEST: Total Mass (M_sun) =    ", npart*mass
 write(iprint,'(a,Es10.3)') "NICIL: SPH TEST: min(dt) (s) =           ", dt*utime
 !
 close(iprintdat)
 close(iprintwarn)
!----------------------------------------------------------------------!
end  program nicil_ex_sph
