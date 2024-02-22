!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine comparing time between dumps
!
! :References: None
!
! :Owner: Esseldeurs Mats
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 use krome_user, only: krome_nmols
 use part,       only: maxp
 implicit none
 character(len=20), parameter, public :: analysistype = 'krome'
 public :: do_analysis
 logical, allocatable :: mask(:)

 private
 real  :: H_init, He_init, C_init, N_init, O_init
 real  :: S_init, Fe_init, Si_init, Mg_init
 real  :: Na_init, P_init, F_init
 real, allocatable    :: abundance(:,:), abundance_prev(:,:)
 character(len=16)    :: abundance_label(krome_nmols)
 integer, allocatable :: iorig_old(:)

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only: isdead_or_accreted, iorig, rhoh, eos_vars, igamma, imu
 use units,      only: utime,unit_density
 use eos,        only: get_temperature, ieos
 use krome_main, only: krome_init, krome
 use krome_user, only: krome_get_names,krome_consistent_x
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character     :: filename
 real, save    :: tprev = 0.
 integer, save :: nprev = 0
 real          :: dt_cgs, rho_cgs, T_gas, gamma, mu
 real          :: abundance_part(krome_nmols)
 integer       :: i, j, iprev(npart)

 if (time == 0.) then
    print*, "initialising KROME"
    call krome_init()
    abundance_label(:) = krome_get_names()
    allocate(abundance(krome_nmols,maxp))
    abundance = 0.
    allocate(abundance_prev(krome_nmols,maxp))
    abundance_prev = 0.
    allocate(iorig_old(maxp))
    iorig_old = 0
    print*, "Initialised KROME"
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev) &
    !$omp shared(abundance,abundance_prev,particlemass,unit_density) &
    !$omp shared(eos_vars,ieos,rho_cgs,T_gas,gamma,mu,j) &
    !$omp private(i,abundance_part)
    do i=1, npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          call chem_init(abundance_part)
          abundance(:,i) = abundance_part
       endif
     enddo
 else
    dt_cgs = (time - tprev)*utime
    print*, "not first step data, timestep = ",dt_cgs, "npart = ",npart, "nprev = ",nprev
   !$omp parallel do default(none) &
   !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev) &
   !$omp shared(abundance,abundance_prev,particlemass,unit_density) &
   !$omp shared(eos_vars,ieos) &
   !$omp private(i,j,abundance_part,rho_cgs,T_gas,gamma,mu)
    outer: do i=1,npart
       print*,i
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          inner: do j=1,nprev
             if (iorig(i) == iorig_old(j)) then
                iprev(i) = j
                exit inner
             endif
          enddo inner
          if (j == iprev(i)) then
             abundance_part(:) = abundance_prev(:,iprev(i))
          else
             call chem_init(abundance_part)
          endif
          rho_cgs = rhoh(xyzh(4,i),particlemass)*unit_density
          !TODO: load in temperature from SPH, as gamma and mu are not available
          gamma = 1.2   ! eos_vars(igamma,i)
          mu    = 2.381 ! eos_vars(imu,i)
          T_gas = get_temperature(ieos,xyzh(1:3, i),rhoh(xyzh(4,i),particlemass),vxyzu(:,i),gamma,mu)
         !  print*, "ieos", ieos
         !  print*, "xyzh(1:3,i)",xyzh(1:3,i)
         !  print*, "rhoh(xyzh(4,i),particlemass)",rhoh(xyzh(4,i),particlemass)
         !  print*, "vxyzu(:,i)",vxyzu(:,i)
         !  print*, "gamma",gamma
         !  print*, "mu",mu
         !  print*, "T_gas = ",T_gas
          T_gas = max(T_gas,2000.0d0)
          call krome_consistent_x(abundance_part)
          print*, "calling krome", rho_cgs, T_gas, dt_cgs
          print*, "abundance_part",abundance_part
          call krome(abundance_part,rho_cgs,T_gas,dt_cgs)
          abundance(:,i) = abundance_part
       endif
    enddo outer
 endif

 nprev = npart
 tprev = time
 iorig_old(1:npart) = iorig(1:npart)
 abundance_prev(:,1:npart) = abundance(:,1:npart)
end subroutine do_analysis

subroutine chem_init(abundance_part)
 use krome_user, only: krome_idx_He,krome_idx_C,krome_idx_N,krome_idx_O,&
       krome_idx_H,krome_idx_S,krome_idx_Fe,krome_idx_Si,krome_idx_Mg,&
       krome_idx_Na,krome_idx_P,krome_idx_F
 real, intent(out) :: abundance_part(krome_nmols)

 ! Initial chemical abundance value for AGB surface
 He_init = 3.11e-1 ! mass fraction
 C_init  = 2.63e-3 ! mass fraction
 N_init  = 1.52e-3 ! mass fraction
 O_init  = 9.60e-3 ! mass fraction

 S_init  = 3.97e-4 ! mass fraction
 Fe_init = 1.17e-3 ! mass fraction
 Si_init = 6.54e-4 ! mass fraction
 Mg_init = 5.16e-4 ! mass fraction

 Na_init = 3.38e-5 ! mass fraction
 P_init  = 8.17e-6 ! mass fraction
 F_init  = 4.06e-7 ! mass fraction

 H_init = 1.0 - He_init - C_init - N_init - O_init - S_init - Fe_init - &
         Si_init - Mg_init - Na_init - P_init - F_init

 abundance_part(:)            = 0.
 abundance_part(krome_idx_He) = He_init
 abundance_part(krome_idx_C)  = C_init
 abundance_part(krome_idx_N)  = N_init
 abundance_part(krome_idx_O)  = O_init
 abundance_part(krome_idx_S)  = S_init
 abundance_part(krome_idx_Fe) = Fe_init
 abundance_part(krome_idx_Si) = Si_init
 abundance_part(krome_idx_Mg) = Mg_init
 abundance_part(krome_idx_Na) = Na_init
 abundance_part(krome_idx_P)  = P_init
 abundance_part(krome_idx_F)  = F_init
 abundance_part(krome_idx_H)  = H_init

!  print*, He_init, C_init, N_init, O_init, S_init, Fe_init, Si_init, Mg_init&
!          , Na_init, P_init, F_init, H_init

end subroutine chem_init

end module analysis