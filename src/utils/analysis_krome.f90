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
 use raytracer,  only: get_all_tau
 implicit none
 character(len=20), parameter, public :: analysistype = 'krome'
 public :: do_analysis
 logical, allocatable :: mask(:)

 private
 real  :: H_init, He_init, C_init, N_init, O_init
 real  :: S_init, Fe_init, Si_init, Mg_init
 real  :: Na_init, P_init, F_init
 real, allocatable    :: abundance(:,:), abundance_prev(:,:), one(:)
 character(len=16)    :: abundance_label(krome_nmols)
 integer, allocatable :: iorig_old(:)
 integer, allocatable :: iprev(:)
 logical :: done_init = .false.

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only: isdead_or_accreted, iorig, rhoh, nptmass, xyzmh_ptmass, iReff
 use linklist,   only: set_linklist
 use units,      only: utime,unit_density,udist
 use eos,        only: get_temperature, ieos, gamma,gmw, init_eos
 use io,         only: fatal
 use krome_main, only: krome_init, krome
 use krome_user, only: krome_get_names,krome_consistent_x
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character     :: filename
 real, save    :: tprev = 0.
 integer, save :: nprev = 0
 real          :: dt_cgs, rho_cgs, T_gas, gammai, mui, AUV, xi
 real          :: abundance_part(krome_nmols), column_density(npart), xyzh_copy(4,npart)
 integer       :: i, j, ierr, completed_iterations, npart_copy = 0

 if (.not.done_init) then
    done_init = .true.
    print*, "initialising KROME"
    call krome_init()
    print*, "Initialised KROME"
    abundance_label(:) = krome_get_names()
    allocate(abundance(krome_nmols,maxp))
    abundance = 0.
    allocate(abundance_prev(krome_nmols,maxp))
    abundance_prev = 0.
    allocate(one(maxp))
    one = 1.
    allocate(iorig_old(maxp))
    iorig_old = 0
    allocate(iprev(maxp))
    iprev = 0
    print*, "setting abundances"
   !$omp parallel do default(none) &
   !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev) &
   !$omp shared(abundance,abundance_prev,particlemass,unit_density) &
   !$omp shared(ieos,rho_cgs,T_gas,j) &
   !$omp private(i,abundance_part)
    do i=1, npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          call chem_init(abundance_part)
          abundance(:,i) = abundance_part
       endif
    enddo
    print*, "abundances set"
    call init_eos(ieos, ierr)
    if (ierr /= 0) call fatal(analysistype, "Failed to initialise EOS")
 else
    dt_cgs = (time - tprev)*utime
    completed_iterations = 0
    print*, "not first step data, timestep = ",dt_cgs, "npart = ",npart, "nprev = ",nprev
    xyzmh_ptmass(iReff,1) = 2.
    print*, "xyzmh_ptmass(iReff,1) = ", xyzmh_ptmass(iReff,1)
    npart_copy = npart
    xyzh_copy = xyzh(:,:npart)
    call set_linklist(npart_copy,npart_copy,xyzh_copy,vxyzu)
    call get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, one, 5, .false., column_density)
    print*, "got column density"
    print*, column_density(1:5)
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev) &
    !$omp shared(abundance,abundance_prev,particlemass,unit_density,udist) &
    !$omp shared(ieos,gamma,gmw,time,completed_iterations,column_density) &
    !$omp private(i,j,abundance_part,rho_cgs,T_gas,gammai,mui,AUV,xi)
    outer: do i=1,npart
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
             !Thermodynamic quantities
             rho_cgs = rhoh(xyzh(4,i),particlemass)*unit_density
             gammai = gamma
             mui    = gmw
             T_gas = get_temperature(ieos,xyzh(1:3, i),rhoh(xyzh(4,i),particlemass),vxyzu(:,i),gammai,mui)
             T_gas = max(T_gas,20.0d0)
             print*, 'T_gas = ', T_gas, ' rho_cgs = ', rho_cgs
             
             !Radiation quantities
             AUV = 4.65 * column_density(i) / 1.87e21
             xi = get_xi(AUV)
             print*, 'AUV = ', AUV, ' xi = ', xi
             
             call krome_consistent_x(abundance_part)
             call krome(abundance_part,rho_cgs,T_gas,dt_cgs)
             abundance(:,i) = abundance_part
       endif
       !$omp atomic
       completed_iterations = completed_iterations + 1
       print*, 'Completed ', completed_iterations, ' of ', npart
    enddo outer
    stop
 endif

 call write_chem(npart, dumpfile)
 nprev = npart
 tprev = time
 iorig_old(1:npart) = iorig(1:npart)
 abundance_prev(:,1:npart) = abundance(:,1:npart)
end subroutine do_analysis

real function get_xi(AUV)
 use physcon, only: pi
 real, intent(in) :: AUV
 real :: xi
 real :: W(6), GA(6), ceta
 integer :: i

 W(1) = 0.17132449
 W(2) = 0.36076157
 W(3) = 0.46791393
 W(4) = W(1)
 W(5) = W(2)
 W(6) = W(3)
 GA(1) = 0.93246951
 GA(2) = 0.66120939
 GA(3) = 0.23861919
 GA(4) = -GA(1)
 GA(5) = -GA(2)
 GA(6) = -GA(3)

 xi = 0.0
 do i=1,6
   ceta = (pi*GA(i)+pi)/2.0
   xi=xi+(W(i)*(sin(ceta)*exp((-AUV*ceta)/sin(ceta))))
 enddo
 xi = (pi/4.0)*xi
 
 get_xi = xi

end function get_xi

subroutine write_chem(npart, dumpfile)
 use krome_user, only: krome_idx_He,krome_idx_C,krome_idx_N,krome_idx_O,&
       krome_idx_H,krome_idx_S,krome_idx_Fe,krome_idx_Si,krome_idx_Mg,&
       krome_idx_Na,krome_idx_P,krome_idx_F
 integer, intent(in)          :: npart
 character(len=*), intent(in) :: dumpfile
 integer :: i, iu

 open(newunit=iu, file=dumpfile//'.comp', status='replace', action='write')
 write(iu, *) '# H, He, C, N, O, S, Fe, Si, Mg, Na, P, F'
 do i=1, npart
    write(iu, *) abundance(krome_idx_H, i),  abundance(krome_idx_He, i), abundance(krome_idx_C, i),  &
                 abundance(krome_idx_N, i),  abundance(krome_idx_O, i),  abundance(krome_idx_S, i),  &
                 abundance(krome_idx_Fe, i), abundance(krome_idx_Si, i), abundance(krome_idx_Mg, i), &
                 abundance(krome_idx_Na, i), abundance(krome_idx_P, i),  abundance(krome_idx_F, i)
 enddo
 close(iu)
 
end subroutine write_chem

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

end subroutine chem_init

end module analysis