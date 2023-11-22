!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module testeos_stratified
!
! Unit tests of the stratified equation of state
!
! :References: None
!
! :Owner: Caitlyn Hardiman
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, physcon, testutils, units
!

 use testutils,     only:checkval,update_test_scores,checkvalbuf,checkvalbuf_end
 !use units,         only:set_units,unit_density

 implicit none
 public :: test_eos_stratified

 integer, parameter :: n = 5, nr = 700, nz = 210
 ! Parameters are found using the fits from Law et al. 2021
 ! Disc order: HD 1632996, IM Lup, GM Aur, AS 209, MWC 480
 real, parameter :: qfacdiscs(n) = (/0.09,0.01,0.005,0.09,0.115/)
 real, parameter :: qfacdisc2s(n) = (/0.305,-0.015,0.275,0.295,0.35/)
 real, parameter :: alpha_zs(n) = (/3.01,4.91,2.57,3.31,2.78/)
 real, parameter :: beta_zs(n) = (/0.42,2.07,0.54,0.02,-0.05/)
 real, parameter :: z0s(n) = (/1.30089579367134,2.1733078802249720E-004,1.0812929024334721, &
          4.5600541967795483,8.8124778825591701/)
 real, parameter :: polyks(n) = (/2./3.*3.222911812370378E-004,2./3.*1.6068568523984949E-004, &
             2./3.*1.2276291046706421E-004, 2./3.*3.3571998045524743E-004, &
             2./3.*4.5645812781352422E-004/)
 real, parameter :: polyk2s(n) = (/4.0858881228052306E-003,1.2253168963394993E-004, &
              2.3614956983147709E-003,2.1885055156599335E-003, &
              6.7732173498498277E-003/)
 real, parameter :: temp_mid0s(n) = (/24,25,20,25,27/)
 real, parameter :: temp_atm0s(n) = (/63,36,48,37,69/)
 real, parameter :: z0_originals(n) = (/9,3,13,5,7/)
 real, parameter :: q_mids(n) = (/-0.18,-0.02,-0.01,-0.18,-0.23/)
 real, parameter :: q_atms(n) = (/-0.61,0.03,-0.55,-0.59,-0.7/)


 private

contains
!----------------------------------------------------------
!+
!  unit tests of stratified equation of state
!+
!----------------------------------------------------------
subroutine test_eos_stratified(ntests,npass)
 use io,        only:master,stdout
 use physcon,   only:solarm,au
 use units,     only:set_units
 integer, intent(inout) :: ntests,npass


 call test_stratified_midplane(ntests,npass)
 call test_stratified_temps(ntests,npass)
 call test_stratified_temps_dartois(ntests,npass)
 !call map_stratified_temps(ntests,npass)


end subroutine test_eos_stratified

!----------------------------------------------------------------------------
!+
!  test ieos=7 has matches ieos=3 in midplane
!+
!----------------------------------------------------------------------------
subroutine test_stratified_midplane(ntests, npass)
 use eos,           only:maxeos,equationofstate,eosinfo,init_eos,qfacdisc, &
                         qfacdisc2,z0,alpha_z,beta_z,polyk,polyk2,istrat
 use units,         only:unit_density
 use io,            only:id,master,stdout
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2),ncheck(2)
 integer :: ierr,ieos,i,j,k
 real    :: rhoi,tempi,xi,yi,zi,ponrhoi,spsoundi,tempi_ref,temp_mid0, &
            temp_atm0,z0_original,q_atm,q_mid,spsoundi_ref
 real    :: errmax


 if (id==master) write(*,"(/,a)") '--> testing stratified disc equation of state'

 ieos = 7
 istrat = 1

 ! ieos=7, stratified eos, requires polyk to be set to avoid undefined
 polyk = 0.1
 polyk2 = 0.1
 tempi = -1
 tempi_ref = -1


 call init_eos(ieos, ierr)
 if (ierr /= 0) then
    write(*,"(/,a)") '--> skipping stratified disc eos test due to init_eos() fail'
    return
 endif

 nfailed = 0
 ncheck  = 0
 errmax = 0

 call eosinfo(ieos,stdout)


 do i=1,5
    call get_disc_params(i,qfacdisc,qfacdisc2,alpha_z,beta_z,z0,polyk,polyk2, &
                         temp_mid0,temp_atm0,z0_original,q_mid,q_atm)
    rhoi = 1e-13/unit_density

    do j=1,1000,10
       xi=j
       do k=1,1
          yi = 0
          zi = 0

          call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi)
          !print*, xi, zi, tempi, spsoundi*unit_velocity, 'cm/s'

          call equationofstate(3,ponrhoi,spsoundi_ref,rhoi,xi,yi,zi,tempi_ref)
          !print*, xi, zi, tempi_ref, spsoundi_ref*unit_velocity, 'cm/s'
          call checkvalbuf(spsoundi,spsoundi_ref,1e-12,'ieos=3 cs in midplane matches ieos=7 cs in midplane', &
                         nfailed(1),ncheck(1),errmax)
          call checkvalbuf(tempi,tempi_ref,1e-12,'ieos=3 temp in midplane matches ieos=7 temp in midplane', &
                         nfailed(2),ncheck(2),errmax)

       enddo
    enddo
    call checkvalbuf_end('ieos=3 cs in midplane matches ieos=7 cs in midplane',ncheck(1),nfailed(1),errmax,1e-12)
    call checkvalbuf_end('ieos=3 temp in midplane matches ieos=7 temp in midplane',ncheck(2),nfailed(2),errmax,1e-12)
    call update_test_scores(ntests,nfailed,npass)
 enddo

end subroutine test_stratified_midplane

!----------------------------------------------------------------------------
!+
!  test ieos=7 produces temperature equal to that from MAPS paper
!+
!----------------------------------------------------------------------------
subroutine test_stratified_temps(ntests, npass)
 use eos,           only:maxeos,equationofstate,eosinfo,init_eos,qfacdisc, &
                         qfacdisc2,z0,alpha_z,beta_z,polyk,polyk2,istrat
 use units,         only:unit_density,set_units
 use physcon,       only:au,solarm
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2),ncheck(2)
 integer :: ieos,ierr,i,j,k,l
 real    :: rhoi,tempi,xi,yi,zi,ponrhoi,spsoundi,temp_ref,temp_mid0, &
            temp_atm0,z0_original,q_atm,q_mid,ri,temp_atm,temp_mid,zq
 real    :: errmax
 integer, parameter :: nstep=20,nmax=1000

 ieos = 7
 istrat = 0

 ! ieos=7, stratified eos, requires polyk to be set to avoid undefined
 polyk = 0.1
 polyk2 = 0.1
 tempi = -1


 call set_units(mass=solarm,dist=au,G=1.d0)

 call init_eos(ieos, ierr)
 if (ierr /= 0) then
    write(*,"(/,a)") '--> skipping stratified disc eos test due to init_eos() fail'
    return
 endif

 nfailed = 0.
 ncheck  = 0.
 errmax = 0.


 do i=1,n
    call get_disc_params(i,qfacdisc,qfacdisc2,alpha_z,beta_z,z0,polyk,polyk2, &
                         temp_mid0,temp_atm0,z0_original,q_mid,q_atm)
    do j=1,nmax,nstep
       xi=j
       do k=1,nmax,nstep
          yi = k
          do l=1,nmax/2,nstep/2
             zi = l
             rhoi = 1e-13/unit_density
             call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi)
             ri = sqrt(xi**2 + yi**2)
             zq = z0_original*(ri/100)**beta_z
             temp_mid = temp_mid0*(ri/100)**q_mid
             temp_atm = temp_atm0*(ri/100)**q_atm
             temp_ref = (temp_mid**4 + 0.5*(1+tanh((abs(zi) - alpha_z*zq)/zq))*temp_atm**4)**(0.25)
             call checkvalbuf(tempi,temp_ref,1e-14,'ieos=7 temp matches temp from Law et al. 2021 equation',&
             nfailed(1),ncheck(1),errmax)
          enddo
       enddo
    enddo
    call checkvalbuf_end('ieos=7 temp matches temp from Law et al. 2021 equation',ncheck(1),nfailed(1),errmax,1e-14)
    call update_test_scores(ntests,nfailed,npass)
 enddo

end subroutine test_stratified_temps

!----------------------------------------------------------------------------
!+
!  test ieos=7 produces temperature equal to that from Dartois paper
!+
!----------------------------------------------------------------------------
subroutine test_stratified_temps_dartois(ntests, npass)
 use eos,           only:maxeos,equationofstate,eosinfo,init_eos,qfacdisc, &
                         qfacdisc2,z0,beta_z,polyk,polyk2,istrat
 use io,            only:master,stdout
 use testutils,     only:checkval,update_test_scores,checkvalbuf,checkvalbuf_end
 use units,         only:unit_density,set_units
 use physcon,       only:au,solarm
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2),ncheck(2)
 integer :: ierr,ieos,j,k,l
 real    :: rhoi,tempi,xi,yi,zi,ponrhoi,spsoundi,temp_ref,temp_mid0, &
            temp_atm0,z0_original,q_atm,q_mid,ri,temp_atm,temp_mid,zq
 real    :: errmax
 integer, parameter :: nstep=20,nmax=1000

 real, parameter :: pi = 4.*atan(1.0)

 ieos = 7

 ! ieos=7, stratified eos, requires polyk to be set to avoid undefined
 polyk = 0.1
 polyk2 = 0.1
 tempi = -1


 call set_units(mass=solarm,dist=au,G=1.d0)

 nfailed = 0.
 ncheck  = 0.
 errmax = 0.

 call init_eos(ieos, ierr)

 qfacdisc = 0.17
 qfacdisc2 = 0.48
 beta_z = 0.07
 z0 = 43.466157604499408
 polyk = 2./3. * 7.7436597566195883E-004
 !polyk = 5.162439837746392E-004
 polyk2 = 2.7824007780848647E-002
 temp_mid0 = 27.6
 temp_atm0 = 85.6
 z0_original = 60
 q_mid = -0.34
 q_atm = -0.96

 rhoi = 1e-13/unit_density

 do j=1,nmax,nstep
    xi=j
    do k=1,nmax,nstep
       yi = k
       do l=1,nmax/2,nstep/2
          zi = l
          istrat = 1
          call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi)
          !print*, 'eos temp', tempi
          ri = sqrt(xi**2 + yi**2)
          zq = z0_original*(ri/100)**beta_z
          temp_mid = temp_mid0*(ri/100)**q_mid
          temp_atm = temp_atm0*(ri/100)**q_atm
          if (zi < zq) then
             temp_ref = temp_atm + (temp_mid - temp_atm)*(cos((pi/2)*(zi/zq)))**2
          else
             temp_ref = temp_atm
          endif
          !print*, 'ref temp', temp_ref
          call checkvalbuf(tempi,temp_ref,1e-14,'ieos=7 temp matches temp from Dartois et al. 2003 equation', &
                          nfailed(1),ncheck(1),errmax)
       enddo
    enddo
 enddo
 call checkvalbuf_end('ieos=7 temp matches temp from Dartois et al. 2003 equation',ncheck(1), &
                         nfailed(1),errmax,1e-14)
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_stratified_temps_dartois

!----------------------------------------------------------------------------
!+
!  print out temperatures for python mapping
!+
!----------------------------------------------------------------------------
subroutine map_stratified_temps(ntests, npass)
 use eos,           only:maxeos,equationofstate,eosinfo,init_eos,qfacdisc, &
                         qfacdisc2,z0,alpha_z,beta_z,polyk,polyk2
 use units,         only:unit_density
 use io,            only:id,master,stdout
 integer, intent(inout) :: ntests,npass
 integer :: ieos,i,j,k,count
 real    :: rhoi,tempi,xi,yi,zi,ponrhoi,spsoundi,temp_mid0, &
            temp_atm0,z0_original,q_atm,q_mid
 real, dimension(nr) :: radius


 if (id==master) write(*,"(/,a)") '--> testing stratified disc temperatures'

 ieos = 7

 call eosinfo(ieos,stdout)

 open(1, file='HD1632996_temps.txt', status = 'replace')
 open(2, file='IMLup_temps.txt', status = 'replace')
 open(3, file='GMAur_temps.txt', status = 'replace')
 open(4, file='AS209_temps.txt', status = 'replace')
 open(5, file='MWC480_temps.txt', status = 'replace')


 do i=1,n
    call get_disc_params(i,qfacdisc,qfacdisc2,alpha_z,beta_z,z0,polyk,polyk2, &
                          temp_mid0,temp_atm0,z0_original,q_mid,q_atm)
    rhoi = 1e-13/unit_density
    do j=0,210
       zi=j
       count = 0
       do k=1,nr
          xi = k
          yi = 0
          tempi = -1
          call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi)
          radius(k) = tempi
       enddo
       write(i,*) radius
    enddo
 enddo

 close(1)
 close(2)
 close(3)
 close(4)
 close(5)

end subroutine map_stratified_temps

subroutine get_disc_params(ndisc,qfacdisc,qfacdisc2,alpha_z,beta_z,z0,polyk, &
                            polyk2,temp_mid0,temp_atm0,z0_original,q_mid,q_atm)
 integer, intent(in) :: ndisc
 real, intent(out) :: qfacdisc,qfacdisc2,alpha_z,beta_z,z0,polyk,polyk2, &
                        temp_mid0,temp_atm0,z0_original,q_mid,q_atm

 qfacdisc = qfacdiscs(ndisc)
 qfacdisc2 = qfacdisc2s(ndisc)
 alpha_z = alpha_zs(ndisc)
 beta_z = beta_zs(ndisc)
 z0 = z0s(ndisc)
 polyk = polyks(ndisc)
 polyk2 = polyk2s(ndisc)
 temp_mid0 = temp_mid0s(ndisc)
 temp_atm0 = temp_atm0s(ndisc)
 z0_original = z0_originals(ndisc)
 q_mid = q_mids(ndisc)
 q_atm = q_atms(ndisc)

end subroutine get_disc_params

end module testeos_stratified
