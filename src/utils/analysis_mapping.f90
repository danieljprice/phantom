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
!  Analysis routine for mapping to kepler
!
!  REFERENCES: None
!
!  OWNER: Nicole Rodrigues
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, dim, part, physcon, units
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'kepler'
 logical :: firstcall = .true.
 public  :: do_analysis

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass
 use physcon,      only: pi,gg,years
 use units,        only: umass,udist,utime
 use part,         only: rhoh
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) ! due to reset center of mass
 real,             intent(in)    :: particlemass,time
 integer, parameter :: nbins = 500                       ! number of bins
 real,    parameter :: rmax  =  1.0                      ! radius of star
 real,    parameter :: g = 5./3.
 logical            :: logr  = .true.
 integer            :: i,j,ii
 integer            :: ibins(nbins)
 real               :: dr,total_mass,ri,ui,hi
 real               :: rbins(nbins),mass(nbins),density(nbins),ubins(nbins),u(nbins),temp(nbins),hbins(nbins),np(nbins)
 logical            :: iexist
 character(len=200) :: fileout
 real :: grid(nbins) = (/(i,i=1,nbins,1)/)
 real :: mu(nbins)   = 1/(8.50479435771475778e-01) ! mean molecular weight (1/Ye assuming full ionisation) in cgs
 real :: mh(nbins)   = 1.6737236e-24               ! mass of hydrogen atom in cgs
 real :: kb(nbins)   = 1.380658e-16                ! boltzmann constant in cgs
 real :: nt1(nbins)  = 0.
 real :: h1(nbins)   = 7.00873167255713247e-01
 real :: he4(nbins)  = 2.83724159091303385e-01
 real :: c12(nbins)  = 1.12052240487196221e-05
 real :: n14(nbins)  = 3.69260191170019125e-03
 real :: o16(nbins)  = 6.61614668185568373e-03
 real :: ne20(nbins) = 1.30599353794390869e-03
 real :: mg24(nbins) = 7.91396084172139185e-04
 real :: si28(nbins) = 8.29395896149018839e-04
 real :: s32(nbins)  = 4.23297905521931848e-04
 real :: ar36(nbins) = 1.13099440383961051e-04
 real :: ca40(nbins) = 7.38096347899203195e-05
 real :: ti44(nbins) = 3.82798105914944417e-06
 real :: cr48(nbins) = 3.42898303339186569e-05
 real :: fe52(nbins) = 0.
 real :: fe54(nbins) = 1.46099277100764633e-03
 real :: ni56(nbins) = 0.
 real :: fe56(nbins) = 0.
 real :: fe(nbins)   = 0.

 !
 !-- Initialise parameters
 !
 rbins = 0.0
 ibins = 0
 ubins = 0.0
 hbins = 0.0
 !
 !--Set bins
 dr   = rmax/float(nbins)    ! radius of each shell
 do i = 1,nbins
    rbins(i) = float(i)*dr   ! radius of each shell from centre of star
 enddo
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 !--Loop over particles putting properties into the correct bin/Bin the data
 do i = 1,npart
    !
    !--Calculate properties of the particle
    !-- i refers to particle, ii refers to bin
    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       ri = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))   ! radius of each particle
       ui = vxyzu(4,i)                                       ! internal energy of each particle
       hi = xyzh(4,i)                                       ! smoothing length
       ii = int((ri-rbins(1))/dr + 1)                    ! binning particles by radius
       if (ii > nbins) cycle
       !
       ibins(ii) = ibins(ii) + 1
       ubins(ii) = ubins(ii) + ui
       hbins(ii) = hbins(ii) + hi
    endif
 enddo
 !
 !--Write results to file
 write(fileout,'(3a)') 'analysisout_',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',26(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'grid',&
       2,'cell mass',&
       2,'outer mass',&
       3,'radius',&
       4,'density',&
       5,'temperature',&
       6,'int energy',&
       7,'nt1',&
       8,'h1',&
       9,'he4',&
       10,'c12',&
       11,'n14',&
       12,'o16',&
       13,'ne20',&
       14,'mg24',&
       15,'si28',&
       16,'s32',&
       17,'ar36',&
       18,'ca40',&
       19,'ti44',&
       20,'cr48',&
       21,'fe52',&
       22,'fe54',&
       23,'ni56',&
       24,'fe56',&
       25,'fe'
 write(iunit,"('#',26(1x,'[',2x,1x,a11,']',2x))") &
       'unit',&
       'g',&
       'g',&
       'cm',&
       'g/cm^3',&
       'K',&
       'erg/g',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac',&
       'mfrac'
 total_mass = 0.0
 do i = 1,nbins
    mass = ibins(i)*particlemass                        ! mass in each shell
    total_mass = total_mass + ibins(i)*particlemass     ! cumulative sum of mass in each shell
    !--density in each shell
    if (mass(i) == 0.) then
       density(i) = 0.
    else
       !density = (1./(((4.*pi)/3.)*(rbins(i)**3-rbins(i-1)**3)))*(ibins(i)*particlemass)
       density = rhoh(hbins(i)/ibins(i),particlemass)*(umass/udist**3)
    endif
    !--internal energy in each shell
    if (ibins(i) == 0.) then
       u = 0.
    else
       u = (ubins(i)/ibins(i))*(1.907e15)
    endif
    !--temperature in each shell
    if (mass(i) == 0.) then
       temp = 0.
    else
       temp = (mu*mh/kb)*(g-1.)*u
    endif
    write(iunit,'(26(1pe18.10,1x))') grid(i),mass(i)*umass,total_mass*umass,rbins(i),density(i),&
                                    temp(i),u(i),nt1(i),h1(i),he4(i),c12(i),n14(i),o16(i),ne20(i),mg24(i),si28(i),&
                                    s32(i),ar36(i),ca40(i),ti44(i),cr48(i),fe52(i),fe54(i),ni56(i),fe56(i),fe(i)
 enddo
 close(iunit)
 !
end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module


