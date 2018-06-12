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
!  Analysis routine for mapping to kepler (bins radially by mass)
!
!  REFERENCES: None
!
!  OWNER: Nicole Rodrigues
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, dim, part, physcon, sortutils, units
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
 use sortutils
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) ! due to reset center of mass
 real,             intent(in)    :: particlemass,time
 integer, parameter :: nbins = 1000                      ! number of bins
 real,    parameter :: mmax = 1.0072                         ! total mass
 real,    parameter :: g = 5./3.                         ! adiabatic index
 logical            :: logr  = .true.
 integer            :: i,j,ii
 integer            :: ibins(nbins)
 real               :: ri,ui,hi,dm,mnew
 real               :: rbins(nbins),mass(nbins),density(nbins),ubins(nbins),u(nbins),temp(nbins),hbins(nbins),m(nbins),mbins(nbins),mpbins(nbins)
 logical            :: iexist
 character(len=200) :: fileout
 real               :: x0(3)
 integer            :: iorder(npart)
 real :: grid(nbins) = (/(i,i=1,nbins,1)/)
 real :: mu(nbins)   = 0.609                       ! mean molecular weight
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
 mbins = 0.0
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 !--Sorting particles
 !
 call set_r2func_origin(x0(1),x0(2),x0(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
 !
 !--Mass of each shell
 dm   = mmax/float(nbins)
 !
 !--Binning data
 ii = 1
 parts: do i = 1,npart
    !
    !--Calculate properties of the particle
    !-- i refers to particle, ii refers to bin
    if (xyzh(4,i)  <  tiny(xyzh)) cycle parts         ! skip dead particles
    j = iorder(i)                                      ! particles in order
    ri = sqrt(dot_product(xyzh(1:3,j),xyzh(1:3,j)))   ! radius of each particle in order
    ui = vxyzu(4,j)                                       ! internal energy of each particle
    hi = xyzh(4,j)                                      ! smoothing length
    mnew  = mass(ii) + particlemass
    !--Adding particles to a bin until they exceed dm
    if (mnew<dm .and. ii<=nbins) then
       ibins(ii) = ibins(ii) + 1
       mass(ii)  = mnew
       rbins(ii) = ri
       ubins(ii) = ubins(ii) + ui
       hbins(ii) = hbins(ii) + hi
       !
       if (ii==1) then
          mbins(ii) = mass(ii)
       else
          mbins(ii) = mbins(ii-1) + mass(ii)
       endif
       !
    else
       ii = ii + 1
    endif
    !
 enddo parts
 !
 !--Write results to file
 write(fileout,'(3a)') 'analysisout_',trim(dumpfile),'.dat'
 fileout=trim(fileout)
 open(iunit,file=fileout)
 write(iunit,"('#',26(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'grid',&
       2,'cell mass',&
       3,'outer mass',&
       4,'radius',&
       5,'density',&
       6,'temperature',&
       7,'int energy',&
       8,'nt1',&
       9,'h1',&
       10,'he4',&
       11,'c12',&
       12,'n14',&
       13,'o16',&
       14,'ne20',&
       15,'mg24',&
       16,'si28',&
       17,'s32',&
       18,'ar36',&
       19,'ca40',&
       20,'ti44',&
       21,'cr48',&
       22,'fe52',&
       23,'fe54',&
       24,'ni56',&
       25,'fe56',&
       26,'fe'
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
 do i = 1,ii
    !--density in each shell
    if (mass(i) == 0.) then
       density(i) = 0.
    else
       density = rhoh(hbins(i)/ibins(i),particlemass)*(umass/udist**3)
    endif
    !--internal energy in each shell
    if (mass(i) == 0.) then
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
    write(iunit,'(26(1pe18.10,1x))') grid(i),mass(i)*umass,mbins(i)*umass,rbins(i)*udist,density(i),&
                                    temp(i),u(i),nt1(i),h1(i),he4(i),c12(i),n14(i),o16(i),ne20(i),mg24(i),si28(i),&
                                    s32(i),ar36(i),ca40(i),ti44(i),cr48(i),fe52(i),fe54(i),ni56(i),fe56(i),fe(i)
 enddo
 close(iunit)
 !
end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
