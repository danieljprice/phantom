!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setsoftenedcore
!
! This module solves for a softened density and pressure profile below
! a softening radius, rcore, of a stellar profile
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: eos, eos_mesa, io, physcon, setcubiccore,
!   setfixedentropycore, table_utils
!
 implicit none
 ! rcore: Radius / Rsun below which we replace the original profile with a
 !        softened profile.
 ! mcore: Mass / Msun of core particle

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that sets a softened core profile
!+
!-----------------------------------------------------------------------
subroutine set_softened_core(eos_type,isoftcore,isofteningopt,regrid_core,rcore,mcore,r,den,pres,m,X,Y,ierr)
 use eos,                 only:X_in,Z_in,init_eos,gmw,get_mean_molecular_weight
 use eos_mesa,            only:init_eos_mesa
 use io,                  only:fatal
 use table_utils,         only:interpolator,yinterp,flip_array
 use setcubiccore,        only:set_cubic_core,find_mcore_given_rcore,&
                               find_rcore_given_mcore,check_rcore_and_mcore
 use setfixedentropycore, only:set_fixedS_softened_core
 use physcon,             only:solarr,solarm
 integer, intent(in) :: eos_type,isoftcore,isofteningopt
 logical, intent(in) :: regrid_core
 real, intent(inout) :: rcore,mcore
 real, intent(inout), allocatable :: r(:),den(:),m(:),pres(:),X(:),Y(:)
 integer             :: core_index,ierr,npts,Ncore
 real                :: Xcore,Zcore,rc
 logical             :: isort_decreasing,iexclude_core_mass
 real, allocatable   :: r1(:),den1(:),pres1(:),m1(:),X1(:),Y1(:)

 ierr = 0
 write(*,'(/,1x,a)') 'Setting softened core profile'
 ! Output data to be sorted from stellar surface to interior?
 isort_decreasing = .true.     ! Needs to be true if to be read by Phantom
 !
 ! Exclude core mass in output mass coordinate?
 iexclude_core_mass = .true.   ! Needs to be true if to be read by Phantom

 ! get values of rcore and mcore
 if (isoftcore == 1) then
    select case (isofteningopt)
    case(1)
       call find_mcore_given_rcore(rcore,r,den,m,mcore,ierr)
       if (ierr==1) call fatal('setup','Cannot find mcore that produces nice profile (mcore/m(h) > 0.98 reached)')
    case(2)
       call find_rcore_given_mcore(mcore,r,den,m,rcore,ierr)
       if (ierr==1) call fatal('setup','Cannot find softening length that produces nice profile (h/r(mcore) < 1.02 reached)')
    case(3) ! Both rcore and mcore are specified, check if values are sensible
       call check_rcore_and_mcore(rcore,mcore,r,den,m,ierr)
       if (ierr==1) call fatal('setup','mcore cannot exceed m(r=h)')
       if (ierr==2) call fatal('setup','softenedrho/rho > tolerance')
       if (ierr==3) call fatal('setup','drho/dr > 0 found in softened profile')
    end select
 endif

 ! set fixed composition (X, Z, mu) of softened region (take on values at the core boundary)
 rc = rcore*solarr
 Xcore = yinterp(X,r,rc)
 Zcore = 1.-Xcore-yinterp(Y,r,rc)
 gmw = get_mean_molecular_weight(Xcore,Zcore)

 write(*,'(1x,3(a,f7.5))') 'Using composition at core boundary: X = ',Xcore,', Z = ',Zcore,', mu = ',gmw
 call interpolator(r,rc,core_index)  ! find index of core
 X(1:core_index) = Xcore
 Y(1:core_index) = yinterp(Y,r,rc)
 if (eos_type==10) then
    X_in = Xcore
    Z_in = Zcore
    if (ierr /= 0) call fatal('set_softened_core','could not initialise equation of state')
 endif
 call init_eos(eos_type,ierr)  ! need to initialise EoS again with newfound composition (also needed for iopacity_type = 1)

 if (regrid_core) then
    ! make copy of original arrays
    npts = size(r)
    allocate(r1(npts),den1(npts),pres1(npts),m1(npts),X1(npts),Y1(npts))
    r1 = r
    den1 = den
    pres1 = pres
    m1 = m
    X1 = X
    Y1 = Y
    Ncore = 5000  ! number of grid points in softened region (hardwired for now)
    call calc_regrid_core(Ncore,rc,core_index,r1,den1,pres1,m1,X1,Y1,r,den,pres,m,X,Y)
    X(:) = X(size(X))
    Y(:) = Y(size(Y))
 endif

 ! call core-softening subroutines
 select case(isoftcore) ! choose type of core-softening
 case(1)
    call set_cubic_core(mcore,rcore,den,r,pres,m)
 case(2)
    call set_fixedS_softened_core(eos_type,mcore,rcore,den,r,pres,m,Xcore,1.-Xcore-Zcore,ierr)
    if (ierr /= 0) call fatal('setup','could not set fixed entropy softened core')
 end select

 ! Reverse arrays so that data is sorted from stellar surface to stellar centre.
 if (isort_decreasing) then
    call flip_array(m)
    call flip_array(pres)
    call flip_array(r)
    call flip_array(den)
    call flip_array(X)
    call flip_array(Y)
 endif

 if (iexclude_core_mass) then
    m = m - mcore*solarm
 endif

end subroutine set_softened_core


!-----------------------------------------------------------------------
!+
!  Increase number of grid points in softened region to help converge in
!  shooting. Currently, use linear grid points
!
!  Ncore: No. of grid points to use in softened region
!+
!-----------------------------------------------------------------------
subroutine calc_regrid_core(Ncore,rcore_cm,icore,r1,den1,pres1,m1,X1,Y1,r2,den2,pres2,m2,X2,Y2)
 integer, intent(in)            :: Ncore
 real, intent(in)               :: rcore_cm
 integer, intent(inout)         :: icore
 real, intent(in), dimension(:) :: r1,den1,pres1,m1,X1,Y1
 real, intent(out), dimension(:), allocatable :: r2,den2,pres2,m2,X2,Y2
 integer                        :: npts,npts_old,i
 real                           :: dr

 npts_old = size(r1)
 npts = npts_old - icore + Ncore

 allocate(r2(npts),den2(npts),pres2(npts),m2(npts),X2(npts),Y2(npts))
 r2(Ncore:npts)    = r1(icore:npts_old)
 den2(Ncore:npts)  = den1(icore:npts_old)
 pres2(Ncore:npts) = pres1(icore:npts_old)
 m2(Ncore:npts)    = m1(icore:npts_old)
 X2(Ncore:npts)    = X1(icore:npts_old)
 Y2(Ncore:npts)    = Y1(icore:npts_old)

 ! Set uniform r grid in softened region
 dr = rcore_cm/real(Ncore)
 do i = 1,Ncore-1
    r2(i) = real(i)*dr
 enddo
 r2(Ncore) = rcore_cm

 icore = Ncore

end subroutine calc_regrid_core

end module setsoftenedcore
