!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
! :Dependencies: eos, io, physcon, setcubiccore, setfixedentropycore,
!   table_utils
!
 implicit none
 real :: rcore,mcore

 ! rcore: Radius / Rsun below which we replace the original profile with a
 !        softened profile.
 ! mcore: Mass / Msun of core particle

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that sets a softened core profile
!+
!-----------------------------------------------------------------------
subroutine set_softened_core(isoftcore,isofteningopt,r,den,pres,m,X,Y,ierr)
 use eos,         only:ieos,X_in,Z_in,init_eos,get_mean_molecular_weight
 use io,          only:fatal
 use table_utils, only:interpolator,yinterp,flip_array
 use setcubiccore,only:set_cubic_core,find_mcore_given_rcore,find_rcore_given_mcore,check_rcore_and_mcore
 use setfixedentropycore,only:set_fixedS_softened_core
 use physcon, only:solarr,solarm
 integer, intent(in) :: isoftcore,isofteningopt
 real, intent(inout) :: r(:),den(:),m(:),pres(:),X(:),Y(:)
 integer             :: core_index,ierr
 real                :: Xcore,Zcore,rc
 logical             :: isort_decreasing,iexclude_core_mass

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

 write(*,'(1x,a,f7.5,a,f7.5,a,f7.5)') 'Using composition at core boundary: X = ',Xcore,', Z = ',Zcore,&
                                      ', mu = ',get_mean_molecular_weight(Xcore,Zcore)
 call interpolator(r,rc,core_index)  ! find index of core
 X(1:core_index) = Xcore
 Y(1:core_index) = yinterp(Y,r,rc)
 if (ieos==10) then
    X_in = Xcore
    Z_in = Zcore
 endif

 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('set_softened_core','could not initialise equation of state')

 ! call core-softening subroutines
 select case(isoftcore) ! choose type of core-softening
 case(1)
    call set_cubic_core(mcore,rcore,den,r,pres,m)
 case(2)
    call set_fixedS_softened_core(mcore,rcore,den,r,pres,m,Xcore,1.-Xcore-Zcore,ierr)
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

end module setsoftenedcore
