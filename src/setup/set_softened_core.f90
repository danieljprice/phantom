!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setsoftenedcore
!
! This module softens the core of a MESA stellar profile given a softening
! radius
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, setcubiccore, setfixedentropycore, table_utils, physcon
!
 implicit none
 real    :: rcore,mcore,hsoft

 ! rcore: Radius below which we replace the original profile with a
 !        softened profile.
 ! mcore: Mass of core particle

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that sets a softened core profile
!+
!-----------------------------------------------------------------------
subroutine set_softened_core(isoftcore,isofteningopt,r,den,pres,m,ene,temp,X,Y,ierr)
 use eos,         only:ieos,X_in,Z_in,gmw,init_eos
 use io,          only:fatal
 use setcubiccore,only:set_cubic_core,find_mcore_given_hsoft,find_hsoft_given_mcore,check_hsoft_and_mcore
 use setfixedentropycore,only:set_fixedS_softened_core
 integer, intent(in)        :: isoftcore,isofteningopt
 real, intent(inout)        :: r(:),den(:),m(:),pres(:),ene(:),temp(:),X(:),Y(:)
 integer                    :: ierr

 ! set mu, X, Z mass fractions to be their values at R/2
 if ((ieos == 10) .or. (ieos == 12)) then
    call get_composition(r,den,pres,temp,X,Y,X_in,Z_in,gmw)
    write(*,'(1x,a,f7.5,a,f7.5,a,f7.5)') 'Using composition at R/2: X = ',X_in,', Z = ',Z_in,', mu = ',gmw
 endif
 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('set_softened_core','could not initialise equation of state')

 ! get values of rcore and mcore
 select case (isofteningopt)
 case(1)
    call find_mcore_given_hsoft(rcore,r,den,m,mcore,ierr)
    if (ierr==1) call fatal('setup','Cannot find mcore that produces nice profile (mcore/m(h) > 0.98 reached)')
 case(2)
    call find_hsoft_given_mcore(mcore,r,den,m,rcore,ierr)
    if (ierr==1) call fatal('setup','Cannot find softening length that produces nice profile (h/r(mcore) < 1.02 reached)')
 case(3) ! Both rcore and mcore are specified, check if values are sensible
    call check_hsoft_and_mcore(rcore,mcore,r,den,m,ierr)
    if (ierr==1) call fatal('setup','mcore cannot exceed m(r=h)')
    if (ierr==2) call fatal('setup','softenedrho/rho > tolerance')
    if (ierr==3) call fatal('setup','drho/dr > 0 found in softened profile')
 end select
 hsoft = 0.5*rcore ! this is set by default so that we have the original unsoftened profile for r > hsoft

 ! call core-softening subroutines
 select case(isoftcore) ! choose type of core-softneing
 case(1)
    call set_cubic_core(mcore,rcore,hsoft,den,r,pres,m,ene,temp,ierr)
    if (ierr /= 0) call fatal('setup','could not set softened core')
 case(2)
    call set_fixedS_softened_core(mcore,rcore,hsoft,den,r,pres,m,ene,temp,ierr)
    !call set_fixedS_surface(mcore,mtab,den,r,pres,en,temp,ierr)
    !call set_fixedS_star(mcore,rcore,mtab,den,r,pres,en,temp,ierr)
    if (ierr /= 0) call fatal('setup','could not set fixed entropy softened core')
 end select

end subroutine set_softened_core


!-----------------------------------------------------------------------
!+
!  Get composition (mean molecular weight, X, Z mass fractions) to be
!  their values at R/2
!+
!-----------------------------------------------------------------------
subroutine get_composition(r,den,pres,temp,X,Y,Xfixed,Zfixed,mu)
 use table_utils, only:interpolator
 use physcon,     only:radconst,kb_on_mh
 real, intent(in)     :: r(:),den(:),pres(:),temp(:),X(:),Y(:)
 real, intent(out)    :: mu,Xfixed,Zfixed
 real                 :: Rstar,pgas
 integer              :: i

 Rstar = r(size(r))
 call interpolator(r, 0.5*Rstar, i)
 pgas = pres(i) - radconst*temp(i)**4/3. ! Assuming pressure due to ideal gas + radiation
 mu = (den(i) * kb_on_mh * temp(i)) / pgas
 Xfixed = X(i)
 Zfixed = 1. - Xfixed - Y(i)

end subroutine get_composition

end module setsoftenedcore
