!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cullendehnen
!
! cullendehnen
!
! :References: None
!
! :Owner: Elisabeth Borchert
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

contains
!-------------------------------------------------------------------------------
!+
!  function to return alphaloc from known values of d(divv)/dt and sound speed
!  for use in Cullen & Dehnen (2010) switch
!+
!-------------------------------------------------------------------------------
pure real function get_alphaloc(divvdti,spsoundi,hi,xi_limiter,alphamin,alphamax)
 !use kernel, only:radkern
 real, intent(in) :: divvdti,spsoundi,hi,xi_limiter,alphamin,alphamax
 real :: source
 real :: temp

 source = 10.*hi**2*xi_limiter*max(-divvdti,0.)
 temp = spsoundi**2 !+ source
 if (temp > epsilon(temp)) then
    get_alphaloc = max(min(source/temp,alphamax),alphamin)
 else
    get_alphaloc = alphamin
 endif

end function get_alphaloc

pure real function xi_limiter(dvdx)
 real(kind=4), intent(in) :: dvdx(9)
 real  :: dvxdx,dvxdy,dvxdz,dvydx,dvydy,dvydz,dvzdx,dvzdy,dvzdz
 real  :: fac,traceS,divv,curlvx,curlvy,curlvz
 
 dvxdx = dvdx(1)
 dvxdy = dvdx(2)
 dvxdz = dvdx(3)
 dvydx = dvdx(4)
 dvydy = dvdx(5)
 dvydz = dvdx(6)
 dvzdx = dvdx(7)
 dvzdy = dvdx(8)
 dvzdz = dvdx(9)
 divv = dvxdx + dvydy + dvzdz
 curlvx = dvzdy - dvydz
 curlvy = dvxdz - dvzdx
 curlvz = dvydx - dvxdy

 fac    = max(-divv,0.)**2 
 traceS = curlvx**2 + curlvy**2 + curlvz**2
 if (fac + traceS > 0.) then
    xi_limiter = fac/(fac + traceS)
 else
    xi_limiter = 1.
 endif

end function xi_limiter

end module cullendehnen