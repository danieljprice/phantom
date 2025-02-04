!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cullendehnen
!
! Utility routines for the Cullen & Dehnen shock detection switch
!
! :References:
!   Cullen & Dehnen (2010), MNRAS 408, 669
!   Price et al. (2018), PASA 35, e031
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

!-------------------------------------------------------------------------------
!+
!  return the xi_limiter function used in the Cullen & Dehnen switch
!  based on the spatial velocity gradients
!+
!-------------------------------------------------------------------------------
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
 if (fac + traceS > epsilon(0.)) then
    xi_limiter = fac/(fac + traceS)
 else
    xi_limiter = 1.
 endif

end function xi_limiter

end module cullendehnen
