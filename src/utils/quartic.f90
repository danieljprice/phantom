!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module quartic
!
! Module for finding solution to quartic equations
!
! :References: Whitehouse & Bate (2004), MNRAS 353, 1078-1094
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 public :: quarticsolve,quarticf

 private

contains
!---------------------------------------------------------
!+
!  generic solver for quartic equation in the form
!
!  x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
!
!  INPUT:
!    array of coefficients a = [a0,a1,a2,a3]
!    previous solution uold
!  OUTPUT:
!    soln: real solution
!    moresweep: whether need to refine boundaries further
!    ierr: error code in case of failure to find root
!
!  returns solution closest to uold
!+
!---------------------------------------------------------
subroutine quarticsolve(a,uold,soln,moresweep,ierr)
 !use cubic, only:cubicsolve
 real, intent(in) :: a(0:3),uold
 integer, intent(out)   :: ierr
 logical, intent(inout) :: moresweep
 real, intent(out)      :: soln
 integer :: rtst
 real :: y1,ub1,uc1,ub2,uc2,ub,uc,a0,a1,a2,a3,a4
 real :: tmin,tmax,quantity1,biggest_term
 real, dimension(2) :: z1,z2,z3,z4
 !real :: x(3)
 !integer :: nreal

 ierr = 0
 a0 = a(0)
 a1 = a(1)
 a2 = a(2)
 a3 = a(3)
 a4 = 1.

 z1(2) = 0.
 z2(2) = 0.
 z3(2) = 0.
 z4(2) = 0.

!  la1 = abs(log10(abs(a1))-10.)
!  test1 = a1

 ! Real root of equation:
 ! y^3 -4 a0 y - a1^2 = 0

 quantity1 = -54.*a2*a1**3*a3 + 12.*a3**2*a0*a2**3 - 3.*a1**2*a3**2*a2**2 - &
             432.*a2*a0*a1**2 + 18.*a1**2*a3**2*a0 + 576.*a1*a3*a0**2 - 432.*a2* &
             a0**2*a3**2 + 240.*a1*a3*a2**2*a0 - 54.*a2*a1*a3**3*a0 + 384.*a2**2* &
             a0**2-48.*a0*a2**4 + 12.*a1**2*a2**3 + 81.*a3**4*a0**2 - 768.*a0**3+ &
             81.*a1**4 + 12.*a1**3*a3**3

 biggest_term = max(abs(-54.*a2*a1**3*a3),abs(12.*a3**2*a0*a2**3), &
                    abs(-3.*a1**2*a3**2*a2**2),abs(-432.*a2*a0*a1**2), &
                    abs(18.*a1**2*a3**2*a0),abs(576.*a1*a3*a0**2), &
                    abs(-432.*a2*a0**2*a3**2),abs(240.*a1*a3*a2**2*a0), &
                    abs(-54.*a2*a1*a3**3*a0),abs(384.*a2**2*a0**2), &
                    abs(-48.*a0*a2**4),abs(12.*a1**2*a2**3),abs(81.*a3**4*a0**2), &
                    abs(-768.*a0**3),abs(81.*a1**4),abs(12.*a1**3*a3**3))

 if (quantity1 < 0. .and. abs(quantity1)/biggest_term < 1.e-12) then
    quantity1 = 0.
    print *,"q1 ",quantity1
 elseif (quantity1 < 0.) then
    print *,"QUARTIC4: Quantity1 is negative. "
    print *,"Quantity1:",quantity1,biggest_term
    print *,"Returning to TRAP with moresweep2=.TRUE."
    moresweep = .true.
    ierr = 1
    return
 endif

 y1 = (-36.*a2*a1*a3 - 288.*a2*a0 + 108.*a1**2 + 108.*a3**2*a0 + 8.*a2**3 + 12.* &
      sqrt(quantity1))**(1./3.)/6. - 6*(a1*a3/3.-4./3.*a0 - &
      a2**2/9.)/(-36.*a2*a1*a3 - 288.*a2*a0 + 108.*a1**2 + 108.*a3**2*a0 + &
      8.*a2**3 + 12.*sqrt(quantity1))**(1./3.)+a2/3.

! call cubicsolve(1.,-a2,(a1*a3-4.*a0),(4.*a2*a0-a1**2-a3**2*a0),x,nreal)
! print*,nreal,x(1:nreal) !,y1
 !if (abs(x(1)-y1) > epsilon(0.)) print*,' ERROR discrepant solutions ',x(1:nreal),y1
! y1 = x(1)

 z1(2) = 0.
 z2(2) = 0.
 z3(2) = 0.
 z4(2) = 0.

 z1(1) = 0.
 z2(1) = 0.
 z3(1) = 0.
 z4(1) = 0.

 ! Solution to quartic
 ! This is solution to two quadratics
 ! v^2+(a2/2... etc)
 ub = (0.25*a3**2 + y1 - a2)
 if (ub < 0.) then
    if (abs(ub) < 5e-15) then ! can recover if close to zero
       ub = 0.
    else
       print*,"quartic4: error, imaginary co-eff b to quadratic"
       print*,'a0,a1=',a0,a1,' y1 = ',y1,ub,epsilon(ub)
       print*,"returning to trapimpl with moresweep2=.true."
       moresweep=.true.
       ierr = 2
       return
    endif
 endif

 uc = ((0.5*y1)**2 - a0)
 if (uc < 0.) then
    print*,"QUARTIC4: Error, imaginary co-eff c to quadratic"
    ierr = 3
    return
 endif

 ub1 = 0.5*a3 + sqrt(ub)
 ub2 = 0.5*a3 - sqrt(ub)
 if (a1 < 0. .and. a0 > 0.) then
    uc1 = 0.5*y1 + sqrt(uc)
    uc2 = 0.5*y1 - sqrt(uc)
 else
    uc1 = 0.5*y1 - sqrt(uc)
    uc2 = 0.5*y1 + sqrt(uc)
 endif

 if (ub2**2 - 4.*uc2 > 0. .and. ub1**2 - 4.*uc1 < 0.) then
    if (abs(a3) >= tiny(0.)) then
       if (abs((a3/2.) - sqrt(ub))/abs(a3/2.) < 1.e-6) then
          print *,"quartic4: error, big - big / big too big for co-eff b"
          ierr = 4
          return
       endif
    endif
    z2(1)=0.5*((-ub2)-sqrt(ub2**2-4.0*uc2))
    z1(1)=0.5*((-ub2)+sqrt(ub2**2-4.0*uc2))
    z2(2)=1
    z1(2)=1
 elseif (ub2**2 - 4.*uc2 < 0. .and. ub1**2 - 4.*uc1 > 0.) then
    if (abs(y1) > tiny(0.) .and. abs(2.*a0/(y1**2 + epsilon(0.))) < 1.e-6 .and. abs(a3) >= tiny(0.)) then
       if (abs((2.*a0/(y1 + epsilon(0.)))/(((a3/2.) + sqrt((a3**2/4.) + y1 - a2))**2)) > 1.e-6) then
          print*,"quartic4: second taylor expansion no longer valid"
          print*,"quartic4: value is ",abs((a0/y1)/((a3/2.)+sqrt((a3**2/4.)+y1-a2)))
          ierr = 5
          return
       endif
       z4(1) = -1.*((a0/y1)/((a3/2.) + sqrt((a3**2/4.) + y1 - a2)))
    else
       z4(1) = 0.5*((-ub1)+sqrt(ub1**2 - 4.*uc1))
    endif

    z3(1) = 0.5*((-ub1) - sqrt(ub1**2 - 4.*uc1))
    z3(2) = 1
    z4(2) = 1
 elseif (ub2**2 - 4.*uc2 < 0. .and. ub1**2 - 4.*uc1 < 0.) then
    !!  NUMERICAL SOLUTION IF ONLY IMAGINARY ARE returnED ANALYTICALLY !!
    print*,"QUARTIC4: All imaginary roots for quartic"
    ierr = 5
    return
 else
    if (abs(a3) >= tiny(0.)) then
       if (abs((a3/2.) - sqrt(ub)) / abs(a3/2.) < 1.e-6) then
          print *,"QUARTIC4: Error, big - big / big too big for co-eff b"
          ierr = 6
          return
       endif
    endif

    if (abs(2.*a0/(y1**2 + epsilon(0.))) < 1.e-6) then
       if (abs(a3) >= tiny(0.)) then
          if ((abs(2.*a0/y1)/(((a3/2.) + sqrt((a3**2/4.) + y1 - a2))**2)) > 1.e-6) then
             print*,"quartic4: second taylor expansion no longer valid"
             print*,"quartic4: value is ", abs((a0/y1)/((a3/2.) + sqrt((a3**2/4.)+y1-a2)))
             ierr = 7
             return
          endif
          z4(1) = -1.*((a0/y1)/((a3/2.) + sqrt((a3**2/4.) + y1 - a2)))
       else
          print *,"quartic4: a3<tiny at position 2 ",a3
       endif
    else
       z4(1)=0.5*((-ub1)+sqrt(ub1**2 - 4.*uc1))
    endif

    z3(1)=0.5*((-ub1)+sqrt(ub1**2 - 4.*uc1))
    z2(1)=0.5*((-ub2)-sqrt(ub2**2 - 4.*uc2))
    z1(1)=0.5*((-ub2)+sqrt(ub2**2 - 4.*uc2))

    z4(2)=1
    z3(2)=1
    z2(2)=1
    z1(2)=1
 endif

 if (z1(2) == 1. .and. z2(2) == 1. .and. z3(2) == 0. .and. z4(2) == 0.) then
    if (z1(1) > 0. .and. z2(1) <= 0.) then
       soln = z1(1)
    elseif (z1(1) <= 0. .and. z2(1) > 0.) then
       soln = z2(1)
    elseif (z1(1) <= 0. .and. z2(1) <= 0.) then
       print*,"failed 1 ",z1(1),z2(1),z3(1),z4(1)
       print*,"coeffs = ",a0,a1,a2,a3,a4
       print*,"     ",quantity1,biggest_term
       print*,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
       return
    elseif (abs(z1(1)-uold) > abs(z2(1)-uold)) then
       soln = z2(1)
       if (abs((soln-uold)/uold) >= 1.0) then
          print*,"change big 1a ",z1(1),z2(1),z3(1),z4(1),uold
          print*,"coeffs = ",a0,a1,a2,a3,a4
          print*,"     ",quantity1,biggest_term
          print*,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
          return
       endif
    else
       soln = z1(1)
       if (abs((soln-uold)/uold) > 1.) then
          print*,"change big 1b ",z1(1),z2(1),z3(1),z4(1),uold
          print*,"coeffs = ",a0,a1,a2,a3,a4
          print*,"     ",quantity1,biggest_term
          print*,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
          return
       endif
    endif

 elseif (z3(2) == 1. .and. z4(2) == 1. .and. z1(2) == 0. .and. z2(2) == 0.) then

    if (z3(1) > 0. .and. z4(1) < 0.) then
       soln = z3(1)
    elseif (z3(1) <= 0. .and. z4(1) >= 0.) then
       soln = z4(1)
    elseif (z3(1) <= 0. .and. z4(1) <= 0.) then
       print*,"failed 2 ",z1(1),z2(1),z3(1),z4(1),uold
       print*,"coeffs = ",a0,a1,a2,a3,a4
       print*,"     ",quantity1,biggest_term
       print*,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
       return
    elseif (abs(z3(1)-uold) >= abs(z4(1)-uold)) then
       soln = z4(1)
       if (abs((soln-uold)/uold) >= 1.) then
          print*,"change big 2a ",z1(1),z2(1),z3(1),z4(1),uold
          print*,"coeffs = ",a0,a1,a2,a3,a4
          print*,"     ",quantity1,biggest_term
          print*,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
          return
       endif
    else
       soln = z3(1)
       if (abs((soln-uold)/uold) >= 1.) then
          print*,"change big 2b ",z1(1),z2(1),z3(1),z4(1),uold
          print*,"coeffs = ",a0,a1,a2,a3,a4
          print*,"     ",quantity1,biggest_term
          print*,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
          return
       endif
    endif

 elseif (z3(2) == 0. .and. z4(2) == 0. .and. z1(2) == 0. .and. z2(2) == 0.) then
    print*,"quartic4: all imaginary in "
    print*,"coeffs = ",a0,a1,a2,a3,a4
    print*,"     ",quantity1,biggest_term
    print*,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
    return
 else                      !four solutions
    rtst = 0
    write (*,*) 'four solutions ',z1(1),z2(1),z3(1),z4(1),a1,a0,uold
    soln = z1(1) ! at least return one of them...
    ierr = 1
    return
    !  if (tsoln1 >= 0. .and. tsoln1 >= tmin .and. tsoln1 <= tmax) then
    !     rtst = rtst + 1
    !     soln = z1(1)
    !  endif
    !  if (tsoln2 >= 0. .and. tsoln2 >= tmin .and. tsoln2 <= tmax) then
    !     rtst = rtst + 1
    !     soln = z2(1)
    !  endif
    !  if (tsoln1 >= 0. .and. tsoln1 >= tmin .and. tsoln1 <= tmax) then
    !     rtst = rtst + 1
    !     soln = z3(1)
    !  endif
    !  if (tsoln1 >= 0. .and. tsoln1 >= tmin .and. tsoln1 <= tmax) then
    !     rtst = rtst + 1
    !     soln = z4(1)
    !  endif

    if (rtst /= 1) then
       print*,"quartic4: there are four solutions and i'm incapable of"
       print*,"picking one. solns are: "
       print*,z1(1),z2(1),z3(1),z4(1)
       !  print*,tsoln1,tsoln2,tsoln3,tsoln4
       print*,"min..max t :",tmin,tmax
       print*,"i'm going back with moresweep2=.true."
       moresweep = .true.
    endif
 endif

end subroutine quarticsolve

real function quarticf(a,x)
 real, intent(in) :: a(0:3),x

 quarticf = a(0) + a(1)*x + a(2)*x**2 + a(3)*x**3 + x**4

end function quarticf

end module quartic
