!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: cubic
!
!  DESCRIPTION:
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module cubic
 implicit none

contains
!-------------------------------------------------------------
! this subroutine finds the real solutions to
! a cubic equation of the form
!
! a*x^3 + b*x^2 + c*x + d
!
! formulae taken from:
! Woan, The Cambridge Handbook of Physics Formulas, 2000, p51
!
! input  : a,b,c,d : coefficients of cubic polynomial
! output : x(3)    : array containing up to 3 real solutions
!                    => x(1:nreal) non-zero, rest set to zero
!          nreal   : number of real solutions
!
! Daniel Price, 22/2/07
! dprice@astro.ex.ac.uk
!-------------------------------------------------------------
subroutine cubicsolve(a,b,c,d,x,nreal,check)
 real,    intent(in)  :: a,b,c,d
 real,    intent(out) :: x(3)
 integer, intent(out) :: nreal
 logical, intent(in), optional :: check
 real :: p,q,det,sqrtdet
 real :: a2,b2,u,v,y1,y2,y3,term,phi
 real, parameter :: eps = 1000.*epsilon(0.)
 real, parameter :: pi = 3.14159265358979323846
 integer :: i

 x = 0.
!
!--handle all trivial cases (quadratic, linear, all zero)
!
 if (abs(a) < eps) then
    det = c**2 - 4.*b*d
    if (det < 0.) then ! no solutions to quadratic
       nreal = 0
    else
       if (abs(b) < eps) then
          !--no solutions if a = 0, b = 0, c = 0
          if (abs(c) < eps) then
             nreal = 0
          else
             !--solve linear equation if a = 0, b = 0
             nreal = 1
             x(1) = -d/c
          endif
       else
          !--solve quadratic for a = 0
          nreal = 2
          sqrtdet = sqrt(det)
          x(1) = 0.5*(-c + sqrtdet)/b
          x(2) = 0.5*(-c - sqrtdet)/b
       endif
    endif
 else
!
!--cubic solution
!
    a2 = a**2
    b2 = b**2
    p = (c/a - b2/(3.*a2))
    q = (2.*b**3/(27.*a2*a) - b*c/(3.*a2) + d/a)
    det = (p**3)/27. + 0.25*q**2
!
!--determine number of solutions
!
    if (det < 0.) then
       !--3 distinct real roots
       nreal = 3
       term = sqrt(abs(p)/3.)
       phi = acos(-0.5*q*term**(-3))

       !--these are the solutions to the reduced cubic
       !  y^3 + py + q = 0
       y1 = 2.*term*cos(phi/3.)
       y2 = -2.*term*cos((phi + pi)/3.)
       y3 = -2.*term*cos((phi - pi)/3.)
    else
       !--1 real, 2 complex roots
       term = -0.5*q + sqrt(det)
       !--must take cube root of positive quantity, then give sign later
       !  (otherwise gives NaNs)
       u = (abs(term))**(1/3.)*SIGN(1.0,term)
       term = -0.5*q - sqrt(det)
       v = (abs(term))**(1/3.)*SIGN(1.0,term)
       nreal = 1
       y1 = u + v
       !--if det=0, 3 real roots, but at least 2 equal, so max of 2 unique roots)
       if (abs(det) < eps) then
          nreal = 2
          y2 = -(u + v)/2.
       endif
       y3 = 0.
    endif
    !--return solutions to original cubic, not reduced cubic
    term = b/(3.*a)
    if (nreal >= 1) x(1) = y1 - term
    if (nreal >= 2) x(2) = y2 - term
    if (nreal >= 3) x(3) = y3 - term

 endif

 if (present(check)) then
    if (check) then
       !--verify the cubic solution
       print*,'verifying: ',a,'x^3 + ',b,'x^2 + ',c,'x + ',d
       do i=1,nreal
          term = a*x(i)**3 + b*x(i)**2 + c*x(i) + d
          if (abs(term) < eps) then
             print*,'root ',i,':',x(i),'f=',term,': OK'
          else
             print*,'root ',i,':',x(i),'f=',term,': FAILED',eps
          endif
       enddo
    endif
 endif
 return

end subroutine cubicsolve

!-------------------------------------------------------------
! this subroutine returns both the real and complex
! solutions to a cubic equation of the form
!
! x^3 + a*x^2 + b*x + c
!
! input  : b,c,d : coefficients of cubic polynomial
! output : x(3)  : array of 3 COMPLEX solutions
!          nreal : number of real solutions
!
! The form of the equation above means that we
! do not need to handle trivial cases (quadratic, etc.)
! and that there will always be 3 solutions.
!
! Daniel Price, daniel.price@monash.edu 21/01/2011
!
!-------------------------------------------------------------
subroutine cubicsolve_complex(a,b,c,x,nreal,check)
 real,    intent(in)  :: a,b,c
 complex, intent(out) :: x(3)
 integer, intent(out), optional :: nreal
 logical, intent(in),  optional :: check
 real :: p,q,det,sqrtdet,sqrtp
 real :: a2,p3,t,sinht,u,v,term,termr,termi,termc
 real, parameter :: eps = 1000.*epsilon(0.)
 real, parameter :: pi = 3.14159265358979323846
 integer :: i

 x = (0.,0.)
!
!--preliminaries
!
 a2 = a*a
 p = (a2 - 3.*b)/9.
 q = a*b/6. - 0.5*c - a2*a/27.
 p3 = p*p*p
 det = p3 - q*q
!
!--cubic solution
!
 if (p < 0) then
    !--one real, two complex roots irrespective of the value of D
    sqrtp = sqrt(-p)
    t = 1./3.*asinh(q/sqrtp**3)
    sinht = sinh(t)
    x(1) = -a/3. + 2.*sqrtp*sinht ! real root
    termr = -a/3. - sqrtp*sinht
    termi = sqrt(-3.*p)*cosh(t)
    x(2) = cmplx(termr,termi)
    x(3) = cmplx(termr,-termi)
    if (present(nreal)) nreal = 1
 elseif (det > 0) then ! p > 0
    !--three real roots
    sqrtp = sqrt(p)
    term = acos(q/sqrt(p3))
    x(1) = -a/3. + 2.*sqrtp*cos(term/3.) ! real
    x(2) = -a/3. + 2.*sqrtp*cos((-2.*pi + term)/3.)
    x(3) = -a/3. + 2.*sqrtp*cos((2.*pi + term)/3.)
    if (present(nreal)) nreal = 3
 else ! p < 0 and d < 0
    !--one real, two complex roots
    sqrtdet = sqrt(-det)
    term = q - sqrtdet
    !--must take cube root of positive quantity, then give sign later
    !  (otherwise gives NaNs)
    u = (abs(term))**(1./3.)*SIGN(1.0,term)
    term = q + sqrtdet
    v = (abs(term))**(1./3.)*SIGN(1.0,term)
    x(1) = -a/3. + u + v  ! real
    termr = -a/3. - 0.5*(u + v)
    termi = sqrt(3.)*0.5*(u - v)
    x(2) = cmplx(termr,termi)
    x(3) = cmplx(termr,-termi)
    if (present(nreal)) then
       nreal = 1
       if (abs(det) < eps) nreal = 2
    endif
 endif

 !--the following lines can be used for debugging
 if (present(check)) then
    if (check) then
       !--verify the cubic solution
       print*,'verifying: x^3 + ',a,'x^2 + ',b,'x + ',c
       do i=1,3
          termc = x(i)**3 + a*x(i)**2 + b*x(i) + c
          if (abs(termc) < eps) then
             print*,'root ',i,':',x(i),'f=',termc,': OK'
          else
             print*,'root ',i,':',x(i),'f=',termc,': FAILED',eps
          endif
       enddo
    endif
 endif

 return
end subroutine cubicsolve_complex

end module cubic
