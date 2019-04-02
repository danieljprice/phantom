!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: mathfunc
!
!  DESCRIPTION:
!  Contains utility routines for mathematical functions
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
module mathfunc
 implicit none

 public :: bessi0_s, bessi1_s, bessk0_s, bessk1_s,IK01A
 public :: poly,gegenbauer_poly,legendre_associated

 private

contains


!--------------------------------------------------------------------
!Gegenbauer polynomials
!--------------------------------------------------------------------
subroutine gegenbauer_poly( n, alpha, x, cx )
!
!Taken from the Polpak http://people.sc.fsu.edu/~jburkardt/
!
!*****************************************************************************80
!
!! GEGENBAUER_POLY computes the Gegenbauer polynomials C(I,ALPHA,X).
!
!  Discussion:
!
!    The Gegenbauer polynomial can be evaluated in Mathematica with
!    the command
!
!      GegenbauerC[n,m,x]
!
!    ALPHA must be greater than -0.5.
!
!    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
!    polynomials of the second kind.
!
!  Differential equation:
!
!    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + N (N + 2 ALPHA) Y = 0
!
!  Recursion:
!
!    C(0,ALPHA,X) = 1,
!    C(1,ALPHA,X) = 2*ALPHA*X
!    C(N,ALPHA,X) = ( (2*N-2+2*ALPHA) * X * C(N-1,ALPHA,X)
!                   + ( -N+2-2*ALPHA)     * C(N-2,ALPHA,X) ) / N
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 )
!      ( 1 - X^2 )^( ALPHA - 0.5 ) * C(N,ALPHA,X)^2 dX
!
!    = PI * 2^( 1 - 2 * ALPHA ) * Gamma ( N + 2 * ALPHA )
!      / ( N! * ( N + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) ALPHA, a parameter which is part of the
!    definition of the Gegenbauer polynomials.  It must be greater than -0.5.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 Gegenbauer
!    polynomials at the point X.
!

 integer, intent(in) :: n
 real(kind=8), intent(in)  :: alpha
 real(kind=8), intent(out) :: cx(0:n)
 real(kind=8), intent(in)  :: x
 integer (kind=4) :: i

 if ( alpha <= -0.5D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEGENBAUER_POLY - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal value of ALPHA = ', alpha
    write ( *, '(a)' ) '  but ALPHA must be greater than -0.5.'
    return
 endif

 if ( n < 0 ) then
    return
 endif

 cx(0) = 1.0D+00

 if ( n == 0 ) then
    return
 endif

 cx(1) = 2.0D+00 * alpha * x

 do i = 2, n
    cx(i) = &
      ( ( real ( 2 * i - 2, kind = 8 ) + 2.0D+00 * alpha ) * x * cx(i-1)   &
      + ( real (   - i + 2, kind = 8 ) - 2.0D+00 * alpha )     * cx(i-2) ) &
      /   real (     i,     kind = 8 )
 enddo
 return
end subroutine

!--------------------------------------------------------------------
!Associated Legendre polynomials
!--------------------------------------------------------------------
subroutine legendre_associated( n, m, x, cx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED evaluates the associated Legendre functions.
!
!  Differential equation:
!
!    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N,X) )
!
!    P00 =    1
!    P10 =    1 X
!    P20 = (  3 X^2 -     1)/2
!    P30 = (  5 X^3 -   3 X)/2
!    P40 = ( 35 X^4 -  30 X^2 +     3)/8
!    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
!    P60 = (231 X^6 - 315 X^4 + 105 X^2 -    5)/16
!    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
!
!    M = 1
!
!    P01 =   0
!    P11 =   1 * sqrt(1-X*X)
!    P21 =   3 * sqrt(1-X*X) * X
!    P31 = 1.5 * sqrt(1-X*X) * (5*X*X-1)
!    P41 = 2.5 * sqrt(1-X*X) * (7*X*X*X-3*X)
!
!    M = 2
!
!    P02 =   0
!    P12 =   0
!    P22 =   3 * (1-X*X)
!    P32 =  15 * (1-X*X) * X
!    P42 = 7.5 * (1-X*X) * (7*X*X-1)
!
!    M = 3
!
!    P03 =   0
!    P13 =   0
!    P23 =   0
!    P33 =  15 * (1-X*X)^1.5
!    P43 = 105 * (1-X*X)^1.5 * X
!
!    M = 4
!
!    P04 =   0
!    P14 =   0
!    P24 =   0
!    P34 =   0
!    P44 = 105 * (1-X*X)^2
!
!  Recursion:
!
!    if N < M:
!      P(N,M) = 0
!    if N = M:
!      P(N,M) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      P(N,M) = X*(2*M+1)*P(M,M)
!    if M+1 < N:
!      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
!
!  Special values:
!
!    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
!    function of the first kind equals the Legendre polynomial of the
!    first kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to be
!    evaluated.  X must satisfy -1 <= X <= 1.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 functions.
!


 integer, intent(in) :: n, m
 real(kind=8), intent(in)  :: x
 real(kind=8), intent(out) :: cx(0:n)
 real(kind=8) :: fact
 integer :: i
 real(kind=8) :: somx2

 if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    return
 endif

 if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    return
 endif

 if ( x < -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no less than -1.'
    return
 endif

 if ( 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no more than 1.'
    return
 endif

 cx(0:m-1) = 0.0D+00

 cx(m) = 1.0D+00
 somx2 = sqrt ( 1.0D+00 - x * x )

 fact = 1.0D+00
 do i = 1, m
    cx(m) = -cx(m) * fact * somx2
    fact = fact + 2.0D+00
 enddo

 if ( m + 1 <= n ) then
    cx(m+1) = x * real ( 2 * m + 1, kind = 8 ) * cx(m)
 endif

 do i = m + 2, n
    cx(i) = ( real ( 2 * i     - 1, kind = 8 ) * x * cx(i-1) &
            + real (   - i - m + 1, kind = 8 ) *     cx(i-2) ) &
            / real (     i - m,     kind = 8 )
 enddo

 return
end subroutine

!--------------------------------------------------------------------
!Below are a set of slightly modified functions from
!Numerical Recipes (Press et al. 1996), Bessel and poly functions.
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! Modified Bessel function I0, in scalar form.
!--------------------------------------------------------------------
real(kind=8) function bessi0_s(x)
 real(kind=8), intent(in) :: x
 !real :: bessi0_s,poly
 real(kind=8) :: ax
 real(kind=8), dimension(7) :: p = (/1.0d0,3.5156229d0,&
     3.0899424d0,1.2067492d0,0.2659732d0,0.360768d-1,&
     0.45813d-2/)
 real(kind=8), dimension(9) :: q = (/0.39894228d0,0.1328592d-1,&
     0.225319d-2,-0.157565d-2,0.916281d-2,&
     -0.2057706d-1,0.2635537d-1,-0.1647633d-1,&
     0.392377d-2/)
 ax=abs(x)
 if (ax < 3.75) then
    bessi0_s=poly(real((x/3.75d0)**2,8),p)
 else
    bessi0_s=(exp(ax)/sqrt(ax))*poly(real(3.75d0/ax,8),q)
 endif

end function bessi0_s

!--------------------------------------------------------------------
! modified bessel function i1, in scalar form
!--------------------------------------------------------------------
real(kind=8) function bessi1_s(x)
 real(kind=8), intent(in) :: x
 !real :: bessi1_s,poly
 real(kind=8) :: ax
 real(kind=8), dimension(7) :: p = (/0.5d0,0.87890594d0,&
     0.51498869d0,0.15084934d0,0.2658733d-1,&
     0.301532d-2,0.32411d-3/)
 real(kind=8), dimension(9) :: q = (/0.39894228d0,-0.3988024d-1,&
     -0.362018d-2,0.163801d-2,-0.1031555d-1,&
     0.2282967d-1,-0.2895312d-1,0.1787654d-1,&
     -0.420059d-2/)

 ax=abs(x)
 if (ax < 3.75) then
    bessi1_s=ax*poly(real((x/3.75d0)**2,8),p)
 else
    bessi1_s=(exp(ax)/sqrt(ax))*poly(real(3.75d0/ax,8),q)
 endif

 if (x < 0.0) bessi1_s=-bessi1_s

end function bessi1_s

!--------------------------------------------------------------------
! modified bessel function k0, in scalar form
!--------------------------------------------------------------------
real(kind=8) function bessk0_s(x)
 real(kind=8), intent(in) :: x
 real(kind=8) :: y
 real(kind=8), dimension(7) :: p = (/-0.57721566d0,0.42278420d0,&
     0.23069756d0,0.3488590d-1,0.262698d-2,0.10750d-3,&
     0.74d-5/)
 real(kind=8), dimension(7) :: q = (/1.25331414d0,-0.7832358d-1,&
     0.2189568d-1,-0.1062446d-1,0.587872d-2,&
     -0.251540d-2,0.53208d-3/)
 !call assert(x > 0.0, 'bessk0_s arg')
 if (x <= 2.0) then
    y=x*x/4.0d0
    bessk0_s=(-log(x/2.0d0)*bessi0_s(x))+poly(y,p)
 else
    y=(2.0d0/x)
    bessk0_s=(exp(-x)/sqrt(x))*poly(y,q)
 endif

end function bessk0_s

!--------------------------------------------------------------------
! modified bessel function k1, in scalar form
!--------------------------------------------------------------------
real(kind=8) function bessk1_s(x)
 real(kind=8), intent(in) :: x
 !real :: bessk1_s,bessi1_s,poly
 real(kind=8) :: y
 real(kind=8), dimension(7) :: p = (/1.0d0,0.15443144d0,&
     -0.67278579d0,-0.18156897d0,-0.1919402d-1,&
     -0.110404d-2,-0.4686d-4/)
 real(kind=8), dimension(7) :: q = (/1.25331414d0,0.23498619d0,&
     -0.3655620d-1,0.1504268d-1,-0.780353d-2,&
     0.325614d-2,-0.68245d-3/)
 !call assert(x > 0.0, 'bessk1_s arg')
 if (x <= 2.0) then
    y=x*x/4.0d0
    bessk1_s=(log(x/2.0d0)*bessi1_s(x))+(1.0d0/x)*poly(y,p)
 else
    y=2.0d0/x
    bessk1_s=(exp(-x)/sqrt(x))*poly(y,q)
 endif

end function bessk1_s


!-----------------------------------------------------------------------------
! Below taken from http://jin.ece.illinois.edu/routines/routines.html
! Provides a higher accuracy than the bessel functions above
!       =========================================================
!       Purpose: Compute modified Bessel functions I0(x), I1(1),
!                K0(x) and K1(x), and their derivatives
!       Input :  x   --- Argument ( x Ú 0 )
!       Output:  BI0 --- I0(x)
!                DI0 --- I0'(x)
!                BI1 --- I1(x)
!                DI1 --- I1'(x)
!                BK0 --- K0(x)
!                DK0 --- K0'(x)
!                BK1 --- K1(x)
!                DK1 --- K1'(x)
!       =========================================================
!
!-----------------------------------------------------------------------------
subroutine ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
 !implicit double precision (a-h,o-z)
 real(kind=8), intent(in)  :: x
 real(kind=8), intent(out) :: bi0,di0,bi1,di1,bk0,dk0,bk1,dk1
 real(kind=8) :: a(12),b(12),a1(8)
 real(kind=8) :: pi,el,x2,r,ca,cb,ct,ww,xr,xr2,w0
 integer :: k,k0

 pi=3.141592653589793d0
 el=0.5772156649015329d0
 x2=x*x
 if (abs(x) < tiny(0.d0)) then
    bi0=1.0d0
    bi1=0.0d0
    bk0=1.0d+300
    bk1=1.0d+300
    di0=0.0d0
    di1=0.5d0
    dk0=-1.0d+300
    dk1=-1.0d+300
    return
 else if (x <= 18.0d0) then
    bi0=1.0d0
    r=1.0d0
    do k=1,50
       r=0.25d0*r*x2/(k*k)
       bi0=bi0+r
       if (dabs(r/bi0) < 1.0d-15) exit
    enddo
    bi1=1.0d0
    r=1.0d0
    do k=1,50
       r=0.25d0*r*x2/(k*(k+1))
       bi1=bi1+r
       if (dabs(r/bi1) < 1.0d-15) exit
    enddo
    bi1=0.5d0*x*bi1
 else
    a = (/0.125d0,7.03125d-2, &
          7.32421875d-2,1.1215209960938d-1,&
          2.2710800170898d-1,5.7250142097473d-1,&
          1.7277275025845d0,6.0740420012735d0,&
          2.4380529699556d01,1.1001714026925d02,&
          5.5133589612202d02,3.0380905109224d03/)
    b = (/-0.375d0,-1.171875d-1,&
          -1.025390625d-1,-1.4419555664063d-1,&
          -2.7757644653320d-1,-6.7659258842468d-1,&
          -1.9935317337513d0,-6.8839142681099d0,&
          -2.7248827311269d01,-1.2159789187654d02,&
          -6.0384407670507d02,-3.3022722944809d03/)
    k0=12
    if (x >= 35.0) k0=9
    if (x >= 50.0) k0=7
    ca=dexp(x)/dsqrt(2.0d0*pi*x)
    bi0=1.0d0
    xr=1.0d0/x
    do k=1,k0
       bi0=bi0+a(k)*xr**k
    enddo
    bi0=ca*bi0
    bi1=1.0d0
    do k=1,k0
       bi1=bi1+b(k)*xr**k
    enddo
    bi1=ca*bi1
 endif
 if (x <= 9.0d0) then
    ct=-(dlog(x/2.0d0)+el)
    bk0=0.0d0
    w0=0.0d0
    ww=bk0   ! added by DJP, fixes compiler warning
    r=1.0d0
    do k=1,50
       w0=w0+1.0d0/k
       r=0.25d0*r/(k*k)*x2
       bk0=bk0+r*(w0+ct)
       if (dabs((bk0-ww)/bk0) < 1.0d-15) exit
       ww=bk0
    enddo
    bk0=bk0+ct
 else
    a1 = (/0.125d0,0.2109375d0,&
           1.0986328125d0,1.1775970458984d01,&
           2.1461706161499d02,5.9511522710323d03,&
           2.3347645606175d05,1.2312234987631d07/)
    cb=0.5d0/x
    xr2=1.0d0/x2
    bk0=1.0d0
    do k=1,8
       bk0=bk0+a1(k)*xr2**k
    enddo
    bk0=cb*bk0/bi0
 endif
 bk1=(1.0d0/x-bi1*bk0)/bi0
 di0=bi1
 di1=bi0-bi1/x
 dk0=-bk1
 dk1=-bk0-bk1/x

 return
end subroutine

!--------------------------------------------------------------------
! polynomial function
!--------------------------------------------------------------------
real(kind=8) function poly(x,coeffs)
 real(kind=8) :: coeffs(5)
 integer :: i,n
 real(kind=8) :: x

 n=size(coeffs)
 if (n <= 0) then
    poly=0.0
 else
    poly=coeffs(n)
    do i=n-1,1,-1
       poly=x*poly+coeffs(i)
    enddo
 endif

end function poly
!--------------------------------------------------------------------

end module mathfunc
