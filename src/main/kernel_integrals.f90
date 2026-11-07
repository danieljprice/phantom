!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module kernel_integrals
!
! Exact overlap integrals of a piecewise polynomial SPH kernel
! with polygonal (2D) or polyhedral (3D) cells, generalising
! Petkova, Laibe & Bonnell (2018), J. Comp. Phys. 353, 300
! to any kernel.
!
! The basic quantity is the integral of the kernel over a
! "vertex pyramid" (3D) or triangle (2D) with apex at the particle,
! bounded by the foot of the perpendicular to an edge and a vertex.
! With all lengths in units of h, the particle is a distance c from
! the wall plane (c=0 in 2D), the foot of the edge lies a distance
! p from the foot of the wall perpendicular, and the vertex lies a
! distance x along the edge. The squared distance to a point x'
! along the edge is r^2 = a + x'^2 with a = c^2 + p^2.
!
! For a kernel with radial antiderivatives (per piece j)
!    2D:  G2(q) = cnorm int_0^q w q' dq'     = sum_n g2(n,j) q^n
!    3D:  A3(s) = int G3(s)/s^2 ds           = sum_n a3(n,j) s^n
!         G3(s) = cnorm int_0^s w s'^2 ds'
! the vertex integrals reduce to sums over the terms
!    S_n = int (r)^n dphi,  phi = atan(x'/p)
! which obey the recurrence
!    S_{n+2} = c^2 S_n + p J_n,   J_n = int_0^x (a+x'^2)^(n/2) dx'
!    J_n = [x (a+x^2)^(n/2) + n a J_{n-2}]/(n+1)
! with J_0 = x, J_{-1} = asinh(x/sqrt(a)), S_0 = phi,
! S_1 = p J_{-1} + c T and T = c S_{-1} = atan(x c/(p r)).
!
! The coefficient tables are computed from the polynomial coefficients
! of the kernel by set_kernel_integral_coeffs
!
! The routines here take everything they need as arguments so
! this module has no dependencies on the rest of the code
!
! Adapted from the kernel_integrals module of SPLASH (D. Price, 2026),
! adding the kernel pieces of Phantom's kernels, the zero-second-moment
! kernel and the integral over a cell (without clipping, for kernels
! with negative parts)
!
! :References: Petkova, Laibe & Bonnell (2018), J. Comp. Phys. 353, 300
!
! :Owner:
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 integer, parameter :: dp = kind(0.d0)
 integer, parameter, public :: maxpieces = 3, maxdeg = 13, maxpow = maxdeg + 2

 public :: set_kernel_integral_coeffs, get_S_terms
 public :: eval_laurent
 public :: edge_int_2d, edge_int_3d
 public :: get_kernel_pieces, zero_moment_pieces, cell_integral_3d

 private

contains

!-----------------------------------------------------------------------
!+
!  Compute S_n for n=0..nmax (and T = c S_{-1}) at distance x along
!  an edge at perpendicular distance p, for a plane at distance c
!+
!-----------------------------------------------------------------------
pure subroutine get_S_terms(x,p,c,nmax,S,T)
 real(dp), intent(in)  :: x,p,c
 integer,  intent(in)  :: nmax
 real(dp), intent(out) :: S(0:maxpow),T
 real(dp) :: a,r,rn,J(-1:maxpow)
 integer  :: n

 a = c*c + p*p
 r = sqrt(a + x*x)
 S(0) = atan2(x,p)
 T = 0._dp
 if (c > 0._dp) T = atan2(x*c,p*r)
 if (nmax < 1) return

 ! J_n = int_0^x (a + x'^2)^(n/2) dx', by forward recurrence
 J(-1) = asinh(x/sqrt(a))
 J(0)  = x
 rn = 1._dp
 do n=1,nmax-2
    rn = rn*r
    J(n) = (x*rn + n*a*J(n-2))/(n+1)
 enddo

 ! S_n by forward recurrence (all terms positive, so stable)
 S(1) = p*J(-1) + c*T
 do n=0,nmax-2
    S(n+2) = c*c*S(n) + p*J(n)
 enddo

end subroutine get_S_terms

!-----------------------------------------------------------------------
!+
!  Sum the coefficient-weighted S_n and T terms over each piece of the
!  kernel crossed along an edge, from x'=0 to each of the two vertices
!  at x'=xv(1) and x'=xv(2). Both vertices share the same breakpoints,
!  so these are evaluated only once. Returns the sum of the S_n terms
!  (n >= 0), the T terms (n = -1) and phi = S_0 for each vertex. The
!  constant shift0 is subtracted from the S_0 coefficient on every piece
!+
!-----------------------------------------------------------------------
pure subroutine sum_over_pieces(xv,p,c,npieces,nmax,qb,coef,shift0,sumS,sumT,phi)
 real(dp), intent(in)  :: xv(2),p,c
 integer,  intent(in)  :: npieces,nmax
 real(dp), intent(in)  :: qb(maxpieces),coef(-1:maxpow,maxpieces+1),shift0
 real(dp), intent(out) :: sumS(2),sumT(2),phi(2)
 real(dp) :: a,xb,Slo(0:maxpow),Shi(0:maxpow),Tlo,Thi,accS,accT
 integer  :: i,j
 logical  :: done(2)

 a = c*c + p*p
 done = .false.

 ! find the piece of the kernel containing the foot of the edge
 j = 1
 do while (j <= npieces)
    if (a < qb(j)**2) exit
    j = j + 1
 enddo

 ! integrate piece by piece, splitting at each kernel breakpoint,
 ! finishing each vertex in the piece where it lies
 Slo  = 0._dp
 Tlo  = 0._dp
 accS = 0._dp
 accT = 0._dp
 do while (j <= npieces)
    xb = sqrt(qb(j)**2 - a)
    do i=1,2
       if (.not.done(i) .and. xv(i) <= xb) then
          call finish_vertex(xv(i),p,c,nmax,coef(:,j),shift0,Slo,Tlo,accS,accT,sumS(i),sumT(i),phi(i))
          done(i) = .true.
       endif
    enddo
    if (all(done)) return
    call get_S_terms(xb,p,c,nmax,Shi,Thi)
    call add_piece(nmax,coef(:,j),shift0,Slo,Shi,Tlo,Thi,accS,accT)
    Slo = Shi
    Tlo = Thi
    j = j + 1
 enddo

 ! vertices beyond the kernel: only S_0 and T are needed
 do i=1,2
    if (.not.done(i)) call finish_vertex(xv(i),p,c,0,coef(:,j),shift0,Slo,Tlo,accS,accT, &
                                         sumS(i),sumT(i),phi(i))
 enddo

end subroutine sum_over_pieces

!-----------------------------------------------------------------------
!+
!  Complete the sums for a vertex at x lying on the current piece
!+
!-----------------------------------------------------------------------
pure subroutine finish_vertex(x,p,c,nmax,coef,shift0,Slo,Tlo,accS,accT,sumS,sumT,phi)
 real(dp), intent(in)  :: x,p,c,coef(-1:maxpow),shift0,Slo(0:maxpow),Tlo,accS,accT
 integer,  intent(in)  :: nmax
 real(dp), intent(out) :: sumS,sumT,phi
 real(dp) :: Shi(0:maxpow),Thi

 sumS = accS
 sumT = accT
 call get_S_terms(x,p,c,nmax,Shi,Thi)
 call add_piece(nmax,coef,shift0,Slo,Shi,Tlo,Thi,sumS,sumT)
 phi = Shi(0)

end subroutine finish_vertex

!-----------------------------------------------------------------------
!+
!  Add the contribution from one sub-interval on a single piece
!+
!-----------------------------------------------------------------------
pure subroutine add_piece(nmax,coef,shift0,Slo,Shi,Tlo,Thi,sumS,sumT)
 integer,  intent(in)    :: nmax
 real(dp), intent(in)    :: coef(-1:maxpow),shift0,Slo(0:maxpow),Shi(0:maxpow),Tlo,Thi
 real(dp), intent(inout) :: sumS,sumT
 integer :: n

 sumS = sumS + (coef(0) - shift0)*(Shi(0) - Slo(0))
 do n=1,nmax
    sumS = sumS + coef(n)*(Shi(n) - Slo(n))
 enddo
 sumT = sumT + coef(-1)*(Thi - Tlo)

end subroutine add_piece

!-----------------------------------------------------------------------
!+
!  Integrals of the normalised 2D kernel over the two triangles formed
!  by the particle, the foot of the perpendicular to an edge at
!  distance q0 and each of the two vertices, at distances x1 and x2
!  along the edge (all in units of h)
!+
!-----------------------------------------------------------------------
pure subroutine edge_int_2d(x1,x2,q0,npieces,nmax,qb,g2,v1,v2)
 real(dp), intent(in)  :: x1,x2,q0
 integer,  intent(in)  :: npieces,nmax
 real(dp), intent(in)  :: qb(maxpieces),g2(-1:maxpow,maxpieces+1)
 real(dp), intent(out) :: v1,v2
 real(dp) :: sumS(2),sumT(2),phi(2)

 v1 = 0._dp
 v2 = 0._dp
 if (q0 <= tiny(q0)) return
 call sum_over_pieces((/x1,x2/),q0,0._dp,npieces,nmax,qb,g2,0._dp,sumS,sumT,phi)
 v1 = sumS(1)
 v2 = sumS(2)

end subroutine edge_int_2d

!-----------------------------------------------------------------------
!+
!  Integrals of the normalised 3D kernel over the two vertex pyramids
!  with apex at the particle, at distance c from the wall plane, whose
!  bases are the right triangles formed by the foot of the perpendicular
!  to the wall, the foot of the perpendicular to an edge at distance p,
!  and each of the two vertices at distances x1 and x2 along the edge
!  (all in units of h)
!+
!-----------------------------------------------------------------------
pure subroutine edge_int_3d(x1,x2,p,c,npieces,nmax,qb,a3,v1,v2)
 real(dp), intent(in)  :: x1,x2,p,c
 integer,  intent(in)  :: npieces,nmax
 real(dp), intent(in)  :: qb(maxpieces),a3(-1:maxpow,maxpieces+1)
 real(dp), intent(out) :: v1,v2
 real(dp) :: sumS(2),sumT(2),phi(2),cA,A0,v(2)

 v1 = 0._dp
 v2 = 0._dp
 if (p <= tiny(p) .or. c <= tiny(c)) return
 ! the constant term of A3(c) is subtracted analytically from each piece
 ! to avoid cancellation of large terms when c >> 1
 call c_times_A3(c,npieces,nmax,qb,a3,cA,A0)
 call sum_over_pieces((/x1,x2/),p,c,npieces,nmax,qb,a3,A0,sumS,sumT,phi)
 v = c*sumS + sumT - cA*phi
 v1 = v(1)
 v2 = v(2)

end subroutine edge_int_3d

!-----------------------------------------------------------------------
!+
!  Evaluate c*A3(c) on the piece of the kernel containing c, giving
!  the constant that sets the radial integral to zero at the wall.
!  Returns c*A3(c) excluding its constant term (cA) and the constant
!  term itself (A0) separately
!+
!-----------------------------------------------------------------------
pure subroutine c_times_A3(c,npieces,nmax,qb,a3,cA,A0)
 real(dp), intent(in)  :: c
 integer,  intent(in)  :: npieces,nmax
 real(dp), intent(in)  :: qb(maxpieces),a3(-1:maxpow,maxpieces+1)
 real(dp), intent(out) :: cA,A0
 integer :: j,n

 j = 1
 do while (j <= npieces)
    if (c < qb(j)) exit
    j = j + 1
 enddo
 cA = 0._dp
 do n=nmax,1,-1
    cA = cA*c + a3(n,j)
 enddo
 cA = cA*c*c + a3(-1,j)
 A0 = a3(0,j)

end subroutine c_times_A3

!-----------------------------------------------------------------------
!+
!  Evaluate sum_n c(n) x^n, where n = -1 is allowed (a 1/x term)
!+
!-----------------------------------------------------------------------
pure real(dp) function eval_laurent(c,x) result(val)
 real(dp), intent(in) :: c(-1:maxpow),x
 integer :: n

 val = 0._dp
 do n=maxpow,0,-1
    val = val*x + c(n)
 enddo
 if (abs(c(-1)) > 0._dp) val = val + c(-1)/x

end function eval_laurent

!-----------------------------------------------------------------------
!+
!  Compute the coefficient tables for the exact kernel integrals from
!  the monomial coefficients w(k,j) of the (unnormalised) kernel on each
!  piece j, which extends from qb(j-1) (or 0) to qb(j). The kernel is
!  normalised here, and the normalisation constants are returned. On
!  each piece, and on the region beyond the kernel (j = npieces+1):
!
!     G2(q) = cnorm2D int_0^q w q' dq'
!     G3(s) = cnorm3D int_0^s w s'^2 ds'
!     A3(s) = int G3(s)/s^2 ds
!
!  with integration constants set so that each is continuous across the
!  breakpoints, with G2(0) = G3(0) = 0 (A3 only needs to be continuous)
!+
!-----------------------------------------------------------------------
pure subroutine set_kernel_integral_coeffs(npieces,qb,w,nmax,g2,a3,cnorm2D,cnorm3D)
 integer,  intent(in)  :: npieces
 real(dp), intent(in)  :: qb(maxpieces),w(0:maxdeg,maxpieces)
 integer,  intent(out) :: nmax
 real(dp), intent(out) :: g2(-1:maxpow,maxpieces+1),a3(-1:maxpow,maxpieces+1)
 real(dp), intent(out) :: cnorm2D,cnorm3D
 real(dp), parameter :: pi = 4._dp*atan(1._dp)
 real(dp) :: g3s(-1:maxpow,maxpieces+1),g3const(maxpieces+1),int2,int3,lo
 integer  :: j,k,n

 ! antiderivatives of w q and w q^2 on each piece, the latter
 ! stored as G3(s)/s to keep the powers within maxpow
 g2  = 0._dp
 g3s = 0._dp
 do j=1,npieces
    do k=0,maxdeg
       g2(k+2,j)  = w(k,j)/(k+2)
       g3s(k+2,j) = w(k,j)/(k+3)
    enddo
 enddo

 ! normalise so that the kernel integrates to unity in 2D and 3D
 int2 = 0._dp
 int3 = 0._dp
 lo = 0._dp
 do j=1,npieces
    int2 = int2 + eval_laurent(g2(:,j),qb(j)) - eval_laurent(g2(:,j),lo)
    int3 = int3 + qb(j)*eval_laurent(g3s(:,j),qb(j)) - lo*eval_laurent(g3s(:,j),lo)
    lo = qb(j)
 enddo
 cnorm2D = 1._dp/(2._dp*pi*int2)
 cnorm3D = 1._dp/(4._dp*pi*int3)
 g2  = cnorm2D*g2
 g3s = cnorm3D*g3s

 ! constants making G2 and G3 continuous across breakpoints, zero at the origin
 call set_continuous(npieces,qb,g2)
 g3const(1) = 0._dp
 do j=2,npieces+1
    g3const(j) = g3const(j-1) + qb(j-1)*(eval_laurent(g3s(:,j-1),qb(j-1)) &
                                       - eval_laurent(g3s(:,j),qb(j-1)))
 enddo

 ! A3 = int G3/s^2 ds: s^(n+1) -> s^n/n, with the constant -> -1/s
 a3 = 0._dp
 do j=1,npieces+1
    a3(-1,j) = -g3const(j)
    do n=2,maxpow
       a3(n,j) = g3s(n,j)/n
    enddo
 enddo
 call set_continuous(npieces,qb,a3)

 ! highest power needed in the recurrences
 nmax = 0
 do n=0,maxpow
    if (any(abs(g2(n,:)) > 0._dp) .or. any(abs(a3(n,:)) > 0._dp)) nmax = n
 enddo

end subroutine set_kernel_integral_coeffs

!-----------------------------------------------------------------------
!+
!  Set the constant term c(0,j) on each piece j > 1 so that the
!  piecewise function is continuous across each breakpoint
!+
!-----------------------------------------------------------------------
pure subroutine set_continuous(npieces,qb,c)
 integer,  intent(in)    :: npieces
 real(dp), intent(in)    :: qb(maxpieces)
 real(dp), intent(inout) :: c(-1:maxpow,maxpieces+1)
 integer :: j

 do j=2,npieces+1
    c(0,j) = 0._dp
    c(0,j) = eval_laurent(c(:,j-1),qb(j-1)) - eval_laurent(c(:,j),qb(j-1))
 enddo

end subroutine set_continuous

!-----------------------------------------------------------------------
!+
!  Monomial coefficients w(k,j) of q^k of each (unnormalised) kernel on
!  each piece j, from qb(j-1) (or zero) to qb(j), for the kernel with
!  the given name; ierr = 1 if the kernel is not known
!+
!-----------------------------------------------------------------------
pure subroutine get_kernel_pieces(name,npieces,qb,w,ierr)
 character(len=*), intent(in)  :: name
 integer,          intent(out) :: npieces,ierr
 real(dp),         intent(out) :: qb(maxpieces),w(0:maxdeg,maxpieces)

 ! integer coefficients over a common denominator, so that the
 ! conversion to double precision is exact
 ierr = 0
 qb = 0
 w  = 0
 select case(trim(name))
 case('Wendland 2/3D C^6') ! (1 - q/2)^8 (1 + 4q + 25/4 q^2 + 4q^3)
    npieces = 1
    qb(1) = 2
    w(0:11,1) = (/1024,0,-2816,0,4224,0,-7392,8448,-4620,1408,-231,16/)/1024._dp
 case('Wendland 2/3D C^4') ! (1 - q/2)^6 (1 + 3q + 35/12 q^2)
    npieces = 1
    qb(1) = 2
    w(0:8,1) = (/768,0,-1792,0,3360,-3584,1680,-384,35/)/768._dp
 case('Wendland 2/3D C^2') ! (1 - q/2)^4 (1 + 2q)
    npieces = 1
    qb(1) = 2
    w(0:5,1) = (/16,0,-40,40,-15,2/)/16._dp
 case('M_6 quintic') ! (3-q)^5 - 6(2-q)^5 + 15(1-q)^5, (3-q)^5 - 6(2-q)^5, (3-q)^5
    npieces = 3
    qb(1:3) = (/1,2,3/)
    w(0:5,1) = (/66,0,-60,0,30,-10/)
    w(0:5,2) = (/51,75,-210,150,-45,5/)
    w(0:5,3) = (/243,-405,270,-90,15,-1/)
 case('M_5 quartic') ! (5/2-q)^4 - 5(3/2-q)^4 + 10(1/2-q)^4, (5/2-q)^4 - 5(3/2-q)^4, (5/2-q)^4
    npieces = 3
    qb(1:3) = (/1,3,5/)/2._dp
    w(0:4,1) = (/230,0,-240,0,96/)/16._dp
    w(0:4,2) = (/220,80,-480,320,-64/)/16._dp
    w(0:4,3) = (/625,-1000,600,-160,16/)/16._dp
 case('M_4 cubic') ! 1 - 3/2 q^2 + 3/4 q^3, 1/4 (2-q)^3
    npieces = 2
    qb(1:2) = (/1,2/)
    w(0:3,1) = (/4,0,-6,3/)/4._dp
    w(0:3,2) = (/8,-12,6,-1/)/4._dp
 case default
    npieces = 0
    ierr = 1
 end select

end subroutine get_kernel_pieces

!-----------------------------------------------------------------------
!+
!  Replace the kernel pieces w by those of w (1 - A q^2), with A chosen
!  so that the second moment of the 3D kernel vanishes, A = <q^2>/<q^4>
!  (the moments of q^2 w q^2 dq are integrated exactly piece by piece).
!  The normalisation is applied later by set_kernel_integral_coeffs
!+
!-----------------------------------------------------------------------
pure subroutine zero_moment_pieces(npieces,qb,w,a)
 integer,  intent(in)    :: npieces
 real(dp), intent(in)    :: qb(maxpieces)
 real(dp), intent(inout) :: w(0:maxdeg,maxpieces)
 real(dp), intent(out)   :: a
 real(dp) :: m4,m6,lo
 integer  :: j,k

 m4 = 0._dp
 m6 = 0._dp
 lo = 0._dp
 do j=1,npieces
    do k=0,maxdeg-2
       m4 = m4 + w(k,j)*(qb(j)**(k+5) - lo**(k+5))/(k+5)
       m6 = m6 + w(k,j)*(qb(j)**(k+7) - lo**(k+7))/(k+7)
    enddo
    lo = qb(j)
 enddo
 a = m4/m6
 do j=1,npieces
    w(2:maxdeg,j) = w(2:maxdeg,j) - a*w(0:maxdeg-2,j)
 enddo

end subroutine zero_moment_pieces

!-----------------------------------------------------------------------
!+
!  Integral of the normalised 3D kernel of a particle at xp with
!  smoothing length h over the box centred on xc with widths dx, by
!  the divergence theorem as the sum of the integrals over its six
!  faces (Petkova, Laibe & Bonnell 2018). Unlike the rendering version
!  the face integrals are not clipped at zero, so that kernels with
!  negative parts are also integrated exactly
!+
!-----------------------------------------------------------------------
pure real(dp) function cell_integral_3d(xp,xc,dx,h,npieces,nmax,qb,a3) result(vint)
 real(dp), intent(in) :: xp(3),xc(3),dx(3),h,qb(maxpieces),a3(-1:maxpow,maxpieces+1)
 integer,  intent(in) :: npieces,nmax
 real(dp) :: d(3)
 integer  :: ia,ib,ic

 d = xc - xp
 vint = 0._dp
 ! faces normal to each axis ic, with in-plane axes ia, ib
 do ic=1,3
    ia = modulo(ic,3) + 1
    ib = modulo(ic+1,3) + 1
    vint = vint + face_integral( d(ic) + 0.5_dp*dx(ic),d(ia),d(ib),dx(ia),dx(ib),h,npieces,nmax,qb,a3) &
                + face_integral(-d(ic) + 0.5_dp*dx(ic),d(ia),d(ib),dx(ia),dx(ib),h,npieces,nmax,qb,a3)
 enddo

end function cell_integral_3d

!-----------------------------------------------------------------------
!+
!  Contribution of one rectangular face at signed distance r0 from the
!  particle (positive if the particle is inside the face's half space),
!  with the face centre offset by (da, db) in the plane and widths wa, wb
!+
!-----------------------------------------------------------------------
pure real(dp) function face_integral(r0,da,db,wa,wb,h,npieces,nmax,qb,a3) result(f)
 real(dp), intent(in) :: r0,da,db,wa,wb,h,qb(maxpieces),a3(-1:maxpow,maxpieces+1)
 integer,  intent(in) :: npieces,nmax

 f = edge_part(r0,0.5_dp*wb + db,0.5_dp*wa - da,0.5_dp*wa + da,h,npieces,nmax,qb,a3) &
   + edge_part(r0,0.5_dp*wb - db,0.5_dp*wa + da,0.5_dp*wa - da,h,npieces,nmax,qb,a3) &
   + edge_part(r0,0.5_dp*wa + da,0.5_dp*wb + db,0.5_dp*wb - db,h,npieces,nmax,qb,a3) &
   + edge_part(r0,0.5_dp*wa - da,0.5_dp*wb - db,0.5_dp*wb + db,h,npieces,nmax,qb,a3)

end function face_integral

!-----------------------------------------------------------------------
!+
!  contribution from one edge of a face: r0 is the distance to the face
!  plane, p0 the distance in the plane to the edge and d1, d2 the
!  distances along the edge to each vertex
!+
!-----------------------------------------------------------------------
pure real(dp) function edge_part(r0,p0,d1,d2,h,npieces,nmax,qb,a3) result(f)
 real(dp), intent(in) :: r0,p0,d1,d2,h,qb(maxpieces),a3(-1:maxpow,maxpieces+1)
 integer,  intent(in) :: npieces,nmax
 real(dp) :: hi1,v1,v2

 f = 0._dp
 if (abs(r0) < tiny(r0)) return
 hi1 = 1._dp/h
 call edge_int_3d(abs(d1)*hi1,abs(d2)*hi1,abs(p0)*hi1,abs(r0)*hi1,npieces,nmax,qb,a3,v1,v2)
 ! combine the two vertex integrals depending on where the foot of the edge lies
 if (d1*d2 >= 0._dp) then
    f = v1 + v2
 elseif (abs(d1) < abs(d2)) then
    f = v2 - v1
 else
    f = v1 - v2
 endif
 f = sign(1._dp,r0)*sign(1._dp,p0)*f

end function edge_part

end module kernel_integrals
