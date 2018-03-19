!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: icosahedron
!
!  DESCRIPTION:
!   Icosahedron package for pixelising the sphere
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io
!+
!--------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       THE ICOSAHEDRON PACKAGE FOR PIXELIZING THE SPHERE        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Written by Max Tegmark, Max-Planck-Institut fuer Physik, Munich
!    April 1996
!    Currently I'm at Univ. of Pennsylvania, max@physics.upenn.edu
!
!    Modifications by Daniel Price daniel.price@monash.edu:
!    19/12/2012 Converted to Fortran 90, added error calls instead
!     of pause statements for incorporation into Phantom
!
!  WHAT IS IT?
!    This FORTRAN package lets the user pixelize the sphere at
!    a wide range of resolutions. It was written primarily for
!    map-making in astronomy and cosmology. It is also useful
!    for doing integrals over the sphere when the integrand is
!    expensive to evaluate, so that one wishes to minimize the
!    number of points used.
!
!  DOCUMENTATION:
!    The package and its purpose is described in detail in
!    a postscript file available from
!      http://www.mpa-garching.mpg.de/~max/icosahedron.html
!    (faster from Europe) and from from
!      http://sns.ias.edu.edu/~max/icosahedron.html
!    (faster from the US). This site also contains the latest
!    version of the source code.
!
!  RULES:
!    The package is public domain, which means that you are
!    allowed to use it for any non-commercial purpose whatsoever
!    free of charge. The only requirement is that you include an
!    appropriate acknowledgement in any publications based on
!    work where the package has been used. Also, if you
!    redistribute the package, you must not remove this text.
!
!  HOW IT WORKS:
!    As a supplement to the above-mentioned postscript file,
!    here is a brief summary of the nitty-gritty details.
!    To use the package, first call the subroutines
!    compute_matrices and compute_corners once and for all,
!    as in the demo routine below. This precomputes some
!    geometric stuff to save time later.
!    The RESOLUTION is a positive integer 1, 2, 3, ... that you
!    can choose freely. It determines the number of pixels, which
!    is given by
!       N = 40*resolution*(resolution-1)+12.
!    For instance, resolution=13 gives N=6252 pixels.
!    The only subroutines you ever need to use are these two:
!    * vector2pixel takes a point on the sphere (specified by a
!      unit vector) and returns the number of the pixel to which
!      this point belongs, i.e., an integer between 0 and N-1.
!    * pixel2vector does the opposite: it takes a pixel number and
!      computes the corresponding unit vector.
!    The subroutine "demo" below illustrates how the routines are used.
!    It produces a text file with the coordinates of all pixels
!    together with their pixel numbers. It also outputs the pixel
!    number reconstructed with vector2pixel, so you can verify that
!    it all works as it should by checking that the first two columns
!    in the file are identical.
!
!  YOU DON'T NEED TO KNOW THIS:
!    The resolution is defined so that the number of pixels along
!    the side of a triangular face is 2*resolution, so there are
!    resolution*(2*resolution+1) pixels on each of the 20 faces of
!    the icosahedron. To avoid counting pixels on edges more than
!    once, the edge pixels are split up half-and-half between the
!    two faces to which the edge belongs, so each face in fact
!    only contains 2*resolution*(resolution-1) pixels if you ignore the corners.
!    The 12 corner pixels aren't lumped in with any face at all,
!    so you can see them listed separately as the last 12 pixels
!    in test.dat if you run the demo.
!    This makes 40*resolution*(resolution-1) + 12 pixels all in all.
!    Thanks to Christopher Weth for catching typos in an earlier version
!    of this documentation!
!
!  FEEDBACK:
!    If you have any questions, comments or suggestions regarding
!    this package, please email them to me at max@ias.edu.
!    Thanks,
!    ;-)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module icosahedron
 use io, only:error
 implicit none
 public :: demo,vector2pixel,pixel2vector,compute_matrices,compute_corners

 private

contains

subroutine demo()
 real               :: vector(3), R(0:19,3,3), v(0:11,3)
 integer            :: resolution, i, j, n, pixel
 character(len=60) :: f
!      resolution = 1      ! 0 pixels per face,  so  12 pixels in total.
!      resolution = 2      ! 4 pixels per face,  so  92 pixels in total.
!      resolution = 3      ! 12 pixels per face, so 252 pixels in total.
!      resolution = 4      ! 24 pixels per face, so 492 pixels in total.
 print *,'Resolution? (1, 2, 3, ... - try say 4)'
 read *,resolution
 n = 2*resolution*(resolution-1)
 n = 20*n + 12
 call compute_matrices(R)
 call compute_corners(v)
 f = 'test.dat'
 open(2,file=f)
 do i=0,n-1
    call pixel2vector(i,resolution,R,v,vector)
    call vector2pixel(vector,resolution,R,v,pixel)
    write(2,'(2i6,3f9.5)') i, pixel, (vector(j),j=1,3)
 enddo
 close(2)
 print *,n,' pixels saved in the file ',f
 return
end subroutine demo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     THESE SUBROUTINES ARE ALL YOU NEED TO CALL FROM OUTSIDE    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      These subroutines convert between unit vectors and         !!!
!!!     pixel numbers, and are the only ones that the user of this !!!
!!!     package calls repeatedly:                           !!!
!!!        subroutine vector2pixel(vector,resolution,R,v,pixel)     !!!
!!!        subroutine pixel2vector(pixel,resolution,R,v,vector)     !!!
!!!                                                                !!!
!!!      These subroutines are called only once, in the beginning,  !!!
!!!     and compute the necessary rotation matrices and corner     !!!
!!!     vectors once and for all:                                  !!!
!!!        subroutine compute_matrices(R)                           !!!
!!!        subroutine compute_corners(v)                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine vector2pixel(vector,resolution,R,v,pixel)
 real,    intent(in)  :: vector(3)
 integer, intent(in)  :: resolution
 integer, intent(out) :: pixel
 real,    intent(in)  :: R(0:19,3,3), v(0:11,3)
 real    :: A(3,3), vec(3), x, y
 integer :: pix, face, pixperface, ifail

 if (resolution < 1) call error('vector2pixel','Resolution must exceed 0')
 pixperface = 2*resolution*(resolution-1)
 call find_face(vector,R,face)
 call getmatrix(face,R,A)
 call vecmatmul2(A,vector,vec)
 x      = vec(1)/vec(3)
 y      = vec(2)/vec(3)
 call adjust(x,y)
 call tangentplanepixel(resolution,x,y,pix,ifail)
 if (ifail > 0) then
    ! Try the runner-up face:
    call find_another_face(vector,R,face)
    call getmatrix(face,R,A)
    call vecmatmul2(A,vector,vec)
    x      = vec(1)/vec(3)
    y      = vec(2)/vec(3)
    call adjust(x,y)
    call tangentplanepixel(resolution,x,y,pix,ifail)
 endif
 pixel = face*pixperface + pix
 if (ifail > 0) then
    ! The pixel wasn't in any of those two faces,
    ! so it must be a corner pixel.
    call find_corner(vector,v,pix)
    pixel = 20*pixperface + pix
 endif

 return
end subroutine vector2pixel

subroutine pixel2vector(pixel,resolution,R,v,vector)
 ! Returns a unit vector pointing towards pixel.
 ! Resolution must be an even, positive integer.
 integer, intent(in)  :: pixel, resolution
 real,    intent(in)  :: R(0:19,3,3), v(0:11,3)
 real,    intent(out) :: vector(3)
 real    :: A(3,3), x, y, norm
 integer :: pix, face, pixperface

 if (resolution < 1) call error('pixel2vector','Resolution must exceed 0')
 pixperface = 2*resolution*(resolution-1)
 if (pixel < 0) call error('pixel2vector','negative pixel number')
 if (pixel >= 20*pixperface+12) call error('pixel2vector','pixel number too large')
 if (pixperface > 0) then
    face = pixel/pixperface
    if (face > 20) face = 20
 else ! There are no pixels at all on the faces - just corners.
    face = 20
 endif
 pix = pixel - face*pixperface
 if (face < 20) then
    ! The pixel is on one of the 20 faces:
    call tangentplanevector(pix,resolution,x,y)
    call unadjust(x,y)
    norm       = sqrt(x*x+y*y+1)
    vector(1)      = x/norm
    vector(2)      = y/norm
    vector(3)      = 1./norm
    call getmatrix(face,R,A)
    call vecmatmul1(A,vector,vector)
 else
    ! This is a corner pixel:
    if (pix > 11) call error('pixel2vector','pixel number too big')
    vector(1) = v(pix,1)
    vector(2) = v(pix,2)
    vector(3) = v(pix,3)
 endif

 return
end subroutine pixel2vector

subroutine compute_matrices(R)
 ! On exit, R will contain the 20 rotation matrices
 ! that rotate the 20 icosahedron faces
 ! into the tangent plane
 ! (the horizontal plane with z=1).
 ! Only called once, so speed is irrelevant.
 real, intent(out) :: R(0:19,3,3)
 real    :: A(3,3), B(3,3), C(3,3), D(3,3), E(3,3)
 real    :: pi, sn, cs, ct, x
 integer :: i,j,n
 do i=1,3
    do j=1,3
       A(i,j) = 0.
       B(i,j) = 0.
       C(i,j) = 0.
       D(i,j) = 0.
    enddo
 enddo
 pi       = 4.*atan(1.)
 x       = 2.*pi/5.
 cs       = cos(x)
 sn      = sin(x)
 A(1,1)      = cs
 A(1,2)      = -sn
 A(2,1)      = sn
 A(2,2)      = cs
 A(3,3)       = 1.
 ! A rotates by 72 degrees around the z-axis.
 x             = pi/5.
 ct             = cos(x)/sin(x)
 cs            = ct/sqrt(3.)
 sn            = sqrt(1-ct*ct/3.)
 C(1,1)      = 1
 C(2,2)      = cs
 C(2,3)      = -sn
 C(3,2)      = sn
 C(3,3)       = cs
 ! C rotates around the x-axis so that the north pole
 ! ends up at the center of face 1.
 cs            = -0.5
 sn            = sqrt(3.)/2
 D(1,1)      = cs
 D(1,2)      = -sn
 D(2,1)      = sn
 D(2,2)      = cs
 D(3,3)      = 1.
 ! D rotates by 120 degrees around z-axis.
 call matmul1(C,D,E)
 call matmul2(E,C,B)      ! B = CDC^t
 ! B rotates face 1 by 120 degrees.
 do i=1,3
    do j=1,3
       E(i,j) = 0.
    enddo
    E(i,i) = 1.
 enddo      ! Now E is the identity matrix.
 call putmatrix(0,R,E)
 call matmul1(B,A,E)
 call matmul1(B,E,E)
 call putmatrix(5,R,E)
 call matmul1(E,A,E)
 call putmatrix(10,R,E)
 call matmul1(E,B,E)
 call matmul1(E,B,E)
 call matmul1(E,A,E)
 call putmatrix(15,R,E)
 do n=0,15,5
    call getmatrix(n,R,E)
    do i=1,4
       call matmul1(A,E,E)
       call putmatrix(n+i,R,E)
    enddo
 enddo
 ! Now the nth matrix in R will rotate
 ! face 1 into face n.
 ! Multiply by C so that they will rotate
 ! the tangent plane into face n instead:
 do n=0,19
    call getmatrix(n,R,E)
    call matmul1(E,C,E)
    call putmatrix(n,R,E)
 enddo
 return
end subroutine compute_matrices

subroutine compute_corners(v)
 ! On exit, v will contain unit vectors pointing toward
 ! the 12 icoshedron corners.
 real, intent(out) :: v(0:11,3)
 real    :: pi, z, rho, dphi
 integer :: i
 pi = 4.*atan(1.)
 dphi = 2.*pi/5.
 ! First corner is at north pole:
 v(0,1) = 0.
 v(0,2) = 0.
 v(0,3) = 1.
 ! The next five lie on a circle, with one on the y-axis:
 z = 0.447213595      ! This is 1/(2 sin^2(pi/5)) - 1
 rho = sqrt(1.-z*z)
 do i=0,4
    v(1+i,1) = -rho*sin(i*dphi)
    v(1+i,2) =  rho*cos(i*dphi)
    v(1+i,3) = z
 enddo
 ! The 2nd half are simply opposite the first half:
 do i=0,5
    v(6+i,1) = -v(i,1)
    v(6+i,2) = -v(i,2)
    v(6+i,3) = -v(i,3)
 enddo
 return
end subroutine compute_corners

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     THE SUBROUTINES BELOW ARE SUBORDINATE TO THOSE ABOVE, AND  !!!
!!!     CAN BE SAFELY IGNORED BY THE GENERAL USER OF THE PACKAGE.  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      These subroutines perform some standard linear algebra:    !!!
!!!        subroutine matmul1(A,B,C)                                !!!
!!!        subroutine matmul2(A,B,C)                                !!!
!!!        subroutine matmul3(A,B,C)                                !!!
!!!        subroutine vecmatmul1(A,b,c)                                 !!!
!!!        subroutine vecmatmul2(A,b,c)                                 !!!
!!!                                                                !!!
!!!      These subroutines copy matrices in and out of storage:     !!!
!!!        subroutine getmatrix(n,R,A)                              !!!
!!!        subroutine putmatrix(n,R,A)                              !!!
!!!                                                                !!!
!!!     These subroutines help vector2pixel reduce the 3D sphere   !!!
!!!     problem to a problem on an equilateral triangle in the     !!!
!!!     z=1 tangent plane (an icosahedron face):                   !!!
!!!        subroutine find_face(vector,R,face)                      !!!
!!!        subroutine find_another_face(vector,R,face)              !!!
!!!        subroutine find_corner(vector,v,corner)                  !!!
!!!                                                                !!!
!!!     These subroutines pixelize this triangle with a regular    !!!
!!!     triangular grid:                                           !!!
!!!        subroutine find_mn(pixel,resolution,m,n)                 !!!
!!!        subroutine tangentplanepixel(resolution,x,y,pix,ifail)   !!!
!!!        subroutine tangentplanevector(pix,resolution,x,y)        !!!
!!!                                                                !!!
!!!     These subroutines reduce the area equalization problem to  !!!
!!!     one on the right triangle in the lower right corner:       !!!
!!!        subroutine find_sixth(x,y,rot,flip)                      !!!
!!!        subroutine rotate_and_flip(rot,flip,x,y)                 !!!
!!!        subroutine adjust(x,y)                                   !!!
!!!        subroutine unadjust(x,y)                                 !!!
!!!                                                                !!!
!!!     These subroutines perform the area equalization mappings   !!!
!!!     on this right triangle:                                    !!!
!!!        subroutine adjust_sixth(x,y)                             !!!
!!!        subroutine unadjust_sixth(x,y)                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matmul1(A,B,C)
 ! Matrix multiplication C = AB.
 ! A, B and C are allowed to be physiclly the same.
 real, intent(in) :: A(3,3), B(3,3)
 real    :: C(3,3)  ! do not declare intent to prevent compiler warnings
 real    :: D(3,3), sum
 integer :: i,j,k

 sum = 0.
 do i=1,3
    do j=1,3
       sum = 0.
       do k=1,3
          sum = sum + A(i,k)*B(k,j)
       enddo
       D(i,j) = sum
    enddo
 enddo
 call copymatrix(D,C)

 return
end subroutine matmul1

subroutine matmul2(A,B,C)
 ! Matrix multiplication C = AB^t
 ! A, B and C are allowed to be physically the same.
 real, intent(in)  :: A(3,3), B(3,3)
 real, intent(out) :: C(3,3)
 real    :: D(3,3), sum
 integer :: i,j,k

 sum = 0.
 do i=1,3
    do j=1,3
       sum = 0.
       do k=1,3
          sum = sum + A(i,k)*B(j,k)
       enddo
       D(i,j) = sum
    enddo
 enddo
 call copymatrix(D,C)

 return
end subroutine matmul2

!subroutine matmul3(A,B,C)
 ! Matrix multiplication C = A^t B
 ! A, B and C are allowed to be physically the same.
! real, intent(in)  :: A(3,3), B(3,3)
! real, intent(out) :: C(3,3)
! real    :: D(3,3), sum
! integer :: i,j,k

! sum = 0.
! do i=1,3
!   do j=1,3
!     sum = 0.
!     do k=1,3
!       sum = sum + A(k,i)*B(k,j)
!     enddo
!     D(i,j) = sum
!   enddo
! enddo
! call copymatrix(D,C)
!
! return
!end subroutine matmul3

subroutine vecmatmul1(A,b,c)
 ! Matrix multiplication c = Ab
 ! b and c are allowed to be physically the same.
 real, intent(in) :: A(3,3), b(3)
 real    :: c(3)  ! intent(out) but do not declare to prevent warnings
 real    :: d(3), sum
 integer :: i,j

 sum = 0.
 do i=1,3
    sum = 0.
    do j=1,3
       sum = sum + A(i,j)*b(j)
    enddo
    d(i) = sum
 enddo
 call copyvector(d,c)

 return
end subroutine vecmatmul1

subroutine vecmatmul2(A,b,c)
 ! Matrix multiplication c = A^tb
 ! b and c are allowed to be physiclly the same.
 real, intent(in)  :: A(3,3), b(3)
 real, intent(out) :: c(3)
 real    :: d(3), sum
 integer :: i,j

 sum = 0.
 do i=1,3
    sum = 0.
    do j=1,3
       sum = sum + A(j,i)*b(j)
    enddo
    d(i) = sum
 enddo
 call copyvector(d,c)

 return
end subroutine vecmatmul2

subroutine copymatrix(A,B)
 ! B = A
 real, intent(in)  :: A(3,3)
 real, intent(out) :: B(3,3)
 integer :: i,j

 do i=1,3
    do j=1,3
       B(i,j) = A(i,j)
    enddo
 enddo

 return
end subroutine copymatrix

subroutine copyvector(a,b)
 ! b = a
 real, intent(in)  :: a(3)
 real, intent(out) :: b(3)
 integer :: i

 do i=1,3
    b(i) = a(i)
 enddo

 return
end subroutine copyvector

subroutine getmatrix(n,R,A)
 ! A = the nth matrix in R
 integer, intent(in)  :: n
 real,    intent(in)  :: R(0:19,3,3)
 real,    intent(out) :: A(3,3)
 integer :: i,j

 do i=1,3
    do j=1,3
       A(i,j) = R(n,i,j)
    enddo
 enddo

 return
end subroutine getmatrix

subroutine putmatrix(n,R,A)
 ! the nth matrix in R = A
 integer, intent(in)    :: n
 real,    intent(inout) :: R(0:19,3,3)
 real,    intent(in)    :: A(3,3)
 integer :: i,j

 do i=1,3
    do j=1,3
       R(n,i,j) = A(i,j)
    enddo
 enddo

 return
end subroutine putmatrix

subroutine find_face(vector,R,face)
 ! Locates the face to which vector points.
 ! Computes the dot product with the vectors
 ! pointing to the center of each face and picks the
 ! largest one.
 ! This simple routine can be substantially accelerated
 ! by adding a bunch of if-statements, to avoid looping
 ! over more than a few faces.
 real,    intent(in)  :: vector(3), R(0:19,3,3)
 integer, intent(out) :: face
 real    :: dot, max
 integer :: n,i

 max  = -17.
 face = 0 ! to avoid compiler warnings
 do n=0,19
    dot = 0.
    do i=1,3
       dot = dot + R(n,i,3)*vector(i)
    enddo
    if (dot > max) then
       face = n
       max = dot
    endif
 enddo

 return
end subroutine find_face

subroutine find_another_face(vector,R,face)
 ! Computes the dot product with the vectors
 ! pointing to the center of each face and picks the
 ! largest one other than face.
 ! This simple routine can be substantially accelerated
 ! by adding a bunch of if-statements, to avoid looping
 ! over more than a few faces.
 real,    intent(in)  :: vector(3), R(0:19,3,3)
 integer, intent(out) :: face
 real    :: dot, max
 integer :: n,facetoavoid,i

 facetoavoid = face
 max = -17.
 do n=0,19
    if (n /= facetoavoid) then
       dot = 0.
       do i=1,3
          dot = dot + R(n,i,3)*vector(i)
       enddo
       if (dot > max) then
          face = n
          max = dot
       endif
    endif
 enddo

 return
end subroutine find_another_face

subroutine find_corner(vector,v,corner)
 ! Locates the corner to which vector points.
 ! Computes the dot product with the vectors
 ! pointing to each corner and picks the
 ! largest one.
 ! This simple routine can be substantially accelerated
 ! by adding a bunch of if-statements, but that's pretty
 ! pointless since it gets called so rarely.
 real,    intent(in)  :: vector(3), v(0:11,3)
 integer, intent(out) :: corner
 real    :: dot, max
 integer :: n,i

 max = -17.
 do n=0,11
    dot = 0.
    do i=1,3
       dot = dot + v(n,i)*vector(i)
    enddo
    if (dot > max) then
       corner = n
       max = dot
    endif
 enddo

 return
end subroutine find_corner

subroutine find_mn(pixel,resolution,m,n)
 ! Computes the integer coordinates (m,n) of the pixel
 ! numbered pix on the basic triangle.
 integer, intent(in)  :: pixel, resolution
 integer, intent(out) :: m,n
 integer :: pix, interiorpix , pixperedge
 pix           = pixel
 interiorpix = (2*resolution-3)*(resolution-1)
 pixperedge  = (resolution)-1
 if (pix < interiorpix) then
    ! The pixel lies in the interior of the triangle.
    m = int((sqrt(1.+8.*pix)-1.)/2. + 0.5/resolution)
    ! 0.5/resolution was added to avoid problems with
    ! rounding errors for the case when n=0.
    ! As long as you don't add more than 2/m, you're OK.
    n = pix - m*(m+1)/2
    m = m + 2
    n = n + 1
    return
 endif
 pix = pix - interiorpix
 if (pix < pixperedge) then
    ! The pixel lies on the bottom edge.
    m = 2*resolution-1
    n = pix+1
    return
 endif
 pix = pix - pixperedge
 if (pix < pixperedge) then
    ! The pixel lies on the right edge.
    m = 2*resolution-(pix+2)
    n = m
    return
 endif
 pix = pix - pixperedge
 ! The pixel lies on the left edge.
 m = pix+1
 n = 0

 return
end subroutine find_mn

subroutine tangentplanepixel(resolution,x,y,pix,ifail)
 ! Finds the hexagon in which the point (x,y) lies
 ! and computes the corresponding pixel number pix.
 ! Returns ifail=0 if (x,y) lies on the face,
 ! otherwise returns ifail=1.
 integer, intent(in)  :: resolution
 real,    intent(in)  :: x,y
 integer, intent(out) :: pix, ifail
 real, parameter :: c = 0.866025404     ! sqrt(3)/2
 real, parameter :: edgelength = 1.3231690765
 ! The edge length of the icosahedron is
 ! sqrt(9 tan^2(pi/5) - 3) when scaled so that
 ! it circumscribes the unit sphere.
 real    :: a, b, d
 integer :: i, j, k, m, n, r2

 r2      = 2*resolution
 a       = 0.5*x
 b       = c*y
 d       = 0.5*edgelength/r2
 i       = int(x/d) + r2
 j       = int((a+b)/d) + r2
 k       = int((a-b)/d) + r2
 m       = (r2+r2-j+k-1)/3
 n       = (i+k+1-r2)/3
 pix     = (m-2)*(m-1)/2 + (n-1)
 ifail = 0
 if (m==r2-1) then            ! On bottom row
    if ((n <= 0).or.(n >= resolution)) then
       ifail=1
    endif
    return                      ! Pix already correct
 endif
 if (n==m) then                  ! On right edge
    k = (r2-1) - m
    if ((k <= 0).or.(k >= resolution)) then
       ifail = 1
    else
       pix = (r2-2)*(resolution-1) + k - 1
    endif
    return
 endif
 if (n==0) then                  ! On left edge
    if ((m <= 0).or.(m >= resolution)) then
       ifail = 1
    else
       pix = (r2-1)*(resolution-1) + m - 1
    endif
 endif

 return
end subroutine tangentplanepixel

subroutine tangentplanevector(pix,resolution,x,y)
 ! Computes the coordinates (x,y) of the pixel
 ! numbered pix on the basic triangle.
 integer, intent(in)  :: pix,resolution
 real,    intent(out) :: x,y
 real, parameter :: c1 = 0.577350269      ! 1/sqrt(3)
 real, parameter :: c2 = 0.866025404      ! sqrt(3)/2
 real, parameter :: edgelength = 1.3231690765
 ! The edge length of the icosahedron is
 ! sqrt(9 tan^2(pi/5) - 3) when scaled so that
 ! it circumscribes the unit sphere.
 integer :: m, n

 call find_mn(pix,resolution,m,n)
 x = edgelength*(n-0.5*m)/(2*resolution-1)
 y = edgelength*(c1-(c2/(2*resolution-1))*m)

 return
end subroutine tangentplanevector

subroutine find_sixth(x,y,rot,flip)
 ! Find out in which sixth of the basic triangle
 ! the point (x,y) lies, identified by the
 ! two integers rot (=0, 1 or 2) and flip = 0 or 1).
 ! rot and flip are defined such that the sixth is
 ! mapped onto the one at the bottom right by
 ! these two steps:
 ! 1. Rotate by 120 degrees anti-clockwise, rot times.
 ! 2. Flip the sign of x if flip = 1, not if flip=0.
 ! The if-statements below go through the six cases
 ! anti-clockwise, starting at the bottom right.
 real,    intent(in)  :: x,y
 integer, intent(out) :: rot,flip
 real, parameter :: c = 1.73205081   ! sqrt(3)
 real :: d

 d = c*y
 if (x >= 0) then
    if (x <= -d) then
       rot  = 0
       flip = 0
    else
       if (x >= d) then
          rot  = 2
          flip = 1
       else
          rot  = 2
          flip = 0
       endif
    endif
 else
    if (x >= -d) then
       rot  = 1
       flip = 1
    else
       if (x <= d) then
          rot  = 1
          flip = 0
       else
          rot  = 0
          flip = 1
       endif
    endif
 endif

 return
end subroutine find_sixth

subroutine rotate_and_flip(rot,flip,x,y)
 real,    intent(inout) :: x,y
 integer, intent(in)    :: rot,flip
 real, parameter :: cs = -0.5
 real, parameter :: c = 0.866025404      ! sqrt(3)/2
 real :: x1, sn

 if (rot > 0) then
    if (rot==1) then
       sn = c      ! Rotate 120 degrees anti-clockwise
    else
       sn = -c      ! Rotate 120 degrees anti-clockwise
    endif
    x1 = x
    x  = cs*x1 - sn*y
    y  = sn*x1 + cs*y
 endif
 if (flip > 0) x = -x

 return
end subroutine rotate_and_flip

subroutine adjust(x,y)
 ! Maps the basic triangle onto itself in such a way
 ! that pixels will have equal area when mapped onto
 ! the sphere.
 real, intent(inout) :: x, y
 integer :: rot, flip

 call find_sixth(x,y,rot,flip)
 call rotate_and_flip(rot,flip,x,y)
 call adjust_sixth(x,y)
 ! Now rotate & flip the sixth back into its
 ! original position:
 if ((flip==0).and.(rot > 0)) then
    call rotate_and_flip(3-rot,flip,x,y)
 else
    call rotate_and_flip(rot,flip,x,y)
 endif

 return
end subroutine adjust

subroutine unadjust(x,y)
 ! Performs the inverse of what adjust does.
 real, intent(inout) :: x, y
 integer :: rot, flip

 call find_sixth(x,y,rot,flip)
 call rotate_and_flip(rot,flip,x,y)
 call unadjust_sixth(x,y)
 ! Now rotate & flip the sixth back into its
 ! original position:
 if ((flip==0).and.(rot > 0)) then
    call rotate_and_flip(3-rot,flip,x,y)
 else
    call rotate_and_flip(rot,flip,x,y)
 endif

 return
end subroutine unadjust

subroutine adjust_sixth(x,y)
 ! Maps the basic right triangle (the sixth of the face that
 ! is in the lower right corner) onto itself in such a way
 ! that pixels will have equal area when mapped onto the sphere.
 real, intent(inout) :: x, y
 real, parameter :: eps = 1.e-14
 real, parameter :: scale = 1.09844
 real, parameter :: g = 1.7320508075689    ! sqrt(3)
 real :: u, v, v2, root, trig

 u     = x  + eps
 v     = -y + eps
 v2    = v*v
 root  = sqrt(1.+4.*v2)
 trig  = atan((g*root-g)/(root+3.))
 y     = sqrt(trig*2./g)
 x     = sqrt((1.+4.*v2)/(1.+u*u+v2))*u*y/v
 x     = scale*x
 y     = -scale*y

 return
end subroutine adjust_sixth

subroutine unadjust_sixth(x,y)
 ! Performs the inverse of what adjust_sixth does.
 real, intent(inout) :: x, y
 real, parameter :: eps=1.e-14
 real, parameter :: scale=1.09844
 real, parameter :: g=1.7320508075689   ! sqrt(3)
 real :: u, v, v2, y2, tmp, trig

 u    =  x/scale + eps
 v    = -y/scale + eps
 v2   = v*v
 trig = tan(g*v2/2.)
 tmp  = (g+3.*trig)/(g-trig)
 y2   = (tmp*tmp-1.)/4.
 y    = sqrt(y2)
 tmp  = v2*(1.+4.*y2) - u*u*y2
 x    = u*y*sqrt((1.+y2)/tmp)
 y    = -y

 return
end subroutine unadjust_sixth

end module icosahedron
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     END OF THE ICOSAHEDRON PACKAGE FOR PIXELIZING THE SPHERE   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
