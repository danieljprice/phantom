!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module hilbert
!
! Peano-Hilbert keys in three dimensions, used to order particles along
! a space filling curve for the MPI domain decomposition.
!
! A key of m bits per dimension has 3m bits. The first 3k bits give the
! cell at level k of an octree of the unit cube, so every octree cell is a
! contiguous range of keys. Since the 8 sub-cubes of a cell are visited in
! Gray code order, a prefix of any number of bits is also an axis-aligned
! box (a half or a quarter of a cell).
!
! The routines p_to_h, h_to_p and their helpers are adapted from RMPCDMD:
!
!   Copyright (c) 2015-2016 Pierre de Buyl and contributors
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are
!   met: (1) Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer. (2)
!   Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
!   (3) Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
!   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! :References: Hamilton (2006), Compact Hilbert indices, Technical Report
!              CS-2006-07, Dalhousie University
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 integer(8), parameter, public :: nbits_key = 21  ! bits per dimension: 63-bit keys

 public :: p_to_h,h_to_p,get_hilbert_keys

 private

 integer(8), parameter :: ndim = 3
 integer(8), parameter :: mask = 2**ndim - 1

contains

!--------------------------------------------------------------------------
!+
!  Hilbert keys of the particles, with nbits_key bits per dimension, in
!  the cube with lower corner xmin and side dxcube. Particles outside the
!  cube are given the key of the nearest cell on its surface
!+
!--------------------------------------------------------------------------
subroutine get_hilbert_keys(np,xyzh,xmin,dxcube,key)
 integer,    intent(in)  :: np
 real,       intent(in)  :: xyzh(:,:),xmin(3),dxcube
 integer(8), intent(out) :: key(:)
 integer(8) :: ip(3),ngrid
 real    :: fac
 integer :: i

 ngrid = 2_8**nbits_key
 fac   = real(ngrid)/dxcube
 !$omp parallel do default(none) schedule(static) shared(np,xyzh,xmin,fac,ngrid,key) private(i,ip)
 do i = 1,np
    ip = int((xyzh(1:3,i) - xmin)*fac,kind=8)
    ip = min(max(ip,0_8),ngrid-1)   ! a particle on the upper face would wrap to 0
    key(i) = p_to_h(ip,nbits_key)
 enddo
 !$omp end parallel do

end subroutine get_hilbert_keys

!--------------------------------------------------------------------------
!+
!  Hilbert key of the integer coordinates p, with m bits per dimension
!+
!--------------------------------------------------------------------------
pure function p_to_h(p,m) result(h)
 integer(8), intent(in) :: p(ndim),m
 integer(8) :: h
 integer(8) :: e,d,i,j,l,w

 h = 0
 e = 0
 d = 2
 do i = m-1,0,-1
    l = 0
    do j = 1,ndim
       l = l + 2**(j-1)*ibits(p(j),i,1)
    enddo
    l = transform(e,d,l)
    w = inverse_gc(l)
    e = ieor(e,rotate_left(entry_point(w),d+1))
    d = modulo(d + intracube_d(w) + 1,ndim)
    h = ior(shiftl(h,ndim),w)
 enddo

end function p_to_h

!--------------------------------------------------------------------------
!+
!  integer coordinates of the Hilbert key h, with m bits per dimension
!+
!--------------------------------------------------------------------------
pure function h_to_p(h,m) result(p)
 integer(8), intent(in) :: h,m
 integer(8) :: p(ndim)
 integer(8) :: e,d,i,j,l,w

 e = 0
 d = 2
 p = 0
 do i = m-1,0,-1
    w = 0
    do j = 0,ndim-1
       w = w + 2**j*ibits(h,i*ndim+j,1)
    enddo
    l = gc(w)
    l = inverse_transform(e,d,l)
    do j = 1,ndim
       p(j) = p(j) + shiftl(ibits(l,j-1,1),i)
    enddo
    e = ieor(e,rotate_left(entry_point(w),d+1))
    d = modulo(d + intracube_d(w) + 1,ndim)
 enddo

end function h_to_p

!--------------------------------------------------------------------------
!+
!  helper functions: rotations within ndim bits, Gray code and its
!  inverse, and the entry point and direction of each sub-cube
!+
!--------------------------------------------------------------------------
pure function rotate_right(x,d) result(r)
 integer(8), intent(in) :: x,d
 integer(8) :: r

 r = iand(ior(shiftr(x,d),shiftl(x,ndim-d)),mask)

end function rotate_right

pure function rotate_left(x,d) result(r)
 integer(8), intent(in) :: x,d
 integer(8) :: r

 r = iand(ior(shiftl(x,d),shiftr(x,ndim-d)),mask)

end function rotate_left

pure function gc(i)
 integer(8), intent(in) :: i
 integer(8) :: gc

 gc = ieor(i,shiftr(i,1))

end function gc

pure function inverse_gc(g)
 integer(8), intent(in) :: g
 integer(8) :: inverse_gc
 integer(8) :: j

 inverse_gc = g
 do j = 1,ndim-1
    inverse_gc = ieor(inverse_gc,shiftr(g,j))
 enddo

end function inverse_gc

pure function entry_point(i)
 integer(8), intent(in) :: i
 integer(8) :: entry_point

 if (i == 0) then
    entry_point = 0
 else
    entry_point = gc(2*((i-1)/2))
 endif

end function entry_point

pure function intercube_g(i) result(g)
 integer(8), intent(in) :: i
 integer(8) :: g

 g = trailz(ieor(gc(i),gc(i+1)))

end function intercube_g

pure function intracube_d(i) result(d)
 integer(8), intent(in) :: i
 integer(8) :: d

 if (i == 0) then
    d = 0
 elseif (modulo(i,2_8) == 0) then
    d = modulo(intercube_g(i-1),ndim)
 else
    d = modulo(intercube_g(i),ndim)
 endif

end function intracube_d

pure function transform(e,d,b) result(t)
 integer(8), intent(in) :: e,d,b
 integer(8) :: t

 t = rotate_right(ieor(b,e),d+1)

end function transform

pure function inverse_transform(e,d,b) result(t)
 integer(8), intent(in) :: e,d,b
 integer(8) :: t

 t = transform(rotate_right(e,d+1),ndim-d-2,b)

end function inverse_transform

end module hilbert
