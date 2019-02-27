! This file is part of RMPCDMD
! Copyright (c) 2015-2016 Pierre de Buyl and contributors
! License: BSD 3-clause (see file LICENSE)

!> Compute compact Hilbert indices
!!
!! Ref. \cite hamilton_compact_hilbert_tr
!
! Modifications by Daniel Price (Feb 2019):
!  Converted precision to int*8 instead of int
!
module hilbert
  implicit none
  private

  public :: p_to_h, h_to_p, compact_p_to_h, compact_h_to_p

  integer(8), parameter :: dim = 3
  integer(8), parameter :: mask=2**dim-1

contains

  pure function rotate_right(x, d)
    integer(8), intent(in) :: x
    integer(8), intent(in) :: d
    integer(8) :: rotate_right
    integer(8) :: tmp

    rotate_right = shiftr(x, d)
    tmp = shiftl(x, dim-d)
    rotate_right = iand(ior(rotate_right, tmp), mask)

  end function rotate_right

  pure function rotate_left(x, d)
    integer(8), intent(in) :: x
    integer(8), intent(in) :: d
    integer(8) :: rotate_left
    integer(8) :: tmp

    rotate_left = shiftl(x, d)
    tmp = shiftr(x, dim-d)
    rotate_left = iand(ior(rotate_left, tmp), mask)

  end function rotate_left

  pure function bin_str(x)
    integer(8), intent(in) :: x

    character(len=dim) :: bin_str
    integer(8) :: i, j, total_size

    total_size = bit_size(x)

    do i=1, dim
       j = dim-(i-1)
       if (ibits(x, i-1, 1).eq.1) then
          bin_str(j:j) = '1'
       else
          bin_str(j:j) = '0'
       end if
    end do

  end function bin_str

  pure function gc(i)
    integer(8), intent(in) :: i
    integer(8) gc

    gc = ieor(i, shiftr(i, 1))

  end function gc

  pure function entry_point(i)
    integer(8), intent(in) :: i
    integer(8) :: entry_point

    if (i .eq. 0) then
       entry_point = 0
    else
       entry_point = gc(2*((i-1)/2))
    end if

  end function entry_point

  pure function exit_point(i)
    integer(8), intent(in) :: i
    integer(8) :: exit_point

    exit_point = ieor(entry_point(mask-i), 2**(dim-1))

  end function exit_point

  pure function inverse_gc(g)
    integer(8), intent(in) :: g
    integer(8) :: inverse_gc
    integer(8) :: j

    inverse_gc = g
    j = 1
    do while ( j .lt. dim )
       inverse_gc = ieor(inverse_gc, shiftr(g, j))
       j = j + 1
    end do

  end function inverse_gc

  pure function intercube_g(i) result(g)
    integer(8), intent(in) :: i
    integer(8) :: g

    g = trailz(ieor(gc(i), gc(i+1)))

  end function intercube_g

  pure function intracube_d(i) result(d)
    integer(8), intent(in) :: i
    integer(8) :: d

    if (i .eq. 0) then
       d = 0
    else if ( modulo(i, 2_8) .eq. 0 ) then
       d = modulo(intercube_g(i-1), dim)
    else
       d = modulo(intercube_g(i), dim)
    end if
  end function intracube_d

  pure function transform(e, d, b) result(t)
    integer(8), intent(in) :: e, d, b
    integer(8) :: t

    t = rotate_right( ieor(b, e), d+1)
  end function transform

  pure function inverse_transform(e, d, b) result(t)
    integer(8), intent(in) :: e, d, b
    integer(8) :: t

    t = transform(rotate_right(e, d+1), dim-d-2, b)

  end function inverse_transform

  pure function p_to_h(p, m) result(h)
    integer(8), intent(in) :: p(dim), m
    integer(8) :: h

    integer(8) :: e, d, i, j, l, w

    h = 0
    e = 0
    d = 2
    do i = m-1, 0, -1
       l = 0
       do j=1, dim
          l = l + 2**(j-1)*ibits(p(j), i, 1)
       end do
       l = transform(e, d, l)
       w = inverse_gc(l)
       e = ieor(e, rotate_left(entry_point(w), d+1))
       d = modulo(d + intracube_d(w) + 1, dim)
       h = ior(shiftl(h, dim), w)
    end do

  end function p_to_h

  pure function h_to_p(h, m) result(p)
    integer(8), intent(in) :: h, m
    integer(8) :: p(dim)

    integer(8) :: e, d, i, j, l, w

    e = 0
    d = 2
    p = 0
    do i=m-1, 0, -1
       w = 0
       do j=0, dim-1
          w = w + 2**j*ibits(h, i*dim+j, 1)
       end do
       l = gc(w)
       l = inverse_transform(e, d, l)
       do j=1, dim
          p(j) = p(j) + shiftl(ibits(l, j-1, 1) , i)
       end do
       e = ieor( e, rotate_left(entry_point(w), d+1))
       d = modulo(d + intracube_d(w) + 1, dim)
    end do

  end function h_to_p

  pure function gcr(i, mu) result(r)
    integer(8), intent(in) :: i, mu
    integer(8) :: r

    integer(8) :: k
    r = 0
    do k=dim-1, 0, -1
       if (ibits(mu, k, 1) .eq. 1) then
          r = ior( shiftl(r, 1), ibits(i, k, 1))
       end if
    end do

  end function gcr

  pure function inverse_gcr(r, mu, pi) result(i)
    integer(8), intent(in) :: r, mu, pi
    integer(8) :: i

    integer(8) :: g, j, k
    i = 0_8
    g = 0_8
    j = -1_8
    do k=0, dim - 1
       j = j + ibits(mu, k, 1_8)
    end do
    do k=dim-1, 0, -1
       if (ibits(mu, k, 1) .eq. 1) then
          i = ior(i, shiftl( ibits(r, j, 1_8), k))
          g = ior(g, shiftl( modulo( ibits(i, k, 1_8)+ibits(i, k+1_8, 1_8), 2_8), k) )
          j = j-1
       else
          g = ior(g, shiftl( ibits(pi, k, 1_8), k))
          i = ior(i, shiftl( modulo( ibits(g, k, 1_8)+ibits(i, k+1, 1_8), 2_8), k) )
       end if
    end do
  end function inverse_gcr

  pure function extract_mask(i, m) result(mu)
    integer(8), intent(in) :: i, m(dim)
    integer(8) :: mu

    integer(8) :: j

    mu = 0
    do j=dim, 1, -1
       mu = shiftl(mu, 1_8)
       if ( m(j) .gt. i ) then
          mu = ior(mu, 1_8)
       end if
    end do

  end function extract_mask

  pure function compact_p_to_h(p, m) result(h)
    integer(8), intent(in) :: p(dim), m(dim)
    integer(8) :: h

    integer(8) :: e, d, max_m, i, j, l, mu, mu_norm, pi, r, w

    h = 0
    e = 0
    d = 2
    max_m = maxval(m)
    do i=max_m-1, 0, -1
       mu = extract_mask(i, m)
       mu_norm = 0
       do j=0, dim-1
          mu_norm = mu_norm + ibits(mu, j, 1)
       end do
       mu = rotate_right(mu, d+1)
       pi = iand(rotate_right(e, d+1), iand(not(mu), mask))
       l = 0
       do j=1, dim
          l = l + 2**(j-1)*ibits(p(j), i, 1)
       end do
       l = transform(e, d, l)
       w = inverse_gc(l)
       r = gcr(w, mu)
       e = ieor(e, rotate_left(entry_point(w), d+1))
       d = modulo(d + intracube_d(w) + 1, dim)
       h = ior( shiftl(h, mu_norm), r)
    end do

  end function compact_p_to_h

  pure function compact_h_to_p(h, m) result(p)
    integer(8), intent(in) :: h, m(dim)
    integer(8) :: p(dim)

    integer(8) :: e, d, i, j, k, max_m, sum_m, mu_norm, mu, pi, r, l, w

    e = 0
    d = 2
    k = 0
    p = 0
    max_m = maxval(m)
    sum_m = sum(m)
    do i=max_m-1, 0, -1
       mu = extract_mask(i, m)
       mu_norm = 0
       do j=0, dim-1
          mu_norm = mu_norm + ibits(mu, j, 1)
       end do
       mu = rotate_right(mu, d+1)
       pi = iand(rotate_right(e, d+1), iand(not(mu), mask))
       r = 0
       do j=0, mu_norm-1
          r = r + 2**j*ibits(h, sum_m - k - mu_norm + j, 1)
       end do
       k = k + mu_norm
       w = inverse_gcr(r, mu, pi)
       l = gc(w)
       l = inverse_transform(e, d, l)
       do j=1, dim
          p(j) = p(j) + shiftl(ibits(l, j-1, 1), i)
       end do
       e = ieor(e, rotate_left(entry_point(w), d+1))
       d = modulo(d + intracube_d(w) + 1, dim)
    end do
  end function compact_h_to_p

end module hilbert
