!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module healpix
!
! This module sets the types used in the Fortran 90 modules (healpix_types.f90)
! of the HEALPIX distribution and follows the example of Numerical Recipes
! Benjamin D. Wandelt October 1997
! Eric Hivon June 1998
! Eric Hivon Oct  2001, edited to be compatible with 'F' compiler
! Eric Hivon July 2002, addition of i8b, i2b, i1b
!                       addition of max_i8b, max_i2b and max_i1b
!            Jan 2005, explicit form of max_i1b because of ifc 8.1.021
!            June 2005, redefine i8b as 16 digit integer because of Nec f90 compiler
!            Mars 2008: i8b same as i4b on machines not supporting 64 bits (NO64BITS flag set)
!            Feb  2009: introduce healpix_version
!
! :References: K. M. Górski et al, 2005, ApJ, 622, 759
!
! :Owner: Mats Esseldeurs
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 character(len=*), parameter, public :: healpix_version = '3.80'
 integer, parameter, public :: i4b = selected_int_kind(9)
 integer, parameter, public :: i8b = selected_int_kind(16)
 integer, parameter, public :: i2b = selected_int_kind(4)
 integer, parameter, public :: i1b = selected_int_kind(2)
 integer, parameter, public :: sp  = selected_real_kind(5,30)
 integer, parameter, public :: dp  = selected_real_kind(12,200)
 integer, parameter, public :: lgt = kind(.TRUE.)
 integer, parameter, public :: spc = kind((1.0_sp, 1.0_sp))
 integer, parameter, public :: dpc = kind((1.0_dp, 1.0_dp))
 !
 integer(I8B),  parameter, public :: max_i8b = huge(1_i8b)
 integer,       parameter, public :: max_i4b = huge(1_i4b)
 integer,       parameter, public :: max_i2b = huge(1_i2b)
 integer,       parameter, public :: max_i1b = 127
 real(kind=sp), parameter, public :: max_sp  = huge(1.0_sp)
 real(kind=dp), parameter, public :: max_dp  = huge(1.0_dp)

 ! Numerical Constant (Double precision)
 real(kind=dp), parameter, public :: QUARTPI=0.785398163397448309615660845819875721049_dp
 real, parameter, public :: HALFPI= 1.570796326794896619231321691639751442099
 real, parameter, public :: PI    = 3.141592653589793238462643383279502884197
 real, parameter, public :: TWOPI = 6.283185307179586476925286766559005768394
 real(kind=dp), parameter, public :: FOURPI=12.56637061435917295385057353311801153679_dp
 real(kind=dp), parameter, public :: SQRT2 = 1.41421356237309504880168872420969807856967_dp
 real(kind=dp), parameter, public :: EULER = 0.5772156649015328606065120900824024310422_dp
 real(kind=dp), parameter, public :: SQ4PI_INV = 0.2820947917738781434740397257803862929220_dp
 real(kind=dp), parameter, public :: TWOTHIRD = 0.6666666666666666666666666666666666666666_dp

 real(kind=DP), parameter, public :: RAD2DEG = 180.0_DP / PI
 real(kind=DP), parameter, public :: DEG2RAD = PI / 180.0_DP
 real(kind=SP), parameter, public :: hpx_sbadval = -1.6375e30_sp
 real(kind=DP), parameter, public :: hpx_dbadval = -1.6375e30_dp

 ! Maximum length of filenames
 integer, parameter :: filenamelen = 1024


 !   ! ---- Normalisation and convention ----
 ! normalisation of spin weighted functions
 real(kind=dp), parameter, public ::  KvS = 1.0_dp ! 1.0 : CMBFAST (Healpix 1.2)
 !   ! sign of Q
 !   real(kind=dp), parameter, public :: sgQ = -1.0_dp ! -1 : CMBFAST (Healpix 1.2)
 !   ! sign of spin weighted function !
 !   real(kind=dp), parameter, public :: SW1 = -1.0_dp ! -1 : Healpix 1.2, bug correction

 ! !  ! normalisation of spin weighted functions
 ! !  real(kind=dp), parameter, public ::  KvS = 2.0_dp ! 2.0 : KKS  (Healpix 1.1)
 ! !  ! sign of Q
 ! !  real(kind=dp), parameter, public :: sgQ = +1.0_dp ! +1 : KKS (Healpix 1.1)
 ! !  ! sign of spin weighted function !
 ! !  real(kind=dp), parameter, public :: SW1 = +1.0_dp ! +1 : Healpix 1.1

 !   real(kind=dp), parameter, public :: iKvS = 1.0_dp / KvS  ! inverse of KvS
 integer(kind=i4b), private, parameter :: ns_max4=8192     ! 2^13
 integer(kind=i4b), private, save, dimension(0:127) :: x2pix1=-1,y2pix1=-1
 integer(kind=i4b), private, save, dimension(0:1023) :: pix2x=-1, pix2y=-1
 integer(i4b), parameter :: oddbits=89478485   ! 2^0 + 2^2 + 2^4+..+2^26
 integer(i4b), parameter :: evenbits=178956970 ! 2^1 + 2^3 + 2^4+..+2^27
 integer(kind=i4b), private, parameter :: ns_max=268435456! 2^28

contains

 !! Returns i with even and odd bit positions interchanged.
function swapLSBMSB(i)
 integer(i4b) :: swapLSBMSB
 integer(i4b), intent(in) :: i

 swapLSBMSB = iand(i,evenbits)/2 + iand(i,oddbits)*2
end function swapLSBMSB

 !! Returns not(i) with even and odd bit positions interchanged.
function invswapLSBMSB(i)
 integer(i4b) :: invswapLSBMSB
 integer(i4b), intent(in) :: i

 invswapLSBMSB = not(swapLSBMSB(i))
end function invswapLSBMSB

 !! Returns i with odd (1,3,5,...) bits inverted.
function invLSB(i)
 integer(i4b) :: invLSB
 integer(i4b), intent(in) :: i

 invLSB = ieor(i,oddbits)
end function invLSB

 !! Returns i with even (0,2,4,...) bits inverted.
function invMSB(i)
 integer(i4b) :: invMSB
 integer(i4b), intent(in) :: i

 invMSB = ieor(i,evenbits)
end function invMSB

 !=======================================================================
 !     vec2pix_nest
 !
 !     renders the pixel number ipix (NESTED scheme) for a pixel which contains
 !     a point on a sphere at coordinate vector (=x,y,z), given the map
 !     resolution parameter nside
 !
 ! 2009-03-10: calculations done directly at nside rather than ns_max
 !=======================================================================
subroutine vec2pix_nest  (nside, vector, ipix)
 integer(i4b), parameter :: MKD = I4B
 integer(kind=I4B), intent(in)                :: nside
 real,              intent(in), dimension(1:) :: vector
 integer(kind=MKD), intent(out)               :: ipix

 integer(kind=MKD) :: ipf,scale,scale_factor
 real(kind=DP)     :: z,za,tt,tp,tmp,dnorm,phi
 integer(kind=I4B) :: jp,jm,ifp,ifm,face_num,ix,iy,ix_low,iy_low,ntt,i,ismax
 character(len=*), parameter :: code = "vec2pix_nest"

 !-----------------------------------------------------------------------
 if (nside<1 .or. nside>ns_max4) call fatal_error(code//"> nside out of range")
 dnorm = sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
 z = vector(3) / dnorm
 phi = 0.0
 if (vector(1) /= 0.0 .or. vector(2) /= 0.0) &
 &     phi = atan2(vector(2),vector(1)) ! phi in ]-pi,pi]

 za = abs(z)
 if (phi < 0.0)    phi = phi + twopi ! phi in [0,2pi[
 tt = phi / halfpi ! in [0,4[
 if (x2pix1(127) <= 0) call mk_xy2pix1()

 if (za <= twothird) then ! equatorial region

    !        (the index of edge lines increase when the longitude=phi goes up)
    jp = int(nside*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
    jm = int(nside*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

    !        finds the face
    ifp = jp / nside  ! in {0,4}
    ifm = jm / nside
    if (ifp == ifm) then          ! faces 4 to 7
       face_num = iand(ifp,3) + 4
    elseif (ifp < ifm) then     ! (half-)faces 0 to 3
       face_num = iand(ifp,3)
    else                            ! (half-)faces 8 to 11
       face_num = iand(ifm,3) + 8
    endif

    ix =         iand(jm, nside-1)
    iy = nside - iand(jp, nside-1) - 1

 else ! polar region, za > 2/3

    ntt = int(tt)
    if (ntt >= 4) ntt = 3
    tp = tt - ntt
    !tmp = sqrt( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]
    tmp = sqrt(vector(1)**2+vector(2)**2) / dnorm ! sin(theta)
    tmp = tmp * sqrt( 3.0_dp / (1.0_dp + za) ) !more accurate

    !        (the index of edge lines increase when distance from the closest pole goes up)
    jp = int( nside * tp          * tmp ) ! line going toward the pole as phi increases
    jm = int( nside * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
    jp = min(nside-1, jp) ! for points too close to the boundary
    jm = min(nside-1, jm)

    !        finds the face and pixel's (x,y)
    if (z >= 0) then
       face_num = ntt  ! in {0,3}
       ix = nside - jm - 1
       iy = nside - jp - 1
    else
       face_num = ntt + 8 ! in {8,11}
       ix =  jp
       iy =  jm
    endif

 endif

 if (nside <= ns_max4) then
    ix_low = iand(ix, 127)
    iy_low = iand(iy, 127)
    ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
    & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
 else
    scale = 1_MKD
    scale_factor = 16384_MKD ! 128*128
    ipf = 0_MKD
    ismax = 1 ! for nside in [2^14, 2^20]
    if (nside >  1048576 ) ismax = 3
    do i=0, ismax
       ix_low = iand(ix, 127) ! last 7 bits
       iy_low = iand(iy, 127) ! last 7 bits
       ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
       scale = scale * scale_factor
       ix  =     ix / 128 ! truncate out last 7 bits
       iy  =     iy / 128
    enddo
    ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
 endif
 ipix = ipf + face_num* int(nside,MKD) * nside    ! in {0, 12*nside**2 - 1}

end subroutine vec2pix_nest

 !=======================================================================
 !     pix2vec_nest
 !
 !     renders vector (x,y,z) coordinates of the nominal pixel center
 !     for the pixel number ipix (NESTED scheme)
 !     given the map resolution parameter nside
 !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
 !     in the order N,W,S,E
 !=======================================================================
subroutine pix2vec_nest  (nside, ipix, vector, vertex)
 integer(i4b), parameter :: MKD = i4b
 integer(kind=I4B), intent(in) :: nside
 integer(kind=MKD), intent(in) :: ipix
 real,              intent(out), dimension(1:) :: vector
 real,     intent(out), dimension(1:,1:), optional :: vertex

 integer(kind=MKD) :: npix, npface, ipf
 integer(kind=I4B) :: ip_low, ip_trunc, ip_med, ip_hi
 integer(kind=I4B) :: face_num, ix, iy, kshift, scale, i, ismax
 integer(kind=I4B) :: jrt, jr, nr, jpt, jp, nl4
 real     :: z, fn, fact1, fact2, sth, phi

 ! coordinate of the lowest corner of each face
 integer(kind=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
 integer(kind=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

 real :: phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn, sin_phi, cos_phi
 real :: z_nv, z_sv, sth_nv, sth_sv
 real :: hdelta_phi
 integer(kind=I4B) :: iphi_mod, iphi_rat
 logical(kind=LGT) :: do_vertex
 integer(kind=i4b) :: diff_phi
 character(len=*), parameter :: code = "pix2vec_nest"

 !-----------------------------------------------------------------------
 if (nside > ns_max4) call fatal_error(code//"> nside out of range")
 npix = nside2npix(nside)       ! total number of points
 if (ipix <0 .or. ipix>npix-1) call fatal_error(code//"> ipix out of range")

 !     initiates the array for the pixel number -> (x,y) mapping
 if (pix2x(1023) <= 0) call mk_pix2xy()

 npface = nside * int(nside, kind=MKD)
 nl4    = 4*nside

 !     finds the face, and the number in the face
 face_num = ipix/npface  ! face number in {0,11}
 ipf = modulo(ipix,npface)  ! pixel number in the face {0,npface-1}

 do_vertex = .false.
 if (present(vertex)) then
    if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
       do_vertex = .true.
    else
       call fatal_error(code//">  vertex array has wrong size ")
    endif
 endif
 fn = real(nside)
 fact1 = 1.0/(3.0*fn*fn)
 fact2 = 2.0/(3.0*fn)

 !     finds the x,y on the face (starting from the lowest corner)
 !     from the pixel number
 if (nside <= ns_max4) then
    ip_low = iand(ipf,1023_MKD)       ! content of the last 10 bits
    ip_trunc =    ipf/1024        ! truncation of the last 10 bits
    ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
    ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
 else
    ix = 0
    iy = 0
    scale = 1
    ismax = 4
    do i=0, ismax
       ip_low = iand(ipf,1023_MKD)
       ix = ix + scale * pix2x(ip_low)
       iy = iy + scale * pix2y(ip_low)
       scale = scale * 32
       ipf   = ipf/1024
    enddo
    ix = ix + scale * pix2x(ipf)
    iy = iy + scale * pix2y(ipf)
 endif

 !     transforms this in (horizontal, vertical) coordinates
 jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
 jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

 !     computes the z coordinate on the sphere
 jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

 z_nv = 0.; z_sv = 0.     ! avoid compiler warnings

 if (jr < nside) then     ! north pole region
    nr = jr
    z = 1. - nr*fact1*nr
    sth = nr * sqrt(fact1 * (1. + z) ) ! more accurate close to pole
    kshift = 0
    if (do_vertex) then
       z_nv = 1. - (nr-1)*fact1*(nr-1)
       z_sv = 1. - (nr+1)*fact1*(nr+1)
    endif

 elseif (jr <= 3*nside) then ! equatorial region
    nr = nside
    z  = (2*nside-jr)*fact2
    sth = sqrt((1.0-z)*(1.0+z)) ! good enough on Equator
    kshift = iand(jr - nside, 1)
    if (do_vertex) then
       z_nv = (2*nside-jr+1)*fact2
       z_sv = (2*nside-jr-1)*fact2
       if (jr == nside) then ! northern transition
          z_nv =  1.0- (nside-1) * fact1 * (nside-1)
       elseif (jr == 3*nside) then  ! southern transition
          z_sv = -1.0 + (nside-1) * fact1 * (nside-1)
       endif
    endif

 elseif (jr > 3*nside) then ! south pole region
    nr = nl4 - jr
    z   = - 1.0 + nr*fact1*nr
    sth = nr * sqrt(fact1 * (1. - z) )
    kshift = 0
    if (do_vertex) then
       z_nv = - 1.0 + (nr+1)*fact1*(nr+1)
       z_sv = - 1.0 + (nr-1)*fact1*(nr-1)
    endif
 endif

 !     computes the phi coordinate on the sphere, in [0,2Pi]
 jp = (jpll(face_num+1)*nr + jpt + 1_MKD + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
 if (jp > nl4) jp = jp - nl4
 if (jp < 1)   jp = jp + nl4

 phi = (jp - (kshift+1)*0.5) * (halfpi / nr)

 ! pixel center
 !
 cos_phi = cos(phi)
 sin_phi = sin(phi)
 vector(1) = sth * cos_phi
 vector(2) = sth * sin_phi
 vector(3) = z

 if (do_vertex) then
    phi_nv = phi
    phi_sv = phi
    diff_phi = 0 ! phi_nv = phi_sv = phisth * 1}
    iphi_rat = (jp-1) / nr      ! in {0,1,2,3}
    iphi_mod = mod(jp-1,nr)
    phi_up   = 0.
    if (nr > 1) phi_up = HALFPI * (iphi_rat +  iphi_mod   /real(nr-1))
    phi_dn             = HALFPI * (iphi_rat + (iphi_mod+1)/real(nr+1))
    if (jr < nside) then            ! North polar cap
       phi_nv = phi_up
       phi_sv = phi_dn
       diff_phi = 3 ! both phi_nv and phi_sv different from phi
    elseif (jr > 3*nside) then     ! South polar cap
       phi_nv = phi_dn
       phi_sv = phi_up
       diff_phi = 3 ! both phi_nv and phi_sv different from phi
    elseif (jr == nside) then      ! North transition
       phi_nv = phi_up
       diff_phi = 1
    elseif (jr == 3*nside) then    ! South transition
       phi_sv = phi_up
       diff_phi = 2
    endif

    hdelta_phi = PI / (4.0*nr)

    ! west vertex
    phi_wv      = phi - hdelta_phi
    vertex(1,2) = sth * cos(phi_wv)
    vertex(2,2) = sth * sin(phi_wv)
    vertex(3,2) = z

    ! east vertex
    phi_ev      = phi + hdelta_phi
    vertex(1,4) = sth * cos(phi_ev)
    vertex(2,4) = sth * sin(phi_ev)
    vertex(3,4) = z

    ! north and south vertices
    sth_nv = sqrt((1.0-z_nv)*(1.0+z_nv))
    sth_sv = sqrt((1.0-z_sv)*(1.0+z_sv))
    if (diff_phi == 0) then
       vertex(1,1) = sth_nv * cos_phi
       vertex(2,1) = sth_nv * sin_phi
       vertex(1,3) = sth_sv * cos_phi
       vertex(2,3) = sth_sv * sin_phi
    else
       vertex(1,1) = sth_nv * cos(phi_nv)
       vertex(2,1) = sth_nv * sin(phi_nv)
       vertex(1,3) = sth_sv * cos(phi_sv)
       vertex(2,3) = sth_sv * sin(phi_sv)
    endif
    vertex(3,1) = z_nv
    vertex(3,3) = z_sv
 endif

end subroutine pix2vec_nest

 !=======================================================================
 !   npix2nside
 !
 ! given npix, returns nside such that npix = 12*nside^2
 !  nside should be a power of 2 smaller than ns_max
 !  if not, -1 is returned
 ! EH, Feb-2000
 ! 2009-03-05, edited, accepts 8-byte npix
 !=======================================================================
function npix2nside  (npix) result(nside_result)
 integer(i4b), parameter :: MKD = I4B
 integer(kind=MKD), parameter  :: npix_max = (12_MKD*ns_max4)*ns_max4
 integer(kind=MKD), intent(in) :: npix
 integer(kind=MKD)             :: npix1, npix2
 integer(kind=I4B)             :: nside_result
 integer(kind=I4B)             :: nside
 character(LEN=*),  parameter  :: code = "npix2nside"
 !=======================================================================

 if (npix < 12 .or. npix > npix_max) then
    print*, code,"> Npix=",npix, &
    & " is out of allowed range: {12,",npix_max,"}"
    nside_result = -1
    return
 endif

 nside = nint( sqrt(npix/12.0_dp) )
 npix1 = (12_MKD*nside)*nside
 if (abs(npix1-npix) > 0) then
    print*, code,"> Npix=",npix, &
    & " is not 12 * Nside * Nside "
    nside_result = -1
    return
 endif

 ! test validity of Nside
 npix2 = nside2npix(nside)
 if (npix2 < 0) then
    nside_result = -1
    return
 endif

 nside_result = nside

end function npix2nside


 !=======================================================================
function nside2npix(nside) result(npix_result)
 !=======================================================================
 ! given nside, returns npix such that npix = 12*nside^2
 !  nside should be a power of 2 smaller than ns_max
 !  if not, -1 is returned
 ! EH, Feb-2000
 ! 2009-03-04: returns i8b result, faster
 !=======================================================================
 integer(kind=I4B)             :: npix_result
 integer(kind=I4B), intent(in) :: nside

 integer(kind=I4B) :: npix
 character(LEN=*), parameter :: code = "nside2npix"
 !=======================================================================

 npix = (12_i4b*nside)*nside
 if (nside < 1 .or. nside > ns_max4 .or. iand(nside-1,nside) /= 0) then
    print*,code,": Nside=",nside," is not a power of 2."
    npix = -1
 endif
 npix_result = npix

end function nside2npix

 !=======================================================================
 ! CHEAP_ISQRT
 !       Returns exact Floor(sqrt(x)) where x is a (64 bit) integer.
 !             y^2 <= x < (y+1)^2         (1)
 !       The double precision floating point operation is not accurate enough
 !        when dealing with 64 bit integers, especially in the vicinity of
 !       perfect squares.
 !=======================================================================
function cheap_isqrt(lin) result (lout)
 integer(i4b), intent(in) :: lin
 integer(i4b) :: lout
 lout = floor(sqrt(dble(lin)), kind=I4B)
 return
end function cheap_isqrt

 !=======================================================================
subroutine mk_pix2xy()
 !=======================================================================
 !     constructs the array giving x and y in the face from pixel number
 !     for the nested (quad-cube like) ordering of pixels
 !
 !     the bits corresponding to x and y are interleaved in the pixel number
 !     one breaks up the pixel number by even and odd bits
 !=======================================================================
 integer(kind=I4B) ::  kpix, jpix, ix, iy, ip, id

 !cc cf block data      data      pix2x(1023) /0/
 !-----------------------------------------------------------------------
 !      print *, 'initiate pix2xy'
 do kpix=0,1023          ! pixel number
    jpix = kpix
    IX = 0
    IY = 0
    IP = 1               ! bit position (in x and y)
    !        do while (jpix/=0) ! go through all the bits
    do
       if (jpix == 0) exit ! go through all the bits
       ID = modulo(jpix,2)  ! bit value (in kpix), goes in ix
       jpix = jpix/2
       IX = ID*IP+IX

       ID = modulo(jpix,2)  ! bit value (in kpix), goes in iy
       jpix = jpix/2
       IY = ID*IP+IY

       IP = 2*IP         ! next bit (in x and y)
    enddo
    pix2x(kpix) = IX     ! in 0,31
    pix2y(kpix) = IY     ! in 0,31
 enddo

end subroutine mk_pix2xy
 !=======================================================================
subroutine mk_xy2pix1()
 !=======================================================================
 !     sets the array giving the number of the pixel lying in (x,y)
 !     x and y are in {1,128}
 !     the pixel number is in {0,128**2-1}
 !
 !     if  i-1 = sum_p=0  b_p * 2^p
 !     then ix = sum_p=0  b_p * 4^p
 !          iy = 2*ix
 !     ix + iy in {0, 128**2 -1}
 !=======================================================================
 integer(kind=I4B):: k,ip,i,j,id
 !=======================================================================

 do i = 0,127           !for converting x,y into
    j  = i           !pixel numbers
    k  = 0
    ip = 1

    do
       if (j==0) then
          x2pix1(i) = k
          y2pix1(i) = 2*k
          exit
       else
          id = modulo(J,2)
          j  = j/2
          k  = ip*id+k
          ip = ip*4
       endif
    enddo
 enddo

end subroutine mk_xy2pix1

subroutine fatal_error (msg)
 character(len=*), intent(in), optional :: msg

 if (present(msg)) then
    print *,'Fatal error: ', trim(msg)
 else
    print *,'Fatal error'
 endif
 call exit_with_status(1)

end subroutine fatal_error

 ! ===========================================================
subroutine exit_with_status (code, msg)
 integer(i4b), intent(in) :: code
 character (len=*), intent(in), optional :: msg

 if (present(msg)) print *,trim(msg)
 print *,'program exits with exit code ', code
 call exit (code)

end subroutine exit_with_status

 !====================================================================
 ! The following is a routine which finds the 7 or 8 neighbours of
 ! any pixel in the nested scheme of the HEALPIX pixelisation.
 !====================================================================
 !  neighbours_nest
 !
 !   Returns list n(8) of neighbours of pixel ipix (in NESTED scheme)
 !   the neighbours are ordered in the following way:
 !   First pixel is the one to the south (the one west of the south
 ! direction is taken
 ! for the pixels which don't have a southern neighbour). From
 ! then on the neighbours are ordered in the clockwise direction
 ! about the pixel with number ipix.
 !
 !   nneigh is the number of neighbours (mostly 8, 8 pixels have 7 neighbours)
 !
 !   Benjamin D. Wandelt October 1997
 !   Added to pix_tools in March 1999
 !   added 'return' for case nside=1, EH, Oct 2005
 !   corrected bugs in case nside=1 and ipix=7, 9 or 11, EH, June 2006
 !   2009-06-16: deals with Nside > 8192
 !====================================================================
subroutine neighbours_nest(nside, ipix, n, nneigh)
 !   use bit_manipulation
 integer(kind=i4b), parameter  ::   MKD = I4B
 !====================================================================
 integer(kind=i4b), intent(in)::  nside
 integer(kind=MKD), intent(in)::  ipix
 integer(kind=MKD), intent(out), dimension(1:):: n
 integer(kind=i4b), intent(out):: nneigh

 integer(kind=i4b) :: ix,ixm,ixp,iy,iym,iyp,ixo,iyo
 integer(kind=i4b) :: face_num,other_face
 integer(kind=i4b) :: ia,ib,ibp,ibm,ib2,icase
 integer(kind=MKD) :: npix,ipf,ipo
 integer(kind=MKD) :: local_magic1,local_magic2,nsidesq
 character(len=*), parameter :: code = "neighbours_nest"

 !     integer(kind=i4b), intrinsic :: IAND

 !--------------------------------------------------------------------
 if (nside <1 .or. nside > ns_max4) call fatal_error(code//"> nside out of range")
 npix = nside2npix(nside) ! total number of points
 nsidesq = npix / 12
 if (ipix <0 .or. ipix>npix-1) call fatal_error(code//"> ipix out of range")

 ! quick and dirty hack for Nside=1

 if (nside == 1) then
    nneigh = 6
    if (ipix==0 ) n(1:6) = (/ 8, 4, 3, 2, 1, 5 /)
    if (ipix==1 ) n(1:6) = (/ 9, 5, 0, 3, 2, 6 /)
    if (ipix==2 ) n(1:6) = (/10, 6, 1, 0, 3, 7 /)
    if (ipix==3 ) n(1:6) = (/11, 7, 2, 1, 0, 4 /)
    if (ipix==4 ) n(1:6) = (/11, 7, 3, 0, 5, 8 /)
    if (ipix==5 ) n(1:6) = (/ 8, 4, 0, 1, 6, 9 /)
    if (ipix==6 ) n(1:6) = (/ 9, 5, 1, 2, 7,10 /)
    if (ipix==7 ) n(1:6) = (/10, 6, 2, 3, 4,11 /)
    if (ipix==8 ) n(1:6) = (/10,11, 4, 0, 5, 9 /)
    if (ipix==9 ) n(1:6) = (/11, 8, 5, 1, 6,10 /)
    if (ipix==10) n(1:6) = (/ 8, 9, 6, 2, 7,11 /)
    if (ipix==11) n(1:6) = (/ 9,10, 7, 3, 4, 8 /)
    return
 endif

 !     initiates array for (x,y)-> pixel number -> (x,y) mapping
 if (x2pix1(127) <= 0) call mk_xy2pix1()

 local_magic1=(nsidesq-1)/3
 local_magic2=2*local_magic1
 face_num=ipix/nsidesq

 ipf=modulo(ipix,nsidesq)   !Pixel number in face

 call pix2xy_nest(nside,ipf,ix,iy)
 ixm=ix-1
 ixp=ix+1
 iym=iy-1
 iyp=iy+1

 nneigh=8                  !Except in special cases below

 !     Exclude corners
 if (ipf==local_magic2)     then !WestCorner
    icase=5
    goto 100
 endif
 if (ipf==(nsidesq-1)) then !NorthCorner
    icase=6
    goto 100
 endif
 if (ipf==0)           then !SouthCorner
    icase=7
    goto 100
 endif
 if (ipf==local_magic1)     then !EastCorner
    icase=8
    goto 100
 endif

 !     Detect edges
 if (iand(ipf,local_magic1)==local_magic1) then !NorthEast
    icase=1
    goto 100
 endif
 if (iand(ipf,local_magic1)==0)      then !SouthWest
    icase=2
    goto 100
 endif
 if (iand(ipf,local_magic2)==local_magic2) then !NorthWest
    icase=3
    goto 100
 endif
 if (iand(ipf,local_magic2)==0)      then !SouthEast
    icase=4
    goto 100
 endif

 !     Inside a face
 call xy2pix_nest(nside, ixm, iym, face_num, n(1))
 call xy2pix_nest(nside, ixm, iy , face_num, n(2))
 call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
 call xy2pix_nest(nside, ix , iyp, face_num, n(4))
 call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
 call xy2pix_nest(nside, ixp, iy , face_num, n(6))
 call xy2pix_nest(nside, ixp, iym, face_num, n(7))
 call xy2pix_nest(nside, ix , iym, face_num, n(8))
 return

100 continue

 ia= face_num/4            !in {0,2}
 ib= modulo(face_num,4)       !in {0,3}
 ibp=modulo(ib+1,4)
 ibm=modulo(ib+4-1,4)
 ib2=modulo(ib+2,4)

 if (ia==0) then          !North Pole region
    select case(icase)
    case(1)              !NorthEast edge
       other_face=0+ibp
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ixm, iym, face_num, n(1))
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo+1 , iyo, other_face, n(5))
       n(6)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo-1, iyo, other_face, n(7))
    case(2)              !SouthWest edge
       other_face=4+ib
       ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
       n(2)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       call xy2pix_nest(nside, ixp, iym, face_num, n(7))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
    case(3)              !NorthWest edge
       other_face=0+ibm
       ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo, iyo-1, other_face, n(3))
       n(4)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo, iyo+1, other_face, n(5))
       call xy2pix_nest(nside, ixm, iym, face_num, n(1))
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ixp, iym, face_num, n(7))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
    case(4)              !SouthEast edge
       other_face=4+ibp
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
       n(8)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
    case(5)              !West corner
       nneigh=7
       other_face=4+ib
       n(2)=other_face*nsidesq+nsidesq-1
       n(1)=n(2)-2
       other_face=0+ibm
       n(3)=other_face*nsidesq+local_magic1
       n(4)=n(3)+2
       n(5)=ipix+1
       n(6)=ipix-1
       n(7)=ipix-2
    case(6)              !North corner
       n(1)=ipix-3
       n(2)=ipix-1
       n(8)=ipix-2
       other_face=0+ibm
       n(4)=other_face*nsidesq+nsidesq-1
       n(3)=n(4)-2
       other_face=0+ib2
       n(5)=other_face*nsidesq+nsidesq-1
       other_face=0+ibp
       n(6)=other_face*nsidesq+nsidesq-1
       n(7)=n(6)-1
    case(7)              !South corner
       other_face=8+ib
       n(1)=other_face*nsidesq+nsidesq-1
       other_face=4+ib
       n(2)=other_face*nsidesq+local_magic1
       n(3)=n(2)+2
       n(4)=ipix+2
       n(5)=ipix+3
       n(6)=ipix+1
       other_face=4+ibp
       n(8)=other_face*nsidesq+local_magic2
       n(7)=n(8)+1
    case(8)              !East corner
       nneigh=7
       n(2)=ipix-1
       n(3)=ipix+1
       n(4)=ipix+2
       other_face=0+ibp
       n(6)=other_face*nsidesq+local_magic2
       n(5)=n(6)+1
       other_face=4+ibp
       n(7)=other_face*nsidesq+nsidesq-1
       n(1)=n(7)-1
    end select ! north

 elseif (ia==1) then      !Equatorial region
    select case(icase)
    case(1)              !NorthEast edge
       other_face=0+ib
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ixm, iym, face_num, n(1))
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
       n(6)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
    case(2)              !SouthWest edge
       other_face=8+ibm
       ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
       n(2)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       call xy2pix_nest(nside, ixp, iym, face_num, n(7))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
    case(3)              !NorthWest edge
       other_face=0+ibm
       ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
       n(4)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
       call xy2pix_nest(nside, ixm, iym, face_num, n(1))
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ixp, iym, face_num, n(7))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
    case(4)              !SouthEast edge
       other_face=8+ib
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
       n(8)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
    case(5)              !West corner
       other_face=8+ibm
       n(2)=other_face*nsidesq+nsidesq-1
       n(1)=n(2)-2
       other_face=4+ibm
       n(3)=other_face*nsidesq+local_magic1
       other_face=0+ibm
       n(4)=other_face*nsidesq
       n(5)=n(4)+1
       n(6)=ipix+1
       n(7)=ipix-1
       n(8)=ipix-2
    case(6)              !North corner
       nneigh=7
       n(1)=ipix-3
       n(2)=ipix-1
       other_face=0+ibm
       n(4)=other_face*nsidesq+local_magic1
       n(3)=n(4)-1
       other_face=0+ib
       n(5)=other_face*nsidesq+local_magic2
       n(6)=n(5)-2
       n(7)=ipix-2
    case(7)              !South corner
       nneigh=7
       other_face=8+ibm
       n(1)=other_face*nsidesq+local_magic1
       n(2)=n(1)+2
       n(3)=ipix+2
       n(4)=ipix+3
       n(5)=ipix+1
       other_face=8+ib
       n(7)=other_face*nsidesq+local_magic2
       n(6)=n(7)+1
    case(8)              !East corner
       other_face=8+ib
       n(8)=other_face*nsidesq+nsidesq-1
       n(1)=n(8)-1
       n(2)=ipix-1
       n(3)=ipix+1
       n(4)=ipix+2
       other_face=0+ib
       n(6)=other_face*nsidesq
       n(5)=n(6)+2
       other_face=4+ibp
       n(7)=other_face*nsidesq+local_magic2
    end select ! equator
 else                    !South Pole region
    select case(icase)
    case(1)              !NorthEast edge
       other_face=4+ibp
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ixm, iym, face_num, n(1))
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
       n(6)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
    case(2)              !SouthWest edge
       other_face=8+ibm
       ipo=modulo(swapLSBMSB(ipf),nsidesq)        !W-E flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       n(2)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo+1, iyo, other_face, n(3))
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       call xy2pix_nest(nside, ixp, iym, face_num, n(7))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
    case(3)              !NorthWest edge
       other_face=4+ib
       ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
       n(4)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
       call xy2pix_nest(nside, ixm, iym, face_num, n(1))
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ix , iym, face_num, n(8))
       call xy2pix_nest(nside, ixp, iym, face_num, n(7))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
    case(4)              !SouthEast edge
       other_face=8+ibp
       call xy2pix_nest(nside, ixm, iy , face_num, n(2))
       call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
       call xy2pix_nest(nside, ix , iyp, face_num, n(4))
       call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       ipo=modulo(swapLSBMSB(ipf),nsidesq) !E-W flip
       call pix2xy_nest(nside,ipo,ixo,iyo)
       call xy2pix_nest(nside, ixo, iyo+1, other_face, n(7))
       n(8)=other_face*nsidesq+ipo
       call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
    case(5)              !West corner
       nneigh=7
       other_face=8+ibm
       n(2)=other_face*nsidesq+local_magic1
       n(1)=n(2)-1
       other_face=4+ib
       n(3)=other_face*nsidesq
       n(4)=n(3)+1
       n(5)=ipix+1
       n(6)=ipix-1
       n(7)=ipix-2
    case(6)              !North corner
       n(1)=ipix-3
       n(2)=ipix-1
       other_face=4+ib
       n(4)=other_face*nsidesq+local_magic1
       n(3)=n(4)-1
       other_face=0+ib
       n(5)=other_face*nsidesq
       other_face=4+ibp
       n(6)=other_face*nsidesq+local_magic2
       n(7)=n(6)-2
       n(8)=ipix-2
    case(7)              !South corner
       other_face=8+ib2
       n(1)=other_face*nsidesq
       other_face=8+ibm
       n(2)=other_face*nsidesq
       n(3)=n(2)+1
       n(4)=ipix+2
       n(5)=ipix+3
       n(6)=ipix+1
       other_face=8+ibp
       n(8)=other_face*nsidesq
       n(7)=n(8)+2
    case(8)              !East corner
       nneigh=7
       other_face=8+ibp
       n(7)=other_face*nsidesq+local_magic2
       n(1)=n(7)-2
       n(2)=ipix-1
       n(3)=ipix+1
       n(4)=ipix+2
       other_face=4+ibp
       n(6)=other_face*nsidesq
       n(5)=n(6)+2
    end select ! south
 endif

end subroutine neighbours_nest


 !=======================================================================
 !  pix2xy_nest
 !     gives the x, y coords in a face from pixel number within the face (NESTED)
 !
 !     Benjamin D. Wandelt 13/10/97
 !
 !     using code from HEALPIX toolkit by K.Gorski and E. Hivon
 !     2009-06-15: deals with Nside > 8192
 !     2012-03-02: test validity of ipf_in instead of undefined ipf
 !                 define ipf as MKD
 !     2012-08-27:  corrected bug on (ix,iy) for Nside > 8192 (MARK)
 !=======================================================================
subroutine pix2xy_nest  (nside, ipf_in, ix, iy)
 integer(kind=i4b), parameter  ::   MKD = I4B
 integer(kind=I4B), intent(in)  :: nside
 integer(kind=MKD), intent(in)  :: ipf_in
 integer(kind=I4B), intent(out) :: ix, iy

 integer(kind=MKD) :: ipf
 integer(kind=I4B) ::  ip_low, ip_trunc, ip_med, ip_hi, scale, i, ismax
 character(len=*), parameter :: code = "pix2xy_nest"

 !-----------------------------------------------------------------------
 if (nside<1 .or. nside>ns_max) call fatal_error(code//"> nside out of range")
 if (ipf_in<0 .or. ipf_in>nside*nside-1) &
 &     call fatal_error(code//"> ipix out of range")
 if (pix2x(1023) <= 0) call mk_pix2xy()

 ipf = ipf_in
 if (nside <= ns_max4) then
    ip_low = iand(ipf,1023_MKD)   ! content of the last 10 bits
    ip_trunc =    ipf/1024        ! truncation of the last 10 bits
    ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
    ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
 else
    ix = 0
    iy = 0
    scale = 1
    ismax = 4
    do i=0, ismax
       ip_low = iand(ipf,1023_MKD)
       ix = ix + scale * pix2x(ip_low)
       iy = iy + scale * pix2y(ip_low) ! corrected 2012-08-27
       scale = scale * 32
       ipf   = ipf/1024
    enddo
    ix = ix + scale * pix2x(ipf)
    iy = iy + scale * pix2y(ipf) ! corrected 2012-08-27
 endif

end subroutine pix2xy_nest

 !=======================================================================
 !     gives the pixel number ipix (NESTED)
 !     corresponding to ix, iy and face_num
 !
 !     Benjamin D. Wandelt 13/10/97
 !     using code from HEALPIX toolkit by K.Gorski and E. Hivon
 !     2009-06-15: deals with Nside > 8192
 !     2012-03-02: test validity of ix_in and iy_in instead of undefined ix and iy
 !=======================================================================
subroutine xy2pix_nest(nside, ix_in, iy_in, face_num, ipix)
 integer(kind=i4b), parameter  ::   MKD = I4B
 !=======================================================================
 integer(kind=I4B), intent(in) ::  nside, ix_in, iy_in, face_num
 integer(kind=MKD), intent(out) :: ipix
 integer(kind=I4B) ::  ix, iy, ix_low, iy_low, i, ismax
 integer(kind=MKD) :: ipf, scale, scale_factor
 character(len=*), parameter :: code = "xy2pix_nest"

 !-----------------------------------------------------------------------
 if (nside<1 .or. nside>ns_max) call fatal_error(code//"> nside out of range")
 if (ix_in<0 .or. ix_in>(nside-1)) call fatal_error(code//"> ix out of range")
 if (iy_in<0 .or. iy_in>(nside-1)) call fatal_error(code//"> iy out of range")
 if (x2pix1(127) <= 0) call mk_xy2pix1()

 ix = ix_in
 iy = iy_in
 if (nside <= ns_max4) then
    ix_low = iand(ix, 127)
    iy_low = iand(iy, 127)
    ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
    & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
 else
    scale = 1_MKD
    scale_factor = 16384_MKD ! 128*128
    ipf = 0_MKD
    ismax = 1 ! for nside in [2^14, 2^20]
    if (nside >  1048576 ) ismax = 3
    do i=0, ismax
       ix_low = iand(ix, 127) ! last 7 bits
       iy_low = iand(iy, 127) ! last 7 bits
       ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
       scale = scale * scale_factor
       ix  =     ix / 128 ! truncate out last 7 bits
       iy  =     iy / 128
    enddo
    ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
 endif
 ipix = ipf + face_num* int(nside,MKD) * nside    ! in {0, 12*nside**2 - 1}

end subroutine xy2pix_nest

end module healpix
