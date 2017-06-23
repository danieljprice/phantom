!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: random
!
!  DESCRIPTION:
!  this module contains a motley collection of random number
!  generator routines, used in various particle setups
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
module random
 implicit none
 public :: ran1,ran2,rayleigh_deviate
 public :: sobseq
 public :: ranset,dran

 private

contains

!!------------------------------------------------------------------------!!
!!
!! Random number generator using the minimal standard generator of
!!  Park & Miller (1988) + shuffling (see Press et al, Numerical Recipes)
!!
!! Period is about 10**8
!!
!! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
!!  endpoints). Call with iseed < 0 to initialise, thereafter do not
!!  alter iseed between calls.
!!
!!------------------------------------------------------------------------!!

real function ran1(iseed)
 integer, intent(inout) :: iseed
 integer, parameter :: ia = 16807, im=2147483647, iq = 127773, ir = 2836
 integer, parameter :: ntab = 32, ndiv = 1+(im-1)/ntab
 integer :: iv(ntab)
 integer :: j,k,iy
 real, parameter :: am = 1./im, eps = 1.2e-7, floatmax = 1.-eps

 save iv,iy
 data iv /ntab*0/, iy /0/
!
!--initialise
!
 if (iseed <= 0 .or. iy==0) then
    iseed = max(-iseed,1)  ! do not allow iseed = 0
    do j = ntab+8,1,-1
       k = iseed/iq
       iseed = ia*(iseed-k*iq) - ir*k
       if (iseed < 0) iseed = iseed + im
       if (j <= ntab) iv(j) = iseed
    enddo
    iy = iv(1)
 endif
!
!--generate random number
!
 k = iseed/iq
 iseed = ia*(iseed-k*iq) - ir*k
 if (iseed < 0) iseed = iseed + im
 j = 1 + iy/ndiv
 iy = iv(j)
 iv(j) = iseed
 ran1 = min(am*iy,floatmax)

 return
end function ran1

!!------------------------------------------------------------------------!!
!!
!! Long period random number generator (see Press et al, Numerical Recipes)
!!
!! Period is about 2 x 10**18
!!
!! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
!!  endpoints). Call with iseed < 0 to initialise, thereafter do not
!!  alter iseed between calls.
!!
!!------------------------------------------------------------------------!!

real function ran2(iseed)
 integer, parameter :: im1=2147483563, im2=2147483399, &
   imm1=im1-1, ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211, &
   ir2=3791, ntab=32,ndiv=1+imm1/ntab
 real, parameter :: am=1./im1, eps=1.2e-7, rnmx=1.-eps
 integer :: iseed
 integer :: iseed2,j,k,iv(ntab),iy
 save iv,iy,iseed2

 data iseed2/123456789/, iv/ntab*0/, iy/0/
!
!--initialise random sequence
!
 if (iseed <= 0) then
    iseed = max(-iseed,1) ! iseed not zero
    iseed2 = iseed
    do j=ntab+8,1,-1
       k = iseed/iq1
       iseed = ia1*(iseed-k*iq1) - k*ir1
       if (iseed < 0) iseed = iseed + im1
       if (j <= ntab) iv(j) = iseed
    enddo
    iy = iv(1)
 endif
 k = iseed/iq1
 iseed = ia1*(iseed-k*iq1) - k*ir1
 if (iseed < 0) iseed = iseed + im1
 k = iseed2/iq2
 iseed2 = ia2*(iseed2-k*iq2) - k*iq2
 if (iseed2 < 0) iseed2 = iseed2 + im2
 j = 1 + iy/ndiv
 iy = iv(j) - iseed2
 iv(j) = iseed
 if (iy < 1) iy = iy + imm1
 ran2 = min(am*iy,rnmx)

 return
end function ran2

!!-------------------------------------------------------------------------
!!
!! Function returns a random number drawn from a Rayleigh distribution
!! P(r) = r*e^(-r^2/(2*s^2))/s^2
!!
!! Useful for drawing amplitudes from a Gaussian distribution,
!! since the modulus is distributed according to a Rayleigh distribution.
!!
!!-------------------------------------------------------------------------
real function rayleigh_deviate(iseed)
 integer :: iseed

 rayleigh_deviate = sqrt(-log(ran1(iseed)))

end function rayleigh_deviate

!!-------------------------------------------------------------------------
!!
!! Quasi Random sobol sequence from Numerical Recipes
!!
!!-------------------------------------------------------------------------
subroutine sobseq(n,x)
 integer, intent(in) :: n
 real    :: x(*)
 integer, parameter :: maxbit=30, maxdim=6
 integer i,im,in,ipp,j,k,l,ip(maxdim),iu(maxdim,maxbit)
 integer iv(maxbit*maxdim),ix(maxdim),mdeg(maxdim)
 real fac
 save ip,mdeg,ix,iv,in,fac
 equivalence (iv,iu)
 data ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
 data iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
 if (n < 0) then
    do k=1,maxdim
       do j=1,mdeg(k)
          iu(k,j)=iu(k,j)*2**(maxbit-j)
       enddo
       do j=mdeg(k)+1,maxbit
          ipp=ip(k)
          i=iu(k,j-mdeg(k))
          i=ieor(i,i/2**mdeg(k))
          do l=mdeg(k)-1,1,-1
             if(iand(ipp,1) /= 0)i=ieor(i,iu(k,j-l))
             ipp=ipp/2
          enddo
          iu(k,j)=i
       enddo
    enddo
    fac=1./2.**maxbit
    in=0
 else
    im=in
    do j=1,maxbit
       if(iand(im,1)==0)goto 1
       im=im/2
    enddo
    stop 'maxbit too small in sobseq'
1   im=(j-1)*maxdim
    do k=1,min(n,maxdim)
       ix(k)=ieor(ix(k),iv(im+k))
       x(k)=ix(k)*fac
    enddo
    in=in+1
 endif
 return
end subroutine sobseq

!--------------------------------------------------
! Black box random number generator routines follow
!--------------------------------------------------

!--------------------------------
subroutine ranset(ir,iseed)
!
!  This routine is used to initialize all of the
!  random number generators in RANPAK.
!
!  The seeding is performed as determined by ir
!  and iseed as follows:
!
!     ir <= 0, initialize all seeds
!     ir > 0, initialize only seed ir
!     iseed < 1, the default seed is used
!     iseed >= 1, use the value passed in iseed
!
!  After initialization, the value of the jth
!  random number seed should be accessed with
!  the RANGET routine.
!
!  CODE DEPENDENCIES: (internal RANPAK COMMON block)
!
!  DATE: DEC. 1, 1994
!  AUTHOR: R.D. STEWART
!
 integer ir,iseed,iset,i
 integer nrg
 parameter (nrg=100)
 integer isd(nrg)
 common /ranpak/ isd

 if (iseed > 0) then
    iset = iseed
 else
    ! use the "demonic" default seed
    iset = 666
 endif

 if (ir < 1) then
    !  initialize all random number generator seeds
    !  to iseed
    do i=1,nrg
       isd(i) = iset
    enddo
 elseif (ir <= nrg) then
    isd(ir) = iset
 endif

 ir = 1

 return
end subroutine ranset

subroutine dran(ir,rand_no)
!
!  PREFERRED DOUBLE PRECISION RANDOM NUMBER GENERATOR
!
!  The random number seed is set initially by a call to the
!  RANSET routine.  ir is a pointer to an internal random
!  number seed and NOT a random number seed.
!
!  COMMENTS: This routine will work correctly on any computer
!            with a maximum integer greater than or equal
!            to 2**31 - 1.  NOTE: The proper function of
!            the routine can verified by checking to see
!            that the random number seed after the generation
!            of 10000 random numbers is 1,043,618,065
!
!            For more information refer to the classic paper
!
!            Park, S.K. and Miller, K.W., "Random number generators:
!            good ones are hard to find."  Communications of the
!            ACM, Vol 31, No 10 (Oct. 1988).
!
!  CODE DEPENDENCIES: (internal RANPAK COMMON block)
!
!  DATE: DEC. 1, 1994
!  AUTHOR: R.D. STEWART
!
 integer nrg,ir
 parameter (nrg=100)
 integer isd(nrg)
 common /ranpak/ isd

 double precision tmp,a,seed,rand_no
 double precision m,minv
 parameter (a=16807.0d+00,m=2147483647.0d+00)
 parameter (minv=1.0d+00/m)

 seed = isd(ir)

 tmp = a*seed
 seed = tmp - m*dint(tmp*minv)
 rand_no = seed*minv

 isd(ir) = int(seed)

 return
end subroutine dran

end module random
