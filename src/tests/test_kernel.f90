!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testkernel
!
!  DESCRIPTION:
!  This module performs unit tests of the kernel routines
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, kernel, testutils
!+
!--------------------------------------------------------------------------
module testkernel
 implicit none
 public :: test_kernel

 private

contains

subroutine test_kernel(ntests,npass)
 use io,        only:id,master
 use kernel,    only:kernelname,get_kernel,wkern,grkern,wab0,gradh0,radkern,radkern2,cnormk, &
                     kernel_softening,get_kernel_grav1,dphidh0
 use testutils, only:checkvalbuf,checkvalbuf_end,checkval
 integer, intent(inout) :: ntests,npass
 integer, parameter :: n = 200
 integer, parameter :: stderr = 0
 integer, parameter :: nktest = 7
 integer            :: nerr(nktest),ncheck(nktest),i
 real :: wval,grkernval,gradhval,wval2,grkernval2,wval3,grkernval3,we,dw,dphidh
 real :: dq,q2,qi,errmax(nktest),eps,tolgrad
 real :: potensoft,fsoft,potensofte,fsofte,dp

 if (id==master) write(stderr,"(a,/)") '--> TESTING '//trim(kernelname)//' kernel'
!
!--check that wab0 and gradh0 are equal to values at q=zero
!
 ntests = ntests + 1
 call get_kernel(0.,0.,wval,grkernval)
 call checkval(wab0,wval,tiny(0.),nerr(1),'wab0 = wab(0)')
 if (nerr(1)==0) npass = npass + 1

 ntests = ntests + 1
 gradhval = -3.*wval
 call checkval(gradh0,gradhval,tiny(0.),nerr(1),'gradh0 = -3.*wab(0)')
 if (nerr(1)==0) npass = npass + 1

 ! test get_kernel_grav1 routine
 ntests = ntests + 1
 call get_kernel_grav1(0.,0.,wval,grkernval,dphidh)
 call checkval(dphidh0,dphidh,tiny(0.),nerr(1),'dphidh0 = dphidh(0)')
 if (nerr(1)==0) npass = npass + 1
!
!--check that radkern2 = radkern**2
!
 ntests = ntests + 1
 call checkval(radkern2,radkern**2,tiny(0.),nerr(1),'radkern2 = radkern*radkern')
 if (nerr(1)==0) npass = npass + 1
!
!--check that all three functions give consistent answers
!  for q=0..radkern2
!
 ntests = ntests + nktest
 errmax = 0.
 dq = radkern/real(n-1)
 !open(unit=1,file='kernelfunc-'//trim(kernelname)//'.out',status='replace')
 ncheck = 0
 nerr = 0
 eps = 1.e-7
 tolgrad = 2.e-5
 do i=1,n
    qi = (i-1)*dq
    q2 = qi*qi
    call get_kernel(q2,qi,wval,grkernval)
    wval2      = wkern(q2,qi)
    grkernval2 = grkern(q2,qi)
    !write(1,*) qi,cnormk*wval,cnormk*grkernval
    call checkvalbuf(wval,wval2,tiny(0.),'W from get_kernel subroutine /= W from wkern function',nerr(1),ncheck(1),errmax(1))
    call checkvalbuf(grkernval,grkernval2,tiny(0.), &
                     'gradW from get_kernel subroutine /= gradW from grkern function',nerr(2),ncheck(2),errmax(2))
    ! check that kernel gradient is OK
    call get_kernel((qi+eps)**2,qi + eps,we,grkernval2)
    dw = (we - wval)/eps
    call checkvalbuf(dw,grkernval,tolgrad, &
                     'gradient of kernel incorrect ',nerr(3),ncheck(3),errmax(3))
    !-----------------
    ! softening tests
    !-----------------
    ! check that fsoft is gradient of potensoft
    call kernel_softening(q2,qi,potensoft,fsoft)
    call kernel_softening((qi+eps)**2,qi+eps,potensofte,fsofte)
    dp = (potensofte - potensoft)/eps
    call checkvalbuf(dp,fsoft,tolgrad,'gradient of potential /= force in kernel',nerr(4),ncheck(4),errmax(4))

    ! test get_kernel_grav1 routine
    call get_kernel_grav1(q2,qi,wval3,grkernval3,dphidh)
    call checkvalbuf(wval,wval3,tiny(0.),&
         'W from get_kernel_grav1 /= W from wkern function',nerr(5),ncheck(5),errmax(5))
    call checkvalbuf(grkernval,grkernval3,tiny(0.),&
         'gradW from get_kernel_grav1 /= gradW from grkern function',nerr(6),ncheck(6),errmax(6))

    ! check that dphidh is gradient of potential w.r.t. h
    call checkvalbuf(dphidh,-potensoft -qi*fsoft,2.e-7,'dphidh /= phi - q*dphi/dq',nerr(7),ncheck(7),errmax(7))
 enddo
 !close(unit=1)
 do i=1,nktest
    if (nerr(i)==0) npass = npass + 1
 enddo

 call checkvalbuf_end('get_kernel == wkern',ncheck(1),nerr(1),errmax(1),tiny(0.))
 call checkvalbuf_end('get_kernel == grkern',ncheck(2),nerr(2),errmax(2),tiny(0.))
 call checkvalbuf_end('w gradient equal to gradw',ncheck(3),nerr(3),errmax(3),tolgrad)
 call checkvalbuf_end('potential gradient = force',ncheck(4),nerr(4),errmax(4),tolgrad)
 call checkvalbuf_end('get_kernel_grav1 == wkern',ncheck(5),nerr(5),errmax(5),tiny(0.))
 call checkvalbuf_end('get_kernel_grav1 == grkern',ncheck(6),nerr(6),errmax(6),tiny(0.))
 call checkvalbuf_end('dphi/dh = phi - q*dphi/dq',ncheck(7),nerr(7),errmax(7),2.e-7)

 if (id==master) write(stderr,"(/,a,/)") '<-- KERNEL TEST COMPLETE'

end subroutine test_kernel

end module testkernel
