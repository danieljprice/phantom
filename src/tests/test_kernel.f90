!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testkernel
!
! This module performs unit tests of the kernel routines
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, kernel, testutils
!
 implicit none
 public :: test_kernel

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of kernel functions
!+
!-----------------------------------------------------------------------
subroutine test_kernel(ntests,npass)
 use io,        only:id,master
 use kernel,    only:kernelname,get_kernel,wkern,grkern,wab0,gradh0,radkern,radkern2,cnormk, &
                     kernel_softening,get_kernel_grav1,dphidh0, &
                     get_kernel_tilde,wab0_tilde,gradh0_tilde
 use testutils, only:checkvalbuf,checkvalbuf_end,checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer, parameter :: n = 200
 integer, parameter :: stderr = 0
 integer, parameter :: nktest = 8
 integer            :: nerr(nktest),ncheck(nktest),i
 real :: wval,grkernval,gradhval,wval2,grkernval2,wval3,grkernval3,we,dw,dphidh
 real :: dq,q2,qi,errmax(nktest),eps,tolgrad
 real :: potensoft,fsoft,potensofte,fsofte,dp
 real :: wtilde,grtilde,wetilde,dwt

 if (id==master) write(*,"(a,/)") '--> TESTING '//trim(kernelname)//' kernel'
!
!--check that wab0 and gradh0 are equal to values at q=zero
!
 call get_kernel(0.,0.,wval,grkernval)
 call checkval(wab0,wval,tiny(0.),nerr(1),'wab0 = wab(0)')
 call update_test_scores(ntests,nerr(1:1),npass)

 gradhval = -3.*wval
 call checkval(gradh0,gradhval,tiny(0.),nerr(1),'gradh0 = -3.*wab(0)')
 call update_test_scores(ntests,nerr(1:1),npass)

 ! test get_kernel_grav1 routine
 call get_kernel_grav1(0.,0.,wval,grkernval,dphidh)
 call checkval(dphidh0,dphidh,tiny(0.),nerr(1),'dphidh0 = dphidh(0)')
 call update_test_scores(ntests,nerr(1:1),npass)
!
!--check that radkern2 = radkern**2
!
 call checkval(radkern2,radkern**2,tiny(0.),nerr(1),'radkern2 = radkern*radkern')
 call update_test_scores(ntests,nerr(1:1),npass)
!
!--number-density (smooth-shift-4) companion: self-value at q=0
!
 call get_kernel_tilde(0.,0.,wtilde,grtilde)
 call checkval(wab0_tilde,wtilde,tiny(0.),nerr(1),'wab0_tilde = Wtilde(0)')
 call update_test_scores(ntests,nerr(1:1),npass)

 gradhval = -3.*wtilde
 call checkval(gradh0_tilde,gradhval,tiny(0.),nerr(1),'gradh0_tilde = -3.*Wtilde(0)')
 call update_test_scores(ntests,nerr(1:1),npass)

 ! flat core and compact support: Wtilde' = 0 at q=0 and at radkern
 call checkval(grtilde,0.,tiny(0.),nerr(1),'dWtilde/dq(0) = 0')
 call update_test_scores(ntests,nerr(1:1),npass)
 call get_kernel_tilde(radkern2,radkern,wtilde,grtilde)
 call checkval(grtilde,0.,2.e-6,nerr(1),'dWtilde/dq(R) = 0')
 call update_test_scores(ntests,nerr(1:1),npass)
!
!--check that all three functions give consistent answers
!  for q=0..radkern2
!
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

    ! number-density companion: finite-difference gradient of Wtilde
    ! Expanded Wtilde is ill-conditioned; compare only where |W'| is not tiny
    call get_kernel_tilde(q2,qi,wtilde,grtilde)
    if (abs(grtilde) > 1.e-2) then
       call get_kernel_tilde((qi+eps)**2,qi+eps,wetilde,grkernval2)
       call get_kernel_tilde(max(qi-eps,0.)**2,max(qi-eps,0.),wval3,grkernval3)
       dwt = 0.5*(wetilde - wval3)/eps
       call checkvalbuf(dwt,grtilde,5.e-3, &
                        'gradient of Wtilde incorrect ',nerr(8),ncheck(8),errmax(8))
    endif
 enddo
 !close(unit=1)
 do i=1,nktest
    call update_test_scores(ntests,nerr(i:i),npass)
 enddo

 call checkvalbuf_end('get_kernel == wkern',ncheck(1),nerr(1),errmax(1),tiny(0.))
 call checkvalbuf_end('get_kernel == grkern',ncheck(2),nerr(2),errmax(2),tiny(0.))
 call checkvalbuf_end('w gradient equal to gradw',ncheck(3),nerr(3),errmax(3),tolgrad)
 call checkvalbuf_end('potential gradient = force',ncheck(4),nerr(4),errmax(4),tolgrad)
 call checkvalbuf_end('get_kernel_grav1 == wkern',ncheck(5),nerr(5),errmax(5),tiny(0.))
 call checkvalbuf_end('get_kernel_grav1 == grkern',ncheck(6),nerr(6),errmax(6),tiny(0.))
 call checkvalbuf_end('dphi/dh = phi - q*dphi/dq',ncheck(7),nerr(7),errmax(7),2.e-7)
 call checkvalbuf_end('Wtilde gradient equal to dWtilde/dq',ncheck(8),nerr(8),errmax(8),5.e-3)

 if (id==master) write(*,"(/,a,/)") '<-- KERNEL TEST COMPLETE'

end subroutine test_kernel

end module testkernel
