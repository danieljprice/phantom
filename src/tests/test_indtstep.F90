!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testindtstep
!
!  DESCRIPTION:
!  test module for individual timestepping utilities
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, testutils, timestep_ind
!+
!--------------------------------------------------------------------------
module testindtstep
 implicit none
 public :: test_indtstep

 private

contains

subroutine test_indtstep(ntests,npass)
 use io,              only:id,master
#ifdef IND_TIMESTEPS
 use testutils,       only:checkvalbuf,checkvalbuf_end,checkval
 use timestep_ind,    only:change_nbinmax,get_newbin,maxbins,istepfrac,nbinmax
 use io,              only:iverbose
#endif
 integer, intent(inout) :: ntests,npass
#ifdef IND_TIMESTEPS
 integer :: istep
 integer(kind=1) :: i,nbinmaxprev,ibini
 real :: dtmax,dt,fracin,fracout,errmax1,errmax2
 integer :: ncheck(2)
 integer :: nfailed(32)
 character(len=10) :: string

 if (id==master) write(*,"(/,a,/)") '--> TESTING IND TIMESTEP UTILS'

 dtmax   = 1.
 !
 !--test that the get_newbin routine returns the correct timestep bin
 !
 ntests = ntests + 1
 nfailed(:) = 0
 dtmax = 1.
 istepfrac = 0
 nbinmax = maxbins
 dt = 3.*dtmax
 do i=1,30
    ibini = 1
    dt = dt/2.
    write(string,"(es10.3)") dt
    call get_newbin(dt,dtmax,ibini)
    call checkval(int(ibini),i-1,0,nfailed(i),'bin for dt = '//trim(string))
 enddo
 !--check dt = huge
 ibini = 0
 call get_newbin(huge(dt),dtmax,ibini)
 call checkval(int(ibini),0,0,nfailed(31),'bin for dt = huge')
 !--check dt = 0.
 call get_newbin(0.,dtmax,ibini)
 call checkval(int(ibini),maxbins,0,nfailed(32),'bin for dt = 0')
 if (all(nfailed==0)) npass = npass + 1

 ncheck  = 0
 nfailed = 0
 errmax1 = 0.
 errmax2 = 0.
 iverbose = 0

 if (id==master) print "(/,a)",' ---> checking change_nbinmax routine'
 ntests = ntests + 2
 !
 !--loop over all possible timestep bins
 !
 do i=1,24  ! 29 is maximum possible, but gets slow
    nbinmax = i
    nbinmaxprev = i
    if (id==master) print "(a,i2)",' checking nbinmax = ',nbinmax
    !
    !--loop over all possible fractional timesteps for this nbins
    !
    do istep = 1,2**i

       nbinmaxprev = i
       istepfrac   = istep
       fracin      = istepfrac/real(2**nbinmaxprev)
       !
       !--check that increasing nbins works
       !
       nbinmax     = nbinmaxprev + 1_1
       call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
       fracout = istepfrac/real(2**nbinmax)
       !print *,istep,'/',2**nbinmaxprev,' fraction before/after = ',fracin,fracout
       call checkvalbuf(fracout,fracin,tiny(0.),'increasing nbins keeps step fraction',nfailed(1),ncheck(1),errmax1)
       !
       !--if stepfrac is a multiple of 2, check that decreasing nbins works
       !
       istepfrac   = istep
       if (mod(istepfrac,2)==0) then
          nbinmaxprev = i
          fracin      = istepfrac/real(2**nbinmaxprev)
          !
          !--check that decreasing nbins works
          !
          nbinmax     = nbinmaxprev - 1_1
          call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
          fracout = istepfrac/real(2**nbinmax)
          !print *,istep,'/',2**nbinmaxprev,' fraction before/after = ',fracin,fracout
          call checkvalbuf(fracout,fracin,tiny(0.),'decreasing nbins keeps step fraction',nfailed(2),ncheck(2),errmax2)

       endif
    enddo
 enddo

 call checkvalbuf_end('increasing nbins keeps step fraction',ncheck(1),nfailed(1),errmax1,tiny(0.))
 call checkvalbuf_end('decreasing nbins keeps step fraction',ncheck(2),nfailed(2),errmax2,tiny(0.))
 if (nfailed(1)==0) npass = npass + 1
 if (nfailed(2)==0) npass = npass + 1

 ! reset nbinmax to avoid problems with subsequent tests
 nbinmax = 0

 if (id==master) write(*,"(/,a)") '<-- IND TIMESTEP UTILS TEST COMPLETE'

#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF IND TIMESTEP UTILS (need -DIND_TIMESTEPS)'

#endif

end subroutine test_indtstep

end module testindtstep
