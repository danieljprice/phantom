!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testindtstep
!
! test module for individual timestepping utilities
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, testutils, timestep_ind
!
 implicit none
 public :: test_indtstep

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests for individual timestep utilities
!+
!-----------------------------------------------------------------------
subroutine test_indtstep(ntests,npass)
 use dim,             only:ind_timesteps
 use io,              only:id,master,iverbose
 use testutils,       only:checkvalbuf,checkvalbuf_end,checkval,update_test_scores
 use timestep_ind,    only:change_nbinmax,get_newbin,maxbins,istepfrac,nbinmax
 integer, intent(inout) :: ntests,npass
 integer :: istep
 integer(kind=1) :: i,nbinmaxprev,ibini
 real :: dtmax,dt,fracin,fracout,errmax1,errmax2
 integer :: ncheck(2)
 integer :: nfailed(32)
 character(len=10) :: string

 if (ind_timesteps) then
    if (id==master) write(*,"(/,a,/)") '--> TESTING IND TIMESTEP UTILS'
 else
    if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF IND TIMESTEP UTILS (need -DIND_TIMESTEPS)'
    return
 endif

 dtmax   = 1.
 !
 !--test that the get_newbin routine returns the correct timestep bin
 !
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
 call update_test_scores(ntests,nfailed,npass)

 ncheck  = 0
 nfailed = 0
 errmax1 = 0.
 errmax2 = 0.
 iverbose = 0

 if (id==master) print "(/,a)",' ---> checking change_nbinmax routine'
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
 call update_test_scores(ntests,nfailed(1:1),npass)
 call update_test_scores(ntests,nfailed(2:2),npass)

 ! reset nbinmax to avoid problems with subsequent tests
 nbinmax = 0

 if (id==master) write(*,"(/,a)") '<-- IND TIMESTEP UTILS TEST COMPLETE'

end subroutine test_indtstep

end module testindtstep
