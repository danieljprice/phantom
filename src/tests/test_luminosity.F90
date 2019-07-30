!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testlum
!
!  DESCRIPTION:
!   Tests lightcurve and timestepping
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: deriv, dim, energies, eos, io, options, part, setdisc,
!    testutils, timestep_ind, timing, viscosity
!+
!--------------------------------------------------------------------------
module testlum
 implicit none
 public :: test_lum

 private

contains

subroutine test_lum(ntests,npass)
 use dim,      only:periodic,lightcurve
 use io,       only:id,master
#ifdef LIGHTCURVE
 use io,       only:iverbose
 use part,     only:npart,npartoftype,massoftype,xyzh,hfact,vxyzu,fxyzu,fext,&
                    igas,divcurlv,iphase,isetphase,maxphase,mhd,dustprop,ddustprop,&
                    Bevol,dBevol,dustfrac,ddustevol,temperature,divcurlB
 use eos,             only:gamma,polyk
 use testutils,       only:checkval,checkvalf
 use energies,        only:compute_energies,ekin,etherm,totlum !etot,eacc,accretedmass
 use setdisc,         only:set_disc
 use deriv,           only:derivs
 use timing,          only:getused
#ifndef DISC_VISCOSITY
 use dim,             only:maxp
 use part,            only:alphaind,maxalpha
 use options,         only:alphau,alphaB
 use viscosity,       only:irealvisc,shearfunc,dt_viscosity,shearparam
#endif
 use options,         only:ieos,alpha,iexternalforce,ipdv_heating,ishock_heating
 !use part,            only:luminosity
#ifdef IND_TIMESTEPS
 use timestep_ind,    only:nactive
#endif
#endif
 integer, intent(inout) :: ntests,npass
#ifdef LIGHTCURVE
 integer                :: i,itest
 real                   :: totlum_saved(2),dtext_dum,etot_saved(2),diff,alpha_in
 real                   :: time
 real(kind=4) :: t1,t2
 integer                :: nfail(1),ii
#endif

!#ifdef DISC_VISCOSITY
!    if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF LIGHTCURVE (cannot have -DDISC_VISCOSITY)'
!    return
!#endif

#ifdef LIGHTCURVE
 if (periodic) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING LUMINOSITY TESTS'
    return
 else
#ifdef IND_TIMESTEPS
    if (id==master) write(*,"(/,a,/)") '--> TESTING LUMINOSITY'
#else
    if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF LIGHTCURVE (need -DIND_TIMESTEPS too)'
    return
#endif
 endif

 npart = min(size(xyzh(1,:)),100000)
 npartoftype(:) = 0
 npartoftype(1) = npart
 gamma = 1.0
 time = 0.0
 hfact = 1.2
 totlum_saved = 7.0
 iexternalforce = 1
 iverbose = 0
 ieos = 2

 alpha_in = 0.1
 alpha = alpha_in

! Run with "alpha_SS" viscosity first
#ifndef DISC_VISCOSITY
 irealvisc = 2
 shearparam = alpha_in
 alpha = 0.
 alphau = 0.
 alphaB = 0.
 if (maxalpha==maxp)  alphaind(:,:) = real(alpha_in,kind=kind(alphaind))
#endif

 call set_disc(id,master=master,&
                   nparttot= npart, &
                   npart   = npart,&
                   rmin    = 0.5, &
                   rmax    = 10.,&
                   p_index = 1.5,    &
                   q_index = 0.75,   &
                   HoverR  = 0.02, &
                   disc_mass = 1.e-4,   &
                   star_mass = 1.0,    &
                   gamma   = gamma,  &
                   particle_mass = massoftype(1), &
                   hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,polyk=polyk,&
                   alpha=alpha,ismooth=.true.,writefile=.false.)

 if (mhd) Bevol(:,:) = 0.

 do ii=1,4
    if (ii == 1) then
       ipdv_heating = 0
       ishock_heating = 1
    elseif (ii == 2) then
       ipdv_heating = 1
       ishock_heating = 0
    elseif (ii==3) then
       ipdv_heating = 1
       ishock_heating = 1
    else
       ipdv_heating = 0
       ishock_heating = 0
    endif

#ifdef IND_TIMESTEPS
    do itest = 1,2
       if (itest == 1) then
          iphase(1:npart) = isetphase(igas,iactive=.true.)
          nactive = npart
       else
          nactive = 10**(itest+1)
          do i=1,npart
             if (i <= nactive) then
                iphase(i) = isetphase(igas,iactive=.true.)
             else
                iphase(i) = isetphase(igas,iactive=.false.)
             endif
          enddo
       endif

       !
       !--calculate derivatives
       !
       !print*,nactive,' particles active'
       call getused(t1)
       fext = 0.
       call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
       call getused(t2)

       !print*,maxalpha,maxp,alphaind(1)

       totlum_saved(itest) = totlum
       call compute_energies(time)
       etot_saved(itest) = ekin + etherm
       !if (itest == 2) print*,totlum_saved
    enddo

! Checking the total energy sum - currently does not pass
    ! do itest=1,2
    !    ntests = ntests + 1
    !    call checkval(etot_saved(itest) + totlum_saved(itest),0.0,1.e-12,nfail,'totlum compared to etot')
    !    if (nfail == 0) npass = npass + 1
    ! enddo

    diff = (totlum_saved(1) - totlum_saved(2))/totlum_saved(1)
    call checkval(totlum_saved(1),totlum_saved(2),0.01,nfail(1),'totlum')
    call update_test_scores(ntests,nfail(1:1),npass)
#endif
    nactive = npart !for later
 enddo

! print*,'COMPARE THIS',totlum_saved,etot,etherm,ekin

!-- Check with regular viscosity
 call getused(t1)
 call derivs(1,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
 call getused(t2)
 totlum_saved(2) = totlum
 diff = (totlum_saved(1) - totlum_saved(2))/totlum_saved(1)
! call checkval(diff,0.0,0.01,nfail,'physical viscosity on/off')
! if (nfail==0) npass = npass + 1

! print*,'COMPARE ',totlum_saved(1),'to ',totlum_saved(2)

 if (id==master) write(*,"(/,a)") '<-- LUMINOSITY TEST COMPLETE'

#else
 if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF LIGHTCURVE (need -DLIGHTCURVE)'
#endif

end subroutine test_lum

end module testlum
