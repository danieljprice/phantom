!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testlum
!
! Tests lightcurve and timestepping
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: deriv, dim, energies, eos, io, options, part, setdisc,
!   testutils, timing, viscosity
!
 implicit none
 public :: test_lum

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of fake lightcurve output
!+
!-----------------------------------------------------------------------
subroutine test_lum(ntests,npass)
 use dim,      only:periodic,disc_viscosity,maxp,ind_timesteps,track_lum
 use io,       only:id,master
 use io,       only:iverbose
 use part,     only:init_part,npart,npartoftype,massoftype,xyzh,hfact,vxyzu,&
                    igas,iphase,isetphase
 use eos,             only:gamma,polyk
 use testutils,       only:checkval,checkvalf,update_test_scores
 use energies,        only:compute_energies,ekin,etherm,totlum !etot,eacc,accretedmass
 use setdisc,         only:set_disc
 use deriv,           only:get_derivs_global
 use part,            only:alphaind,maxalpha
 use options,         only:alphau,alphaB
 use viscosity,       only:irealvisc,shearfunc,dt_viscosity,shearparam
 use options,         only:ieos,alpha,iexternalforce,ipdv_heating,ishock_heating
 use units,           only:unit_luminosity,set_units
 use physcon,         only:solarm,au
 !use part,            only:luminosity
 integer, intent(inout) :: ntests,npass
 integer                :: i,itest
 real                   :: totlum_saved(2),etot_saved(2),diff,alpha_in
 real                   :: time
 real, parameter        :: tol = 2.e-5
 integer                :: nfail(1),ii
 integer :: nactive

!#ifdef DISC_VISCOSITY
!    if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF LIGHTCURVE (cannot have -DDISC_VISCOSITY)'
!    return
!#endif

 if (.not.track_lum) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF LIGHTCURVE (need track_lum=.true.)'
    return
 endif

 if (periodic) then
    if (id==master) write(*,"(/,a)") '--> SKIPPING LUMINOSITY TESTS (need PERIODIC=no)'
    return
 else
    if (ind_timesteps) then
       if (id==master) write(*,"(/,a,/)") '--> TESTING LUMINOSITY'
    else
       if (id==master) write(*,"(/,a)") '--> SKIPPING TEST OF LIGHTCURVE (need -DIND_TIMESTEPS too)'
       return
    endif
 endif

 call init_part()
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

 call set_units(mass=solarm,dist=au,G=1.d0)

 alpha_in = 0.1
 alpha = alpha_in

! Run with "alpha_SS" viscosity first
 if (.not.disc_viscosity) then
    irealvisc = 2
    shearparam = alpha_in
    alpha = 0.
    alphau = 0.
    alphaB = 0.
    if (maxalpha==maxp) alphaind(:,:) = real(alpha_in,kind=kind(alphaind))
 endif

 call set_disc(id,master=master,nparttot=npart,npart=npart,rmin=0.5,rmax=10.,p_index=1.5,q_index=0.75, &
               HoverR=0.02,disc_mass=1.e-4,star_mass=1.0,gamma=gamma,particle_mass=massoftype(1), &
               hfact=hfact,xyzh=xyzh,vxyzu=vxyzu,polyk=polyk,alpha=alpha,ismooth=.true.,writefile=.false.)

 do ii=1,4
    select case(ii)
    case(1)
       ipdv_heating = 0
       ishock_heating = 1
    case(2)
       ipdv_heating = 1
       ishock_heating = 0
    case(3)
       ipdv_heating = 1
       ishock_heating = 1
    case default
       ipdv_heating = 0
       ishock_heating = 0
    end select

    if (ind_timesteps) then
       do itest = 1,2
          if (itest == 1) then
             ! during test 1, all particles are active
             iphase(1:npart) = isetphase(igas,iactive=.true.)
             nactive = npart
          else
             ! during test 2, only a fraction of particles are active
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
          call get_derivs_global()
          call compute_energies(time)

          totlum_saved(itest) = totlum*unit_luminosity
          etot_saved(itest) = ekin + etherm
       enddo
       !
       ! check that the total luminosity is the same regardless of whether particles are active or not
       ! i.e. that the luminosity is updated on active particles but preserved on non-active particles 
       ! between force updates
       !
       diff = (totlum_saved(1) - totlum_saved(2))/totlum_saved(1)
       call checkval(totlum_saved(1),totlum_saved(2),tol,nfail(1),'total luminosity')
       call update_test_scores(ntests,nfail(1:1),npass)
    endif
    nactive = npart !for later
 enddo

 if (id==master) write(*,"(/,a)") '<-- LUMINOSITY TEST COMPLETE'

end subroutine test_lum

end module testlum
