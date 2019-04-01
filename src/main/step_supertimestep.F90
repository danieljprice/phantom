!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: supertimestep
!
!  DESCRIPTION: This is the module to control super-timestepping.
!               If use_sts=true, then this is called from evolve.f rather
!               than step.
!               We determine if super-timestepping is required, and if
!               so, we implement it.  This method weakens the stability
!               inherent in the predictor-corrector method, thus megasteps
!               are implemented when necessary; a megastep in dt/Nmega, then
!               superstepping is calculated on this revised timestep.
!               Nmega is predicted by limiting the growth of B in any given
!               step.  Thus, the options are
!               1) no super-timestepping required
!               2) can use N ~ sqrt(dt/dtdiff)
!               3) increased N using Nmegatseps
!               4) using Nreal ~ dt/didiff (since Nsupersteps = Nreal for small Nreal)
!               5) using Nreal since Nsts > nnu, or Nsts*Nmega>Nreal
!               The subroutines called are in utils_supertimestep.f
!
!  REFERENCES: Alexiades V., Amiez G., Gremaud P.A., 1996, Commun. Numer. Meth. Eng., 12, 31
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: evwrite, io, io_summary, part, step_lf_global, timestep,
!    timestep_ind, timestep_sts
!+
!--------------------------------------------------------------------------
module supertimestep
 use io_summary,     only: summary_counter,iosum_nsts
 use timestep_sts,   only: iNosts,iNsts,iNmegasts,iNostsSml,iNostsBig,icase_sts, &
                           nbinmaxsts,Nmegasts_now,Nmegasts_next,Nreal,Nsts

 implicit none
 !--Local variables required to be saved
 integer         :: ists_loopctr
 !--Subroutines
 public          :: step_sts

contains
!--------------------------------------------------------------------------
subroutine step_sts(npart,nactive,time,dt,dtextforce,dtnew,iprint)
#ifdef IND_TIMESTEPS
 use timestep_sts,   only: sts_init_step
 use timestep_ind,   only: nbinmax
 use part,           only: ibin,iphase
#else
 use timestep,       only: dtcourant,dtforce,dterr
 use timestep_sts,   only: sts_get_dtau_array
#endif
 use timestep,       only: dtdiff
 use timestep_sts,   only: sts_it_n,dtau,Nmegasts_done,bigdt, &
                           sts_initialise_activity,sts_set_active_particles,nnu
 use io,             only: fatal
 use step_lf_global, only: step
 integer, intent(inout) :: npart,nactive
 integer, intent(in)    :: iprint
 real,    intent(in)    :: time,dt
 real,    intent(inout) :: dtextforce
 real,    intent(out)   :: dtnew
 real, parameter        :: time_tol = 1.0d-8
 integer                :: nactive_sts
#ifndef IND_TIMESTEPS
 integer(kind=1), parameter :: nbinmax = 0  ! The simplest way to keep everything clean
#endif
 real                   :: dtau0,dtdiff_in,dtsum,timei
 real                   :: dttemp3(3),dttemp2(2)
 !
 ! Initialise values
 !
 dtdiff_in   = dtdiff
 dtdiff      = bigdt
 dttemp2     = bigdt
 dttemp3     = bigdt
 timei       = time
 dtsum       = 0.0d0
 sts_it_n    = .true.
#ifdef IND_TIMESTEPS
 ! Determine activity of particles
 ! This is required even if sts is not being used since
 ! isactive(:) is always being checked when use_sts=.true.
 ! This does not affect the activity of particles, just determines who
 ! will be and if they're normally active or sts-active
 call sts_initialise_activity(nactive_sts,npart,ibin,iphase)
#else
 Nmegasts_done = 0
 call sts_get_dtau_array(Nmegasts_next,dt,dtdiff_in)
 nactive_sts = npart
#endif
 !
 ! No superstepping required; just call step and exit subroutine
 if (icase_sts == iNosts) then
    call step(npart,nactive,time,dt,dtextforce,dtnew)
    call summary_counter(iosum_nsts)
    return
 endif
 !
 ! Print statements to inform user status of super-timestepping
 call sts_print_output(nactive_sts,time,nbinmax,iprint)
 !
 ! Set only sts-active particles to active.  The logical input is to set only
 ! sts-active particles to active.
#ifdef IND_TIMESTEPS
 call sts_set_active_particles(npart,nactive,.false.)
#endif
 !
 ! Perform super-timestepping
 !
 sts_it_n     = .false.
 Nmegasts_now = Nmegasts_next
 ists_loopctr = 0
 do while (ists_loopctr < Nmegasts_now)
    Nmegasts_done = Nmegasts_done + 1
    ists_loopctr  = ists_loopctr  + 1
    dtau0         = dtau(ists_loopctr)
    if ( ists_loopctr==Nmegasts_now ) then
#ifdef IND_TIMESTEPS
       ! For the final iteration, move all active particles that were defined
       ! by set_active_particles in evolve.  The logical input is to set all
       ! active particle to active.
       call sts_set_active_particles(npart,nactive,.true.)
#endif
       sts_it_n = .true.
    endif
    !
    ! Call step
    call step(npart,nactive,timei,dtau0,dtextforce,dtnew)
    call summary_counter(iosum_nsts)
    dtsum      = dtsum + dtau0
    timei      = timei + dtau0

    dttemp2(1) = min(dttemp2(1),dtnew     )
    dttemp2(2) = min(dttemp2(2),dtextforce)
#ifndef IND_TIMESTEPS
    dttemp3(1) = min(dttemp3(1),dtcourant )
    dttemp3(2) = min(dttemp3(2),dtforce   )
    dttemp3(3) = min(dttemp3(3),dterr     )
#endif
 enddo

 dtnew      = dttemp2(1)
 dtextforce = dttemp2(2)
#ifndef IND_TIMESTEPS
 dtcourant  = dttemp3(1)
 dtforce    = dttemp3(2)
 dterr      = dttemp3(3)
#endif
 if (abs(1.0-dtsum/dt) >= (sqrt(real(Nmegasts_now))*time_tol)) then
    write(iprint,'(a,I4,3Es16.7)') 'Super-timestepping: Nmegasts,dt,dt_sum,abs(1.0-dt_sum/dt)  : ' &
    ,Nmegasts_now,dt,dtsum,abs(1.0-dtsum/dt)
    call fatal ('step','Super-timestepping: superstep not reaching the correct time')
 endif

end subroutine step_sts
!--------------------------------------------------------------------------
!+
!  A routine to handle printing to output and the calls to summary
!+
!--------------------------------------------------------------------------
subroutine sts_print_output(nactive_sts,time,nbinmax,iprint)
 use io,             only: fatal,iverbose
 use io_summary,     only: summary_variable, &
                           iosumstse,iosumsts,iosumstsi,iosumstsm,iosumstsr,iosumstsri, &
                           iosumstsd,iosumstsdi,iosumstso,iosumstsoi, &
                           iosumstsnn,iosumstsnm,iosumstsns,iosumstsnl
 integer,         intent(in)  :: nactive_sts,iprint
 integer(kind=1), intent(in)  :: nbinmax
 real,            intent(in)  :: time
 integer                      :: Nviaibin
 !
 ! Print to screen, if requested
 if (iverbose>=1) then
    if (icase_sts==iNsts) then
       write(iprint,10) 'Super-timestepping: Enabled with Nsts = ',Nmegasts_next,' at time = ',time,' with Nreal = ',Nreal
    else if (icase_sts==iNmegasts) then
       write(iprint,10) 'Super-timestepping: Enabled with Nsts*Nmega = ',Nmegasts_next,' at time = ',time,' with Nreal = ',Nreal
    else if (icase_sts==iNostsSml .or. icase_sts==iNostsBig) then
       write(iprint,20)      &
       'Super-timestepping: Disabled.  Using Nreal = ',Nreal,' at time = ',time
    else
       ! safe since this call should never happen
       call fatal('Super-timestepping','Illegal case.  You should not be able to trigger this.')
    endif
 endif
10 format(a,I4,a,Es16.7,a,I4)
20 format(a,I4,a,Es16.7)

 Nviaibin = 2**(nbinmaxsts-nbinmax)
 ! Update summary
 if (icase_sts==iNsts) then
    call summary_variable('sts',iosumstse ,0,real(Nmegasts_next))
    call summary_variable('sts',iosumsts  ,0,real(Nreal)        )
    call summary_variable('sts',iosumstsnn,0,real(nactive_sts)  )
    call summary_variable('sts',iosumstsi ,0,real(Nviaibin)     )
 else if (icase_sts==iNmegasts) then
    call summary_variable('sts',iosumstsm ,0,real(Nmegasts_next))
    call summary_variable('sts',iosumstsr ,0,real(Nreal)        )
    call summary_variable('sts',iosumstsnm,0,real(nactive_sts)  )
    call summary_variable('sts',iosumstsri,0,real(Nviaibin)     )
 else if (icase_sts==iNostsSml) then
    call summary_variable('sts',iosumstsd ,0,real(Nreal)        )
    call summary_variable('sts',iosumstsns,0,real(nactive_sts)  )
    call summary_variable('sts',iosumstsdi,0,real(Nviaibin)     )
 else if (icase_sts==iNostsBig) then
    call summary_variable('sts',iosumstso ,0,real(Nreal)        )
    call summary_variable('sts',iosumstsnl,0,real(nactive_sts)  )
    call summary_variable('sts',iosumstsoi,0,real(Nviaibin)     )
 else
    call fatal('Super-timestepping', 'BROKEN.  Illegal case')
 endif

end subroutine sts_print_output
!--------------------------------------------------------------------------
end module supertimestep
