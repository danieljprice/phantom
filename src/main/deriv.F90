!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module deriv
!
! this module is a wrapper for the main derivative evaluation
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: cons2prim, densityforce, derivutils, dim, externalforces,
!   forces, forcing, growth, io, linklist, metric_tools, options, part,
!   porosity, ptmass, ptmass_radiation, radiation_implicit, timestep,
!   timestep_ind, timing
!
 implicit none

 public :: derivs, get_derivs_global
 real, private :: stressmax

 private

contains

!-------------------------------------------------------------
!+
!  calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!+
!-------------------------------------------------------------
subroutine derivs(icall,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                  Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
                  dustevol,ddustevol,filfac,dustfrac,eos_vars,time,dt,dtnew,pxyzu,&
                  dens,metrics,apr_level)
 use dim,            only:maxvxyzu,mhd,fast_divcurlB,gr,periodic,do_radiation,&
                          sink_radiation,use_dustgrowth,ind_timesteps
 use io,             only:iprint,fatal,error
 use linklist,       only:set_linklist
 use densityforce,   only:densityiterate
 use ptmass,         only:ipart_rhomax,ptmass_calc_enclosed_mass,ptmass_boundary_crossing
 use externalforces, only:externalforce
 use part,           only:dustgasprop,dvdx,Bxyz,set_boundaries_to_active,&
                          nptmass,xyzmh_ptmass,sinks_have_heating,dust_temp,VrelVf,fxyz_drag
 use timestep_ind,   only:nbinmax
 use timestep,       only:dtmax,dtcourant,dtforce,dtrad
#ifdef DRIVING
 use forcing,        only:forceit
#endif
 use growth,           only:get_growth_rate
 use porosity,         only:get_disruption,get_probastick
 use ptmass_radiation, only:get_dust_temperature
 use timing,         only:get_timings
 use forces,         only:force
 use part,           only:mhd,gradh,alphaind,igas,iradxi,ifluxx,ifluxy,ifluxz,ithick
 use derivutils,     only:do_timing
 use cons2prim,      only:cons2primall,cons2prim_everything
 use metric_tools,   only:init_metric
 use radiation_implicit, only:do_radiation_implicit,ierr_failed_to_converge
 use options,        only:implicit_radiation,implicit_radiation_store_drad,use_porosity
 integer,      intent(in)    :: icall
 integer,      intent(inout) :: npart
 integer,      intent(in)    :: nactive
 real,         intent(inout) :: xyzh(:,:)
 real,         intent(inout) :: vxyzu(:,:)
 real,         intent(inout) :: fxyzu(:,:)
 real,         intent(in)    :: fext(:,:)
 real(kind=4), intent(out)   :: divcurlv(:,:)
 real(kind=4), intent(out)   :: divcurlB(:,:)
 real,         intent(in)    :: Bevol(:,:)
 real,         intent(out)   :: dBevol(:,:)
 real,         intent(inout) :: rad(:,:)
 real,         intent(out)   :: eos_vars(:,:)
 real,         intent(out)   :: drad(:,:)
 real,         intent(inout) :: radprop(:,:)
 real,         intent(in)    :: dustevol(:,:)
 real,         intent(inout) :: dustprop(:,:)
 real,         intent(out)   :: dustfrac(:,:)
 real,         intent(out)   :: ddustevol(:,:),ddustprop(:,:)
 real,         intent(inout) :: filfac(:)
 real,         intent(in)    :: time,dt
 real,         intent(out)   :: dtnew
 real,         intent(inout) :: pxyzu(:,:), dens(:)
 real,         intent(inout) :: metrics(:,:,:,:)
 integer(kind=1), intent(in) :: apr_level(:)
 integer                     :: ierr,i
 real(kind=4)                :: t1,tcpu1,tlast,tcpulast

 t1    = 0.
 tcpu1 = 0.
 call get_timings(t1,tcpu1)
 tlast    = t1
 tcpulast = tcpu1
!
!--check for errors in input options
!
 if (icall < 0 .or. icall > 2) call fatal('deriv','invalid icall on input')
!
! icall is a flag to say whether or not positions have changed
! since the last call to derivs.
!
! icall = 1 is the "standard" call to derivs: calculates all derivatives
! icall = 2 does not remake the link list and does not recalculate density
!           (ie. only re-evaluates the SPH force term using updated values
!            of the input variables)
!
! call link list to find neighbours
!
 if (icall==1 .or. icall==0) then
    call set_linklist(npart,nactive,xyzh,vxyzu)

    if (gr) then
       ! Recalculate the metric after moving particles to their new tasks
       call init_metric(npart,xyzh,metrics)
    endif

    if (nptmass > 0 .and. periodic) call ptmass_boundary_crossing(nptmass,xyzmh_ptmass)
 endif

 call do_timing('link',tlast,tcpulast,start=.true.)


 !
 ! compute disruption of dust particles
 !
 if (use_dustgrowth .and. use_porosity) call get_disruption(npart,xyzh,filfac,dustprop,dustgasprop)
!
! calculate density by direct summation
!

 if (icall==1) then
    call densityiterate(1,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                        stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level)
    if (.not. fast_divcurlB) then
       ! Repeat the call to calculate all the non-density-related quantities in densityiterate.
       ! This needs to be separate for an accurate calculation of divcurlB which requires an up-to-date rho.
       ! if fast_divcurlB = .false., then all additional quantities are calculated during the previous call
       call densityiterate(3,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                           stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level)
    endif
    set_boundaries_to_active = .false.     ! boundary particles are no longer treated as active
    call do_timing('dens',tlast,tcpulast)
 endif

 if (gr) then
    call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
 else
    call cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,Bevol,Bxyz,dustevol,dustfrac,alphaind)
 endif
 call do_timing('cons2prim',tlast,tcpulast)

 !
 ! implicit radiation update
 !
 if (do_radiation .and. implicit_radiation .and. dt > 0.) then
    call do_radiation_implicit(dt,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
    if (ierr /= 0 .and. ierr /= ierr_failed_to_converge) call fatal('radiation','Failed in radiation')
    call do_timing('radiation',tlast,tcpulast)
 endif

!
! compute forces
!
#ifdef DRIVING
 ! forced turbulence -- call driving routine
 call forceit(time,npart,xyzh,vxyzu,fxyzu)
 call do_timing('driving',tlast,tcpulast)
#endif
 stressmax = 0.
 if (sinks_have_heating(nptmass,xyzmh_ptmass)) call ptmass_calc_enclosed_mass(nptmass,npart,xyzh)
 call force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&
            rad,drad,radprop,dustprop,dustgasprop,dustfrac,ddustevol,fext,fxyz_drag,&
            ipart_rhomax,dt,stressmax,eos_vars,dens,metrics,apr_level)
 call do_timing('force',tlast,tcpulast)

 if (use_dustgrowth) then ! compute growth rate of dust particles
    call get_growth_rate(npart,xyzh,vxyzu,dustgasprop,VrelVf,dustprop,filfac,ddustprop(1,:))!--we only get dm/dt (i.e 1st dimension of ddustprop)
    ! compute growth rate and probability of sticking/bouncing of porous dust
    if (use_porosity) call get_probastick(npart,xyzh,ddustprop(1,:),dustprop,dustgasprop,filfac)
 endif

!
! compute dust temperature
!
 if (sink_radiation .and. maxvxyzu == 4) then
    call get_dust_temperature(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp)
 endif

 if (do_radiation .and. implicit_radiation .and. .not.implicit_radiation_store_drad) then
    !$omp parallel do shared(drad,fxyzu,npart) private(i)
    do i=1,npart
       drad(:,i) = 0.
       fxyzu(4,i) = 0.
    enddo
    !$omp end parallel do
 endif
!
! set new timestep from Courant/forces condition
!
 if (ind_timesteps) then
    dtnew = dtmax/2.**nbinmax  ! minimum timestep over all particles
 else
    dtnew = min(dtforce,dtcourant,dtrad,dtmax)
 endif

 call do_timing('total',t1,tcpu1,lunit=iprint)

end subroutine derivs

!--------------------------------------
!+
!  wrapper for the call to derivs
!  so only one line needs changing
!  if interface changes
!
!  this should NOT be called during timestepping, it is useful
!  for when one requires just a single call to evaluate derivatives
!  and store them in the global shared arrays
!  does not work for sink GR yet
!+
!--------------------------------------
subroutine get_derivs_global(tused,dt_new,dt,icall)
 use part,         only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                        Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,filfac,&
                        dustfrac,ddustevol,eos_vars,pxyzu,dens,metrics,dustevol,gr,&
                        apr_level
 use timing,       only:printused,getused
 use io,           only:id,master
 use cons2prim,    only:prim2consall
 use metric_tools, only:init_metric
 real(kind=4), intent(out), optional :: tused
 real,         intent(out), optional :: dt_new
 real,         intent(in),  optional :: dt  ! optional argument needed to test implicit radiation routine
 integer     , intent(in),  optional :: icall
 real(kind=4) :: t1,t2
 real    :: dtnew
 real    :: time,dti
 integer :: icalli

 time = 0.
 dti = 0.
 icalli = 1
 if (present(dt)) dti = dt
 if (present(icall)) icalli = icall
 call getused(t1)
 ! update conserved quantities in the GR code
 if (gr) then
    call init_metric(npart,xyzh,metrics)
    call prim2consall(npart,xyzh,metrics,vxyzu,pxyzu,use_dens=.false.,dens=dens)
 endif

 ! evaluate derivatives
 call derivs(icalli,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
             rad,drad,radprop,dustprop,ddustprop,dustevol,ddustevol,filfac,dustfrac,&
             eos_vars,time,dti,dtnew,pxyzu,dens,metrics,apr_level)
 call getused(t2)
 if (id==master .and. present(tused)) call printused(t1)
 if (present(tused)) tused = t2 - t1
 if (present(dt_new)) dt_new = dtnew

end subroutine get_derivs_global

end module deriv
