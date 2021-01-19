!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!   forces, forcing, growth, io, linklist, part, photoevap, ptmass,
!   ptmass_radiation, timestep, timestep_ind, timing
!
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

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
                  dustevol,ddustevol,dustfrac,eos_vars,time,dt,dtnew,pxyzu,dens,metrics)
 use dim,            only:maxvxyzu
 use io,             only:iprint,fatal
 use linklist,       only:set_linklist
 use densityforce,   only:densityiterate
 use ptmass,         only:ipart_rhomax
 use externalforces, only:externalforce
 use part,           only:dustgasprop,gamma_chem,dvdx,Bxyz
#ifdef IND_TIMESTEPS
 use timestep_ind,   only:nbinmax
#else
 use timestep,       only:dtcourant,dtforce,dtrad
#endif
 use timestep,       only:dtmax
#ifdef DRIVING
 use forcing,        only:forceit
#endif
#ifdef PHOTO
 use photoevap,      only:find_ionfront,photo_ionize
 use part,           only:massoftype
#endif
#ifdef DUSTGROWTH
 use growth,         only:get_growth_rate
 use part,           only:VrelVf
#endif
#ifdef SINK_RADIATION
 use ptmass_radiation, only:get_dust_temperature_from_ptmass
 use part,             only:dust_temp,nptmass,xyzmh_ptmass
#endif
#ifdef PERIODIC
 use ptmass,         only:ptmass_boundary_crossing
 use part,           only:nptmass,xyzmh_ptmass
#endif
 use part,           only:mhd,gradh,alphaind,igas
 use timing,         only:get_timings
 use forces,         only:force
 use part,           only:iradxi,ifluxx,ifluxy,ifluxz,ithick
 use derivutils,     only:do_timing
#ifdef GR
 use cons2prim,      only:cons2primall
#endif
 use cons2prim,      only:cons2prim_everything
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
 real,         intent(in)    :: rad(:,:)
 real,         intent(out)   :: eos_vars(:,:)
 real,         intent(out)   :: drad(:,:)
 real,         intent(inout) :: radprop(:,:)
 real,         intent(in)    :: dustevol(:,:)
 real,         intent(inout) :: dustprop(:,:)
 real,         intent(out)   :: dustfrac(:,:)
 real,         intent(out)   :: ddustevol(:,:),ddustprop(:,:)
 real,         intent(in)    :: time,dt
 real,         intent(out)   :: dtnew
 real,         intent(inout) :: pxyzu(:,:), dens(:)
 real,         intent(in)    :: metrics(:,:,:,:)
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
#ifdef PERIODIC
    if (nptmass > 0) call ptmass_boundary_crossing(nptmass,xyzmh_ptmass)
#endif
 endif

 call do_timing('link',tlast,tcpulast,start=.true.)

#ifdef PHOTO
 !
 ! update location of particles on grid and calculate the location of the ionization front
 !
 call find_ionfront(time,npart,xyzh,massoftype(igas))
 !
 ! update the temperatures of the particles depending on whether ionized or not
 !
 call photo_ionize(vxyzu,npart)
#endif
!
! calculate density by direct summation
!
 if (icall==1) then
    call densityiterate(1,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                        stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx)
    call do_timing('dens',tlast,tcpulast)
 endif

#ifdef GR
 call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
#else
 call cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,gamma_chem,Bevol,Bxyz,dustevol,dustfrac,alphaind)
#endif

!
! compute forces
!
#ifdef DRIVING
 ! forced turbulence -- call driving routine
 call forceit(time,npart,xyzh,vxyzu,fxyzu)
 call do_timing('driving',tlast,tcpulast)
#endif
 stressmax = 0.
 call force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&
            rad,drad,radprop,dustprop,dustgasprop,dustfrac,ddustevol,&
            ipart_rhomax,dt,stressmax,eos_vars,dens,metrics)
 call do_timing('force',tlast,tcpulast)

#ifdef DUSTGROWTH
 ! compute growth rate of dust particles
 call get_growth_rate(npart,xyzh,vxyzu,dustgasprop,VrelVf,dustprop,ddustprop(1,:))!--we only get ds/dt (i.e 1st dimension of ddustprop)
#endif

#ifdef SINK_RADIATION
 !compute dust temperature
 if (maxvxyzu >= 4) call get_dust_temperature_from_ptmass(npart,xyzh,nptmass,xyzmh_ptmass,dust_temp)
#endif
!
! set new timestep from Courant/forces condition
!
#ifdef IND_TIMESTEPS
 dtnew = dtmax/2**nbinmax  ! minimum timestep over all particles
#else
 dtnew = min(dtforce,dtcourant,dtrad,dtmax)
#endif


 call do_timing('total',t1,tcpu1,lunit=iprint)

 return

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
!+
!--------------------------------------
subroutine get_derivs_global(tused,dt_new)
 use part,   only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
                dustfrac,ddustevol,eos_vars,pxyzu,dens,metrics,dustevol
 use timing, only:printused,getused
 use io,     only:id,master
 real(kind=4), intent(out), optional :: tused
 real,         intent(out), optional :: dt_new
 real(kind=4) :: t1,t2
 real :: dtnew
 real :: time,dt

 time = 0.
 dt = 0.
 call getused(t1)
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
             rad,drad,radprop,dustprop,ddustprop,dustevol,ddustevol,dustfrac,eos_vars,&
             time,dt,dtnew,pxyzu,dens,metrics)
 call getused(t2)
 if (id==master .and. present(tused)) call printused(t1)
 if (present(tused)) tused = t2 - t1
 if (present(dt_new)) dt_new = dtnew

end subroutine get_derivs_global

end module deriv
