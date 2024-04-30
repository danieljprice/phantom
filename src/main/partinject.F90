!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module partinject
!
! This module contains routines to inject/move particles in the simulation.
!
!  These routines were previously in the part module, but it had to be
!  separated because of a circular dependency with the timestep_ind module.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: cons2prim, cooling_ism, dim, eos, extern_gr, io,
!   metric_tools, options, part, timestep_ind
!
 implicit none

 public :: add_or_update_particle, add_or_update_sink
 public :: update_injected_particles

 !
 ! Use this flag if particles are updated rather than injected (e.g. inject_sne)
 ! see inject_sne for use; currently only valid for gas particles
 !
 logical, public :: updated_particle = .false.

 private

contains

!-----------------------------------------------------------------------
!+
!  Inject or update a particle into the simulation.
!+
!-----------------------------------------------------------------------
subroutine add_or_update_particle(itype,position,velocity,h,u,particle_number,npart,npartoftype,xyzh,vxyzu,JKmuS)
 use part, only:maxp,iamtype,iphase,maxvxyzu,iboundary,nucleation,eos_vars,abundance
 use part, only:maxalpha,alphaind,maxgradh,gradh,fxyzu,fext,set_particle_type
 use part, only:mhd,Bevol,dBevol,Bxyz,divBsymm!,dust_temp
 use part, only:divcurlv,divcurlB,ndivcurlv,ndivcurlB,ntot,ibin,imu,igamma
 use part, only:iorig,norig
 use io,   only:fatal
 use eos,  only:gamma,gmw
 use dim,  only:ind_timesteps,update_muGamma,h2chemistry
 use timestep_ind, only:nbinmax
 use cooling_ism,  only:abund_default
 integer, intent(in)    :: itype
 real,    intent(in)    :: position(3), velocity(3), h, u
 real,    intent(in), optional :: JKmuS(:)
 integer, intent(in)    :: particle_number
 integer, intent(inout) :: npart, npartoftype(:)
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 logical :: new_particle
 integer :: itype_old

 new_particle = .false.
 if (particle_number == npart+1) then
! This particle doesn't already exist. Create it.
    new_particle = .true.
    npart = npart + 1
    ntot = npart ! reduce_mpi('+',npart)
    if (npart  >  maxp) then
       call fatal('Add particle','npart > maxp')
    endif
    npartoftype(itype) = npartoftype(itype) + 1
    ! add particle ID
    norig                  = norig + 1
    iorig(particle_number) = norig
 elseif (particle_number  >  npart + 1) then
    call fatal('Add particle', 'Incorrect particle number (> npart + 1).')
 elseif (particle_number <= npart) then
    ! This particle already exists. Update it.
    itype_old = iamtype(iphase(particle_number))
    npartoftype(itype_old) = npartoftype(itype_old) - 1
    npartoftype(itype) = npartoftype(itype)+1
 endif

 call set_particle_type(particle_number,itype)

 ! Update particle type, position, size, velocity and energy
 xyzh(1,particle_number) = position(1)
 xyzh(2,particle_number) = position(2)
 xyzh(3,particle_number) = position(3)
 if (itype /= iboundary .or. new_particle) xyzh(4,particle_number) = h

 vxyzu(1,particle_number) = velocity(1)
 vxyzu(2,particle_number) = velocity(2)
 vxyzu(3,particle_number) = velocity(3)
 if (maxvxyzu>=4) vxyzu(4,particle_number) = u

 fxyzu(:,particle_number) = 0.
 fext(:,particle_number) = 0.

 if (mhd) then
    Bevol(:,particle_number) = 0.
    dBevol(:,particle_number) = 0.
    Bxyz(:,particle_number) = 0.
    divBsymm(particle_number) = 0.
 endif

 if (ndivcurlv > 0) divcurlv(:,particle_number) = 0.
 if (ndivcurlB > 0) divcurlB(:,particle_number) = 0.
 if (maxalpha==maxp) alphaind(:,particle_number) = 0.
 if (maxgradh==maxp) gradh(:,particle_number) = 0.
 !if (store_dust_temperature) dust_temp(:,particle_number) = 0.

 if (ind_timesteps) ibin(particle_number) = nbinmax
 if (present(jKmuS)) nucleation(:,particle_number) = JKmuS(:)
 if (update_muGamma) then
    eos_vars(imu,particle_number) = gmw
    eos_vars(igamma,particle_number) = gamma
 endif
 if (h2chemistry) abundance(:,particle_number) = abund_default

end subroutine add_or_update_particle

!-----------------------------------------------------------------------
!+
!  Inject or update a sink particle into the simulation.
!+
!-----------------------------------------------------------------------
subroutine add_or_update_sink(position,velocity,radius,mass,sink_number)
 use io,   only:fatal
 use part, only:nptmass,maxptmass,xyzmh_ptmass,vxyz_ptmass
 real,    intent(in) :: position(3), velocity(3), radius, mass
 integer, intent(in) :: sink_number

 if (sink_number == nptmass+1) then
    ! This sink particle doesn't already exists. Create it.
    nptmass = nptmass + 1
    if (nptmass  >  maxptmass) then
       call fatal('Add sink','nptmass > maxptmass')
    endif
 elseif (sink_number  >  nptmass+1) then
    call fatal('Add sink', 'Incorrect sink number (> maxptmass + 1).')
 endif

 ! Update sink position, size, velocity and mass
 xyzmh_ptmass(1,sink_number) = position(1)
 xyzmh_ptmass(2,sink_number) = position(2)
 xyzmh_ptmass(3,sink_number) = position(3)
 xyzmh_ptmass(4,sink_number) = mass
 xyzmh_ptmass(5,sink_number) = radius
 vxyz_ptmass(1,sink_number) = velocity(1)
 vxyz_ptmass(2,sink_number) = velocity(2)
 vxyz_ptmass(3,sink_number) = velocity(3)

end subroutine add_or_update_sink

!-----------------------------------------------------------------------
!+
!  Update additional quantities needed on injected particles
!  In GR code this is a prim2cons solve
!  otherwise we set the twas variable correctly for individual timesteps
!+
!-----------------------------------------------------------------------
subroutine update_injected_particles(npartold,npart,istepfrac,nbinmax,time,dtmax,dt,dtinject)
 use dim,          only:ind_timesteps
 use timestep_ind, only:get_newbin,change_nbinmax,get_dt
 use part,         only:twas,ibin,ibin_old,iphase,igas,iunknown
#ifdef GR
 use part,         only:xyzh,vxyzu,pxyzu,dens,metrics,metricderivs,fext
 use cons2prim,    only:prim2consall
 use metric_tools, only:init_metric,imet_minkowski,imetric
 use extern_gr,    only:get_grforce_all
 use options,      only:iexternalforce
#endif
 integer,         intent(in)    :: npartold,npart
 integer,         intent(inout) :: istepfrac
 integer(kind=1), intent(inout) :: nbinmax
 real,            intent(inout) :: dt
 real,            intent(in)    :: time,dtmax,dtinject
 integer                        :: i
 integer(kind=1)                :: nbinmaxprev
#ifdef GR
 real                           :: dtext_dum
#endif
 !
 !--Exit if particles not added or updated
 !
 if (npartold==npart .and. .not.updated_particle) return

#ifdef GR
 !
 ! after injecting particles, reinitialise metrics on all particles
 !
 call init_metric(npart,xyzh,metrics,metricderivs)
 call prim2consall(npart,xyzh,metrics,vxyzu,dens,pxyzu,use_dens=.false.)
 if (iexternalforce > 0 .and. imetric /= imet_minkowski) then
    call get_grforce_all(npart,xyzh,metrics,metricderivs,vxyzu,dens,fext,dtext_dum) ! Not 100% sure if this is needed here
 endif
#endif

 if (ind_timesteps) then
    ! find timestep bin associated with dtinject
    nbinmaxprev = nbinmax
    call get_newbin(dtinject,dtmax,nbinmax,allow_decrease=.false.)
    if (nbinmax > nbinmaxprev) then ! update number of bins if needed
       call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
    endif
    ! put all injected particles on shortest bin
    do i=npartold+1,npart
       ibin(i)     = nbinmax
       ibin_old(i) = nbinmax ! for particle waking to ensure that neighbouring particles are promptly woken
       twas(i)     = time + 0.5*get_dt(dtmax,ibin(i))
    enddo
 else
    ! For global timestepping, reset the timestep, since this is otherwise
    ! not updated until after the call to step.
    dt = min(dt,dtinject)
 endif

 ! if a particle was updated rather than added, reset iphase & set timestep (if individual timestepping)
 if (updated_particle) then
    do i=1,npart
       if (iphase(i) == iunknown) then
          iphase(i) = igas
          if (ind_timesteps) then
             ibin(i)     = nbinmax
             ibin_old(i) = nbinmax
             twas(i)     = time + 0.5*get_dt(dtmax,ibin(i))
          endif
       endif
    enddo
 endif

end subroutine update_injected_particles

end module partinject
