!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: partinject
!
!  DESCRIPTION:
!  This module contains routines to inject/move particles in the simulation.
!
!  These routines were previously in the part module, but it had to be
!  separated because of a circular dependency with the timestep_ind module.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: cons2prim, extern_gr, io, metric_tools, options, part,
!    timestep_ind
!+
!--------------------------------------------------------------------------
module partinject
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 public :: add_or_update_particle, add_or_update_sink
 public :: update_injected_particles
 private

contains

!-----------------------------------------------------------------------
!+
!  Inject or update a particle into the simulation.
!+
!-----------------------------------------------------------------------
subroutine add_or_update_particle(itype,position,velocity,h,u,particle_number,npart,npartoftype,xyzh,vxyzu)
 use part, only:maxp,iamtype,iphase,maxvxyzu
 use part, only:maxalpha,alphaind,maxgradh,gradh,fxyzu,fext,set_particle_type
 use part, only:mhd,Bevol,dBevol,Bxyz,divBsymm
 use part, only:divcurlv,divcurlB,ndivcurlv,ndivcurlB,ntot
 use io,   only:fatal
#ifdef IND_TIMESTEPS
 use part,         only:ibin
 use timestep_ind, only:nbinmax
#endif
 integer, intent(in)    :: itype
 real,    intent(in)    :: position(3), velocity(3), h, u
 integer, intent(in)    :: particle_number
 integer, intent(inout) :: npart, npartoftype(:)
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)

 integer :: itype_old

 if (particle_number == npart+1) then
    ! This particle doesn't already exist. Create it.
    npart = npart + 1
    ntot = npart ! reduce_mpi('+',npart)
    if (npart  >  maxp) then
       call fatal('Add particle','npart > maxp')
    endif
    npartoftype(itype) = npartoftype(itype) + 1
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
 xyzh(4,particle_number) = h
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
#ifdef IND_TIMESTEPS
 ibin(particle_number) = nbinmax
#endif
end subroutine

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
end subroutine

!-----------------------------------------------------------------------
!+
!  Update additional quantities needed on injected particles
!  In GR code this is a prim2cons solve
!  otherwise we set the twas variable correctly for individual timesteps
!+
!-----------------------------------------------------------------------
subroutine update_injected_particles(npartold,npart,istepfrac,nbinmax,time,dtmax,dt,dtinject)
#ifdef IND_TIMESTEPS
 use timestep_ind, only:get_newbin,change_nbinmax,get_dt
 use part,         only:twas,ibin
#endif
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
#ifdef IND_TIMESTEPS
 integer(kind=1) :: nbinmaxprev
 integer :: i
#endif
#ifdef GR
 real :: dtext_dum
#endif

 if (npartold==npart) return
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
#ifdef IND_TIMESTEPS
 ! find timestep bin associated with dtinject
 nbinmaxprev = nbinmax
 call get_newbin(dtinject,dtmax,nbinmax,allow_decrease=.false.)
 if (nbinmax > nbinmaxprev) then ! update number of bins if needed
    call change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
 endif
 ! put all injected particles on shortest bin
 do i=npartold+1,npart
    ibin(i) = nbinmax
    twas(i) = time + 0.5*get_dt(dtmax,ibin(i))
 enddo
#endif

end subroutine update_injected_particles

end module partinject
