!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: deriv
!
!  DESCRIPTION:
!  this module is a wrapper for the main derivative evaluation
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: densityforce, dim, externalforces, forces, forcing,
!    growth, io, krome_interface, linklist, mpiutils, part, photoevap,
!    ptmass, timestep, timing
!+
!--------------------------------------------------------------------------
module deriv
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 public :: derivs
 real, private :: stressmax

 private

contains

!-------------------------------------------------------------
!+
!  calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!+
!-------------------------------------------------------------
subroutine derivs(icall,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
                  dustfrac,ddustevol,temperature,time,dt,dtnew)
 use dim,            only:maxvxyzu
 use io,             only:iprint,fatal
 use linklist,       only:set_linklist
 use densityforce,   only:densityiterate
 use timestep,       only:dtcourant,dtforce,dtmax
 use ptmass,         only:ipart_rhomax
 use externalforces, only:externalforce
#ifdef DRIVING
 use forcing,        only:forceit
#endif
#ifdef PHOTO
 use photoevap,      only:find_ionfront,photo_ionize
 use part,           only:massoftype
#endif
#ifdef DUSTGROWTH
 use growth,         only:get_growth_rate
 use part,           only:St,csound
#endif
 use part,         only:mhd,gradh,alphaind,igas
 use timing,       only:get_timings
 use forces,       only:force
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
 real,         intent(in)    :: dustfrac(:,:)
 real,         intent(inout) :: dustprop(:,:)
 real,         intent(out)   :: ddustevol(:,:),ddustprop(:,:)
 real,         intent(inout) :: temperature(:)
 real,         intent(in)    :: time,dt
 real,         intent(out)   :: dtnew
 logical, parameter :: itiming = .true.
 real(kind=4)       :: t1,tcpu1,tlast,tcpulast

 t1 = 0.
 tcpu1 = 0.
 if (itiming) call get_timings(t1,tcpu1)
 tlast = t1
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
 if (icall==1 .or. icall==0) call set_linklist(npart,nactive,xyzh,vxyzu)

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
                        stressmax,fxyzu,fext,alphaind,gradh)
    call do_timing('dens',tlast,tcpulast)
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
 call force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,dustprop,ddustprop,&
            dustfrac,ddustevol,ipart_rhomax,dt,stressmax,temperature)
 call do_timing('force',tlast,tcpulast)
#ifdef DUSTGROWTH
 !
 ! compute growth rate of dust particles with respect to their positions
 !
 call get_growth_rate(npart,xyzh,vxyzu,csound,St,dustprop,ddustprop(1,:))!--we only get ds/dt (i.e 1st dimension of ddustprop)
#endif
!
! set new timestep from Courant/forces condition
!
 dtnew = min(dtforce,dtcourant,dtmax)

 call do_timing('total',t1,tcpu1,lunit=iprint)

 return

contains

!-------------------------------------------------------------
!+
!  interface to timing routines
!+
!-------------------------------------------------------------
subroutine do_timing(label,tlast,tcpulast,start,lunit)
 use io,       only:iverbose,id,master
 use mpiutils, only:reduce_mpi
 use timing,   only:timer,increment_timer,log_timing,timer_dens,timer_force,timer_link
 character(len=*), intent(in)    :: label
 real(kind=4),     intent(inout) :: tlast,tcpulast
 logical,          intent(in), optional :: start
 integer,          intent(in), optional :: lunit
 real(kind=4) :: t2,tcpu2,tcpu,twall

 call get_timings(t2,tcpu2)
 twall = reduce_mpi('+',t2-tlast)
 tcpu  = reduce_mpi('+',tcpu2-tcpulast)

 if (label=='dens') then
    call increment_timer(timer_dens,t2-tlast,tcpu2-tcpulast)
 else if (label=='force') then
    call increment_timer(timer_force,t2-tlast,tcpu2-tcpulast)
 else if (label=='link') then
    call increment_timer(timer_link,t2-tlast,tcpu2-tcpulast)
 endif

 if (itiming .and. iverbose >= 2) then
    if (id==master) then
       if (present(start)) then
          call log_timing(label,t2-tlast,tcpu,start=.true.)
       elseif (present(lunit)) then
          call log_timing(label,t2-tlast,tcpu,iunit=lunit)
       else
          call log_timing(label,t2-tlast,tcpu)
       endif
    endif
 endif
 tlast = t2
 tcpulast = tcpu2

end subroutine do_timing

end subroutine derivs

end module deriv
