!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  This module crudely initialises two galaxies.  It simply reads in
!  data from a run made using Hydra.  Units being read in are
!  1e5M_sun, kpc, 1e5yr.  The former two remain the same here, while
!  the time-unit is modified such that G == 1.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, datafiles, dim, io, mpiutils, part, physcon,
!    timestep, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxp,maxvxyzu
 use io,           only:master,fatal,warning
 use mpiutils,     only:bcast_mpi
 use timestep,     only:tmax,dtmax
 use physcon,      only:solarm,years,kpc
 use units,        only:set_units,udist,utime
 use part,         only:set_particle_type,igas,istar,idarkmatter,iamtype,iphase
 use boundary,     only:set_boundary
 use datafiles,    only:find_phantom_datafile
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(in)    :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=120)               :: filename
 integer                          :: i,ndark,nstar,ngas,itype,ctrd,ctrs,ctrg,ierr,lu
 real                             :: massdark,massstar,massgas
 real                             :: polykset
 real                             :: utmp(maxp)
 logical                          :: iexist
 !
 ! Open file and read data
 !
 filename = find_phantom_datafile('galaxiesP.dat','galaxy_merger')

 open(newunit=lu,file=filename,status='old',action='read',iostat=ierr)
 if (ierr /= 0) call fatal('setup','unable to open '//trim(filename))
 read(lu,*) npart,ndark,nstar,ngas,massdark,massstar,massgas,time
 ctrd = 0
 ctrs = 0
 ctrg = 0
 do i = 1,npart
    read(lu,*) itype,xyzh(1:3,i),vxyzu(1:3,i),utmp(i),xyzh(4,i)
    if (itype== 0) then
       call set_particle_type(i,idarkmatter)
       ctrd = ctrd + 1
    else if (itype==-1) then
       call set_particle_type(i,istar)
       ctrs = ctrs + 1
    else if (itype== 1) then
       call set_particle_type(i,igas)
       ctrg = ctrg + 1
    endif
 enddo
 close(lu)
 if (ctrd/=ndark) call fatal('setup','read in incorrect number of dark matter particles')
 if (ctrs/=nstar) call fatal('setup','read in incorrect number of star particles')
 if (ctrg/=ngas ) call fatal('setup','read in incorrect number of gas particles')
 npartoftype = 0
 npartoftype(idarkmatter) = ndark
 npartoftype(istar) = nstar
 npartoftype(igas) = ngas
 massoftype = 0.0
 massoftype(idarkmatter) = massdark
 massoftype(istar) = massstar
 massoftype(igas) = massgas
 !
 ! set units
 !
 call set_units(dist=kpc,mass=1.0d5*solarm,G = 1.0d0)
 !
 ! reset velocity with new time unit (read in with t = 1d5years)
 !
 vxyzu(1:3,:) = vxyzu(1:3,:)*utime/(1.0d5*years)
 !
 ! set energies (if not isothermal)
 !
 polykset = 3.0d5*utime/udist
 polyk = polykset**2
 gamma = 5./3.
 if (maxvxyzu >= 4) then
    do i = 1,npart
       if (iamtype(iphase(i))==igas) then
          if (utmp(i)==0) then
             vxyzu(4,i) = polyk/(gamma * (gamma-1.0))
          else
             vxyzu(4,i) = utmp(i)*(utime/(1.0d5*years))**2
          endif
       else
          vxyzu(4,i) = 0.0
       endif
    enddo
 endif
 print*,' polyk = ',polyk
 call bcast_mpi(polykset)
 !
 ! set general parameters (only if not already done so)
 !
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 time        = time*(1.0d5*years)/utime
 if (.not. iexist) then
    tmax      = (0.1500d10*years)/utime
    dtmax     = (0.0005d10*years)/utime
 endif

 print*, 'n_total,n_dark,n_star,n_gas: ',npart,ndark,nstar,ngas
 print*, 'm_dark,m_star,m_gas: ',massoftype(idarkmatter), massoftype(istar), massoftype(igas)

end subroutine setpart

end module setup
