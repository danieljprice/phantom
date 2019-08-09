!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testsedov
!
!  DESCRIPTION:
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, deriv, dim, energies, eos, evolve, evwrite,
!    initial_params, io, io_summary, mpiutils, options, part, physcon,
!    testutils, timestep, unifdis, viscosity
!+
!--------------------------------------------------------------------------
module testsedov
 implicit none

 public :: test_sedov

contains

subroutine test_sedov(ntests,npass)
 use dim,      only:maxp,maxvxyzu,maxalpha,use_dust
 use io,       only:id,master,iprint,ievfile,iverbose,real4
 use boundary, only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use unifdis,  only:set_unifdis
 use part,     only:mhd,npart,npartoftype,massoftype,xyzh,vxyzu,hfact,fxyzu,fext,ntot, &
                    divcurlv,divcurlB,Bevol,dBevol,Bextx,Bexty,Bextz,alphaind,&
                    dustfrac,ddustevol,dustevol,dustprop,ddustprop,temperature
 use part,     only:iphase,maxphase,igas,isetphase
 use eos,      only:gamma,polyk
 use options,  only:ieos,tolh,alpha,alphau,alphaB,beta
 use physcon,  only:pi
 use deriv,    only:derivs
 use timestep, only:time,tmax,dtmax,C_cour,C_force,dt,tolv
#ifndef IND_TIMESTEPS
 use timestep, only:dtcourant,dtforce
#endif
 use testutils, only:checkval,update_test_scores
 use evwrite,   only:init_evfile,write_evfile
 use energies,  only:etot,totmom,angtot,mdust
 use evolve,    only:evol
 use viscosity, only:irealvisc
 use io_summary,only:summary_reset
 use initial_params, only:etot_in,angtot_in,totmom_in,mdust_in
 use mpiutils,  only:reduceall_mpi
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2)
 integer :: i,itmp,ierr,iu
 real    :: psep,denszero,enblast,rblast,prblast,gam1,dtext_dum
 real    :: totmass,etotin,momtotin,etotend,momtotend
 character(len=20) :: logfile,evfile,dumpfile

#ifndef PERIODIC
 if (id==master) write(*,"(/,a)") '--> SKIPPING Sedov blast wave (needs -DPERIODIC)'
 return
#endif
#ifdef DISC_VISCOSITY
 if (id==master) write(*,"(/,a)") '--> SKIPPING Sedov blast wave (cannot use -DDISC_VISCOSITY)'
 return
#endif

 testsedv: if (maxvxyzu >= 4) then
    if (id==master) write(*,"(/,a)") '--> testing Sedov blast wave'
    call summary_reset ! reset since summary will be written by evol if there are warnings; want only warnings from this test
    call set_boundary(-0.5,0.5,-0.5,0.5,-0.5,0.5)
    time      = 0.
    hfact     = 1.2
    ieos      = 2
    iverbose  = 1 !max(iverbose,1)
    alpha     = 1.
    alphau    = 1.
    alphaB    = 0.
    beta      = 2.
    tolh      = 1.e-5
    if (maxalpha==maxp) alphaind(1,:) = real4(alpha)
    irealvisc = 0
    tolv      = 1.e-3
    iu        = 4
!
!--setup particles
!
    npart = 16
    psep  = dxbound/npart

    denszero = 1.0
    polyk    = 0.
    enblast  = 1.0
    rblast   = 2.*hfact*psep
    gamma    = 5./3.
    gam1     =  gamma - 1.
    prblast  = gam1*enblast/(4./3.*pi*rblast**3)
    npart    = 0

    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,hfact,npart,xyzh)

    npartoftype(:) = 0
    npartoftype(1) = npart
    ntot           = npart

    totmass = denszero*dxbound*dybound*dzbound
    massoftype(:) = 0.
    massoftype(igas) = totmass/reduceall_mpi('+',npart)
    print*,' npart = ',npart,' particle mass = ',massoftype(igas)

    do i=1,npart
       if (maxphase==maxp) iphase(i) = isetphase(igas,iactive=.true.)
       vxyzu(:,i) = 0.
       if (use_dust) then
          dustfrac(:,i) = 0.
          dustevol(:,i) = 0.
       endif

       if ((xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2) < rblast*rblast) then
          vxyzu(iu,i) = prblast/(gam1*denszero)
       else
          vxyzu(iu,i) = 0.
       endif
       if (mhd) then
          Bevol(:,i) = 0.
          Bextx = 0.
          Bexty = 0.
          Bextz = 0.
       endif
    enddo
    tmax = 0.1
    dtmax = tmax
!
!--call derivs the first time around
!
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtext_dum)
!
!--now call evolve
!
#ifndef IND_TIMESTEPS
    dt = min(C_cour*dtcourant,C_force*dtforce)
#endif
    C_cour = 0.15
    C_force = 0.25
    iprint = 6
    logfile  = 'test01.log'
    evfile   = 'test01.ev'
    dumpfile = 'test001'

    call init_evfile(ievfile,evfile,.true.)
    call write_evfile(time,dt)
    etotin    = etot
    momtotin  = totmom
    etot_in   = etot
    angtot_in = angtot
    totmom_in = totmom
    mdust_in  = mdust
    call evol('test.in',logfile,evfile,dumpfile)
    call write_evfile(time,dt)
    etotend   = etot
    momtotend = totmom

    nfailed(:) = 0
    call checkval(etotend,etotin,4.7e-4,nfailed(1),'total energy')
    call checkval(momtotend,momtotin,7.e-15,nfailed(2),'linear momentum')

    ! delete temporary files
    close(unit=ievfile,status='delete',iostat=ierr)

    itmp = 201
    open(unit=itmp,file='test002',status='old',iostat=ierr)
    close(unit=itmp,status='delete',iostat=ierr)

    open(unit=itmp,file='test.in',status='old',iostat=ierr)
    close(unit=itmp,status='delete',iostat=ierr)

    call update_test_scores(ntests,nfailed,npass)
 else
    if (id==master) write(*,"(/,a)") '--> SKIPPING Sedov blast wave (needs thermal energy: maxvxyzu=4)'

 endif testsedv

end subroutine test_sedov

end module testsedov
