!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module testsedov
!
! No description
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, deriv, dim, domain, energies, eos, evolve,
!   evwrite, io, io_summary, mpiutils, options, part,
!   physcon, radiation_utils, testutils, timestep, unifdis, units,
!   viscosity
!
 implicit none

 public :: test_sedov

contains
!-----------------------------------------------------------------------
!+
!   Unit test of the "complete" code, performing sedov test
!+
!-----------------------------------------------------------------------
subroutine test_sedov(ntests,npass)
 use dim,      only:maxp,maxvxyzu,maxalpha,use_dust,periodic,do_radiation
 use io,       only:id,master,iprint,ievfile,iverbose,real4
 use boundary, only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use unifdis,  only:set_unifdis
 use part,     only:init_part,npart,npartoftype,massoftype,xyzh,vxyzu,hfact,ntot, &
                    alphaind,rad,radprop,ikappa
 use part,     only:iphase,maxphase,igas,isetphase
 use eos,      only:gamma,polyk
 use options,  only:ieos,tolh,alpha,alphau,alphaB,beta
 use physcon,  only:pi,au,solarm
 use deriv,    only:get_derivs_global
 use timestep, only:time,tmax,dtmax,C_cour,C_force,dt,tolv,bignumber
 use units,    only:set_units
#ifndef IND_TIMESTEPS
 use timestep, only:dtcourant,dtforce,dtrad
#endif
 use testutils, only:checkval,update_test_scores
 use evwrite,   only:init_evfile,write_evfile
 use energies,  only:etot,totmom,angtot,mdust,etot_in,angtot_in,totmom_in,mdust_in
 use evolve,    only:evol
 use viscosity, only:irealvisc
 use io_summary,only:summary_reset
 use mpiutils,  only:reduceall_mpi
 use domain,    only:i_belong
 use radiation_utils, only:set_radiation_and_gas_temperature_equal
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(2)
 integer :: i,itmp,ierr,iu
 real    :: psep,denszero,enblast,rblast,prblast,gam1
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
 if (do_radiation) call set_units(dist=au,mass=solarm,G=1.d0)

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
    call init_part()
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

    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,hfact,&
                     npart,xyzh,periodic,mask=i_belong)

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
       if ((xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2) < rblast*rblast) then
          vxyzu(iu,i) = prblast/(gam1*denszero)
       else
          vxyzu(iu,i) = 0.
       endif
    enddo
    if (do_radiation) then
       call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad)
       radprop(ikappa,1:npart) = bignumber
    endif
    tmax = 0.1
    dtmax = tmax
    C_cour = 0.15
    C_force = 0.25
!
!--call derivs the first time around
!
    call get_derivs_global()
!
!--now call evolve
!
#ifndef IND_TIMESTEPS
    dt = min(dtcourant,dtforce,dtrad)
#endif
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
