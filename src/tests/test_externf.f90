!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testexternf
!
!  DESCRIPTION:
!  Unit tests of the externalforces module
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: extern_corotate, externalforces, io, part, physcon,
!    testutils, unifdis, units
!+
!--------------------------------------------------------------------------
module testexternf
 implicit none
 public :: test_externf

 private

contains

subroutine test_externf(ntests,npass)
 use io,       only:id,master
 use part,     only:npart,xyzh,hfact,massoftype,igas
 use testutils,only:checkval,checkvalf,checkvalbuf_start,checkvalbuf,checkvalbuf_end
 use externalforces, only:externalforcetype,externalforce,accrete_particles, &
                          was_accreted,iexternalforce_max,initialise_externalforces,&
                          accradius1,update_externalforce,is_velocity_dependent,&
                          externalforce_vdependent,update_vdependent_extforce_leapfrog,&
                          iext_lensethirring,iext_prdrag,iext_einsteinprec,iext_spiral,&
                          iext_neutronstar,iext_staticsine,iext_gwinspiral
 use extern_corotate, only:omega_corotate
 use unifdis,  only:set_unifdis
 use units,    only:set_units
 use physcon,  only:pc,solarm
 integer, intent(inout) :: ntests,npass
 integer                :: i,iextf,nfail1,ierr
 logical                :: dotest1,dotest2,dotest3,accreted
 integer :: nfailed(7),ncheck(7),ierrmax
 real :: psep,fxi,fyi,fzi,dtf,time,pmassi,dhi
 real :: fextxi,fextyi,fextzi,dumx,dumy,dumz,pot1,pot2
 real :: xerrmax,yerrmax,zerrmax,ferrmaxx,ferrmaxy,ferrmaxz
 real :: xi(4),v1(3),fext_iteration(3),fexti(3),vhalfx,vhalfy,vhalfz,dt
 real :: xmini(3),xmaxi(3),poti
 real, parameter :: tolf = 1.5e-3
 real, parameter :: tolfold = 1.e-10

 if (id==master) write(*,"(/,a,/)") '--> TESTING EXTERNAL FORCES MODULE'

 dotest1 = .true.
 dotest2 = .true.
 dotest3 = .true.
!
!--setup 1000 particles, placed randomly
!
 xmini(:) = -100.
 xmaxi(:) = 100.
 hfact      = 1.2
 accradius1 = 100.  ! should be >6 for Lense-Thirring to pass
 psep  = (xmaxi(1) - xmini(1))/10.
 npart = 0
 call set_unifdis('random',id,master,xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3),&
                  psep,hfact,npart,xyzh)
 !dhi   = 0.001*hfact*psep
 dhi   = 1.e-8*psep
 massoftype(igas) = 1./real(npart)
!
!--Test 1: check that external force is the derivative of potential
!
 test1: if (dotest1) then

    time = 0.1 ! just something non-zero
    nfailed(:) = 0
    ncheck(:) = 0
    omega_corotate = 0.5
    do iextf=1,iexternalforce_max
       if (externalforcetype(iextf) /= 'none') then
          select case(iextf)
          case(iext_spiral)
             call set_units(dist=100.*pc,mass=1.d05*solarm,G=1.d0)
          case(iext_lensethirring,iext_prdrag,iext_einsteinprec)
             call set_units(c=1.d0)
          case default
             call set_units(G=1.d0)
          end select
          call initialise_externalforces(iextf,ierr)

          select case(iextf)
          case(iext_neutronstar,iext_gwinspiral)
             !--OK to fail initialisation due to not finding file
             call checkval(ierr,ierr,0,nfailed(1),trim(externalforcetype(iextf))//' external force initialisation')
          case default
             call checkval(ierr,0,0,nfailed(1),trim(externalforcetype(iextf))//' external force initialisation')
          end select

          call update_externalforce(iextf,0.,0.)

          nfailed(2:) = 0
          ncheck(2:) = 0
          xerrmax = 0.
          yerrmax = 0.
          zerrmax = 0.
          ferrmaxx = 0.
          ferrmaxy = 0.
          ferrmaxz = 0.
          do i=1,npart
             xi(:) = xyzh(:,i)
             call externalforce(iextf,xi(1),xi(2),xi(3),xi(4),time, &
                                fxi,fyi,fzi,pot1,dtf)
             !--get derivatives of potential
             call externalforce(iextf,xi(1)+dhi,xi(2),xi(3),xi(4),time, &
                                dumx,dumy,dumz,pot2,dtf)
             fextxi = -(pot2 - pot1)/dhi
             call externalforce(iextf,xi(1),xi(2)+dhi,xi(3),xi(4),time, &
                                dumx,dumy,dumz,pot2,dtf)
             fextyi = -(pot2 - pot1)/dhi
             call externalforce(iextf,xi(1),xi(2),xi(3)+dhi,xi(4),time, &
                                dumx,dumy,dumz,pot2,dtf)
             fextzi = -(pot2 - pot1)/dhi

             call checkvalbuf(fxi,fextxi,tolf,'fextx = -grad phi',nfailed(2),ncheck(2),xerrmax)
             call checkvalbuf(fyi,fextyi,tolf,'fexty = -grad phi',nfailed(3),ncheck(3),yerrmax)
             call checkvalbuf(fzi,fextzi,tolf,'fextz = -grad phi',nfailed(4),ncheck(4),zerrmax)
          enddo
          call checkvalbuf_end('fextx = -grad phi',ncheck(2),nfailed(2),xerrmax,tolf)
          call checkvalbuf_end('fexty = -grad phi',ncheck(3),nfailed(3),yerrmax,tolf)
          call checkvalbuf_end('fextz = -grad phi',ncheck(4),nfailed(4),zerrmax,tolf)
          ntests = ntests + 2
          if (nfailed(1)==0) npass = npass + 1
          if (all(nfailed(2:4)==0)) npass = npass + 1
       endif

       nfailed(:) = 0
    enddo

 endif test1
!
!--Test 2: check that the result of accrete_particles and the was_accreted function
!          agree with each other
!
 test2: if (dotest2) then
    if (id==master) write(*,"(/,a)") '--> testing accrete_particles routine'
    xi(:) = 0.
    pmassi = 1./1000.
    nfail1 = 0
    ncheck(1) = 0
    !call checkvalbuf_start('accreted=was_accreted')
    do iextf=1,iexternalforce_max
       xi(4) = 1.
       call accrete_particles(iextf,xi(1),xi(2),xi(3),xi(4),pmassi,time,accreted)
       call checkvalbuf(was_accreted(iextf,xi(4)),accreted,'accrete/=was_accreted',nfail1,ncheck(1))
    enddo
    ierrmax = 0
    call checkvalbuf_end('accreted=was_accreted for all externf',ncheck(1),nfail1,ierrmax,0)
    ntests = ntests + 1
    if (nfail1==0) npass = npass + 1
 endif test2
!
!--Test 3: check that the update_leapfrog_vdependent routines are correct
!
 test3: if (dotest3) then
    if (id==master) write(*,"(/,a)") '--> testing velocity-dependent external force solvers'
    omega_corotate = 0.5 ! not too high as causes iterations in unit test to fail
    do iextf=1,iexternalforce_max
       if (is_velocity_dependent(iextf)) then
          select case(iextf)
          case(iext_lensethirring,iext_prdrag,iext_einsteinprec)
             call set_units(c=1.d0)
          case default
             call set_units(G=1.d0)
          end select
          call initialise_externalforces(iextf,ierr)
          call checkval(ierr,0,0,nfailed(1),trim(externalforcetype(iextf))//' external force initialisation')
          ntests = ntests + 1
          if (nfailed(1)==0) npass = npass + 1

          call update_externalforce(iextf,0.,0.)
          !
          ! take an arbitrary position and velocity
          !
          xi(1:3) = (/1.14,1.02,1.12/)  ! close to R=1 so black hole spin important
          vhalfx = 0.03   ! keep fairly small, so v-dependent part makes big change
          vhalfy = 0.043
          vhalfz = 0.025
          dt = 0.325     ! something reasonable
          fxi = -0.0789  ! non-zero, but small so that v-dependent
          fyi = 0.036    ! part is dominant component of the force
          fzi = -0.01462
          !
          ! get an explicit evaluation of the external force
          ! and solve v^1 = v^1/2 + dt/2*[f1(x^1) + f1(x^1,v^1)]
          ! by iterating 20 times
          !
          v1 = (/vhalfx + 0.5*dt*fxi,vhalfy + 0.5*dt*fyi,vhalfz + 0.5*dt*fzi/)
          do i=1,30
             call externalforce_vdependent(iextf,xi(1:3),v1,fext_iteration,poti)
             v1(1) = vhalfx + 0.5*dt*(fxi + fext_iteration(1))
             v1(2) = vhalfy + 0.5*dt*(fyi + fext_iteration(2))
             v1(3) = vhalfz + 0.5*dt*(fzi + fext_iteration(3))
             !print*,'fext_iteration = ',fext_iteration
          enddo
          !
          ! call update_leapfrog routine to get analytic solution
          !
          call update_vdependent_extforce_leapfrog(iextf,vhalfx,vhalfy,vhalfz,&
                                 fxi,fyi,fzi,fexti,dt,xi(1),xi(2),xi(3))
          !
          ! check that these agree with each other
          !
          call checkval(fexti(1),fext_iteration(1),1.e-14,nfailed(1),'fx(x1,v1)=v12 + dt/2*[f(x1) + f(x1,v1)]')
          call checkval(fexti(2),fext_iteration(2),1.e-14,nfailed(2),'fy(x1,v1)=v12 + dt/2*[f(x1) + f(x1,v1)]')
          call checkval(fexti(3),fext_iteration(3),1.e-14,nfailed(3),'fz(x1,v1)=v12 + dt/2*[f(x1) + f(x1,v1)]')
          !print*,'fext_iteration = ',fext_iteration
          !print*,'fext_exact     = ',fexti
          ntests = ntests + 1
          if (all(nfailed(1:3)==0)) npass = npass + 1

       endif
    enddo
 endif test3

 if (id==master) write(*,"(/,a)") '<-- EXTERNAL FORCE TESTS COMPLETE'

end subroutine test_externf

end module testexternf
