!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testgnewton
!
!  DESCRIPTION:
!  unit tests of the gnewton external force module
!
!  REFERENCES: Tejeda E., Rosswog S., 2013, MNRAS, 433, 1930
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: extern_gnewton, externalforces, io, options, part,
!    physcon, testutils, timestep, units
!+
!--------------------------------------------------------------------------
module testgnewton
 implicit none
 public :: test_gnewton

 private

contains

subroutine test_gnewton(ntests,npass)
 use io,              only:id,master
 use part,            only:xyzh,vxyzu
 use timestep,        only:C_force
 use testutils,       only:checkval
 use physcon,         only:solarm,solarr,c,pi
 use options,         only:iexternalforce
 use units,           only:set_units,udist,utime
 use extern_gnewton,  only:get_gnewton_energy
 use externalforces,  only:externalforce,externalforce_vdependent,mass1,iext_gnewton
 integer, intent(inout) :: ntests,npass
 logical                :: test_relatorbit
 real    :: t,dtnew
 integer :: nfailed(3)
 real    :: ccode
 real    :: e,ra,rp,rela,relp,rg,va,period,a
 real    :: fextx,fexty,fextz,phi,dtf,dt
 real    :: energyin,angmomin,angmomxin,angmomyin,angmomzin
 real    :: energy,angmom,angmomx,angmomy,angmomz
 real    :: errenergy,errangmom
 real    :: ck,k,pign,precang,rb,vdotr

 if (id==master) write(*,"(/,a,/)") '--> TESTING GNEWTON MODULE'

 test_relatorbit = .true.
 !
 !--Test : eccentric orbit of a test particle
 !
 testrelatorbit: if (test_relatorbit) then
    if (id==master) write(*,"(/,a)") '--> testing relativistic orbit'
    !
    !--setup a test particle on an eccentric orbit
    !
    call set_units(dist=solarr,mass=solarm,G=1.d0)

    mass1 = 1e6 ! mass of the central object
    ccode = c/(udist/utime)
    rg = (mass1)/(ccode)**2
    ra=180.
    rp=20.
    e=(ra-rp)/(ra+rp)
    a=0.5*(ra+rp)
    rela=1-2*rg/ra
    relp=1-2*rg/rp
    va=sqrt( (2*mass1*(1/ra-1/rp)*rela) / ( 1-(relp/rela)*(ra/rp)**2 ) ) ! velocity at apocenter

    xyzh(1,1) = ra
    xyzh(2,1) = 0.
    xyzh(3,1) = 0.
    xyzh(4,1) = 1.e-2

    vxyzu(1,1) = 0.
    vxyzu(2,1) = va
    vxyzu(3,1) = 0.

    iexternalforce = iext_gnewton ! choose the gnewton force

    call get_gnewton_energy(xyzh(1:3,1),vxyzu(1:3,1),mass1,energyin,angmomxin,angmomyin,angmomzin)
    angmomin = sqrt(angmomxin**2 + angmomyin**2 + angmomzin**2)

    if (id==master) then
       print "(/,2x,a)",'---------- orbital parameters ---------- '
       print "(6(2x,a,f15.3,/),2x,a,5x,es10.3)", &
             'apocenter               :',ra, &
             'pericenter              :',rp, &
             'semi-major axis         :',a, &
             'eccentricity            :',e, &
             'initial velocity        :',va, &
             'gravitational radius    :',rg, &
             'central mass            :',mass1
       print "(2x,40('-'),/)"
    endif
    !
    !--evolve this
    !
    C_force = 1e-4
    if (id==master) print*,'integrating the orbit for one azimuthal period... ( C_force = ',C_force,')'

    t = 0.
    call externalforce(iexternalforce,xyzh(1,1),xyzh(2,1),xyzh(3,1),xyzh(4,1), &
                             t,fextx,fexty,fextz,phi,dtf)
    dt = C_force*dtf ! initial timestep
    period = 2.*pi*sqrt(a**3/mass1)

    errenergy = 0.
    errangmom = 0.
    vdotr = 0.

    evolve : do while (t < 0.75*period .or. vdotr > 0.)
       call get_gnewton_energy(xyzh(1:3,1),vxyzu(1:3,1),mass1,energy,angmomx,angmomy,angmomz)
       vdotr = xyzh(1,1)*vxyzu(1,1) + xyzh(2,1)*vxyzu(2,1) + xyzh(3,1)*vxyzu(3,1)
       angmom = sqrt(angmomx**2 + angmomy**2 + angmomz**2)
       errenergy = max(errenergy,abs(energy - energyin))
       errangmom = max(errangmom,abs(angmom - angmomin))
       call step_lf(t,dt,dtnew)
       t = t + dt
       dt = dtnew
    enddo evolve

    call checkval(energyin+errenergy,energyin,1.e-3,nfailed(1),'energy')
    call checkval(angmomin+errangmom,angmomin,1.e-3,nfailed(2),'angular momentum')
    !
    !--compute the precession angle (eq. (3.23) of Tejeda and Rosswog (2013))
    !
    precang=atan2(xyzh(2,1),xyzh(1,1))
    if (precang < 0.) precang=2.*pi+precang ! the precession angle must be between 0 and 2 pi
    rb=-mass1/energyin-ra-rp
    k=sqrt((rb*(ra-rp))/(rp*(ra-rb)))
    ck=rf(0.,1.-k**2,1.)
    pign=(2.*angmomin*ck)/sqrt(-2*energyin*rp*(ra-rb))

    call checkval(precang,2.*(pign-real(pi)),1.e-5,nfailed(3),'precession angle')

    ntests = ntests + 1
    if (all(nfailed(1:3)==0)) npass = npass + 1

 endif testrelatorbit

 if (id==master) write(*,"(/,a)") '<-- GNEWTON TEST COMPLETE'

end subroutine test_gnewton

!-----------------------------------------------------------------------
!+
!  Leapfrog algorithm for a test particle
!+
!-----------------------------------------------------------------------
subroutine step_lf(t,dt,dtnew)
 use externalforces, only:externalforce,update_vdependent_extforce_leapfrog,externalforce_vdependent
 use timestep,       only:C_force
 use part,           only:xyzh,vxyzu
 use options,        only:iexternalforce
 real,    intent(in)    :: t,dt
 real,    intent(out)   :: dtnew
 real    :: dtf,hdt,phi
 real    :: fextx,fexty,fextz,fextv(3),fx,fy,fz,vxhalf,vyhalf,vzhalf

 hdt = 0.5*dt

 call externalforce(iexternalforce,xyzh(1,1),xyzh(2,1),xyzh(3,1),xyzh(4,1), &
                          t,fextx,fexty,fextz,phi)
 call externalforce_vdependent(iexternalforce,xyzh(1:3,1),vxyzu(1:3,1),fextv,phi)

 fx = fextx + fextv(1)
 fy = fexty + fextv(2)
 fz = fextz + fextv(3)

 vxhalf = vxyzu(1,1) + hdt*fx
 vyhalf = vxyzu(2,1) + hdt*fy
 vzhalf = vxyzu(3,1) + hdt*fz

 xyzh(1,1) = xyzh(1,1) + dt*vxhalf
 xyzh(2,1) = xyzh(2,1) + dt*vyhalf
 xyzh(3,1) = xyzh(3,1) + dt*vzhalf

 call externalforce(iexternalforce,xyzh(1,1),xyzh(2,1),xyzh(3,1),xyzh(4,1), &
                          t,fextx,fexty,fextz,phi)

 fx = fextx
 fy = fexty
 fz = fextz

 call update_vdependent_extforce_leapfrog(iexternalforce,&
        vxhalf,vyhalf,vzhalf, &
        fx,fy,fz,fextv,dt,xyzh(1,1),xyzh(2,1),xyzh(3,1))

 vxyzu(1,1) = vxhalf + hdt*fx
 vxyzu(2,1) = vyhalf + hdt*fy
 vxyzu(3,1) = vzhalf + hdt*fz

 call externalforce(iexternalforce,xyzh(1,1),xyzh(2,1),xyzh(3,1),xyzh(4,1), &
                          t,fextx,fexty,fextz,phi,dtf)
 dtnew = C_force*dtf ! new timestep

end subroutine step_lf

!-----------------------------------------------------------------------
!+
!  Computes the elliptic integral
!+
!-----------------------------------------------------------------------
real function rf(x,y,z)
 real, intent(in) :: x,y,z
 real, parameter :: errtol=.08,THIRD=1./3.,C1=1./24.,C2=.1,C3=3./44.,C4=1./14.
 real :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
 real :: errmax

 if(min(x,y,z) < 0..or.min(x+y,x+z,y+z) < tiny(x) .or.max(x,y,z) > huge(x)) print*,"invalid arguments in rf"
 xt=x
 yt=y
 zt=z
 errmax = huge(errmax)
 do while (errmax > errtol)
    sqrtx=sqrt(xt)
    sqrty=sqrt(yt)
    sqrtz=sqrt(zt)
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
    xt=.25*(xt+alamb)
    yt=.25*(yt+alamb)
    zt=.25*(zt+alamb)
    ave=THIRD*(xt+yt+zt)
    delx=(ave-xt)/ave
    dely=(ave-yt)/ave
    delz=(ave-zt)/ave
    errmax = max(abs(delx),abs(dely),abs(delz))
 enddo
 e2=delx*dely-delz**2
 e3=delx*dely*delz
 rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
 return

end function rf

end module testgnewton
