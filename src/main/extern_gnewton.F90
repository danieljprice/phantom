!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: extern_gnewton
!
!  DESCRIPTION:
! This module contains routines relating to the computation
! of the force in a generalized Newtonian potential (Tejeda
! & Rosswog 2013) which reproduced accurately several features
! of the Schwarzschild space-time.
!
!  REFERENCES: Tejeda E., Rosswog S., 2013, MNRAS, 433, 1930
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: fastmath, io, physcon, units
!+
!--------------------------------------------------------------------------
module extern_gnewton
 implicit none
 public  :: get_gnewton_spatial_force, get_gnewton_vdependent_force
 public  :: update_gnewton_leapfrog
 public  :: get_gnewton_energy

 private

contains

!------------------------------------------------
!+
!  compute the spatial part of the acceleration
!+
!------------------------------------------------
subroutine get_gnewton_spatial_force(xi,yi,zi,mass,fextxi,fextyi,fextzi,phi)
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use units, only: udist, utime
 use physcon, only: c
 real, intent(in)    :: xi,yi,zi,mass
 real, intent(inout) :: fextxi,fextyi,fextzi
 real, intent(out)   :: phi
 real                :: r2,dr,dr3,ccode,rg

 ccode = c/(udist/utime)
 rg = (mass)/(ccode)**2

 r2 = xi*xi + yi*yi + zi*zi

 if (r2 > epsilon(r2)) then
#ifdef FINVSQRT
    dr  = finvsqrt(r2)
#else
    dr = 1./sqrt(r2)
#endif

    dr3 = dr**3
    fextxi = fextxi - mass*xi*dr3*(1.-2.*rg*dr)**2
    fextyi = fextyi - mass*yi*dr3*(1.-2.*rg*dr)**2
    fextzi = fextzi - mass*zi*dr3*(1.-2.*rg*dr)**2
    phi    = -mass*dr

 endif

end subroutine get_gnewton_spatial_force

!-----------------------------------------------------------------------
!+
!  Routine to return velocity-dependent part
!+
!-----------------------------------------------------------------------
subroutine get_gnewton_vdependent_force(xyzi,veli,mass,fexti)
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use physcon, only:c
 use units,   only:utime,udist
 real, intent(in)  :: xyzi(3), veli(3)
 real, intent(in)  :: mass
 real, intent(out) :: fexti(3)
 real              :: r2,dr,dr3,dr5,dr_rel,A,B,ccode,rg
 real              :: xi,yi,zi,vxi,vyi,vzi

 ccode = c/(udist/utime)
 rg = (mass)/(ccode)**2

 xi=xyzi(1)
 yi=xyzi(2)
 zi=xyzi(3)
 vxi=veli(1)
 vyi=veli(2)
 vzi=veli(3)

 r2 = xi*xi + yi*yi + zi*zi

 if (r2 > epsilon(r2)) then
#ifdef FINVSQRT
    dr  = finvsqrt(r2)
#else
    dr = 1./sqrt(r2)
#endif

    dr_rel = 1./(1.-2.*rg*dr)
    dr3 = dr**3
    dr5 = dr**5

    A=xi*vxi+yi*vyi+zi*vzi
    B=(xi*vyi-yi*vxi)**2+(xi*vzi-zi*vxi)**2+(zi*vyi-yi*vzi)**2

    fexti(1) =  2*rg*dr3*dr_rel*vxi*A - 3*rg*dr5*xi*B
    fexti(2) =  2*rg*dr3*dr_rel*vyi*A - 3*rg*dr5*yi*B
    fexti(3) =  2*rg*dr3*dr_rel*vzi*A - 3*rg*dr5*zi*B

 endif

end subroutine get_gnewton_vdependent_force

subroutine update_gnewton_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dt,xi,yi,zi,mass)
 use io,             only:fatal
 real, intent(in)    :: dt,xi,yi,zi, mass
 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(inout) :: fexti(3)
 real                :: fextv(3)
 real                :: v1x, v1y, v1z, v1xold, v1yold, v1zold, vhalf2, erri, dton2
 logical             :: converged
 integer             :: its, itsmax
 integer, parameter  :: maxitsext = 50 ! maximum number of iterations on external force
 character(len=30), parameter :: label = 'update_gnewton_leapfrog'
 real, parameter :: tolv = 1.e-2
 real, parameter :: tolv2 = tolv*tolv
 real,dimension(3) :: pos,vel

 itsmax = maxitsext
 its = 0
 converged = .false.
 dton2 = 0.5*dt

 v1x = vhalfx
 v1y = vhalfy
 v1z = vhalfz
 vhalf2 = vhalfx*vhalfx + vhalfy*vhalfy + vhalfz*vhalfz
 fextv = 0. ! to avoid compiler warning

 iterations : do while (its < itsmax .and. .not.converged)
    its = its + 1
    erri = 0.
    v1xold = v1x
    v1yold = v1y
    v1zold = v1z
    pos = (/xi,yi,zi/)
    vel = (/v1x,v1y,v1z/)
    call get_gnewton_vdependent_force(pos,vel,mass,fextv)
!    xi = pos(1)
!    yi = pos(2)
!    zi = pos(3)
    v1x = vel(1)
    v1y = vel(2)
    v1z = vel(3)

    v1x = vhalfx + dton2*(fxi + fextv(1))
    v1y = vhalfy + dton2*(fyi + fextv(2))
    v1z = vhalfz + dton2*(fzi + fextv(3))

    erri = (v1x - v1xold)**2 + (v1y - v1yold)**2 + (v1z - v1zold)**2
    erri = erri / vhalf2
    converged = (erri < tolv2)

 enddo iterations

 if (its >= maxitsext) call fatal(label,'VELOCITY ITERATIONS ON EXTERNAL FORCE NOT CONVERGED!!')

 fexti(1) = fextv(1)
 fexti(2) = fextv(2)
 fexti(3) = fextv(3)

 fxi = fxi + fexti(1)
 fyi = fyi + fexti(2)
 fzi = fzi + fexti(3)

end subroutine update_gnewton_leapfrog


!-----------------------------------------------------------------------
!+
!  Routine to return energy and angular momentum
!+
!-----------------------------------------------------------------------
subroutine get_gnewton_energy(xyzi,veli,mass,energy,angmomx,angmomy,angmomz)
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use physcon, only:c
 use units,   only:utime,udist
 implicit none
 real, dimension(3), intent(in)    :: xyzi, veli
 real,               intent(in)    :: mass
 real, intent(inout)   :: energy, angmomx, angmomy, angmomz
 real                              :: r2,rg,ccode
 real                              :: xi,yi,zi,vxi,vyi,vzi
 real                              :: dr,dr2,dr_rel,dr_rel2
 real                              :: A,B,A2

 xi=xyzi(1)
 yi=xyzi(2)
 zi=xyzi(3)
 vxi=veli(1)
 vyi=veli(2)
 vzi=veli(3)

 ccode = c/(udist/utime)
 rg = (mass)/(ccode)**2

 r2 = xi*xi + yi*yi + zi*zi

#ifdef FINVSQRT
 dr = finvsqrt(r2)
#else
 dr = 1./sqrt(r2)
#endif

 dr_rel = 1./(1.-2.*rg*dr)
 dr_rel2 = dr_rel**2
 dr2 = dr**2

 A=xi*vxi+yi*vyi+zi*vzi
 B=(xi*vyi-yi*vxi)**2+(xi*vzi-zi*vxi)**2+(zi*vyi-yi*vzi)**2
 A2=A**2

 energy = 0.5*(A2*dr_rel2*dr2 + B*dr2*dr_rel) - mass*dr ! energy
 angmomx = (yi*vzi-zi*vyi)*dr_rel !
 angmomy = (zi*vxi-xi*vzi)*dr_rel ! angular momentum
 angmomz = (xi*vyi-yi*vxi)*dr_rel !

end subroutine get_gnewton_energy

end module extern_gnewton
