!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: dustywaves
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
!  DEPENDENCIES: cubic
!+
!--------------------------------------------------------------------------

! ----------------------------------------------------------------------
! compute exact solution for a linear wave
! plots a sine function with a given amplitude, period and wavelength
! ----------------------------------------------------------------------
module dustywaves
 implicit none

contains

subroutine exact_dustywave(time,ampl,cs,Kdragin,lambda,x0,ymean_gas,ymean_dust,xplot,vgaso,vdusto,rhogaso,rhodusto,ierr)
 use cubic,   only:cubicsolve_complex
 !use plotlib, only:plot_line,plot_sls
 integer :: i
 real, parameter :: pi = 3.1415926536
 real,    intent(in)  :: time, ampl, cs, Kdragin, lambda, x0, ymean_gas, ymean_dust
 real,    intent(in)  :: xplot(:)
 real,    intent(out) :: vgaso(size(xplot)),vdusto(size(xplot)),rhogaso(size(xplot)),rhodusto(size(xplot))
 integer, intent(out) :: ierr
 real :: rhodeq,rhogeq,rhodsol,rhogsol,vdeq,vgeq,vgsol,vdsol
 real :: aa,bb,cc,w1r,w2r,w3r,w1i,w2i,w3i
 real :: k,xk,arg1,arg2,arg3,vgas,vdust,rhogas,rhodust
 real :: vd1r,vd1i,vd2r,vd2i,vd3r,vd3i
 real :: vg1r,vg1i,vg2r,vg2i,vg3r,vg3i
 real :: rhod1r,rhod1i,rhod2r,rhod2i,rhod3r,rhod3i
 real :: rhog1r,rhog1i,rhog2r,rhog2i,rhog3r,rhog3i
 real :: tgas1,tdust1,Kdrag
 complex :: xc(3)

 Kdrag = Kdragin

 print*,'plotting two-fluid gas/dust linear wave solution ... '
 print*,' lambda = ',lambda,' ampl = ',ampl,' cs = ',cs,' Kdrag = ',Kdrag
!
! check for errors
!
 ierr = 0
 if (ampl < 0.) then
    print*,'error: amplitude < 0 on input'
    ierr = 1
    return
 endif
 if (lambda <= 0.) then
    print*,'error: lambda <= 0 on input'
    ierr = 2
    return
 endif
 if (cs <= 0) then
    print*,'error: sound speed <= 0 on input'
    ierr = 3
    return
 endif
 if (ymean_gas <= 0) then
    print*,'error: mean density <= 0 on input'
    ierr = 4
    return
 endif
 if (ymean_dust <= 0) then
    print*,'error: mean density <= 0 on input'
    ierr = 5
    return
 endif
 if (Kdrag < 0) then
    print*,'error: drag coefficient < 0 on input'
    ierr = 6
    return
 elseif (abs(Kdrag) < 1.e-8) then
    print*,' WARNING: Kdrag = 0 on input; using tiny to avoid divergence '
    Kdrag = 0.
 endif

 rhodeq  = ymean_dust ! initial dust density
 rhogeq  = ymean_gas ! initial gas density
 print*,' rho(dust),0 = ',ymean_dust,' rho(gas),0 = ',ymean_gas
 rhodsol = ampl*ymean_dust  ! amplitude of dust density perturbation
 rhogsol = ampl*ymean_gas ! amplitude of gas density perturbation
 vdeq    = 0.
 vgeq    = 0.
 vgsol   = ampl  ! amplitude of gas velocity perturbation
 vdsol   = ampl  ! amplitude of dust velocity perturbation
 !cs      = 1.0
 !Kdrag   = 1.0
 k       = 2.*pi/lambda ! wavenumber

 vd1r = 0.
 vd1i = 0.
 vd2r = 0.
 vd2i = 0.
 vd3r = 0.
 vd3i = 0.
 rhod1r = 0.
 rhod1i = 0.
 rhod2r = 0.
 rhod2i = 0.
 rhod3r = 0.
 rhod3i = 0.
 !
 !--solve cubic to get the 3 solutions for omega
 !  (these each have both real and imaginary components,
 !   labelled w1r, w1i etc.)
 !
 tdust1 = Kdrag/rhodeq
 tgas1  = Kdrag/rhogeq
 aa = (tdust1 + tgas1)
 bb = k**2*cs**2
 cc = bb*tdust1

 call cubicsolve_complex(aa,bb,cc,xc)
 !--get solutions for (w = iy instead of y)
 xc = xc*cmplx(0.,1.)
 print*,' roots are ',xc

 w1r = real(xc(1))
 w2r = real(xc(2))
 w3r = real(xc(3))

 w1i = aimag(xc(1))
 w2i = aimag(xc(2))
 w3i = aimag(xc(3))

!-------------------------------
! G A S  V E L O C I T I E S
!-------------------------------
 vg3r =(k*Kdrag*vdsol*w3r**2*w2r + k*Kdrag*vdsol*w3r**2*w1r - k*Kdrag*vdsol*w3r*w2r*w1r - k*Kdrag*vdsol*w3r*w3i**2 +&
        k*Kdrag*vdsol*w2i*w1i*w3r - k*Kdrag*vdsol*w2r*w1i*w3i + k*Kdrag*vdsol*w2r*w3i**2 - k*Kdrag*vdsol*w2i*w3i*w1r +&
        k*Kdrag*vdsol*w3i**2*w1r - k*Kdrag*vgsol*w3r**2*w2r - k*Kdrag*vgsol*w3r**2*w1r + k*Kdrag*vgsol*w3r*w2r*w1r +&
        k*Kdrag*vgsol*w3r*w3i**2 - k*Kdrag*vgsol*w2i*w1i*w3r + k*Kdrag*vgsol*w2r*w1i*w3i - k*Kdrag*vgsol*w2r*w3i**2 +&
        k*Kdrag*vgsol*w2i*w3i*w1r - k*Kdrag*vgsol*w3i**2*w1r - rhogsol*w3r*w3i**2*w2r*w1i - rhogsol*w3r*w3i**2*w2i*w1r&
        + k*rhogeq*vgsol*w3r**3*w1i + k*rhogeq*vgsol*w3r**3*w2i - k*rhogeq*vgsol*w3r**2*w2r*w3i -&
        k*rhogeq*vgsol*w3r**2*w3i*w1r - k*rhogeq*vgsol*w3r*w2r**2*w1i - k*rhogeq*vgsol*w3r*w1i*w2i**2 +&
        k*rhogeq*vgsol*w3r*w3i**2*w2i - k*rhogeq*vgsol*w2r*w3i**3 - rhogsol*w3r**3*w1i*w2r - rhogsol*w3r**3*w2i*w1r +&
        rhogsol*w3r**2*w2r**2*w1i + rhogsol*w3r**2*w1i*w2i**2 + rhogsol*w3r**2*w2i*w1r**2 + rhogsol*w3r**2*w2i*w1i**2 -&
        rhogsol*w2r**2*w1i**2*w3i + rhogsol*w2r**2*w1i*w3i**2 - rhogsol*w3i*w1r**2*w2r**2 + rhogsol*w3i**3*w1r*w2r +&
        rhogsol*w3i**2*w2i*w1r**2 + rhogsol*w2i*w1i**2*w3i**2 - rhogsol*w1i**2*w2i**2*w3i - rhogsol*w2i**2*w1r**2*w3i +&
        rhogsol*w2i**2*w3i**2*w1i - rhogsol*w2i*w3i**3*w1i - rhogsol*k**2*cs**2*w3r**2*w2i -&
        rhogsol*k**2*cs**2*w3r**2*w1i + rhogsol*k**2*cs**2*w3r*w1i*w2r + rhogsol*k**2*cs**2*w3r*w2i*w1r -&
        rhogsol*k**2*cs**2*w3i*w2r*w1r + rhogsol*k**2*cs**2*w2i*w3i*w1i - rhogsol*k**2*cs**2*w3i**2*w1i -&
        rhogsol*k**2*cs**2*w3i**2*w2i + rhogsol*k**2*cs**2*w3r**2*w3i - k*rhogeq*vgsol*w3r*w2i*w1r**2 +&
        k*rhogeq*vgsol*w3r*w3i**2*w1i - k*rhogeq*vgsol*w3r*w2i*w1i**2 + k*rhogeq*vgsol*w3i*w1r*w2r**2 +&
        k*rhogeq*vgsol*w3i*w2r*w1r**2 + k*rhogeq*vgsol*w3i*w2r*w1i**2 + k*rhogeq*vgsol*w3i*w1r*w2i**2 -&
        k*rhogeq*vgsol*w3i**3*w1r - k*Kdrag*vdsol*w3r**3 + k*Kdrag*vgsol*w3r**3 + rhogsol*k**2*cs**2*w3i**3 +&
        rhogsol*w3r**2*w3i*w2r*w1r - rhogsol*w3r**2*w2i*w3i*w1i)/rhogeq/k/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 -&
        2*w2i*w3i + w3r**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)

 vg3i = - 1/rhogeq/k*( - w3r*w3i**2*cs**2*k**2*rhogsol - w3r**3*rhogsol*k**2*cs**2 + w3r**3*w2i*rhogsol*w1i +&
        w3r**2*w3i*w2i*rhogeq*k*vgsol - w3r**2*w2i*w3i*w1r*rhogsol - w2r*w3r**2*w3i*rhogsol*w1i -&
        w2r*w1i*w3i**3*rhogsol - w2r*w3r**3*w1r*rhogsol - w3r**2*w1i*k*Kdrag*vgsol + w3r**2*w2i*k*Kdrag*vdsol +&
        w3r*w3i**2*k*rhogeq*w1r*vgsol + w3r*w1i*w3i**2*rhogsol*w2i - w3r**2*w2i*k*Kdrag*vgsol +&
        w3r**2*w1i*k*Kdrag*vdsol - w3r**2*w1i**2*k*rhogeq*vgsol + w3r**3*k*rhogeq*w1r*vgsol +&
        w3r**2*w3i*w1i*rhogeq*k*vgsol - w3r**2*w1r**2*rhogeq*k*vgsol + w3r**2*w1r*cs**2*k**2*rhogsol +&
        w3r**2*w3i*k*Kdrag*vgsol - w3r**2*w3i*k*Kdrag*vdsol + w2r*w3r*rhogeq*k*vgsol*w3i**2 +&
        w2r*w3r**2*cs**2*k**2*rhogsol + w2r*w3r**3*vgsol*k*rhogeq + w2r*w3i**2*cs**2*k**2*rhogsol -&
        w3i**2*k*rhogeq*w1r**2*vgsol + w3i**2*cs**2*k**2*rhogsol*w1r - 2*w3r**2*k*w2i*rhogeq*w1i*vgsol -&
        w1i**2*rhogeq*k*vgsol*w3i**2 - 2*w3i**2*k*w2i*rhogeq*w1i*vgsol - w3r**2*w2i**2*rhogeq*k*vgsol +&
        w2i*w1i**2*w3i*rhogeq*k*vgsol - w2i*w3i*w1i*k*Kdrag*vdsol + w2i*w3i*w1i*k*Kdrag*vgsol -&
        w2i*w3i*w1r*cs**2*k**2*rhogsol + w2i*rhogeq*k*vgsol*w3i*w1r**2 + w2r*w3i*w1r*k*Kdrag*vdsol -&
        w2i*w3i**2*k*Kdrag*vgsol - w2r*w3i*w1r*k*Kdrag*vgsol + w2i*w3i**2*k*Kdrag*vdsol -&
        2*w2r*w3i**2*k*rhogeq*w1r*vgsol + w2i*w3i**3*rhogeq*k*vgsol + w3r*w2i**2*w1r*rhogeq*k*vgsol +&
        w2r**2*w3i**2*w1r*rhogsol - w1i*w3i**2*k*Kdrag*vgsol + w1i*rhogeq*k*vgsol*w3i**3 + w1i*w3i**2*k*Kdrag*vdsol -&
        w2r*w3r*w3i**2*w1r*rhogsol - w3i**3*k*Kdrag*vdsol + w3i**3*k*Kdrag*vgsol - w2r*w1i*w3i*cs**2*k**2*rhogsol -&
        w1r*w2i*w3i**3*rhogsol + w3i*w2r**2*w1i*rhogeq*k*vgsol + w1i**2*w3i**2*w2r*rhogsol -&
        w3i**2*w2i**2*rhogeq*k*vgsol + w2i**2*w3i*w1i*rhogeq*k*vgsol + w3r*w2i*w1i*cs**2*k**2*rhogsol -&
        w3r*w1i**2*w2r**2*rhogsol - w3r*w2i**2*w1i**2*rhogsol + w3r**2*w1i**2*w2r*rhogsol + w3r**2*w2i**2*w1r*rhogsol -&
        w2r**2*w3i**2*k*rhogeq*vgsol + w3r**2*w2r*w1r**2*rhogsol - w3r**2*w2r**2*rhogeq*k*vgsol +&
        w3r**2*w2r**2*w1r*rhogsol - w3r*w1r**2*w2i**2*rhogsol - w3r*w2r**2*w1r**2*rhogsol - w3r*w2i*w1r*k*Kdrag*vdsol +&
        w3r*w2i*w1r*k*Kdrag*vgsol - w3r*w2r*w1r*cs**2*k**2*rhogsol - 2*w3r**2*w2r*k*rhogeq*w1r*vgsol +&
        w3r*w2r**2*k*rhogeq*w1r*vgsol + w3r*w2r*w1i*k*Kdrag*vgsol + w3r*w2r*w1i**2*k*rhogeq*vgsol -&
        w3r*w2r*w1i*k*Kdrag*vdsol + w3r*w2r*w1r**2*rhogeq*k*vgsol + w2i**2*w3i**2*w1r*rhogsol +&
        w2r*w3i**2*w1r**2*rhogsol)/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i + w3r**2)/(w1i**2 - 2*w3i*w1i +&
        w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)

 vg2r = - (w2r**2*rhogsol*w2i*w3i*w1i + w2r**2*rhogsol*k**2*cs**2*w1i + k*Kdrag*vdsol*w3r*w2r*w1r +&
        k*Kdrag*vdsol*w2i*w1i*w3r - k*Kdrag*vdsol*w2r*w1i*w3i + k*Kdrag*vdsol*w2i*w3i*w1r - k*Kdrag*vgsol*w3r*w2r*w1r -&
        k*Kdrag*vgsol*w2i*w1i*w3r + k*Kdrag*vgsol*w2r*w1i*w3i - k*Kdrag*vgsol*w2i*w3i*w1r - rhogsol*w3r**2*w2r**2*w1i -&
        rhogsol*w3r**2*w1i*w2i**2 + rhogsol*w3r**2*w2i*w1r**2 + rhogsol*w3r**2*w2i*w1i**2 - rhogsol*w2r**2*w1i**2*w3i -&
        rhogsol*w2r**2*w1i*w3i**2 - rhogsol*w3i*w1r**2*w2r**2 + rhogsol*w3i**2*w2i*w1r**2 + rhogsol*w2i*w1i**2*w3i**2 -&
        rhogsol*w1i**2*w2i**2*w3i - rhogsol*w2i**2*w1r**2*w3i - rhogsol*w2i**2*w3i**2*w1i + w2r**3*rhogsol*w3i*w1r -&
        rhogsol*k**2*cs**2*w3r*w1i*w2r + rhogsol*k**2*cs**2*w3r*w2i*w1r - rhogsol*k**2*cs**2*w3i*w2r*w1r -&
        rhogsol*k**2*cs**2*w2i*w3i*w1i - k*rhogeq*vgsol*w3r*w2i*w1r**2 - k*rhogeq*vgsol*w3r*w2i*w1i**2 +&
        k*rhogeq*vgsol*w3i*w2r*w1r**2 + k*rhogeq*vgsol*w3i*w2r*w1i**2 + w2r**2*rhogsol*k**2*cs**2*w3i +&
        w2i**3*rhogsol*w3i*w1i - w2r**3*k*Kdrag*vgsol + w2r**3*k*Kdrag*vdsol - w2r**3*vgsol*k*rhogeq*w1i +&
        w2i*vgsol*k*rhogeq*w1r*w2r**2 + w2i**3*vgsol*k*rhogeq*w1r - w2i**2*w2r*vgsol*k*rhogeq*w1i -&
        w3r*w2i**3*rhogsol*w1r + w3r*w2r**3*rhogsol*w1i - w3r*w2r**2*rhogsol*w2i*w1r - w3r*w2i**2*k*Kdrag*vdsol +&
        w3r*w2i**2*rhogsol*w1i*w2r + w3r*w2i**2*k*Kdrag*vgsol - w3r**2*vgsol*k*rhogeq*w2i*w1r +&
        w3r**2*w2r*k*rhogeq*vgsol*w1i + w3r*w2r**2*k*Kdrag*vgsol - w3r*w2r**2*k*Kdrag*vdsol +&
        w3r*w2r**2*k*rhogeq*vgsol*w2i + w3r*w2i**3*k*rhogeq*vgsol - w2r**3*k*rhogeq*vgsol*w3i -&
        w2i**3*rhogsol*k**2*cs**2 - vgsol*k*rhogeq*w2i*w1r*w3i**2 + w2i**2*k*Kdrag*vgsol*w1r - w2i**2*k*Kdrag*vgsol*w2r&
        - w2i**2*k*Kdrag*vdsol*w1r + w2i**2*rhogsol*k**2*cs**2*w3i - w2i**2*k*rhogeq*vgsol*w2r*w3i +&
        w2i**2*rhogsol*k**2*cs**2*w1i + w2i**2*rhogsol*w3i*w2r*w1r - w2r**2*k*Kdrag*vdsol*w1r +&
        w2r**2*k*Kdrag*vgsol*w1r + w2r*vgsol*k*rhogeq*w1i*w3i**2 + w2i**2*k*Kdrag*vdsol*w2r -&
        w2r**2*rhogsol*k**2*cs**2*w2i)/rhogeq/k/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i + w3r**2)/(w2r**2 +&
        w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)

 vg2i =(w1r**2*w2i**2*k*rhogeq*vgsol - w2r**2*w2i*k*Kdrag*vgsol + w2i**2*w1i**2*k*rhogeq*vgsol -&
        w1i*w2i**3*rhogeq*k*vgsol - w2i**3*w3i*rhogeq*k*vgsol - w3r**2*k*w2i*rhogeq*w1i*vgsol -&
        w3i**2*k*w2i*rhogeq*w1i*vgsol + w3r*w2i*w1i*w2r**2*rhogsol + w3r**2*w2i**2*rhogeq*k*vgsol -&
        w2i**2*w1r*cs**2*k**2*rhogsol - w2i*w1i**2*w3i*rhogeq*k*vgsol + w2i*w3i*w1i*k*Kdrag*vdsol -&
        w2i*w3i*w1i*k*Kdrag*vgsol - w2r**2*w3i*w2i*rhogeq*k*vgsol + w2i*w3i*w1r*cs**2*k**2*rhogsol -&
        w2i*rhogeq*k*vgsol*w3i*w1r**2 + w2r**3*cs**2*k**2*rhogsol - w2r**2*w1i*w2i*rhogeq*k*vgsol +&
        w2r**2*w2i*w3i*w1r*rhogsol + w2r*w3i*w1r*k*Kdrag*vdsol - w2i**2*w1i*k*Kdrag*vdsol - w2r*w3i*w1r*k*Kdrag*vgsol -&
        w2r*w3i**2*k*rhogeq*w1r*vgsol - w2r*w1r*w2i**2*rhogeq*k*vgsol + 2*w3r*w2i**2*w1r*rhogeq*k*vgsol +&
        w2r**2*w1i*k*Kdrag*vgsol - w2r**2*w3i**2*w1r*rhogsol + w3r*w2i**3*w1i*rhogsol + w2i**2*w1i*k*Kdrag*vgsol -&
        w2r*w1i*w3i*cs**2*k**2*rhogsol + 2*w3i*w2r**2*w1i*rhogeq*k*vgsol + w1i**2*w3i**2*w2r*rhogsol +&
        w3i**2*w2i**2*rhogeq*k*vgsol + 2*w2i**2*w3i*w1i*rhogeq*k*vgsol + w3r*w2i*w1i*cs**2*k**2*rhogsol -&
        w3r*w1i**2*w2r**2*rhogsol - w3r*w2i**2*w1i**2*rhogsol + w3r**2*w1i**2*w2r*rhogsol - w3r**2*w2i**2*w1r*rhogsol +&
        w2r**2*w2i*k*Kdrag*vdsol - w2r**3*k*rhogeq*w1r*vgsol + w2i**2*w3i*k*Kdrag*vgsol - w2i**2*w3i*k*Kdrag*vdsol +&
        w2r*w2i**2*cs**2*k**2*rhogsol - w2r*w2i**2*w3i*rhogsol*w1i + w2r**2*w3i**2*k*rhogeq*vgsol +&
        w2r**2*w1r**2*rhogeq*k*vgsol - w2r**2*w1r*cs**2*k**2*rhogsol + w3r**2*w2r*w1r**2*rhogsol +&
        w2r**2*w3i*k*Kdrag*vgsol - w2r**2*w3i*k*Kdrag*vdsol - w2r**2*w1i*k*Kdrag*vdsol + w2r**2*w1i**2*k*rhogeq*vgsol -&
        w3r*w2r**2*cs**2*k**2*rhogsol + w3r*w2r*w2i**2*w1r*rhogsol - w3r*w2i**2*cs**2*k**2*rhogsol -&
        w3r*w2r**3*rhogeq*k*vgsol + w3r**2*w2r**2*rhogeq*k*vgsol - w3r**2*w2r**2*w1r*rhogsol - w2r**3*w3i*rhogsol*w1i -&
        w3r*w1r**2*w2i**2*rhogsol - w3r*w2r**2*w1r**2*rhogsol - w3r*w2i*w1r*k*Kdrag*vdsol + w3r*w2i*w1r*k*Kdrag*vgsol +&
        w3r*w2r**3*w1r*rhogsol + w3r*w2r*w1r*cs**2*k**2*rhogsol - w3r**2*w2r*k*rhogeq*w1r*vgsol +&
        2*w3r*w2r**2*k*rhogeq*w1r*vgsol - w3r*w2r*w2i**2*rhogeq*k*vgsol - w3r*w2r*w1i*k*Kdrag*vgsol -&
        w3r*w2r*w1i**2*k*rhogeq*vgsol + w3r*w2r*w1i*k*Kdrag*vdsol - w3r*w2r*w1r**2*rhogeq*k*vgsol +&
        w2i**3*w3i*w1r*rhogsol - w2i**3*k*Kdrag*vgsol - w2i**2*w3i**2*w1r*rhogsol + w2i**3*k*Kdrag*vdsol +&
        w2r*w3i**2*w1r**2*rhogsol)/rhogeq/k/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i + w3r**2)/(w2r**2 +&
        w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)

 vg1r =( - rhogsol*w1i**3*w3i*w2i + rhogsol*k**2*cs**2*w1i**3 - rhogsol*w3i*w2r*w1r**3 - k*Kdrag*vgsol*w2r*w1r**2 +&
        k*Kdrag*vdsol*w2r*w1r**2 - k*Kdrag*vdsol*w1r*w1i**2 - k*Kdrag*vgsol*w2r*w1i**2 + k*Kdrag*vdsol*w2r*w1i**2 +&
        k*Kdrag*vgsol*w1r*w1i**2 - k*Kdrag*vdsol*w1r**3 + k*Kdrag*vgsol*w1r**3 - w3r*rhogsol*w2i*w1r**3 -&
        w3r*Kdrag*w1i**2*k*vgsol + w3r*k*Kdrag*vdsol*w1r**2 + w3r*k*Kdrag*vdsol*w1i**2 - w3r*Kdrag*w1r**2*k*vgsol -&
        k*Kdrag*vdsol*w3r*w2r*w1r - k*Kdrag*vdsol*w2i*w1i*w3r - k*Kdrag*vdsol*w2r*w1i*w3i + k*Kdrag*vdsol*w2i*w3i*w1r +&
        k*Kdrag*vgsol*w3r*w2r*w1r + k*Kdrag*vgsol*w2i*w1i*w3r + k*Kdrag*vgsol*w2r*w1i*w3i - k*Kdrag*vgsol*w2i*w3i*w1r +&
        k*rhogeq*vgsol*w3r*w2r**2*w1i + k*rhogeq*vgsol*w3r*w1i*w2i**2 - rhogsol*w3r**2*w2r**2*w1i -&
        rhogsol*w3r**2*w1i*w2i**2 + rhogsol*w3r**2*w2i*w1r**2 + rhogsol*w3r**2*w2i*w1i**2 + rhogsol*w2r**2*w1i**2*w3i -&
        rhogsol*w2r**2*w1i*w3i**2 + rhogsol*w3i*w1r**2*w2r**2 + rhogsol*w3i**2*w2i*w1r**2 + rhogsol*w2i*w1i**2*w3i**2 +&
        rhogsol*w1i**2*w2i**2*w3i + rhogsol*w2i**2*w1r**2*w3i - rhogsol*w2i**2*w3i**2*w1i -&
        rhogsol*k**2*cs**2*w3r*w1i*w2r + rhogsol*k**2*cs**2*w3r*w2i*w1r + rhogsol*k**2*cs**2*w3i*w2r*w1r +&
        rhogsol*k**2*cs**2*w2i*w3i*w1i - k*rhogeq*vgsol*w3i*w1r*w2r**2 - k*rhogeq*vgsol*w3i*w1r*w2i**2 -&
        w3r**2*vgsol*k*rhogeq*w2i*w1r + w3r**2*w2r*k*rhogeq*vgsol*w1i + w3r*rhogsol*w1i**3*w2r -&
        vgsol*k*rhogeq*w2i*w1r*w3i**2 - w3r*w1r**2*k*rhogeq*vgsol*w1i - w3r*w1i**3*k*rhogeq*vgsol -&
        w3r*w1r*rhogsol*w2i*w1i**2 + w3r*rhogsol*w1i*w2r*w1r**2 + w2r*vgsol*k*rhogeq*w1i*w3i**2 -&
        rhogsol*w3i*w2r*w1r*w1i**2 - rhogeq*w2r*w1i*k*vgsol*w1r**2 - rhogeq*w2r*w1i**3*k*vgsol +&
        k*rhogeq*vgsol*w2i*w1r**3 - rhogsol*k**2*cs**2*w2i*w1r**2 + k*rhogeq*vgsol*w3i*w1r*w1i**2 +&
        rhogsol*k**2*cs**2*w1i*w1r**2 - rhogsol*k**2*cs**2*w3i*w1i**2 + k*rhogeq*vgsol*w3i*w1r**3 -&
        rhogsol*k**2*cs**2*w1i**2*w2i - rhogsol*k**2*cs**2*w3i*w1r**2 - rhogsol*w3i*w1r**2*w2i*w1i +&
        w1r*k*rhogeq*vgsol*w2i*w1i**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)/(w2r**2 + w1r**2 +&
        w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)/rhogeq/k

 vg1i = - ( - w1r**2*w2i**2*k*rhogeq*vgsol - w2i**2*w1i**2*k*rhogeq*vgsol - w3r**2*w1i**2*k*rhogeq*vgsol -&
        w3r**2*w1r**2*rhogeq*k*vgsol - w3r*w2r*w1r**3*rhogsol - w3r*w1i**2*w2r*rhogsol*w1r -&
        w3i**2*k*rhogeq*w1r**2*vgsol + w3r**2*k*w2i*rhogeq*w1i*vgsol - w1i**2*rhogeq*k*vgsol*w3i**2 +&
        w3i**2*k*w2i*rhogeq*w1i*vgsol + w1i**3*k*rhogeq*vgsol*w3i - w3r*w2i*w1i*rhogsol*w1r**2 - w1i**3*k*Kdrag*vdsol -&
        2*w2i*w1i**2*w3i*rhogeq*k*vgsol + w2i*w3i*w1r*rhogsol*w1i**2 - w1i**2*k*Kdrag*vgsol*w3i -&
        w1i*k*Kdrag*vdsol*w1r**2 + w1i**2*k*rhogeq*vgsol*w3r*w1r + w1i*k*Kdrag*vgsol*w1r**2 - w1r**3*cs**2*k**2*rhogsol&
        + cs**2*k**2*rhogsol*w3r*w1r**2 + w3i*k*Kdrag*vdsol*w1r**2 - w2i*w3i*w1i*k*Kdrag*vdsol +&
        w2i*w3i*w1i*k*Kdrag*vgsol + cs**2*k**2*rhogsol*w3r*w1i**2 + w3i*k*Kdrag*vdsol*w1i**2 +&
        w1i**3*w2i*rhogeq*k*vgsol + w2r*k*rhogeq*w1r*vgsol*w1i**2 - w2i*k*Kdrag*vgsol*w1i**2 + w2i*k*Kdrag*vdsol*w1i**2&
        + w1r**2*rhogeq*k*vgsol*w3i*w1i + w1r**3*rhogeq*k*vgsol*w3r + w2i*k*Kdrag*vdsol*w1r**2 +&
        w2i*w3i*w1r*cs**2*k**2*rhogsol - 2*w2i*rhogeq*k*vgsol*w3i*w1r**2 - w3i*k*Kdrag*vgsol*w1r**2 -&
        w2r*w3i*w1r*k*Kdrag*vdsol + w2r*w3i*w1r*k*Kdrag*vgsol + w2r*w3i**2*k*rhogeq*w1r*vgsol +&
        w3r*w2i**2*w1r*rhogeq*k*vgsol - w2r**2*w3i**2*w1r*rhogsol - w1i**3*w2r*rhogsol*w3i - w2r*w1r**2*rhogsol*w3i*w1i&
        - w2i*k*Kdrag*vgsol*w1r**2 - w2r*w1i*w3i*cs**2*k**2*rhogsol + w2i*w3i*w1r**3*rhogsol +&
        w3i*w2r**2*w1i*rhogeq*k*vgsol + w1i**2*w3i**2*w2r*rhogsol + w1i*w2i*rhogeq*k*vgsol*w1r**2 +&
        w2i**2*w3i*w1i*rhogeq*k*vgsol - w3r*w2i*w1i*cs**2*k**2*rhogsol + w3r*w1i**2*w2r**2*rhogsol +&
        w3r*w2i**2*w1i**2*rhogsol + w3r**2*w1i**2*w2r*rhogsol - w3r**2*w2i**2*w1r*rhogsol -&
        w2r**2*w1r**2*rhogeq*k*vgsol + w3r**2*w2r*w1r**2*rhogsol - w2r**2*w1i**2*k*rhogeq*vgsol -&
        w1r*cs**2*k**2*rhogsol*w1i**2 - w3r**2*w2r**2*w1r*rhogsol + w3r*w1r**2*w2i**2*rhogsol +&
        w3r*w2r**2*w1r**2*rhogsol - w3r*w2i*w1i**3*rhogsol + w1i**3*k*Kdrag*vgsol - w3r*w2i*w1r*k*Kdrag*vdsol +&
        w3r*w2i*w1r*k*Kdrag*vgsol - w3r*w2r*w1r*cs**2*k**2*rhogsol + w3r**2*w2r*k*rhogeq*w1r*vgsol +&
        w3r*w2r**2*k*rhogeq*w1r*vgsol - w3r*w2r*w1i*k*Kdrag*vgsol - 2*w3r*w2r*w1i**2*k*rhogeq*vgsol +&
        w3r*w2r*w1i*k*Kdrag*vdsol - 2*w3r*w2r*w1r**2*rhogeq*k*vgsol + w2r*k*rhogeq*w1r**3*vgsol +&
        w2r*w1r**2*cs**2*k**2*rhogsol - w2i**2*w3i**2*w1r*rhogsol + rhogsol*k**2*cs**2*w1i**2*w2r +&
        w2r*w3i**2*w1r**2*rhogsol)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)/(w2r**2 + w1r**2 +&
        w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)/rhogeq/k

!-------------------------------
! D U S T  V E L O C I T I E S
!-------------------------------
 if (Kdrag > 0.) then

    vd3r = - (rhogeq*cs**2*k**2*w2r**2*w1r**2*rhogsol + rhogeq**2*w3i**4*k*w1r*vgsol -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r - w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r + w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r -&
        rhogeq**2*cs**2*k**3*w2i**2*w1r*vgsol + rhogeq*cs**4*k**4*w2r*w1r*rhogsol -&
        rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol - rhogeq*cs**2*k**3*w2r*w1i*Kdrag*vgsol -&
        rhogeq**2*cs**2*k**3*w2r*w1i**2*vgsol + rhogeq*cs**2*k**3*w2r*w1i*Kdrag*vdsol -&
        rhogeq**2*cs**2*k**3*w2r*w1r**2*vgsol - rhogeq*w3i**4*cs**2*k**2*rhogsol - rhogeq*cs**4*k**4*w2i*w1i*rhogsol +&
        rhogeq*cs**2*k**2*w1i**2*w2r**2*rhogsol + rhogeq**2*w3i**4*w2r*k*vgsol - rhogeq*w3i**3*k*Kdrag*vdsol*w1r +&
        2*rhogeq*w3i**3*k*Kdrag*vgsol*w2r + rhogeq*w2i**2*w1i**2*cs**2*k**2*rhogsol +&
        rhogeq*w2i**2*w1r**2*cs**2*k**2*rhogsol - rhogeq*w3i**4*w2r*w1r*rhogsol + rhogeq*w3i**4*w1i*rhogsol*w2i +&
        w3i*cs**4*k**4*rhogeq*rhogsol*w2i + w3i*cs**4*k**4*rhogeq*rhogsol*w1i + w3i**2*Kdrag**2*k*vgsol*w1r -&
        w3i**3*Kdrag*rhogsol*k**2*cs**2 - w3i**3*Kdrag*rhogsol*w2r*w1r + w3i**3*Kdrag*rhogsol*w2i*w1i -&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1i - w3i*cs**2*k**2*rhogeq*rhogsol*w1i*w2i**2 +&
        rhogeq*w3i**3*rhogsol*k**2*cs**2*w2i + rhogeq*w3i**3*rhogsol*k**2*cs**2*w1i - rhogeq*w3i**3*rhogsol*w2r**2*w1i&
        - rhogeq*w3i**3*rhogsol*w1i*w2i**2 - rhogeq*w3i**3*rhogsol*w2i*w1r**2 -&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w1r**2 - w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w1i**2 -&
        rhogeq*w3i**3*rhogsol*w2i*w1i**2 - rhogeq*w3i**3*k*Kdrag*vdsol*w2r + 2*rhogeq*w3i**3*k*Kdrag*vgsol*w1r +&
        w3r*rhogeq**2*cs**2*k**3*w2i**2*vgsol + w3r*rhogeq*w2r*w3i**2*cs**2*k**2*rhogsol -&
        w3r*rhogeq**2*w3i**2*k*w1r**2*vgsol + w3r*rhogeq*w3i**2*cs**2*k**2*rhogsol*w1r +&
        2*w3r*rhogeq*cs**2*k**3*w3i*Kdrag*vdsol - w3r*rhogeq*cs**4*k**4*w2r*rhogsol +&
        w3r*rhogeq**2*cs**2*k**3*w2r**2*vgsol + w3r**4*rhogeq**2*w2r*vgsol*k - w3r**4*rhogeq*w2r*w1r*rhogsol +&
        w3r**4*rhogeq**2*k*w1r*vgsol - w3r**4*rhogeq*rhogsol*k**2*cs**2 + w3r**4*rhogeq*w2i*rhogsol*w1i +&
        2*w3r*rhogeq*w2i*w3i*w1i*k*Kdrag*vgsol + 2*w3r*rhogeq**2*w2i*k*vgsol*w3i*w1r**2 +&
        2*w3r*rhogeq*w2r*w3i*w1r*k*Kdrag*vdsol - 2*w3r*rhogeq*w2i*w3i**2*k*Kdrag*vgsol - w3r*w3i**2*Kdrag**2*k*vgsol +&
        w3r*w3i**2*Kdrag**2*k*vdsol - w3r*rhogeq**2*w1i**2*k*vgsol*w3i**2 - 2*w3r*rhogeq**2*w3i**2*k*w2i*w1i*vgsol +&
        2*w3r*rhogeq**2*w2i*w1i**2*w3i*k*vgsol + w3r*w3i**2*Kdrag*rhogsol*w1i*w2r + w3r*w3i**2*Kdrag*rhogsol*w2i*w1r -&
        w3r*Kdrag*rhogsol*k**2*cs**2*w1i*w2r - w3r*Kdrag*rhogsol*k**2*cs**2*w2i*w1r +&
        w3r*Kdrag*k*rhogeq*vgsol*w2i*w1r**2 + w3r*Kdrag*k*rhogeq*vgsol*w2i*w1i**2 + w3r*Kdrag**2*k*vdsol*w2r*w1r -&
        w3r*Kdrag**2*k*vdsol*w2i*w1i - w3r*Kdrag**2*k*vgsol*w2r*w1r + w3r*Kdrag**2*k*vgsol*w2i*w1i +&
        w3r*Kdrag*k*rhogeq*vgsol*w2r**2*w1i + w3r*Kdrag*k*rhogeq*vgsol*w1i*w2i**2 - w3r*rhogeq**2*w2r**2*w3i**2*k*vgsol&
        - 2*w3r*rhogeq*w2r*w3i*w1r*k*Kdrag*vgsol + w3r*rhogeq*w2i*w3i**2*k*Kdrag*vdsol -&
        2*w3r*rhogeq**2*w2r*w3i**2*k*w1r*vgsol + w3r*rhogeq*w2r**2*w3i**2*w1r*rhogsol -&
        2*w3r*rhogeq*w1i*w3i**2*k*Kdrag*vgsol + w3r*rhogeq*w1i*w3i**2*k*Kdrag*vdsol +&
        2*w3r*rhogeq**2*w3i*w2r**2*w1i*k*vgsol + w3r*rhogeq*w1i**2*w3i**2*w2r*rhogsol -&
        w3r*rhogeq**2*w3i**2*w2i**2*k*vgsol - w3r*rhogeq*cs**2*k**2*w1i**2*w2r*rhogsol -&
        w3r*rhogeq*cs**2*k**2*w2r*w1r**2*rhogsol - w3r*rhogeq*cs**2*k**2*w2r**2*w1r*rhogsol +&
        w3r*rhogeq*w2i**2*w3i**2*w1r*rhogsol + w3r*rhogeq*w2r*w3i**2*w1r**2*rhogsol -&
        2*w3r*rhogeq**2*cs**2*k**3*w3i*w2i*vgsol - w3r*rhogeq*w2i**2*w1r*cs**2*k**2*rhogsol +&
        2*w3r*rhogeq**2*w2i**2*w3i*w1i*k*vgsol + w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol -&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol + w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol -&
        w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol + w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol -&
        2*w3r*rhogeq**2*cs**2*k**3*w3i*w1i*vgsol + w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol -&
        w3r*rhogeq*cs**4*k**4*w1r*rhogsol - 2*w3r*rhogeq*cs**2*k**3*w3i*Kdrag*vgsol +&
        2*w3r*rhogeq**2*cs**2*k**3*w2r*w1r*vgsol - 2*w3r*rhogeq*w2i*w3i*w1i*k*Kdrag*vdsol +&
        2*w3r*rhogeq**2*cs**2*k**3*w2i*w1i*vgsol + w3r**2*w3i*rhogeq*rhogsol*k**2*cs**2*w1i -&
        w3r**2*w3i*rhogeq*rhogsol*w2r**2*w1i - w3r**2*w3i*rhogeq*rhogsol*w1i*w2i**2 -&
        w3r**2*w3i*rhogeq*rhogsol*w2i*w1r**2 - w3r**2*w3i*rhogeq*rhogsol*w2i*w1i**2 -&
        w3r**2*w3i*rhogeq*k*Kdrag*vdsol*w2r - w3r**2*Kdrag*rhogsol*w3i*w2r*w1r + w3r**2*Kdrag*rhogsol*w2i*w3i*w1i -&
        w3r**2*w3i*rhogeq*k*Kdrag*vdsol*w1r + w3r**2*w3i*rhogeq*rhogsol*k**2*cs**2*w2i -&
        w3r**2*Kdrag*rhogsol*k**2*cs**2*w3i + w3r**2*Kdrag*rhogsol*k**2*cs**2*w2i + w3r**2*Kdrag*rhogsol*k**2*cs**2*w1i&
        + 2*w3r**2*Kdrag*k*rhogeq*vgsol*w2r*w3i + 2*w3r**2*Kdrag*k*rhogeq*vgsol*w3i*w1r +&
        2*w3r**2*rhogeq**2*w2r*k*vgsol*w3i**2 - w3r**2*rhogeq**2*cs**2*k**3*w2r*vgsol +&
        w3r**2*rhogeq*w2r*w1i*k*Kdrag*vgsol + w3r**2*rhogeq**2*w2r*w1i**2*k*vgsol - w3r**2*rhogeq*w2r*w1i*k*Kdrag*vdsol&
        + w3r**2*rhogeq**2*w2r*w1r**2*k*vgsol - w3r**2*rhogeq**2*cs**2*k**3*w1r*vgsol - w3r**2*Kdrag**2*k*vdsol*w2r -&
        w3r**2*rhogeq*w2i**2*w1i**2*rhogsol - w3r**2*rhogeq*w1r**2*w2i**2*rhogsol - w3r**2*rhogeq*w2r**2*w1r**2*rhogsol&
        + w3r**2*rhogeq*cs**4*k**4*rhogsol - w3r**2*rhogeq*w1i**2*w2r**2*rhogsol - w3r**2*Kdrag*rhogsol*w1i*w2i**2 -&
        w3r**2*Kdrag*rhogsol*w2i*w1r**2 - w3r**2*Kdrag*rhogsol*w2i*w1i**2 - w3r**2*Kdrag*rhogsol*w2r**2*w1i -&
        w3r**2*Kdrag**2*k*vdsol*w1r + w3r**2*Kdrag**2*k*vgsol*w2r + w3r**2*Kdrag**2*k*vgsol*w1r +&
        w3r**3*Kdrag*rhogsol*w2i*w1r + w3r**3*rhogeq*w2i*k*Kdrag*vdsol - w3r**3*rhogeq**2*w2i**2*k*vgsol +&
        w3r**2*rhogeq**2*w2i**2*w1r*k*vgsol - 2*w3r**2*rhogeq*w2r*w3i**2*w1r*rhogsol +&
        2*w3r**2*rhogeq**2*w3i**2*k*w1r*vgsol + 2*w3r**2*rhogeq*w1i*w3i**2*rhogsol*w2i -&
        w3r**2*rhogeq*w2i*w1r*k*Kdrag*vdsol + w3r**2*rhogeq*w2i*w1r*k*Kdrag*vgsol - w3r**3*rhogeq**2*w1i**2*k*vgsol -&
        2*w3r**3*rhogeq**2*w2r*k*w1r*vgsol + w3r**3*rhogeq*w2r**2*w1r*rhogsol + w3r**3*rhogeq*w2r*cs**2*k**2*rhogsol -&
        2*w3r**3*rhogeq*w2i*k*Kdrag*vgsol - w3r**3*rhogeq**2*w1r**2*k*vgsol + w3r**3*rhogeq*w1r*cs**2*k**2*rhogsol -&
        w3r**3*rhogeq**2*w2r**2*k*vgsol + w3r**3*rhogeq*w1i*k*Kdrag*vdsol - 2*w3r**3*rhogeq**2*k*w2i*w1i*vgsol +&
        w3r**3*rhogeq*w2r*w1r**2*rhogsol + w3r**3*rhogeq*w1i**2*w2r*rhogsol + w3r**3*rhogeq*w2i**2*w1r*rhogsol -&
        2*w3r**3*rhogeq*w1i*k*Kdrag*vgsol + w3r**3*Kdrag*rhogsol*w1i*w2r + w3r**2*rhogeq**2*w2r**2*k*w1r*vgsol -&
        2*w3r**2*rhogeq*w3i**2*cs**2*k**2*rhogsol + w3r**3*Kdrag**2*k*vdsol + w3i**2*rhogeq**2*cs**2*k**3*w2r*vgsol -&
        w3r**3*Kdrag**2*k*vgsol + w3i**2*rhogeq*w2i*w1r*k*Kdrag*vdsol - w3i**2*Kdrag**2*k*vdsol*w2r +&
        w3i**2*Kdrag**2*k*vgsol*w2r - w3i**2*Kdrag*rhogsol*w2r**2*w1i - w3i**2*Kdrag*rhogsol*w1i*w2i**2 -&
        w3i**2*Kdrag*rhogsol*w2i*w1r**2 - w3i**2*Kdrag*rhogsol*w2i*w1i**2 + w3i**2*Kdrag*rhogsol*k**2*cs**2*w2i +&
        w3i**2*Kdrag*rhogsol*k**2*cs**2*w1i - Kdrag*k*rhogeq*vgsol*w3i*w1r*w2r**2 - Kdrag*k*rhogeq*vgsol*w3i*w2r*w1r**2&
        - Kdrag*k*rhogeq*vgsol*w3i*w2r*w1i**2 - Kdrag*k*rhogeq*vgsol*w3i*w1r*w2i**2 +&
        Kdrag*rhogsol*k**2*cs**2*w3i*w2r*w1r - Kdrag*rhogsol*k**2*cs**2*w2i*w3i*w1i + Kdrag*rhogsol*w2r**2*w1i**2*w3i +&
        Kdrag*rhogsol*w3i*w1r**2*w2r**2 + Kdrag*rhogsol*w1i**2*w2i**2*w3i + Kdrag*rhogsol*w2i**2*w1r**2*w3i +&
        Kdrag**2*k*vdsol*w2r*w1i*w3i + Kdrag**2*k*vdsol*w2i*w3i*w1r - Kdrag**2*k*vgsol*w2r*w1i*w3i -&
        Kdrag**2*k*vgsol*w2i*w3i*w1r - w3i**2*Kdrag**2*k*vdsol*w1r + rhogeq*cs**2*k**3*w2i*w1r*Kdrag*vdsol -&
        rhogeq*cs**2*k**3*w2i*w1r*Kdrag*vgsol - w3i**2*rhogeq**2*w2i**2*w1r*k*vgsol -&
        w3i**2*rhogeq**2*w2r**2*k*w1r*vgsol - w3i**2*rhogeq*w2i*w1r*k*Kdrag*vgsol - w3i**2*rhogeq*w2r*w1i*k*Kdrag*vgsol&
        - w3i**2*rhogeq**2*w2r*w1i**2*k*vgsol + w3i**2*rhogeq*w2r*w1i*k*Kdrag*vdsol -&
        w3i**2*rhogeq**2*w2r*w1r**2*k*vgsol + w3i**2*rhogeq**2*cs**2*k**3*w1r*vgsol +&
        w3i**2*rhogeq*w2i**2*w1i**2*rhogsol + w3i**2*rhogeq*w1r**2*w2i**2*rhogsol + w3i**2*rhogeq*w2r**2*w1r**2*rhogsol&
        - w3i**2*rhogeq*cs**4*k**4*rhogsol + w3i**2*rhogeq*w1i**2*w2r**2*rhogsol)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2&
        + w3i**2 - 2*w3r*w1r)/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i + w3r**2)/k/rhogeq/Kdrag

    vd3i = - (cs**2*k**3*rhogeq*w3i**2*Kdrag*vgsol - cs**2*k**3*rhogeq*w3i**2*Kdrag*vdsol -&
        2*cs**4*k**4*rhogeq*w3r*w3i*rhogsol - 4*cs**2*k**2*rhogeq*w3r*w3i*w2r*w1r*rhogsol +&
        4*cs**2*k**2*rhogeq*w3r*w3i*w2i*rhogsol*w1i - cs**2*k**3*rhogeq*w3r**2*Kdrag*vgsol +&
        cs**2*k**3*rhogeq*w3r**2*Kdrag*vdsol - 2*cs**2*k**3*rhogeq**2*w3i*w2r*w1r*vgsol +&
        2*cs**2*k**3*rhogeq**2*w3r*w3i*w2r*vgsol + 2*cs**2*k**3*rhogeq**2*w3r*w3i*w1r*vgsol -&
        w3r**2*w1i*k*Kdrag**2*vgsol + w3r**2*w1i*k*Kdrag**2*vdsol + rhogeq*w3i*w3r**2*w1i*k*Kdrag*vdsol -&
        rhogeq**2*w3i*w3r**2*w1i**2*k*vgsol + 2*rhogeq**2*w3i**2*w3r**2*w1i*k*vgsol -&
        rhogeq**2*w3i*w3r**2*w1r**2*k*vgsol + rhogeq*w3i*w3r**2*w1r*cs**2*k**2*rhogsol +&
        2*rhogeq*w3i**2*w3r**2*k*Kdrag*vgsol - 2*rhogeq*w3i**2*w3r**2*k*Kdrag*vdsol +&
        rhogeq*w3i*w2r*w3r**2*cs**2*k**2*rhogsol + rhogeq*w3i**3*w2r*cs**2*k**2*rhogsol -&
        rhogeq**2*w3i**3*k*w1r**2*vgsol + rhogeq*w3i**3*cs**2*k**2*rhogsol*w1r - 2*rhogeq**2*w3i*w3r**2*k*w2i*w1i*vgsol&
        - rhogeq**2*w3i**3*w1i**2*k*vgsol - 2*rhogeq**2*w3i**3*k*w2i*w1i*vgsol - rhogeq**2*w3i*w3r**2*w2i**2*k*vgsol +&
        rhogeq**2*w3i**2*w2i*w1i**2*k*vgsol - rhogeq*w3i**2*w2i*w1i*k*Kdrag*vdsol +&
        2*rhogeq**2*w3i**2*w3r**2*w2i*k*vgsol - 2*rhogeq*w3i**2*w3r**2*w2i*w1r*rhogsol -&
        2*rhogeq*w3i**2*w2r*w3r**2*rhogsol*w1i - rhogeq*w3i**4*w2r*w1i*rhogsol + rhogeq*w3i*w3r**2*w2i*k*Kdrag*vdsol +&
        w3r**2*w3i*k*Kdrag**2*vgsol - w3r**2*w3i*k*Kdrag**2*vdsol - w2i*w3i**2*k*Kdrag**2*vgsol +&
        w2i*w3i**2*k*Kdrag**2*vdsol - w1i*w3i**2*k*Kdrag**2*vgsol + w1i*w3i**2*k*Kdrag**2*vdsol +&
        rhogeq**2*w3i**2*w2i*k*vgsol*w1r**2 + rhogeq*w3i**2*w2r*w1r*k*Kdrag*vdsol + rhogeq*w3i**3*w2i*k*Kdrag*vdsol -&
        2*rhogeq**2*w3i**3*w2r*k*w1r*vgsol + rhogeq**2*w3i**4*w2i*k*vgsol + 2*rhogeq**2*w3i*w3r*w2i**2*w1r*k*vgsol +&
        rhogeq*w3i**3*w1i**2*w2r*rhogsol - rhogeq**2*w3i**3*w2i**2*k*vgsol - 2*rhogeq*w3i*w3r*w1i**2*w2r**2*rhogsol -&
        2*rhogeq*w3i*w3r*w2i**2*w1i**2*rhogsol + rhogeq*w3i*w3r**2*w1i**2*w2r*rhogsol +&
        rhogeq*w3i*w3r**2*w2i**2*w1r*rhogsol - rhogeq**2*w3i**3*w2r**2*k*vgsol + rhogeq*w3i*w3r**2*w2r*w1r**2*rhogsol -&
        2*rhogeq**2*w3i*w3r**2*w2r*k*w1r*vgsol + rhogeq*w3i**3*w2r**2*w1r*rhogsol + rhogeq**2*w3i**4*w1i*k*vgsol +&
        rhogeq*w3i**3*w1i*k*Kdrag*vdsol - rhogeq*w3i**4*k*Kdrag*vdsol + rhogeq*w3i**4*k*Kdrag*vgsol -&
        rhogeq*w3i**4*w1r*w2i*rhogsol - Kdrag*w2r*w3r*w3i**2*w1r*rhogsol - Kdrag*w1r*w2i*w3i**3*rhogsol +&
        Kdrag*w2r*w3i**2*cs**2*k**2*rhogsol - Kdrag*w3i**2*k*rhogeq*w1r**2*vgsol + Kdrag*w3i**2*cs**2*k**2*rhogsol*w1r&
        - Kdrag*w1i**2*rhogeq*k*vgsol*w3i**2 + Kdrag*w3r*w1i*w3i**2*rhogsol*w2i - Kdrag*w3r**2*w1i**2*k*rhogeq*vgsol -&
        Kdrag*w3r**2*w1r**2*rhogeq*k*vgsol + Kdrag*w3r**2*w1r*cs**2*k**2*rhogsol + Kdrag*w2r*w3r**2*cs**2*k**2*rhogsol&
        - 2*rhogeq*w3i*w3r*w2r*w1i*k*Kdrag*vdsol + 2*rhogeq**2*w3i*w3r*w2r*w1r**2*k*vgsol +&
        rhogeq*w3i**3*w2i**2*w1r*rhogsol + rhogeq*w3i**3*w2r*w1r**2*rhogsol - Kdrag*w3r*w3i**2*cs**2*k**2*rhogsol -&
        Kdrag*w3r**3*rhogsol*k**2*cs**2 + Kdrag*w3r**3*w2i*rhogsol*w1i - Kdrag*w3r**2*w2i*w3i*w1r*rhogsol -&
        Kdrag*w2r*w3r**2*w3i*rhogsol*w1i - Kdrag*w2r*w1i*w3i**3*rhogsol + 2*rhogeq**2*w3i*w3r*w2r*w1i**2*k*vgsol -&
        Kdrag*w2r*w3r**3*w1r*rhogsol + rhogeq*w3r**3*k*Kdrag*vdsol*w2r - rhogeq**2*w3i*w3r**2*w2r**2*k*vgsol +&
        rhogeq*w3i*w3r**2*w2r**2*w1r*rhogsol - 2*rhogeq*w3i*w3r*w1r**2*w2i**2*rhogsol -&
        2*rhogeq*w3i*w3r*w2r**2*w1r**2*rhogsol - 2*rhogeq*w3i*w3r*w2i*w1r*k*Kdrag*vdsol +&
        2*rhogeq*w3i*w3r*w2i*w1r*k*Kdrag*vgsol + 2*rhogeq**2*w3i*w3r*w2r**2*k*w1r*vgsol +&
        2*rhogeq*w3i*w3r*w2r*w1i*k*Kdrag*vgsol + rhogeq**2*w3r**4*k*vgsol*w2i - rhogeq*w3r**4*rhogsol*w1i*w2r -&
        rhogeq*w3r**4*rhogsol*w2i*w1r + rhogeq*w3r**3*rhogsol*w2r**2*w1i + rhogeq*w3r**3*rhogsol*w1i*w2i**2 +&
        rhogeq*w3r**3*rhogsol*w2i*w1r**2 + rhogeq*w3r**3*rhogsol*w2i*w1i**2 + rhogeq*w3r*rhogsol*w2r**2*w1i*w3i**2 -&
        rhogeq*w3r**4*k*Kdrag*vdsol - rhogeq*w3r**2*k*Kdrag*vdsol*w2r*w1r + rhogeq*w3r**2*k*Kdrag*vdsol*w2i*w1i +&
        rhogeq*w3r*k*Kdrag*vdsol*w2r*w3i**2 + rhogeq*w3r*k*Kdrag*vdsol*w3i**2*w1r -&
        rhogeq*w3r*rhogsol*k**2*cs**2*w3i**2*w1i - rhogeq*w3r*rhogsol*k**2*cs**2*w3i**2*w2i -&
        rhogeq**2*w3r**2*k*vgsol*w2i*w1r**2 - rhogeq**2*w3r**2*k*vgsol*w2i*w1i**2 + rhogeq**2*w3r**4*k*vgsol*w1i +&
        rhogeq*w3r**4*k*Kdrag*vgsol + rhogeq*w3r*rhogsol*w3i**2*w2i*w1r**2 + rhogeq*w3r*rhogsol*w2i*w1i**2*w3i**2 +&
        rhogeq*w3r*rhogsol*w2i**2*w3i**2*w1i - rhogeq*w3r**3*rhogsol*k**2*cs**2*w2i -&
        rhogeq*w3r**3*rhogsol*k**2*cs**2*w1i + rhogeq*w3r**3*k*Kdrag*vdsol*w1r + w3r**2*w2i*k*Kdrag**2*vdsol -&
        w3r**2*w2i*k*Kdrag**2*vgsol - w3i**3*k*Kdrag**2*vdsol + w3i**3*k*Kdrag**2*vgsol -&
        w3r*Kdrag*rhogsol*w2r**2*w1r**2 - w3r*Kdrag*rhogsol*w2r**2*w1i**2 + w3r*Kdrag*k*rhogeq*vgsol*w2r*w1i**2 +&
        w3r*Kdrag*k*rhogeq*vgsol*w2r*w1r**2 - w3r*rhogeq*cs**2*k**2*rhogsol*w2r**2*w1i +&
        w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r - w3r*Kdrag*cs**2*k**2*rhogsol*w2r*w1r -&
        w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r + w3r**2*Kdrag*rhogsol*w2r*w1r**2 +&
        2*w3r**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i - w3r**2*rhogeq**2*k*vgsol*w2r**2*w1i -&
        w3r**2*Kdrag*rhogeq*k*vgsol*w2r*w1r - w3r**2*Kdrag*rhogeq*k*vgsol*w2r**2 - w3r**2*Kdrag*rhogeq*k*vgsol*w2i**2 -&
        w3r**2*rhogeq**2*cs**2*k**3*vgsol*w1i + w3r*Kdrag**2*k*vgsol*w2r*w1i + w3i**2*rhogeq**2*k*vgsol*w2r**2*w1i +&
        rhogeq*cs**2*k**3*vdsol*Kdrag*w2r*w1r + w3i*Kdrag*k*rhogeq*vgsol*w2r**2*w1i +&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1r + w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w1r**2 +&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w1i**2 - 3*w3i**2*Kdrag*rhogeq*k*vgsol*w2r*w1r -&
        w3i**2*Kdrag*rhogeq*k*vgsol*w2i**2 - w3i*rhogeq**2*cs**2*k**3*vgsol*w2r**2 -&
        2*w3i**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i - w3i**2*Kdrag*rhogeq*k*vgsol*w2r**2 -&
        rhogeq*cs**4*k**4*rhogsol*w2r*w1i + w3i*rhogeq*cs**4*k**4*rhogsol*w2r - w3i*Kdrag*cs**2*k**2*rhogsol*w2r*w1i -&
        cs**2*k**3*rhogeq*Kdrag*vgsol*w2r*w1r + w3i**2*cs**2*k**3*rhogeq**2*vgsol*w1i + w3i**2*Kdrag*rhogsol*w2r**2*w1r&
        + w3i*Kdrag**2*k*vdsol*w2r*w1r + w3i**2*Kdrag*rhogsol*w2r*w1r**2 + w3i**2*Kdrag*rhogsol*w2r*w1i**2 -&
        w3i*Kdrag**2*k*vgsol*w2r*w1r + w3r**2*Kdrag*rhogsol*w2r*w1i**2 + w3r**2*Kdrag*rhogsol*w2r**2*w1r +&
        rhogeq**2*cs**2*k**3*vgsol*w1i*w2r**2 + w3r*Kdrag*k*rhogeq*vgsol*w2r**2*w1r - w3r*Kdrag**2*k*vdsol*w2r*w1i -&
        w3r*rhogeq*cs**2*k**2*rhogsol*w2i*w1r**2 + w3r*cs**4*k**4*rhogeq*rhogsol*w1i +&
        w3r*Kdrag*cs**2*k**2*rhogsol*w1i*w2i - w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r - w3i*Kdrag**2*k*vdsol*w1i*w2i +&
        w3i*Kdrag**2*k*vgsol*w1i*w2i + w3r*Kdrag**2*k*vgsol*w2i*w1r - cs**2*k**3*rhogeq*Kdrag*vdsol*w1i*w2i -&
        w3i**2*Kdrag*rhogeq*k*vgsol*w1i*w2i + w3r*Kdrag*k*rhogeq*vgsol*w2i**2*w1r + w3i*rhogeq*k*Kdrag*vgsol*w2i**2*w1i&
        + w3i*rhogeq*k*Kdrag*vgsol*w2i*w1r**2 + w3i*rhogeq*k*Kdrag*vgsol*w1i**2*w2i -&
        w3i*rhogeq**2*cs**2*k**3*vgsol*w1r**2 - w3i*rhogeq**2*cs**2*k**3*vgsol*w1i**2 -&
        3*w3r**2*Kdrag*k*rhogeq*vgsol*w1i*w2i + w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r -&
        w3i*Kdrag*cs**2*k**2*rhogsol*w2i*w1r - rhogeq*cs**4*k**4*rhogsol*w2i*w1r +&
        2*w3r**2*rhogeq*cs**2*k**2*rhogsol*w2i*w1r - w3r**2*vgsol*rhogeq**2*k**3*cs**2*w2i +&
        w3i*rhogeq*cs**2*k**3*vdsol*Kdrag*w2i + w3i*rhogeq*cs**2*k**3*vdsol*Kdrag*w1i +&
        cs**2*k**3*rhogeq*Kdrag*vgsol*w1i*w2i - w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1i -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i + w3i*rhogeq*cs**4*k**4*rhogsol*w1r +&
        w3i**2*vgsol*rhogeq**2*k**3*cs**2*w2i + w3i**2*vgsol*k*rhogeq**2*w2i**2*w1i -&
        2*w3i**2*rhogeq*cs**2*k**2*rhogsol*w2i*w1r - w3r**2*vgsol*k*rhogeq**2*w2i**2*w1i -&
        w3i*rhogeq**2*cs**2*k**3*vgsol*w2i**2 - 2*w3i*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i -&
        w3r*Kdrag*rhogsol*w2i**2*w1r**2 + rhogeq**2*cs**2*k**3*vgsol*w1i**2*w2i + rhogeq**2*cs**2*k**3*vgsol*w2i*w1r**2&
        - w3r*Kdrag**2*k*vdsol*w2i*w1r + w3r*cs**4*k**4*rhogeq*rhogsol*w2i - w3r*rhogeq*cs**2*k**2*rhogsol*w1i**2*w2i -&
        w3r*rhogeq*cs**2*k**2*rhogsol*w2i**2*w1i - w3r*Kdrag*rhogsol*w2i**2*w1i**2 +&
        w3i*rhogeq*cs**2*k**2*rhogsol*w2i**2*w1r + w3i**2*Kdrag*rhogsol*w2i**2*w1r + w3r**2*Kdrag*rhogsol*w2i**2*w1r +&
        rhogeq**2*cs**2*k**3*vgsol*w1i*w2i**2)/Kdrag/rhogeq/k/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i +&
        w3r**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)

    vd2r = - ( - rhogeq*w1i*rhogsol*k**2*cs**2*w3i**2*w2i + rhogeq*w1i**2*rhogsol*k**2*cs**2*w3r**2 +&
        rhogeq*w1i**2*rhogsol*k**2*cs**2*w3i**2 - rhogeq*w1i*rhogsol*k**2*cs**2*w3r**2*w2i +&
        rhogeq**2*w1i**2*k*vgsol*w3r*w2r**2 - rhogeq**2*w1i**2*k*vgsol*w3r*w2i**2 - rhogeq**2*w2r**3*k*vgsol*w3i**2 +&
        rhogeq*w2r**3*w3i**2*w1r*rhogsol + 2*rhogeq**2*w1r*w2r**2*w2i**2*k*vgsol - rhogeq**2*w1r**2*w2r*k*vgsol*w2i**2&
        - 2*cs**2*k**3*rhogeq**2*w2r*vgsol*w2i*w1i + cs**2*k**2*rhogeq*w2r*w1r*rhogsol*w2i**2 +&
        cs**2*k**2*rhogeq*w2r**3*w1r*rhogsol - w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r + w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r + rhogeq**2*cs**2*k**3*w2i**2*w1r*vgsol -&
        rhogeq*cs**4*k**4*w2r*w1r*rhogsol - rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol +&
        rhogeq*cs**2*k**3*w2r*w1i*Kdrag*vgsol + rhogeq**2*cs**2*k**3*w2r*w1i**2*vgsol -&
        rhogeq*cs**2*k**3*w2r*w1i*Kdrag*vdsol + rhogeq**2*cs**2*k**3*w2r*w1r**2*vgsol +&
        rhogeq*cs**4*k**4*w2i*w1i*rhogsol + w3i*cs**4*k**4*rhogeq*rhogsol*w2i - w3i*cs**4*k**4*rhogeq*rhogsol*w1i -&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w1r**2 - w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w1i**2 +&
        w3r*rhogeq**2*cs**2*k**3*w2i**2*vgsol - w3r*rhogeq*cs**4*k**4*w2r*rhogsol -&
        w3r*rhogeq**2*cs**2*k**3*w2r**2*vgsol - w3r*Kdrag*rhogsol*k**2*cs**2*w1i*w2r +&
        w3r*Kdrag*rhogsol*k**2*cs**2*w2i*w1r - w3r*Kdrag*k*rhogeq*vgsol*w2i*w1r**2 -&
        w3r*Kdrag*k*rhogeq*vgsol*w2i*w1i**2 + w3r*Kdrag**2*k*vdsol*w2r*w1r + w3r*Kdrag**2*k*vdsol*w2i*w1i -&
        w3r*Kdrag**2*k*vgsol*w2r*w1r - w3r*Kdrag**2*k*vgsol*w2i*w1i + w3r*Kdrag*k*rhogeq*vgsol*w2r**2*w1i -&
        w3r*Kdrag*k*rhogeq*vgsol*w1i*w2i**2 - w3r*rhogeq*cs**2*k**2*w1i**2*w2r*rhogsol -&
        w3r*rhogeq*cs**2*k**2*w2r*w1r**2*rhogsol - w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol -&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol + w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol +&
        w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol - w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol -&
        w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol + w3r*rhogeq*cs**4*k**4*w1r*rhogsol +&
        2*w3r*rhogeq**2*cs**2*k**3*w2r*w1r*vgsol + w3r**2*rhogeq**2*cs**2*k**3*w2r*vgsol +&
        w3r**2*rhogeq*w2r*w1i*k*Kdrag*vgsol - w3r**2*rhogeq**2*cs**2*k**3*w1r*vgsol +&
        w3r**2*rhogeq*w2i**2*w1i**2*rhogsol + w3r**2*rhogeq*w1r**2*w2i**2*rhogsol - w3r**2*rhogeq*w2r**2*w1r**2*rhogsol&
        - w3r**2*rhogeq*w1i**2*w2r**2*rhogsol - w3r**2*Kdrag*rhogsol*w1i*w2i**2 + w3r**2*Kdrag*rhogsol*w2i*w1r**2 +&
        w3r**2*Kdrag*rhogsol*w2i*w1i**2 - w3r**2*Kdrag*rhogsol*w2r**2*w1i - w3r**2*rhogeq**2*w2i**2*w1r*k*vgsol -&
        w3r**2*rhogeq*w2i*w1r*k*Kdrag*vgsol + w3r**2*rhogeq**2*w2r**2*k*w1r*vgsol +&
        w3i**2*rhogeq**2*cs**2*k**3*w2r*vgsol - w3i**2*Kdrag*rhogsol*w2r**2*w1i - w3i**2*Kdrag*rhogsol*w1i*w2i**2 +&
        w3i**2*Kdrag*rhogsol*w2i*w1r**2 + w3i**2*Kdrag*rhogsol*w2i*w1i**2 + Kdrag*k*rhogeq*vgsol*w3i*w1r*w2r**2 +&
        Kdrag*k*rhogeq*vgsol*w3i*w2r*w1r**2 + Kdrag*k*rhogeq*vgsol*w3i*w2r*w1i**2 - Kdrag*k*rhogeq*vgsol*w3i*w1r*w2i**2&
        - Kdrag*rhogsol*k**2*cs**2*w3i*w2r*w1r - Kdrag*rhogsol*k**2*cs**2*w2i*w3i*w1i - Kdrag*rhogsol*w2r**2*w1i**2*w3i&
        - Kdrag*rhogsol*w3i*w1r**2*w2r**2 - Kdrag*rhogsol*w1i**2*w2i**2*w3i - Kdrag*rhogsol*w2i**2*w1r**2*w3i -&
        Kdrag**2*k*vdsol*w2r*w1i*w3i + Kdrag**2*k*vdsol*w2i*w3i*w1r + Kdrag**2*k*vgsol*w2r*w1i*w3i -&
        Kdrag**2*k*vgsol*w2i*w3i*w1r - rhogeq*cs**2*k**3*w2i*w1r*Kdrag*vdsol + rhogeq*cs**2*k**3*w2i*w1r*Kdrag*vgsol -&
        w3i**2*rhogeq**2*w2i**2*w1r*k*vgsol + w3i**2*rhogeq**2*w2r**2*k*w1r*vgsol - w3i**2*rhogeq*w2i*w1r*k*Kdrag*vgsol&
        + w3i**2*rhogeq*w2r*w1i*k*Kdrag*vgsol - w3i**2*rhogeq**2*cs**2*k**3*w1r*vgsol +&
        w3i**2*rhogeq*w2i**2*w1i**2*rhogsol + w3i**2*rhogeq*w1r**2*w2i**2*rhogsol - w3i**2*rhogeq*w2r**2*w1r**2*rhogsol&
        - w3i**2*rhogeq*w1i**2*w2r**2*rhogsol - 2*w3i*rhogeq*k*Kdrag*vdsol*w2r*w2i*w1i +&
        w3i*rhogeq*rhogsol*k**2*cs**2*w2i*w2r**2 + 2*rhogeq**2*w1i**2*k*vgsol*w2r*w3i*w2i -&
        rhogeq*w2i*w1r*k*Kdrag*vdsol*w2r**2 - 2*vgsol*Kdrag*rhogeq*k**3*cs**2*w2r*w2i -&
        2*cs**2*k**2*rhogeq*rhogsol*w2i**2*w2r**2 + 2*vdsol*Kdrag*rhogeq*k**3*cs**2*w2r*w2i -&
        cs**2*k**2*rhogeq*rhogsol*w2r**4 - rhogeq**2*w2r*k*vgsol*w3i**2*w2i**2 - w3i*rhogeq*k*Kdrag*vdsol*w1r*w2r**2 +&
        w3i*rhogeq*k*Kdrag*vdsol*w2r*w2i**2 - w3i*rhogeq*rhogsol*w2i*w1r**2*w2r**2 +&
        2*w3i*rhogeq*rhogsol*w2r**2*w1i*w2i**2 + rhogeq*cs**4*k**4*rhogsol*w2r**2 + w3i*rhogeq*k*Kdrag*vdsol*w2r**3 +&
        w3i*rhogeq*rhogsol*w2r**4*w1i - rhogeq*w1i*w3i**2*rhogsol*w2i*w2r**2 + rhogeq*w2r*w3i**2*w1r*rhogsol*w2i**2 +&
        2*w3i*rhogeq**2*k*vgsol*w2i*w2r*w1r**2 - 2*w3i*cs**2*k**3*rhogeq**2*vgsol*w2r*w2i +&
        2*w3i*cs**2*k**3*rhogeq**2*vgsol*w2r*w1i - rhogeq*w3i**2*cs**2*k**2*rhogsol*w2r*w1r +&
        2*rhogeq**2*w2r*k*vgsol*w3i**2*w2i*w1i + 2*rhogeq*w2r*w1i*k*Kdrag*vgsol*w2i*w3i +&
        w3r*rhogeq*w2r**3*cs**2*k**2*rhogsol + w3r*rhogeq*w2r**3*w1r**2*rhogsol + w3r*rhogeq**2*w2r**4*k*vgsol +&
        w3r*rhogeq**2*k*vgsol*w2i**4 - w3r*vdsol*Kdrag*k*rhogeq*w2i**3 - w3i*rhogeq*rhogsol*w2i**3*w1r**2 +&
        w3i*rhogeq*rhogsol*w2i**4*w1i - w3r**2*rhogeq*rhogsol*w2i**3*w1i - w3i**2*rhogeq*rhogsol*w2i**3*w1i -&
        w3r*rhogeq**2*k*vgsol*w2i**2*w1r**2 + w3i*vdsol*Kdrag*k*rhogeq*w2i**2*w1r - rhogeq*rhogsol*k**2*cs**2*w2i**4 +&
        w3r**2*cs**2*k**2*rhogeq*rhogsol*w1r**2 - rhogeq*k*Kdrag*vdsol*w2i**3*w1r +&
        w3i**2*cs**2*k**2*rhogeq*rhogsol*w1r**2 + w3i*cs**2*k**2*rhogeq*rhogsol*w2i**3 - w3r*rhogeq*rhogsol*w2i**4*w1r&
        - rhogeq**2*w1i**2*w2r**3*vgsol*k - rhogeq*w1i**2*w2i**3*rhogsol*w3i - cs**4*k**4*rhogeq*rhogsol*w2i**2 +&
        w2i**2*k*Kdrag**2*vdsol*w2r + w2r**2*k*Kdrag**2*vgsol*w1r - w2r**2*k*Kdrag**2*vdsol*w1r -&
        w2i**2*k*Kdrag**2*vgsol*w2r - w2i**2*k*Kdrag**2*vdsol*w1r + w2r**3*k*Kdrag**2*vdsol +&
        w3r*w2r**2*k*Kdrag**2*vgsol - w3r*w2r**2*k*Kdrag**2*vdsol + w2i**2*k*Kdrag**2*vgsol*w1r +&
        w3r*w2i**2*k*Kdrag**2*vgsol - w3r*w2i**2*k*Kdrag**2*vdsol - w2r**3*k*Kdrag**2*vgsol -&
        2*Kdrag*w2r**3*k*rhogeq*vgsol*w3i + Kdrag*w2i**2*rhogsol*k**2*cs**2*w3i - 2*Kdrag*w2i**2*k*rhogeq*vgsol*w2r*w3i&
        - 2*Kdrag*w2r**3*vgsol*k*rhogeq*w1i + 2*Kdrag*w2i*vgsol*k*rhogeq*w1r*w2r**2 + 2*Kdrag*w2i**3*vgsol*k*rhogeq*w1r&
        - 2*Kdrag*w2i**2*w2r*vgsol*k*rhogeq*w1i + Kdrag*w2i**3*rhogsol*w3i*w1i - Kdrag*w3r*w2i**3*rhogsol*w1r +&
        Kdrag*w3r*w2r**3*rhogsol*w1i - Kdrag*w3r*w2r**2*rhogsol*w2i*w1r + Kdrag*w3r*w2i**2*rhogsol*w1i*w2r +&
        Kdrag*w2r**3*rhogsol*w3i*w1r + rhogeq*w1i*w2i**2*k*Kdrag*vdsol*w2r + rhogeq*w1i*w2r**2*rhogsol*k**2*cs**2*w2i +&
        Kdrag*w2r**2*rhogsol*w2i*w3i*w1i + Kdrag*w2r**2*rhogsol*k**2*cs**2*w1i + Kdrag*w2r**2*rhogsol*k**2*cs**2*w3i -&
        2*rhogeq**2*w1i*w2i**2*k*vgsol*w2r*w3i - rhogeq*w1i*w3r*w2r**2*k*Kdrag*vdsol -&
        2*rhogeq**2*w1i*w2r**3*k*vgsol*w3i + rhogeq*w1i*w2i**3*rhogsol*k**2*cs**2 + rhogeq*w1i*w3r*w2i**2*k*Kdrag*vdsol&
        + rhogeq*w1i**2*w3r*w2i**2*rhogsol*w2r + rhogeq*w1i*w2r**3*k*Kdrag*vdsol + rhogeq*w1i**2*w3r*w2r**3*rhogsol -&
        rhogeq*w1i**2*w2r**2*rhogsol*w2i*w3i - rhogeq**2*w1i**2*w2i**2*w2r*vgsol*k +&
        Kdrag*w2i**2*rhogsol*k**2*cs**2*w1i + Kdrag*w2i**2*rhogsol*w3i*w2r*w1r - Kdrag*w2r**2*rhogsol*k**2*cs**2*w2i -&
        Kdrag*w2i**3*rhogsol*k**2*cs**2 + 2*Kdrag*w3r*w2r**2*k*rhogeq*vgsol*w2i + 2*Kdrag*w3r*w2i**3*k*rhogeq*vgsol -&
        2*w3r*w2r*Kdrag*k*rhogeq*vgsol*w2i*w1r + 2*w3r*rhogeq*w2i*k*Kdrag*vdsol*w2r*w1r -&
        2*w3r*w2r*rhogeq**2*vgsol*k*w2i**2*w1r + w3r*rhogeq**2*w1r**2*k*vgsol*w2r**2 +&
        w3r*rhogeq*w2r*cs**2*k**2*rhogsol*w2i**2 - w3r*rhogeq*w2r**4*w1r*rhogsol + w3r*rhogeq*w2r*w1r**2*rhogsol*w2i**2&
        - 2*w3r*rhogeq**2*w2r**3*k*w1r*vgsol - 2*w3r*rhogeq*w2r**2*w1r*rhogsol*w2i**2 -&
        w3r*rhogeq*w2i*k*Kdrag*vdsol*w2r**2 + 2*w3r*rhogeq**2*w2i**2*k*vgsol*w2r**2 +&
        w3r**2*rhogeq*w2r*w1r*rhogsol*w2i**2 - w3r**2*rhogeq**2*w2r*vgsol*k*w2i**2 -&
        w3r**2*rhogeq*w2i*rhogsol*w1i*w2r**2 - w3r**2*rhogeq**2*w2r**3*vgsol*k + 2*w3r**2*rhogeq**2*w2r*vgsol*k*w2i*w1i&
        - w3r**2*rhogeq*rhogsol*k**2*cs**2*w2r*w1r + w3r**2*rhogeq*w2r**3*w1r*rhogsol + rhogeq**2*w1r*w2i**4*k*vgsol +&
        rhogeq**2*w1r*w2r**4*k*vgsol - rhogeq**2*w1r**2*w2r**3*k*vgsol)/k/rhogeq/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2&
        - 2*w2i*w3i + w3r**2)/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)/Kdrag

    vd2i =1/k*(w3r*rhogeq*k*Kdrag*vdsol*w2r**2*w1r - w3r*rhogeq*k*Kdrag*vdsol*w2r*w2i**2 -&
        w3r*rhogeq*k*Kdrag*vdsol*w2r**3 - w3r*Kdrag*rhogsol*w2r**2*w1r**2 - w3r*Kdrag*rhogsol*w2r**2*w1i**2 +&
        w3r*rhogeq*rhogsol*w2r**4*w1i - w3r*Kdrag*k*rhogeq*vgsol*w2r*w1i**2 - w3r*Kdrag*k*rhogeq*vgsol*w2r*w1r**2 -&
        2*w3r*rhogeq*cs**2*k**2*rhogsol*w2r**2*w1i - w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r -&
        2*w3r*rhogeq**2*k*vgsol*w2r*w2i*w1i**2 - 2*w3r*rhogeq**2*k*vgsol*w2r*w2i*w1r**2 +&
        2*w3r*rhogeq**2*k*vgsol*w2r**2*w2i*w1r - w3r*rhogeq*cs**2*k**2*rhogsol*w2i*w2r**2 +&
        w3r*Kdrag*cs**2*k**2*rhogsol*w2r*w1r - w3r*Kdrag*cs**2*k**2*rhogsol*w2r**2 + w3r*Kdrag*rhogsol*w2r**2*w1i*w2i -&
        2*w3r*cs**2*k**3*rhogeq**2*vgsol*w2i*w2r + w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r +&
        4*w3r*rhogeq*cs**2*k**2*rhogsol*w2r*w2i*w1r + w3r*Kdrag*rhogsol*w2r*w1r*w2i**2 +&
        2*w3r*rhogeq*rhogsol*w2i**2*w2r**2*w1i + w3r**2*Kdrag*rhogsol*w2r*w1r**2 +&
        2*w3r**2*rhogeq*rhogsol*w2r*w2i*w1i**2 - w3r*rhogeq*rhogsol*w1i**2*w2r**2*w2i -&
        w3r*rhogeq*rhogsol*w2i*w1r**2*w2r**2 - 2*w3r*w2r*Kdrag*rhogeq*k*vgsol*w2i*w1i +&
        w3r**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i + w3r**2*rhogeq**2*k*vgsol*w2i*w2r**2 +&
        w3r**2*rhogeq**2*k*vgsol*w2r**2*w1i - w3r**2*Kdrag*rhogeq*k*vgsol*w2r*w1r + w3r**2*Kdrag*rhogeq*k*vgsol*w2r**2&
        - w3r**2*rhogeq*rhogsol*w2r*w2i**2*w1i - w3r**2*rhogeq*rhogsol*w2r**2*w2i*w1r +&
        w3r**2*Kdrag*rhogeq*k*vgsol*w2i**2 - w3r**2*rhogeq**2*cs**2*k**3*vgsol*w1i +&
        2*w3r**2*rhogeq*rhogsol*w2r*w2i*w1r**2 - 2*w3r**2*rhogeq**2*k*vgsol*w2r*w2i*w1r - w3r*Kdrag**2*k*vgsol*w2r*w1i&
        + rhogeq*k*Kdrag*vdsol*w2r**4 - w3i*rhogeq*rhogsol*w2r**3*w1i**2 + Kdrag**2*k*vdsol*w2i*w2r**2 -&
        2*vgsol*k*rhogeq**2*w2i**2*w2r**2*w1i + vgsol*k*rhogeq**2*w2i*w1r**2*w2r**2 -&
        2*vgsol*rhogeq**2*k**3*cs**2*w2r*w2i*w1r - rhogeq*cs**2*k**2*rhogsol*w2r**2*w2i*w1r +&
        rhogeq*cs**2*k**2*rhogsol*w2r**3*w1i + rhogeq*cs**2*k**2*rhogsol*w2r*w2i**2*w1i +&
        w3i*rhogeq**2*k*vgsol*w2r**2*w1r**2 - w3i*rhogeq**2*k*vgsol*w2r**4 + w3i**2*rhogeq**2*k*vgsol*w2i*w2r**2 -&
        2*w3i**2*rhogeq**2*k*vgsol*w2r*w2i*w1r + w3i**2*rhogeq**2*k*vgsol*w2r**2*w1i -&
        rhogeq*cs**2*k**3*vdsol*Kdrag*w2r**2 + rhogeq*cs**2*k**3*vdsol*Kdrag*w2r*w1r - Kdrag**2*k*vdsol*w2r**2*w1i +&
        3*w3i*Kdrag*k*rhogeq*vgsol*w2r**2*w1i - 2*w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1r +&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w2i**2 + w3i*cs**2*k**2*rhogeq*rhogsol*w2r**3 +&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w1r**2 - 4*w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w2i*w1i +&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w1i**2 + w3i*rhogeq**2*k*vgsol*w2r**2*w1i**2 -&
        2*w3i*rhogeq**2*k*vgsol*w2i**2*w2r**2 - w3i**2*Kdrag*rhogeq*k*vgsol*w2r*w1r +&
        w3i**2*Kdrag*rhogeq*k*vgsol*w2i**2 + Kdrag**2*k*vgsol*w2r**2*w1i + w3i*rhogeq**2*cs**2*k**3*vgsol*w2r**2 +&
        Kdrag*k*rhogeq*vgsol*w2r**2*w1r**2 + Kdrag*k*rhogeq*vgsol*w2r**2*w1i**2 +&
        w3i**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i - Kdrag**2*k*vgsol*w2i*w2r**2 + 2*rhogeq*cs**4*k**4*rhogsol*w2i*w2r +&
        2*w3i*rhogeq*rhogsol*w2i**2*w2r**2*w1r + w3i**2*Kdrag*rhogeq*k*vgsol*w2r**2 + w3i*rhogeq*rhogsol*w2r**4*w1r -&
        w3i*rhogeq*rhogsol*w2r**3*w1r**2 - rhogeq*cs**4*k**4*rhogsol*w2r*w1i - w3i*rhogeq*cs**4*k**4*rhogsol*w2r -&
        rhogeq*k*Kdrag*vdsol*w2r**2*w1i*w2i + Kdrag*cs**2*k**2*rhogsol*w2r**3 + 2*rhogeq*k*Kdrag*vdsol*w2i**2*w2r**2 -&
        rhogeq*k*Kdrag*vdsol*w2r**3*w1r + Kdrag*cs**2*k**2*rhogsol*w2r*w2i**2 - rhogeq*k*Kdrag*vdsol*w2r*w1r*w2i**2 -&
        Kdrag*cs**2*k**2*rhogsol*w2r**2*w1r + w3i*Kdrag*rhogsol*w2r**2*w2i*w1r - w3i*Kdrag*rhogsol*w2r*w2i**2*w1i -&
        w3i*Kdrag*rhogsol*w2r**3*w1i - w3i*Kdrag*cs**2*k**2*rhogsol*w2r*w1i - w3i**2*rhogeq*rhogsol*w2r**2*w2i*w1r -&
        w3i**2*rhogeq*rhogsol*w2r*w2i**2*w1i - w3i*rhogeq*rhogsol*w2r*w2i**2*w1r**2 +&
        2*w3i**2*rhogeq*rhogsol*w2r*w2i*w1r**2 - w3i**2*rhogeq*rhogsol*w2r**3*w1i +&
        2*w3i**2*rhogeq*rhogsol*w2r*w2i*w1i**2 - w3i*rhogeq*rhogsol*w2r*w2i**2*w1i**2 -&
        cs**2*k**3*rhogeq*Kdrag*vgsol*w2r*w1r + cs**2*k**3*rhogeq*Kdrag*vgsol*w2r**2 -&
        w3i**2*cs**2*k**3*rhogeq**2*vgsol*w1i - w3i**2*Kdrag*rhogsol*w2r**2*w1r + w3i*Kdrag**2*k*vdsol*w2r*w1r +&
        w3i**2*Kdrag*rhogsol*w2r*w1r**2 + w3i**2*Kdrag*rhogsol*w2r*w1i**2 - w3i*Kdrag**2*k*vdsol*w2r**2 +&
        w3i*Kdrag**2*k*vgsol*w2r**2 - w3i*Kdrag**2*k*vgsol*w2r*w1r + 2*w3i*rhogeq*k*Kdrag*vdsol*w2r*w2i*w1r -&
        w3i*rhogeq*k*Kdrag*vdsol*w2r**2*w1i - w3i*rhogeq*k*Kdrag*vdsol*w2i*w2r**2 + w3r**2*Kdrag*rhogsol*w2r*w1i**2 -&
        w3r**2*Kdrag*rhogsol*w2r**2*w1r + w3r*Kdrag*rhogsol*w2r**3*w1r - w3r**2*rhogeq*rhogsol*w2r**3*w1i +&
        Kdrag**2*k*vdsol*w2i**3 - Kdrag**2*k*vgsol*w2i**3 - Kdrag*rhogeq*k*vgsol*w2r**4 +&
        rhogeq**2*cs**2*k**3*vgsol*w1i*w2r**2 - 2*Kdrag*rhogeq*k*vgsol*w2i**2*w2r**2 -&
        2*Kdrag*rhogeq*k*vgsol*w2r*w1r*w2i*w3i + 2*rhogeq**2*k*vgsol*w2r**2*w1i*w2i*w3i +&
        2*w3r*rhogeq*k*Kdrag*vdsol*w2r*w2i*w1i + w3r*Kdrag*k*rhogeq*vgsol*w2r**2*w1r + w3r*Kdrag**2*k*vdsol*w2r*w1i -&
        rhogeq**2*k*vgsol*w2r**4*w1i - w3r*rhogeq*k*Kdrag*vdsol*w2i**2*w1r - w3r*rhogeq*cs**2*k**2*rhogsol*w2i*w1r**2 -&
        w3r*rhogeq*cs**2*k**2*rhogsol*w2i**3 + 2*w3r*vgsol*k**3*cs**2*rhogeq**2*w2i*w1r +&
        w3r*cs**4*k**4*rhogeq*rhogsol*w1i + w3r*Kdrag*cs**2*k**2*rhogsol*w1i*w2i -&
        w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r + w3i*Kdrag**2*k*vdsol*w1i*w2i - w3i*Kdrag**2*k*vdsol*w2i**2 +&
        rhogeq*k*Kdrag*vdsol*w2i**4 - w3i*Kdrag**2*k*vgsol*w1i*w2i + w3i*Kdrag**2*k*vgsol*w2i**2 +&
        w3r*Kdrag**2*k*vgsol*w2i*w1r + 2*w3r*rhogeq**2*k*vgsol*w2i**3*w1r - cs**2*k**3*rhogeq*Kdrag*vdsol*w1i*w2i -&
        w3i**2*Kdrag*rhogeq*k*vgsol*w1i*w2i - Kdrag*rhogeq*k*vgsol*w2i**4 + Kdrag*rhogeq*k*vgsol*w2i**2*w1r**2 +&
        Kdrag*rhogeq*k*vgsol*w2i**2*w1i**2 + 3*w3r*Kdrag*k*rhogeq*vgsol*w2i**2*w1r +&
        w3i*rhogeq*k*Kdrag*vgsol*w2i**2*w1i - w3i*rhogeq*k*Kdrag*vgsol*w2i*w1r**2 - w3i*rhogeq*k*Kdrag*vgsol*w1i**2*w2i&
        - w3i*rhogeq*k*Kdrag*vdsol*w2i**3 + w3i*rhogeq*k*Kdrag*vdsol*w2i**2*w1i - w3i*rhogeq**2*k*vgsol*w2i**2*w1i**2 +&
        2*w3i*rhogeq**2*k*vgsol*w2i**3*w1i - w3i*rhogeq**2*k*vgsol*w2i**4 - w3i*rhogeq**2*cs**2*k**3*vgsol*w1r**2 -&
        w3i*rhogeq**2*cs**2*k**3*vgsol*w1i**2 - w3r**2*Kdrag*k*rhogeq*vgsol*w1i*w2i +&
        w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r - rhogeq*cs**2*k**2*rhogsol*w2i**3*w1r +&
        w3i*Kdrag*cs**2*k**2*rhogsol*w2i*w1r - rhogeq*cs**4*k**4*rhogsol*w2i*w1r -&
        w3r**2*rhogeq*cs**2*k**2*rhogsol*w2i*w1r + w3r**2*vgsol*rhogeq**2*k**3*cs**2*w2i -&
        w3i*rhogeq*cs**2*k**3*vdsol*Kdrag*w2i + w3i*rhogeq*cs**2*k**3*vdsol*Kdrag*w1i +&
        cs**2*k**3*rhogeq*Kdrag*vgsol*w1i*w2i - w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1i +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i + w3i*rhogeq*cs**4*k**4*rhogsol*w1r +&
        w3i**2*vgsol*rhogeq**2*k**3*cs**2*w2i - w3i**2*vgsol*k*rhogeq**2*w2i**2*w1i -&
        w3i**2*rhogeq*cs**2*k**2*rhogsol*w2i*w1r + w3i**2*vgsol*k*rhogeq**2*w2i**3 -&
        w3r**2*vgsol*k*rhogeq**2*w2i**2*w1i + w3r**2*vgsol*k*rhogeq**2*w2i**3 - w3i*rhogeq**2*cs**2*k**3*vgsol*w2i**2 +&
        2*w3i*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i - cs**2*k**3*rhogeq*Kdrag*vgsol*w2i**2 +&
        cs**2*k**3*rhogeq*Kdrag*vdsol*w2i**2 - w3i*rhogeq**2*k*vgsol*w2i**2*w1r**2 - rhogeq*k*Kdrag*vdsol*w2i**3*w1i -&
        Kdrag*cs**2*k**2*rhogsol*w2i**2*w1r + w3r*Kdrag*rhogsol*w2i**3*w1i - w3r*Kdrag*rhogsol*w2i**2*w1r**2 +&
        w3i*Kdrag*rhogsol*w2i**3*w1r + rhogeq**2*cs**2*k**3*vgsol*w1i**2*w2i + rhogeq**2*cs**2*k**3*vgsol*w2i*w1r**2 -&
        w3r*Kdrag**2*k*vdsol*w2i*w1r + vgsol*k*rhogeq**2*w1i**2*w2r**2*w2i - w3r*cs**4*k**4*rhogeq*rhogsol*w2i -&
        w3r*rhogeq*cs**2*k**2*rhogsol*w1i**2*w2i + 2*w3r*rhogeq*cs**2*k**2*rhogsol*w2i**2*w1i -&
        w3r*Kdrag*rhogsol*w2i**2*w1i**2 + 2*w3i*rhogeq*cs**2*k**2*rhogsol*w2i**2*w1r - w3r**2*rhogeq*rhogsol*w2i**3*w1r&
        - w3i**2*rhogeq*rhogsol*w2i**3*w1r - w3i**2*Kdrag*rhogsol*w2i**2*w1r - w3r**2*Kdrag*rhogsol*w2i**2*w1r -&
        rhogeq**2*k*vgsol*w2i**4*w1i + rhogeq**2*k*vgsol*w2i**3*w1r**2 - w3r*rhogeq*rhogsol*w2i**3*w1i**2 +&
        w3r*rhogeq*rhogsol*w2i**4*w1i - w3r*rhogeq*rhogsol*w2i**3*w1r**2 + rhogeq**2*k*vgsol*w2i**3*w1i**2 -&
        Kdrag**2*k*vdsol*w2i**2*w1i + w3i*rhogeq*rhogsol*w2i**4*w1r + Kdrag**2*k*vgsol*w2i**2*w1i -&
        w3r*Kdrag*cs**2*k**2*rhogsol*w2i**2 - rhogeq**2*cs**2*k**3*vgsol*w1i*w2i**2)/(w2r**2 - 2*w3r*w2r + w2i**2 +&
        w3i**2 - 2*w2i*w3i + w3r**2)/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)/rhogeq/Kdrag

    vd1r =(w3r*Kdrag*rhogsol*w1i**3*w2r - w3r*w1r**2*k*Kdrag**2*vgsol - w3r*w1r**3*rhogeq*rhogsol*w2i**2 -&
        w3r*rhogeq*w2r**2*w1r**3*rhogsol - w3r*w1i**2*k*Kdrag**2*vgsol + w3r*rhogeq*w1i**4*w2r*rhogsol +&
        w3r*Kdrag*rhogsol*w1i*w2r*w1r**2 - 2*w3r*rhogeq*w1i*k*Kdrag*vdsol*w2r*w1r +&
        2*w3r*rhogeq*w2r*w1i*k*Kdrag*vgsol*w1r - w3r*rhogeq*w1i**2*w2r**2*rhogsol*w1r + w3r*w1i**2*k*Kdrag**2*vdsol +&
        w3r*rhogeq*w2i*w1r**2*k*Kdrag*vdsol + w3r*w1r**2*k*Kdrag**2*vdsol - w3r*Kdrag*rhogsol*w2i*w1r**3 +&
        w3r**2*w1i**2*rhogeq**2*k*w1r*vgsol + w3r**2*w1r**3*rhogeq**2*k*vgsol - 2*w3r**2*w2i*w1i*rhogeq**2*k*w1r*vgsol&
        - w3r**2*rhogeq*rhogsol*k**2*cs**2*w2r**2 - w3r**2*rhogeq*rhogsol*k**2*cs**2*w2i**2 -&
        w3r**2*rhogeq*w2r*w1r*rhogsol*w1i**2 + w3r**2*w2i*w1i*rhogeq*w1r**2*rhogsol +&
        2*w3r*w2r*w1i**2*rhogeq**2*k*w1r*vgsol + 2*w3r*w2r*w1r**3*rhogeq**2*k*vgsol -&
        w3r*rhogeq*w1i**2*rhogsol*k**2*cs**2*w1r - w3r*cs**2*k**2*rhogeq*rhogsol*w1r**3 -&
        2*w3r*w1i**2*rhogeq**2*w1r**2*k*vgsol + w3r*w1i**3*rhogeq*k*Kdrag*vdsol - 2*w3r*w1i**3*Kdrag*k*rhogeq*vgsol -&
        w3r*w1i**2*rhogeq*w1r*rhogsol*w2i**2 + w3r*w1r**2*rhogeq*w1i*k*Kdrag*vdsol - w3r*w1r*Kdrag*rhogsol*w2i*w1i**2 -&
        w3r*w2i*w1i**2*rhogeq*k*Kdrag*vdsol - 2*w3r*w1r**2*Kdrag*k*rhogeq*vgsol*w1i +&
        2*w3r*rhogeq*w1i**2*w2r*rhogsol*w1r**2 + rhogeq*w1i*rhogsol*k**2*cs**2*w3i**2*w2i +&
        rhogeq*w1i*rhogsol*k**2*cs**2*w3r**2*w2i + rhogeq**2*w1i**2*k*vgsol*w3r*w2r**2 +&
        rhogeq**2*w1i**2*k*vgsol*w3r*w2i**2 - w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r + w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r - rhogeq**2*cs**2*k**3*w2i**2*w1r*vgsol +&
        rhogeq*cs**4*k**4*w2r*w1r*rhogsol - rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol -&
        rhogeq*cs**2*k**3*w2r*w1i*Kdrag*vgsol - rhogeq**2*cs**2*k**3*w2r*w1i**2*vgsol +&
        rhogeq*cs**2*k**3*w2r*w1i*Kdrag*vdsol + rhogeq**2*cs**2*k**3*w2r*w1r**2*vgsol -&
        rhogeq*cs**4*k**4*w2i*w1i*rhogsol + w3i*cs**4*k**4*rhogeq*rhogsol*w2i - w3i*cs**4*k**4*rhogeq*rhogsol*w1i +&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1i + w3i*cs**2*k**2*rhogeq*rhogsol*w1i*w2i**2 +&
        w3r*rhogeq**2*cs**2*k**3*w2i**2*vgsol - w3r*rhogeq*cs**4*k**4*w2r*rhogsol +&
        w3r*rhogeq**2*cs**2*k**3*w2r**2*vgsol - w3r*Kdrag*rhogsol*k**2*cs**2*w1i*w2r +&
        w3r*Kdrag*rhogsol*k**2*cs**2*w2i*w1r - w3r*Kdrag*k*rhogeq*vgsol*w2i*w1r**2 +&
        w3r*Kdrag*k*rhogeq*vgsol*w2i*w1i**2 - w3r*Kdrag**2*k*vdsol*w2r*w1r - w3r*Kdrag**2*k*vdsol*w2i*w1i +&
        w3r*Kdrag**2*k*vgsol*w2r*w1r + w3r*Kdrag**2*k*vgsol*w2i*w1i + w3r*Kdrag*k*rhogeq*vgsol*w2r**2*w1i +&
        w3r*Kdrag*k*rhogeq*vgsol*w1i*w2i**2 + w3r*rhogeq*cs**2*k**2*w2r**2*w1r*rhogsol +&
        w3r*rhogeq*w2i**2*w1r*cs**2*k**2*rhogsol - w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol -&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol + w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol +&
        w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol - w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol +&
        w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol + w3r*rhogeq*cs**4*k**4*w1r*rhogsol -&
        2*w3r*rhogeq**2*cs**2*k**3*w2r*w1r*vgsol + w3r**2*rhogeq**2*cs**2*k**3*w2r*vgsol +&
        w3r**2*rhogeq*w2r*w1i*k*Kdrag*vgsol + w3r**2*rhogeq**2*w2r*w1i**2*k*vgsol - w3r**2*rhogeq**2*w2r*w1r**2*k*vgsol&
        - w3r**2*rhogeq**2*cs**2*k**3*w1r*vgsol - w3r**2*rhogeq*w2i**2*w1i**2*rhogsol +&
        w3r**2*rhogeq*w1r**2*w2i**2*rhogsol + w3r**2*rhogeq*w2r**2*w1r**2*rhogsol - w3r**2*rhogeq*w1i**2*w2r**2*rhogsol&
        - w3r**2*Kdrag*rhogsol*w1i*w2i**2 + w3r**2*Kdrag*rhogsol*w2i*w1r**2 + w3r**2*Kdrag*rhogsol*w2i*w1i**2 -&
        w3r**2*Kdrag*rhogsol*w2r**2*w1i - w3r**2*rhogeq*w2i*w1r*k*Kdrag*vgsol + w3i**2*rhogeq**2*cs**2*k**3*w2r*vgsol -&
        w3i**2*Kdrag*rhogsol*w2r**2*w1i - w3i**2*Kdrag*rhogsol*w1i*w2i**2 + w3i**2*Kdrag*rhogsol*w2i*w1r**2 +&
        w3i**2*Kdrag*rhogsol*w2i*w1i**2 - Kdrag*k*rhogeq*vgsol*w3i*w1r*w2r**2 - Kdrag*k*rhogeq*vgsol*w3i*w2r*w1r**2 +&
        Kdrag*k*rhogeq*vgsol*w3i*w2r*w1i**2 - Kdrag*k*rhogeq*vgsol*w3i*w1r*w2i**2 +&
        Kdrag*rhogsol*k**2*cs**2*w3i*w2r*w1r + Kdrag*rhogsol*k**2*cs**2*w2i*w3i*w1i + Kdrag*rhogsol*w2r**2*w1i**2*w3i +&
        Kdrag*rhogsol*w3i*w1r**2*w2r**2 + Kdrag*rhogsol*w1i**2*w2i**2*w3i + Kdrag*rhogsol*w2i**2*w1r**2*w3i -&
        Kdrag**2*k*vdsol*w2r*w1i*w3i + Kdrag**2*k*vdsol*w2i*w3i*w1r + Kdrag**2*k*vgsol*w2r*w1i*w3i -&
        Kdrag**2*k*vgsol*w2i*w3i*w1r + rhogeq*cs**2*k**3*w2i*w1r*Kdrag*vdsol - rhogeq*cs**2*k**3*w2i*w1r*Kdrag*vgsol -&
        w1r**3*vdsol*Kdrag*k*rhogeq*w2i - Kdrag*rhogsol*k**2*cs**2*w2i*w1r**2 + 2*Kdrag*k*rhogeq*vgsol*w2i*w1r**3 -&
        w3i**2*rhogeq*w2i*w1r*k*Kdrag*vgsol + w3i**2*rhogeq*w2r*w1i*k*Kdrag*vgsol + w3i**2*rhogeq**2*w2r*w1i**2*k*vgsol&
        - w3i**2*rhogeq**2*w2r*w1r**2*k*vgsol - w3i**2*rhogeq**2*cs**2*k**3*w1r*vgsol -&
        w3i**2*rhogeq*w2i**2*w1i**2*rhogsol + w3i**2*rhogeq*w1r**2*w2i**2*rhogsol + w3i**2*rhogeq*w2r**2*w1r**2*rhogsol&
        - w3i**2*rhogeq*w1i**2*w2r**2*rhogsol + rhogeq*w1i**3*w2r**2*rhogsol*w3i - 2*rhogeq*w2r*w1i**3*k*Kdrag*vgsol +&
        rhogeq*w2r**2*w1r**2*rhogsol*w3i*w1i + rhogeq**2*w2r**2*k*w1r*vgsol*w1i**2 + Kdrag**2*k*vdsol*w2r*w1r**2 -&
        Kdrag**2*k*vgsol*w2r*w1r**2 - rhogeq*cs**2*k**2*w2r*w1r**3*rhogsol + Kdrag**2*k*vgsol*w1r**3 -&
        Kdrag**2*k*vdsol*w1r**3 + rhogeq**2*w2r**2*k*w1r**3*vgsol - 2*rhogeq**2*w2r**2*k*w1r*vgsol*w3i*w1i +&
        rhogeq*w3i**2*cs**2*k**2*rhogsol*w2r*w1r - 2*rhogeq*w2r*w1i*k*Kdrag*vgsol*w1r**2 -&
        rhogeq*cs**2*k**2*w1i**2*w2r*rhogsol*w1r - w3r*rhogeq**2*k*vgsol*w2i**2*w1r**2 - rhogeq**2*w2r*w1r**4*k*vgsol -&
        rhogeq*w2r*w3i**2*w1r**3*rhogsol + w3i*rhogeq*k*Kdrag*vdsol*w2r*w1r**2 + rhogeq*w2r*w1i*k*Kdrag*vdsol*w1r**2 -&
        Kdrag*rhogsol*w3i*w2r*w1r**3 - 2*rhogeq**2*w2r*w1i**2*k*vgsol*w1r**2 - rhogeq*w3i**2*cs**2*k**2*rhogsol*w2r**2&
        - rhogeq*w3i**2*cs**2*k**2*rhogsol*w2i**2 - Kdrag*rhogsol*w3i*w2r*w1r*w1i**2 -&
        w3i*rhogeq*k*Kdrag*vdsol*w2r*w1i**2 + Kdrag**2*k*vdsol*w2r*w1i**2 - rhogeq*w2r*w3i**2*w1r*rhogsol*w1i**2 -&
        Kdrag**2*k*vgsol*w2r*w1i**2 + rhogeq*w2r*w1i**3*k*Kdrag*vdsol - rhogeq**2*w2r*w1i**4*k*vgsol +&
        cs**2*k**2*rhogeq*rhogsol*w1r**4 + w1r**2*w3i*rhogeq*rhogsol*w2i**2*w1i +&
        2*rhogeq*w1i**2*rhogsol*k**2*cs**2*w1r**2 + rhogeq*w1i**4*rhogsol*k**2*cs**2 -&
        w2i*w1i**3*rhogeq*rhogsol*k**2*cs**2 - w1r**4*w3i*rhogeq*rhogsol*w2i - w1r**2*rhogeq*cs**4*k**4*rhogsol +&
        w2i*w1i**3*rhogeq*rhogsol*w3i**2 - w3i*cs**2*k**2*rhogeq*rhogsol*w1i*w1r**2 +&
        w2i*w1i*rhogeq*w1r**2*rhogsol*w3i**2 - w3i*cs**2*k**2*rhogeq*rhogsol*w1i**3 + w1i**3*w3i*rhogeq*rhogsol*w2i**2&
        - w2i*w1i*cs**2*k**2*rhogeq*rhogsol*w1r**2 + w1i**2*rhogeq**2*k*vgsol*w2i**2*w1r -&
        2*w2i*w1i*rhogeq**2*k*w1r*vgsol*w3i**2 - 2*vgsol*k*rhogeq*Kdrag*w2i*w1r*w3i*w1i +&
        w1r**3*rhogeq**2*k*vgsol*w2i**2 - w1i**4*w3i*rhogeq*rhogsol*w2i + 2*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w1r +&
        rhogeq*cs**4*k**4*rhogsol*w1i**2 - 2*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w1r + Kdrag*rhogsol*k**2*cs**2*w1i**3 -&
        Kdrag**2*k*vdsol*w1r*w1i**2 + rhogeq**2*w3i**2*k*w1r**3*vgsol + Kdrag**2*k*vgsol*w1r*w1i**2 -&
        Kdrag*rhogsol*w1i**3*w3i*w2i + 2*w1r*Kdrag*k*rhogeq*vgsol*w2i*w1i**2 + 2*w2i*w1i*rhogeq**2*cs**2*k**3*w1r*vgsol&
        - Kdrag*rhogsol*w3i*w1r**2*w2i*w1i + rhogeq**2*w3i**2*k*w1r*vgsol*w1i**2 +&
        2*w3i*rhogeq*k*Kdrag*vdsol*w1r*w2i*w1i - Kdrag*rhogsol*k**2*cs**2*w3i*w1r**2 - w3i*rhogeq*k*Kdrag*vdsol*w1r**3&
        + 2*rhogeq**2*cs**2*k**3*w1r*vgsol*w3i*w1i + 2*Kdrag*k*rhogeq*vgsol*w3i*w1r**3 -&
        Kdrag*rhogsol*k**2*cs**2*w1i**2*w2i - Kdrag*rhogsol*k**2*cs**2*w3i*w1i**2 + Kdrag*rhogsol*k**2*cs**2*w1i*w1r**2&
        - 2*rhogeq**2*k*vgsol*w2i**2*w1r*w3i*w1i - 2*w1r**2*w3i*rhogeq*rhogsol*w1i**2*w2i -&
        w1r*w1i**2*vdsol*Kdrag*k*rhogeq*w2i - 2*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i*w3i -&
        w3i*rhogeq*k*Kdrag*vdsol*w1r*w1i**2 + 2*Kdrag*k*rhogeq*vgsol*w3i*w1r*w1i**2 - w3r**2*rhogeq*w2r*w1r**3*rhogsol&
        - w3r*w1i**4*rhogeq**2*k*vgsol + w3r*rhogeq*w2r*w1r**4*rhogsol - w3r*rhogeq**2*w1r**4*k*vgsol +&
        2*w1i**2*rhogeq**2*k*w1r*vgsol*w2i*w3i + w3r**2*w2i*w1i**3*rhogeq*rhogsol - w3r*rhogeq**2*w1r**2*k*vgsol*w2r**2&
        + 2*w1r**3*rhogeq**2*k*vgsol*w2i*w3i + w3r**2*rhogeq*rhogsol*k**2*cs**2*w2r*w1r)/(w1i**2 - 2*w3i*w1i + w3r**2 +&
        w1r**2 + w3i**2 - 2*w3r*w1r)/Kdrag/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)/rhogeq/k

    vd1i = - (w3r**2*rhogeq*rhogsol*w1i**3*w2r + w3r*Kdrag*rhogsol*k**2*cs**2*w1i**2 +&
        w3r*rhogeq*rhogsol*k**2*cs**2*w1i**3 - 2*w3r*rhogeq*rhogsol*w2i*w1i**2*w1r**2 +&
        w3r*rhogeq*rhogsol*w1i*w2i**2*w1r**2 + w3r*rhogeq*rhogsol*w2r**2*w1i*w1r**2 -&
        2*w3r*w2r*rhogeq**2*k*vgsol*w1i*w1r**2 - 2*w3r*w2r*rhogeq**2*cs**2*k**3*vgsol*w1i -&
        2*w3r*w2r*rhogeq**2*k*vgsol*w1i**3 - w3r*Kdrag*rhogsol*w2r*w1i**2*w1r +&
        w3r*rhogeq*rhogsol*k**2*cs**2*w1i*w1r**2 + w3r*rhogeq*k*Kdrag*vdsol*w2r*w1i**2 +&
        w3r*Kdrag*rhogsol*k**2*cs**2*w1r**2 - w3r*Kdrag*w2i*rhogsol*w1i*w1r**2 - w3r*rhogeq*k*Kdrag*vdsol*w2r*w1r**2 +&
        w3r*rhogeq*k*Kdrag*vdsol*w1r**3 - 2*w3r*rhogeq*k*Kdrag*vdsol*w1r*w2i*w1i + w3r*rhogeq*k*Kdrag*vdsol*w1r*w1i**2&
        + 2*w3r*vgsol*k*rhogeq**2*w2i**2*w1i*w1r + 2*w3r*Kdrag*k*rhogeq*vgsol*w1i*w2i*w1r +&
        2*w3r*rhogeq**2*cs**2*k**3*vgsol*w1i*w1r - 2*w3r**2*rhogeq*rhogsol*w1i*w2r**2*w1r -&
        w3r**2*rhogeq**2*k*vgsol*w1i*w1r**2 - 2*w3r**2*rhogeq*rhogsol*w2i**2*w1r*w1i +&
        2*w3r**2*rhogeq**2*k*vgsol*w1i*w2r*w1r + w3r**2*rhogeq*rhogsol*w1i*w2r*w1r**2 +&
        rhogeq**2*w3i**2*w2i*w1i**2*k*vgsol - rhogeq**2*w3i**2*w2i*k*vgsol*w1r**2 - Kdrag*w3i**2*k*rhogeq*w1r**2*vgsol&
        - Kdrag*w1i**2*rhogeq*k*vgsol*w3i**2 - Kdrag*w3r**2*w1i**2*k*rhogeq*vgsol - Kdrag*w3r**2*w1r**2*rhogeq*k*vgsol&
        - rhogeq**2*w3r**2*k*vgsol*w2i*w1r**2 + rhogeq**2*w3r**2*k*vgsol*w2i*w1i**2 + w3r*Kdrag*rhogsol*w2r**2*w1r**2 +&
        w3r*Kdrag*rhogsol*w2r**2*w1i**2 - 3*w3r*Kdrag*k*rhogeq*vgsol*w2r*w1i**2 - w3r*Kdrag*k*rhogeq*vgsol*w2r*w1r**2 +&
        w3r*rhogeq*cs**2*k**2*rhogsol*w2r**2*w1i - w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r -&
        w3r*Kdrag*cs**2*k**2*rhogsol*w2r*w1r + w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r + w3r**2*Kdrag*rhogsol*w2r*w1r**2&
        + w3r**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i + w3r**2*Kdrag*rhogeq*k*vgsol*w2r*w1r -&
        w3r**2*rhogeq**2*cs**2*k**3*vgsol*w1i - w3r*Kdrag**2*k*vgsol*w2r*w1i -&
        4*w3r*cs**2*k**2*rhogeq*rhogsol*w2r*w1i*w1r + 2*w3r*rhogeq**2*k*vgsol*w2r**2*w1i*w1r -&
        Kdrag*rhogsol*w2r*w1i**3*w3i + cs**2*k**2*rhogeq*rhogsol*w2r*w1i**3 + 2*Kdrag*rhogeq*k*vgsol*w2r*w1r*w3i*w1i -&
        rhogeq**2*k*vgsol*w2r**2*w1i**3 - 2*rhogeq*k*Kdrag*vdsol*w1r**2*w1i**2 - rhogeq*k*Kdrag*vdsol*w1r**4 -&
        Kdrag*rhogsol*w2r*w1r**2*w3i*w1i + Kdrag*cs**2*k**2*rhogsol*w2r*w1r**2 +&
        cs**2*k**2*rhogeq*rhogsol*w2r*w1i*w1r**2 - rhogeq**2*k*vgsol*w2r**2*w1i*w1r**2 -&
        vgsol*k*rhogeq**2*w2i**2*w1i**3 - rhogeq*k*Kdrag*vdsol*w1i**4 - 2*rhogeq**2*k*vgsol*w1i**3*w2i*w3i -&
        cs**2*k**3*rhogeq*Kdrag*vgsol*w1r**2 + 4*rhogeq*cs**2*k**2*rhogsol*w2i*w1r*w3i*w1i -&
        rhogeq*cs**2*k**2*rhogsol*w2i*w1r*w1i**2 + w3r**2*rhogeq*rhogsol*w2i*w1r*w1i**2 -&
        rhogeq*cs**2*k**2*rhogsol*w2i*w1r**3 + cs**2*k**3*rhogeq*Kdrag*vdsol*w1r**2 -&
        2*cs**4*k**4*rhogeq*rhogsol*w1i*w1r + w1i**3*k*Kdrag**2*vgsol - vgsol*k*rhogeq**2*w2i**2*w1i*w1r**2 -&
        Kdrag**2*k*vgsol*w2i*w1r**2 + rhogeq*k*Kdrag*vdsol*w2r*w1r**3 + rhogeq*k*Kdrag*vdsol*w2i*w1i*w1r**2 +&
        cs**2*k**3*rhogeq*Kdrag*vgsol*w1i**2 - cs**2*k**3*rhogeq*Kdrag*vdsol*w1i**2 + w1i*k*Kdrag**2*vgsol*w1r**2 +&
        rhogeq*w3i*w1i*k*Kdrag*vdsol*w1r**2 + 2*rhogeq**2*w3i*w1i**2*k*vgsol*w1r**2 -&
        rhogeq**2*w3i**2*w1i*k*vgsol*w1r**2 - rhogeq*w3i*w1r**3*cs**2*k**2*rhogsol +&
        rhogeq*w3i*w2i*k*Kdrag*vdsol*w1r**2 + rhogeq*w3i**2*w2r*rhogsol*w1i*w1r**2 + rhogeq*w3i**2*w2i*w1r**3*rhogsol +&
        rhogeq*k*Kdrag*vdsol*w2r*w1r*w1i**2 + rhogeq*w3i*w2i**2*w1r*rhogsol*w1i**2 -&
        Kdrag*w1r*cs**2*k**2*rhogsol*w1i**2 + Kdrag*w2i*w3i*w1r*rhogsol*w1i**2 + Kdrag*w2r*cs**2*k**2*rhogsol*w1i**2 +&
        rhogeq*w3i**2*w2i*w1r*rhogsol*w1i**2 + rhogeq*w3i*w1i**2*w2r**2*rhogsol*w1r +&
        2*rhogeq**2*w3i**2*w1i*k*vgsol*w2r*w1r - 2*rhogeq*w3i*w1i*k*Kdrag*vdsol*w2r*w1r +&
        rhogeq*w3i*w1i**3*k*Kdrag*vdsol - rhogeq*w3i*w1i**2*k*Kdrag*vdsol*w2i -&
        rhogeq*w3i*w1r*cs**2*k**2*rhogsol*w1i**2 - 2*rhogeq*w3i**2*w2r**2*rhogsol*w1i*w1r +&
        rhogeq*w3i**2*w2r*rhogsol*w1i**3 - 2*rhogeq*w3i**2*w2i**2*w1r*rhogsol*w1i -&
        2*rhogeq*w3i*w1i**2*w2r*rhogsol*w1r**2 + 2*Kdrag*w1i**2*k*rhogeq*vgsol*w1r**2 - w1i**3*k*Kdrag**2*vdsol +&
        Kdrag*w1r**4*rhogeq*k*vgsol + rhogeq*k*Kdrag*vdsol*w2i*w1i**3 + 2*rhogeq**2*k*vgsol*w2i*w1i**2*w1r**2 +&
        2*rhogeq**2*cs**2*k**3*vgsol*w1i*w2r*w1r + Kdrag**2*k*vdsol*w2i*w1r**2 - 2*rhogeq**2*k*vgsol*w1i*w1r**2*w2i*w3i&
        + w3i*k*Kdrag**2*vdsol*w1r**2 + rhogeq*w3i*w2i**2*w1r**3*rhogsol - w3i*k*Kdrag**2*vgsol*w1r**2 +&
        rhogeq**2*w3i*w1r**4*k*vgsol - w1i*k*Kdrag**2*vdsol*w1r**2 - w3i*rhogeq**2*k*vgsol*w2r**2*w1r**2 -&
        rhogeq*cs**2*k**3*vdsol*Kdrag*w2r*w1r + w3i*Kdrag*k*rhogeq*vgsol*w2r**2*w1i -&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1r + 2*w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w1r**2 -&
        2*w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w1i**2 + w3i*rhogeq**2*k*vgsol*w2r**2*w1i**2 +&
        w3i**2*Kdrag*rhogeq*k*vgsol*w2r*w1r - w1i**2*k*Kdrag**2*vgsol*w2i + w3i*rhogeq**2*cs**2*k**3*vgsol*w2r**2 -&
        Kdrag*k*rhogeq*vgsol*w2r**2*w1r**2 - Kdrag*k*rhogeq*vgsol*w2r**2*w1i**2 +&
        w3i**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i + rhogeq*cs**4*k**4*rhogsol*w2r*w1i -&
        w3i*rhogeq*cs**4*k**4*rhogsol*w2r - w3i*Kdrag*cs**2*k**2*rhogsol*w2r*w1i + w1i**2*k*Kdrag**2*vdsol*w2i -&
        rhogeq**2*w3i**2*w1i**3*k*vgsol + rhogeq**2*w3i*w1i**4*k*vgsol - rhogeq*w3i*w2r*w1r**4*rhogsol -&
        Kdrag*w1r**3*cs**2*k**2*rhogsol + Kdrag*w2i*w3i*w1r**3*rhogsol + rhogeq*w3i*w2r**2*w1r**3*rhogsol -&
        w3i*k*Kdrag**2*vgsol*w1i**2 + w3i*k*Kdrag**2*vdsol*w1i**2 + rhogeq**2*k*vgsol*w2i*w1r**4 -&
        rhogeq*w3i*w1i**4*w2r*rhogsol + Kdrag*w1i**4*k*rhogeq*vgsol + cs**2*k**3*rhogeq*Kdrag*vgsol*w2r*w1r -&
        w3i**2*cs**2*k**3*rhogeq**2*vgsol*w1i - w3i**2*Kdrag*rhogsol*w2r**2*w1r - w3i*Kdrag**2*k*vdsol*w2r*w1r +&
        w3i**2*Kdrag*rhogsol*w2r*w1r**2 + w3i**2*Kdrag*rhogsol*w2r*w1i**2 + w3i*Kdrag**2*k*vgsol*w2r*w1r +&
        rhogeq**2*k*vgsol*w2i*w1i**4 - w3r*Kdrag*rhogsol*w2r*w1r**3 + w3r*rhogeq*rhogsol*w2r**2*w1i**3 -&
        w3r**2*rhogeq**2*k*vgsol*w1i**3 + w3r**2*rhogeq*rhogsol*w2i*w1r**3 + w3r**2*Kdrag*rhogsol*w2r*w1i**2 -&
        w3r**2*Kdrag*rhogsol*w2r**2*w1r - w3r*rhogeq*rhogsol*w2i*w1i**4 + w3r*rhogeq*rhogsol*w1i**3*w2i**2 -&
        w3r*Kdrag*w2i*rhogsol*w1i**3 - w3r*rhogeq*rhogsol*w2i*w1r**4 - rhogeq**2*cs**2*k**3*vgsol*w1i*w2r**2 +&
        w3r*Kdrag*k*rhogeq*vgsol*w2r**2*w1r + w3r*Kdrag**2*k*vdsol*w2r*w1i + 2*w3r*rhogeq*cs**2*k**2*rhogsol*w2i*w1r**2&
        + w3r*cs**4*k**4*rhogeq*rhogsol*w1i - w3r*Kdrag*cs**2*k**2*rhogsol*w1i*w2i -&
        w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r - w3i*Kdrag**2*k*vdsol*w1i*w2i + w3i*Kdrag**2*k*vgsol*w1i*w2i +&
        w3r*Kdrag**2*k*vgsol*w2i*w1r + cs**2*k**3*rhogeq*Kdrag*vdsol*w1i*w2i + w3i**2*Kdrag*rhogeq*k*vgsol*w1i*w2i -&
        Kdrag*rhogeq*k*vgsol*w2i**2*w1r**2 - Kdrag*rhogeq*k*vgsol*w2i**2*w1i**2 + w3r*Kdrag*k*rhogeq*vgsol*w2i**2*w1r +&
        w3i*rhogeq*k*Kdrag*vgsol*w2i**2*w1i - 3*w3i*rhogeq*k*Kdrag*vgsol*w2i*w1r**2 -&
        w3i*rhogeq*k*Kdrag*vgsol*w1i**2*w2i + w3i*rhogeq**2*k*vgsol*w2i**2*w1i**2 -&
        w3i*rhogeq**2*cs**2*k**3*vgsol*w1r**2 + w3i*rhogeq**2*cs**2*k**3*vgsol*w1i**2 +&
        w3r**2*Kdrag*k*rhogeq*vgsol*w1i*w2i + w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r +&
        w3i*Kdrag*cs**2*k**2*rhogsol*w2i*w1r + rhogeq*cs**4*k**4*rhogsol*w2i*w1r -&
        w3r**2*rhogeq*cs**2*k**2*rhogsol*w2i*w1r + w3r**2*vgsol*rhogeq**2*k**3*cs**2*w2i -&
        w3i*rhogeq*cs**2*k**3*vdsol*Kdrag*w2i + w3i*rhogeq*cs**2*k**3*vdsol*Kdrag*w1i -&
        cs**2*k**3*rhogeq*Kdrag*vgsol*w1i*w2i - w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1i +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i + w3i*rhogeq*cs**4*k**4*rhogsol*w1r +&
        w3i**2*vgsol*rhogeq**2*k**3*cs**2*w2i - w3i**2*rhogeq*cs**2*k**2*rhogsol*w2i*w1r +&
        w3i*rhogeq**2*cs**2*k**3*vgsol*w2i**2 - 2*w3i*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i -&
        w3i*rhogeq**2*k*vgsol*w2i**2*w1r**2 + w3r*Kdrag*rhogsol*w2i**2*w1r**2 + rhogeq**2*cs**2*k**3*vgsol*w1i**2*w2i -&
        rhogeq**2*cs**2*k**3*vgsol*w2i*w1r**2 - w3r*Kdrag**2*k*vdsol*w2i*w1r - w3r*cs**4*k**4*rhogeq*rhogsol*w2i -&
        2*w3r*rhogeq*cs**2*k**2*rhogsol*w1i**2*w2i + w3r*rhogeq*cs**2*k**2*rhogsol*w2i**2*w1i +&
        w3r*Kdrag*rhogsol*w2i**2*w1i**2 - w3i*rhogeq*cs**2*k**2*rhogsol*w2i**2*w1r - w3i**2*Kdrag*rhogsol*w2i**2*w1r -&
        w3r**2*Kdrag*rhogsol*w2i**2*w1r - rhogeq**2*cs**2*k**3*vgsol*w1i*w2i**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2&
        + w3i**2 - 2*w3r*w1r)/Kdrag/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)/rhogeq/k
 endif

!-------------------------------
! G A S  D E N S I T I E S
!-------------------------------
 rhog3r =(w3r**2*k*rhogeq*vgsol*w1i + w3r**2*w2i*rhogeq*k*vgsol - w3r**2*w2r*w1i*rhogsol - w3i**2*w2i*rhogeq*k*vgsol +&
        w3i**2*w2r*w1i*rhogsol - w3i*w1i**2*w2r*rhogsol - w3i*w2i**2*w1r*rhogsol - w3i*w2r**2*w1r*rhogsol -&
        w3i*w2r*w1r**2*rhogsol + w3i**2*rhogsol*w2i*w1r - w3i*w1i*k*Kdrag*vdsol + k*Kdrag*vdsol*w2i*w1i -&
        k*Kdrag*vgsol*w2i*w1i + w3i*w1i**2*k*rhogeq*vgsol + w3i*w1r**2*rhogeq*k*vgsol - w3i*w1r*cs**2*k**2*rhogsol -&
        w3i*w2i*k*Kdrag*vdsol + w3i*w1i*k*Kdrag*vgsol + w3i*w2i*k*Kdrag*vgsol + w3i*w2r**2*rhogeq*k*vgsol -&
        w3i**2*k*Kdrag*vgsol - k*Kdrag*vdsol*w2r*w1r + w3i**2*k*Kdrag*vdsol - k*rhogeq*vgsol*w2r**2*w1i -&
        w3i*w2r*cs**2*k**2*rhogsol + w3i*w2i**2*rhogeq*k*vgsol + k*Kdrag*vgsol*w2r*w1r - k*rhogeq*vgsol*w2i*w1r**2 +&
        rhogsol*k**2*cs**2*w2i*w1r + rhogsol*k**2*cs**2*w1i*w2r - k*rhogeq*vgsol*w1i*w2i**2 - w3i**2*k*rhogeq*vgsol*w1i&
        - k*rhogeq*vgsol*w2i*w1i**2 + w3r*k*Kdrag*vdsol*w2r + 2*w3r*w3i*rhogsol*k**2*cs**2 - w3r*rhogsol*k**2*cs**2*w1i&
        + w3r*k*Kdrag*vdsol*w1r - w3r*rhogsol*k**2*cs**2*w2i + 2*w3r*w3i*w2r*w1r*rhogsol - 2*w3r*w3i*w2i*rhogsol*w1i -&
        w3r*k*Kdrag*vgsol*w2r - w3r*k*Kdrag*vgsol*w1r + w3r*rhogsol*w2i*w1i**2 + w3r*rhogsol*w2i*w1r**2 +&
        w3r*rhogsol*w2r**2*w1i + w3r*rhogsol*w1i*w2i**2 + w3r**2*k*Kdrag*vgsol - w3r**2*k*Kdrag*vdsol -&
        w3r**2*rhogsol*w2i*w1r + 2*w3i*w2r*k*rhogeq*w1r*vgsol + 2*w3i*k*w2i*rhogeq*w1i*vgsol -&
        2*w3r*w3i*w2r*vgsol*k*rhogeq - 2*w3r*w3i*k*rhogeq*w1r*vgsol)/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i&
        + w3r**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)

 rhog3i = - ( - rhogsol*w1i**2*w2r**2 - rhogsol*w1i**2*w2i**2 - rhogsol*w1r**2*w2r**2 - rhogsol*w1r**2*w2i**2 -&
        2*w3r*rhogeq*k*vgsol*w2i*w1i + rhogsol*w2i**2*w3i*w1i + rhogsol*k**2*cs**2*w1r*w3r - rhogsol*k**2*cs**2*w3i*w1i&
        - vgsol*k*rhogeq*w3r*w2i**2 - vgsol*Kdrag*k*w3r*w1i - vgsol*Kdrag*k*w3r*w2i + 2*w3r*rhogeq*k*vgsol*w3i*w1i +&
        rhogsol*w1i**2*w2i*w3i + rhogsol*w1i*w2i*w3r**2 + rhogsol*w2i*w3i*w1r**2 - rhogsol*w1i*w2i*w3i**2 +&
        rhogsol*w2i**2*w3r*w1r - rhogsol*k**2*cs**2*w2i*w3i - vgsol*k*rhogeq*w3r*w1i**2 + w3r*rhogsol*k**2*cs**2*w2r -&
        2*w3r*w2r*k*rhogeq*w1r*vgsol + w3r**2*w2r*vgsol*k*rhogeq - w3r*w2r**2*rhogeq*k*vgsol -&
        vgsol*k*rhogeq*w3r*w1r**2 + vgsol*k*rhogeq*w3r**2*w1r + vdsol*Kdrag*k*w3r*w1i - vdsol*Kdrag*k*w2i*w1r +&
        vdsol*Kdrag*k*w3r*w2i + vdsol*Kdrag*k*w3i*w1r + rhogsol*k**2*cs**2*w1i*w2i + vgsol*Kdrag*k*w2i*w1r -&
        vgsol*Kdrag*k*w3i*w1r - vgsol*k*rhogeq*w1r*w3i**2 + vgsol*k*rhogeq*w2i**2*w1r - rhogsol*k**2*cs**2*w2r*w1r -&
        w3r**2*w2r*w1r*rhogsol + w3r*w2r*rhogsol*w1i**2 + w3r*w2r*w1r**2*rhogsol + w3r*w2r**2*w1r*rhogsol -&
        2*w3r*w1r*rhogsol*w2i*w3i + 2*vgsol*Kdrag*k*w3i*w3r - 2*vdsol*Kdrag*k*w3i*w3r - w3r**2*rhogsol*k**2*cs**2 +&
        2*w3r*rhogeq*k*vgsol*w2i*w3i + rhogsol*k**2*cs**2*w3i**2 - 2*w2r*rhogsol*w3i*w1i*w3r +&
        vgsol*k*rhogeq*w2r*w1r**2 + vgsol*k*rhogeq*w1r*w2r**2 + vgsol*k*rhogeq*w2r*w1i**2 - vgsol*k*rhogeq*w2r*w3i**2 +&
        vgsol*Kdrag*k*w1i*w2r - vgsol*Kdrag*k*w2r*w3i + rhogsol*w2r**2*w3i*w1i + rhogsol*w2r*w1r*w3i**2 +&
        vdsol*Kdrag*k*w2r*w3i - vdsol*Kdrag*k*w1i*w2r)/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i +&
        w3r**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)

 rhog2r =( - w3r**2*k*rhogeq*vgsol*w1i + w3r**2*w2i*rhogeq*k*vgsol + w3r**2*w2r*w1i*rhogsol +&
        w3i**2*w2i*rhogeq*k*vgsol + w3i**2*w2r*w1i*rhogsol + w3i*w1i**2*w2r*rhogsol + w3i*w2i**2*w1r*rhogsol -&
        w3i*w2r**2*w1r*rhogsol + w3i*w2r*w1r**2*rhogsol + 2*rhogsol*w3r*w2i*w1r*w2r + 2*cs**2*k**2*rhogsol*w2r*w2i -&
        k*Kdrag*vgsol*w2i**2 + k*Kdrag*vgsol*w2r**2 - 2*rhogsol*w2r*w2i*w3i*w1i + k*Kdrag*vdsol*w2i**2 -&
        k*Kdrag*vdsol*w2r**2 - 2*vgsol*rhogeq*k*w3r*w2r*w2i - 2*vgsol*rhogeq*k*w2i*w1r*w2r +&
        2*vgsol*rhogeq*k*w3r*w2i*w1r - w3i**2*rhogsol*w2i*w1r + w3i*w1i*k*Kdrag*vdsol - k*Kdrag*vdsol*w2i*w1i +&
        k*Kdrag*vgsol*w2i*w1i - w3i*w1i**2*k*rhogeq*vgsol - w3i*w1r**2*rhogeq*k*vgsol + w3i*w1r*cs**2*k**2*rhogsol -&
        w3i*w2i*k*Kdrag*vdsol - w3i*w1i*k*Kdrag*vgsol + w3i*w2i*k*Kdrag*vgsol + w3i*w2r**2*rhogeq*k*vgsol +&
        k*Kdrag*vdsol*w2r*w1r + k*rhogeq*vgsol*w2r**2*w1i - w3i*w2r*cs**2*k**2*rhogsol - w3i*w2i**2*rhogeq*k*vgsol -&
        k*Kdrag*vgsol*w2r*w1r + k*rhogeq*vgsol*w2i*w1r**2 - rhogsol*k**2*cs**2*w2i*w1r - rhogsol*k**2*cs**2*w1i*w2r -&
        k*rhogeq*vgsol*w1i*w2i**2 - w3i**2*k*rhogeq*vgsol*w1i + k*rhogeq*vgsol*w2i*w1i**2 + w3r*k*Kdrag*vdsol*w2r +&
        w3r*rhogsol*k**2*cs**2*w1i - w3r*k*Kdrag*vdsol*w1r - w3r*rhogsol*k**2*cs**2*w2i - w3r*k*Kdrag*vgsol*w2r +&
        w3r*k*Kdrag*vgsol*w1r - w3r*rhogsol*w2i*w1i**2 - w3r*rhogsol*w2i*w1r**2 - w3r*rhogsol*w2r**2*w1i +&
        w3r*rhogsol*w1i*w2i**2 - w3r**2*rhogsol*w2i*w1r + 2*w3i*k*w2i*rhogeq*w1i*vgsol)/(w2r**2 - 2*w3r*w2r + w2i**2 +&
        w3i**2 - 2*w2i*w3i + w3r**2)/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)

 rhog2i =( - 2*w2i*w2r*vgsol*k*rhogeq*w1i + 2*w3r*w2i*rhogsol*w1i*w2r - 2*w2i*k*rhogeq*vgsol*w2r*w3i +&
        rhogsol*w2i**2*w3i*w1i + rhogsol*k**2*cs**2*w1r*w3r - rhogsol*k**2*cs**2*w3i*w1i + vgsol*k*rhogeq*w3r*w2i**2 -&
        vgsol*Kdrag*k*w3r*w1i + vgsol*Kdrag*k*w3r*w2i - rhogsol*w1i**2*w2i*w3i - rhogsol*w1i*w2i*w3r**2 -&
        rhogsol*w2i*w3i*w1r**2 - rhogsol*w1i*w2i*w3i**2 - rhogsol*w2i**2*w3r*w1r + rhogsol*k**2*cs**2*w2i*w3i -&
        vgsol*k*rhogeq*w3r*w1i**2 - w3r*rhogsol*k**2*cs**2*w2r + 2*w3r*w2r*k*rhogeq*w1r*vgsol +&
        w3r**2*w2r*vgsol*k*rhogeq - w3r*w2r**2*rhogeq*k*vgsol - vgsol*k*rhogeq*w3r*w1r**2 - vgsol*k*rhogeq*w3r**2*w1r +&
        vdsol*Kdrag*k*w3r*w1i - vdsol*Kdrag*k*w2i*w1r - vdsol*Kdrag*k*w3r*w2i + vdsol*Kdrag*k*w3i*w1r +&
        rhogsol*k**2*cs**2*w1i*w2i + vgsol*Kdrag*k*w2i*w1r - vgsol*Kdrag*k*w3i*w1r - vgsol*k*rhogeq*w1r*w3i**2 +&
        vgsol*k*rhogeq*w2i**2*w1r - rhogsol*k**2*cs**2*w2r*w1r - w3r**2*w2r*w1r*rhogsol - w3r*w2r*rhogsol*w1i**2 -&
        w3r*w2r*w1r**2*rhogsol + w3r*w2r**2*w1r*rhogsol + rhogsol*w3r**2*w1i**2 + rhogsol*w3r**2*w1r**2 +&
        rhogsol*w1i**2*w3i**2 + vgsol*k*rhogeq*w2r*w1r**2 - vgsol*k*rhogeq*w1r*w2r**2 + vgsol*k*rhogeq*w2r*w1i**2 +&
        rhogsol*w3i**2*w1r**2 - w2i**2*rhogsol*k**2*cs**2 + vgsol*k*rhogeq*w2r*w3i**2 + vgsol*Kdrag*k*w1i*w2r +&
        vgsol*Kdrag*k*w2r*w3i + w2r**2*rhogsol*k**2*cs**2 + 2*w2i*k*Kdrag*vdsol*w2r + 2*w2i*rhogsol*w3i*w2r*w1r -&
        2*w2i*k*Kdrag*vgsol*w2r + 2*w2r*w3i*k*rhogeq*w1i*vgsol - rhogsol*w2r**2*w3i*w1i - rhogsol*w2r*w1r*w3i**2 -&
        vdsol*Kdrag*k*w2r*w3i - vdsol*Kdrag*k*w1i*w2r)/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i +&
        w3r**2)/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)

 rhog1r = - ( - w3r**2*k*rhogeq*vgsol*w1i - 2*rhogsol*k**2*cs**2*w1i*w1r + 2*k*rhogeq*vgsol*w1i*w2r*w1r +&
        2*rhogsol*w2i*w1r*w3i*w1i + 2*w3r*k*rhogeq*vgsol*w1i*w1r - 2*w3r*w2r*k*rhogeq*vgsol*w1i -&
        2*w3r*w2r*w1i*rhogsol*w1r + w3r**2*w2i*rhogeq*k*vgsol + w3r**2*w2r*w1i*rhogsol + w3i**2*w2i*rhogeq*k*vgsol +&
        w3i**2*w2r*w1i*rhogsol - w3i*w1i**2*w2r*rhogsol - w3i*w2i**2*w1r*rhogsol - w3i*w2r**2*w1r*rhogsol +&
        w3i*w2r*w1r**2*rhogsol - w3i**2*rhogsol*w2i*w1r + w3i*w1i*k*Kdrag*vdsol + k*Kdrag*vdsol*w2i*w1i -&
        k*Kdrag*vgsol*w2i*w1i + w3i*w1i**2*k*rhogeq*vgsol - w3i*w1r**2*rhogeq*k*vgsol + w3i*w1r*cs**2*k**2*rhogsol -&
        w3i*w2i*k*Kdrag*vdsol - w3i*w1i*k*Kdrag*vgsol + w3i*w2i*k*Kdrag*vgsol + w3i*w2r**2*rhogeq*k*vgsol -&
        k*Kdrag*vdsol*w2r*w1r - k*rhogeq*vgsol*w2r**2*w1i - w3i*w2r*cs**2*k**2*rhogsol + w3i*w2i**2*rhogeq*k*vgsol +&
        k*Kdrag*vgsol*w2r*w1r - k*rhogeq*vgsol*w2i*w1r**2 + rhogsol*k**2*cs**2*w2i*w1r + rhogsol*k**2*cs**2*w1i*w2r -&
        k*rhogeq*vgsol*w1i*w2i**2 - w3i**2*k*rhogeq*vgsol*w1i + k*rhogeq*vgsol*w2i*w1i**2 + w3r*k*Kdrag*vdsol*w2r +&
        w3r*rhogsol*k**2*cs**2*w1i - w3r*k*Kdrag*vdsol*w1r - w3r*rhogsol*k**2*cs**2*w2i - w3r*k*Kdrag*vgsol*w2r +&
        w3r*k*Kdrag*vgsol*w1r - Kdrag*w1r**2*k*vgsol - w3r*rhogsol*w2i*w1i**2 + w3r*rhogsol*w2i*w1r**2 +&
        Kdrag*w1i**2*k*vgsol + w3r*rhogsol*w2r**2*w1i + w3r*rhogsol*w1i*w2i**2 + k*Kdrag*vdsol*w1r**2 -&
        w3r**2*rhogsol*w2i*w1r - k*Kdrag*vdsol*w1i**2 - 2*w3i*k*w2i*rhogeq*w1i*vgsol)/(w1i**2 - 2*w3i*w1i + w3r**2 +&
        w1r**2 + w3i**2 - 2*w3r*w1r)/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)

 rhog1i = - ( - 2*vdsol*Kdrag*k*w1i*w1r + rhogsol*w2i**2*w3i*w1i + rhogsol*k**2*cs**2*w1i**2 -&
        rhogsol*k**2*cs**2*w1r**2 + rhogsol*k**2*cs**2*w1r*w3r - rhogsol*k**2*cs**2*w3i*w1i + vgsol*k*rhogeq*w3r*w2i**2&
        + 2*vgsol*Kdrag*k*w1i*w1r - vgsol*Kdrag*k*w3r*w1i + vgsol*Kdrag*k*w3r*w2i - rhogsol*w1i**2*w2i*w3i +&
        rhogsol*w1i*w2i*w3r**2 - 2*rhogsol*w1i*w2i*w3r*w1r + rhogsol*w2i*w3i*w1r**2 + rhogsol*w1i*w2i*w3i**2 +&
        rhogsol*w2i**2*w3r*w1r + rhogsol*k**2*cs**2*w2i*w3i - vgsol*k*rhogeq*w3r*w1i**2 - w3r*rhogsol*k**2*cs**2*w2r -&
        2*w3r*w2r*k*rhogeq*w1r*vgsol + w3r**2*w2r*vgsol*k*rhogeq + w3r*w2r**2*rhogeq*k*vgsol +&
        vgsol*k*rhogeq*w3r*w1r**2 + 2*vgsol*k*rhogeq*w1r*w2i*w1i - vgsol*k*rhogeq*w3r**2*w1r + vdsol*Kdrag*k*w3r*w1i +&
        vdsol*Kdrag*k*w2i*w1r - vdsol*Kdrag*k*w3r*w2i + vdsol*Kdrag*k*w3i*w1r - w3r**2*rhogsol*w2r**2 -&
        rhogsol*k**2*cs**2*w1i*w2i - vgsol*Kdrag*k*w2i*w1r - vgsol*Kdrag*k*w3i*w1r + 2*vgsol*k*rhogeq*w1i*w1r*w3i -&
        vgsol*k*rhogeq*w1r*w3i**2 - 2*vgsol*k*rhogeq*w2i*w1r*w3i - vgsol*k*rhogeq*w2i**2*w1r +&
        rhogsol*k**2*cs**2*w2r*w1r + w3r**2*w2r*w1r*rhogsol + w3r*w2r*rhogsol*w1i**2 - w3r*w2r*w1r**2*rhogsol +&
        w3r*w2r**2*w1r*rhogsol + vgsol*k*rhogeq*w2r*w1r**2 - vgsol*k*rhogeq*w1r*w2r**2 - vgsol*k*rhogeq*w2r*w1i**2 -&
        rhogsol*w2i**2*w3i**2 - rhogsol*w2i**2*w3r**2 + vgsol*k*rhogeq*w2r*w3i**2 - vgsol*Kdrag*k*w1i*w2r +&
        vgsol*Kdrag*k*w2r*w3i - 2*rhogsol*w2r*w1r*w3i*w1i + rhogsol*w2r**2*w3i*w1i + rhogsol*w2r*w1r*w3i**2 -&
        vdsol*Kdrag*k*w2r*w3i + vdsol*Kdrag*k*w1i*w2r - rhogsol*w2r**2*w3i**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 +&
        w3i**2 - 2*w3r*w1r)/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)

!-------------------------------
! D U S T  D E N S I T I E S
!-------------------------------
 if (Kdrag > 0.) then

    rhod3r = - rhodeq*( - w3r**3*rhogeq*w1i**2*w2r**2*rhogsol - w3r**2*Kdrag**2*k*vgsol*w2r*w1r +&
        w3r**2*Kdrag**2*k*vgsol*w2i*w1i - w3r**2*Kdrag**2*k*vdsol*w2i*w1i + w3r**2*Kdrag*k*rhogeq*vgsol*w2i*w1i**2 +&
        w3r**2*Kdrag**2*k*vdsol*w2r*w1r - w3r**2*Kdrag*rhogsol*k**2*cs**2*w2i*w1r +&
        w3r**2*Kdrag*k*rhogeq*vgsol*w2i*w1r**2 - w3r**2*Kdrag*rhogsol*k**2*cs**2*w1i*w2r +&
        w3r**2*rhogeq**2*w2i*w1i**2*w3i*k*vgsol - 4*w3r**2*rhogeq**2*w3i**2*k*w2i*w1i*vgsol -&
        2*w3r**2*rhogeq*w2i*w3i**2*k*Kdrag*vgsol + w3r**2*rhogeq*w2r*w3i*w1r*k*Kdrag*vdsol -&
        2*w3r**2*rhogeq**2*w1i**2*k*vgsol*w3i**2 - w3r**2*rhogeq*w2i*w3i*w1i*k*Kdrag*vgsol +&
        w3r**5*rhogeq*w2i*rhogsol*w1i + w3r**2*rhogeq**2*w2i*k*vgsol*w3i*w1r**2 - w3r**5*rhogeq*rhogsol*k**2*cs**2 -&
        w3r**5*rhogeq*w2r*w1r*rhogsol + w3r**5*rhogeq**2*k*w1r*vgsol + 3*w3r**2*rhogeq*cs**2*k**3*w3i*Kdrag*vdsol +&
        w3r**5*rhogeq**2*w2r*vgsol*k + 2*w3r**2*rhogeq*w3i**2*cs**2*k**2*rhogsol*w1r -&
        w3r**2*rhogeq*cs**4*k**4*w2r*rhogsol - 2*w3r**2*rhogeq**2*w3i**2*k*w1r**2*vgsol +&
        2*w3r**2*rhogeq*w2r*w3i**2*cs**2*k**2*rhogsol - w3r*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i**2 -&
        w3r**3*Kdrag*rhogsol*w2i*w1i**2 - w3r**3*Kdrag*rhogsol*w2r**2*w1i - w3r**3*Kdrag*rhogsol*w1i*w2i**2 -&
        w3r**3*Kdrag*rhogsol*w2i*w1r**2 - w3r**3*rhogeq*w2r*w1i*k*Kdrag*vdsol + w3r**3*rhogeq**2*w2r*w1r**2*k*vgsol +&
        w3r**3*rhogeq**2*w2r*w1i**2*k*vgsol - w3r**3*rhogeq**2*cs**2*k**3*w1r*vgsol +&
        w3r**3*rhogeq*w2r*w1i*k*Kdrag*vgsol - w3r**3*Kdrag**2*k*vdsol*w2r + 2*w3r**3*Kdrag*k*rhogeq*vgsol*w3i*w1r +&
        2*w3r**3*rhogeq**2*w2r*k*vgsol*w3i**2 + 2*w3r**3*Kdrag*k*rhogeq*vgsol*w2r*w3i +&
        w3r**3*Kdrag*rhogsol*k**2*cs**2*w1i - w3r**3*rhogeq**2*cs**2*k**3*w2r*vgsol -&
        2*w3r**3*Kdrag*rhogsol*k**2*cs**2*w3i + w3r**3*Kdrag*rhogsol*k**2*cs**2*w2i +&
        2*w3r**3*Kdrag*rhogsol*w2i*w3i*w1i - w3r**3*rhogeq*w2i**2*w1i**2*rhogsol - 2*w3r**3*Kdrag*rhogsol*w3i*w2r*w1r -&
        w3r**2*rhogeq*w2i*w3i*w1i*k*Kdrag*vdsol - w3r**3*rhogeq*w1r**2*w2i**2*rhogsol -&
        w3r**2*rhogeq*cs**4*k**4*w1r*rhogsol - 3*w3r**2*rhogeq**2*cs**2*k**3*w3i*w1i*vgsol -&
        3*w3r**2*rhogeq*cs**2*k**3*w3i*Kdrag*vgsol - w3r**2*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol -&
        w3r**2*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol + w3r**2*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol +&
        w3r**2*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol - w3r**3*rhogeq*w2r**2*w1r**2*rhogsol +&
        w3r**2*rhogeq**2*w2i**2*w3i*w1i*k*vgsol + w3r**2*rhogeq**2*cs**2*k**3*w1i**2*vgsol +&
        w3r**3*rhogeq*cs**4*k**4*rhogsol - 3*w3r**2*rhogeq**2*cs**2*k**3*w3i*w2i*vgsol +&
        2*w3r**2*rhogeq*w2r*w3i**2*w1r**2*rhogsol + w3r**2*rhogeq**2*cs**2*k**3*w1r**2*vgsol -&
        2*w3r**2*rhogeq**2*w3i**2*w2i**2*k*vgsol + 2*w3r**2*rhogeq*w2i**2*w3i**2*w1r*rhogsol +&
        w3r**2*rhogeq**2*w3i*w2r**2*w1i*k*vgsol + 2*w3r**2*rhogeq*w1i*w3i**2*k*Kdrag*vdsol +&
        2*w3r**2*rhogeq*w1i**2*w3i**2*w2r*rhogsol - 2*w3r**2*rhogeq*w1i*w3i**2*k*Kdrag*vgsol +&
        2*w3r**2*rhogeq*w2r**2*w3i**2*w1r*rhogsol - 4*w3r**2*rhogeq**2*w2r*w3i**2*k*w1r*vgsol -&
        3*w3r**2*rhogeq*w2r*w3i*w1r*k*Kdrag*vgsol - 2*w3r**2*rhogeq**2*w2r**2*w3i**2*k*vgsol +&
        2*w3r**2*rhogeq*w2i*w3i**2*k*Kdrag*vdsol + w3r**3*rhogeq**2*w2r**2*k*w1r*vgsol -&
        2*w3r**3*rhogeq*w3i**2*cs**2*k**2*rhogsol + w3r**4*rhogeq*w2i**2*w1r*rhogsol + w3r**4*Kdrag*rhogsol*w1i*w2r +&
        w3r**4*rhogeq*w2r*w1r**2*rhogsol + w3r**4*rhogeq*w1i**2*w2r*rhogsol - 2*w3r**4*rhogeq**2*k*w2i*w1i*vgsol -&
        2*w3r**4*rhogeq*w1i*k*Kdrag*vgsol - w3r**4*rhogeq**2*w2r**2*k*vgsol + w3r**4*rhogeq*w1i*k*Kdrag*vdsol -&
        w3r**4*rhogeq**2*w1r**2*k*vgsol + w3r**4*rhogeq*w1r*cs**2*k**2*rhogsol + w3r**4*rhogeq*w2r*cs**2*k**2*rhogsol -&
        2*w3r**4*rhogeq*w2i*k*Kdrag*vgsol - w3r**4*rhogeq**2*w1i**2*k*vgsol + w3r**4*rhogeq*w2r**2*w1r*rhogsol +&
        w3r**3*rhogeq*w2i*w1r*k*Kdrag*vgsol + w3r**4*Kdrag**2*k*vdsol + 2*w3r**3*rhogeq*w1i*w3i**2*rhogsol*w2i -&
        2*w3r**4*rhogeq**2*w2r*k*w1r*vgsol + 2*w3r**3*rhogeq**2*w3i**2*k*w1r*vgsol -&
        w3r**3*rhogeq*w2i*w1r*k*Kdrag*vdsol + w3r**3*rhogeq**2*w2i**2*w1r*k*vgsol -&
        2*w3r**3*rhogeq*w2r*w3i**2*w1r*rhogsol + w3r**4*rhogeq*w2i*k*Kdrag*vdsol - w3r**4*rhogeq**2*w2i**2*k*vgsol +&
        w3r**3*Kdrag**2*k*vgsol*w1r + w3r**4*Kdrag*rhogsol*w2i*w1r - w3r**4*Kdrag**2*k*vgsol -&
        w3r**3*Kdrag**2*k*vdsol*w1r + w3r**3*Kdrag**2*k*vgsol*w2r + w1i*w3i**3*k*Kdrag**2*vdsol +&
        2*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r + w3r*rhogeq**2*w3i**4*w2r*k*vgsol +&
        2*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r - 2*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r -&
        2*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r + w3r*rhogeq**2*w3i**4*k*w1r*vgsol -&
        w3r*rhogeq*w3i**4*w2r*w1r*rhogsol + 2*w3r*w3i*cs**4*k**4*rhogeq*rhogsol*w2i + w3r*rhogeq*w3i**4*w1i*rhogsol*w2i&
        - w2i*w3i**3*k*Kdrag**2*vgsol + 2*w3r*rhogeq*w3i**3*k*Kdrag*vgsol*w2r + w2i*w3i**3*k*Kdrag**2*vdsol -&
        w1i*w3i**3*k*Kdrag**2*vgsol - w3r*rhogeq*w3i**4*cs**2*k**2*rhogsol -&
        2*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w1i*w2i**2 - 2*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1i +&
        2*w3r*w3i**3*Kdrag*rhogsol*w2i*w1i - 2*w3r*w3i**3*Kdrag*rhogsol*w2r*w1r - 2*w3r*w3i**3*Kdrag*rhogsol*k**2*cs**2&
        + w3r*w3i**2*Kdrag**2*k*vgsol*w1r + 2*w3r*w3i*cs**4*k**4*rhogeq*rhogsol*w1i +&
        w3r*rhogeq*cs**2*k**2*w1i**2*w2r**2*rhogsol + w3r*rhogeq*cs**2*k**2*w2r**2*w1r**2*rhogsol -&
        w3r*rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol + w3r**2*rhogeq**2*cs**2*k**3*w2r**2*vgsol +&
        w3r**2*rhogeq*w2r**2*w1i*k*Kdrag*vgsol - 2*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w1i**2 -&
        2*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w1r**2 + w1r*w3i**2*cs**2*k**2*rhogeq*rhogsol*w2i**2 -&
        w1r*w3r**2*cs**2*k**2*rhogeq*rhogsol*w2i**2 - w2r*rhogeq*w1i**2*rhogsol*k**2*cs**2*w3r**2 -&
        w2r*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol + w2r*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol -&
        w2r*w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol - w2r*w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol +&
        w2r*w3r*rhogeq*cs**4*k**4*w1r*rhogsol + 2*w2r*w3r**2*rhogeq**2*cs**2*k**3*w1r*vgsol +&
        rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r*w1r - w3r*w3i**2*rhogeq*w2i*w1r*k*Kdrag*vdsol +&
        3*w3r*w3i**2*rhogeq**2*cs**2*k**3*w2r*vgsol + w3r*w3i**2*Kdrag*rhogsol*k**2*cs**2*w1i +&
        2*w3r*rhogeq*w3i**3*k*Kdrag*vgsol*w1r + rhogeq*w3i**2*cs**2*k**2*rhogsol*w1r*w2r**2 -&
        w3r**2*rhogeq*rhogsol*k**2*cs**2*w2r**2*w1r - w2r*w3r**2*cs**2*k**2*rhogeq*rhogsol*w1r**2 +&
        rhogeq*w2r*w3i**2*cs**2*k**2*rhogsol*w1r**2 - w3r*w3i**2*Kdrag*rhogsol*w2r**2*w1i -&
        w3r*w3i**2*Kdrag**2*k*vdsol*w2r + w3r*w3i**2*Kdrag**2*k*vgsol*w2r - rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r*w1r +&
        rhogeq**2*cs**2*k**3*w3i*w1i*vgsol*w2r**2 - w3r*w3i**2*Kdrag*rhogsol*w2i*w1r**2 +&
        w3r*w3i**2*Kdrag*rhogsol*k**2*cs**2*w2i - w3r*w3i**2*Kdrag*rhogsol*w1i*w2i**2 -&
        2*w3i**2*cs**2*k**3*rhogeq**2*vgsol*w2r*w1r - w3i**2*cs**2*k**3*rhogeq**2*vgsol*w2r**2 -&
        w3r*w3i**2*Kdrag*rhogsol*w2i*w1i**2 - w3i*rhogsol*k**4*cs**4*rhogeq*w2r*w1i +&
        rhogeq*cs**2*k**2*w1i**2*w2r*rhogsol*w3i**2 + Kdrag*k*rhogeq*vgsol*w2r**2*w1i*w3i**2 -&
        w3r*w3i**2*Kdrag**2*k*vdsol*w1r - rhogeq*w3i**5*w2r*w1i*rhogsol + w3r*w3i**2*rhogeq*w2r*w1i*k*Kdrag*vgsol +&
        w3r*w3i**2*rhogeq**2*w2r*w1i**2*k*vgsol + w3r*w3i**2*rhogeq*w2i*w1r*k*Kdrag*vgsol -&
        w3r*w3i**2*rhogeq*w2i**2*w1i**2*rhogsol + w3r*w3i**2*rhogeq**2*w2r**2*k*w1r*vgsol -&
        w3r*w3i**2*rhogeq*w1r**2*w2i**2*rhogsol - w3r*w3i**2*rhogeq*w2r**2*w1r**2*rhogsol -&
        3*w3r*w3i**2*rhogeq*cs**4*k**4*rhogsol + w3r*w3i**2*rhogeq**2*w2i**2*w1r*k*vgsol -&
        4*cs**2*k**2*rhogeq*w3r*w3i**2*w2r*w1r*rhogsol - rhogeq**2*w3i**4*w1i**2*k*vgsol -&
        cs**2*k**3*rhogeq*w3i**3*Kdrag*vdsol - w3r*w3i**2*rhogeq*w1i**2*w2r**2*rhogsol +&
        cs**2*k**3*rhogeq*w3i**3*Kdrag*vgsol + 3*w3r*w3i**2*rhogeq**2*cs**2*k**3*w1r*vgsol -&
        w3r*w3i**2*rhogeq*w2r*w1i*k*Kdrag*vdsol + w3r*w3i**2*rhogeq**2*w2r*w1r**2*k*vgsol +&
        4*cs**2*k**2*rhogeq*w3r*w3i**2*w2i*rhogsol*w1i - rhogeq**2*w3i**4*k*w1r**2*vgsol +&
        2*rhogeq*w3i**3*w3r**2*k*Kdrag*vgsol + 2*rhogeq**2*w3i**3*w3r**2*w1i*k*vgsol +&
        rhogeq*w3i**4*cs**2*k**2*rhogsol*w1r - 2*rhogeq*w3i**3*w3r**2*k*Kdrag*vdsol +&
        rhogeq*w3i**4*w2r*cs**2*k**2*rhogsol - 2*rhogeq**2*w3i**4*k*w2i*w1i*vgsol - rhogeq*w3i**3*w2i*w1i*k*Kdrag*vdsol&
        + rhogeq**2*w3i**3*w2i*w1i**2*k*vgsol + Kdrag*w2i**2*vgsol*k*rhogeq*w1i*w3r**2 -&
        w3i**2*rhogeq**2*cs**2*k**3*vgsol*w2i**2 + 2*w3r**2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i*w1i - rhogeq**2*w3i**4*w2r**2*k*vgsol +&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w1i + w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1r**2 +&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i**2 + rhogeq*w3i**4*w1i**2*w2r*rhogsol - rhogeq**2*w3i**4*w2i**2*k*vgsol&
        + rhogeq**2*w3i**5*w2i*k*vgsol + rhogeq*w3i**3*w2r*w1r*k*Kdrag*vdsol + rhogeq*w3i**4*w2i*k*Kdrag*vdsol +&
        rhogeq**2*w3i**3*w2i*k*vgsol*w1r**2 - 2*rhogeq**2*w3i**4*w2r*k*w1r*vgsol +&
        w3r*rhogeq*cs**2*k**2*w1r**2*rhogsol*w2i**2 + w3r*rhogeq*cs**2*k**2*w1i**2*rhogsol*w2i**2 -&
        2*rhogeq*w3i**3*w3r**2*w2i*w1r*rhogsol - 2*rhogeq*w3i**3*w2r*w3r**2*rhogsol*w1i +&
        w3r**2*rhogeq**2*cs**2*k**3*vgsol*w2i**2 - w3i*cs**4*k**4*rhogeq*rhogsol*w2i*w1r -&
        w3r*rhogeq*cs**4*k**4*rhogsol*w2i*w1i + 2*rhogeq**2*w3i**3*w3r**2*w2i*k*vgsol +&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w1r - w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w1r +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i*w1i + Kdrag*w3i**3*cs**2*k**2*rhogsol*w1r +&
        rhogeq*w3i**4*w2i**2*w1r*rhogsol - Kdrag*w3i**3*k*rhogeq*w1r**2*vgsol + Kdrag*w2r*w3i**3*cs**2*k**2*rhogsol -&
        rhogeq*w3i**5*w1r*w2i*rhogsol - Kdrag*w1r*w2i*w3i**4*rhogsol - rhogeq*w3i**5*k*Kdrag*vdsol +&
        rhogeq*w3i**5*k*Kdrag*vgsol - 2*w3i**2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i + rhogeq**2*w3i**5*w1i*k*vgsol +&
        rhogeq*w3i**4*w2r**2*w1r*rhogsol + rhogeq*w3i**4*w1i*k*Kdrag*vdsol - Kdrag*w1i**2*rhogeq*k*vgsol*w3i**3 +&
        rhogeq*w3i**4*w2r*w1r**2*rhogsol - Kdrag*w2r*w1i*w3i**4*rhogsol + w3i**4*k*Kdrag**2*vgsol -&
        w3i**4*k*Kdrag**2*vdsol + w3i**3*rhogeq**2*k*vgsol*w2r**2*w1i - w3i**2*Kdrag**2*k*vgsol*w2r*w1r +&
        w3i**3*Kdrag*rhogsol*w2r*w1r**2 + w3i**3*Kdrag*rhogsol*w2r*w1i**2 + w3i**2*Kdrag**2*k*vdsol*w2r*w1r +&
        w3i**3*cs**2*k**3*rhogeq**2*vgsol*w1i - w3i**2*Kdrag*cs**2*k**2*rhogsol*w2r*w1i +&
        w3i**3*Kdrag*rhogsol*w2r**2*w1r + w3i**2*rhogeq*cs**4*k**4*rhogsol*w2r -&
        2*w3i**3*cs**2*k**2*rhogeq*rhogsol*w2r*w1i - w3i**3*Kdrag*rhogeq*k*vgsol*w2r**2 -&
        3*w3i**3*Kdrag*rhogeq*k*vgsol*w2r*w1r - w3i**3*Kdrag*rhogeq*k*vgsol*w2i**2 -&
        w3i**2*rhogeq**2*cs**2*k**3*vgsol*w1i**2 + w3i**2*rhogeq*k*Kdrag*vgsol*w1i**2*w2i -&
        w3i**2*rhogeq**2*cs**2*k**3*vgsol*w1r**2 + w3i**2*rhogeq*k*Kdrag*vgsol*w2i*w1r**2 -&
        w3i**3*Kdrag*rhogeq*k*vgsol*w1i*w2i - w3i**2*Kdrag**2*k*vdsol*w1i*w2i + w3i**2*Kdrag**2*k*vgsol*w1i*w2i -&
        w3i**2*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i - w3i**2*cs**2*k**3*rhogeq*Kdrag*vgsol*w1i +&
        w3i**3*Kdrag*rhogsol*w2i**2*w1r + w3i**2*rhogeq*cs**2*k**3*vdsol*Kdrag*w1i +&
        w3i**2*rhogeq*cs**2*k**3*vdsol*Kdrag*w2i - w3i**2*Kdrag*cs**2*k**2*rhogsol*w2i*w1r -&
        2*w3i**3*rhogeq*cs**2*k**2*rhogsol*w2i*w1r + w3i**3*vgsol*k*rhogeq**2*w2i**2*w1i +&
        w3i**3*vgsol*rhogeq**2*k**3*cs**2*w2i + w3i**2*rhogeq*cs**4*k**4*rhogsol*w1r -&
        w3i*Kdrag*w3r**2*w1r**2*rhogeq*k*vgsol - w3i*Kdrag*w3r**2*w1i**2*k*rhogeq*vgsol -&
        w3i*w3r**2*w1i*k*Kdrag**2*vgsol + w3i*w3r**2*w1i*k*Kdrag**2*vdsol + w3i*Kdrag*w3r**2*w1r*cs**2*k**2*rhogsol +&
        w3i*Kdrag*w2r*w3r**2*cs**2*k**2*rhogsol + Kdrag*w2i**2*k*rhogeq*vgsol*w3i**2*w1i -&
        w3i*rhogeq*w3r**4*k*Kdrag*vdsol - w3i*rhogeq*w3r**4*rhogsol*w2i*w1r - w3i*rhogeq*w3r**4*rhogsol*w1i*w2r +&
        w3i*rhogeq**2*w3r**4*k*vgsol*w2i - w3i*w3r**2*Kdrag*rhogeq*k*vgsol*w2i**2 + w3i*w3r**2*Kdrag*rhogsol*w2r*w1r**2&
        - w3i*w3r**2*w2i*k*Kdrag**2*vgsol + w3i*w3r**2*w2i*k*Kdrag**2*vdsol + w3i*rhogeq*w3r**4*k*Kdrag*vgsol +&
        w3i*rhogeq**2*w3r**4*k*vgsol*w1i + w3i*w3r**2*Kdrag*rhogsol*w2r**2*w1r + w3i*w3r**2*Kdrag*rhogsol*w2r*w1i**2 +&
        w3i*w3r**2*Kdrag*rhogsol*w2i**2*w1r + 2*w3i*w3r**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i -&
        w3i*w3r**2*Kdrag*rhogeq*k*vgsol*w2r**2 + 2*w3i*w3r**2*rhogeq*cs**2*k**2*rhogsol*w2i*w1r)/(w3r**2 +&
        w3i**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 -&
        2*w2i*w3i + w3r**2)/rhogeq/Kdrag

    rhod3i = - rhodeq*(w3r*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w1i + w2r*w3i*cs**2*k**3*rhogeq**2*vgsol*w1i**2 +&
        w2r*w3i*cs**2*k**3*rhogeq**2*vgsol*w1r**2 - w3r**2*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1i +&
        rhogeq*w2r**2*w3i**2*cs**2*k**2*rhogsol*w1i - w3i**2*Kdrag*rhogsol*w2i**2*w1r**2 -&
        w3i**2*Kdrag*rhogsol*w2i**2*w1i**2 - w3r**2*Kdrag*rhogsol*w2i**2*w1i**2 +&
        rhogeq*w1i*rhogsol*k**2*cs**2*w3i**2*w2i**2 - w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1i**2 -&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1r**2 - w2r*w3i*cs**4*k**4*rhogeq*rhogsol*w1r -&
        w2r*w3r*rhogeq*cs**4*k**4*rhogsol*w1i + w2r*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1i +&
        rhogeq*w1i**2*rhogsol*k**2*cs**2*w3i**2*w2i - rhogeq*w1i**2*rhogsol*k**2*cs**2*w3r**2*w2i +&
        w2r*w3r*rhogeq*cs**2*k**3*Kdrag*vdsol*w1r - w2r*w3r*rhogeq*cs**2*k**3*Kdrag*vgsol*w1r +&
        w3i*cs**4*k**4*rhogeq*rhogsol*w1i*w2i - w3i*cs**2*k**2*rhogeq*rhogsol*w2i**2*w1r**2 +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r*w2i - w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r*w2i -&
        rhogeq*w1i*rhogsol*k**2*cs**2*w3r**2*w2i**2 - w3i*cs**2*k**2*rhogeq*rhogsol*w2i**2*w1i**2 +&
        w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2i + w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol*w2i +&
        w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol*w2i - w3r*rhogeq*cs**4*k**4*w1r*rhogsol*w2i -&
        w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2i + w3r**2*rhogeq*w2i**2*w1r*k*Kdrag*vgsol +&
        w3r**2*rhogeq*w1r*k*Kdrag*vgsol*w2r**2 - w3r**2*Kdrag*rhogsol*w1i**2*w2r**2 -&
        w3r**2*Kdrag*rhogsol*w2i**2*w1r**2 - w2r*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1i -&
        w3i**2*Kdrag*rhogsol*w1r**2*w2r**2 - w3i**2*Kdrag*rhogsol*w1i**2*w2r**2 - w3r**2*Kdrag*rhogsol*w1r**2*w2r**2 +&
        w3i**2*rhogeq*w2i**2*w1r*k*Kdrag*vgsol + w3i**2*rhogeq*w1r*k*Kdrag*vgsol*w2r**2 -&
        w3r**2*cs**2*k**2*rhogeq*rhogsol*w1r**2*w2i + w3i**2*cs**2*k**2*rhogeq*rhogsol*w1r**2*w2i +&
        rhogeq*w3r**4*rhogsol*w2r**2*w1i + rhogeq*w3r**4*rhogsol*w1i*w2i**2 + rhogeq*w3r**4*rhogsol*w2i*w1i**2 +&
        rhogeq*w3r**4*rhogsol*w2i*w1r**2 - w3r*rhogeq*w2r*w3i**2*w1r*k*Kdrag*vgsol - rhogeq**2*w3i**5*w2r*k*vgsol -&
        rhogeq**2*w3i**5*k*w1r*vgsol + rhogeq*w3i**5*w2r*w1r*rhogsol - w3i**2*cs**4*k**4*rhogeq*rhogsol*w2i +&
        2*w3r**2*rhogeq*w3i**3*cs**2*k**2*rhogsol - 3*w3r*rhogeq*w2i*w3i**2*w1i*k*Kdrag*vgsol -&
        w3r*rhogeq*w2r*w3i**2*w1r*k*Kdrag*vdsol - 2*w3r**2*Kdrag*k*rhogeq*vgsol*w3i**2*w1r -&
        2*w3r**2*Kdrag*k*rhogeq*vgsol*w2r*w3i**2 - 2*w3r**2*rhogeq**2*w2r*k*vgsol*w3i**3 +&
        3*w3r*rhogeq**2*cs**2*k**3*w3i**2*w1i*vgsol + 3*w3r*rhogeq*cs**2*k**3*w3i**2*Kdrag*vgsol -&
        w3r*rhogeq**2*w2i**2*w3i**2*w1i*k*vgsol + 3*w3r*rhogeq**2*cs**2*k**3*w3i**2*w2i*vgsol -&
        w3r*rhogeq**2*w3i**2*w2r**2*w1i*k*vgsol + 2*w3r*rhogeq*w1i*w3i**3*k*Kdrag*vgsol -&
        w3r*rhogeq**2*w2i*w1i**2*w3i**2*k*vgsol + 2*w3r*rhogeq*w2i*w3i**3*k*Kdrag*vgsol -&
        w3r*rhogeq**2*w2i*k*vgsol*w3i**2*w1r**2 - 3*w3r*rhogeq*cs**2*k**3*w3i**2*Kdrag*vdsol -&
        w3i**3*rhogeq*w2i**2*w1i**2*rhogsol + w3i**2*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r -&
        w3i**3*rhogeq*w1r**2*w2i**2*rhogsol + w3i**3*rhogeq**2*w2r**2*k*w1r*vgsol + w3i**3*rhogeq**2*w2i**2*w1r*k*vgsol&
        + w3i**3*rhogeq*cs**4*k**4*rhogsol - rhogeq*w3i**5*w1i*rhogsol*w2i - w3i**4*Kdrag*rhogsol*w2i*w1i +&
        rhogeq*w3i**5*cs**2*k**2*rhogsol + w3i**4*Kdrag*rhogsol*w2r*w1r + w3i**2*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r +&
        w3i**4*Kdrag*rhogsol*k**2*cs**2 - w3i**2*cs**4*k**4*rhogeq*rhogsol*w1i - w3i**3*Kdrag**2*k*vgsol*w1r -&
        w3i**3*rhogeq**2*cs**2*k**3*w2r*vgsol + w3i**3*Kdrag*rhogsol*w2r**2*w1i - w3i**3*Kdrag*rhogsol*k**2*cs**2*w1i -&
        2*rhogeq*w3i**4*k*Kdrag*vgsol*w1r - 2*w3r**2*rhogeq*w1i*w3i**3*rhogsol*w2i -&
        2*w3r**2*rhogeq**2*w3i**3*k*w1r*vgsol + 2*w3r**2*rhogeq*w2r*w3i**3*w1r*rhogsol +&
        w3r*rhogeq*w2i*w3i**2*w1i*k*Kdrag*vdsol - w3i**3*rhogeq*w2i*w1r*k*Kdrag*vdsol -&
        w3i**3*Kdrag*rhogsol*k**2*cs**2*w2i - w3i**2*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r +&
        w3i**3*Kdrag*rhogsol*w2i*w1i**2 - w3i**2*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r + w3i**3*Kdrag**2*k*vdsol*w1r +&
        w3i**3*rhogeq**2*w2r*w1i**2*k*vgsol - w3i**3*rhogeq*w2r**2*w1r**2*rhogsol - 2*rhogeq*w3i**4*k*Kdrag*vgsol*w2r -&
        Kdrag*w3r**4*rhogsol*k**2*cs**2 - w3i**3*Kdrag**2*k*vgsol*w2r - rhogeq*w3r**5*rhogsol*w2i*w1r +&
        w3i**3*Kdrag*rhogsol*w1i*w2i**2 + w3i**3*rhogeq*w2r*w1i*k*Kdrag*vgsol + Kdrag*w3r**4*w2i*rhogsol*w1i +&
        w3i**3*rhogeq*w2i*w1r*k*Kdrag*vgsol + rhogeq**2*w3i**4*w3r*w1i*k*vgsol - rhogeq*w3i**4*w3r*k*Kdrag*vdsol +&
        w3i**3*rhogeq**2*w2r*w1r**2*k*vgsol + w3i**3*Kdrag**2*k*vdsol*w2r + rhogeq*w3i**4*w3r*k*Kdrag*vgsol -&
        w3i**3*rhogeq**2*cs**2*k**3*w1r*vgsol - rhogeq*w3i**4*w2r*w3r*rhogsol*w1i - w3i**3*rhogeq*w1i**2*w2r**2*rhogsol&
        - rhogeq*w3i**4*w3r*w2i*w1r*rhogsol - w3i**3*rhogeq*w2r*w1i*k*Kdrag*vdsol - 2*w3r*w3i**3*k*Kdrag**2*vdsol +&
        2*w3i**2*rhogeq**2*w3r**3*k*vgsol*w1i + 2*w3i**2*rhogeq*w3r**3*k*Kdrag*vgsol +&
        w3i**2*w3r*Kdrag*rhogsol*w2r*w1i**2 + w3i**2*w3r*Kdrag*rhogsol*w2r**2*w1r -&
        w3i**2*w3r*Kdrag*rhogeq*k*vgsol*w2r**2 + w3i**2*w3r*Kdrag*rhogsol*w2i**2*w1r -&
        2*w3i**2*w3r*cs**2*k**2*rhogeq*rhogsol*w2r*w1i - 2*w3i**2*w3r*rhogeq*cs**2*k**2*rhogsol*w2i*w1r -&
        w3i*w3r**2*rhogeq*w1i**2*w2r**2*rhogsol + 2*w3i*w3r*rhogeq*cs**4*k**4*w2r*rhogsol -&
        w3i*w3r**4*rhogeq*w2i*rhogsol*w1i + w3i*w3r**4*rhogeq*rhogsol*k**2*cs**2 + w3i*w3r**4*rhogeq*w2r*w1r*rhogsol -&
        w3i*w3r**4*rhogeq**2*w2r*vgsol*k + rhogeq**2*w3i**4*w3r*w2i*k*vgsol - w3i*w3r**4*rhogeq**2*k*w1r*vgsol +&
        w3i*w3r**2*Kdrag*rhogsol*w2r**2*w1i + w3i*w3r**2*Kdrag*rhogsol*w2i*w1r**2 + w3i*w3r**2*Kdrag*rhogsol*w1i*w2i**2&
        + w3i*w3r**2*Kdrag**2*k*vdsol*w2r + w3i*w3r**2*rhogeq**2*w2r*w1r**2*k*vgsol +&
        w3i*w3r**2*rhogeq**2*w2r*w1i**2*k*vgsol + 3*w3i*w3r**2*rhogeq**2*cs**2*k**3*w1r*vgsol +&
        2*w3i*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol - w3i*w3r**2*Kdrag*rhogsol*k**2*cs**2*w1i +&
        3*w3i*w3r**2*rhogeq**2*cs**2*k**3*w2r*vgsol + 2*w3i*w3r*rhogeq*cs**4*k**4*w1r*rhogsol +&
        2*w3i*w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol - w3i*w3r**2*Kdrag*rhogsol*k**2*cs**2*w2i -&
        w3i*w3r**2*rhogeq*w1r**2*w2i**2*rhogsol - w3i*w3r**2*rhogeq*w2i**2*w1i**2*rhogsol +&
        w3i**2*w3r*w1i*k*Kdrag**2*vdsol - w3i**2*Kdrag*w3r*w1r**2*rhogeq*k*vgsol -&
        w3i**2*Kdrag*w3r*w1i**2*k*rhogeq*vgsol + w3i**2*Kdrag*w3r*w1r*cs**2*k**2*rhogsol +&
        w3i**2*Kdrag*w2r*w3r*cs**2*k**2*rhogsol - w3i**2*w3r*w1i*k*Kdrag**2*vgsol -&
        2*w3i**2*rhogeq*w3r**3*k*Kdrag*vdsol - 2*w3i**2*rhogeq*w3r**3*rhogsol*w2i*w1r -&
        2*w3i**2*rhogeq*w3r**3*rhogsol*w1i*w2r + w3i**2*w3r*Kdrag*rhogsol*w2r*w1r**2 -&
        w3i**2*w3r*Kdrag*rhogeq*k*vgsol*w2i**2 - 2*w3i*w3r**3*Kdrag**2*k*vdsol - w3i**2*w3r*w2i*k*Kdrag**2*vgsol +&
        w3i**2*w3r*w2i*k*Kdrag**2*vdsol - w3i*w3r**2*Kdrag**2*k*vgsol*w1r + 2*w3i*w3r**3*rhogeq*w2i*k*Kdrag*vgsol -&
        2*w3i*w3r**3*Kdrag*rhogsol*w2i*w1r + w3i*w3r**2*Kdrag**2*k*vdsol*w1r -&
        2*w3i*w3r*rhogeq**2*cs**2*k**3*vgsol*w2i**2 - 2*w3i*w3r*rhogeq**2*cs**2*k**3*w2r**2*vgsol -&
        w3i*w3r**2*rhogeq*w2r**2*w1r**2*rhogsol - 2*w3i*w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol -&
        2*w3i*w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol - 3*w3i*w3r**2*rhogeq*cs**4*k**4*rhogsol +&
        w3i*w3r**2*rhogeq**2*w2r**2*k*w1r*vgsol - w3i*w3r**2*rhogeq*w2r*w1i*k*Kdrag*vdsol -&
        2*w3i*w3r**3*Kdrag*rhogsol*w1i*w2r + 2*w3i*w3r**3*rhogeq*w1i*k*Kdrag*vgsol +&
        w3i*w3r**2*rhogeq*w2r*w1i*k*Kdrag*vgsol + w3i*w3r**2*Kdrag*rhogsol*w2i*w1i**2 +&
        rhogeq*rhogsol*w2r**2*w1i*w3i**4 - rhogeq*w3r**5*k*Kdrag*vdsol + rhogeq*w3r**3*k*Kdrag*vdsol*w2i*w1i -&
        rhogeq*w3r**3*k*Kdrag*vdsol*w2r*w1r + 2*rhogeq*w3r**2*rhogsol*w2r**2*w1i*w3i**2 +&
        w3i*w3r**2*rhogeq*w2i*w1r*k*Kdrag*vgsol + 2*rhogeq*w3r**2*rhogsol*w1i*w2i**2*w3i**2 +&
        2*rhogeq*w3r**2*rhogsol*w2i*w1r**2*w3i**2 + 2*rhogeq*w3r**2*rhogsol*w2i*w1i**2*w3i**2 -&
        2*Kdrag*w2r*w3r*w3i**3*rhogsol*w1i + rhogeq*w3r**4*k*Kdrag*vdsol*w2r + 2*rhogeq*w3r**2*k*Kdrag*vdsol*w2r*w3i**2&
        - Kdrag*w3r**3*w1r**2*rhogeq*k*vgsol - Kdrag*w3r**3*w1i**2*k*rhogeq*vgsol + Kdrag*w3r**3*w1r*cs**2*k**2*rhogsol&
        + Kdrag*w2r*w3r**3*cs**2*k**2*rhogsol + cs**2*k**3*rhogeq*w3r**3*Kdrag*vdsol +&
        w3i*w3r**2*rhogeq**2*w2i**2*w1r*k*vgsol - rhogeq**2*w3r**3*k*vgsol*w2i*w1r**2 -&
        rhogeq*rhogsol*k**2*cs**2*w3i**4*w2i + rhogeq**2*w3r**5*k*vgsol*w1i + rhogeq*w3r**5*k*Kdrag*vgsol -&
        w3i*w3r**2*rhogeq*w2i*w1r*k*Kdrag*vdsol - w3i*w3r**2*Kdrag**2*k*vgsol*w2r + rhogeq*rhogsol*w3i**4*w2i*w1r**2 +&
        2*w3i*w3r**3*Kdrag**2*k*vgsol + w3i*rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol + rhogeq*rhogsol*w2i*w1i**2*w3i**4 +&
        rhogeq*k*Kdrag*vdsol*w2r*w3i**4 + 2*w3i*w1r*w3r*cs**2*k**2*rhogeq*rhogsol*w2i**2 + rhogeq**2*w3r**5*k*vgsol*w2i&
        + 2*rhogeq*w3r**2*k*Kdrag*vdsol*w3i**2*w1r - 2*w3i*w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol -&
        rhogeq*w3r**5*rhogsol*w1i*w2r + rhogeq*k*Kdrag*vdsol*w3i**4*w1r - 2*rhogeq*w3r**2*rhogsol*k**2*cs**2*w3i**2*w1i&
        - 2*w3i*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol - Kdrag*w2r*w3r**4*w1r*rhogsol + w3r**3*Kdrag*rhogsol*w2r**2*w1r&
        + Kdrag**2*k*vgsol*w2r*w1i*w3i**2 + w3r**2*Kdrag*k*rhogeq*vgsol*w2r*w1i**2 +&
        w3r**2*Kdrag*k*rhogeq*vgsol*w2r*w1r**2 + w3r**2*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r -&
        w3r**2*Kdrag*cs**2*k**2*rhogsol*w2r*w1r - w3r**2*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r +&
        rhogeq**2*cs**2*k**3*w1r*vgsol*w2i**2*w3i + 2*w3r**3*cs**2*k**2*rhogeq*rhogsol*w2r*w1i +&
        2*w3i**2*rhogeq**2*w3r**3*k*vgsol*w2i + w3r**3*Kdrag*rhogsol*w2r*w1r**2 - w3r**3*Kdrag*rhogeq*k*vgsol*w2i**2 -&
        w3r**3*rhogeq**2*cs**2*k**3*vgsol*w1i + 2*w3i*w2r*w3r*cs**2*k**2*rhogeq*rhogsol*w1r**2 -&
        w3r**3*rhogeq**2*k*vgsol*w2r**2*w1i - w3r**3*Kdrag*rhogeq*k*vgsol*w2r**2 - rhogeq*w3r**4*rhogsol*k**2*cs**2*w2i&
        + rhogeq*rhogsol*w2i**2*w3i**4*w1i - rhogeq*w3r**4*rhogsol*k**2*cs**2*w1i -&
        4*w3i*w3r*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i + w3r**3*w1i*k*Kdrag**2*vdsol - w3r**3*w2i*k*Kdrag**2*vgsol +&
        rhogeq*w3r**4*k*Kdrag*vdsol*w1r + w3r**3*w2i*k*Kdrag**2*vdsol - rhogeq*rhogsol*k**2*cs**2*w3i**4*w1i -&
        rhogeq**2*w3r**3*k*vgsol*w2i*w1i**2 - 2*rhogeq*w3r**2*rhogsol*k**2*cs**2*w3i**2*w2i -&
        w3r**3*w1i*k*Kdrag**2*vgsol - 4*w3i*w2r*w3r*rhogeq**2*cs**2*k**3*w1r*vgsol +&
        w3r**2*cs**4*k**4*rhogeq*rhogsol*w2i + 2*w3i*w3r*rhogeq*rhogsol*k**2*cs**2*w2r**2*w1r +&
        2*w3r*w3i**3*k*Kdrag**2*vgsol - w3r**2*Kdrag**2*k*vdsol*w2i*w1r - Kdrag**2*k*vdsol*w2i*w1r*w3i**2 -&
        w3r**3*vgsol*k*rhogeq**2*w2i**2*w1i - w3r**3*Kdrag*rhogeq*k*vgsol*w2r*w1r +&
        2*w3r**3*rhogeq*cs**2*k**2*rhogsol*w2i*w1r - 3*w3r**3*Kdrag*k*rhogeq*vgsol*w1i*w2i +&
        w3r**2*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r - w3r**2*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r +&
        Kdrag**2*k*vgsol*w2i*w1r*w3i**2 + rhogeq**2*cs**2*k**3*vgsol*w1i*w2r**2*w3r - w3r**2*Kdrag**2*k*vdsol*w2r*w1i -&
        Kdrag**2*k*vdsol*w2r*w1i*w3i**2 + w3r**2*cs**4*k**4*rhogeq*rhogsol*w1i +&
        w3r**2*Kdrag*cs**2*k**2*rhogsol*w1i*w2i + w3r**3*Kdrag*rhogsol*w2r*w1i**2 -&
        cs**2*k**3*rhogeq*w3r**3*Kdrag*vgsol + 2*w3i*w2r*rhogeq*w1i**2*rhogsol*k**2*cs**2*w3r +&
        4*cs**2*k**2*rhogeq*w3r**2*w3i*w2i*rhogsol*w1i - 4*cs**2*k**2*rhogeq*w3r**2*w3i*w2r*w1r*rhogsol -&
        2*Kdrag*w1r*w2i*w3i**3*rhogsol*w3r + Kdrag*k*rhogeq*vgsol*w2r*w1r**2*w3i**2 + w3r**2*Kdrag**2*k*vgsol*w2r*w1i -&
        Kdrag*cs**2*k**2*rhogsol*w2r*w1r*w3i**2 + w3i**3*Kdrag*rhogsol*w2i*w1r**2 +&
        Kdrag*k*rhogeq*vgsol*w2r*w1i**2*w3i**2 + w3r**2*Kdrag**2*k*vgsol*w2i*w1r +&
        Kdrag*cs**2*k**2*rhogsol*w1i*w2i*w3i**2 + w3r**3*Kdrag*rhogsol*w2i**2*w1r -&
        w3r**3*vgsol*rhogeq**2*k**3*cs**2*w2i)/(w3r**2 + w3i**2)/(w1i**2 - 2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 -&
        2*w3r*w1r)/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i + w3r**2)/rhogeq/Kdrag

    rhod2r = - rhodeq*(k*Kdrag**2*vdsol*w2i**3*w1i - k*Kdrag**2*vgsol*w2i**3*w1i -&
        2*w3r*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i**2 + w3r*Kdrag*rhogsol*w2i**3*w1i**2 -&
        Kdrag*k*rhogeq*vgsol*w2i**3*w1i**2 + w3r*Kdrag*rhogsol*w2i**3*w1r**2 - Kdrag*k*rhogeq*vgsol*w2i**3*w1r**2 -&
        2*w2i*w3r**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i + w2i*w3r**2*rhogeq**2*k*vgsol*w2r**2*w1i -&
        w2i*w3r**2*Kdrag*rhogeq*k*vgsol*w2r**2 - w3r*rhogeq**2*cs**2*k**3*w2r**3*vgsol -&
        w3r*Kdrag*rhogsol*k**2*cs**2*w1i*w2r**2 + w3r*Kdrag**2*k*vdsol*w2r**2*w1r - w3r*Kdrag**2*k*vgsol*w2r**2*w1r +&
        w3r*Kdrag*k*rhogeq*vgsol*w2r**3*w1i - w3r*rhogeq*cs**2*k**2*w1i**2*w2r**2*rhogsol -&
        w3r*rhogeq*cs**2*k**2*w2r**2*w1r**2*rhogsol + 2*w3r*rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol +&
        w3r**2*rhogeq**2*cs**2*k**3*w2r**2*vgsol + w3r**2*rhogeq*w2r**2*w1i*k*Kdrag*vgsol -&
        w3r**2*Kdrag*rhogsol*w2r**3*w1i + w3r**2*rhogeq**2*w2r**3*k*w1r*vgsol + w2i*w3r*rhogeq*k*Kdrag*vdsol*w2r**2*w1r&
        + w2i*w3r*Kdrag*rhogsol*w2r**2*w1r**2 + w2i*w3r*Kdrag*rhogsol*w2r**2*w1i**2 - w2i*w3r*rhogeq*rhogsol*w2r**4*w1i&
        + 2*w2i*w3r*rhogeq*cs**2*k**2*rhogsol*w2r**2*w1i + 2*w2i*w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r +&
        w2i*w3r*Kdrag*cs**2*k**2*rhogsol*w2r**2 - 2*w2i*w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r -&
        4*w3r*rhogeq*cs**2*k**2*rhogsol*w2r*w2i**2*w1r - w3r*rhogeq*k*Kdrag*vdsol*w2r*w2i**2*w1i +&
        rhogeq**2*k*vgsol*w2i**5*w1i - w3r*rhogeq*rhogsol*w2i**5*w1i + 4*w3i*cs**2*k**2*rhogeq*rhogsol*w2r*w2i**2*w1i +&
        2*w3r**2*rhogeq*rhogsol*w2r**2*w2i**2*w1r - w3r**2*rhogeq*rhogsol*w2r*w2i**2*w1r**2 +&
        w3r**2*rhogeq**2*k*vgsol*w2r*w2i**2*w1r + 2*vgsol*k*rhogeq**2*w2i**3*w2r**2*w1i +&
        2*w3r*rhogeq*cs**2*k**2*rhogsol*w2i**2*w2r**2 + 3*w3r*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w2r -&
        2*w3r*Kdrag*rhogsol*w2r*w1r*w2i**3 - 2*w3r*rhogeq*rhogsol*w2i**3*w2r**2*w1i -&
        w3r**2*rhogeq*rhogsol*w2r*w2i**2*w1i**2 + 2*w3r*rhogeq*rhogsol*w1i**2*w2r**2*w2i**2 +&
        2*w3r*rhogeq*rhogsol*w2i**2*w1r**2*w2r**2 + w3r*w2r*Kdrag*rhogeq*k*vgsol*w2i**2*w1i -&
        2*w3r**2*rhogeq**2*k*vgsol*w2i**2*w2r**2 + w3r*rhogeq**2*k*vgsol*w2r*w2i**2*w1i**2 +&
        w3r*rhogeq**2*k*vgsol*w2r*w2i**2*w1r**2 - 4*w3r*rhogeq**2*k*vgsol*w2r**2*w2i**2*w1r -&
        w3r*rhogeq*cs**4*k**4*w2r**2*rhogsol - w3r**2*rhogeq*w2r**3*w1r**2*rhogsol +&
        w1r*Kdrag*rhogsol*k**2*cs**2*w2i**3 + w2i*w3r**2*Kdrag*rhogsol*w2r**2*w1r - 2*w2i*w3r*Kdrag*rhogsol*w2r**3*w1r&
        - 3*w2i*w3r*Kdrag*k*rhogeq*vgsol*w2r**2*w1r + w2i*rhogeq**2*k*vgsol*w2r**4*w1i +&
        rhogeq**2*w1i**2*k*vgsol*w3r*w2r**3 + 2*w2i*w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1r +&
        w1r*w3i**2*cs**2*k**2*rhogeq*rhogsol*w2i**2 + w1r*w3r**2*cs**2*k**2*rhogeq*rhogsol*w2i**2 +&
        2*Kdrag*w3r*w2r**3*k*rhogeq*vgsol*w2i + w3r*rhogeq**2*w1r**2*k*vgsol*w2r**3 - w3r*rhogeq*w2r**5*w1r*rhogsol -&
        rhogeq*w1i*w3r*w2r**3*k*Kdrag*vdsol + rhogeq*w1i**2*w3r*w2r**4*rhogsol + w3r*rhogeq*w2r**4*w1r**2*rhogsol +&
        w3r*rhogeq**2*w2r**5*k*vgsol + w3r*w2r**3*k*Kdrag**2*vgsol - w3r*w2r**3*k*Kdrag**2*vdsol +&
        Kdrag*w3r*w2r**4*rhogsol*w1i + w3r*rhogeq*w2r**4*cs**2*k**2*rhogsol +&
        w2r*rhogeq*w1i**2*rhogsol*k**2*cs**2*w3r**2 - w2r*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol +&
        w2r*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol - w2r*w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol -&
        w2r*w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol + w2r*w3r*rhogeq*cs**4*k**4*w1r*rhogsol -&
        w2r*w3r**2*rhogeq**2*cs**2*k**3*w1r*vgsol - w2r*w3r**2*Kdrag*rhogsol*w1i*w2i**2 -&
        rhogeq*w2r*w3i*w1r*k*Kdrag*vdsol*w2i**2 + rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r*w1r -&
        rhogeq*w2i*w3i*w1i*k*Kdrag*vdsol*w2r**2 - rhogeq*w3i**2*cs**2*k**2*rhogsol*w1r*w2r**2 -&
        2*w3r*rhogeq**2*w2r**4*k*w1r*vgsol - 2*w3r*rhogeq*w2r**3*w1r*rhogsol*w2i**2 +&
        2*w3r*rhogeq**2*w2i**2*k*vgsol*w2r**3 - w3r**2*rhogeq**2*w2r**4*vgsol*k -&
        w3r**2*rhogeq*rhogsol*k**2*cs**2*w2r**2*w1r + w3r**2*rhogeq*w2r**4*w1r*rhogsol +&
        w2r*w3r**2*cs**2*k**2*rhogeq*rhogsol*w1r**2 - w2r*w3r*rhogeq*rhogsol*w2i**4*w1r +&
        Kdrag*rhogsol*k**2*cs**2*w1i*w2r**3 - rhogeq**2*w2r**4*w3i**2*k*vgsol +&
        rhogeq*w2r*w3i**2*cs**2*k**2*rhogsol*w1r**2 + rhogeq*cs**4*k**4*w2r**3*rhogsol +&
        w2r*w3r*rhogeq**2*k*vgsol*w2i**4 + rhogeq*w2r**4*w3i**2*w1r*rhogsol + rhogeq**2*cs**2*k**3*w2r**2*vgsol*w1r**2&
        - 3*rhogeq*cs**4*k**4*w2r*rhogsol*w2i**2 - 2*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w2r*w1r +&
        2*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w2r*w1r - rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r**2 -&
        2*rhogeq*w2r*w3i**2*cs**2*k**2*rhogsol*w2i*w1i - rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r*w1r -&
        rhogeq**2*w1i**2*k*vgsol*w2r**4 + rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r**2 +&
        2*rhogeq**2*cs**2*k**3*w3i*w1i*vgsol*w2r**2 - 3*rhogeq**2*cs**2*k**3*w2r**2*vgsol*w2i*w1i +&
        2*rhogeq*cs**4*k**4*w2r*rhogsol*w2i*w1i - rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2r**2 -&
        3*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w2r**2 - 2*rhogeq**2*w1i**2*k*vgsol*w2r**2*w2i**2 +&
        rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2r**2 + rhogeq**2*cs**2*k**3*w1i**2*vgsol*w2r**2 +&
        3*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w2r**2 - 4*rhogeq**2*w3i*w2r**2*w1i*k*vgsol*w2i**2 -&
        2*rhogeq**2*w3i*w2r**4*w1i*k*vgsol + 2*rhogeq*w2r**2*w3i**2*w1r*rhogsol*w2i**2 +&
        2*w2r*Kdrag*w3r*w2i**3*k*rhogeq*vgsol + rhogeq*cs**2*k**2*w2r**4*w1r*rhogsol +&
        rhogeq**2*w2r**3*w3i**2*k*w1r*vgsol - rhogeq**2*cs**2*k**3*w2r**3*w1r*vgsol -&
        rhogeq*cs**4*k**4*w1r*rhogsol*w2r**2 + w2r*w3r*w2i**2*k*Kdrag**2*vgsol - w2r*w3r*w2i**2*k*Kdrag**2*vdsol -&
        rhogeq*w2r*w3i**2*w1r**2*rhogsol*w2i**2 + 2*rhogeq*cs**2*k**2*w2r**2*w1r*rhogsol*w2i**2 +&
        rhogeq**2*w2i*k*vgsol*w3i*w1r**2*w2r**2 - 3*rhogeq**2*cs**2*k**3*w3i*w2i*vgsol*w2r**2 -&
        2*Kdrag*k*rhogeq*vgsol*w2r**4*w1i - rhogeq*w2r**3*w3i**2*w1r**2*rhogsol + w3i*vdsol*Kdrag*k*rhogeq*w2r**4 -&
        vdsol*Kdrag*k*rhogeq*w2i*w2r**4 - 2*w3i*rhogeq*rhogsol*w2r**2*w2i**3*w1r + w3i*rhogeq*rhogsol*w2r*w2i**4*w1i -&
        w3i*rhogeq*rhogsol*w2i**5*w1r + w3i*rhogeq*rhogsol*w2r**5*w1i - w3i*rhogeq*rhogsol*w2i*w2r**4*w1r +&
        2*w3i*rhogeq*rhogsol*w2i**2*w2r**3*w1i - w3i**2*cs**2*k**3*rhogeq**2*vgsol*w2r*w1r +&
        w3i**2*cs**2*k**3*rhogeq**2*vgsol*w2r**2 - 2*rhogsol*k**2*cs**2*Kdrag*w2i*w2r**3 +&
        w3i*rhogsol*k**2*cs**2*Kdrag*w2r**3 + w3i*rhogsol*k**2*cs**2*Kdrag*w2r*w2i**2 -&
        2*rhogsol*k**2*cs**2*Kdrag*w2r*w2i**3 - w3i*vgsol*k*Kdrag**2*w2i*w2r**2 + vdsol*k*Kdrag**2*w2r**4 -&
        cs**2*k**2*rhogeq*rhogsol*w2r*w2i**4 + 2*w3i*rhogeq**2*k*vgsol*w2r**2*w2i**3 -&
        2*cs**2*k**2*rhogeq*rhogsol*w2r**3*w2i**2 + w3i*rhogeq**2*k*vgsol*w2i*w2r**4 +&
        2*w3i*Kdrag*rhogsol*w2r*w2i**3*w1i - cs**2*k**2*rhogeq*rhogsol*w2r**5 + w3i*rhogeq**2*k*vgsol*w2i**5 +&
        2*rhogeq**2*k*vgsol*w2i**2*w2r**3*w1r + w3i*Kdrag*rhogsol*w2r**4*w1r + 2*w3i*Kdrag*rhogsol*w2i*w2r**3*w1i +&
        rhogeq*w2r**2*w3i*w1r**2*k*Kdrag*vgsol + w3i*vdsol*k*Kdrag**2*w2i*w2r**2 - vgsol*k*Kdrag**2*w2r**4 +&
        2*vdsol*Kdrag*k*rhogeq*w2i**2*w2r**2*w1i - vdsol*Kdrag*k*rhogeq*w2i**5 + vdsol*Kdrag*k*rhogeq*w2r**4*w1i -&
        2*vdsol*Kdrag*k*rhogeq*w2r**2*w2i**3 - w3i*rhogsol*k**4*cs**4*rhogeq*w2r*w1i +&
        2*vgsol*k*rhogeq*Kdrag*w2i*w2r**3*w1r + 2*vgsol*k*rhogeq*Kdrag*w2r**2*w2i**3 +&
        2*vgsol*k*rhogeq*Kdrag*w2r*w2i**3*w1r + vgsol*k*rhogeq*Kdrag*w2i*w2r**4 - 2*w3i*vgsol*k*rhogeq*Kdrag*w2r**4 +&
        vgsol*k*rhogeq*Kdrag*w2i**5 - Kdrag*w2r**3*rhogsol*w1i**2*w3i - Kdrag*w2r**2*k*rhogeq*vgsol*w2i*w3i**2 -&
        Kdrag*w2i**2*rhogsol*w3i*w2r*w1r**2 + Kdrag*w2r**2*rhogsol*k**2*cs**2*w2i*w1r -&
        Kdrag*w2i**2*rhogsol*w1i**2*w2r*w3i + w2i**2*k*Kdrag**2*vgsol*w2r*w1r - w2r**3*k*Kdrag**2*vdsol*w1r -&
        w2i**2*k*Kdrag**2*vdsol*w2r*w1r + rhogeq**2*w1i**2*w2r**2*k*vgsol*w2i*w3i +&
        rhogeq*cs**2*k**2*w1i**2*w2r*rhogsol*w3i**2 + Kdrag*k*rhogeq*vgsol*w2r**2*w1i*w3i**2 +&
        Kdrag*k*rhogeq*vgsol*w2r**2*w1i**2*w3i + w2r**2*k*Kdrag**2*vgsol*w3i*w1i - w2r**2*k*Kdrag**2*vdsol*w3i*w1i +&
        rhogeq**2*k*vgsol*w2r**5*w1r - 2*w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w2r*w1i**2 -&
        2*w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w2r*w1r**2 - 2*rhogeq**2*k*vgsol*w2i**2*w2r**2*w1r**2 +&
        rhogeq**2*k*vgsol*w2r*w2i**4*w1r - Kdrag*w2r**2*k*rhogeq*vgsol*w2i*w3i*w1i -&
        Kdrag*w2r**2*k*rhogeq*vgsol*w2i*w1i**2 + Kdrag*w2r**3*k*rhogeq*vgsol*w3i*w1r -&
        Kdrag*w2i*vgsol*k*rhogeq*w1r**2*w2r**2 + Kdrag*w2i**2*k*rhogeq*vgsol*w2r*w3i*w1r + w2r**3*k*Kdrag**2*vgsol*w1r&
        - Kdrag*w2r**3*rhogsol*w1i*w3i**2 - Kdrag*w2r**2*rhogsol*k**2*cs**2*w3i*w1r +&
        rhogeq**2*w1i*w2r**2*k*vgsol*w2i*w3i**2 + Kdrag*w2r**2*rhogsol*w2i*w1r*w3i**2 -&
        rhogeq*w1i**2*w2r**3*rhogsol*w3i**2 - Kdrag*w2r**3*rhogsol*w3i*w1r**2 - rhogeq*w1i**2*w2i**2*rhogsol*w2r*w3i**2&
        - Kdrag*w2i**2*rhogsol*w1i*w2r*w3i**2 - w3r**2*rhogeq*w1i**2*w2r**3*rhogsol +&
        2*rhogsol*k**4*cs**4*rhogeq*w2r*w2i*w3i + 3*cs**2*k**3*rhogeq**2*vgsol*w2r*w1r*w2i**2 +&
        2*rhogeq*w2i**2*k*Kdrag*vdsol*w2r**2*w3i - rhogeq*w2r**3*w3i*w1r*k*Kdrag*vdsol +&
        Kdrag*rhogsol*k**2*cs**2*w1i*w2r*w2i**2 - 2*Kdrag*k*rhogeq*vgsol*w2r**2*w1i*w2i**2 -&
        2*rhogeq**2*w2r**2*w3i**2*k*vgsol*w2i**2 - Kdrag**2*k*vgsol*w2i*w1i*w2r**2 + Kdrag**2*k*vdsol*w2i*w1i*w2r**2 +&
        w2r*rhogeq**2*vgsol*k*w2i**2*w1r*w3i**2 - 2*Kdrag*w2r**2*k*rhogeq*vgsol*w2i**2*w3i - k*Kdrag**2*vdsol*w2i**4 +&
        k*Kdrag**2*vgsol*w2i**4 + Kdrag*w2i**2*vgsol*k*rhogeq*w1i*w3r**2 - w3i**2*rhogeq**2*cs**2*k**3*vgsol*w2i**2 +&
        w2i**2*k*Kdrag**2*vdsol*w3r*w1r - w2i**2*k*Kdrag**2*vdsol*w3i*w1i + w2i**2*k*Kdrag**2*vgsol*w3i*w1i -&
        Kdrag*rhogsol*k**2*cs**2*w3i*w1r*w2i**2 - w2i**2*k*Kdrag**2*vgsol*w3r*w1r +&
        w3r**2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i - w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i*w1i -&
        rhogeq**2*vgsol*k*w2i**4*w1r**2 - rhogeq**2*vgsol*k*w2i**4*w1i**2 - w3r**2*vgsol*k*rhogeq*Kdrag*w2i**3 -&
        w3r*Kdrag*rhogsol*w2i**4*w1i + cs**2*k**3*rhogeq*vgsol*Kdrag*w2i**3 + w3r**2*rhogeq**2*vgsol*k*w2i**3*w1i +&
        cs**2*k**3*rhogeq**2*vgsol*w2i**3*w1i + w3r*vdsol*Kdrag*k*rhogeq*w2i**3*w1r - w3i*vgsol*k*Kdrag**2*w2i**3 +&
        w3i*vdsol*k*Kdrag**2*w2i**3 - w3i**2*vgsol*k*rhogeq*Kdrag*w2i**3 - w3i*vgsol*k*rhogeq*Kdrag*w2i**3*w1i -&
        3*w3r*vgsol*k*rhogeq*Kdrag*w2i**3*w1r - w3i*Kdrag*rhogsol*w2i**4*w1r + w3r*rhogsol*k**2*cs**2*Kdrag*w2i**3 +&
        w3r**2*rhogeq*rhogsol*w2i**4*w1r + w3i**2*rhogeq*rhogsol*w2i**4*w1r + w3i**2*Kdrag*rhogsol*w2i**3*w1r +&
        w3r**2*Kdrag*rhogsol*w2i**3*w1r - 2*w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w1i +&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**3 + w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1r**2 +&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i**2 - cs**2*k**3*rhogeq*vdsol*Kdrag*w2i**3 -&
        2*w3r*rhogeq**2*k*vgsol*w2i**4*w1r - w3r**2*rhogeq**2*vgsol*k*w2i**4 + w3r*cs**2*k**2*rhogeq*rhogsol*w2i**4 -&
        w3i**2*rhogeq**2*vgsol*k*w2i**4 - w3i*vdsol*Kdrag*k*rhogeq*w2i**3*w1i + w3i*vdsol*Kdrag*k*rhogeq*w2i**4 +&
        rhogeq*rhogsol*k**2*cs**2*w2i**4*w1r + rhogeq*k*Kdrag*vdsol*w2i**4*w1i - 2*w3i*rhogeq**2*k*vgsol*w2i**4*w1i +&
        w3i*rhogeq**2*k*vgsol*w2i**3*w1i**2 + w3i*rhogeq**2*k*vgsol*w2i**3*w1r**2 -&
        2*w3r*cs**2*k**2*rhogeq*rhogsol*w2i**3*w1i + w3i**2*rhogeq**2*vgsol*k*w2i**3*w1i +&
        w3r*rhogeq*cs**2*k**2*w1r**2*rhogsol*w2i**2 + w3r*rhogeq*cs**2*k**2*w1i**2*rhogsol*w2i**2 +&
        w3r*rhogeq*cs**4*k**4*rhogsol*w2i**2 - w3r*Kdrag*rhogsol*k**2*cs**2*w1i*w2i**2 -&
        w3r**2*rhogeq**2*cs**2*k**3*vgsol*w2i**2 - w3i*cs**4*k**4*rhogeq*rhogsol*w2i*w1r -&
        w3r*rhogeq*cs**4*k**4*rhogsol*w2i*w1i + rhogeq*cs**4*k**4*w1r*rhogsol*w2i**2 -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i**2 + w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i**2 -&
        rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2i**2 + rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2i**2 +&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w1r - rhogeq**2*cs**2*k**3*w1r**2*vgsol*w2i**2 -&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w1r - rhogeq**2*cs**2*k**3*w1i**2*vgsol*w2i**2 +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i*w1i + Kdrag*w2i**2*k*rhogeq*vgsol*w3i*w1i**2 +&
        w3i**2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i - 2*w3i*cs**2*k**2*rhogeq*rhogsol*w2i**3*w1r +&
        w3r*rhogeq*rhogsol*w2i**4*w1r**2 + w3r*rhogeq*rhogsol*w2i**4*w1i**2 - rhogeq**2*k*vgsol*w2r**4*w1r**2 +&
        Kdrag*w2i**2*k*rhogeq*vgsol*w3i**2*w1i + Kdrag*w2i**2*k*rhogeq*vgsol*w3i*w1r**2)/(w2i**2 +&
        w2r**2)/rhogeq/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i + w3r**2)/(w2r**2 + w1r**2 + w2i**2 -&
        2*w2i*w1i - 2*w2r*w1r + w1i**2)/Kdrag

    rhod2i =rhodeq*( - w3r**2*Kdrag*rhogsol*w1i*w2i**3 - w2r*w3i*vdsol*k*Kdrag**2*w2i**2 +&
        w2r*w3i**2*vgsol*k*rhogeq*Kdrag*w2i**2 + 3*w2r*w3i*vgsol*k*rhogeq*Kdrag*w2i**2*w1i +&
        w2r*w3r*vgsol*k*rhogeq*Kdrag*w2i**2*w1r + w2r*w3r**2*vgsol*k*rhogeq*Kdrag*w2i**2 -&
        3*w2r*cs**2*k**3*rhogeq*vgsol*Kdrag*w2i**2 - w2r*w3r*rhogsol*k**2*cs**2*Kdrag*w2i**2 +&
        2*w2r*w3i*Kdrag*rhogsol*w2i**3*w1r - w2r*w3i*cs**2*k**3*rhogeq**2*vgsol*w1i**2 +&
        3*w2r*cs**2*k**3*rhogeq*vdsol*Kdrag*w2i**2 - w2r*w3r**2*Kdrag*rhogsol*w2i**2*w1r -&
        3*w2r*w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**2 - 2*Kdrag*k*rhogeq*vgsol*w2r**3*w1i*w2i +&
        4*w2r*w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i - w2r*w3i**2*Kdrag*rhogsol*w2i**2*w1r -&
        w2r*w3i*cs**2*k**3*rhogeq**2*vgsol*w1r**2 - w2r*w3i*vdsol*Kdrag*k*rhogeq*w2i**2*w1i +&
        w2r*w3i*vgsol*k*Kdrag**2*w2i**2 + Kdrag*vgsol*k*rhogeq*w1r**2*w2r**3 - Kdrag*w2i*k*rhogeq*vgsol*w2r**2*w3i*w1r&
        - Kdrag*w2r**3*rhogsol*w1r*w3i**2 - 3*cs**2*k**3*rhogeq**2*vgsol*w2r**2*w1r*w2i +&
        Kdrag*rhogsol*k**2*cs**2*w1i*w2r**2*w2i - w3i**2*Kdrag*rhogsol*w1i*w2i**3 +&
        4*w2r*w3r*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i - w2r*w3r*Kdrag*rhogsol*w2i**2*w1i**2 -&
        Kdrag**2*k*vdsol*w1i*w2r**3 - w2r**2*rhogeq**2*vgsol*k*w2i*w1r*w3i**2 -&
        2*w2r*w1r*w3r**2*cs**2*k**2*rhogeq*rhogsol*w2i - w2r*k*Kdrag**2*vdsol*w2i**2*w1i +&
        w2r*w3r*rhogeq*rhogsol*w2i**4*w1i - 2*w2r*w1r*w3i**2*cs**2*k**2*rhogeq*rhogsol*w2i +&
        w2r*Kdrag*k*rhogeq*vgsol*w2i**2*w1i**2 - w2r*rhogeq**2*k*vgsol*w2i**4*w1i - w2r*w3r*Kdrag*rhogsol*w2i**2*w1r**2&
        + w2r*k*Kdrag**2*vgsol*w2i**2*w1i - 2*Kdrag*w2r**3*k*rhogeq*vgsol*w2i*w3i + Kdrag**2*k*vgsol*w1i*w2r**3 +&
        rhogeq*w1i**2*w2i*rhogsol*w2r**2*w3i**2 - Kdrag*w2i*rhogsol*w1i*w2r**2*w3i**2 -&
        rhogsol*k**4*cs**4*rhogeq*w2r**2*w3i + rhogeq**2*w1i*w2r**3*k*vgsol*w3i**2 +&
        w3r**2*rhogeq**2*k*vgsol*w2r**3*w1i + w3r**2*Kdrag*rhogeq*k*vgsol*w2r**3 + w2r*w3i*rhogeq*rhogsol*w2i**4*w1r -&
        w2r*w3i*rhogeq**2*k*vgsol*w2i**4 + w3r**2*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1i +&
        w2r*vdsol*Kdrag*k*rhogeq*w2i**4 - w2r*vgsol*k*rhogeq*Kdrag*w2i**4 + 2*w2r*w3i**2*rhogeq**2*cs**2*k**3*vgsol*w2i&
        + w2r*Kdrag*k*rhogeq*vgsol*w2i**2*w1r**2 - w2r*w1r*Kdrag*rhogsol*k**2*cs**2*w2i**2 -&
        2*vgsol*k*rhogeq**2*w2i**2*w2r**3*w1i - 3*w3r*cs**2*k**3*rhogeq**2*vgsol*w2i*w2r**2 +&
        w3r*rhogeq*k*Kdrag*vdsol*w2r**3*w1r - w3r*Kdrag*rhogsol*w2r**3*w1r**2 - w3r*Kdrag*rhogsol*w2r**3*w1i**2 +&
        w3r*rhogeq*rhogsol*w2r**5*w1i - 4*w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w2i*w1i +&
        w3r**2*rhogeq*rhogsol*w2r**2*w2i*w1r**2 - w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r**2 +&
        w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r**2 + w3r*rhogeq*k*Kdrag*vdsol*w2r**2*w2i*w1i -&
        w3r*Kdrag*cs**2*k**2*rhogsol*w2r**3 - 2*w3r*rhogeq*cs**2*k**2*rhogsol*w2r**3*w1i +&
        4*w3r*rhogeq*cs**2*k**2*rhogsol*w2r**2*w2i*w1r - w3r**2*rhogeq**2*k*vgsol*w2r**2*w2i*w1r -&
        w3r*rhogeq**2*k*vgsol*w2r**2*w2i*w1r**2 - w3r**2*Kdrag*rhogsol*w2r**3*w1r + w3r*Kdrag*rhogsol*w2r**4*w1r +&
        w3r*Kdrag*k*rhogeq*vgsol*w2r**3*w1r - 2*w2r*k*Kdrag**2*vgsol*w2i**3 + rhogeq*w2r**2*w3i*w1r*k*Kdrag*vdsol*w2i -&
        rhogeq*w3i*w1i*k*Kdrag*vdsol*w2r**3 + 2*w2r**2*w3r*rhogeq**2*k*vgsol*w2i**3 -&
        2*w3i*cs**2*k**2*rhogeq*rhogsol*w2r**3*w1r - w3r*rhogeq**2*k*vgsol*w2r**2*w2i*w1i**2 -&
        rhogeq**2*k*vgsol*w2r**5*w1i - w2r**2*w3r**2*Kdrag*rhogsol*w1i*w2i + 2*w3r*rhogeq*rhogsol*w2i**2*w2r**3*w1i +&
        w3r**2*rhogeq*rhogsol*w2r**2*w2i*w1i**2 - w3r*w2r**2*Kdrag*rhogeq*k*vgsol*w2i*w1i +&
        rhogeq**2*cs**2*k**3*w2r**3*vgsol*w1i + 3*rhogeq*cs**4*k**4*w2r**2*rhogsol*w2i +&
        rhogeq*cs**2*k**3*Kdrag*vdsol*w2r**2*w1r - rhogeq*cs**2*k**3*Kdrag*vgsol*w2r**2*w1r -&
        2*w2r**2*w3r*rhogeq*rhogsol*w2i**3*w1r - w3r*rhogeq*w2r**4*w1r*rhogsol*w2i + w3r*rhogeq**2*w2i*k*vgsol*w2r**4 +&
        2*w2r*k*Kdrag**2*vdsol*w2i**3 - rhogeq*cs**2*k**3*Kdrag*vdsol*w2r**3 - w2r**2*w3r*w2i*k*Kdrag**2*vdsol +&
        2*w2r**2*Kdrag*w3r*w2i**2*k*rhogeq*vgsol + rhogeq*w2r**2*w3i**2*w1r**2*rhogsol*w2i +&
        rhogeq**2*k*vgsol*w3i*w1r**2*w2r**3 - rhogeq*cs**4*k**4*w2r**2*rhogsol*w1i +&
        rhogeq*cs**2*k**3*Kdrag*vgsol*w2r**3 + rhogeq*w2r**2*w3i**2*cs**2*k**2*rhogsol*w1i +&
        w3i**2*Kdrag*rhogsol*w2i**2*w1r**2 + rhogeq**2*cs**2*k**3*w3i*vgsol*w2r**3 +&
        2*w3i*rhogeq*rhogsol*w2r**2*w2i**3*w1i + w3i*rhogeq*rhogsol*w2r**5*w1r + w3i*rhogeq*rhogsol*w2i*w2r**4*w1i +&
        w2r**2*w3r*w2i*k*Kdrag**2*vgsol + rhogsol*k**2*cs**2*Kdrag*w2r**4 - 2*w3i*rhogeq**2*k*vgsol*w2r**3*w2i**2 -&
        cs**2*k**2*rhogeq*rhogsol*w2r**4*w2i - 2*cs**2*k**2*rhogeq*rhogsol*w2r**2*w2i**3 - w3i*rhogeq**2*k*vgsol*w2r**5&
        + rhogeq**2*k*vgsol*w2i*w2r**4*w1r + vdsol*Kdrag*k*rhogeq*w2r**5 + 2*w3i*rhogeq*rhogsol*w2r**3*w2i**2*w1r +&
        w3i**2*Kdrag*rhogsol*w2i**2*w1i**2 + 2*vdsol*Kdrag*k*rhogeq*w2r**3*w2i**2 +&
        2*vgsol*k*rhogeq*Kdrag*w2r**2*w2i**2*w1r + w3i*rhogsol*k**2*cs**2*Kdrag*w2r**2*w2i +&
        w3i*vgsol*k*Kdrag**2*w2r**3 - w2i*k*Kdrag**2*vdsol*w2r**2*w1r + rhogeq**2*w1i**2*w2r**3*k*vgsol*w3i -&
        w3i*vdsol*k*Kdrag**2*w2r**3 - w3i*Kdrag*rhogsol*w2r**4*w1i - 2*vgsol*k*rhogeq*Kdrag*w2r**3*w2i**2 -&
        Kdrag*w2i*rhogsol*w1i**2*w2r**2*w3i + w2i*k*Kdrag**2*vgsol*w2r**2*w1r + 3*Kdrag*w2r**3*k*rhogeq*vgsol*w3i*w1i +&
        Kdrag*w2r**3*k*rhogeq*vgsol*w1i**2 + w3r**2*Kdrag*rhogsol*w2i**2*w1i**2 -&
        rhogeq*w1i*rhogsol*k**2*cs**2*w3i**2*w2i**2 + w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1i**2 +&
        w3i*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1r**2 + 2*rhogeq**2*k*vgsol*w2r**2*w2i**3*w1r -&
        vgsol*k*rhogeq*Kdrag*w2r**5 + Kdrag*w2r**3*k*rhogeq*vgsol*w3i**2 - Kdrag*w2i*rhogsol*w3i*w2r**2*w1r**2 -&
        Kdrag*w2r**3*rhogsol*k**2*cs**2*w1r + 2*Kdrag*w2i**4*vgsol*k*rhogeq*w1r - 2*Kdrag*w2i**3*w2r*vgsol*k*rhogeq*w1i&
        + Kdrag*w2i**4*rhogsol*w3i*w1i + rhogeq*w1i*w2r**4*rhogsol*k**2*cs**2 +&
        2*rhogeq*w1i*w2r**2*rhogsol*k**2*cs**2*w2i**2 - Kdrag*w3r*w2i**4*rhogsol*w1r +&
        w2r*w3i*cs**4*k**4*rhogeq*rhogsol*w1r + rhogeq*w1i*w3r*w2i**3*k*Kdrag*vdsol +&
        2*Kdrag*w2r**3*rhogsol*w3i*w1r*w2i - 2*w2r*w3r*rhogeq*cs**4*k**4*rhogsol*w2i +&
        2*w2r*w3r**2*rhogeq**2*cs**2*k**3*vgsol*w2i + rhogeq*w1i*w2i**4*rhogsol*k**2*cs**2 +&
        w2r*w3i*rhogeq**2*k*vgsol*w2i**2*w1r**2 + 2*w2r*w3r*cs**2*k**2*rhogeq*rhogsol*w2i**2*w1i -&
        2*w2r*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i + 2*w2r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2i -&
        2*w2r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2i + w2r*w3r*rhogeq*cs**4*k**4*rhogsol*w1i +&
        2*w2r*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i - w2r*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1i -&
        2*w2r*rhogeq*cs**4*k**4*w1r*rhogsol*w2i + w2r*w3i*rhogeq**2*k*vgsol*w2i**2*w1i**2 +&
        w2r*w3i**2*rhogeq**2*vgsol*k*w2i**2*w1i - 2*w2r*w3r*rhogeq*cs**2*k**2*w1r**2*rhogsol*w2i -&
        2*w2r*w3r*rhogeq*cs**2*k**2*w1i**2*rhogsol*w2i - w2r*w3i**2*rhogeq**2*cs**2*k**3*vgsol*w1i +&
        2*w2r*w3i*cs**2*k**2*rhogeq*rhogsol*w2i**2*w1r - rhogeq**2*w1i**2*k*vgsol*w3r*w2i**3 +&
        rhogeq*w1i**2*rhogsol*k**2*cs**2*w3i**2*w2i + rhogeq*w1i**2*rhogsol*k**2*cs**2*w3r**2*w2i -&
        w2r*w3r*rhogeq*cs**2*k**3*Kdrag*vdsol*w1r + 2*w2r*rhogeq**2*cs**2*k**3*w1r**2*vgsol*w2i +&
        w2r*w3r*rhogeq*cs**2*k**3*Kdrag*vgsol*w1r + 2*w2r*rhogeq**2*cs**2*k**3*w1i**2*vgsol*w2i -&
        w3i*cs**4*k**4*rhogeq*rhogsol*w1i*w2i - w3i*cs**2*k**2*rhogeq*rhogsol*w2i**2*w1r**2 +&
        w3r*rhogeq**2*cs**2*k**3*w2i**3*vgsol - w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r*w2i +&
        rhogeq**2*cs**2*k**3*w2i**3*w1r*vgsol + w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r*w2i -&
        rhogeq*w1i*rhogsol*k**2*cs**2*w3r**2*w2i**2 - w3r*w2i**3*k*Kdrag**2*vdsol +&
        w3r*Kdrag*rhogsol*k**2*cs**2*w1r*w2r**2 - w3r*Kdrag*k*rhogeq*vgsol*w2i**2*w1r**2 -&
        w3r*Kdrag*k*rhogeq*vgsol*w1r**2*w2r**2 - w3i*cs**2*k**2*rhogeq*rhogsol*w2i**2*w1i**2 +&
        w3i*cs**4*k**4*rhogeq*rhogsol*w2i**2 + rhogeq*cs**4*k**4*w2i**2*w1i*rhogsol -&
        w3r*Kdrag*k*rhogeq*vgsol*w1i**2*w2r**2 + w3r*Kdrag**2*k*vdsol*w2i**2*w1i + w3r*Kdrag**2*k*vdsol*w1i*w2r**2 +&
        w3r*Kdrag*rhogsol*k**2*cs**2*w2i**2*w1r - w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2i -&
        w3r*rhogeq*cs**2*k**3*w2i**2*Kdrag*vdsol + w3r**2*rhogeq*w2i**3*w1i**2*rhogsol - 2*w2i*k*Kdrag**2*vgsol*w2r**3&
        - w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol*w2i - w3r*Kdrag**2*k*vgsol*w2i**2*w1i - w3r*Kdrag**2*k*vgsol*w1i*w2r**2&
        - w3r*Kdrag*k*rhogeq*vgsol*w1i*w2i**3 - w3r*Kdrag*k*rhogeq*vgsol*w2i**2*w1i**2 +&
        w3r**2*rhogeq*w1r**2*w2i**3*rhogsol - w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol*w2i +&
        w3r*rhogeq*cs**4*k**4*w1r*rhogsol*w2i - w3r**2*rhogeq**2*w2i**3*w1r*k*vgsol -&
        w3r**2*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i + w3r*rhogeq*cs**2*k**3*w2i**2*Kdrag*vgsol +&
        w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2i - w3r**2*rhogeq*w2i**2*w1r*k*Kdrag*vgsol -&
        w3r**2*rhogeq*w1r*k*Kdrag*vgsol*w2r**2 + w3r**2*Kdrag*rhogsol*w1i**2*w2r**2 +&
        w3r**2*Kdrag*rhogsol*w2i**2*w1r**2 - Kdrag*rhogsol*k**2*cs**2*w3i*w1i*w2r**2 -&
        Kdrag*k*rhogeq*vgsol*w3i*w1r*w2i**3 - Kdrag*rhogsol*w1i**2*w2i**3*w3i +&
        w2r*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1i - Kdrag*rhogsol*k**2*cs**2*w2i**2*w3i*w1i -&
        Kdrag**2*k*vgsol*w2i**2*w3i*w1r + w3i**2*Kdrag*rhogsol*w1r**2*w2r**2 + w3i**2*Kdrag*rhogsol*w1i**2*w2r**2 +&
        w3r**2*Kdrag*rhogsol*w1r**2*w2r**2 + w2r*w3r**2*rhogeq**2*vgsol*k*w2i**2*w1i -&
        3*w2r*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w1i + w2r*w3r*vdsol*Kdrag*k*rhogeq*w2i**2*w1r -&
        w2r*w3r**2*rhogeq**2*cs**2*k**3*vgsol*w1i + 2*w2r*w3r*Kdrag*rhogsol*w2i**3*w1i -&
        w3i**2*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i - rhogeq*cs**2*k**3*w2i**2*w1r*Kdrag*vdsol -&
        Kdrag*rhogsol*w2i**3*w1r**2*w3i - w3i**2*rhogeq*w2i**2*w1r*k*Kdrag*vgsol -&
        w3i**2*rhogeq*w1r*k*Kdrag*vgsol*w2r**2 - rhogeq*w1r*k*Kdrag*vdsol*w2r**4 -&
        2*rhogeq*w2i**2*w1r*k*Kdrag*vdsol*w2r**2 - 2*rhogeq*w1i*w3i**2*rhogsol*w2i**2*w2r**2 -&
        rhogeq*w1i*w3i**2*rhogsol*w2r**4 - 2*w3i*rhogeq*rhogsol*w2i**2*w1r**2*w2r**2 -&
        w3r*rhogeq**2*k*vgsol*w2i**3*w1r**2 - w3i*rhogeq*rhogsol*w1r**2*w2r**4 +&
        2*w3i*rhogeq*rhogsol*k**2*cs**2*w2i**2*w2r**2 + w3i*rhogeq*rhogsol*k**2*cs**2*w2r**4 +&
        w3i**2*rhogeq*w1r**2*w2i**3*rhogsol - w3r**2*rhogeq*rhogsol*w2i**4*w1i -&
        2*w3r**2*rhogeq*rhogsol*w2i**2*w1i*w2r**2 - w3i**2*rhogeq*rhogsol*w2i**4*w1i +&
        w3i**2*rhogeq*w2i**3*w1i**2*rhogsol - rhogeq*rhogsol*k**2*cs**2*w2i**5 +&
        w3r**2*cs**2*k**2*rhogeq*rhogsol*w1r**2*w2i - w3i*rhogeq*rhogsol*w2i**4*w1r**2 +&
        w3i**2*cs**2*k**2*rhogeq*rhogsol*w1r**2*w2i + w3i*cs**2*k**2*rhogeq*rhogsol*w2i**4 +&
        2*w2i*k*Kdrag**2*vdsol*w2r**3 - w3r*rhogeq*rhogsol*w2i**5*w1r - rhogeq*w1i**2*w2i**4*rhogsol*w3i -&
        cs**4*k**4*rhogeq*rhogsol*w2i**3 + w3r*w2i**3*k*Kdrag**2*vgsol - rhogeq*k*Kdrag*vdsol*w2i**4*w1r +&
        w3i*vdsol*Kdrag*k*rhogeq*w2i**3*w1r + w3i*rhogeq*rhogsol*w2i**5*w1i + w3r*rhogeq**2*k*vgsol*w2i**5 -&
        w3r*vdsol*Kdrag*k*rhogeq*w2i**4 - 2*w3r*vdsol*Kdrag*k*rhogeq*w2i**2*w2r**2 -&
        2*Kdrag*w2i**3*k*rhogeq*vgsol*w2r*w3i - w2i**3*k*Kdrag**2*vdsol*w1r + w2i**3*k*Kdrag**2*vgsol*w1r +&
        2*Kdrag*w3r*w2r**3*rhogsol*w1i*w2i + Kdrag*w2i**3*rhogsol*k**2*cs**2*w3i -&
        2*rhogeq*w1i**2*w2i**2*rhogsol*w3i*w2r**2 - rhogeq*w1i**2*w2r**4*rhogsol*w3i - Kdrag*w2i**4*rhogsol*k**2*cs**2&
        + Kdrag*w2i**3*rhogsol*k**2*cs**2*w1i - w3r**2*rhogeq*rhogsol*w1i*w2r**4 - w3r*rhogeq*k*Kdrag*vdsol*w2r**4 +&
        2*Kdrag*w3r*w2i**4*k*rhogeq*vgsol - w3i**2*rhogeq**2*w2i**3*w1r*k*vgsol + Kdrag**2*k*vdsol*w2i**2*w3i*w1r +&
        Kdrag**2*k*vdsol*w3i*w1r*w2r**2 - Kdrag**2*k*vgsol*w3i*w1r*w2r**2 + rhogeq*cs**2*k**3*w2i**2*w1r*Kdrag*vgsol +&
        rhogeq**2*w1r*w2i**5*k*vgsol)/(w2i**2 + w2r**2)/rhogeq/(w2r**2 - 2*w3r*w2r + w2i**2 + w3i**2 - 2*w2i*w3i +&
        w3r**2)/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r + w1i**2)/Kdrag

    rhod1r =(2*w2i*w3r*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r*w3i**2 - w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w1r*w3i**2 +&
        rhogeq*w2r*w3i*w1r*k*Kdrag*vdsol*w2i**2*w3r**2 + w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w1r*w3i**2 -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i*w1i*w3r**2 - Kdrag*w2i**2*k*rhogeq*vgsol*w2r*w3i*w1r*w3r**2 +&
        Kdrag*w2r**2*k*rhogeq*vgsol*w2i*w3i*w1i*w3r**2 + 4*w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w2r*w1i**2*w3r**2 -&
        rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r*w1r*w3r**2 + rhogeq*w2i*w3i*w1i*k*Kdrag*vdsol*w2r**2*w3r**2 +&
        w2i*w3r*rhogeq*k*Kdrag*vdsol*w2r**2*w1r*w3i**2 - 4*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w2r*w1r*w3r**2 +&
        4*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w2r*w1r*w3r**2 + w2r*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w3i**2 +&
        2*rhogeq*cs**2*k**3*w3i**2*Kdrag*vgsol*w2r*w1r*w1i - 2*w2i*w3i**2*cs**2*k**2*rhogeq*rhogsol*w2r**2*w1r*w1i +&
        w1i**2*rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r*w1r - w1i**2*rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r*w1r +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i*w1i*w3r**2 - w3r*rhogeq*k*Kdrag*vdsol*w2r*w2i**2*w1i*w3i**2 +&
        4*w3r*w2i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r*w3i*w1i - 2*rhogeq*cs**2*k**3*w3i**2*Kdrag*vdsol*w2r*w1r*w1i -&
        2*w2r*rhogeq**2*cs**2*k**3*w1r*vgsol*w3i*w1i*w2i**2 - 4*w2i*w3r**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i*w3i**2 +&
        3*w3r*w2r*Kdrag*rhogeq*k*vgsol*w2i**2*w1i*w3i**2 - w2i*w3r*Kdrag*k*rhogeq*vgsol*w2r**2*w1r*w3i**2 -&
        w2r*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w3i**2 + rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r*w1r*w3r**2 +&
        2*w3r*w3i*cs**4*k**4*rhogeq*rhogsol*w2i*w2r**2 - 2*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r*w2i**2 -&
        2*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r**3 - 4*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w1r*w2r**2 +&
        2*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r**3 + 2*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r*w2i**2 +&
        4*w3r*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w1r*w2r**2 + 2*w3r*w3i*cs**4*k**4*rhogeq*rhogsol*w2i**3 +&
        4*w2r*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i*w1r*w3i**2 - 4*rhogeq**2*k*vgsol*w2r*w2i**2*w1r*w3i**3*w1i -&
        w3r*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i**4 - 2*w3r*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i**2*w2r**2 -&
        2*w3r**3*Kdrag*rhogsol*w2r**2*w1i*w2i**2 - rhogeq*w3i**3*w2i*w1r*rhogsol*w2r**2*w1i**2 +&
        2*w1r**2*rhogsol*k**4*cs**4*rhogeq*w2r*w2i*w3i - rhogeq*w3i**3*k*Kdrag*vdsol*w2r**2*w1i**2 +&
        w1i**2*rhogeq**2*k*vgsol*w2r*w2i**2*w1r*w3i**2 + w1i**3*w2i*rhogeq**2*k*vgsol*w2r**2*w3i**2 -&
        2*w2r*rhogeq*w1i**3*rhogsol*k**2*cs**2*w3i**2*w2i + w1i**2*rhogeq**2*w2r**3*k*w1r*vgsol*w3i**2 -&
        2*w2r*cs**2*k**2*rhogeq*rhogsol*w1r**2*w3i**2*w1i*w2i - 2*Kdrag*k*rhogeq*vgsol*w2r**3*w1i*w1r*w3i**2 -&
        2*w2r*Kdrag*rhogeq*k*vgsol*w2i**2*w1i*w1r*w3i**2 + w1r**3*rhogeq**2*k*vgsol*w2r*w2i**2*w3i**2 -&
        2*w1r**3*w2r*rhogeq**2*cs**2*k**3*vgsol*w2i*w3i - 2*w3i**2*cs**2*k**3*rhogeq**2*vgsol*w2r**2*w2i**2 -&
        w1r**3*w2r*rhogeq**2*cs**2*k**3*vgsol*w3i**2 - w3i*rhogsol*k**4*cs**4*rhogeq*w2r*w1i*w2i**2 -&
        w3i*rhogsol*k**4*cs**4*rhogeq*w2r**3*w1i - w2i*rhogeq*k*Kdrag*vdsol*w2r**2*w1r**2*w3i**2 +&
        w1r**3*rhogeq**2*w2r**3*k*vgsol*w3i**2 - w3i**2*rhogeq*k*Kdrag*vdsol*w2r**2*w1i**2*w2i -&
        w3i**2*cs**2*k**3*rhogeq**2*vgsol*w2r**4 + 2*rhogeq**2*cs**2*k**3*w3i*w1i*vgsol*w2r**2*w2i**2 -&
        rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r**3*w1r - 2*w3r*w3i**2*Kdrag*rhogsol*w2r**2*w1i*w2i**2 -&
        rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r*w1r*w2i**2 + rhogeq**2*cs**2*k**3*w3i*w1i*vgsol*w2r**4 +&
        2*w2r*w3r**2*cs**2*k**2*rhogeq*rhogsol*w1r**2*w2i**2 + 2*w2r**3*w3r**2*cs**2*k**2*rhogeq*rhogsol*w1r**2 -&
        2*rhogeq*w2r*w3i**2*cs**2*k**2*rhogsol*w1r**2*w2i**2 - 2*rhogeq*w2r**3*w3i**2*cs**2*k**2*rhogsol*w1r**2 +&
        rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r*w1r*w2i**2 + rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r**3*w1r +&
        w2r*w3r*rhogeq*cs**4*k**4*w1r*rhogsol*w2i**2 + w2r**3*w3r*rhogeq*cs**4*k**4*w1r*rhogsol +&
        w2r*w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol*w2i**2 + w2r**3*w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol +&
        w2r*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2i**2 + w2r**3*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol -&
        w2r**3*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol - w2r*w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol*w2i**2 -&
        w2r**3*w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol - w2r*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2i**2 +&
        4*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w2i**3*w1i**2 + 4*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w2i*w1i**2*w2r**2 -&
        w3r*rhogeq**2*cs**2*k**3*w2r**4*w1r*vgsol + 2*w3r**2*rhogeq**2*cs**2*k**3*w2r**2*vgsol*w2i**2 +&
        w3r**2*rhogeq**2*cs**2*k**3*w2r**4*vgsol + w3r*rhogeq*cs**2*k**2*w2r**4*w1r**2*rhogsol +&
        2*w3r*rhogeq*cs**2*k**2*w1i**2*w2r**2*rhogsol*w2i**2 - 4*w3r*w3i*cs**4*k**4*rhogeq*rhogsol*w1i*w2i**2 +&
        w3r*rhogeq*cs**2*k**2*w1i**2*w2r**4*rhogsol + 2*w3r*rhogeq*cs**2*k**2*w2r**2*w1r**2*rhogsol*w2i**2 -&
        2*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w2r**4*w1i - 2*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w1i*w2i**4 -&
        4*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w1i*w2i**2*w2r**2 - w3i**3*Kdrag*rhogeq*k*vgsol*w2i**4 +&
        2*w3i**3*Kdrag*rhogsol*w2r**2*w1r*w2i**2 + w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i**3*w1i +&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i*w1i*w2r**2 - w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w1r*w2r**2 -&
        w3r*rhogeq*cs**4*k**4*rhogsol*w2i**3*w1i + w3r*rhogeq*cs**2*k**3*w2i**3*Kdrag*vdsol*w1r +&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w1r*w2r**2 - w3r*rhogeq*cs**2*k**3*w2i**3*Kdrag*vgsol*w1r +&
        w3r*rhogeq*cs**2*k**2*w1i**2*rhogsol*w2i**4 + w3r**2*rhogeq**2*cs**2*k**3*vgsol*w2i**4 -&
        w3i*cs**4*k**4*rhogeq*rhogsol*w2i**3*w1r - w3i*cs**4*k**4*rhogeq*rhogsol*w2i*w1r*w2r**2 -&
        w3r*rhogeq*cs**4*k**4*rhogsol*w2i*w1i*w2r**2 + w3r*rhogeq*cs**2*k**2*w1r**2*rhogsol*w2i**4 +&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**4*w1i + w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**3*w1r**2 +&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1r**2*w2r**2 - w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**3*w1i**2 -&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i**2*w2r**2 - w3i**2*rhogeq**2*cs**2*k**3*vgsol*w2i**4 +&
        Kdrag*w2i**4*vgsol*k*rhogeq*w1i*w3r**2 - w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i**3*w1i -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i*w1i*w2r**2 - 2*w3i**3*cs**2*k**3*rhogeq**2*vgsol*w2r*w1r*w2i -&
        rhogeq*w3i**3*w2i*w1r**3*rhogsol*w2r**2 + 2*w3i*w3r**2*Kdrag*rhogsol*w2r**2*w1r*w2i**2 -&
        w3i*w3r**2*Kdrag*rhogeq*k*vgsol*w2i**4 + Kdrag*w2i**4*k*rhogeq*vgsol*w3i**2*w1i +&
        2*w1r**2*rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r**2 - 2*w1r**2*rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r**2 +&
        w1r**4*w2r*cs**2*k**2*rhogeq*rhogsol*w3i**2 + w3r*rhogeq**2*cs**2*k**3*w2r**3*vgsol*w3i**2 +&
        2*w2i*w3r**2*rhogeq**2*k*vgsol*w2r**2*w1i*w3i**2 - w2i*w3r**4*Kdrag*rhogeq*k*vgsol*w2r**2 +&
        w2i*w3r**4*rhogeq**2*k*vgsol*w2r**2*w1i - 2*w2i*w3r**4*cs**2*k**2*rhogeq*rhogsol*w2r*w1i -&
        w3r*Kdrag*rhogsol*w2i**3*w1r**2*w3i**2 - Kdrag*w3r**3*w2r**4*rhogsol*w1i +&
        w3r*Kdrag*rhogsol*w2i**3*w1i**2*w3i**2 + 2*Kdrag*k*rhogeq*vgsol*w2i**3*w1r**2*w3i**2 -&
        2*w1i**2*w2r*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i*w3i - w1i**2*w2r*rhogeq**2*cs**2*k**3*w1r*vgsol*w3i**2 -&
        w3r**3*w2r**3*k*Kdrag**2*vdsol + 2*Kdrag*k*rhogeq*vgsol*w2i**3*w1r**2*w3r**2 +&
        w1r**2*w2i*rhogeq**2*k*vgsol*w2r**2*w1i*w3i**2 + k*Kdrag**2*vgsol*w2i**3*w1i*w3i**2 +&
        k*Kdrag**2*vgsol*w2i**3*w1i*w3r**2 + w3r**3*w2r**3*k*Kdrag**2*vgsol +&
        2*w1i**2*w2r*cs**2*k**2*rhogeq*rhogsol*w1r**2*w3i**2 - w3r**3*rhogeq*w2r**4*w1r**2*rhogsol -&
        rhogeq*w1i**2*w3r**3*w2r**4*rhogsol - w3r**4*rhogeq*w2r**3*w1r**2*rhogsol - k*Kdrag**2*vdsol*w2i**3*w1i*w3i**2&
        - w3r**4*Kdrag*rhogsol*w2r**3*w1i - w3r**3*Kdrag*rhogsol*w2i**3*w1r**2 + w3r**3*Kdrag*rhogsol*w2i**3*w1i**2 -&
        k*Kdrag**2*vdsol*w2i**3*w1i*w3r**2 - rhogeq**2*w3i**2*k*w1r**4*vgsol*w2r**2 +&
        Kdrag*w2r**3*rhogsol*w1i**2*w3i**3 + 2*w3r**2*rhogeq**2*w2r**3*k*w1r*vgsol*w3i**2 +&
        w2i*w3r**3*Kdrag*rhogsol*w2r**2*w1i**2 + w3r**4*rhogeq*w2r**2*w1i*k*Kdrag*vgsol -&
        w2i*w3r**3*Kdrag*rhogsol*w2r**2*w1r**2 + w3r**4*rhogeq**2*w2r**3*k*w1r*vgsol -&
        2*w3r*rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol*w3i**2 - 2*w3r**2*Kdrag*rhogsol*w2r**3*w1i*w3i**2 +&
        2*w2r*rhogeq**2*cs**2*k**3*w1r*vgsol*w3i**3*w1i + w3r**4*rhogeq**2*cs**2*k**3*w2r**2*vgsol +&
        2*w3r*rhogeq*cs**2*k**2*w2r**2*w1r**2*rhogsol*w3i**2 + 2*w3r**2*rhogeq**2*cs**2*k**3*w2r**2*vgsol*w3i**2 +&
        rhogeq**2*cs**2*k**3*w2r**2*vgsol*w1r**2*w3i*w1i - w1i**3*w3i*rhogsol*k**4*cs**4*rhogeq*w2r -&
        2*w2r*cs**2*k**2*rhogeq*rhogsol*w1r**2*w3i**3*w1i + 3*w3r*Kdrag*k*rhogeq*vgsol*w2r**3*w1i*w3i**2 +&
        2*rhogeq*w3i**3*cs**2*k**2*rhogsol*w1r*w2r**2*w1i + w3i**3*Kdrag*rhogsol*w2r**4*w1r -&
        2*w3r**3*rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol + 2*w3r**3*rhogeq*cs**2*k**2*w2r**2*w1r**2*rhogsol -&
        rhogeq*w2r**3*w3i**4*w1r**2*rhogsol - w3r**3*Kdrag**2*k*vgsol*w2r**2*w1r + w3r**3*Kdrag**2*k*vdsol*w2r**2*w1r +&
        3*w3r**3*Kdrag*k*rhogeq*vgsol*w2r**3*w1i + 2*rhogeq*cs**4*k**4*w2r*rhogsol*w2i*w1i**2*w3i -&
        w3r*Kdrag**2*k*vgsol*w2r**2*w1r*w3i**2 + rhogeq*w2r**4*w3i**4*w1r*rhogsol - rhogeq**2*w2r**4*w3i**4*k*vgsol +&
        w3r**4*rhogeq*w2r**4*w1r*rhogsol + w3r*Kdrag**2*k*vdsol*w2r**2*w1r*w3i**2 -&
        w3r*Kdrag*rhogsol*k**2*cs**2*w1i*w2r**2*w3i**2 - w3r**4*rhogeq**2*w2r**4*vgsol*k +&
        w3r**3*rhogeq**2*cs**2*k**3*w2r**3*vgsol - w3r**3*Kdrag*rhogsol*k**2*cs**2*w1i*w2r**2 -&
        2*w2i*w3r**2*Kdrag*rhogeq*k*vgsol*w2r**2*w3i**2 + 2*rhogeq*cs**4*k**4*w1r*rhogsol*w2r**2*w3i*w1i -&
        w3r**4*rhogeq*w1i**2*w2r**3*rhogsol - 2*w2i*w3r**3*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r +&
        4*w3r**2*rhogeq*rhogsol*w2r**2*w2i**2*w1r*w3i**2 + w2i*w3r**3*Kdrag*cs**2*k**2*rhogsol*w2r**2 +&
        w2i*w3r*Kdrag*cs**2*k**2*rhogsol*w2r**2*w3i**2 - 2*rhogeq**2*cs**2*k**3*w2r**3*w1r*vgsol*w3i*w1i -&
        Kdrag*w2r**3*rhogsol*w3i**3*w1r**2 - w1r**3*rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r -&
        w3r**4*rhogeq*rhogsol*w2r*w2i**2*w1r**2 + w2i*w3r*Kdrag*rhogsol*w2r**2*w1i**2*w3i**2 +&
        2*w3r**4*rhogeq*rhogsol*w2r**2*w2i**2*w1r + 2*w2i*w3r**3*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r -&
        rhogeq*w1i**2*w2r**3*rhogsol*w3i**4 - Kdrag*w2r**3*rhogsol*w1i*w3i**4 -&
        w2i*w3r*Kdrag*rhogsol*w2r**2*w1r**2*w3i**2 + w2i*w3r**3*rhogeq*k*Kdrag*vdsol*w2r**2*w1r +&
        2*w3r**2*rhogeq*w2r**2*w1i*k*Kdrag*vgsol*w3i**2 - 4*w3r**2*rhogeq**2*k*vgsol*w2i**2*w2r**2*w3i**2 +&
        2*w3r**2*rhogeq**2*k*vgsol*w2r*w2i**2*w1r*w3i**2 - 2*w3r*rhogeq*rhogsol*w2i**2*w1r**2*w2r**2*w3i**2 -&
        w3r**4*rhogeq*rhogsol*w2r*w2i**2*w1i**2 + w3i**3*Kdrag*rhogsol*w2i**4*w1r -&
        2*w3r*rhogeq*rhogsol*w1i**2*w2r**2*w2i**2*w3i**2 - 2*w3r**2*rhogeq*rhogsol*w2r*w2i**2*w1i**2*w3i**2 +&
        w3r**4*rhogeq*rhogsol*w2i**4*w1r + w3i**3*vdsol*k*Kdrag**2*w2i**3 +&
        rhogeq*w1i*w3i**2*k*Kdrag*vdsol*w2r**2*w1r**2 + w3r**3*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w2r -&
        w3i**3*vgsol*k*Kdrag**2*w2i**3 + w3r**4*rhogeq**2*k*vgsol*w2r*w2i**2*w1r - w3r**3*Kdrag*rhogsol*w2i**4*w1i -&
        2*w3r**2*rhogeq*rhogsol*w2r*w2i**2*w1r**2*w3i**2 - w3r**3*rhogeq*k*Kdrag*vdsol*w2r*w2i**2*w1i -&
        w3r**3*rhogeq*rhogsol*w2i**4*w1i**2 + 2*w3i**2*rhogsol*k**4*cs**4*rhogeq*w2r*w1i**2 +&
        2*w2i*w3r**2*Kdrag*rhogsol*w2r**2*w1r*w3i**2 - w1r*Kdrag*rhogsol*k**2*cs**2*w2i**3*w3i**2 +&
        2*w3r*rhogeq**2*k*vgsol*w2r**2*w2i**2*w1r*w3i**2 - w3r**3*rhogeq*rhogsol*w2i**4*w1r**2 -&
        w1r*Kdrag*rhogsol*k**2*cs**2*w2i**3*w3r**2 - 2*w3r*rhogeq**2*k*vgsol*w2r*w2i**2*w1r**2*w3i**2 -&
        w3i**4*rhogeq**2*vgsol*k*w2i**4 + w3r**4*Kdrag*rhogsol*w2i**3*w1r - w3r**4*rhogeq**2*vgsol*k*w2i**4 +&
        w3i**4*Kdrag*rhogsol*w2i**3*w1r - w3r*rhogeq*cs**4*k**4*w2r**2*rhogsol*w3i**2 +&
        2*w3r*rhogeq**2*k*vgsol*w2r*w2i**2*w1i**2*w3i**2 + w3i**4*rhogeq*rhogsol*w2i**4*w1r +&
        3*w3r**3*w2r*Kdrag*rhogeq*k*vgsol*w2i**2*w1i + 2*w3r**3*rhogeq**2*k*vgsol*w2r**2*w2i**2*w1r -&
        2*w3r**3*rhogeq**2*k*vgsol*w2r*w2i**2*w1r**2 + 2*w3r**3*rhogeq**2*k*vgsol*w2r*w2i**2*w1i**2 -&
        2*w3r**4*rhogeq**2*k*vgsol*w2i**2*w2r**2 - 2*w3r**3*rhogeq*rhogsol*w2i**2*w1r**2*w2r**2 +&
        w3r*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w2r*w3i**2 - 2*w3r**3*rhogeq*rhogsol*w1i**2*w2r**2*w2i**2 +&
        w1r*w3r**4*cs**2*k**2*rhogeq*rhogsol*w2i**2 + w1r*w3i**4*cs**2*k**2*rhogeq*rhogsol*w2i**2 +&
        2*rhogeq**2*w1i**2*k*vgsol*w3r*w2r**3*w3i**2 + w2i*w3r**4*Kdrag*rhogsol*w2r**2*w1r -&
        2*w3r**2*rhogeq*w2r**3*w1r**2*rhogsol*w3i**2 - w3r**3*rhogeq*cs**4*k**4*w2r**2*rhogsol +&
        w1r**3*rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r - 2*w3r**3*rhogeq**2*w1r**2*k*vgsol*w2r**3 +&
        2*rhogeq**2*w1i**2*k*vgsol*w3r**3*w2r**3 - w2i*w3r**3*Kdrag*k*rhogeq*vgsol*w2r**2*w1r -&
        rhogeq*w1i*w3r**3*w2r**3*k*Kdrag*vdsol - 2*w3r*rhogeq**2*w1r**2*k*vgsol*w2r**3*w3i**2 +&
        w2r*w3r**3*rhogeq**2*cs**2*k**3*w1r**2*vgsol - w2r*w3r**3*rhogeq**2*cs**2*k**3*w1i**2*vgsol -&
        w1i**3*rhogeq*w2r**2*k*Kdrag*vgsol*w3i**2 + 2*w1r*w3i**2*cs**2*k**2*rhogeq*rhogsol*w2i**2*w3r**2 +&
        w2r*rhogeq*w1i**2*rhogsol*k**2*cs**2*w3r**4 - rhogeq*w3i**4*cs**2*k**2*rhogsol*w1r*w2r**2 -&
        2*w2r*w3r**2*Kdrag*rhogsol*w1i*w2i**2*w3i**2 - w2r*w3r**4*rhogeq**2*cs**2*k**3*w1r*vgsol -&
        rhogeq*w1i*w3r*w2r**3*k*Kdrag*vdsol*w3i**2 + w2r*w3r**3*rhogeq*cs**4*k**4*w1r*rhogsol -&
        Kdrag*w3r*w2r**4*rhogsol*w1i*w3i**2 - w3r*w2r**3*k*Kdrag**2*vdsol*w3i**2 + w3r*w2r**3*k*Kdrag**2*vgsol*w3i**2 +&
        rhogeq*rhogsol*w2r*w2i**2*w1i**3*w3i**3 - 4*rhogeq**2*w2r**3*k*w1r*vgsol*w3i**3*w1i +&
        w2r**2*k*Kdrag**2*vdsol*w3i**2*w1i**2 - w3r*rhogeq*w2r**4*w1r**2*rhogsol*w3i**2 -&
        rhogeq*w1i**2*w3r*w2r**4*rhogsol*w3i**2 + w3r*rhogeq**2*w2r**4*k*w1r*vgsol*w3i**2 -&
        w2r**2*k*Kdrag**2*vgsol*w3i**2*w1i**2 - Kdrag*rhogsol*k**2*cs**2*w1i*w2r**3*w3r**2 +&
        w2r*w3r**3*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol - w2r*w3r**3*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol +&
        w2r*w3r**4*cs**2*k**2*rhogeq*rhogsol*w1r**2 + 2*w2r*rhogeq*w1i**2*rhogsol*k**2*cs**2*w3r**2*w3i**2 -&
        w3r**4*rhogeq*rhogsol*k**2*cs**2*w2r**2*w1r - rhogeq**2*w1i**4*k*vgsol*w3i**2*w2r**2 -&
        Kdrag*rhogsol*k**2*cs**2*w1i*w2r**3*w3i**2 + w2r*w3r*rhogeq**2*cs**2*k**3*w1r**2*vgsol*w3i**2 -&
        w2r*w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol*w3i**2 - w2r*w3r**4*Kdrag*rhogsol*w1i*w2i**2 +&
        rhogeq*w2r*w3i**3*w1r*k*Kdrag*vdsol*w2i**2 - 2*w2r*w3r**2*rhogeq**2*cs**2*k**3*w1r*vgsol*w3i**2 +&
        rhogeq**2*w3i**3*w1i**3*k*vgsol*w2r**2 + rhogeq*w2r*w3i**4*cs**2*k**2*rhogsol*w1r**2 +&
        w2r*w3r*rhogeq*cs**4*k**4*w1r*rhogsol*w3i**2 - 2*w2r*rhogeq*w1i**3*rhogsol*k**2*cs**2*w3i**3 -&
        rhogeq*cs**2*k**3*w3i**3*Kdrag*vdsol*w2r**2 + rhogeq*cs**4*k**4*w2r*rhogsol*w2i**2*w3i**2 +&
        w3r**3*rhogeq**2*w2r**4*k*w1r*vgsol - w1r**2*rhogeq*w2r**2*w1i*k*Kdrag*vgsol*w3i**2 -&
        rhogeq*cs**4*k**4*w2r*rhogsol*w2i**2*w3r**2 + 2*rhogeq**2*cs**2*k**3*w2r**2*vgsol*w1r**2*w3i**2 +&
        rhogeq*cs**2*k**3*w3i**3*Kdrag*vdsol*w2r*w1r + 2*rhogeq**2*w3i**3*w2r**2*w1i*k*vgsol*w2i**2 -&
        rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r**2*w3r**2 - w1r**2*w3i*rhogsol*k**4*cs**4*rhogeq*w2r*w1i -&
        2*rhogeq**2*cs**2*k**3*w1i**2*vgsol*w2r**2*w3i**2 + rhogeq*cs**4*k**4*w2r**3*rhogsol*w3i**2 -&
        rhogeq*cs**4*k**4*w2r**3*rhogsol*w3r**2 + 2*w2r*w3r**2*cs**2*k**2*rhogeq*rhogsol*w1r**2*w3i**2 +&
        2*w3r**2*rhogeq*w2r**4*w1r*rhogsol*w3i**2 - 2*w3r**2*rhogeq**2*w2r**4*vgsol*k*w3i**2 -&
        2*rhogeq*w3i**2*cs**2*k**2*rhogsol*w1r*w2r**2*w3r**2 + rhogeq*w2i*w3i**3*w1i*k*Kdrag*vdsol*w2r**2 +&
        rhogeq*cs**2*k**3*w3i**3*Kdrag*vgsol*w2r**2 + rhogeq*cs**2*k**2*w2r**4*w1r*rhogsol*w3i**2 -&
        rhogeq*cs**2*k**2*w2r**4*w1r*rhogsol*w3r**2 + 2*rhogeq**2*cs**2*k**3*w2r**2*vgsol*w2i*w1i*w3i**2 -&
        2*w2i*w3r*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r*w3i**2 + rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r**2*w3r**2 -&
        rhogeq*cs**2*k**3*w3i**3*Kdrag*vgsol*w2r*w1r + w2r*w3r**3*w2i**2*k*Kdrag**2*vgsol +&
        rhogeq**2*w3i*w2r**4*w1i*k*vgsol*w3r**2 + rhogeq**2*w2r**3*w3i**4*k*w1r*vgsol +&
        2*rhogeq*w2r**2*w3i**4*w1r*rhogsol*w2i**2 + rhogeq**2*w3i**3*w2r**4*w1i*k*vgsol -&
        2*rhogeq*w2r*w3i**4*cs**2*k**2*rhogsol*w2i*w1i - w2r*w3r**3*w2i**2*k*Kdrag**2*vdsol -&
        rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2r**2*w3r**2 - rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w2r**2*w3i**2 +&
        rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w2r**2*w3r**2 + rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2r**2*w3i**2 +&
        w2r*w3r*w2i**2*k*Kdrag**2*vgsol*w3i**2 - rhogeq*cs**4*k**4*w1r*rhogsol*w2r**2*w3i**2 +&
        3*rhogeq*cs**4*k**4*w1r*rhogsol*w2r**2*w3r**2 + rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2r**2*w3r**2 -&
        rhogeq*w2r*w3i**4*w1r**2*rhogsol*w2i**2 - 4*rhogeq*cs**4*k**4*w2r*rhogsol*w2i*w1i*w3i**2 -&
        2*rhogeq**2*cs**2*k**3*w2r**3*w1r*vgsol*w3r**2 + rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w2r**2*w3i**2 +&
        Kdrag*k*rhogeq*vgsol*w2r**4*w1i*w3i**2 + Kdrag*k*rhogeq*vgsol*w2r**4*w1i*w3r**2 -&
        rhogeq**2*cs**2*k**3*w3i**3*w2i*vgsol*w2r**2 + w1i**3*w3r**3*rhogeq*rhogsol*w2i**3 +&
        2*rhogeq**2*w2i*k*vgsol*w3i**3*w1r**2*w2r**2 - w2r*w3r*w2i**2*k*Kdrag**2*vdsol*w3i**2 -&
        rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w2r**2*w3r**2 - rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2r**2*w3i**2 -&
        w3i**4*cs**2*k**3*rhogeq**2*vgsol*w2r*w1r + 2*rhogeq*cs**2*k**2*w2r**2*w1r*rhogsol*w2i**2*w3i**2 -&
        2*rhogeq*cs**2*k**2*w2r**2*w1r*rhogsol*w2i**2*w3r**2 + 2*rhogeq**2*w3i*w2r**2*w1i*k*vgsol*w2i**2*w3r**2 -&
        w3i**3*vgsol*k*Kdrag**2*w2i*w2r**2 + w3i**3*rhogsol*k**2*cs**2*Kdrag*w2r*w2i**2 +&
        w3i*rhogsol*k**2*cs**2*Kdrag*w2r**3*w3r**2 + w3i**3*rhogsol*k**2*cs**2*Kdrag*w2r**3 +&
        w3i**4*cs**2*k**3*rhogeq**2*vgsol*w2r**2 - rhogeq**2*cs**2*k**3*w3i*w2i*vgsol*w2r**2*w3r**2 +&
        2*rhogeq**2*w2i*k*vgsol*w3i*w1r**2*w2r**2*w3r**2 - w3i*vgsol*k*Kdrag**2*w2i*w2r**2*w3r**2 +&
        w3i*Kdrag*rhogsol*w2r**4*w1r*w3r**2 - 2*w3r*rhogeq*cs**2*k**2*w1i**3*w2r**2*rhogsol*w3i +&
        2*w3r*rhogeq**2*cs**2*k**3*w2r**3*vgsol*w3i*w1i - 2*w3r*w1r**2*w2i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r -&
        2*w3r*w1i**3*w2i*rhogeq*cs**2*k**2*rhogsol*w2r**2 + 2*w3r*w1r**2*w2i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r +&
        2*w3r*w1i**2*rhogeq*cs**2*k**2*w2r**2*w1r**2*rhogsol + w3r*rhogeq*w1i**3*w3i**2*rhogsol*w2i*w2r**2 +&
        w3i*rhogsol*k**2*cs**2*Kdrag*w2r*w2i**2*w3r**2 + 2*w3r*w2r*w1i**3*rhogeq**2*cs**2*k**3*vgsol*w2i +&
        2*w3r*w2r*rhogeq**2*cs**2*k**3*w1i**3*vgsol*w3i + w3r*w1i**4*rhogeq*cs**2*k**2*w2r**2*rhogsol +&
        2*w3r*w2r*w1r**2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i + w3r*w1i**2*rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol +&
        w3r*rhogeq*w2r**3*w1r**3*rhogsol*w3i**2 - 2*w3r*w2r*rhogeq*cs**2*k**3*w1i**2*Kdrag*vdsol*w3i +&
        2*w3r*w2r*rhogeq*cs**2*k**3*w1i**2*Kdrag*vgsol*w3i - 2*w3r*w2r*cs**2*k**2*rhogeq*rhogsol*w1r**3*w3i**2 +&
        2*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2r**2*w1r - 2*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2r**2*w1r +&
        2*w3r*rhogeq**2*cs**2*k**3*w2r**2*vgsol*w2i*w1i*w1r - 2*w3r*w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**3*w1r +&
        4*w3r*rhogeq*rhogsol*k**2*cs**2*w2r**3*w1r*w3i*w1i + 2*w3r*w2r*Kdrag*rhogsol*w1i*w2i**2*w3i**2*w1r -&
        2*w3r*w1r**2*rhogeq*cs**4*k**4*w2r**2*rhogsol + w3r*rhogeq*rhogsol*w2r*w2i**2*w1i**2*w3i**2*w1r -&
        w3r*w1r**2*w2r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol + w3r*w1r**2*w2r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol +&
        2*w3r*w2r*rhogeq**2*cs**2*k**3*w1r**2*vgsol*w3i*w1i - 2*w3r*w2r*rhogeq*w1i**2*rhogsol*k**2*cs**2*w3i**2*w1r +&
        w3r*rhogeq*rhogsol*w2r*w2i**2*w1r**3*w3i**2 - 4*w3r*rhogsol*k**4*cs**4*rhogeq*w2r*w2i*w3i*w1r +&
        w3r*w1r**3*rhogeq**2*cs**2*k**3*w2r**2*vgsol + w3r*rhogeq**2*w3i**2*k*w1r**3*vgsol*w2r**2 -&
        4*w3r*rhogeq**2*cs**2*k**3*w2r**2*w1r*vgsol*w3i*w1i - 2*w3r*rhogeq*cs**2*k**2*w2r**2*w1r**2*rhogsol*w3i*w1i +&
        2*w3r*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w2r*w3i*w1i + 4*w3r*w2i*cs**2*k**2*rhogeq*rhogsol*w2r*w1i*w3i**2*w1r -&
        2*w3r*w1i**2*rhogeq*cs**2*k**2*rhogsol*w2r*w2i**2*w1r - 2*w3r*w1i**2*w2i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2r +&
        2*w3r*w1i**2*w2i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r + w3r*rhogeq*w1i**2*w2r**3*rhogsol*w3i**2*w1r -&
        2*w3r*rhogeq**2*cs**2*k**3*w3i*w2i*vgsol*w2r**2*w1r + w3r**3*rhogeq*w2r**3*w1r**3*rhogsol -&
        2*w3r**3*rhogeq*w2r**2*w1i*k*Kdrag*vgsol*w1r - w3i*rhogsol*k**4*cs**4*rhogeq*w2r*w1i*w3r**2 +&
        w3i*vdsol*k*Kdrag**2*w2i*w2r**2*w3r**2 + 2*rhogeq*w2r**2*w3i**3*w1r**2*k*Kdrag*vgsol -&
        4*w3r*rhogeq**2*cs**2*k**3*w2r*w1i**2*vgsol*w2i*w3i - w3r*w1r**4*w2r*rhogeq**2*cs**2*k**3*vgsol +&
        4*w3r*w2i**2*cs**2*k**2*rhogeq*rhogsol*w2r*w1i*w1r*w3i + w3r*w1r**3*w2r*rhogeq*cs**4*k**4*rhogsol -&
        w3i**3*vgsol*k*rhogeq*Kdrag*w2r**4 - 2*w3r*rhogeq*w2r**2*w1i*k*Kdrag*vgsol*w3i**2*w1r +&
        w3r*rhogeq**2*w3i**2*k*w1r*vgsol*w2r**2*w1i**2 + 4*w3r*rhogeq**2*cs**2*k**3*w2r*w1r**2*vgsol*w2i*w3i +&
        2*w3r*cs**2*k**3*rhogeq**2*w2r*vgsol*w2i*w1i*w3i**2 - 2*w3r*w2r**3*rhogeq*w1i**2*rhogsol*k**2*cs**2*w1r -&
        2*w3r*w2r**3*cs**2*k**2*rhogeq*rhogsol*w1r**3 + rhogeq*w2r**3*w1r**2*rhogsol*w3i**3*w1i -&
        2*w3r*rhogeq*cs**2*k**3*w3i*Kdrag*vdsol*w2r*w1r**2 - 4*w3r*w2i*rhogeq**2*k*vgsol*w2r**2*w1i*w3i**2*w1r -&
        2*w3r*w1r**3*rhogeq*cs**2*k**2*rhogsol*w2r*w2i**2 - 2*w3r*w2i*w1i*rhogeq**2*cs**2*k**3*w1r*vgsol*w3i**2 +&
        w3i**3*vdsol*k*Kdrag**2*w2i*w2r**2 + 2*rhogeq*w2r**2*w3i*w1r**2*k*Kdrag*vgsol*w3r**2 -&
        2*w3r*w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1r**3 - 2*w3r*w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i**2*w1r +&
        4*w3r*rhogeq**2*cs**2*k**3*w3i*w1i*vgsol*w2i**2*w1r + 2*w3r*rhogeq*cs**2*k**3*w3i*Kdrag*vgsol*w2r*w1r**2 +&
        2*w3r*Kdrag*rhogsol*w2r**3*w1i*w3i**2*w1r - 2*w3r*w1i**2*w2r*rhogeq**2*cs**2*k**3*w1r**2*vgsol +&
        w3r*w1i**2*w2r*rhogeq*cs**4*k**4*w1r*rhogsol + w3r*w1i**3*w2r*rhogeq*cs**2*k**3*Kdrag*vdsol -&
        w3r*w1i**3*w2r*rhogeq*cs**2*k**3*Kdrag*vgsol - w3r*w1i**4*w2r*rhogeq**2*cs**2*k**3*vgsol +&
        w3r*rhogeq*w1i*w3i**2*rhogsol*w2i*w2r**2*w1r**2 - 2*w3r*w1r**2*w2i*rhogeq*cs**2*k**2*rhogsol*w2r**2*w1i +&
        Kdrag*w2i**2*rhogsol*w1i**2*w2r*w3i**3 + Kdrag*w2i**2*rhogsol*w1i**2*w2r*w3i*w3r**2 -&
        Kdrag*w2i**2*rhogsol*w3i*w2r*w1r**2*w3r**2 - Kdrag*w2i**2*rhogsol*w3i**3*w2r*w1r**2 -&
        Kdrag*w2r**2*k*rhogeq*vgsol*w2i*w3i**4 - w3i*vgsol*k*rhogeq*Kdrag*w2r**4*w3r**2 +&
        2*w1r**3*rhogeq*cs**2*k**2*rhogsol*w2i**2*w3i**2 + Kdrag*w2r**3*rhogsol*w1i**2*w3i*w3r**2 -&
        w1r**3*w3i**3*rhogeq*rhogsol*w2i**3 - 4*w3r**3*w2i*rhogeq**2*k*vgsol*w2r**2*w1i*w1r -&
        w3i**3*rhogsol*k**4*cs**4*rhogeq*w2r*w1i + rhogeq*cs**2*k**2*w1i**2*w2r*rhogsol*w3i**4 -&
        w2r**3*k*Kdrag**2*vgsol*w1r*w3r**2 - 2*w1i**2*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w3i**2 -&
        2*rhogeq**2*w1i**2*w2r**2*k*vgsol*w2i*w3i**3 - w2r**2*k*Kdrag**2*vdsol*w3i**3*w1i +&
        w2i**2*k*Kdrag**2*vdsol*w2r*w1r*w3i**2 + w1i**4*w2r*rhogeq*rhogsol*k**2*cs**2*w3i**2 +&
        w2r**2*k*Kdrag**2*vgsol*w3i**3*w1i + w2i**2*k*Kdrag**2*vdsol*w2r*w1r*w3r**2 -&
        2*rhogeq**2*w1i**2*w2r**2*k*vgsol*w2i*w3i*w3r**2 + w2r**3*k*Kdrag**2*vdsol*w1r*w3i**2 -&
        w1r**2*vdsol*Kdrag*k*rhogeq*w2i**3*w3i**2 - w1r**2*vdsol*Kdrag*k*rhogeq*w2i**3*w3r**2 +&
        w2r**3*k*Kdrag**2*vdsol*w1r*w3r**2 - Kdrag*w2r**2*rhogsol*k**2*cs**2*w2i*w1r*w3i**2 -&
        w1r**2*Kdrag*k*rhogeq*vgsol*w1i*w2i**2*w3i**2 - w2i**2*k*Kdrag**2*vgsol*w2r*w1r*w3i**2 -&
        w2i**2*k*Kdrag**2*vgsol*w2r*w1r*w3r**2 - w1r**2*Kdrag*k*rhogeq*vgsol*w1i*w2i**2*w3r**2 -&
        Kdrag*w2r**2*rhogsol*k**2*cs**2*w2i*w1r*w3r**2 + w2i**2*k*Kdrag**2*vgsol*w3i**3*w1i +&
        2*w3r*w3i*cs**4*k**4*rhogeq*rhogsol*w2i*w1r**2 - w2i**2*k*Kdrag**2*vdsol*w3i**3*w1i +&
        w2i**2*k*Kdrag**2*vdsol*w3r**3*w1r - 2*w3r**3*w2r*cs**2*k**2*rhogeq*rhogsol*w1r**3 +&
        2*w3r**3*w2r*Kdrag*rhogsol*w1i*w2i**2*w1r - w3i**4*rhogeq**2*cs**2*k**3*vgsol*w2i**2 -&
        Kdrag*w2i**2*k*rhogeq*vgsol*w2r*w3i**3*w1r + 2*Kdrag*w2i*vgsol*k*rhogeq*w1r**2*w2r**2*w3i**2 +&
        w3r**3*rhogeq*w2i*rhogsol*w1i*w2r**2*w1r**2 - 2*w3r**3*w2r*rhogeq*w1i**2*rhogsol*k**2*cs**2*w1r +&
        2*Kdrag*w2i*vgsol*k*rhogeq*w1r**2*w2r**2*w3r**2 - Kdrag*w2r**3*k*rhogeq*vgsol*w3i*w1r*w3r**2 -&
        2*rhogeq**2*w2r**2*w3i**4*k*vgsol*w2i**2 - Kdrag*w2i**2*rhogsol*w1i*w2r*w3i**4 -&
        rhogeq*w1i**2*w2i**2*rhogsol*w2r*w3i**4 - Kdrag*w2r**3*rhogsol*w3i*w1r**2*w3r**2 +&
        Kdrag*w2r**2*k*rhogeq*vgsol*w2i*w3i**3*w1i - 2*rhogeq*w1i**2*w2r**3*rhogsol*w3i**2*w3r**2 -&
        w2r**2*k*Kdrag**2*vdsol*w3i*w1i*w3r**2 + Kdrag*w2r**2*rhogsol*w2i*w1r*w3i**4 +&
        w2r**2*k*Kdrag**2*vgsol*w3i*w1i*w3r**2 - w2r**3*k*Kdrag**2*vgsol*w1r*w3i**2 +&
        4*w3i**3*cs**2*k**2*rhogeq*rhogsol*w2i*w2r*w1i**2 + Kdrag*k*rhogeq*vgsol*w2r**2*w1i*w3i**4 +&
        2*w3r**3*cs**2*k**3*rhogeq**2*w2r*vgsol*w2i*w1i + w3r**3*rhogeq*w2i*rhogsol*w1i**3*w2r**2 +&
        2*w3r**3*Kdrag*rhogsol*w2r**3*w1i*w1r - Kdrag**2*k*vdsol*w2i*w1i*w2r**2*w3r**2 +&
        w3r**3*rhogeq**2*k*w1r**3*vgsol*w2r**2 + w3r*rhogeq*cs**2*k**2*w1r**4*rhogsol*w2i**2 +&
        Kdrag**2*k*vgsol*w2i*w1i*w2r**2*w3i**2 - 2*cs**2*k**3*rhogeq**2*vgsol*w2r*w1r*w2i**2*w3r**2 +&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w1r**3 + w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i*w1i*w1r**2 -&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w1r**3 + Kdrag**2*k*vgsol*w2i*w1i*w2r**2*w3r**2 +&
        w3r**3*rhogeq**2*k*w1r*vgsol*w2r**2*w1i**2 + rhogeq*w2r**3*w3i**3*w1r*k*Kdrag*vdsol +&
        2*rhogsol*k**4*cs**4*rhogeq*w2r*w2i*w3i**3 + 2*rhogsol*k**4*cs**4*rhogeq*w2r*w2i*w3i*w3r**2 -&
        w2i**2*k*Kdrag**2*vgsol*w3r**3*w1r - Kdrag*w2r**2*rhogsol*k**2*cs**2*w3i*w1r*w3r**2 +&
        rhogeq**2*w1i*w2r**2*k*vgsol*w2i*w3i**4 - Kdrag*w2r**2*rhogsol*k**2*cs**2*w3i**3*w1r -&
        Kdrag*w2r**3*k*rhogeq*vgsol*w3i**3*w1r - w3i**4*vgsol*k*rhogeq*Kdrag*w2i**3 +&
        w3i*vdsol*k*Kdrag**2*w2i**3*w3r**2 - w3i*vgsol*k*Kdrag**2*w2i**3*w3r**2 +&
        2*w2i**3*w1i*w3r*rhogeq**2*cs**2*k**3*w1r*vgsol + w3r**4*rhogeq**2*vgsol*k*w2i**3*w1i -&
        w3r*Kdrag*rhogsol*w2i**4*w1i*w3i**2 - Kdrag*rhogsol*k**2*cs**2*w1i*w2r*w2i**2*w3i**2 -&
        w2i**2*k*Kdrag**2*vgsol*w3r*w1r*w3i**2 - Kdrag*rhogsol*k**2*cs**2*w3i**3*w1r*w2i**2 -&
        Kdrag*rhogsol*k**2*cs**2*w1i*w2r*w2i**2*w3r**2 + rhogeq*w2r**3*w3i*w1r*k*Kdrag*vdsol*w3r**2 -&
        w3r**4*vgsol*k*rhogeq*Kdrag*w2i**3 + w2i**2*k*Kdrag**2*vgsol*w3i*w1i*w3r**2 +&
        2*w3r*rhogeq*cs**2*k**2*w1i**2*rhogsol*w2i**2*w1r**2 - w2i**2*k*Kdrag**2*vdsol*w3i*w1i*w3r**2 +&
        w1i**3*rhogeq**2*cs**2*k**3*w3i*vgsol*w2r**2 + w2i**2*k*Kdrag**2*vdsol*w3r*w1r*w3i**2 -&
        2*w3i**2*rhogeq**2*cs**2*k**3*vgsol*w2i**2*w3r**2 + 4*w3r**3*w2i*cs**2*k**2*rhogeq*rhogsol*w2r*w1i*w1r +&
        Kdrag*w2i**2*vgsol*k*rhogeq*w1i*w3r**4 + w3r**3*rhogeq*rhogsol*w2r*w2i**2*w1i**2*w1r -&
        2*Kdrag*w2r**2*k*rhogeq*vgsol*w2i**2*w3i**3 + w2r*rhogeq**2*vgsol*k*w2i**2*w1r*w3i**4 -&
        Kdrag**2*k*vdsol*w2i*w1i*w2r**2*w3i**2 - 2*w2i**3*w1i**3*w3r*rhogeq*cs**2*k**2*rhogsol +&
        w3r**3*vdsol*Kdrag*k*rhogeq*w2i**3*w1r + 2*cs**2*k**3*rhogeq**2*vgsol*w2i**3*w1i*w3i**2 +&
        2*w3i**2*Kdrag*rhogsol*w2i**3*w1r*w3r**2 - Kdrag*rhogsol*k**2*cs**2*w3i*w1r*w2i**2*w3r**2 +&
        2*Kdrag*w2i**2*vgsol*k*rhogeq*w1i*w3r**2*w3i**2 - w3r*rhogeq*cs**4*k**4*rhogsol*w2i*w1i*w1r**2 -&
        w3i*cs**4*k**4*rhogeq*rhogsol*w2i*w1r**3 + 2*w3r**2*rhogeq**2*vgsol*k*w2i**3*w1i*w3i**2 -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i*w1i*w1r**2 + w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1r**4 -&
        2*Kdrag*w2r**2*k*rhogeq*vgsol*w2i**2*w3i*w3r**2 - cs**2*k**3*rhogeq*vgsol*Kdrag*w2i**3*w3i**2 +&
        2*w3r**2*rhogeq*rhogsol*w2i**4*w1r*w3i**2 + 2*Kdrag*k*rhogeq*vgsol*w2r**2*w1i*w2i**2*w3i**2 +&
        2*Kdrag*k*rhogeq*vgsol*w2r**2*w1i*w2i**2*w3r**2 + w3r**3*rhogsol*k**2*cs**2*Kdrag*w2i**3 +&
        w3i*Kdrag*rhogsol*w2i**4*w1r*w3r**2 + cs**2*k**3*rhogeq*vgsol*Kdrag*w2i**3*w3r**2 +&
        w1r**2*w3r**3*rhogeq*rhogsol*w2i**3*w1i + 2*w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i**2*w1r**2 -&
        2*w3r**2*vgsol*k*rhogeq*Kdrag*w2i**3*w3i**2 + w3r**4*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i +&
        w1r**2*vgsol*k*rhogeq**2*w2i**3*w1i*w3r**2 + w1r**2*vgsol*k*rhogeq**2*w2i**3*w1i*w3i**2 +&
        w1r**2*w3r*rhogeq*rhogsol*w2i**3*w1i*w3i**2 - 2*w2i*Kdrag*rhogsol*w2r**2*w1r*w3i**3*w1i +&
        rhogeq*rhogsol*w2r*w2i**2*w1r**2*w3i**3*w1i + w3r*rhogeq**2*k*vgsol*w2i**4*w1r*w3i**2 -&
        2*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2i**2*w1r + 2*w3r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2i**2*w1r +&
        w3r**3*rhogeq**2*k*vgsol*w2i**4*w1r + cs**2*k**3*rhogeq*vdsol*Kdrag*w2i**3*w3i**2 -&
        cs**2*k**3*rhogeq*vdsol*Kdrag*w2i**3*w3r**2 - w3i**3*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i**2 -&
        w3r*vgsol*k*rhogeq*Kdrag*w2i**3*w1r*w3i**2 + w3i*vgsol*k*rhogeq*Kdrag*w2i**3*w1i*w3r**2 -&
        w3i**3*cs**2*k**3*rhogeq**2*vgsol*w2i**3 + w3i**3*cs**2*k**3*rhogeq**2*vgsol*w2i*w1r**2 -&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**3*w3r**2 + w3r*vdsol*Kdrag*k*rhogeq*w2i**3*w1r*w3i**2 -&
        2*w2i**3*w1i*w3i**3*Kdrag*rhogsol*w1r + w3r*rhogeq*cs**2*k**2*w1i**4*rhogsol*w2i**2 -&
        2*w2i**3*w1i*w3r*rhogeq*cs**2*k**2*w1r**2*rhogsol - 2*w2i**2*w1i**2*w3i*cs**2*k**3*rhogeq*Kdrag*vgsol +&
        2*w2i**2*w1i**2*w3i*cs**2*k**3*rhogeq*Kdrag*vdsol + 2*w3i**3*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w1i +&
        2*w2i**2*w1i**2*w3r*rhogeq*cs**4*k**4*rhogsol - w3i**3*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i*w1i +&
        w3r*rhogsol*k**2*cs**2*Kdrag*w2i**3*w3i**2 - w3r**3*vgsol*k*rhogeq*Kdrag*w2i**3*w1r +&
        2*w3r**2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i*w3i**2 + w3i**3*vgsol*k*rhogeq*Kdrag*w2i**3*w1i +&
        2*w2i*w1i**2*w3r*w3i*cs**4*k**4*rhogeq*rhogsol + w3i**3*rhogeq**2*k*vgsol*w2i**4*w1i -&
        2*w2i**3*w1i*w3i*w3r**2*Kdrag*rhogsol*w1r + w3r**3*rhogeq*w1i**2*w2r**3*rhogsol*w1r +&
        w3r**3*rhogeq*rhogsol*w2r*w2i**2*w1r**3 + rhogeq*rhogsol*k**2*cs**2*w2i**4*w1r*w3i**2 -&
        rhogeq*rhogsol*k**2*cs**2*w2i**4*w1r*w3r**2 + w3i*vdsol*Kdrag*k*rhogeq*w2i**3*w1i*w3r**2 -&
        w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w1r*w1i**2 + w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i*w1i**3 +&
        w3i**3*vdsol*Kdrag*k*rhogeq*w2i**3*w1i + w3r*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w1r*w1i**2 -&
        2*w3r**2*rhogeq**2*vgsol*k*w2i**4*w3i**2 - w3r*rhogeq*cs**4*k**4*rhogsol*w2i*w1i**3 -&
        w3i*cs**4*k**4*rhogeq*rhogsol*w2i*w1r*w1i**2 + w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i**4 -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i*w1i**3 - w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1i**2*w3r**2 +&
        w3i*cs**2*k**3*rhogeq**2*vgsol*w2i*w1r**2*w3r**2 + 2*w3i*cs**2*k**3*rhogeq**2*vgsol*w2i**2*w1i*w3r**2 +&
        2*w3i*rhogsol*k**4*cs**4*rhogeq*w1i*w2i**2*w1r + 2*w3i*rhogeq**2*k*vgsol*w2i**3*w1r**2*w3r**2 -&
        2*w3i*rhogeq**2*k*vgsol*w2i**3*w1i**2*w3r**2 - 2*w3r*rhogeq*cs**2*k**2*w1r**2*rhogsol*w2i**2*w3i**2 +&
        w3i*rhogeq**2*k*vgsol*w2i**4*w1i*w3r**2 + w3i**4*rhogeq**2*vgsol*k*w2i**3*w1i +&
        2*w3i**3*rhogeq**2*k*vgsol*w2i**3*w1r**2 - 2*w3i**3*rhogeq**2*k*vgsol*w2i**3*w1i**2 -&
        w3i*cs**4*k**4*rhogeq*rhogsol*w2i*w1r*w3r**2 - w3r**4*rhogeq**2*cs**2*k**3*vgsol*w2i**2 -&
        w3r*Kdrag*rhogsol*k**2*cs**2*w1i*w2i**2*w3i**2 + w3r**3*rhogeq*cs**4*k**4*rhogsol*w2i**2 -&
        w3r*rhogeq*rhogsol*w2i**4*w1r**2*w3i**2 + rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2i**2*w3r**2 -&
        rhogeq*cs**4*k**4*w1r*rhogsol*w2i**2*w3i**2 - rhogeq*cs**4*k**4*w1r*rhogsol*w2i**2*w3r**2 -&
        w3r**3*rhogeq*cs**4*k**4*rhogsol*w2i*w1i + 3*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2i**2*w3i**2 +&
        2*w2i*w1i*rhogeq*cs**4*k**4*w1r*rhogsol*w3i**2 + 2*w2i*w1i*rhogeq*cs**4*k**4*w1r*rhogsol*w3r**2 -&
        2*w2i**3*w1i*rhogeq*cs**2*k**2*w1r*rhogsol*w3i**2 + 2*w2i**3*w1i*rhogeq*cs**2*k**2*w1r*rhogsol*w3r**2 -&
        rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w2i**2*w3r**2 + w3i*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i**2*w3r**2 -&
        w3i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i**2*w3r**2 - w3i**3*cs**4*k**4*rhogeq*rhogsol*w2i*w1r -&
        w3r**3*Kdrag*rhogsol*k**2*cs**2*w1i*w2i**2 - w3r*rhogeq*cs**4*k**4*rhogsol*w2i*w1i*w3i**2 +&
        w3r*rhogeq*cs**4*k**4*rhogsol*w2i**2*w3i**2 - 2*w3r**2*rhogeq**2*w1i**2*k*vgsol*w2r**2*w1r**2 +&
        w3r**2*w1r**3*rhogeq**2*w2r**3*k*vgsol - 2*w3r**3*rhogeq*cs**2*k**2*w1r**2*rhogsol*w2i**2 -&
        w3r**2*w3i*rhogeq*rhogsol*w2i*w1r*w2r**2*w1i**2 - w3r**2*w3i*rhogeq*k*Kdrag*vdsol*w2r**2*w1i**2 -&
        w3r*rhogeq*rhogsol*w2i**4*w1i**2*w3i**2 - 4*w3r**2*rhogeq**2*k*vgsol*w2r*w2i**2*w1r*w3i*w1i -&
        4*w3r**2*w2r*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i*w1r - w2i**2*k*Kdrag**2*vdsol*w1r**2*w3r**2 +&
        w2i**2*k*Kdrag**2*vgsol*w1r**2*w3i**2 + w2i**2*k*Kdrag**2*vgsol*w1r**2*w3r**2 -&
        4*vgsol*k*rhogeq**2*w2i**3*w1i*w3i**2*w3r*w1r - 4*vgsol*k*rhogeq**2*w2i**3*w1i*w3r**3*w1r -&
        2*rhogeq**2*cs**2*k**3*w1i**2*vgsol*w2i**2*w3r**2 + 2*w2i**2*w1i*Kdrag*rhogsol*k**2*cs**2*w1r*w3r**2 +&
        2*w2i**2*w1i*Kdrag*rhogsol*k**2*cs**2*w1r*w3i**2 + w2i**2*w1i**2*Kdrag**2*k*vdsol*w3r**2 -&
        3*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w2i**2*w3i**2 + 2*rhogeq**2*cs**2*k**3*w1r**2*vgsol*w2i**2*w3r**2 +&
        w3i**3*cs**2*k**3*rhogeq*Kdrag*vdsol*w2i**2 - w2i**2*w1i**2*Kdrag**2*k*vgsol*w3r**2 +&
        w2i**2*w1i**2*Kdrag**2*k*vdsol*w3i**2 - w2i**2*w1i**2*Kdrag**2*k*vgsol*w3i**2 -&
        w3i**3*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i**2 + w3i**4*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i +&
        w3i**3*cs**2*k**3*rhogeq*Kdrag*vgsol*w2i*w1i - w3r**3*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w1r +&
        w3r**3*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w1r - w2i**2*k*Kdrag**2*vdsol*w1r**2*w3i**2 +&
        2*Kdrag*w2i**2*k*rhogeq*vgsol*w3i*w1r**2*w3r**2 + 2*w1r**2*rhogeq*cs**2*k**3*w2i*Kdrag*vgsol*w3r**2 -&
        w1r**2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i*w3i**2 + w1r**2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i*w3r**2 +&
        2*Kdrag*w2i**2*k*rhogeq*vgsol*w3i**3*w1r**2 - 2*w1r**2*rhogeq**2*w1i**2*k*vgsol*w2i**2*w3i**2 -&
        2*w1r**2*rhogeq**2*w1i**2*k*vgsol*w2i**2*w3r**2 + Kdrag*w2i**2*k*rhogeq*vgsol*w3i**4*w1i -&
        w1i**2*vdsol*Kdrag*k*rhogeq*w2i**3*w3i**2 - w1i**2*vdsol*Kdrag*k*rhogeq*w2i**3*w3r**2 +&
        w1i**3*vdsol*Kdrag*k*rhogeq*w2i**2*w3r**2 + w1i**3*vdsol*Kdrag*k*rhogeq*w2i**2*w3i**2 -&
        w1i**2*w3i**3*rhogeq*rhogsol*w2i**3*w1r - w1i**2*w3i*rhogeq*rhogsol*w2i**3*w1r*w3r**2 -&
        Kdrag**2*k*vdsol*w2r**2*w1r**2*w3i**2 + 2*w1i**2*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w3i**2 +&
        2*w1i**2*rhogeq*cs**2*k**2*w1r*rhogsol*w2i**2*w3i**2 - 2*rhogeq**2*cs**2*k**3*vgsol*w2i*w1i*w3r**3*w1r -&
        2*rhogeq*cs**2*k**2*w1r*rhogsol*w2i**2*w3i**3*w1i + Kdrag**2*k*vgsol*w2r**2*w1r**2*w3i**2 -&
        2*rhogeq*cs**2*k**2*w1r*rhogsol*w2i**2*w3r**2*w3i*w1i + w1i**3*vgsol*k*rhogeq**2*w2i**3*w3r**2 +&
        w1i**3*vgsol*k*rhogeq**2*w2i**3*w3i**2 - 2*rhogeq**2*w1i**2*k*vgsol*w3i**2*w2r**2*w1r**2 +&
        2*Kdrag*rhogsol*k**2*cs**2*w1i*w2r**2*w3i**2*w1r + 2*w3r**2*Kdrag*rhogsol*k**2*cs**2*w1i*w2r**2*w1r -&
        2*w3r**2*rhogeq*cs**2*k**2*w2r**2*w1r**3*rhogsol + w1i**3*w3r*rhogeq*rhogsol*w2i**3*w3i**2 +&
        w3r**2*w1r**3*w2r*rhogeq**2*cs**2*k**3*vgsol - w1r**4*rhogeq**2*k*vgsol*w2i**2*w3i**2 -&
        w1r**4*rhogeq**2*k*vgsol*w2i**2*w3r**2 - w1i**4*rhogeq**2*k*vgsol*w2i**2*w3r**2 -&
        w1i**4*rhogeq**2*k*vgsol*w2i**2*w3i**2 + w1r**2*vdsol*Kdrag*k*rhogeq*w2i**2*w1i*w3r**2 -&
        w1r**3*w3i*rhogeq*rhogsol*w2i**3*w3r**2 + w1r**2*vdsol*Kdrag*k*rhogeq*w2i**2*w1i*w3i**2 -&
        2*w1r**2*rhogeq*cs**2*k**3*w2i*Kdrag*vdsol*w3r**2 - w3r**2*w2i*rhogeq*k*Kdrag*vdsol*w2r**2*w1r**2 -&
        w3r**2*rhogeq**2*w1r**4*k*vgsol*w2r**2 + w3r**2*w1r**4*w2r*cs**2*k**2*rhogeq*rhogsol -&
        w1i**3*Kdrag*k*rhogeq*vgsol*w2i**2*w3r**2 - w1i**3*Kdrag*k*rhogeq*vgsol*w2i**2*w3i**2 +&
        rhogeq**2*w3i**3*w1i*k*vgsol*w2r**2*w1r**2 + w3r**2*rhogeq*w1i*k*Kdrag*vdsol*w2r**2*w1r**2 +&
        w3r**2*rhogeq*rhogsol*w2r*w2i**2*w1r**2*w3i*w1i + w3r**2*w1i**2*rhogeq**2*k*vgsol*w2r*w2i**2*w1r -&
        4*w3r**2*rhogeq**2*w2r**3*k*w1r*vgsol*w3i*w1i - 2*w3r**2*w1r**2*w2i*cs**2*k**2*rhogeq*rhogsol*w2r*w1i +&
        w3r**2*rhogeq*w1i**3*k*Kdrag*vdsol*w2r**2 + w3r**2*rhogeq*rhogsol*w2r*w2i**2*w1i**3*w3i -&
        w3r**2*rhogeq**2*w1i**4*k*vgsol*w2r**2 - w3r**2*Kdrag**2*k*vgsol*w2r**2*w1i**2 +&
        w3r**2*Kdrag**2*k*vdsol*w2r**2*w1i**2 - 2*w3r**2*w1i**3*w2i*cs**2*k**2*rhogeq*rhogsol*w2r +&
        w3r**2*w1i**3*w2i*rhogeq**2*k*vgsol*w2r**2 + w3r**2*rhogeq*w2r**3*w1r**2*rhogsol*w3i*w1i +&
        2*w3r**2*rhogeq*rhogsol*k**2*cs**2*w2r**2*w1r*w3i*w1i - 2*w3r**2*w2r*cs**2*k**2*rhogeq*rhogsol*w1r**2*w3i*w1i +&
        w3r**2*w1i**2*rhogeq**2*w2r**3*k*w1r*vgsol - 2*w3r**2*w2r*rhogeq*w1i**3*rhogsol*k**2*cs**2*w3i -&
        2*w3r**2*w2i*Kdrag*rhogsol*w2r**2*w1r*w3i*w1i - 2*w3r**2*w2r*Kdrag*rhogeq*k*vgsol*w2i**2*w1i*w1r -&
        2*w3r**2*Kdrag*k*rhogeq*vgsol*w2r**3*w1i*w1r - w3r**2*w3i*rhogeq*rhogsol*w2i*w1r**3*w2r**2 -&
        w3r**2*w1i**3*rhogeq*w2r**2*k*Kdrag*vgsol + 2*w3r**2*w2r*rhogeq**2*cs**2*k**3*w1r*vgsol*w3i*w1i -&
        2*w3r**2*w2r*rhogeq*cs**2*k**3*w1i*Kdrag*vdsol*w1r + 2*w3r**2*w2r*rhogeq*cs**2*k**3*w1i*Kdrag*vgsol*w1r +&
        w3r**2*w1r**3*rhogeq**2*k*vgsol*w2r*w2i**2 - 2*w3r**2*w2r*rhogeq*cs**4*k**4*w1r**2*rhogsol +&
        w1i**3*rhogeq**2*cs**2*k**3*vgsol*w2i*w3r**2 - w1i**3*rhogeq**2*cs**2*k**3*vgsol*w2i*w3i**2 -&
        w3r**2*w1r**2*rhogeq*w2r**2*w1i*k*Kdrag*vgsol + rhogeq*w1i**3*w3i**2*k*Kdrag*vdsol*w2r**2 +&
        w3r**2*w1r**2*w2i*rhogeq**2*k*vgsol*w2r**2*w1i + w2i**2*w1r**3*w3r**3*rhogeq**2*k*vgsol +&
        w3r**2*Kdrag**2*k*vgsol*w2r**2*w1r**2 - w3r**2*Kdrag**2*k*vdsol*w2r**2*w1r**2 +&
        w3r**2*w1i**4*w2r*rhogeq*rhogsol*k**2*cs**2 - rhogeq*w3i**3*k*Kdrag*vdsol*w2r**2*w1r**2 -&
        w2i**2*rhogeq**2*cs**2*k**3*w3i*w1i**3*vgsol + w3r**2*rhogeq*w1i**3*w2r**3*rhogsol*w3i -&
        w2i**2*w3r*rhogeq**2*cs**2*k**3*w1i**2*vgsol*w1r + 2*w3r**2*w2i*rhogeq*cs**2*k**2*rhogsol*w2r**2*w1i*w1r +&
        w3r**2*w3i*rhogeq**2*k*vgsol*w1i*w2r**2*w1r**2 - w3r**2*w3i*rhogeq*k*Kdrag*vdsol*w2r**2*w1r**2 -&
        w2i**2*w3r*rhogeq**2*cs**2*k**3*w1r**3*vgsol - w3r**2*w2i*w1i**2*rhogeq*k*Kdrag*vdsol*w2r**2 -&
        w2i**2*rhogeq**2*cs**2*k**3*w3i*w1i*vgsol*w1r**2 - 2*w2i**2*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w1i*w1r**2 -&
        2*w2i**2*Kdrag*k*rhogeq*vgsol*w1i*w3r**3*w1r - 2*w2i**2*w3r*w3i*cs**2*k**2*rhogeq*rhogsol*w1i**3 -&
        w2i**2*w1r**2*w3i*vdsol*Kdrag*k*rhogeq*w3r**2 - 2*w3r**2*w1i**2*rhogeq*rhogsol*k**2*cs**2*w2r**2*w1r +&
        w3r**2*w3i*rhogeq**2*k*vgsol*w1i**3*w2r**2 - 2*w3r**2*w2r*rhogeq**2*cs**2*k**3*w1r*vgsol*w2i*w3i +&
        2*w3r**2*w1i**2*w2r*cs**2*k**2*rhogeq*rhogsol*w1r**2 + w3r**2*w1i**2*w2r*rhogeq**2*cs**2*k**3*w1r*vgsol +&
        w2i**2*w1i**3*rhogeq**2*w3i**3*k*vgsol + w2i**2*w1i**2*w3r*rhogeq**2*k*w1r*vgsol*w3i**2 +&
        w2i**2*w1i**3*rhogeq**2*w3i*k*vgsol*w3r**2 + w2i**2*w1i**2*w3r**3*rhogeq**2*k*w1r*vgsol -&
        w2i**2*w1r**2*w3i**3*vdsol*Kdrag*k*rhogeq - 2*w2i**2*Kdrag*k*rhogeq*vgsol*w1i*w3i**2*w3r*w1r -&
        w2i**2*w1i**2*w3i**3*vdsol*Kdrag*k*rhogeq - w2i**2*w1i**2*w3i*vdsol*Kdrag*k*rhogeq*w3r**2 +&
        w2i**2*w1r**3*w3r*rhogeq**2*k*vgsol*w3i**2 + w2i**2*w1r**2*rhogeq**2*w3i**3*w1i*k*vgsol +&
        w2i**2*w1r**2*rhogeq**2*w3i*w1i*k*vgsol*w3r**2 - 4*w3r*w2i*cs**2*k**3*rhogeq*Kdrag*vgsol*w2r*w3i*w1i +&
        w3r*w1r**4*rhogeq*cs**2*k**2*w2r**2*rhogsol + rhogeq*w1i**3*w2r**3*rhogsol*w3i**3)*rhodeq/(w1i**2 - 2*w3i*w1i +&
        w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)/(w3r**2 + w3i**2)/Kdrag/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i - 2*w2r*w1r&
        + w1i**2)/rhogeq/(w2i**2 + w2r**2)

    rhod1i = - ( - 2*rhogsol*Kdrag*w1i**3*rhodeq*w3r**2*w2i**3 - 2*rhodsol*rhogeq*Kdrag*w3i**4*w2i**2*w2r**2 -&
        rhodsol*rhogeq*Kdrag*w2r**4*w3i**4 - rhodsol*rhogeq*Kdrag*w3i**4*w2i**4 -&
        2*rhogsol*Kdrag*w1i**3*rhodeq*w3r**2*w2i*w2r**2 - 2*rhodsol*rhogeq*Kdrag*w2r**4*w3i**2*w3r**2 -&
        2*rhogsol*Kdrag*w1i**3*rhodeq*w3r**2*w2r**2*w3i - 2*rhogsol*Kdrag*w1i**3*rhodeq*w3r**2*w2i**2*w3i -&
        2*rhogsol*Kdrag*w1i**3*rhodeq*w2r**2*w3i**3 - 2*rhogsol*Kdrag*w1i**3*rhodeq*w2i**2*w3i**3 -&
        2*rhogsol*Kdrag*w1i**3*rhodeq*w2i**3*w3i**2 - rhodsol*rhogeq*Kdrag*w3r**4*w2i**4 -&
        rhodsol*rhogeq*Kdrag*w3r**4*w2r**4 - 2*rhodsol*rhogeq*Kdrag*w3r**4*w2i**2*w2r**2 -&
        4*rhodsol*rhogeq*Kdrag*w3r**2*w2i**2*w2r**2*w3i**2 - 2*rhodsol*rhogeq*Kdrag*w3i**2*w2i**4*w3r**2 -&
        2*rhogsol*Kdrag*w1i**3*rhodeq*w2i*w2r**2*w3i**2 + rhogsol*rhogeq*k**4*cs**4*w1i**3*rhodeq*w3r*w2r -&
        rhogsol*rhogeq*k**4*cs**4*w1i**3*rhodeq*w2i*w3i + w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r**2*w2i**2 -&
        w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r*w3r*w3i**2 -&
        w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r*w2r*w2i**2 - w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r*w2r**3 -&
        w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i**2*w2i**4 - w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i**4*w2i**2 +&
        w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i**4*w2r**2 - vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w3r*w2i**3 +&
        w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**2*w2i**4 - w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w2r**4*w3i**2 -&
        2*w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i**2*w2i**2*w2r**2 -&
        vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w3r**3*w2i - w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**4*w2i**2 +&
        w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**2*w2r**4 - 2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w2r**2*w3i**3 +&
        cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3i**4*w2i +&
        2*w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**2*w2i**2*w2r**2 -&
        2*w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**2*w2i**2*w3i**2 +&
        w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**4*w2r**2 + rhogeq*k*Kdrag*vdsol*w1i**2*rhodeq*w3r*w2i**2*w3i**2 +&
        rhogeq*k*Kdrag*vdsol*w1i**2*rhodeq*w3r*w2r**2*w3i**2 + rhogeq*k*Kdrag*vdsol*w1i**2*rhodeq*w2r*w2i**2*w3i**2 +&
        rhogeq*k*Kdrag*vdsol*w1i**2*rhodeq*w2r**3*w3r**2 + rhogeq*k*Kdrag*vdsol*w1i**2*rhodeq*w2r*w2i**2*w3r**2 +&
        rhogeq*k*Kdrag*vdsol*w1i**2*rhodeq*w2r**3*w3i**2 + rhogeq*k*Kdrag*vdsol*w1i**2*rhodeq*w3r**3*w2i**2 +&
        rhogeq*k*Kdrag*vdsol*w1i**2*rhodeq*w3r**3*w2r**2 + 2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w2i**2*w3i**3 +&
        2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w2i*w2r**2*w3i**2 +&
        cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3i*w2i**4 + 2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w2i**3*w3i**2 +&
        2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3i*w2r**2*w2i**2 +&
        cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3i*w2r**4 +&
        2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3i**2*w2i*w3r**2 +&
        2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3r**2*w2i**2*w3i -&
        2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3r**2*w2r**2*w3i -&
        2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3r**2*w2i*w2r**2 -&
        2*cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3r**2*w2i**3 -&
        2*w1i*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w2i**3*w3i -&
        2*w1i*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i**2*w2i**2 -&
        2*w1i*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i*w2i*w3r**2 -&
        2*w1i*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w2i*w3i*w2r**2 +&
        2*w1i*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**2*w2r**2 -&
        2*w1i*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i**3*w2i +&
        w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r**2*w3i**2 + cs**2*k**2*rhogsol*w1i**2*rhogeq*rhodeq*w3r**4*w2i +&
        w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i**3*w3i + w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3i**2*w2i**2 +&
        w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3i**3*w2i + w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i*w3i*w2r**2 -&
        3*w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r**2*w2r**2 +&
        4*w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r*w2r*w2i*w3i -&
        w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r*w3r**3 + w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3i*w2i*w3r**2 +&
        2*w1r*w1i*rhogsol*Kdrag*rhodeq*w3r*w2r**2*w2i*w3i**2 + w1r*vgsol*k*w1i**2*rhogeq**2*rhodeq*w3r**2*w2r**2*w3i +&
        w1r*vgsol*k*w1i**2*rhogeq**2*rhodeq*w3r**2*w2i**2*w3i + w1r**2*rhogsol*rhogeq*rhodeq*w2r**4*w3i**3 -&
        2*rhogsol*rhogeq*k**2*cs**2*w1r*rhodeq*w2r*w3i**4*w2i + 2*w1r*w1i*rhogsol*Kdrag*rhodeq*w3r**3*w2i*w2r**2 +&
        2*w1r*w1i*rhogsol*Kdrag*rhodeq*w3i*w2i**2*w3r**2*w2r + 2*w1r*w1i*rhogsol*Kdrag*rhodeq*w3i*w2r**3*w3r**2 +&
        2*w1r*w1i*rhogsol*Kdrag*rhodeq*w3r*w2i**3*w3i**2 - 4*vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w2r*w2i*w3i**2 -&
        vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w2r*w3i**3 - 2*rhogsol*rhogeq*k**2*cs**2*w1r*rhodeq*w2r*w2i*w3r**4 -&
        4*rhogsol*rhogeq*k**2*cs**2*w1r*rhodeq*w2r*w2i*w3i**2*w3r**2 -&
        2*rhogsol*rhogeq*k**2*cs**2*w1r*rhodeq*w3i*w3r*w2i**4 - 2*rhogsol*rhogeq*k**2*cs**2*w1r*rhodeq*w2r**4*w3i*w3r -&
        4*rhogsol*rhogeq*k**2*cs**2*w1r*rhodeq*w3r*w3i*w2r**2*w2i**2 + 2*w1r*w1i*rhogsol*Kdrag*rhodeq*w2i**3*w3r**3 -&
        vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w2r*w3i*w3r**2 -&
        vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w3r*w2i*w3i**2 -&
        vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w3r*w2i*w2r**2 -&
        4*vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w3i*w3r*w2i**2 - vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w2r**3*w3i&
        - vgsol*Kdrag*w1i*k**3*cs**2*rhogeq*rhodeq*w2r*w3i*w2i**2 +&
        w1r*vgsol*k*w1i**2*rhogeq**2*rhodeq*w2i*w2r**2*w3i**2 + 4*w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2r*w2i*w3r**2 +&
        w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2r*w3i*w3r**2 + w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w3r*w2i*w3i**2 +&
        w1r**2*rhogsol*rhogeq*rhodeq*w2r**2*w2i*w3r**4 + w1r*vgsol*k*w1i**2*rhogeq**2*rhodeq*w2r**2*w3i**3 +&
        w1r*vgsol*k*w1i**2*rhogeq**2*rhodeq*w2i**2*w3i**3 + w1r*vgsol*k*w1i**2*rhogeq**2*rhodeq*w3r**2*w2i*w2r**2 +&
        w1r*vgsol*k*w1i**2*rhogeq**2*rhodeq*w3r**2*w2i**3 + 2*w1r*w1i*rhogsol*Kdrag*rhodeq*w2r*w2i**2*w3i**3 +&
        2*w1r*w1i*rhogsol*Kdrag*rhodeq*w3i**3*w2r**3 - 2*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w2r**3*w3r**2 -&
        2*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w2r*w2i**2*w3r**2 + w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2r**3*w3i&
        + w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2r*w3i*w2i**2 + w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2r*w3i**3 +&
        vgsol*Kdrag*k*rhogeq*w1r**3*rhodeq*w3r**2*w2i**2 + w1r*vgsol*k*w1i**2*rhogeq**2*rhodeq*w2i**3*w3i**2 +&
        w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w3r**3*w2i + w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w3r*w2i**3 +&
        w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w3r*w2i*w2r**2 + vgsol*Kdrag*k*rhogeq*w1r**3*rhodeq*w3r**2*w2r**2 +&
        vgsol*Kdrag*k*rhogeq*w1r**3*rhodeq*w3i**2*w2i**2 + vgsol*Kdrag*k*rhogeq*w1r**3*rhodeq*w2r**2*w3i**2 +&
        2*w1i*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**2*w2r**2*w3i**2 -&
        2*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w3r**3*w2r**2 +&
        2*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w3r**3*w2i**2 +&
        4*w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w3i*w3r*w2r**2 -&
        2*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w3r*w2r**2*w3i**2 +&
        2*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w3r*w2i**2*w3i**2 +&
        4*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w2r*w2i*w3r**2*w3i +&
        4*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w2i*w3i*w3r*w2r**2 +&
        4*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w2i**3*w3i*w3r - 2*w1r*vdsol*Kdrag**2*k*w1i*rhodeq*w3r**2*w2r**2 -&
        2*w1r*vdsol*Kdrag**2*k*w1i*rhodeq*w3i**2*w2i**2 - 2*w1r*vdsol*Kdrag**2*k*w1i*rhodeq*w2r**2*w3i**2 -&
        w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r**4*w2i + 4*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w2r*w2i*w3i**3 +&
        2*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w2r*w2i**2*w3i**2 - 2*w1r*vdsol*Kdrag**2*k*w1i*rhodeq*w3r**2*w2i**2&
        - 2*w1r**3*rhogsol*Kdrag*rhodeq*w3r**3*w2r**2 - 2*w1r**3*rhogsol*Kdrag*rhodeq*w2r**3*w3i**2 +&
        w1r*cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r**3*w2i**2 - w1i**4*Kdrag*rhogeq*rhodsol*w2r**2*w3i**2 -&
        2*w1r**3*rhogsol*Kdrag*rhodeq*w3r**3*w2i**2 - 2*w1r**3*rhogsol*Kdrag*rhodeq*w2r*w2i**2*w3r**2 +&
        w1r**2*rhogsol*Kdrag*rhodeq*w2i*w2r**2*w3i**3 + w1r**2*rhogsol*Kdrag*rhodeq*w2i**3*w3i**3 +&
        w1r**2*rhogsol*Kdrag*rhodeq*w3i**4*w2r**2 - w1i**4*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2 -&
        w1i**4*Kdrag*rhogeq*rhodsol*w3r**2*w2r**2 - w1i**4*Kdrag*rhogeq*rhodsol*w3i**2*w2i**2 +&
        2*w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**2*w3i**2 + 2*w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2r**2*w3i**2 +&
        3*w1r**2*rhogsol*Kdrag*rhodeq*w3r*w2r**3*w3i**2 + w1r**2*rhogsol*Kdrag*rhodeq*w2r**4*w3i**2 +&
        w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2r**4 + w1r**2*rhogsol*Kdrag*rhodeq*w2i**3*w3i*w3r**2 +&
        w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**4 + w1r**2*rhogsol*Kdrag*rhodeq*w2i*w3i*w2r**2*w3r**2 +&
        3*w1r**2*rhogsol*Kdrag*rhodeq*w2r*w3r**3*w2i**2 + 2*w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**2*w2r**2 +&
        3*w1r**2*rhogsol*Kdrag*rhodeq*w2r*w3r*w2i**2*w3i**2 + 2*w1r**2*rhogsol*Kdrag*rhodeq*w3i**2*w2i**2*w2r**2 +&
        w1r**2*rhogsol*Kdrag*rhodeq*w3i**2*w2i**4 + w1r**2*rhogsol*Kdrag*rhodeq*w3i**4*w2i**2 +&
        w1r**2*rhogeq*k*Kdrag*vdsol*rhodeq*w3r*w2i**2*w3i**2 + w1r**2*rhogeq*k*Kdrag*vdsol*rhodeq*w3r*w2r**2*w3i**2 +&
        w1r**2*rhogeq*k*Kdrag*vdsol*rhodeq*w2r**3*w3i**2 + w1r**2*rhogeq*k*Kdrag*vdsol*rhodeq*w2r*w2i**2*w3r**2 +&
        w1r**2*rhogeq*k*Kdrag*vdsol*rhodeq*w2r*w2i**2*w3i**2 + w1r**2*rhogsol*Kdrag*rhodeq*w3r**4*w2i**2 +&
        w1r**2*rhogsol*Kdrag*rhodeq*w3r**4*w2r**2 + 3*w1r**2*rhogsol*Kdrag*rhodeq*w3r**3*w2r**3 -&
        2*rhogsol*rhogeq*k**2*cs**2*w1i**3*rhodeq*w3i*w2i*w3r**2 - 2*rhogsol*rhogeq*k**2*cs**2*w1i**3*rhodeq*w3i**3*w2i&
        - 2*rhogsol*rhogeq*k**2*cs**2*w1i**3*rhodeq*w2i*w3i*w2r**2 -&
        2*rhogsol*rhogeq*k**2*cs**2*w1i**3*rhodeq*w2i**3*w3i + 2*rhogsol*rhogeq*k**2*cs**2*w1i**3*rhodeq*w3r**2*w2r**2&
        - 2*rhogsol*rhogeq*k**2*cs**2*w1i**3*rhodeq*w3i**2*w2i**2 + w1r**2*rhogeq*k*Kdrag*vdsol*rhodeq*w3r**3*w2i**2 +&
        w1r**2*rhogeq*k*Kdrag*vdsol*rhodeq*w3r**3*w2r**2 + w1r**2*rhogeq*k*Kdrag*vdsol*rhodeq*w2r**3*w3r**2 +&
        rhogsol*rhogeq*w1i**3*rhodeq*w2r*w3r*w2i**2*w3i**2 - rhogsol*rhogeq*w1i**3*rhodeq*w2i*w2r**2*w3i**3 -&
        rhogsol*rhogeq*w1i**3*rhodeq*w2i**3*w3i**3 + w1i**4*rhogsol*Kdrag*rhodeq*w3r**2*w2i**2 +&
        rhogsol*rhogeq*w1i**3*rhodeq*w2r*w3r**3*w2i**2 - rhogsol*rhogeq*w1i**3*rhodeq*w2i**3*w3i*w3r**2 -&
        rhogsol*rhogeq*w1i**3*rhodeq*w2i*w3i*w2r**2*w3r**2 + rhogsol*rhogeq*w1i**3*rhodeq*w3r**3*w2r**3 +&
        rhogsol*rhogeq*w1i**3*rhodeq*w3r*w2r**3*w3i**2 + w1i**4*rhogsol*Kdrag*rhodeq*w3r**2*w2r**2 +&
        w1i**4*rhogsol*Kdrag*rhodeq*w3i**2*w2i**2 + w1i**4*rhogsol*Kdrag*rhodeq*w2r**2*w3i**2 +&
        w1r*Kdrag*rhogeq*k*vgsol*w1i**2*rhodeq*w3r**2*w2i**2 + w1r*Kdrag*rhogeq*k*vgsol*w1i**2*rhodeq*w3r**2*w2r**2 +&
        w1r*Kdrag*rhogeq*k*vgsol*w1i**2*rhodeq*w3i**2*w2i**2 + w1r*Kdrag*rhogeq*k*vgsol*w1i**2*rhodeq*w2r**2*w3i**2 +&
        w1i*vdsol*k*Kdrag**2*rhodeq*w3r**3*w2r**2 + w1i*vdsol*k*Kdrag**2*rhodeq*w2r**3*w3r**2 +&
        w1i*vdsol*k*Kdrag**2*rhodeq*w2r*w2i**2*w3r**2 + w1i*vdsol*k*Kdrag**2*rhodeq*w3r*w2i**2*w3i**2 +&
        w1i*vdsol*k*Kdrag**2*rhodeq*w3r*w2r**2*w3i**2 + w1i*vdsol*k*Kdrag**2*rhodeq*w2r**3*w3i**2 +&
        w1i*vdsol*k*Kdrag**2*rhodeq*w2r*w2i**2*w3i**2 + Kdrag*cs**2*k**2*rhogsol*w1i**2*rhodeq*w2r**2*w3i**2 -&
        rhogsol*w1r**3*rhogeq*rhodeq*w2r*w2i**2*w3i**3 - rhogsol*w1r**3*rhogeq*rhodeq*w3i**3*w2r**3 +&
        w1i*vdsol*k*Kdrag**2*rhodeq*w3r**3*w2i**2 - rhogsol*w1r**3*rhogeq*rhodeq*w3r*w2r**2*w2i*w3i**2 +&
        Kdrag*cs**2*k**2*rhogsol*w1i**2*rhodeq*w3r**2*w2i**2 + Kdrag*cs**2*k**2*rhogsol*w1i**2*rhodeq*w3r**2*w2r**2 +&
        Kdrag*cs**2*k**2*rhogsol*w1i**2*rhodeq*w3i**2*w2i**2 - w1i**2*Kdrag*rhogeq*rhodsol*w3i**4*w2i**2 -&
        w1i**2*Kdrag*rhogeq*rhodsol*w3i**4*w2r**2 - 4*w1i**2*Kdrag*rhogeq*rhodsol*w2i**3*w3i**3 -&
        rhogsol*w1r**3*rhogeq*rhodeq*w3r**3*w2i*w2r**2 - 2*w1i**2*Kdrag*rhogeq*rhodsol*w3r**2*w2r**2*w3i**2 -&
        w1i**2*Kdrag*rhogeq*rhodsol*w2r**4*w3i**2 - 2*w1i**2*Kdrag*rhogeq*rhodsol*w3i**2*w2i**2*w2r**2 -&
        4*w1i**2*Kdrag*rhogeq*rhodsol*w2i*w2r**2*w3i**3 - 4*w1i**2*Kdrag*rhogeq*rhodsol*w2i*w3i*w2r**2*w3r**2 -&
        w1i**2*Kdrag*rhogeq*rhodsol*w3i**2*w2i**4 - rhogsol*w1r**3*rhogeq*rhodeq*w2i**3*w3r**3 -&
        rhogsol*w1r**3*rhogeq*rhodeq*w3i*w2i**2*w3r**2*w2r - rhogsol*w1r**3*rhogeq*rhodeq*w3i*w2r**3*w3r**2 -&
        rhogsol*w1r**3*rhogeq*rhodeq*w3r*w2i**3*w3i**2 - 2*w1i**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2*w2r**2 -&
        4*w1i**2*Kdrag*rhogeq*rhodsol*w2i**3*w3i*w3r**2 - w1i**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**4 -&
        2*w1i**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2*w3i**2 - w1i**2*Kdrag*rhogeq*rhodsol*w3r**4*w2r**2 -&
        w1i**2*Kdrag*rhogeq*rhodsol*w3r**2*w2r**4 + 2*w1i*Kdrag*rhogeq*rhodsol*w2r**4*w3r**2*w3i +&
        4*w1i*Kdrag*rhogeq*rhodsol*w2r**2*w3i*w3r**2*w2i**2 + 2*w1i*Kdrag*rhogeq*rhodsol*w2r**4*w3i**3 +&
        2*w1i*Kdrag*rhogeq*rhodsol*w2i**4*w3i**3 + 2*w1i*Kdrag*rhogeq*rhodsol*w2r**2*w2i*w3r**4 +&
        4*w1i*Kdrag*rhogeq*rhodsol*w3i**2*w2r**2*w2i*w3r**2 + 4*w1i*Kdrag*rhogeq*rhodsol*w2i**3*w3i**2*w3r**2 -&
        2*w1i*rhogsol*rhogeq*rhodeq*w3i**4*w2i**2*w2r**2 - w1i*rhogsol*rhogeq*rhodeq*w2r**4*w3i**4 -&
        w1i*rhogsol*rhogeq*rhodeq*w3i**4*w2i**4 + 2*w1i*Kdrag*rhogeq*rhodsol*w2i**3*w3r**4 -&
        2*w1i*rhogsol*rhogeq*rhodeq*w2r**4*w3i**2*w3r**2 + 2*w1i*Kdrag*rhogeq*rhodsol*w2i**4*w3r**2*w3i +&
        2*w1i*Kdrag*rhogeq*rhodsol*w2i*w2r**2*w3i**4 + 2*w1i*Kdrag*rhogeq*rhodsol*w3i**4*w2i**3 +&
        4*w1i*Kdrag*rhogeq*rhodsol*w2r**2*w2i**2*w3i**3 - w1i**2*Kdrag*rhogeq*rhodsol*w3r**4*w2i**2 +&
        2*w1r**2*rhogsol*rhogeq*rhodeq*w2r**2*w2i**2*w3i**3 - w1i*rhogsol*rhogeq*rhodeq*w3r**4*w2i**4 +&
        w1r**2*rhogsol*rhogeq*rhodeq*w2i**4*w3r**2*w3i + w1r**2*rhogsol*rhogeq*rhodeq*w2r**4*w3r**2*w3i +&
        w1r**2*rhogsol*rhogeq*rhodeq*w2i*w2r**2*w3i**4 + w1r**2*rhogsol*rhogeq*rhodeq*w2i**4*w3i**3 +&
        2*w1r**2*rhogsol*rhogeq*rhodeq*w3i**2*w2r**2*w2i*w3r**2 +&
        2*w1r**2*rhogsol*rhogeq*rhodeq*w2r**2*w3i*w3r**2*w2i**2 -&
        2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i*w2r**2*w2i**2 + 2*w1r**2*rhogsol*rhogeq*rhodeq*w2i**3*w3i**2*w3r**2&
        + w1r**2*rhogsol*rhogeq*rhodeq*w3i**4*w2i**3 - w1i*rhogsol*rhogeq*rhodeq*w3r**4*w2r**4 -&
        2*w1i*rhogsol*rhogeq*rhodeq*w3r**4*w2i**2*w2r**2 - 4*w1i*rhogsol*rhogeq*rhodeq*w3r**2*w2i**2*w2r**2*w3i**2 -&
        2*w1i*rhogsol*rhogeq*rhodeq*w3i**2*w2i**4*w3r**2 - w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i*w2i**4 -&
        2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2i**3*w3i**2 - 2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2i**2*w3i**3 -&
        2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3r*w2i*w3i**2 - w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i**4*w2i -&
        w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i*w2r**4 - 2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i*w3r*w2i**2*w2r -&
        2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i*w3r*w2r**3 -&
        2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2i*w2r**2*w3i**2 -&
        2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i**2*w2i*w3r**2 -&
        2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r**2*w2i**2*w3i + w1r**2*rhogsol*rhogeq*rhodeq*w2i**3*w3r**4 -&
        2*w1r*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3r**3*w2i - vgsol*Kdrag*rhogeq*k**3*cs**2*w1i**3*rhodeq*w2r*w3i -&
        w1r*rhogsol*rhogeq*w1i**2*rhodeq*w3i**3*w2r**3 + rhogsol*w1r*k**4*cs**4*rhogeq*w1i**2*rhodeq*w2i*w3r -&
        w1r*rhogsol*rhogeq*w1i**2*rhodeq*w3i*w2r**3*w3r**2 - w1r*rhogsol*rhogeq*w1i**2*rhodeq*w3r*w2i**3*w3i**2 -&
        w1r*rhogsol*rhogeq*w1i**2*rhodeq*w3r*w2r**2*w2i*w3i**2 - w1r*rhogsol*rhogeq*w1i**2*rhodeq*w3r**3*w2i*w2r**2 -&
        w1r*rhogsol*rhogeq*w1i**2*rhodeq*w2i**3*w3r**3 + rhogsol*rhogeq*k**2*cs**2*w1i**4*rhodeq*w2r**2*w3i -&
        w1r*rhogsol*rhogeq*w1i**2*rhodeq*w3i*w2i**2*w3r**2*w2r - w1r*rhogsol*rhogeq*w1i**2*rhodeq*w2r*w2i**2*w3i**3 +&
        rhogsol*w1r*k**4*cs**4*rhogeq*w1i**2*rhodeq*w2r*w3i - vgsol*Kdrag*rhogeq*k**3*cs**2*w1i**3*rhodeq*w2i*w3r -&
        w1i*rhogsol*Kdrag*rhodeq*w2i*w2r**2*w3i**4 - w1i*rhogsol*Kdrag*rhodeq*w2i**4*w3i**3 -&
        w1i*rhogsol*Kdrag*rhodeq*w3i**4*w2i**3 + rhogsol*rhogeq*k**2*cs**2*w1i**4*rhodeq*w2i*w3r**2 -&
        2*w1i*rhogsol*Kdrag*rhodeq*w2i**3*w3i**2*w3r**2 - w1i*rhogsol*Kdrag*rhodeq*w2i**4*w3r**2*w3i -&
        w1i*rhogsol*Kdrag*rhodeq*w2r**4*w3r**2*w3i - w1i*rhogsol*Kdrag*rhodeq*w2r**4*w3i**3 -&
        2*w1i*rhogsol*Kdrag*rhodeq*w3i**2*w2r**2*w2i*w3r**2 - 2*w1i*rhogsol*Kdrag*rhodeq*w2r**2*w3i*w3r**2*w2i**2 -&
        2*w1i*rhogsol*Kdrag*rhodeq*w2r**2*w2i**2*w3i**3 + rhogsol*rhogeq*k**2*cs**2*w1i**4*rhodeq*w2i**2*w3i +&
        rhogsol*rhogeq*k**2*cs**2*w1i**4*rhodeq*w2i*w3i**2 + rhogsol*rhogeq*k**4*cs**4*w1r**3*rhodeq*w2i*w3r +&
        rhogsol*rhogeq*k**4*cs**4*w1r**3*rhodeq*w2r*w3i + rhogsol*Kdrag*w1r**4*rhodeq*w3r**2*w2i**2 +&
        rhogsol*Kdrag*w1r**4*rhodeq*w3r**2*w2r**2 + rhogsol*Kdrag*w1r**4*rhodeq*w3i**2*w2i**2 +&
        rhogsol*Kdrag*w1r**4*rhodeq*w2r**2*w3i**2 - w1i*rhogsol*Kdrag*rhodeq*w2i**3*w3r**4 -&
        w1i*rhogsol*Kdrag*rhodeq*w2r**2*w2i*w3r**4 - vgsol*k**3*cs**2*rhogeq**2*w1r**4*rhodeq*w2i*w3r -&
        vgsol*k**3*cs**2*rhogeq**2*w1r**4*rhodeq*w2r*w3i + 2*w1r**3*Kdrag*rhogeq*rhodsol*w3r*w2i**2*w3i**2 +&
        2*w1r**3*Kdrag*rhogeq*rhodsol*w3r*w2r**2*w3i**2 + 2*w1r**3*Kdrag*rhogeq*rhodsol*w2r**3*w3i**2 +&
        2*w1r**3*Kdrag*rhogeq*rhodsol*w2r*w2i**2*w3i**2 + 2*w1r**3*Kdrag*rhogeq*rhodsol*w3r**3*w2i**2 +&
        2*w1r**3*Kdrag*rhogeq*rhodsol*w3r**3*w2r**2 + 2*w1r**3*Kdrag*rhogeq*rhodsol*w2r**3*w3r**2 +&
        2*w1r**3*Kdrag*rhogeq*rhodsol*w2r*w2i**2*w3r**2 + 2*w1r*w1i**2*Kdrag*rhogeq*rhodsol*w3r**3*w2r**2 +&
        2*w1r*w1i**2*Kdrag*rhogeq*rhodsol*w2r**3*w3r**2 + 2*w1r*w1i**2*Kdrag*rhogeq*rhodsol*w2r*w2i**2*w3r**2 +&
        2*w1i**3*Kdrag*rhogeq*rhodsol*w2i**2*w3i**3 + 2*w1i**3*Kdrag*rhogeq*rhodsol*w2i**3*w3i**2 +&
        2*w1r*w1i**2*Kdrag*rhogeq*rhodsol*w3r**3*w2i**2 + 2*w1i**3*Kdrag*rhogeq*rhodsol*w3r**2*w2i**3 +&
        2*w1i**3*Kdrag*rhogeq*rhodsol*w3r**2*w2r**2*w3i + 2*w1i**3*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2*w3i +&
        2*w1i**3*Kdrag*rhogeq*rhodsol*w2i*w2r**2*w3i**2 + 2*w1i**3*Kdrag*rhogeq*rhodsol*w3r**2*w2i*w2r**2 +&
        2*w1i**3*Kdrag*rhogeq*rhodsol*w2r**2*w3i**3 + 2*w1r*w1i**2*Kdrag*rhogeq*rhodsol*w3r*w2i**2*w3i**2 +&
        2*w1r*w1i**2*Kdrag*rhogeq*rhodsol*w3r*w2r**2*w3i**2 + 2*w1r*w1i**2*Kdrag*rhogeq*rhodsol*w2r**3*w3i**2 +&
        2*w1r*w1i**2*Kdrag*rhogeq*rhodsol*w2r*w2i**2*w3i**2 + 2*w1i*w1r**2*Kdrag*rhogeq*rhodsol*w2i**2*w3i**3 +&
        2*w1i*w1r**2*Kdrag*rhogeq*rhodsol*w2i*w2r**2*w3i**2 + 2*w1i*w1r**2*Kdrag*rhogeq*rhodsol*w2i**3*w3i**2 +&
        2*w1i*w1r**2*Kdrag*rhogeq*rhodsol*w2r**2*w3i**3 - 2*w1r**2*w1i**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2 -&
        2*w1r**2*w1i**2*Kdrag*rhogeq*rhodsol*w3r**2*w2r**2 - 2*w1r**2*w1i**2*Kdrag*rhogeq*rhodsol*w3i**2*w2i**2 -&
        2*w1r**2*w1i**2*Kdrag*rhogeq*rhodsol*w2r**2*w3i**2 + 2*w1i*w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i*w2r**2 +&
        2*w1i*w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**3 + 2*w1i*w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2r**2*w3i +&
        2*w1i*w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2*w3i - 2*w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2*w3i**2 -&
        2*w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2r**2*w3i**2 - 4*w1r**2*Kdrag*rhogeq*rhodsol*w3r*w2r**3*w3i**2 -&
        w1r**2*Kdrag*rhogeq*rhodsol*w2r**4*w3i**2 - 2*w1r**2*Kdrag*rhogeq*rhodsol*w3i**2*w2i**2*w2r**2 -&
        w1r**2*Kdrag*rhogeq*rhodsol*w3i**2*w2i**4 - w1r**2*Kdrag*rhogeq*rhodsol*w3i**4*w2i**2 -&
        w1r**2*Kdrag*rhogeq*rhodsol*w3i**4*w2r**2 - w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2r**4 -&
        2*w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2*w2r**2 - w1r**2*Kdrag*rhogeq*rhodsol*w3r**2*w2i**4 -&
        4*w1r**2*Kdrag*rhogeq*rhodsol*w2r*w3r*w2i**2*w3i**2 - w1r**2*Kdrag*rhogeq*rhodsol*w3r**4*w2i**2 -&
        w1r**2*Kdrag*rhogeq*rhodsol*w3r**4*w2r**2 - 4*w1r**2*Kdrag*rhogeq*rhodsol*w2r*w3r**3*w2i**2 -&
        4*w1r**2*Kdrag*rhogeq*rhodsol*w3r**3*w2r**3 + 2*w1r*Kdrag*rhogeq*rhodsol*w2r*w3i**4*w2i**2 +&
        4*w1r*Kdrag*rhogeq*rhodsol*w3r**3*w2r**2*w2i**2 + 4*w1r*Kdrag*rhogeq*rhodsol*w2r**3*w3i**2*w3r**2 +&
        4*w1r*Kdrag*rhogeq*rhodsol*w2r*w3i**2*w2i**2*w3r**2 + 2*w1r*Kdrag*rhogeq*rhodsol*w2r*w2i**2*w3r**4 +&
        2*w1r*Kdrag*rhogeq*rhodsol*w3r**3*w2i**4 - 4*w1i*w1r*Kdrag*rhogeq*rhodsol*w3r*w2r**2*w2i*w3i**2 -&
        4*w1i*w1r*Kdrag*rhogeq*rhodsol*w2r*w2i**2*w3i**3 - 4*w1i*w1r*Kdrag*rhogeq*rhodsol*w3i**3*w2r**3 -&
        4*w1i*w1r*Kdrag*rhogeq*rhodsol*w3r*w2i**3*w3i**2 + 2*w1r*Kdrag*rhogeq*rhodsol*w2r**3*w3r**4 +&
        2*w1r*Kdrag*rhogeq*rhodsol*w3r**3*w2r**4 + 2*w1r*Kdrag*rhogeq*rhodsol*w3r*w2r**4*w3i**2 +&
        2*w1r*Kdrag*rhogeq*rhodsol*w3r*w2i**4*w3i**2 + 4*w1r*Kdrag*rhogeq*rhodsol*w3r*w2r**2*w3i**2*w2i**2 +&
        2*w1r*Kdrag*rhogeq*rhodsol*w2r**3*w3i**4 - 4*w1i*w1r*Kdrag*rhogeq*rhodsol*w2i**3*w3r**3 -&
        4*w1i*w1r*Kdrag*rhogeq*rhodsol*w3i*w2i**2*w3r**2*w2r - 4*w1i*w1r*Kdrag*rhogeq*rhodsol*w3i*w2r**3*w3r**2 -&
        4*w1i*w1r*Kdrag*rhogeq*rhodsol*w3r**3*w2i*w2r**2 - vdsol*Kdrag*k*rhogeq*w1r**3*rhodeq*w3i**2*w2i**2 -&
        vdsol*Kdrag*k*rhogeq*w1r**3*rhodeq*w2r**2*w3i**2 + vdsol*Kdrag*k**3*cs**2*rhogeq*w1r**3*rhodeq*w2i*w3i -&
        vdsol*Kdrag*k*rhogeq*w1r**3*rhodeq*w3r**2*w2i**2 - vdsol*Kdrag*k**3*cs**2*rhogeq*w1r**3*rhodeq*w3r*w2r -&
        vdsol*Kdrag*k*rhogeq*w1r**3*rhodeq*w3r**2*w2r**2 - cs**2*k**2*rhogsol*Kdrag*w1r**2*rhodeq*w3r**2*w2i**2 -&
        cs**2*k**2*rhogsol*Kdrag*w1r**2*rhodeq*w3r**2*w2r**2 - cs**2*k**2*rhogsol*Kdrag*w1r**2*rhodeq*w3i**2*w2i**2 -&
        cs**2*k**2*rhogsol*Kdrag*w1r**2*rhodeq*w2r**2*w3i**2 - 2*rhogsol*rhogeq*k**2*cs**2*w1r**3*rhodeq*w3r*w2i*w3i**2&
        - 2*rhogsol*rhogeq*k**2*cs**2*w1r**3*rhodeq*w3i*w3r*w2i**2 -&
        2*rhogsol*rhogeq*k**2*cs**2*w1r**3*rhodeq*w2r**3*w3i - 2*rhogsol*rhogeq*k**2*cs**2*w1r**3*rhodeq*w3r**3*w2i +&
        2*vgsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2r*w2i*w3i -&
        2*rhogsol*rhogeq*k**2*cs**2*w1r**3*rhodeq*w2r*w2i*w3r**2 -&
        2*rhogsol*rhogeq*k**2*cs**2*w1r**3*rhodeq*w3i*w3r*w2r**2 -&
        2*rhogsol*rhogeq*k**2*cs**2*w1r**3*rhodeq*w2r*w3i*w2i**2 -&
        2*rhogsol*rhogeq*k**2*cs**2*w1r**3*rhodeq*w2r*w2i*w3i**2 -&
        2*vdsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w3r*w2i*w3i -&
        2*vdsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2r*w2i*w3i -&
        2*vgsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2r*w3r**2 -&
        2*vgsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2r**2*w3r +&
        2*vgsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w3r*w2i*w3i +&
        2*vdsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2r*w3r**2 +&
        2*vdsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2r**2*w3r -&
        2*rhogsol*k**4*cs**4*rhogeq*w1r**2*rhodeq*w3r*w2r*w2i - 2*rhogsol*k**4*cs**4*rhogeq*w1r**2*rhodeq*w3r*w2r*w3i -&
        2*rhogsol*k**4*cs**4*rhogeq*w1r**2*rhodeq*w2r**2*w3i - vgsol*k**3*cs**2*rhogeq**2*w1r**3*rhodeq*w2i*w3i**2 -&
        2*rhogsol*k**4*cs**4*rhogeq*w1r**2*rhodeq*w2i*w3r**2 + 2*vgsol*k**3*cs**2*rhogeq**2*w1r**3*rhodeq*w3r*w2r*w3i +&
        2*vgsol*k**3*cs**2*rhogeq**2*w1r**3*rhodeq*w3r*w2r*w2i + vgsol*k**3*cs**2*rhogeq**2*w1r**3*rhodeq*w2r**2*w3i +&
        w1i*vdsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2r*w3i +&
        w1i*vdsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2i*w3r + vgsol*k**3*cs**2*rhogeq**2*w1r**3*rhodeq*w2i*w3r**2 -&
        vgsol*k**3*cs**2*rhogeq**2*w1r**3*rhodeq*w2i**2*w3i - w1i*vgsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2i*w3r -&
        w1i*vgsol*Kdrag*rhogeq*k**3*cs**2*w1r**2*rhodeq*w2r*w3i - w1r*rhogsol*Kdrag*rhodeq*w2r**3*w3i**4 -&
        w1r*rhogsol*Kdrag*rhodeq*w2r*w3i**4*w2i**2 - w1r*rhogsol*Kdrag*rhodeq*w3r*w2r**4*w3i**2 -&
        w1r*rhogsol*Kdrag*rhodeq*w3r*w2i**4*w3i**2 - 2*w1r*rhogsol*Kdrag*rhodeq*w3r*w2r**2*w3i**2*w2i**2 -&
        2*w1r*rhogsol*Kdrag*rhodeq*w3r**3*w2r**2*w2i**2 - 2*w1r*rhogsol*Kdrag*rhodeq*w2r**3*w3i**2*w3r**2 -&
        2*w1r*rhogsol*Kdrag*rhodeq*w2r*w3i**2*w2i**2*w3r**2 - w1r*rhogsol*Kdrag*rhodeq*w2r**3*w3r**4 -&
        w1r*rhogsol*Kdrag*rhodeq*w2r*w2i**2*w3r**4 - w1r*rhogsol*Kdrag*rhodeq*w3r**3*w2i**4 -&
        w1r*rhogsol*Kdrag*rhodeq*w3r**3*w2r**4 - 2*w1r*w1i**2*rhogsol*Kdrag*rhodeq*w2r*w2i**2*w3r**2 -&
        2*w1r*w1i**2*rhogsol*Kdrag*rhodeq*w3r*w2i**2*w3i**2 - 2*w1r*w1i**2*rhogsol*Kdrag*rhodeq*w2r**3*w3i**2 -&
        2*w1r*w1i**2*rhogsol*Kdrag*rhodeq*w3r**3*w2i**2 - 2*w1r*w1i**2*rhogsol*Kdrag*rhodeq*w3r**3*w2r**2 +&
        w1r*cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r*w2r**2*w3i**2 + w1r*cs**2*k**2*rhogsol*Kdrag*rhodeq*w2r**3*w3i**2 +&
        w1r*cs**2*k**2*rhogsol*Kdrag*rhodeq*w2r*w2i**2*w3i**2 + w1r*cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r*w2i**2*w3i**2 -&
        2*w1r*w1i**2*rhogsol*Kdrag*rhodeq*w2r**3*w3r**2 - 2*w1r*w1i**2*rhogsol*Kdrag*rhodeq*w3r*w2r**2*w3i**2 -&
        2*w1r*w1i**2*rhogsol*Kdrag*rhodeq*w2r*w2i**2*w3i**2 + w1r*cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r**3*w2r**2 +&
        w1r*cs**2*k**2*rhogsol*Kdrag*rhodeq*w2r**3*w3r**2 + w1r*cs**2*k**2*rhogsol*Kdrag*rhodeq*w2r*w2i**2*w3r**2 -&
        2*w1r**3*rhogsol*Kdrag*rhodeq*w2r**3*w3r**2 - 2*w1r**3*rhogsol*Kdrag*rhodeq*w3r*w2i**2*w3i**2 -&
        2*w1r**3*rhogsol*Kdrag*rhodeq*w3r*w2r**2*w3i**2 - 2*w1r**3*rhogsol*Kdrag*rhodeq*w2r*w2i**2*w3i**2 -&
        2*w1i**2*vgsol*rhogeq*k*Kdrag*rhodeq*w3r**3*w2i**2 - 2*w1i**2*vgsol*rhogeq*k*Kdrag*rhodeq*w3r**3*w2r**2 -&
        2*w1i**2*vgsol*rhogeq*k*Kdrag*rhodeq*w2r**3*w3r**2 - 2*w1i**2*vgsol*rhogeq*k*Kdrag*rhodeq*w2r*w2i**2*w3r**2 -&
        2*w1i**2*vgsol*rhogeq*k*Kdrag*rhodeq*w3r*w2i**2*w3i**2 - 2*w1i**2*vgsol*rhogeq*k*Kdrag*rhodeq*w3r*w2r**2*w3i**2&
        - 2*w1i**2*vgsol*rhogeq*k*Kdrag*rhodeq*w2r**3*w3i**2 - 2*w1i**2*vgsol*rhogeq*k*Kdrag*rhodeq*w2r*w2i**2*w3i**2 +&
        w1i*k*vgsol*rhogeq**2*rhodeq*w2r*w3i**4*w2i**2 + w1i*k*vgsol*rhogeq**2*rhodeq*w2r**3*w3i**4 +&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*w1i**2*rhodeq*w3r*w2r -&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*w1i**2*rhodeq*w2i*w3i - vgsol*k**3*cs**2*rhogeq**2*w1i**4*rhodeq*w2i*w3r -&
        vgsol*k**3*cs**2*rhogeq**2*w1i**4*rhodeq*w2r*w3i + w1i*k*vgsol*rhogeq**2*rhodeq*w3r*w2r**4*w3i**2 +&
        w1i*k*vgsol*rhogeq**2*rhodeq*w3r*w2i**4*w3i**2 + 2*w1i*k*vgsol*rhogeq**2*rhodeq*w3r*w2r**2*w3i**2*w2i**2 +&
        w1i*k*vgsol*rhogeq**2*rhodeq*w2r**3*w3r**4 + w1i*k*vgsol*rhogeq**2*rhodeq*w2r*w2i**2*w3r**4 +&
        w1i*k*vgsol*rhogeq**2*rhodeq*w3r**3*w2i**4 - vdsol*Kdrag*k*w1r*rhogeq*w1i**2*rhodeq*w3i**2*w2i**2 -&
        vdsol*Kdrag*k*w1r*rhogeq*w1i**2*rhodeq*w2r**2*w3i**2 - vdsol*Kdrag*k*w1r*rhogeq*w1i**2*rhodeq*w3r**2*w2r**2 +&
        w1i*k*vgsol*rhogeq**2*rhodeq*w3r**3*w2r**4 + 2*w1i*k*vgsol*rhogeq**2*rhodeq*w3r**3*w2r**2*w2i**2 +&
        2*w1i*k*vgsol*rhogeq**2*rhodeq*w2r**3*w3i**2*w3r**2 + 2*w1i*k*vgsol*rhogeq**2*rhodeq*w2r*w3i**2*w2i**2*w3r**2 -&
        2*w1i*w1r**2*rhogsol*Kdrag*rhodeq*w2i**2*w3i**3 - 2*w1i*w1r**2*rhogsol*Kdrag*rhodeq*w2i*w2r**2*w3i**2 -&
        2*w1i*w1r**2*rhogsol*Kdrag*rhodeq*w2i**3*w3i**2 - 2*w1i*w1r**2*rhogsol*Kdrag*rhodeq*w2r**2*w3i**3 -&
        2*vgsol*w1r**2*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w2i*w3r -&
        2*vgsol*w1r**2*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w2r*w3i -&
        vdsol*Kdrag*k*w1r*rhogeq*w1i**2*rhodeq*w3r**2*w2i**2 - 2*w1i*w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i*w2r**2 -&
        2*w1i*w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**3 - 2*w1i*w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2r**2*w3i -&
        2*w1i*w1r**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**2*w3i - vgsol*k*w1r*rhogeq**2*rhodeq*w3i**4*w2i**3 -&
        2*vgsol*k*w1r*rhogeq**2*rhodeq*w2r**2*w2i**2*w3i**3 - 2*vgsol*k*w1r*rhogeq**2*rhodeq*w2r**2*w3i*w3r**2*w2i**2 -&
        vgsol*k*w1r*rhogeq**2*rhodeq*w2r**4*w3i**3 - vgsol*k*w1r*rhogeq**2*rhodeq*w2i**4*w3i**3 -&
        vgsol*k*w1r*rhogeq**2*rhodeq*w2i**4*w3r**2*w3i - vgsol*k*w1r*rhogeq**2*rhodeq*w2r**4*w3r**2*w3i -&
        vgsol*k*w1r*rhogeq**2*rhodeq*w2i*w2r**2*w3i**4 - vgsol*k*w1r*rhogeq**2*rhodeq*w2r**2*w2i*w3r**4 -&
        2*vgsol*k*w1r*rhogeq**2*rhodeq*w3i**2*w2r**2*w2i*w3r**2 - 2*vgsol*k*w1r*rhogeq**2*rhodeq*w2i**3*w3i**2*w3r**2 +&
        rhogsol*w1i*k**4*cs**4*rhogeq*w1r**2*rhodeq*w3r*w2r - rhogsol*w1i*k**4*cs**4*rhogeq*w1r**2*rhodeq*w2i*w3i +&
        2*w1i**2*rhogsol*Kdrag*rhodeq*w3i**2*w2i**2*w2r**2 - vgsol*k*w1r*rhogeq**2*rhodeq*w2i**3*w3r**4 +&
        2*w1i**2*rhogsol*Kdrag*rhodeq*w3r**2*w2r**2*w3i**2 + w1i**2*rhogsol*Kdrag*rhodeq*w3r**2*w2r**4 +&
        2*w1i**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**2*w2r**2 + w1i**2*rhogsol*Kdrag*rhodeq*w3r*w2r**3*w3i**2 +&
        w1i**2*rhogsol*Kdrag*rhodeq*w2r*w3r*w2i**2*w3i**2 + w1i**2*rhogsol*Kdrag*rhodeq*w3r**3*w2r**3 +&
        3*w1i**2*rhogsol*Kdrag*rhodeq*w2i**3*w3i*w3r**2 + w1i**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**4 +&
        3*w1i**2*rhogsol*Kdrag*rhodeq*w2i*w3i*w2r**2*w3r**2 + w1i**2*rhogsol*Kdrag*rhodeq*w3r**4*w2r**2 +&
        2*w1i**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**2*w3i**2 + w1i**2*rhogsol*Kdrag*rhodeq*w2r**4*w3i**2 +&
        w1i**2*rhogsol*Kdrag*rhodeq*w3i**2*w2i**4 + w1i**2*rhogsol*Kdrag*rhodeq*w3i**4*w2i**2 +&
        w1i**2*rhogsol*Kdrag*rhodeq*w3i**4*w2r**2 + 3*w1i**2*rhogsol*Kdrag*rhodeq*w2i*w2r**2*w3i**3 +&
        3*w1i**2*rhogsol*Kdrag*rhodeq*w2i**3*w3i**3 + 2*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w2i**3*w3i +&
        2*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w2i*w3i*w2r**2 +&
        4*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w3i**2*w2i**2 -&
        8*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2r*w2i*w3i + w1i**2*rhogsol*Kdrag*rhodeq*w2r*w3r**3*w2i**2 +&
        2*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3r**3 +&
        2*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w3i*w2i*w3r**2 +&
        4*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w3r**2*w2r**2 - w1r*vgsol*k*Kdrag**2*rhodeq*w2i**2*w3i**3 -&
        w1r*vgsol*k*Kdrag**2*rhodeq*w2i**3*w3i**2 - w1r*vgsol*k*Kdrag**2*rhodeq*w2r**2*w3i**3 +&
        2*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3r*w3i**2 +&
        2*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2r*w2i**2 +&
        2*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2r**3 + 2*w1i*vgsol*w1r*k**3*cs**2*rhogeq**2*rhodeq*w3i**3*w2i&
        - w1i*rhogsol*w1r**2*rhogeq*rhodeq*w2i**3*w3i**3 + w1i**2*rhogsol*Kdrag*rhodeq*w3r**4*w2i**2 -&
        w1i*rhogsol*w1r**2*rhogeq*rhodeq*w2i*w3i*w2r**2*w3r**2 + w1i*rhogsol*w1r**2*rhogeq*rhodeq*w3r*w2r**3*w3i**2 +&
        w1i*rhogsol*w1r**2*rhogeq*rhodeq*w2r*w3r*w2i**2*w3i**2 - w1i*rhogsol*w1r**2*rhogeq*rhodeq*w2i**3*w3i*w3r**2 -&
        w1i*rhogsol*w1r**2*rhogeq*rhodeq*w2i*w2r**2*w3i**3 - w1r*vgsol*k*Kdrag**2*rhodeq*w3r**2*w2i*w2r**2 -&
        w1r*vgsol*k*Kdrag**2*rhodeq*w3r**2*w2i**3 - w1r*vgsol*k*Kdrag**2*rhodeq*w3r**2*w2r**2*w3i -&
        w1r*vgsol*k*Kdrag**2*rhodeq*w3r**2*w2i**2*w3i - w1r*vgsol*k*Kdrag**2*rhodeq*w2i*w2r**2*w3i**2 +&
        w1i*rhogsol*w1r**2*rhogeq*rhodeq*w3r**3*w2r**3 + w1i*rhogsol*w1r**2*rhogeq*rhodeq*w2r*w3r**3*w2i**2 +&
        vgsol*Kdrag**2*k*rhodeq*w3r*w2r**2*w2i*w3i**2 + vgsol*Kdrag**2*k*rhodeq*w2r*w2i**2*w3i**3 +&
        vgsol*Kdrag**2*k*rhodeq*w3i**3*w2r**3 + vgsol*Kdrag**2*k*rhodeq*w3r**3*w2i*w2r**2 +&
        vgsol*Kdrag**2*k*rhodeq*w3i*w2r**3*w3r**2 + vgsol*Kdrag**2*k*rhodeq*w3r*w2i**3*w3i**2 -&
        2*w1r*rhogsol*k**2*cs**2*rhogeq*w1i**2*rhodeq*w3i*w3r*w2i**2 -&
        2*w1r*rhogsol*k**2*cs**2*rhogeq*w1i**2*rhodeq*w3i*w3r*w2r**2 + vgsol*Kdrag**2*k*rhodeq*w2i**3*w3r**3 -&
        2*w1r*rhogsol*k**2*cs**2*rhogeq*w1i**2*rhodeq*w2r*w2i*w3r**2 + vgsol*Kdrag**2*k*rhodeq*w3i*w2i**2*w3r**2*w2r -&
        2*w1r*rhogsol*k**2*cs**2*rhogeq*w1i**2*rhodeq*w3r**3*w2i -&
        2*w1r*rhogsol*k**2*cs**2*rhogeq*w1i**2*rhodeq*w3r*w2i*w3i**2 -&
        2*w1r*rhogsol*k**2*cs**2*rhogeq*w1i**2*rhodeq*w2r**3*w3i -&
        2*w1r*rhogsol*k**2*cs**2*rhogeq*w1i**2*rhodeq*w2r*w3i*w2i**2 -&
        2*w1r*rhogsol*k**2*cs**2*rhogeq*w1i**2*rhodeq*w2r*w2i*w3i**2 + 2*w1i*w1r*vgsol*k*Kdrag**2*rhodeq*w3r**2*w2r**2&
        + 2*w1i*w1r*vgsol*k*Kdrag**2*rhodeq*w3i**2*w2i**2 + 2*w1i*w1r*vgsol*k*Kdrag**2*rhodeq*w3r**2*w2i**2 +&
        2*w1i*w1r*vgsol*k*Kdrag**2*rhodeq*w2r**2*w3i**2 + vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w2r*w3i**3 -&
        2*vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w2r*w2i*w3i**2 +&
        2*vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w3i*w3r*w2r**2 - vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w2r**3*w3i&
        - vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w2r*w3i*w2i**2 + vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w3r*w2i**3&
        - 2*vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w3i*w3r*w2i**2 -&
        vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w3r**3*w2i + vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w3r*w2i*w2r**2 -&
        vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w3r*w2i*w3i**2 - w1i*rhogeq*k*Kdrag*vdsol*rhodeq*w3r*w2r**2*w2i*w3i**2&
        - w1i*rhogeq*k*Kdrag*vdsol*rhodeq*w2r*w2i**2*w3i**3 + vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w2r*w3i*w3r**2 -&
        w1i*rhogeq*k*Kdrag*vdsol*rhodeq*w3r*w2i**3*w3i**2 - w1i*rhogeq*k*Kdrag*vdsol*rhodeq*w3i**3*w2r**3 +&
        2*vgsol*k**3*cs**2*rhogeq**2*w1i**2*rhodeq*w2r*w2i*w3r**2 - w1i*rhogeq*k*Kdrag*vdsol*rhodeq*w3r**3*w2i*w2r**2 -&
        w1i*rhogeq*k*Kdrag*vdsol*rhodeq*w2i**3*w3r**3 - w1i*rhogeq*k*Kdrag*vdsol*rhodeq*w3i*w2i**2*w3r**2*w2r -&
        w1i*rhogeq*k*Kdrag*vdsol*rhodeq*w3i*w2r**3*w3r**2 + vdsol*Kdrag*k**3*cs**2*rhogeq*w1i**3*rhodeq*w2i*w3r +&
        vdsol*Kdrag*k**3*cs**2*rhogeq*w1i**3*rhodeq*w2r*w3i +&
        2*w1r*vgsol*rhogeq**2*k**3*cs**2*w1i**2*rhodeq*w3r*w2r*w2i -&
        w1r*vgsol*rhogeq**2*k**3*cs**2*w1i**2*rhodeq*w2i**2*w3i -&
        w1r*vgsol*rhogeq**2*k**3*cs**2*w1i**2*rhodeq*w2i*w3i**2 +&
        w1r*vgsol*rhogeq**2*k**3*cs**2*w1i**2*rhodeq*w2i*w3r**2 +&
        2*w1r*vgsol*rhogeq**2*k**3*cs**2*w1i**2*rhodeq*w3r*w2r*w3i +&
        w1r*vgsol*rhogeq**2*k**3*cs**2*w1i**2*rhodeq*w2r**2*w3i - vgsol*k**3*cs**2*rhogeq**2*w1i**3*rhodeq*w2r*w3r**2 -&
        vgsol*k**3*cs**2*rhogeq**2*w1i**3*rhodeq*w2r**2*w3r - w1i*cs**2*k**2*rhogsol*Kdrag*rhodeq*w2i**3*w3i**2 +&
        vgsol*k**3*cs**2*rhogeq**2*w1i**3*rhodeq*w2i**2*w3r + 2*vgsol*k**3*cs**2*rhogeq**2*w1i**3*rhodeq*w2r*w2i*w3i +&
        vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r**3*w2r**2 - vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r**3*w2i**2 +&
        vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r**3*w3r**2 - w1i*cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**2*w3i -&
        w1i*cs**2*k**2*rhogsol*Kdrag*rhodeq*w2r**2*w3i**3 - w1i*cs**2*k**2*rhogsol*Kdrag*rhodeq*w2i**2*w3i**3 -&
        w1i*cs**2*k**2*rhogsol*Kdrag*rhodeq*w2i*w2r**2*w3i**2 +&
        4*w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w2i*w3i**2 -&
        w1i*cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i*w2r**2 - w1i*cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r**2*w2r**2*w3i +&
        w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r*w2i*w2r**2 +&
        4*w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3i*w3r*w2i**2 +&
        w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w3i*w3r**2 +&
        w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r*w2i*w3i**2 - w1i*cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r**2*w2i**3 +&
        vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i**2*w3i**3 + 2*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3i**4*w2i +&
        w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r**3*w2i + w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r*w2i**3 +&
        w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r**3*w3i + w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w3i*w2i**2 +&
        w1i*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w3i**3 + vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2i**3*w3i**2 +&
        vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2r**2*w2i*w3i**2 + 2*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i*w3r*w2i**4 +&
        vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r**3*w2i*w2r**2 + vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i*w2r**3*w3r**2 +&
        4*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i*w3i**2*w3r**2 + 2*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i*w3r**4&
        + vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r*w2i**2*w3i**2 + vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2i**3*w3r**3 +&
        vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i*w2i**2*w3r**2*w2r +&
        4*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r*w3i*w2r**2*w2i**2 + 2*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r**4*w3i*w3r&
        + vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3i**3*w2r**3 + 2*vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r*w2i*w3r**2*w3i -&
        vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r*w2r**2*w3i**2 - vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r*w2i**2*w3r**2&
        + vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r*w2i**2*w3i**2 +&
        2*vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i*w3i*w3r*w2r**2 +&
        2*vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i**3*w3i*w3r + vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r**3*w3i**2 +&
        2*vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r*w2i*w3i**3 - vdsol*Kdrag**2*k*rhodeq*w3r*w2r**2*w2i*w3i**2 -&
        vdsol*Kdrag**2*k*rhodeq*w2r*w2i**2*w3i**3 - vdsol*Kdrag**2*k*rhodeq*w3r**3*w2i*w2r**2 -&
        vdsol*Kdrag**2*k*rhodeq*w3i*w2i**2*w3r**2*w2r - vdsol*Kdrag**2*k*rhodeq*w3i*w2r**3*w3r**2 -&
        vdsol*Kdrag**2*k*rhodeq*w3r*w2i**3*w3i**2 - vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r**3*w2r**2 +&
        vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w3r**3*w2i**2 - vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r**3*w3r**2 -&
        w1i*w1r**2*vgsol*k*rhogeq**2*rhodeq*w3r**3*w2r**2 - w1i*w1r**2*vgsol*k*rhogeq**2*rhodeq*w2r**3*w3r**2 -&
        w1i*w1r**2*vgsol*k*rhogeq**2*rhodeq*w2r*w2i**2*w3r**2 +&
        2*w1i*w1r*vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i*w3i**2 +&
        2*w1i*w1r*vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r**2*w3i +&
        2*w1i*w1r*vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i**2*w3i -&
        w1i*w1r**2*vgsol*k*rhogeq**2*rhodeq*w3r*w2i**2*w3i**2 - w1i*w1r**2*vgsol*k*rhogeq**2*rhodeq*w3r*w2r**2*w3i**2 -&
        w1i*w1r**2*vgsol*k*rhogeq**2*rhodeq*w2r**3*w3i**2 - w1i*w1r**2*vgsol*k*rhogeq**2*rhodeq*w2r*w2i**2*w3i**2 -&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2i*w3i*w2r**2 - w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2i**3*w3i -&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3i**2*w2i**2 - w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r**2*w3i**2&
        + 2*w1i*w1r*vgsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i*w3r**2 -&
        4*w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r*w2r*w2i*w3i -&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3i**3*w2i - w1i*w1r**2*vgsol*k*rhogeq**2*rhodeq*w3r**3*w2i**2 -&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3i*w2i*w3r**2 -&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r**2*w2i**2 +&
        3*w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r**2*w2r**2 - vgsol*k*rhogeq**2*w1i**3*rhodeq*w2r*w2i**2*w3i**2 +&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w3r**3 - vgsol*k*rhogeq**2*w1i**3*rhodeq*w2r**3*w3r**2 -&
        vgsol*k*rhogeq**2*w1i**3*rhodeq*w2r*w2i**2*w3r**2 - vgsol*k*rhogeq**2*w1i**3*rhodeq*w3r*w2i**2*w3i**2 -&
        vgsol*k*rhogeq**2*w1i**3*rhodeq*w3r*w2r**2*w3i**2 - vdsol*Kdrag**2*k*rhodeq*w2i**3*w3r**3 -&
        vgsol*k*rhogeq**2*w1i**3*rhodeq*w3r**3*w2r**2 - vgsol*k*rhogeq**2*w1i**3*rhodeq*w2r**3*w3i**2 +&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w3r*w3i**2 +&
        w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r*w2r*w2i**2 + w1r*vgsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r*w2r**3 -&
        vgsol*k*rhogeq**2*w1i**3*rhodeq*w3r**3*w2i**2 - vdsol*Kdrag**2*k*rhodeq*w3i**3*w2r**3 +&
        w1r*vdsol*k*Kdrag**2*rhodeq*w2i**2*w3i**3 + w1r*vdsol*k*Kdrag**2*rhodeq*w2i*w2r**2*w3i**2 +&
        w1r*vdsol*k*Kdrag**2*rhodeq*w2i**3*w3i**2 + w1r*vdsol*k*Kdrag**2*rhodeq*w3r**2*w2r**2*w3i +&
        w1r*vdsol*k*Kdrag**2*rhodeq*w3r**2*w2i**2*w3i + w1r*vdsol*k*Kdrag**2*rhodeq*w2r**2*w3i**3 +&
        w1r*vdsol*k*Kdrag**2*rhodeq*w3r**2*w2i*w2r**2 + w1r*vdsol*k*Kdrag**2*rhodeq*w3r**2*w2i**3 +&
        cs**2*k**2*rhogsol*Kdrag*rhodeq*w2i**3*w3i*w3r**2 - cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r*w2r**3*w3i**2 -&
        cs**2*k**2*rhogsol*Kdrag*rhodeq*w2r*w3r*w2i**2*w3i**2 - cs**2*k**2*rhogsol*Kdrag*rhodeq*w3r**3*w2r**3 +&
        cs**2*k**2*rhogsol*Kdrag*rhodeq*w2i*w3i*w2r**2*w3r**2 - w1r**4*Kdrag*rhogeq*rhodsol*w2r**2*w3i**2 +&
        vgsol*Kdrag*k**3*cs**2*rhogeq*w1r**3*rhodeq*w3r*w2r - w1r**4*Kdrag*rhogeq*rhodsol*w3i**2*w2i**2 -&
        vgsol*Kdrag*k**3*cs**2*rhogeq*w1r**3*rhodeq*w2i*w3i - cs**2*k**2*rhogsol*Kdrag*rhodeq*w2r*w3r**3*w2i**2 +&
        cs**2*k**2*rhogsol*Kdrag*rhodeq*w2i*w2r**2*w3i**3 + cs**2*k**2*rhogsol*Kdrag*rhodeq*w2i**3*w3i**3 +&
        2*w1i**2*rhogsol*rhogeq*rhodeq*w2i**3*w3i**2*w3r**2 + w1i**2*rhogsol*rhogeq*rhodeq*w2i**4*w3r**2*w3i +&
        w1i**2*rhogsol*rhogeq*rhodeq*w2r**4*w3r**2*w3i + w1i**2*rhogsol*rhogeq*rhodeq*w2r**4*w3i**3 +&
        2*w1i**2*rhogsol*rhogeq*rhodeq*w3i**2*w2r**2*w2i*w3r**2 +&
        2*w1i**2*rhogsol*rhogeq*rhodeq*w2r**2*w3i*w3r**2*w2i**2 + w1i**2*rhogsol*rhogeq*rhodeq*w2i*w2r**2*w3i**4 +&
        w1i**2*rhogsol*rhogeq*rhodeq*w2i**4*w3i**3 + w1i**2*rhogsol*rhogeq*rhodeq*w3i**4*w2i**3 +&
        2*w1i**2*rhogsol*rhogeq*rhodeq*w2r**2*w2i**2*w3i**3 - w1r**4*Kdrag*rhogeq*rhodsol*w3r**2*w2i**2 -&
        w1r**4*Kdrag*rhogeq*rhodsol*w3r**2*w2r**2 + w1r**3*vgsol*k*rhogeq**2*rhodeq*w2i**2*w3i**3 +&
        w1r**3*vgsol*k*rhogeq**2*rhodeq*w2i*w2r**2*w3i**2 + w1r**3*vgsol*k*rhogeq**2*rhodeq*w2i**3*w3i**2 +&
        w1r**3*vgsol*k*rhogeq**2*rhodeq*w2r**2*w3i**3 + w1i**2*rhogsol*rhogeq*rhodeq*w2i**3*w3r**4 +&
        w1i**2*rhogsol*rhogeq*rhodeq*w2r**2*w2i*w3r**4 + w1r**3*vgsol*k*rhogeq**2*rhodeq*w3r**2*w2i*w2r**2 +&
        w1r**3*vgsol*k*rhogeq**2*rhodeq*w3r**2*w2r**2*w3i + w1r**3*vgsol*k*rhogeq**2*rhodeq*w3r**2*w2i**2*w3i +&
        w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*w1i**2*rhodeq*w2i*w3i +&
        2*w1r**2*rhogsol*rhogeq*k**2*cs**2*w1i**2*rhodeq*w2i*w3r**2 +&
        2*w1r**2*rhogsol*rhogeq*k**2*cs**2*w1i**2*rhodeq*w2i**2*w3i -&
        w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*w1i**2*rhodeq*w3r*w2r - 2*rhogsol*rhogeq*k**4*cs**4*w1i**2*rhodeq*w3r*w2r*w3i&
        + 2*rhogsol*rhogeq*k**4*cs**4*w1i**2*rhodeq*w2i**2*w3i + 2*rhogsol*rhogeq*k**4*cs**4*w1i**2*rhodeq*w2i*w3i**2 +&
        2*w1i*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2i**2*w3i**3 + 2*w1i*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2i*w2r**2*w3i**2 +&
        2*w1i*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2i**3*w3i**2 + 2*w1i*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2i**2*w3i -&
        2*rhogsol*rhogeq*k**4*cs**4*w1i**2*rhodeq*w3r*w2r*w2i +&
        2*w1r**2*rhogsol*rhogeq*k**2*cs**2*w1i**2*rhodeq*w2i*w3i**2 +&
        2*w1r**2*rhogsol*rhogeq*k**2*cs**2*w1i**2*rhodeq*w2r**2*w3i +&
        2*w1i*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2i**3 + 2*w1i*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2r**2*w3i +&
        2*w1i*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2r**2*w3i**3 + 2*Kdrag*vgsol*rhogeq*k**3*cs**2*w1i**2*rhodeq*w2r*w2i*w3i&
        + 2*Kdrag*vgsol*rhogeq*k**3*cs**2*w1i**2*rhodeq*w2r*w3i**2 -&
        2*Kdrag*vdsol*rhogeq*k**3*cs**2*w1i**2*rhodeq*w2r*w3i**2 +&
        2*Kdrag*vgsol*rhogeq*k**3*cs**2*w1i**2*rhodeq*w2i**2*w3r -&
        2*Kdrag*vdsol*rhogeq*k**3*cs**2*w1i**2*rhodeq*w2i**2*w3r -&
        2*Kdrag*vdsol*rhogeq*k**3*cs**2*w1i**2*rhodeq*w3r*w2i*w3i + w1r**3*vgsol*k*rhogeq**2*rhodeq*w3r**2*w2i**3 -&
        2*Kdrag*vdsol*rhogeq*k**3*cs**2*w1i**2*rhodeq*w2r*w2i*w3i +&
        2*Kdrag*vgsol*rhogeq*k**3*cs**2*w1i**2*rhodeq*w3r*w2i*w3i +&
        2*w1i*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2i*w2r**2 - w1i*vgsol*Kdrag**2*k*rhodeq*w3r**3*w2i**2 -&
        w1i*vgsol*Kdrag**2*k*rhodeq*w3r**3*w2r**2 - w1i*vgsol*Kdrag**2*k*rhodeq*w2r**3*w3r**2 -&
        w1i*vgsol*Kdrag**2*k*rhodeq*w2r*w2i**2*w3r**2 - w1i*vgsol*Kdrag**2*k*rhodeq*w2r*w2i**2*w3i**2 -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r**3*w2r**2 -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i**2*w3r**2 -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i*w3r**2*w3i -&
        2*w1i*w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2r*w3r**2 - 2*w1i*w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2i**2*w3r -&
        2*w1i*w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2r**2*w3r - 2*w1i*w1r*rhogsol*rhogeq*k**4*cs**4*rhodeq*w2r*w3i**2 -&
        w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w2i**3*w3i - 3*w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3i**2*w2i**2 +&
        w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w2r**2*w3i**2 + 2*w1r**2*rhogsol*Kdrag*w1i**2*rhodeq*w3r**2*w2r**2 +&
        2*w1r**2*rhogsol*Kdrag*w1i**2*rhodeq*w3i**2*w2i**2 + 2*w1r**2*rhogsol*Kdrag*w1i**2*rhodeq*w2r**2*w3i**2 -&
        w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3i*w2i*w3r**2 + w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r**2*w2r**2 +&
        w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w2r*w3r*w3i**2 + w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r**2*w2i**2 +&
        4*w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r*w2r*w2i*w3i + w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r*w2r**3 +&
        w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r*w2r*w2i**2 - w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3i**3*w2i -&
        w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w2i*w3i*w2r**2 - 2*w1i*w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i*w3r**2&
        - 2*w1i*w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i**2*w3i -&
        2*w1i*w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2i*w3i**2 -&
        2*w1i*w1r*vdsol*Kdrag*rhogeq*k**3*cs**2*rhodeq*w2r**2*w3i + 2*w1r**2*rhogsol*Kdrag*w1i**2*rhodeq*w3r**2*w2i**2&
        + rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r**2*w2i**2*w3i - 2*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3i*w3r*w2r**3 -&
        2*rhogsol*k**4*cs**4*rhogeq*rhodeq*w2r*w3r*w2i*w3i**2 + rhogsol*k**4*cs**4*rhogeq*rhodeq*w2i**3*w3i**2 +&
        rhogsol*k**4*cs**4*rhogeq*rhodeq*w2i*w2r**2*w3i**2 - rhogsol*k**4*cs**4*rhogeq*rhodeq*w2r**2*w3i**3 +&
        w1i*rhogsol*k**4*cs**4*rhogeq*rhodeq*w2r*w3r**3 + w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i**4*w2i -&
        2*rhogsol*k**4*cs**4*rhogeq*rhodeq*w2r*w3r**3*w2i - rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r**2*w2i*w2r**2 -&
        rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r**2*w2r**2*w3i - 2*rhogsol*k**4*cs**4*rhogeq*rhodeq*w3i*w3r*w2i**2*w2r +&
        w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i*w2r**4 + w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i*w2i**4 +&
        2*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i*w2r**2*w2i**2 +&
        w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3r**4*w2i + 4*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w2r*w3r**3*w2i +&
        2*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i**2*w2i*w3r**2 + vdsol*Kdrag*k*w1r*rhogeq*rhodeq*w2i*w2r**2*w3i**3&
        + rhogsol*k**4*cs**4*rhogeq*rhodeq*w2i**2*w3i**3 - vdsol*Kdrag*k*w1r*rhogeq*rhodeq*w2r*w3r*w2i**2*w3i**2 +&
        4*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i*w3r*w2i**2*w2r +&
        4*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w3i*w3r*w2r**3 +&
        4*w1r**2*rhogsol*rhogeq*k**2*cs**2*rhodeq*w2r*w3r*w2i*w3i**2 - rhogsol*k**4*cs**4*rhogeq*rhodeq*w3r**2*w2i**3 -&
        vdsol*Kdrag*k*w1r*rhogeq*rhodeq*w3r**3*w2r**3 + w1i*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3i**2 -&
        vdsol*Kdrag*k*w1r*rhogeq*rhodeq*w2r*w3r**3*w2i**2 + vdsol*Kdrag*k*w1r*rhogeq*rhodeq*w2i**3*w3i*w3r**2 +&
        vdsol*Kdrag*k*w1r*rhogeq*rhodeq*w2i*w3i*w2r**2*w3r**2 - vdsol*Kdrag*k*w1r*rhogeq*rhodeq*w3r*w2r**3*w3i**2 +&
        vdsol*Kdrag*k*w1r*rhogeq*rhodeq*w2i**3*w3i**3 + Kdrag*rhogeq*k*vgsol*rhodeq*w2r**3*w3i**4 +&
        Kdrag*rhogeq*k*vgsol*rhodeq*w2r*w3i**4*w2i**2 + 2*Kdrag*rhogeq*k*vgsol*rhodeq*w3r*w2r**2*w3i**2*w2i**2 -&
        w1i*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3r**2 -&
        w1i*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r**2*w3r +&
        w1i*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2i**2*w3r +&
        2*w1i*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2i*w3i +&
        2*w1i*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i*w3i +&
        2*Kdrag*rhogeq*k*vgsol*rhodeq*w2r*w3i**2*w2i**2*w3r**2 + Kdrag*rhogeq*k*vgsol*rhodeq*w3r*w2r**4*w3i**2 +&
        Kdrag*rhogeq*k*vgsol*rhodeq*w3r*w2i**4*w3i**2 + Kdrag*rhogeq*k*vgsol*rhodeq*w3r**3*w2i**4 +&
        Kdrag*rhogeq*k*vgsol*rhodeq*w3r**3*w2r**4 + 2*Kdrag*rhogeq*k*vgsol*rhodeq*w2r**3*w3i**2*w3r**2 +&
        Kdrag*rhogeq*k*vgsol*rhodeq*w2r*w2i**2*w3r**4 + 2*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**3*w2r**2*w2i**2 -&
        vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3i**3 + vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3i*w2i**2 +&
        2*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i*w3i**2 -&
        2*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w3i*w3r*w2r**2 + vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r**3*w3i&
        - vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2i**3 +&
        2*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w3i*w3r*w2i**2 +&
        2*w1r*cs**2*k**2*rhogsol*w1i*rhogeq*rhodeq*w2r**3*w3i**2 -&
        vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3i*w3r**2 -&
        vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2i*w2r**2 + vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w3r**3*w2i -&
        2*vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i*w3r**2 -&
        vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w2i**2*w3i**2 + vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r*w2r**2*w3i**2&
        - vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w3r*w2i**2*w3i**2 -&
        2*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2i*w3i*w3r*w2r**2 -&
        2*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w2i*w3r**2*w3i -&
        2*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2i**3*w3i*w3r - vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r**3*w3i**2 -&
        2*vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w2i*w3i**3 + vdsol*Kdrag*k**3*cs**2*rhogeq*rhodeq*w2r*w2i**2*w3r**2&
        + 2*vgsol*k**3*cs**2*rhogeq**2*w1i**3*rhodeq*w3r*w2i*w3i + vgsol*k**3*cs**2*rhogeq**2*w1i**3*rhodeq*w2r*w3i**2&
        + w1i*Kdrag*rhogeq*k*vgsol*rhodeq*w3i**3*w2r**3 + w1i*Kdrag*rhogeq*k*vgsol*rhodeq*w3i*w2i**2*w3r**2*w2r +&
        w1i*Kdrag*rhogeq*k*vgsol*rhodeq*w3r*w2i**3*w3i**2 + w1i*Kdrag*rhogeq*k*vgsol*rhodeq*w3r*w2r**2*w2i*w3i**2 +&
        w1i*Kdrag*rhogeq*k*vgsol*rhodeq*w2r*w2i**2*w3i**3 - 2*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3i**2*w2i**2*w2r**2 +&
        w1i*Kdrag*rhogeq*k*vgsol*rhodeq*w2i**3*w3r**3 + w1i*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**3*w2i*w2r**2 +&
        w1i*Kdrag*rhogeq*k*vgsol*rhodeq*w3i*w2r**3*w3r**2 - w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3i**2*w2i**4 -&
        w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3i**4*w2i**2 - 3*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2i*w3i*w2r**2*w3r**2 -&
        2*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2i**2*w2r**2 - w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2r**4*w3i**2 -&
        w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3i**4*w2r**2 - 3*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2i*w2r**2*w3i**3 -&
        3*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2i**3*w3i**3 + vgsol*w1r**2*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2i*w3i**2 +&
        Kdrag*rhogeq*k*vgsol*rhodeq*w2r**3*w3r**4 - 2*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2i**2*w3i**2 -&
        2*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2r**2*w3i**2 - w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2r*w3r*w2i**2*w3i**2&
        - 3*w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2i**3*w3i*w3r**2 - w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2i**4 -&
        w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**2*w2r**4 - w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**4*w2r**2 -&
        w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**3*w2r**3 - w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w2r*w3r**3*w2i**2 +&
        rhogsol*rhogeq*k**2*cs**2*w1r**4*rhodeq*w2i*w3r**2 + rhogsol*rhogeq*k**2*cs**2*w1r**4*rhodeq*w2i**2*w3i +&
        rhogsol*rhogeq*k**2*cs**2*w1r**4*rhodeq*w2i*w3i**2 + rhogsol*rhogeq*k**2*cs**2*w1r**4*rhodeq*w2r**2*w3i -&
        w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r**4*w2i**2 - w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r**4*w3r -&
        w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2i**4 - 2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2i**3*w3i*w3r -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2i**2*w3r*w2r**2 -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2i*w3i*w3r*w2r**2 -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w3r*w2r**2*w3i**2 - w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3i**4 -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w2i*w3i**3 -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3i**2*w3r**2 -&
        2*w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r**3*w3r**2 - w1i*vgsol*k**3*cs**2*rhogeq**2*rhodeq*w2r*w3r**4 -&
        w1i*vgsol*Kdrag**2*k*rhodeq*w3r*w2i**2*w3i**2 - w1i*vgsol*Kdrag**2*k*rhodeq*w3r*w2r**2*w3i**2 -&
        w1i*vgsol*Kdrag**2*k*rhodeq*w2r**3*w3i**2 - w1r*Kdrag*rhogeq*k*vgsol*rhodeq*w3r*w2r**3*w3i**2)/(w1i**2 -&
        2*w3i*w1i + w3r**2 + w1r**2 + w3i**2 - 2*w3r*w1r)/(w3r**2 + w3i**2)/Kdrag/(w2r**2 + w1r**2 + w2i**2 - 2*w2i*w1i&
        - 2*w2r*w1r + w1i**2)/rhogeq/(w2i**2 + w2r**2)
 endif

 print*,'w1 = ',w1r,w1i
 print*,'w2 = ',w2r,w2i
 print*,'w3 = ',w3r,w3i
 print*,' vgas1 = ',vg1r,vg1i
 print*,' vgas2 = ',vg2r,vg2i
 print*,' vgas3 = ',vg3r,vg3i
 if (Kdrag > 0.) then
    print*,' vdust1 = ',vd1r,vd1i
    print*,' vdust2 = ',vd2r,vd2i
    print*,' vdust3 = ',vd3r,vd3i
 endif
 print*,' rhogas1 = ',rhog1r,rhog1i
 print*,' rhogas2 = ',rhog2r,rhog2i
 print*,' rhogas3 = ',rhog3r,rhog3i
 if (Kdrag > 0.) then
    print*,' rhodust1 = ',rhod1r,rhod1i
    print*,' rhodust2 = ',rhod2r,rhod2i
    print*,' rhodust3 = ',rhod3r,rhod3i
 endif
!-------------------------------
! F I N A L  S O L U T I O N
!-------------------------------
 do i=1,size(xplot)
    xk =  2.*pi/lambda*(xplot(i)-x0)
    arg1 = xk - w1r*time
    arg2 = xk - w2r*time
    arg3 = xk - w3r*time
    vgas = vgeq &
          + vg1r*exp(w1i*time)*cos(arg1) - vg1i*exp(w1i*time)*sin(arg1) &
          + vg2r*exp(w2i*time)*cos(arg2) - vg2i*exp(w2i*time)*sin(arg2) &
          + vg3r*exp(w3i*time)*cos(arg3) - vg3i*exp(w3i*time)*sin(arg3)

    vdust = vdeq &
            + vd1r*exp(w1i*time)*cos(arg1) - vd1i*exp(w1i*time)*sin(arg1) &
            + vd2r*exp(w2i*time)*cos(arg2) - vd2i*exp(w2i*time)*sin(arg2) &
            + vd3r*exp(w3i*time)*cos(arg3) - vd3i*exp(w3i*time)*sin(arg3)

    rhogas = rhogeq &
            + rhog1r*exp(w1i*time)*cos(arg1) - rhog1i*exp(w1i*time)*sin(arg1) &
            + rhog2r*exp(w2i*time)*cos(arg2) - rhog2i*exp(w2i*time)*sin(arg2) &
            + rhog3r*exp(w3i*time)*cos(arg3) - rhog3i*exp(w3i*time)*sin(arg3)

    rhodust = rhodeq &
             + rhod1r*exp(w1i*time)*cos(arg1) - rhod1i*exp(w1i*time)*sin(arg1) &
             + rhod2r*exp(w2i*time)*cos(arg2) - rhod2i*exp(w2i*time)*sin(arg2) &
             + rhod3r*exp(w3i*time)*cos(arg3) - rhod3i*exp(w3i*time)*sin(arg3)

    vdusto(i) = vdust
    vgaso(i)  = vgas
    rhodusto(i) = rhodust
    rhogaso(i)  = rhogas
 enddo

 return
end subroutine exact_dustywave

end module dustywaves
