!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: eos_helmholtz
!
!  DESCRIPTION:
!   The equation state based on the Helmholtz free energy. 
!   "Perfect thermodynamic consistency."
!   Primarily used to model degenerate matter in white dwarfs.
!
!  REFERENCES:
!   Timmes & Swesty (2000), ApJS, 126, 501-516.
!
!  OWNER: Terrence Tricco
!
!
!  RUNTIME PARAMETERS: None
!
!+
!--------------------------------------------------------------------------
module eos_helmholtz
 implicit none

! subroutines to read/initialise tables, and get pressure/sound speed
 public :: eos_helmholtz_init
 public :: eos_helmholtz_write_inopt
 public :: eos_helmholtz_set_relaxflag
 public :: eos_helmholtz_pres_sound          ! performs iterations, called by eos.F90
 public :: eos_helmholtz_compute_pres_sound  ! the actual eos calculation
 public :: eos_helmholtz_cv_dpresdt

 private :: psi0, dpsi0, ddpsi0
 private :: psi1, dpsi1, ddpsi1
 private :: psi2, dpsi2, ddpsi2
 private :: xpsi0, xdpsi0
 private :: xpsi1, xdpsi1
 private :: h3, h5

 integer, public :: relaxflag = 1

 integer, parameter :: imax = 271
 integer, parameter :: jmax = 101

 ! these set the mixture of species
 ! currently hard-coded to 50/50 carbon-oxygen
 ! abar = [sum (mass fraction / # nucleons) ]^(-1)
 real, parameter :: abar   = 13.7142857
 ! zbar = [sum (mass fraction / # protons) ]^(-1)
 real, parameter :: zbar   = 6.857142857

! sizes of the tables
! normal table, big table, bigger table, denser bigger table


! original
!      parameter        (imax = 211, jmax = 71)

! standard
!      parameter        (imax = 271, jmax = 101)

! twice as dense
!      parameter        (imax = 541, jmax = 201)

! half as dense
!      parameter        (imax = 136, jmax = 51)



! for the electrons
! density and temperature
 real :: tlo, tstpi, dlo, dstpi
 real :: d(imax), t(jmax)

! for storing the differences
 real :: dt_sav(jmax),  dt2_sav(jmax), &
         dti_sav(jmax), dt2i_sav(jmax), dt3i_sav(jmax), &
         dd_sav(imax),  dd2_sav(imax), &
         ddi_sav(imax), dd2i_sav(imax), dd3i_sav(imax)

! for the helmholtz free energy tables
 real :: f(imax,jmax),    fd(imax,jmax),   ft(imax,jmax), &
         fdd(imax,jmax),  ftt(imax,jmax),  fdt(imax,jmax), &
         fddt(imax,jmax), fdtt(imax,jmax), fddtt(imax,jmax)

! for the pressure derivative with density tables
 real :: dpdf(imax,jmax), dpdfd(imax,jmax), dpdft(imax,jmax), dpdfdt(imax,jmax)

! for chemical potential tables
 real :: ef(imax,jmax), efd(imax,jmax), eft(imax,jmax), efdt(imax,jmax)

! for the number density tables
 real :: xf(imax,jmax), xfd(imax,jmax), xft(imax,jmax), xfdt(imax,jmax)

 contains


!----------------------------------------------------------------
!+
!  initialise the tables used by the Helmholtz equation of state
!+
!----------------------------------------------------------------
 subroutine eos_helmholtz_init()
   use io, only:fatal
   use datafiles, only:find_phantom_datafile
   character(len=80) :: filename
   integer :: i, j
   integer :: ierr
   real    :: tsav, dsav, dth, dt2, dti, dt2i, dt3i, &
                       dd, dd2, ddi, dd2i, dd3i
   real    :: thi, tstp, dhi, dstp

   ! check that the relaxflag is sensible, set to relax if not
   if (relaxflag /= 0 .and. relaxflag /= 1) then
      call eos_helmholtz_set_relaxflag(1)
   endif

   ! find the table datafile
   filename = find_phantom_datafile('helm_table.dat', 'eos/helmholtz')
   open(unit=19,file=trim(filename),status='old',iostat=ierr)
   if (ierr /= 0) then
       call fatal('eos_helmholtz', 'could not find helm_table.dat')
   endif

   ! for standard table limits
   tlo   = 3.0d0
   thi   = 13.0d0
   tstp  = (thi - tlo)/float(jmax-1)
   tstpi = 1.0d0/tstp
   dlo   = -12.0d0
   dhi   = 15.0d0
   dstp  = (dhi - dlo)/float(imax-1)
   dstpi = 1.0d0/dstp

   ! read the helmholtz free energy and its derivatives
   do j=1,jmax
      tsav = tlo + (j-1)*tstp
      t(j) = 10.0d0**(tsav)
      do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                  fddt(i,j),fdtt(i,j),fddtt(i,j)
      enddo
   enddo


   ! read the pressure derivative with density table
   do j=1,jmax
      do i=1,imax
         read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
      enddo
   enddo

   ! read the electron chemical potential table
   do j=1,jmax
      do i=1,imax
         read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
      enddo
   enddo

   ! read the number density table
   do j=1,jmax
      do i=1,imax
         read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
      enddo
   enddo

   ! close the file and write a summary message
   close(unit=19)

   ! construct the temperature and density deltas and their inverses
   do j=1,jmax-1
      dth          = t(j+1) - t(j)
      dt2         = dth * dth
      dti         = 1.0d0/dth
      dt2i        = 1.0d0/dt2
      dt3i        = dt2i*dti
      dt_sav(j)   = dth
      dt2_sav(j)  = dt2
      dti_sav(j)  = dti
      dt2i_sav(j) = dt2i
      dt3i_sav(j) = dt3i
   end do
   do i=1,imax-1
      dd          = d(i+1) - d(i)
      dd2         = dd * dd
      ddi         = 1.0d0/dd
      dd2i        = 1.0d0/dd2
      dd3i        = dd2i*ddi
      dd_sav(i)   = dd
      dd2_sav(i)  = dd2
      ddi_sav(i)  = ddi
      dd2i_sav(i) = dd2i
      dd3i_sav(i) = dd3i
   enddo

 end subroutine eos_helmholtz_init



!----------------------------------------------------------------
!+
!  write options to the input file (currently only relaxflag)
!+
!----------------------------------------------------------------
 subroutine eos_helmholtz_write_inopt(iunit)
   use infile_utils, only:write_inopt
   integer, intent(in) :: iunit

   call write_inopt(relaxflag, 'relaxflag', '0=evolve, 1=relaxation on (keep T const)', iunit)

 end subroutine eos_helmholtz_write_inopt

 
!----------------------------------------------------------------
!+
!  set the relaxflag based on input file read
!  
!  called by eos_read_inopt in eos.F90
!+
!----------------------------------------------------------------
 subroutine eos_helmholtz_set_relaxflag(tmp)
   use io, only:fatal
   integer, intent(in) :: tmp
   character(len=30), parameter  :: label = 'read_options_eos_helmholtz'

   relaxflag = tmp

   if (relaxflag /= 0 .and. relaxflag /= 1) call fatal(label, 'relax flag incorrect, try using 0 (evolve) or 1 (relaxation)')

 end subroutine eos_helmholtz_set_relaxflag



!----------------------------------------------------------------
!+
!  This is called by eos.F90 to get P/rho and sound speed,
!  performing iterations until ue and T converge.
!  
!  The density passed in is in code units, and it returns values
!  also in code units.
!+
!----------------------------------------------------------------
 subroutine eos_helmholtz_pres_sound(tempi,rhoi,ponrhoi,spsoundi,eni)
   use units,   only:unit_density,unit_pressure,unit_ergg,unit_velocity
   real, intent(inout) :: tempi
   real, intent(in)    :: rhoi
   real, intent(out)   :: ponrhoi
   real, intent(out)   :: spsoundi
   real, intent(inout) :: eni
   integer, parameter  :: maxiter = 10
   real,    parameter  :: tol = 0.05  ! temperature convergence
   integer :: itercount
   real :: cgsrhoi, cgspresi, cgsspsoundi, cgseni, cgseni_eos, cgsdendti
   real :: errorp, rel
   real :: tnew, tprev

   cgsrhoi = rhoi * unit_density
   
   call eos_helmholtz_compute_pres_sound(tempi, cgsrhoi, cgspresi, cgsspsoundi, cgseni_eos, cgsdendti)

   ! relaxation: 
   ! constant temperature, set internal energy of particles to result from eos
   if (relaxflag == 1) then
       eni = cgseni_eos / unit_ergg
  
   ! dynamical evolution:
   ! ue is evolved in time, iterate eos to solve for temperature when eos ue converges with particle ue
   elseif (relaxflag == 0) then
       
       cgseni = eni * unit_ergg

       ! Newton-Raphson iterations
       tnew = tempi - (cgseni_eos - cgseni) / cgsdendti

       ! disallow large temperature changes 
       if (tnew > 2.0 * tempi) then
           tnew = 2.0 * tempi
       endif
       if (tnew < 0.5 * tempi) then
           tnew = 0.5 * tempi
       endif

       ! temperature and density limits are given in section 2.3 of Timmes & Swesty (2000)
       if (tnew > 1.0e11) then
           tnew = 1.0e11
       endif
       if (tnew < 1.0e4) then
           tnew = 1.0e4
       endif

       itercount = 0
       iterations: do while ((errorp > tol) .and. (itercount < maxiter))

           itercount = itercount + 1

           ! get new pressure, sound speed, energy for this temperature and density
           call eos_helmholtz_compute_pres_sound(tnew, cgsrhoi, cgspresi, cgsspsoundi, cgseni_eos, cgsdendti)
  
           ! store temperature of previous iteration
           tprev = tnew

           ! iterate to new temperature
           tnew = tnew - (cgseni_eos - cgseni) / cgsdendti

           ! disallow large temperature changes 
           if (tnew > 2.0 * tprev) then
               tnew = 2.0 * tprev
           endif
           if (tnew < 0.5 * tprev) then
               tnew = 0.5 * tprev
           endif

           ! exit if gas is too cold or too hot
           ! temperature and density limits are given in section 2.3 of Timmes & Swesty (2000)
           if (tnew > 1.0e11) then
               tnew = 1.0e11
               errorp = 1.0
           endif
           if (tnew < 1.0e4) then
               tnew = 1.0e4
               errorp = 1.0
           endif
       enddo iterations

       rel = abs(tnew - tempi) / tempi

       ! store new temperature, or energy if iterations fail
       ! TODO: currently we just use the final temperature from the eos and assume we have converged
!       if ((itercount > maxiter) .or. (rel > 0.05)) then
!           eni = cgseni_eos / unit_ergg
!       else
           tempi = tnew
!       endif

       if (itercount > maxiter) then
           print *, 'Helmholtz eos fail to converge'
       endif
   else
       print *, 'error in relaxflag in Helmholtz equation of state'
   endif
 
   ! convert cgs values to code units and return these values 
   ponrhoi  = cgspresi / (unit_pressure * rhoi)
   spsoundi = cgsspsoundi / unit_velocity
 
 end subroutine eos_helmholtz_pres_sound


! psi0 and its derivatives
real function psi0(z)
  real, intent(in) :: z
  psi0   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
end function psi0

real function dpsi0(z)
  real, intent(in) :: z
  dpsi0  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
end function dpsi0

real function ddpsi0(z)
  real, intent(in) :: z
  ddpsi0 = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)
end function ddpsi0


! psi1 and its derivatives
real function psi1(z)
  real, intent(in) :: z
  psi1   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
end function psi1

real function dpsi1(z)
  real, intent(in) :: z
  dpsi1  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
end function dpsi1

real function ddpsi1(z)
  real, intent(in) :: z
  ddpsi1 = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
end function ddpsi1


! psi2  and its derivatives
real function psi2(z)
  real, intent(in) :: z
  psi2   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
end function psi2

real function dpsi2(z)
  real, intent(in) :: z
  dpsi2  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
end function dpsi2

real function ddpsi2(z)
  real, intent(in) :: z
  ddpsi2 = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)
end function ddpsi2


! cubic hermite polynomial statement functions
! psi0 & derivatives
real function xpsi0(z)
  real, intent(in) :: z
  xpsi0  = z * z * (2.0d0*z - 3.0d0) + 1.0
end function xpsi0

real function xdpsi0(z)
  real, intent(in) :: z
  xdpsi0 = z * (6.0d0*z - 6.0d0)
end function xdpsi0


! psi1 & derivatives
real function xpsi1(z)
  real, intent(in) :: z
  xpsi1  = z * ( z * (z - 2.0d0) + 1.0d0)
end function xpsi1

real function xdpsi1(z)
  real, intent(in) :: z
  xdpsi1 = z * (3.0d0*z - 4.0d0) + 1.0d0
end function xdpsi1


! bicubic hermite polynomial statement function
real function h3(i,j,fi,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md)
  integer, intent(in) :: i, j
  real, intent(in)    :: fi(:)
  real, intent(in)    :: w0t, w1t, w0mt, w1mt, w0d, w1d, w0md, w1md

  h3 =   fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
       + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
       + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
       + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
       + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
       + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
       + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
       + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt

end function h3


! biquintic hermite polynomial statement function
real function h5(i,j,fi,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)
   integer, intent(in) :: i, j
   real, intent(in)    :: fi(:)
   real, intent(in)    :: w0t, w1t, w2t
   real, intent(in)    :: w0mt, w1mt, w2mt
   real, intent(in)    :: w0d, w1d, w2d
   real, intent(in)    :: w0md, w1md, w2md

   h5 =  fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
       + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
       + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
       + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
       + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
       + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
       + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
       + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
       + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
       + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
       + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
       + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
       + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
       + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
       + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
       + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
       + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
       + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt
 
end function h5



! given a temperature temp [K], density den [g/cm**3], and a composition
! characterized by abar and zbar, this routine returns most of the other
! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
! their derivatives with respect to temperature, density, abar, and zbar.
! other quantites such the normalized chemical potential eta (plus its
! derivatives), number density of electrons and positron pair (along
! with their derivatives), adiabatic indices, specific heats, and
! relativistically correct sound speed are also returned.
!
! this routine assumes planckian photons, an ideal gas of ions,
! and an electron-positron gas with an arbitrary degree of relativity
! and degeneracy. interpolation in a table of the helmholtz free energy
! is used to return the electron-positron thermodynamic quantities.
! all other derivatives are analytic.
!
! references: cox & giuli chapter 24 ; timmes & swesty apj 1999
 
 subroutine eos_helmholtz_compute_pres_sound(temp,den,pres,sound,ener,denerdt)
  use physcon, only:mass_proton_cgs,kboltz,c,planckh,steboltz,qe,avogadro,pi,fourpi,atomic_mass_unit

  real, intent(in)  :: temp, den
  real, intent(out) :: pres, sound, ener, denerdt
  
  real, parameter :: sioncon = (2.0d0 * pi * atomic_mass_unit * kboltz)/(planckh * planckh)
  real, parameter :: forth   = 4.0d0/3.0d0
  real, parameter :: kergavo = kboltz * avogadro
  real, parameter :: asol    = steboltz*4.0/c    
  real, parameter :: asoli3  = asol/3.0d0
  real, parameter :: light2  = c * c

  integer :: i, j
  real    :: ytot1,ye, &
                   x,y,zz,zzi,deni,tempi,xni,dxnidd, &
                   dpepdt,dpepdd,deepdt,dsepdt, &
                   dpraddd,dpraddt,deraddt,dpiondd,dpiondt, &
                   deiondt, &
                   kt,ktinv,prad,erad,srad,pion,eion, &
                   sion,pele,eele,sele,dpresdd, &
                   dpresdt,cv, &
                   gam1,chit,chid, &
                   s

  real :: pgas,dpgasdd,dpgasdt, &
                   egas,degasdt

! for the interpolations
  integer :: iat,jat
  real :: free,df_d,df_t,df_dd,df_tt,df_dt
  real :: xt,xd,mxt,mxd, &
                   si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                   si0d,si1d,si2d,si0md,si1md,si2md, &
                   dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                   dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                   ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                   ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
                   z, din
  real :: fi(36)

! for the uniform background coulomb correction
  real :: dsdd,lami,inv_lami,lamidd, &
                   plasg,plasgdd,plasgdt, &
                   ecoul,decouldd,decouldt, &
                   pcoul,dpcouldd,dpcouldt
  real, parameter :: a1    = -0.898004d0
  real, parameter :: b1    =  0.96786d0
  real, parameter :: c1    =  0.220703d0
  real, parameter :: d1    = -0.86097d0
  real, parameter :: e1    =  2.5269d0
  real, parameter :: a2    =  0.29561d0
  real, parameter :: b2    =  1.9885d0
  real, parameter :: c2    =  0.288675d0
  real, parameter :: third =  1.0d0/3.0d0
  real, parameter :: esqu  =  qe * qe


! start of pipeline loop, normal execution starts here

    if (temp < 1.0e4) then
        print *, 'WARNING: ', temp, ' eos_helmholtz: temperature below range available (10^4 K - 10^11 K)'
    endif
    if (temp > 1.0e11) then
        print *, 'WARNING: eos_helmholtz: temperature exceeds range available (10^4 - 10^11 K)'
    endif

    if (den < 1.0e-6) then
        print *, 'WARNING: eos_helmholtz: density below range available (10-^6 - 10^11 g/cm^3)'
    endif
    if (den > 1.0e11) then
        print *, 'WARNING: eos_helmholtz: density exceeds range available (10-^6 - 10^11 g/cm^3)'
    endif

    if (temp  <=  0.0) stop 'temp less than 0 in helmeos'
    if (den   <=  0.0) then
        print*, 'den less than 0 in helmeos'
        stop
    end if

    ytot1 = 1.0d0/abar
    ye    = max(1.0d-16,ytot1 * zbar)


! initialize
   deni    = 1.0d0/den
   tempi   = 1.0d0/temp
   kt      = kboltz * temp
   ktinv   = 1.0d0/kt


! radiation section:
   prad    = asoli3 * temp * temp * temp * temp
   dpraddd = 0.0d0
   dpraddt = 4.0d0 * prad*tempi

   erad    = 3.0d0 * prad*deni

   deraddt = 3.0d0 * dpraddt*deni

   srad    = (prad*deni + erad)*tempi


! ion section:
    xni     = avogadro * ytot1 * den

    dxnidd  = avogadro * ytot1

    pion    = xni * kt
    dpiondd = dxnidd * kt
    dpiondt = xni * kboltz

    eion    = 1.5d0 * pion*deni

    deiondt = 1.5d0 * dpiondt*deni


! sackur-tetrode equation for the ion entropy of
! a single ideal gas characterized by abar
    x       = abar*abar*sqrt(abar) * deni/avogadro
    s       = sioncon * temp
    z       = x * s * sqrt(s)
    y       = log(z)

!        y       = 1.0d0/(abar*kt)
!        yy      = y * sqrt(y)
!        z       = xni * sifac * yy
!        etaion  = log(z)

    sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y



! electron-positron section:

! assume complete ionization

! enter the table with ye*den
    din = ye*den


! bomb proof the input

! hash locate this temperature and density
    jat = int((log10(temp) - tlo)*tstpi) + 1
    jat = max(1,min(jat,jmax-1))
    iat = int((log10(din) - dlo)*dstpi) + 1
    iat = max(1,min(iat,imax-1))


! access the table locations only once
    fi(1)  = f(iat,jat)
    fi(2)  = f(iat+1,jat)
    fi(3)  = f(iat,jat+1)
    fi(4)  = f(iat+1,jat+1)
    fi(5)  = ft(iat,jat)
    fi(6)  = ft(iat+1,jat)
    fi(7)  = ft(iat,jat+1)
    fi(8)  = ft(iat+1,jat+1)
    fi(9)  = ftt(iat,jat)
    fi(10) = ftt(iat+1,jat)
    fi(11) = ftt(iat,jat+1)
    fi(12) = ftt(iat+1,jat+1)
    fi(13) = fd(iat,jat)
    fi(14) = fd(iat+1,jat)
    fi(15) = fd(iat,jat+1)
    fi(16) = fd(iat+1,jat+1)
    fi(17) = fdd(iat,jat)
    fi(18) = fdd(iat+1,jat)
    fi(19) = fdd(iat,jat+1)
    fi(20) = fdd(iat+1,jat+1)
    fi(21) = fdt(iat,jat)
    fi(22) = fdt(iat+1,jat)
    fi(23) = fdt(iat,jat+1)
    fi(24) = fdt(iat+1,jat+1)
    fi(25) = fddt(iat,jat)
    fi(26) = fddt(iat+1,jat)
    fi(27) = fddt(iat,jat+1)
    fi(28) = fddt(iat+1,jat+1)
    fi(29) = fdtt(iat,jat)
    fi(30) = fdtt(iat+1,jat)
    fi(31) = fdtt(iat,jat+1)
    fi(32) = fdtt(iat+1,jat+1)
    fi(33) = fddtt(iat,jat)
    fi(34) = fddtt(iat+1,jat)
    fi(35) = fddtt(iat,jat+1)
    fi(36) = fddtt(iat+1,jat+1)


! various differences
    xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
    xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
    mxt = 1.0d0 - xt
    mxd = 1.0d0 - xd

! the six density and six temperature basis functions
    si0t =   psi0(xt)
    si1t =   psi1(xt)*dt_sav(jat)
    si2t =   psi2(xt)*dt2_sav(jat)

    si0mt =  psi0(mxt)
    si1mt = -psi1(mxt)*dt_sav(jat)
    si2mt =  psi2(mxt)*dt2_sav(jat)

    si0d =   psi0(xd)
    si1d =   psi1(xd)*dd_sav(iat)
    si2d =   psi2(xd)*dd2_sav(iat)

    si0md =  psi0(mxd)
    si1md = -psi1(mxd)*dd_sav(iat)
    si2md =  psi2(mxd)*dd2_sav(iat)

! derivatives of the weight functions
    dsi0t =   dpsi0(xt)*dti_sav(jat)
    dsi1t =   dpsi1(xt)
    dsi2t =   dpsi2(xt)*dt_sav(jat)

    dsi0mt = -dpsi0(mxt)*dti_sav(jat)
    dsi1mt =  dpsi1(mxt)
    dsi2mt = -dpsi2(mxt)*dt_sav(jat)

    dsi0d =   dpsi0(xd)*ddi_sav(iat)
    dsi1d =   dpsi1(xd)
    dsi2d =   dpsi2(xd)*dd_sav(iat)

    dsi0md = -dpsi0(mxd)*ddi_sav(iat)
    dsi1md =  dpsi1(mxd)
    dsi2md = -dpsi2(mxd)*dd_sav(iat)

! second derivatives of the weight functions
    ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
    ddsi1t =   ddpsi1(xt)*dti_sav(jat)
    ddsi2t =   ddpsi2(xt)

    ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
    ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
    ddsi2mt =  ddpsi2(mxt)

!        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
!        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
!        ddsi2d =   ddpsi2(xd)

!        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
!        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
!        ddsi2md =  ddpsi2(mxd)


! the free energy
    free  = h5(iat,jat,fi, &
            si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density
    df_d  = h5(iat,jat,fi, &
            si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
            dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


! derivative with respect to temperature
    df_t = h5(iat,jat,fi, &
            dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density**2
!        df_dd = h5(iat,jat,
!     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
!     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

! derivative with respect to temperature**2

    df_tt = h5(iat,jat,fi, &
          ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to temperature and density
    df_dt = h5(iat,jat,fi, &
            dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
            dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



! now get the pressure derivative with density, chemical potential, and
! electron positron number densities
! get the interpolation weight functions
    si0t   =  xpsi0(xt)
    si1t   =  xpsi1(xt)*dt_sav(jat)

    si0mt  =  xpsi0(mxt)
    si1mt  =  -xpsi1(mxt)*dt_sav(jat)

    si0d   =  xpsi0(xd)
    si1d   =  xpsi1(xd)*dd_sav(iat)

    si0md  =  xpsi0(mxd)
    si1md  =  -xpsi1(mxd)*dd_sav(iat)


! derivatives of weight functions
    dsi0t  = xdpsi0(xt)*dti_sav(jat)
    dsi1t  = xdpsi1(xt)

    dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
    dsi1mt = xdpsi1(mxt)

    dsi0d  = xdpsi0(xd)*ddi_sav(iat)
    dsi1d  = xdpsi1(xd)

    dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
    dsi1md = xdpsi1(mxd)


! look in the pressure derivative only once
    fi(1)  = dpdf(iat,jat)
    fi(2)  = dpdf(iat+1,jat)
    fi(3)  = dpdf(iat,jat+1)
    fi(4)  = dpdf(iat+1,jat+1)
    fi(5)  = dpdft(iat,jat)
    fi(6)  = dpdft(iat+1,jat)
    fi(7)  = dpdft(iat,jat+1)
    fi(8)  = dpdft(iat+1,jat+1)
    fi(9)  = dpdfd(iat,jat)
    fi(10) = dpdfd(iat+1,jat)
    fi(11) = dpdfd(iat,jat+1)
    fi(12) = dpdfd(iat+1,jat+1)
    fi(13) = dpdfdt(iat,jat)
    fi(14) = dpdfdt(iat+1,jat)
    fi(15) = dpdfdt(iat,jat+1)
    fi(16) = dpdfdt(iat+1,jat+1)

! pressure derivative with density
    dpepdd  = h3(iat,jat,fi, &
                   si0t,   si1t,   si0mt,   si1mt, &
                   si0d,   si1d,   si0md,   si1md)
    dpepdd  = max(ye * dpepdd,1.0d-30)



! look in the electron chemical potential table only once
    fi(1)  = ef(iat,jat)
    fi(2)  = ef(iat+1,jat)
    fi(3)  = ef(iat,jat+1)
    fi(4)  = ef(iat+1,jat+1)
    fi(5)  = eft(iat,jat)
    fi(6)  = eft(iat+1,jat)
    fi(7)  = eft(iat,jat+1)
    fi(8)  = eft(iat+1,jat+1)
    fi(9)  = efd(iat,jat)
    fi(10) = efd(iat+1,jat)
    fi(11) = efd(iat,jat+1)
    fi(12) = efd(iat+1,jat+1)
    fi(13) = efdt(iat,jat)
    fi(14) = efdt(iat+1,jat)
    fi(15) = efdt(iat,jat+1)
    fi(16) = efdt(iat+1,jat+1)


! electron chemical potential etaele


! derivative with respect to density



! derivative with respect to temperature



! look in the number density table only once
    fi(1)  = xf(iat,jat)
    fi(2)  = xf(iat+1,jat)
    fi(3)  = xf(iat,jat+1)
    fi(4)  = xf(iat+1,jat+1)
    fi(5)  = xft(iat,jat)
    fi(6)  = xft(iat+1,jat)
    fi(7)  = xft(iat,jat+1)
    fi(8)  = xft(iat+1,jat+1)
    fi(9)  = xfd(iat,jat)
    fi(10) = xfd(iat+1,jat)
    fi(11) = xfd(iat,jat+1)
    fi(12) = xfd(iat+1,jat+1)
    fi(13) = xfdt(iat,jat)
    fi(14) = xfdt(iat+1,jat)
    fi(15) = xfdt(iat,jat+1)
    fi(16) = xfdt(iat+1,jat+1)

! electron + positron number densities


! derivative with respect to density


! derivative with respect to temperature


! derivative with respect to abar and zbar


! the desired electron-positron thermodynamic quantities

! dpepdd at high temperatures and low densities is below the
! floating point limit of the subtraction of two large terms.
! since dpresdd doesn't enter the maxwell relations at all, use the
! bicubic interpolation done above instead of the formally correct expression
    x       = din * din
    pele    = x * df_d
    dpepdt  = x * df_dt
!        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)


    x       = ye * ye
    sele    = -df_t * ye

    dsepdt  = -df_tt * ye
    eele    = ye*free + temp * sele

    deepdt  = temp * dsepdt


! coulomb section:

! uniform background corrections only
! from yakovlev & shalybkov 1989
! lami is the average ion seperation
! plasg is the plasma coupling parameter

    z        = forth * pi
    s        = z * xni
    dsdd     = z * dxnidd


    lami     = 1.0d0/s**third
    inv_lami = 1.0d0/lami
    z        = -third * lami
    lamidd   = z * dsdd/s

    plasg    = zbar*zbar*esqu*ktinv*inv_lami
    z        = -plasg * inv_lami
    plasgdd  = z * lamidd
    plasgdt  = -plasg*ktinv * kboltz


! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    if (plasg  >=  1.0) then
     x        = plasg**(0.25d0)
     y        = avogadro * ytot1 * kboltz
     ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
     pcoul    = third * den * ecoul


     y        = avogadro*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
     decouldd = y * plasgdd
     decouldt = y * plasgdt + ecoul/temp


     y        = third * den
     dpcouldd = third * ecoul + y*decouldd
     dpcouldt = y * decouldt



! yakovlev & shalybkov 1989 equations 102, 103, 104
    else if (plasg  <  1.0) then
     x        = plasg*sqrt(plasg)
     y        = plasg**b2
     z        = c2 * x - third * a2 * y
     pcoul    = -pion * z
     ecoul    = 3.0d0 * pcoul/den


     s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
     dpcouldd = -dpiondd*z - pion*s*plasgdd
     dpcouldt = -dpiondt*z - pion*s*plasgdt


     s        = 3.0d0/den
     decouldt = s * dpcouldt

    end if


! bomb proof
    x   = prad + pion + pele + pcoul
    y   = erad + eion + eele + ecoul
    z   = srad + sion + sele 

!        write(6,*) x,y,z
    if (x  <=  0.0 .or. y  <=  0.0 .or. z  <=  0.0) then
!        if (x  <=  0.0 .or. y  <=  0.0) then
!        if (x  <=  0.0) then

!         write(6,*)
!         write(6,*) 'coulomb corrections are causing a negative pressure'
!         write(6,*) 'setting all coulomb corrections to zero'
!         write(6,*)

     pcoul    = 0.0d0
     dpcouldd = 0.0d0
     dpcouldt = 0.0d0
     ecoul    = 0.0d0
     decouldd = 0.0d0
     decouldt = 0.0d0
    end if


! sum all the gas components
   pgas    = pion + pele + pcoul

   egas    = eion + eele + ecoul

   !print *, dpiondt,dpepdt,dpcouldt
   dpgasdd = dpiondd + dpepdd + dpcouldd
   dpgasdt = dpiondt + dpepdt + dpcouldt


   degasdt = deiondt + deepdt + decouldt



! add in radiation to get the total
   pres    = prad + pgas
   ener    = erad + egas


   dpresdd = dpraddd + dpgasdd
   
   dpresdt = dpraddt + dpgasdt
!    print *, deraddt, degasdt
   denerdt = deraddt + degasdt

! for the totals
   zz    = pres*deni
   zzi   = den/pres

   chit  = temp/pres * dpresdt
   chid  = dpresdd*zzi
   cv    = denerdt
   x     = zz * chit/(temp * cv)

   gam1  = chit*x + chid


   z     = 1.0d0 + (ener + light2)*zzi

   sound = c * sqrt(gam1/z)
   


! end of pipeline loop
 end subroutine eos_helmholtz_compute_pres_sound



!----------------------------------------------------------------
!+
!  Broken off from the main helmeos routine to return just c_v and dpres/dt
!+
!----------------------------------------------------------------
subroutine eos_helmholtz_cv_dpresdt(temp,den,cv,dpresdt)
  use physcon, only:mass_proton_cgs,kboltz,c,planckh,steboltz,qe,avogadro,pi,fourpi,atomic_mass_unit
  real, intent(in)  :: temp, den
  real, intent(out) :: cv, dpresdt
  integer :: i,j
  real    :: ytot1,ye, &
             x,y,zz,zzi,deni,tempi,xni, &
             dpepdt,deepdt,dsepdt, &
             dpraddt,deraddt,dpiondt, &
             deiondt, &
             kt,ktinv,prad,erad,srad,pion,eion, &
             pele,eele,sele, &
             denerdt,s 
  real    :: dpgasdt, degasdt

  real, parameter :: sioncon = (2.0d0 * pi * atomic_mass_unit * kboltz)/(planckh * planckh)
  real, parameter :: forth   = 4.0d0/3.0d0
  real, parameter :: asol    = steboltz*4.0/c    
  real, parameter :: asoli3  = asol/3.0d0
  real, parameter :: light2  = c * c

  integer :: iat, jat
  real    :: free,df_d,df_t,df_dd,df_tt,df_dt
  real    :: xt,xd,mxt,mxd, &
                     si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                     si0d,si1d,si2d,si0md,si1md,si2md, &
                     dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                     dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                     ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                     ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
                     z,din
  real    :: fi(36)

    ! for the uniform background coulomb correction
  real    :: lami,inv_lami, &
             plasg,plasgdt, &
             ecoul,decouldt, &
             pcoul,dpcouldt

  real, parameter :: a1    = -0.898004d0
  real, parameter :: b1    =  0.96786d0
  real, parameter :: c1    =  0.220703d0
  real, parameter :: d1    = -0.86097d0
  real, parameter :: e1    =  2.5269d0
  real, parameter :: a2    =  0.29561d0
  real, parameter :: b2    =  1.9885d0
  real, parameter :: c2    =  0.288675d0
  real, parameter :: third =  1.0d0/3.0d0
  real, parameter :: esqu  =  qe * qe


    ! start of pipeline loop, normal execution starts here

    if (temp  <=  0.0) stop 'temp less than 0 in helmeos'
    if (den   <=  0.0) then
        print*, 'den less than 0 in helmeos'
        stop
    end if

    ytot1 = 1.0d0/abar
    ye    = max(1.0d-16,ytot1 * zbar)

    ! initialize
    deni    = 1.0d0/den
    tempi   = 1.0d0/temp

    kt      = kboltz * temp
    ktinv   = 1.0d0/kt

    ! radiation section:
    prad    = asoli3 * temp * temp * temp * temp
    dpraddt = 4.0d0 * prad*tempi

    erad    = 3.0d0 * prad*deni
    deraddt = 3.0d0 * dpraddt*deni

    srad    = (prad*deni + erad)*tempi

    ! ion section:
    xni     = avogadro * ytot1 * den
    pion    = xni * kt
    dpiondt = xni * kboltz
    deiondt = 1.5d0 * dpiondt*deni

    eion    = 1.5d0 * pion*deni


    ! assume complete ionization

    ! enter the table with ye*den
    din = ye*den
    ! bomb proof the input

    ! hash locate this temperature and density
    jat = int((log10(temp) - tlo)*tstpi) + 1
    jat = max(1,min(jat,jmax-1))
    iat = int((log10(din) - dlo)*dstpi) + 1
    iat = max(1,min(iat,imax-1))


! access the table locations only once
    fi(1)  = f(iat,jat)
    fi(2)  = f(iat+1,jat)
    fi(3)  = f(iat,jat+1)
    fi(4)  = f(iat+1,jat+1)
    fi(5)  = ft(iat,jat)
    fi(6)  = ft(iat+1,jat)
    fi(7)  = ft(iat,jat+1)
    fi(8)  = ft(iat+1,jat+1)
    fi(9)  = ftt(iat,jat)
    fi(10) = ftt(iat+1,jat)
    fi(11) = ftt(iat,jat+1)
    fi(12) = ftt(iat+1,jat+1)
    fi(13) = fd(iat,jat)
    fi(14) = fd(iat+1,jat)
    fi(15) = fd(iat,jat+1)
    fi(16) = fd(iat+1,jat+1)
    fi(17) = fdd(iat,jat)
    fi(18) = fdd(iat+1,jat)
    fi(19) = fdd(iat,jat+1)
    fi(20) = fdd(iat+1,jat+1)
    fi(21) = fdt(iat,jat)
    fi(22) = fdt(iat+1,jat)
    fi(23) = fdt(iat,jat+1)
    fi(24) = fdt(iat+1,jat+1)
    fi(25) = fddt(iat,jat)
    fi(26) = fddt(iat+1,jat)
    fi(27) = fddt(iat,jat+1)
    fi(28) = fddt(iat+1,jat+1)
    fi(29) = fdtt(iat,jat)
    fi(30) = fdtt(iat+1,jat)
    fi(31) = fdtt(iat,jat+1)
    fi(32) = fdtt(iat+1,jat+1)
    fi(33) = fddtt(iat,jat)
    fi(34) = fddtt(iat+1,jat)
    fi(35) = fddtt(iat,jat+1)
    fi(36) = fddtt(iat+1,jat+1)


! various differences
    xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
    xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
    mxt = 1.0d0 - xt
    mxd = 1.0d0 - xd

! the six density and six temperature basis functions
    si0t =   psi0(xt)
    si1t =   psi1(xt)*dt_sav(jat)
    si2t =   psi2(xt)*dt2_sav(jat)

    si0mt =  psi0(mxt)
    si1mt = -psi1(mxt)*dt_sav(jat)
    si2mt =  psi2(mxt)*dt2_sav(jat)

    si0d =   psi0(xd)
    si1d =   psi1(xd)*dd_sav(iat)
    si2d =   psi2(xd)*dd2_sav(iat)

    si0md =  psi0(mxd)
    si1md = -psi1(mxd)*dd_sav(iat)
    si2md =  psi2(mxd)*dd2_sav(iat)

! derivatives of the weight functions
    dsi0t =   dpsi0(xt)*dti_sav(jat)
    dsi1t =   dpsi1(xt)
    dsi2t =   dpsi2(xt)*dt_sav(jat)

    dsi0mt = -dpsi0(mxt)*dti_sav(jat)
    dsi1mt =  dpsi1(mxt)
    dsi2mt = -dpsi2(mxt)*dt_sav(jat)

    dsi0d =   dpsi0(xd)*ddi_sav(iat)
    dsi1d =   dpsi1(xd)
    dsi2d =   dpsi2(xd)*dd_sav(iat)

    dsi0md = -dpsi0(mxd)*ddi_sav(iat)
    dsi1md =  dpsi1(mxd)
    dsi2md = -dpsi2(mxd)*dd_sav(iat)

! second derivatives of the weight functions
    ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
    ddsi1t =   ddpsi1(xt)*dti_sav(jat)
    ddsi2t =   ddpsi2(xt)

    ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
    ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
    ddsi2mt =  ddpsi2(mxt)

!        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
!        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
!        ddsi2d =   ddpsi2(xd)

!        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
!        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
!        ddsi2md =  ddpsi2(mxd)


! the free energy
    free  = h5(iat,jat,fi, &
            si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density
    df_d  = h5(iat,jat,fi, &
            si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
            dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


! derivative with respect to temperature
    df_t = h5(iat,jat,fi, &
            dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
    ! derivative with respect to temperature**2

    df_tt = h5(iat,jat,fi, &
          ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
    ! derivative with respect to temperature and density
    df_dt = h5(iat,jat,fi, &
            dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
            dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

    ! now get the pressure derivative with density, chemical potential, and
    ! electron positron number densities
    ! get the interpolation weight functions

    ! derivatives of weight functions


    ! look in the pressure derivative only once
    fi(1)  = dpdf(iat,jat)
    fi(2)  = dpdf(iat+1,jat)
    fi(3)  = dpdf(iat,jat+1)
    fi(4)  = dpdf(iat+1,jat+1)
    fi(5)  = dpdft(iat,jat)
    fi(6)  = dpdft(iat+1,jat)
    fi(7)  = dpdft(iat,jat+1)
    fi(8)  = dpdft(iat+1,jat+1)
    fi(9)  = dpdfd(iat,jat)
    fi(10) = dpdfd(iat+1,jat)
    fi(11) = dpdfd(iat,jat+1)
    fi(12) = dpdfd(iat+1,jat+1)
    fi(13) = dpdfdt(iat,jat)
    fi(14) = dpdfdt(iat+1,jat)
    fi(15) = dpdfdt(iat,jat+1)
    fi(16) = dpdfdt(iat+1,jat+1)
    ! pressure derivative with density


    ! look in the electron chemical potential table only once
    fi(1)  = ef(iat,jat)
    fi(2)  = ef(iat+1,jat)
    fi(3)  = ef(iat,jat+1)
    fi(4)  = ef(iat+1,jat+1)
    fi(5)  = eft(iat,jat)
    fi(6)  = eft(iat+1,jat)
    fi(7)  = eft(iat,jat+1)
    fi(8)  = eft(iat+1,jat+1)
    fi(9)  = efd(iat,jat)
    fi(10) = efd(iat+1,jat)
    fi(11) = efd(iat,jat+1)
    fi(12) = efd(iat+1,jat+1)
    fi(13) = efdt(iat,jat)
    fi(14) = efdt(iat+1,jat)
    fi(15) = efdt(iat,jat+1)
    fi(16) = efdt(iat+1,jat+1)
    ! electron chemical potential etaele


    ! look in the number density table only once
    fi(1)  = xf(iat,jat)
    fi(2)  = xf(iat+1,jat)
    fi(3)  = xf(iat,jat+1)
    fi(4)  = xf(iat+1,jat+1)
    fi(5)  = xft(iat,jat)
    fi(6)  = xft(iat+1,jat)
    fi(7)  = xft(iat,jat+1)
    fi(8)  = xft(iat+1,jat+1)
    fi(9)  = xfd(iat,jat)
    fi(10) = xfd(iat+1,jat)
    fi(11) = xfd(iat,jat+1)
    fi(12) = xfd(iat+1,jat+1)
    fi(13) = xfdt(iat,jat)
    fi(14) = xfdt(iat+1,jat)
    fi(15) = xfdt(iat,jat+1)
    fi(16) = xfdt(iat+1,jat+1)


    ! electron + positron number densities

    ! derivative with respect to density

    ! derivative with respect to temperature


    ! derivative with respect to abar and zbar


    ! the desired electron-positron thermodynamic quantities

    ! dpepdd at high temperatures and low densities is below the
    ! floating point limit of the subtraction of two large terms.
    ! since dpresdd doesn't enter the maxwell relations at all, use the
    ! bicubic interpolation done above instead of the formally correct expression
    x       = din * din
    pele    = x * df_d
    dpepdt  = x * df_dt

    x       = ye * ye
    sele    = -df_t * ye

    dsepdt  = -df_tt * ye

    eele    = ye*free + temp * sele

    deepdt  = temp * dsepdt

    ! coulomb section:
    z        = forth * pi
    s        = z * xni

    lami     = 1.0d0/s**third
    inv_lami = 1.0d0/lami

    plasg    = zbar*zbar*esqu*ktinv*inv_lami
    plasgdt  = -plasg*ktinv * kboltz


    ! uniform background corrections only
    ! from yakovlev & shalybkov 1989
    ! lami is the average ion seperation
    ! plasg is the plasma coupling parameter


    ! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    if (plasg  >=  1.0) then
        x        = plasg**(0.25d0)
        y        = avogadro * ytot1 * kboltz
        ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)


        y        = avogadro*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
        decouldt = y * plasgdt + ecoul/temp

        y        = third * den
        dpcouldt = y * decouldt


    ! yakovlev & shalybkov 1989 equations 102, 103, 104
    else if (plasg  <  1.0) then
        x        = plasg*sqrt(plasg)

        y        = plasg**b2
        z        = c2 * x - third * a2 * y
        pcoul    = -pion * z
        ecoul    = 3.0d0 * pcoul/den
        s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
        dpcouldt = -dpiondt*z - pion*s*plasgdt

        s        = 3.0d0/den
        decouldt = s * dpcouldt

    end if

   
    ! bomb proof
    x   = prad + pion + pele + pcoul
    y   = erad + eion + eele + ecoul
    z   = srad + sele 

!        write(6,*) x,y,z
    if (x  <=  0.0 .or. y  <=  0.0 .or. z  <=  0.0) then
!        if (x  <=  0.0 .or. y  <=  0.0) then
!        if (x  <=  0.0) then

!         write(6,*)
!         write(6,*) 'coulomb corrections are causing a negative pressure'
!         write(6,*) 'setting all coulomb corrections to zero'
!         write(6,*)
     dpcouldt = 0.0d0
     ecoul    = 0.0d0
     decouldt = 0.0d0
    end if
  
    ! sum all the gas components
    dpgasdt = dpiondt + dpepdt + dpcouldt
    degasdt = deiondt + deepdt + decouldt

    ! add in radiation to get the total
    dpresdt = dpraddt + dpgasdt
    denerdt = deraddt + degasdt

    ! for the totals
    cv    = denerdt

end subroutine eos_helmholtz_cv_dpresdt

end module eos_helmholtz
