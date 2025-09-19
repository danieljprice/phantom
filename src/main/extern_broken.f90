!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_broken
!
! This module contains routines relating to the computation
!    of a gravitational force/potential that breaks at a radial location
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - eps_soft1  : *Plummer softening of primary*
!   - m1         : *m1 inside of rbreak*
!   - m2         : *m2 outside of rbreak*
!   - rbreak     : *break radius*
!   - delta_r    : *radial distance to smooth between regions*
!
! :Dependencies: dump_utils, infile_utils, io, physcon
!
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !

 real, public :: bmass1 = 1.0
 real, public :: bmass2 = 2.0
 real, public :: baccradius1  = 1.0
 real, public :: eps_soft1 = 0.0
 real, public :: rbreak = 3.
 real, public :: delta_r = 0.2
 real, public :: eps2_soft = 0.0
 real, public :: omega1 = 1.0
 real, public :: nratio = 1.5
 logical, public :: surface_force = .false.

 public :: brokenstar_potential, broken_omegapotential, mass_eff
 public :: write_options_extern_brokenstar, read_options_extern_brokenstar
 private

 real, private :: x1,y1,x2,y2

 ! Required for HDF5 compatibility
 real, public :: a0 = 0.
 real, public :: direction = 0.

contains


subroutine brokenstar_potential(xi,yi,zi,fextxi,fextyi,fextzi,phi)
 real, intent(in)  :: xi,yi,zi
 real, intent(out) :: fextxi,fextyi,fextzi,phi
 real             :: r2, r, rinv, rinv3, meff
 real             :: rlow, rhigh, h, a, Cshift, half, q
 real             :: rsoft2, rinv_soft, rinv_soft3, rinv_soft_high

fextxi = 0.0
fextyi = 0.0
fextzi = 0.0
phi    = 0.0

half  = 0.5*delta_r
rlow  = rbreak - half
rhigh = rbreak + half
q     = bmass2 - bmass1      

r2       = xi*xi + yi*yi + zi*zi
r        = sqrt(r2)                         
rsoft2   = r2 + eps2_soft
if (rsoft2 <= epsilon(rsoft2)) return        ! particle at the origin

rinv_soft  = 1.0 / sqrt(rsoft2)             
rinv_soft3 = rinv_soft*rinv_soft*rinv_soft   

! ---------- enclosed mass & force ----------
meff      = mass_eff(r, bmass1, bmass2, rlow, delta_r)  
rinv_soft3 = meff*rinv_soft3      
fextxi = fextxi - xi*rinv_soft3
fextyi = fextyi - yi*rinv_soft3
fextzi = fextzi - zi*rinv_soft3

rinv_soft_high = 1.0 / sqrt(rhigh*rhigh + eps2_soft)
Cshift         = -(bmass1 - bmass2)*rinv_soft_high

if (r <= rlow) then                     ! inner zone
    phi = phi - bmass1*rinv_soft

elseif (r >= rhigh) then                ! outer zone
    phi = phi - bmass2*rinv_soft + Cshift

else                                    ! inside the linear break layer
    phi = - ( bmass2/rhigh                                   &
              + bmass1*(rinv_soft - 1.0/rhigh)               &
              + q/delta_r*( log(rhigh/r) + rlow*(1.0/rhigh - rinv_soft) ) ) &
          + Cshift
endif


end subroutine brokenstar_potential


! ---------------------------------------------------------------------
! + 
! A potential where Omega is constant between two regions
! + 
! ---------------------------------------------------------------------

subroutine broken_omegapotential(xi,yi,zi,fextxi,fextyi,fextzi,phi)
   real, intent(in)    :: xi,yi,zi
   real, intent(inout) :: fextxi,fextyi,fextzi, phi
   real :: r2, r, half, rlow, rhigh, s, d
   real :: omega_eff, force_coeff
   real :: a, b, Dc, Cc, phi_layer, omega2

   !========== geometry =========================================
   half  = 0.5*delta_r
   rlow  = rbreak - half
   rhigh = rbreak + half
   d     = delta_r
   omega2 = nratio*omega1
   s     = omega2 - omega1        ! Ω-step

   !========== radius ===========================================
   r2 = xi*xi + yi*yi + zi*zi
   if (r2 <= epsilon(r2)) return   ! at the centre => no force
   r  = sqrt(r2)

   !========== Ω(r) and force ===================================
   if (r <= rlow) then
      omega_eff = omega1
   elseif (r >= rhigh) then
      omega_eff = omega2
   else
      omega_eff = omega1 + s*(r - rlow)/d   ! linear ramp
   end if

   force_coeff = - omega_eff*omega_eff      ! -Ω²

   fextxi = fextxi + force_coeff*xi
   fextyi = fextyi + force_coeff*yi
   fextzi = fextzi + force_coeff*zi

   !========== potential ========================================
   if (r <= rlow) then
      phi = phi + 0.5*omega1*omega1*r*r
      return
   end if

   ! ---- constants needed for continuity (computed once) -------
   a = 2.0*omega1*s/d
   b = (s*s)/(d*d)
   Dc = -a*rlow**3/6.0 + b*rlow**4/12.0   ! shift at rlow

   ! φ inside the layer
   if (r < rhigh) then
      phi_layer = 0.5*omega1*omega1*r*r                              &
                + a*( r**3/3.0 - rlow*r*r/2.0 )                   &
                + b*( r**4/4.0 - 2.0*r**3*rlow/3.0             &
                       + r**2*rlow*rlow/2.0 )                        &
                - Dc
      phi = phi + phi_layer
      return
   end if

   ! ---- outer zone (compute C the first time we need it) ------
   Cc = ( 0.5*omega1*omega1*rhigh*rhigh                               &
         + a*( rhigh**3/3.0 - rlow*rhigh*rhigh/2.0 )              &
         + b*( rhigh**4/4.0 - 2.0*rhigh**3*rlow/3.0            &
                + rhigh**2*rlow*rlow/2.0 )                           &
         - Dc )                                                          &
       - 0.5*omega2*omega2*rhigh*rhigh

   phi = phi + 0.5*omega2*omega2*r*r + Cc

end subroutine broken_omegapotential



!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_extern_brokenstar(iunit)
 use infile_utils,   only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(bmass1,'mass1','mass of object inside rbreak',iunit)
 call write_inopt(bmass2,'mass2','mass of object outside rbreak',iunit)
 call write_inopt(rbreak,'rbreak','break radius',iunit)
 call write_inopt(delta_r,'delta_r','radial distance to smooth between regions',iunit)
 call write_inopt(eps_soft1,'eps_soft1','Plummer softening of primary',iunit)


end subroutine write_options_extern_brokenstar


subroutine write_options_extern_broken_solidbody(iunit)
 use infile_utils,   only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(omega1,'omega1','angular velocity of object inside rbreak',iunit)
 call write_inopt(nratio,'nratio','mass of object outside rbreak',iunit)
 call write_inopt(rbreak,'rbreak','break radius',iunit)
 call write_inopt(delta_r,'delta_r','radial distance to smooth between regions',iunit)

end subroutine write_options_extern_broken_solidbody

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_extern_brokenstar(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,error
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: where = 'read_options_externbinary'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('mass1')
    read(valstring,*,iostat=ierr) bmass1
    ngot = ngot + 1
    if (bmass1 < 0.) then
       call fatal(where,'invalid setting for mass1 (<0)')
    endif
 case('mass2')
    read(valstring,*,iostat=ierr) bmass2
    ngot = ngot + 1
    if (bmass2 < epsilon(bmass2)) &
       call fatal(where,'invalid setting for mass2 (<0)')
    if (bmass2/bmass1 > 1.e10)  call error(where,'binary mass ratio is huge!!!')
 case('rbreak')
    read(valstring,*,iostat=ierr) rbreak
    ngot = ngot + 1
 case('delta_r')
    read(valstring,*,iostat=ierr) delta_r
    ngot = ngot + 1
 case('eps_soft1')
    read(valstring,*,iostat=ierr) eps_soft1
    if (eps_soft1 < 0.)  call fatal(where,'negative eps_soft1')
    eps2_soft = eps_soft1*eps_soft1

 case default
    imatch = .false.
 end select

 igotall = (ngot >= 5)

end subroutine read_options_extern_brokenstar

subroutine read_options_extern_broken_omegapot(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,error
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: where = 'read_options_externbinary'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('omega1')
    read(valstring,*,iostat=ierr) omega1
    ngot = ngot + 1
    if (omega1 < 0.) then
       call fatal(where,'invalid setting for omega1 (<0)')
    endif
 case('nratio')
    read(valstring,*,iostat=ierr) nratio
    ngot = ngot + 1
    if (nratio < epsilon(nratio)) &
       call fatal(where,'invalid setting for nratio (<0)')
 case('rbreak')
    read(valstring,*,iostat=ierr) rbreak
    ngot = ngot + 1
 case('delta_r')
    read(valstring,*,iostat=ierr) delta_r
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 4)

end subroutine read_options_extern_broken_omegapot


pure function mass_eff(r, m1, m2, rlow, dR) result(m)
   real, intent(in) :: r,m1,m2,rlow,dR
   real             :: m, x
   if (r <= rlow) then
       m = m1
   elseif (r >= rlow + dR) then
       m = m2
   else
       x = (r-rlow)/dR       
       m = m1 + (m2-m1)*x     
   endif
end function mass_eff


end module extern_broken
