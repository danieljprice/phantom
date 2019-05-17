!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: extern_prdrag
!
!  DESCRIPTION:
! This module contains routines relating to the computation
! of radial and transverse radiation forces from a central sink particle
! centered at the origin.
!
! This module is intended to represent radiation drag as
! generically as possible. Drag is parametrized using "beta"
! which is defined as the ratio of radiation force to gravitational
! force. You need to specify your own routine for calculating
! beta, and include it in this module. This is the exhaustive list
! of lines you need to change:
!
! subroutine get_prdrag_spatial_force-- use beta_module, only:beta
! subroutine get_prdrag_vdependent_force-- use beta_module, only:beta
! subroutine update_prdrag_leapfrog-- use beta_module, only:beta
! subroutine write_options_prdrag-- use beta_module, only:write_options_beta
! subroutine read_options_prdrag-- use beta_module, only:read_options_beta
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, fastmath, infile_utils, io, lumin_nsdisc, physcon,
!    units
!+
!--------------------------------------------------------------------------
module extern_prdrag
 use eos, only:qfacdisc

 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 real, private    :: k2 = 1.        ! transverse drag
 real, private    :: k0 = 1.        ! radiation pressure
 real, private    :: k1 = 1.        ! redshift

 public  :: get_prdrag_spatial_force, get_prdrag_vdependent_force
 public  :: update_prdrag_leapfrog
 public  :: read_options_prdrag, write_options_prdrag

 private

contains

!------------------------------------------------
!+
!  compute the spatial part of the acceleration
!+
!------------------------------------------------
subroutine get_prdrag_spatial_force(xi,yi,zi,MStar,fextxi,fextyi,fextzi,phi)
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use lumin_nsdisc, only:beta
 use physcon, only: gg
 use units, only: udist, umass, utime
 real, intent(in)    :: xi,yi,zi,Mstar
 real, intent(inout) :: fextxi,fextyi,fextzi
 real, intent(out)   :: phi
 real                   :: r2,dr,dr3,betai,Mbdr3,rbetai,gcode

 gcode = gg / (udist**3/(utime**2*umass))

 r2 = xi*xi + yi*yi + zi*zi
 betai = beta(xi,yi,zi)
 rbetai = k0*betai
 if (r2 > epsilon(r2)) then
#ifdef FINVSQRT
    dr  = finvsqrt(r2)
#else
    dr = 1./sqrt(r2)
#endif
    dr3 = dr**3
    Mbdr3 = Mstar*gcode*(1.-rbetai)*dr3
    fextxi = fextxi - xi*Mbdr3
    fextyi = fextyi - yi*Mbdr3
    fextzi = fextzi - zi*Mbdr3
    phi    = -Mstar*gcode*dr*(1.-rbetai)
 endif

end subroutine get_prdrag_spatial_force

!-----------------------------------------------------------------------
!+
!  Routine to return velocity-dependent part of PR drag force
!+
!-----------------------------------------------------------------------
subroutine get_prdrag_vdependent_force(xyzi,vel,Mstar,fexti)
 use lumin_nsdisc, only: beta      !Change your Poynting-Robertson here.
 use physcon, only: gg, c
 use units, only: udist, umass, utime
 real, intent(in)  :: xyzi(3), vel(3)
 real, intent(in)  :: Mstar
 real, intent(out) :: fexti(3)
 real :: rhat(3)
 real                              :: betai, r, r2, r3, vr, gcode, ccode

 ccode = c  / (udist/utime)
 gcode = gg / (udist**3/(utime**2*umass))

 r2     = dot_product( xyzi, xyzi )
 r      = sqrt(r2)
 r3     = r*r2

 rhat = xyzi/r
 vr = dot_product(vel, rhat)

 betai  = beta( xyzi(1), xyzi(2), xyzi(3) )

 fexti = (-betai*Mstar*gcode/ccode)* &
            ( (vr/r3)*xyzi*k1 + vel/r2*k2 )

end subroutine get_prdrag_vdependent_force

subroutine update_prdrag_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dt,xi,yi,zi,Mstar)
 use lumin_nsdisc,  only:beta
 use physcon, only: c
 use units, only: udist, utime
 use io,            only:warn
 real, intent(in)    :: dt,xi,yi,zi, Mstar
 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(inout) :: fexti(3)
 real                              :: r, r2, r3, Q, betai
 real                              :: Tx, Ty, Tz, vonex, voney, vonez
 real                              :: denominator, vrhalf, vrone, twoQondt
 real                              :: xi2, yi2, zi2, ccode, kd
 character(len=30), parameter :: label = 'update_prdrag_leapfrog'

 ccode = c  / (udist/utime)

 xi2    = xi*xi
 yi2    = yi*yi
 zi2    = zi*zi
 kd     = k1 - k2

 r2     = (xi2 + yi2 + zi2)
 r      = sqrt(r2)
 r3     = r*r2
 vrhalf = vhalfx*xi + vhalfy*yi + vhalfz*zi

 betai       = beta( xi, yi, zi )
 Q           = Mstar*betai*dt/(2.*ccode*r*r)
 twoQondt    = 2.*Q/dt
 denominator = -r2*( k2*kd*Q*Q + (kd-k2)*Q - 1 )

 Tx = vhalfx + 0.5*dt*fxi
 Ty = vhalfy + 0.5*dt*fyi
 Tz = vhalfz + 0.5*dt*fzi

 vonex = (-(Q*k1*xi)*(Ty*yi+Tz*zi)+Q*kd*Tx*r2-Tx*(r2+Q*k1*xi2))/denominator
 voney = (-(Q*k1*yi)*(Tx*xi+Tz*zi)+Q*kd*Ty*r2-Ty*(r2+Q*k1*yi2))/denominator
 vonez = (-(Q*k1*zi)*(Tx*xi+Ty*yi)+Q*kd*Tz*r2-Tz*(r2+Q*k1*zi2))/denominator

 vrone = (vonex*xi + voney*yi + vonez*zi)/r      ! vr = rhat dot v

 fexti(1) = twoQondt * (vonex*k2 + k1*vrone*xi/r)
 fexti(2) = twoQondt * (voney*k2 + k1*vrone*yi/r)
 fexti(3) = twoQondt * (vonez*k2 + k1*vrone*zi/r)

 fxi      = fxi + fexti(1)
 fyi      = fyi + fexti(2)
 fzi      = fzi + fexti(3)

end subroutine update_prdrag_leapfrog

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_prdrag(iunit)
 use infile_utils,         only:write_inopt
 use lumin_nsdisc,         only:write_options_lumin_nsdisc
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to Poynting-Robertson drag'

 call write_inopt(k0, 'RadiationPressure', &
                  'Radiation pressure multiplier', iunit)
 call write_inopt(k2, 'TransverseDrag', &
                  'Transverse multiplier', iunit)
 call write_inopt(k1, 'Redshift', &
                  'Redshift multiplier', iunit)

 call write_options_lumin_nsdisc(iunit)

end subroutine write_options_prdrag

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_prdrag(name,valstring,imatch,igotall,ierr)
 use io,                   only:fatal, warning
 use lumin_nsdisc,         only:read_options_lumin_nsdisc
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_prdrag'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('RadiationPressure')
    read(valstring,*,iostat=ierr) k0
    ngot = ngot + 1
 case('TransverseDrag')
    read(valstring,*,iostat=ierr) k2
    ngot = ngot + 1
 case('Redshift')
    read(valstring,*,iostat=ierr) k1
    ngot = ngot + 1
 case default
    imatch = .false.
    call read_options_lumin_nsdisc(name,valstring,imatch,igotall,ierr)
 end select

 igotall = (ngot >= 1)

end subroutine read_options_prdrag

end module extern_prdrag
