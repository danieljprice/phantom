!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module damping
!
! Various implementations for artificial damping of velocities
! either to relax particles into equilibrium or enforce boundary conditions
!
! :References:
!   Gingold & Monaghan (1977) MNRAS 181, 375
!   Reichardt et al. (2019) MNRAS 484, 631 (used idamp=2)
!   González-Bolívar et al. (2022) MNRAS 517, 3181 (making idamp=2 obsolete)
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - damp   : *damping timescale as fraction of orbital timescale*
!   - idamp  : *artificial damping of velocities (0=off, 1=constant, 2=star, 3=disc)*
!   - r1in   : *inner boundary of inner disc damping zone*
!   - r1out  : *inner boundary of outer disc damping zone*
!   - r2in   : *outer boundary of inner disc damping zone*
!   - r2out  : *outer boundary of outer disc damping zone*
!   - tdyn_s : *dynamical timescale of star in seconds - damping is dependent on it*
!
! :Dependencies: infile_utils, physcon, units
!
 implicit none

 public  :: calc_damp,apply_damp,get_damp_fac_disc
 public  :: write_options_damping,read_options_damping

 private

 integer, public    :: idamp     = 0
 integer, parameter :: idamp_max = 3  ! maximum allowed value of idamp
 real, public :: tdyn_s      = 0.
 real, public :: damp        = 0.
 real, public :: r1in = 0.3
 real, public :: r2in = 0.357
 real, public :: r1out = 2.52
 real, public :: r2out = 3.0

contains

!-----------------------------------------------------------------------
!+
!  calculates damping factor
!+
!-----------------------------------------------------------------------
subroutine calc_damp(time, damp_fac)
 use units,   only:utime
 use physcon, only:pi
 real, intent(out)   :: damp_fac
 real, intent(in)    :: time
 real                :: tau1, tau2, tdyn_star, orbital_period

 select case(idamp)
 case(3)
    orbital_period = 2.*pi*sqrt(r1in**3)  ! G=M=1
    damp_fac = damp/orbital_period ! fraction of orbital time at r=r1in with G=M=1
 case(2)
    tdyn_star = tdyn_s / utime
    tau1 = tdyn_star * 0.1
    tau2 = tdyn_star
    if (time > 5. * tdyn_star) then
       damp_fac = 0.
    elseif (time > 2. * tdyn_star) then
       damp_fac = (tau1 * (tau2 / tau1)**((time - 2. * tdyn_star) / (3. * tdyn_star)))**(-1)
    else
       damp_fac = tau1**(-1)
    endif
 case(1)
    damp_fac = damp  ! fraction per timestep
 case default
    damp_fac = 0.
 end select

end subroutine calc_damp

!-----------------------------------------------------------------------
!+
!  apply damping term to velocity evolution (accelerations)
!+
!-----------------------------------------------------------------------
subroutine apply_damp(fextx, fexty, fextz, vxyz, xyz, damp_fac)
 real, intent(inout) :: fextx, fexty, fextz
 real, intent(in)    :: vxyz(3), xyz(3), damp_fac
 real :: v0(3),fac

 v0 = 0.
 fac = 1.
 !
 ! for idamp=3 damping is only applied in the boundary zones,
 ! hence damping factor depends on spatial location
 ! also in this case we relax to a prescribed velocity, not zero
 !
 if (idamp==3) fac = get_damp_fac_disc(xyz,v0)

 fextx = fextx - damp_fac*(vxyz(1)-v0(1))*fac
 fexty = fexty - damp_fac*(vxyz(2)-v0(2))*fac
 fextz = fextz - damp_fac*(vxyz(3)-v0(3))*fac

end subroutine apply_damp

!-----------------------------------------------------------------------
!+
!  radial damping zones for inner and outer boundary of a disc
!+
!-----------------------------------------------------------------------
real function get_damp_fac_disc(xyz,v0) result(fac)
 use physcon, only:pi
 real, intent(in) :: xyz(3)
 real, intent(out) :: v0(3)
 real :: rcyl,omega,vphi

 rcyl = sqrt(xyz(1)**2 + xyz(2)**2)

 if (rcyl < r2in) then
    fac = 1. - (sin(0.5*pi*(rcyl - r1in)/(r2in - r1in)))**2
 elseif (rcyl > r1out) then
    fac = (sin(0.5*pi*(rcyl - r1out)/(r2out - r1out)))**2*sqrt((r1in/r2out)**3)
 else
    fac = 0.
 endif

 omega = sqrt(1./rcyl**3)
 vphi = rcyl*omega

 v0(1) = -vphi*xyz(2)/rcyl   ! sin(phi) = y/R
 v0(2) =  vphi*xyz(1)/rcyl   ! cos(phi) = x/R
 v0(3) = 0.

end function get_damp_fac_disc

!-----------------------------------------------------------------------
!+
!  applies damping to external force
!+
!-----------------------------------------------------------------------
subroutine write_options_damping(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 ! do not write damping options if idamp == 0 (i.e. it is a hidden option)
 ! this is because it is mostly superseded by the relax-o-matic routines
 ! and can lead to confusion
 if (idamp <= 0) return

 write(iunit,"(/,a)") '# options controlling damping'
 call write_inopt(idamp,'idamp','artificial damping of velocities (0=off, 1=constant, 2=star, 3=disc)',iunit)

 select case(idamp)
 case(1)
    call write_inopt(damp,'damp','fraction to damp velocity each timestep (if > 0, v=0 initially)',iunit)
 case(2)
    call write_inopt(tdyn_s,'tdyn_s','dynamical timescale of star in seconds - damping is dependent on it',iunit)
 case(3)
    call write_inopt(damp,'damp','damping timescale as fraction of orbital timescale',iunit)
    call write_inopt(r1in,'r1in','inner boundary of inner disc damping zone',iunit)
    call write_inopt(r2in,'r2in','outer boundary of inner disc damping zone',iunit)
    call write_inopt(r1out,'r1out','inner boundary of outer disc damping zone',iunit)
    call write_inopt(r2out,'r2out','outer boundary of outer disc damping zone',iunit)
 end select

end subroutine write_options_damping

!-----------------------------------------------------------------------
!+
!  reads damping options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_damping(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(idamp,'idamp',db,errcount=nerr,min=0,max=idamp_max,default=0)
 if (idamp > 0) call read_inopt(damp,'damp',db,errcount=nerr,min=0.,max=1.,default=0.)
 select case(idamp)
 case(3)
    call read_inopt(r1in,'r1in',db,errcount=nerr,min=0.)
    call read_inopt(r2in,'r2in',db,errcount=nerr,min=r1in)
    call read_inopt(r1out,'r1out',db,errcount=nerr,min=r2in)
    call read_inopt(r2out,'r2out',db,errcount=nerr,min=r1out)
 case(2)
    call read_inopt(tdyn_s,'tdyn_s',db,errcount=nerr,min=0.)
 end select

end subroutine read_options_damping

end module damping
