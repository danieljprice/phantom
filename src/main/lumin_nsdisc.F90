!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: lumin_nsdisc
!
!  DESCRIPTION:
! This module contains routines for calculating beta, the
! ratio of radiation to gravitational force, for an accretion disc
! surrounding a neutron star. It contains associated functions
! for calculating opacity, accretion luminosity, etc.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, fastmath, infile_utils, io, physcon, units
!+
!--------------------------------------------------------------------------

module lumin_nsdisc
 use physcon, only: pi
 implicit none

 real :: LumAcc                   = 0.0                 ! luminosity from accretion
 integer, public  :: BurstProfile = 0                   ! Burst luminosity / profile
 real, private :: Lstar           = 0.0                 ! total luminosity of star as a fraction of Eddington
 real, private :: LEdd            = 0.0                 ! Eddington luminosity in code units
 real, private :: frac_X          = 0.7                 ! Hydrogen fraction
 integer, private :: rad_trans    = 0                   ! Radiation transport prescription
 real          :: AccLumEff       = 1.0 ! Fraction of accretion luminosity that
 ! participates in radiation force. This is
 ! mostly only relevant for compact objects.
 ! number of gridpoints for the rho, tau, beta grids
 integer, parameter :: nr = 63, nth=16 , nph=24, nz=30
 real :: ringrid(0:nr), cyl_ringrid(0:nr), thmingrid(0:nth), phmingrid(0:nph), zmingrid(0:nz)
 integer :: npgrid(0:nr-1,0:nth-1,0:nph-1), w92npart( 0:nr-1 ), cyl_npgrid(0:nr-1,0:nz,0:nph-1)
 real :: densitygrid(0:nr-1,0:nth-1,0:nph-1), ldensitygrid(0:nr-1,0:nth-1,0:nph-1)
 real :: tauradgrid(0:nr-1,0:nth-1,0:nph-1), ltauradgrid(0:nr-1,0:nth-1,0:nph-1)
 real :: betagrid(0:nr-1,0:nth-1,0:nph-1), lbetagrid(0:nr-1,0:nth-1,0:nph-1) !the grids

 real :: cyl_densitygrid(0:nr-1,0:nz,0:nph-1), cyl_ldensitygrid(0:nr-1,0:nz,0:nph-1)
 real :: cyl_tauradgrid(0:nr-1,0:nz,0:nph-1), cyl_ltauradgrid(0:nr-1,0:nz,0:nph-1)
 real :: cyl_betagrid(0:nr-1,0:nz,0:nph-1), cyl_lbetagrid(0:nr-1,0:nz,0:nph-1) !cylindrical grids

 real :: w92betagrid( 0:nr-1 ), w92sumbeta( 0:nr-1 )                         !cylindrical grid for comparison with W92
 real :: thetamin = 0.0, thetamax = pi, rmin = 1.0, rmax = 1001.0
 real :: zmin = -200, zmax=200
 real :: phimax = 2*pi, phimin=0
 real :: Lstar_burst
 real, parameter :: eps = 1.e-6
 integer, private :: made_grid_points = 0

 public :: beta, set_Lstar, calc_sigma, calc_scaleheight
 public :: read_options_lumin_nsdisc, write_options_lumin_nsdisc
 public :: LumAcc, Lstar_burst, AccLumEff, ringrid, thmingrid, thetamin, thetamax, rmin, rmax
 public :: nr, nth, nph, densitygrid, tauradgrid, betagrid, lbetagrid, make_beta_grids
 public :: get_grid_points, bilin_interp, get_grid_bins, sphere_segment, get_bracket_grid_points
 public :: ldensitygrid, ltauradgrid, careful_log, phmingrid, phimin, phimax, w92betagrid
 public :: cyl_ringrid, zmingrid, cyl_npgrid, cyl_densitygrid, cyl_ldensitygrid, cyl_tauradgrid
 public :: cyl_ltauradgrid, cyl_betagrid, cyl_lbetagrid, zmin, zmax, nz

 private :: calc_kappa, eps

 private

contains

!----------------------------------------------------------------
!+
!  Sets the location of the grid points
!+
!----------------------------------------------------------------

subroutine make_grid_points()
 use physcon, only:pi, twopi
 integer :: rbin, thbin, phbin, zbin
 real    :: A, B, C, tempr
 A = (rmax-rmin)/(nr*nr)
 B = 2.*(thetamin-thetamax)/nth
 C = 2.*(zmin-zmax)/nz

 do rbin=0, nr-1
    tempr = rmin + A*rbin**2
    ringrid(rbin) = tempr
    w92betagrid(rbin) = tempr
    cyl_ringrid(rbin) = tempr
 enddo

 thmingrid(0)=thetamin
 do thbin=1,nth/2-1
    thmingrid(thbin) = B*(1.0*thbin*thbin/nth-thbin)+thetamin
    thmingrid(nth-thbin) = thetamax - thmingrid(thbin)
 enddo
 thmingrid(nth/2) = (thetamin+thetamax)/2;

 zmingrid(0)=zmin
 do zbin=1,nz/2-1
    zmingrid(zbin) = C*(1.0*zbin*zbin/nz-zbin)+zmin
    zmingrid(nz-zbin) = (zmax+zmin)/2. - zmingrid(zbin)
 enddo
 zmingrid(nz/2) = (zmin+zmax)/2.

 do phbin=0,nph-1
    phmingrid(phbin) = phbin * twopi/nph
 enddo

 made_grid_points = 1

end subroutine make_grid_points

!----------------------------------------------------------------
!+
!  Given a set of coordinates r, theta, phi, finds the cell
!  those coords are in
!+
!----------------------------------------------------------------

subroutine get_grid_bins( r, zt, rbin, ztbin, phi, phibin )
 use physcon, only:pi, twopi
 use io, only : fatal
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 real,    intent(in)  :: r, phi, zt
 integer, intent(out) :: rbin, ztbin, phibin
 real :: B, C, ztnew
#ifdef FINVSQRT
 rbin  = int( nr*finvsqrt( (rmax-rmin)/((r-rmin))))    !optimized for speed not readability
#else
 rbin  = int( nr*sqrt( (r-rmin)/(rmax-rmin)))
#endif

 B = 2.*(thetamin-thetamax)/(nth)
 C = 2.*(zmin-zmax)/nz

 select case(rad_trans)
 case(0,1)
    if( zt < (thetamin+thetamax)/2. ) then
       ztbin = int( (sqrt( (nth*B)**2 + 4*nth*B*(zt-thetamin) ) + nth*B)/(2.*B) )
    else
       ztnew = thetamin + thetamax - zt
       ztbin = int( (sqrt( (nth*B)**2 + 4*nth*B*(ztnew-thetamin) ) + nth*B)/(2.*B) )
       ztbin = nth-ztbin-1
    endif
    if( ztbin < 0 .or. ztbin>nth-1 ) then
       call fatal( 'lumin_nsdisc', 'Array out of bounds error in get_grid_bins (theta)' )
    endif
 case(2)
    if( zt < (zmin+zmax)/2. ) then
       ztbin = int( (sqrt( (nz*C)**2 + 4*nz*C*(zt-zmin) ) + nz*C)/(2.*C) )
    else
       ztnew = zmin + zmax - zt
       ztbin = int( (sqrt( (nz*C)**2 + 4*nz*C*(ztnew-zmin) ) + nz*C)/(2.*C) )
       ztbin = nz-ztbin-1
    endif
    if( ztbin>nz ) ztbin = nz
    if( ztbin<0 ) ztbin = 0

 end select
 phibin = int( phi*nph/twopi)
 if( rbin>nr-1 ) rbin=nr-1        ! Avoids segfaults for distant particles
 if( rbin<0 ) rbin = 0            ! Avoids segfaults for accreted particles

end subroutine get_grid_bins

!----------------------------------------------------------------
!+
!  Given a bin in an array, finds the inner and outer edges, and the
!  midpoint
!+
!----------------------------------------------------------------

subroutine get_grid_points( array, ix, nx, maxx, xin, xout, xmid )
 use io, only : fatal
 integer, intent(in)  :: nx, ix
 real,    intent(in)  :: array(0:nx-1), maxx
 real,    intent(out) :: xin, xout, xmid
 if( ix<0.or.ix>=nx ) then
    call fatal( 'lumin_nsdisc', 'Array out of bounds error in get_grid_points' )
 endif
 xin = array(ix)

 if( ix==nx-1 ) then
    xout = maxx
 else
    xout = array(ix+1)
 endif

 xmid = ( xin + xout )/2.

end subroutine get_grid_points

!----------------------------------------------------------------
!+
!  Given a bin in an array, finds the midpoints of that bin and the next
!  one out. This is called in preparation for bilin_interp
!+
!----------------------------------------------------------------

subroutine get_bracket_grid_points( array, ix, nx, maxx, x1, x2 )
 use io, only : fatal
 integer, intent(in)  :: nx, ix
 real,    intent(in)  :: array(0:nx-1), maxx
 real,    intent(out) :: x1, x2
 real :: minimum, maximum, boundary

 if( ix<0.or.ix>=nx-1 ) then
    call fatal( 'lumin_nsdisc', 'Array out of bounds error in get_bracket_grid_points' )
    !this should never happen. Checks in the calling function should avoid passing bad ix values.
 endif

 minimum = array(ix)
 boundary = array(ix+1)

 if( ix==nx-2 ) then
    maximum = maxx
 else
    maximum = array(ix+2)
 endif

 x1 = minimum  + (boundary-minimum)/2.
 x2 = boundary + (maximum-boundary)/2.

end subroutine

!----------------------------------------------------------------
!+
!  Generates a set of grids containing rho, tau, and beta
!  Calculates density by counting the number of particles in
!  each spherical r,theta bin
!  Calculates tau by integrating radially from NS surface
!  Calculates beta = exp(-tau). You still need to multiply beta
!  by L*/LEdd to get the true beta used to calculate PR drag
!+
!----------------------------------------------------------------

subroutine make_beta_grids(xyzh,particlemass,npart)
 use units, only: udist, umass

 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 real,    intent(in) :: particlemass
 real :: theta, r, phi, kappa, r_cyl
 real :: dr, dz, cell_phin, cell_phout, cell_phmid
 real :: cell_rin, cell_rout, cell_rmid, cell_thin, cell_thout, cell_thmid
 real :: cell_zin, cell_zout, cell_zmid
 real :: cell_volume, logzero, x, y, z
 integer :: rbin, thbin, phbin, zbin, ipart, totpart

 kappa = calc_kappa( frac_X, 2. ) * (umass/(udist*udist))   !kappa in code units

 logzero = careful_log( 0.0 )

 call make_grid_points()
 do rbin=0, nr-1   !clear the density grid, since it is retained after every call
    do phbin=0, nph-1
       do thbin=0, nth-1
          ldensitygrid( rbin, thbin, phbin ) = logzero
          densitygrid( rbin, thbin, phbin ) = 0.0
          npgrid( rbin, thbin, phbin ) = 0
       enddo
       do zbin=0, nz
          cyl_ldensitygrid(rbin, zbin, phbin) = logzero
          cyl_densitygrid(rbin, zbin, phbin) = 0.0
          cyl_npgrid(rbin, zbin, phbin) = 0
       enddo
    enddo
 enddo

 totpart = 0
 do ipart=1, npart    !fills the density grid by counting particles and assigning
    !each one to a cell
    x = xyzh(1, ipart)
    y = xyzh(2, ipart)
    z = xyzh(3, ipart)

    if( isnan(x).or.isnan(y).or.isnan(z) ) then
       x=0.
       y=0.
       z=0.
    endif

    r = sqrt( x**2 + y**2 + z**2 )
    r_cyl = sqrt(x**2 + y**2)

    if( r>rmin.and.r<rmax ) then

       theta = acos( z/r )
       if( theta > thetamax ) theta = thetamax
       if( theta < thetamin ) theta = thetamin

       phi = pi + atan2( y, x )
    endif
    if (phbin >= nph ) then   !phi is cyclic
       phbin = phbin - nph
    endif
    cell_phin  = phmingrid( phbin )
    if( phbin == nph-1 ) then
       cell_phout = phimax
    else
       cell_phout = phmingrid( phbin+1 )
    endif
    if( r>rmin.and.r<rmax ) then
       select case(rad_trans)
       case(0,1)
          call get_grid_bins( r, theta, rbin, thbin, phi, phbin )

          cell_rin   = ringrid( rbin )
          cell_rout  = ringrid( rbin + 1 )

          if( rbin == nr-1 ) then
             cell_rout = rmax
          endif

          cell_thin  = thmingrid( thbin )
          cell_thout = thmingrid( thbin + 1 )

          if( thbin == nth-1 ) then
             cell_thout = thetamax
          endif
          cell_rmid  = (cell_rout+cell_rin)/2.
          cell_thmid = (cell_thout+cell_thin)/2.
          cell_phmid = (cell_phout+cell_phin)/2.
          cell_volume = sphere_segment( cell_rin, cell_rout, cell_thin, cell_thout, cell_phin, cell_phout )
          if ( (rbin>=0).and.(rbin<nr).and.&
              (thbin>=0).and.(thbin<nth).and.&
              (phbin>=0).and.(phbin<nph) ) then

             npgrid( rbin, thbin, phbin ) = npgrid( rbin, thbin, phbin ) + 1
             densitygrid( rbin, thbin, phbin ) = densitygrid( rbin, thbin, phbin ) + particlemass/(cell_volume)
             totpart = totpart + 1
          endif

       case(2)
          if( z>zmax ) z=zmax-eps
          if( z<zmin ) z=zmin+eps
          call get_grid_bins( r_cyl, z, rbin, zbin, phi, phbin )

          cell_rin   = cyl_ringrid( rbin )
          cell_rout  = cyl_ringrid( rbin + 1 )
          if( rbin == nr-1 ) then
             cell_rout = rmax
          endif
          cell_zin = zmingrid( zbin )
          cell_zout = zmingrid( zbin+1)
          if( zbin == nz ) then
             cell_zout = zmax
          endif
          cell_volume = cylinder_segment( cell_rin, cell_rout, cell_zin, cell_zout, cell_phin, cell_phout )
          if ( (rbin>=0).and.(rbin<nr).and.&
              (zbin>=0).and.(zbin<nz).and.&
              (phbin>=0).and.(phbin<nph) ) then

             cyl_npgrid( rbin, zbin, phbin ) = cyl_npgrid( rbin, zbin, phbin ) + 1
             cyl_densitygrid( rbin, zbin, phbin ) = cyl_densitygrid( rbin, zbin, phbin ) + particlemass/(cell_volume)
             totpart = totpart + 1
          endif
       end select
    endif
 enddo

 if (rad_trans==0 .or. rad_trans==1 ) then
    !Now calculate tau and beta
    do phbin=0, nph-1
       do thbin=0, nth-1
          call get_grid_points( ringrid, 0, nr, rmax, cell_rin, cell_rout, cell_rmid )
          call get_grid_points( thmingrid, thbin, nth, thetamax, cell_thin, cell_thout, cell_thmid )
          call get_grid_points( phmingrid, phbin, nph, phimax, cell_phin, cell_phout, cell_phmid )
          dr = (cell_rout-cell_rin)/2.

          ldensitygrid( 0, thbin, phbin ) = careful_log(densitygrid( 0, thbin, phbin ) )
          tauradgrid(0, thbin, phbin) = 0.
          ltauradgrid(0, thbin, phbin) = careful_log(tauradgrid(0, thbin, phbin))
          betagrid( 0, thbin, phbin ) = tau_to_beta(tauradgrid(0, thbin, phbin))
          lbetagrid(0, thbin, phbin ) = careful_log(betagrid(0, thbin, phbin))
          do rbin=1, nr-2
             call get_grid_points( ringrid, rbin, nr, rmax, cell_rin, cell_rout, cell_rmid )
             ldensitygrid( rbin, thbin, phbin ) = careful_log(densitygrid( rbin, thbin, phbin ) )
             dr = (cell_rout-cell_rin)/2.0
             tauradgrid( rbin, thbin, phbin ) = tauradgrid( rbin-1, thbin, phbin ) + &
                   kappa * dr * ( densitygrid(rbin,thbin,phbin)+densitygrid(rbin+1,thbin,phbin) )

             ltauradgrid( rbin, thbin, phbin ) = careful_log(tauradgrid( rbin, thbin, phbin ))
             betagrid( rbin, thbin, phbin ) = tau_to_beta( tauradgrid(rbin, thbin, phbin ) )
             lbetagrid( rbin, thbin, phbin ) = careful_log(betagrid( rbin, thbin, phbin ))
          enddo
          dr = (rmax - ringrid(nr-1))/2.0
          call get_grid_points( ringrid, nr-1, nr, rmax, cell_rin, cell_rout, cell_rmid )
          ldensitygrid( nr-1, thbin, phbin ) = careful_log(densitygrid( nr-1, thbin, phbin ))
          tauradgrid( nr-1, thbin, phbin ) = tauradgrid( nr-2, thbin, phbin )&
         + kappa * dr * densitygrid( nr-1, thbin, phbin)/2.
          ltauradgrid( nr-1, thbin, phbin ) = careful_log(tauradgrid( nr-1, thbin, phbin ))
          betagrid( nr-1, thbin, phbin ) = tau_to_beta( tauradgrid(nr-1, thbin, phbin ) )
          lbetagrid( nr-1, thbin, phbin ) = careful_log(betagrid( nr-1, thbin, phbin ))
       enddo
    enddo
 endif
 if (rad_trans==1) then
    do rbin=0, nr-1
       w92sumbeta(rbin) = 0.0
       w92npart(rbin) = 0
    enddo
    do ipart=0, npart
       x = xyzh(1, ipart)
       y = xyzh(2, ipart)
       z = xyzh(3, ipart)
       if( isnan(x).or.isnan(y).or.isnan(z) ) then
          x=0
          y=0
          z=0
       endif
       r = sqrt( x**2 + y**2 + z**2 )
       r_cyl = sqrt( x**2 + y**2)
       if (r==0.0) then
          theta=0.0
       else
          theta = acos( z/r )
       endif
       if( theta > thetamax ) theta = thetamax
       if( theta < thetamin ) theta = thetamin
       phi = pi + atan2( y, x )
       if( r_cyl>= 0) then
          call get_grid_bins( r_cyl, 0., rbin, thbin, 0., phbin )
       else
          rbin=-1
          thbin=-1
          phbin=-1
       endif
       if( rbin>=0.and.rbin<nr) then
          w92npart(rbin) = w92npart(rbin)+1
          w92sumbeta(rbin) = w92sumbeta(rbin) + beta_by_interp(r, theta, phi)
       endif
    enddo

    do rbin=0, nr-1
       w92betagrid(rbin) = w92sumbeta(rbin)/w92npart(rbin)
    enddo

 endif
 if (rad_trans==2) then        !Case for vertical integration on a cylindrical grid
    do rbin=0, nr-1
       do phbin=0, nph-1
          !We want to integrate vertically from zmin to the midplane, then from zmax to midplane.
          !Start at zbin=0
          call get_grid_points( cyl_ringrid, rbin, nr, rmax, cell_rin, cell_rout, cell_rmid )
          call get_grid_points( phmingrid, phbin, nph, phimax, cell_phin, cell_phout, cell_phmid )
          call get_grid_points( zmingrid, 0, nz, zmax, cell_zin, cell_zout, cell_zmid )

          cyl_ldensitygrid( rbin, 0, phbin ) = careful_log(cyl_densitygrid( rbin, 0, phbin ) )
          cyl_tauradgrid(rbin, 0, phbin) = 0.
          cyl_ltauradgrid(rbin, 0, phbin) = careful_log(cyl_tauradgrid(rbin, 0, phbin))
          cyl_betagrid( rbin, 0, phbin ) = tau_to_beta(cyl_tauradgrid(rbin, 0, phbin))
          cyl_lbetagrid(rbin, 0, phbin ) = careful_log(cyl_betagrid(rbin, 0, phbin))
          do zbin=1, nz/2-1
             call get_grid_points( zmingrid, zbin, nz, zmax, cell_zin, cell_zout, cell_zmid )
             cyl_ldensitygrid( rbin, zbin, phbin ) = careful_log(cyl_densitygrid( rbin, zbin, phbin ) )
             dz = (cell_zout-cell_zin)/2.0
             cyl_tauradgrid( rbin, zbin, phbin ) = cyl_tauradgrid( rbin, zbin-1, phbin ) + &
                   kappa * dz * ( cyl_densitygrid(rbin,zbin,phbin)+cyl_densitygrid(rbin,zbin+1,phbin) )

             cyl_ltauradgrid( rbin, zbin, phbin ) = careful_log(cyl_tauradgrid( rbin, zbin, phbin ))
             cyl_betagrid( rbin, zbin, phbin ) = tau_to_beta( cyl_tauradgrid(rbin, zbin, phbin ) )
             cyl_lbetagrid( rbin, zbin, phbin ) = careful_log(cyl_betagrid( rbin, zbin, phbin ))

          enddo
          !Now go from zmax to midplane
          call get_grid_points( cyl_ringrid, rbin, nr, rmax, cell_rin, cell_rout, cell_rmid )
          call get_grid_points( phmingrid, phbin, nph, phimax, cell_phin, cell_phout, cell_phmid )
          call get_grid_points( zmingrid, nz-1, nz, zmax, cell_zin, cell_zout, cell_zmid )
          dz=(cell_zout-cell_zin)/2.

          cyl_ldensitygrid( rbin, nz-1, phbin ) = careful_log(cyl_densitygrid( rbin, nz-1, phbin ) )
          cyl_tauradgrid(rbin, nz-1, phbin) = 0.
          cyl_ltauradgrid(rbin, nz-1, phbin) = careful_log(cyl_tauradgrid(rbin, nz-1, phbin))
          cyl_betagrid( rbin, nz-1, phbin ) = tau_to_beta(cyl_tauradgrid(rbin, nz-1, phbin))
          cyl_lbetagrid(rbin, nz-1, phbin ) = careful_log(cyl_betagrid(rbin, nz-1, phbin))
          do zbin=nz-2, nz/2, -1
             call get_grid_points( zmingrid, zbin, nz, zmax, cell_zin, cell_zout, cell_zmid )
             cyl_ldensitygrid( rbin, zbin, phbin ) = careful_log(cyl_densitygrid( rbin, zbin, phbin ) )
             dz = (cell_zout-cell_zin)/2.0
             cyl_tauradgrid( rbin, zbin, phbin ) = cyl_tauradgrid( rbin, zbin+1, phbin ) + &
                   kappa * dz * ( cyl_densitygrid(rbin,zbin,phbin)+cyl_densitygrid(rbin,zbin-1,phbin) )
             cyl_ltauradgrid( rbin, zbin, phbin ) = careful_log(cyl_tauradgrid( rbin, zbin, phbin ))
             cyl_betagrid( rbin, zbin, phbin ) = tau_to_beta( cyl_tauradgrid(rbin, zbin, phbin ) )
             cyl_lbetagrid( rbin, zbin, phbin ) = careful_log(cyl_betagrid( rbin, zbin, phbin ))
          enddo
       enddo
    enddo
 endif
end subroutine make_beta_grids

real function sphere_segment( rin, rout, thin, thout, phin, phout )
 !returns the volume of a segment of a sphere with given
 !inner and outer radii, and limiting spherical angles
 real, intent(in) :: rin, rout, thin, thout, phin, phout
 sphere_segment = (rout**3 - rin**3)*(cos(thin) - cos(thout))*(phout-phin)
end function sphere_segment

real function cylinder_segment( rin, rout, zin, zout, phin, phout )
 real, intent(in) :: rin, rout, zin, zout, phin, phout
 cylinder_segment = abs( zout-zin )*(phout-phin)*(rout**2-rin**2)/2.
end function cylinder_segment

real function tau_to_beta( tau )
 !Returns exp(-tau), but with error handling
 !for very high tau
 real, intent(in) :: tau
 real, parameter  :: taubig = range(tau)*log(10.)

 if( tau >= taubig ) then
    tau_to_beta = 0.
 else
    tau_to_beta = exp(-tau)
 endif

end function tau_to_beta

!------------------------------------------------
!+
! Calculates the total luminosity of the star,
! including perhaps luminosity from accretion
!+
!------------------------------------------------

subroutine set_Lstar( BurstProfile, time, dmdt, Mstar )
 use units, only:utime, umass, udist
 use physcon, only:fourpi
 real,    intent(in) :: time, dmdt, Mstar
 integer, intent(in) :: BurstProfile
 real             :: ptime, ptime2

!this assumes c=G=1.
 LEdd = fourpi*Mstar/(calc_kappa( frac_X, 2. ) / ( udist*udist / umass ))
 ptime  = time*utime
 ptime2 = ptime*ptime

 select case( BurstProfile )

 case(-1)        ! Test case. I will modify this one frequently.
    if( time < 9999.) then
       Lstar_burst=0.0
       AccLumEff=0.0
    else
       Lstar_burst=1.0
       AccLumEff=0.0
    endif
 case(0)         ! No luminosity, either from burning or from accretion feedback
    Lstar_burst = 0.0
    AccLumEff   = 0.0
 case(1)         !Time-variable luminosity profile.
    ! 00.00 - 00.25 = no luminosity at all to remove initial transient
    ! 00.25 - 01.00 = no L*, linearly ramp up accretion feedback from 0 to 1
    ! 01.00 - 01.50 = Linear rise in L* from 0 to LEdd/2
    ! 01.50 - 10.00 = Quadratic decay back to zero
    ! 10.00 -       = L* = 0, AccLumEff = 1.
    if( ptime < 0.25 ) then
       Lstar_burst = 0.0
       AccLumEff   = 0.0
    else if( ptime < 1.0  ) then
       Lstar_burst = 0.0
       AccLumEff   = (ptime*4. - 1.)/3.
    else if( ptime < 1.5  ) then
       Lstar_burst = (ptime - 1.)
       AccLumEff   = 1.
    else if( ptime < 10. ) then
       Lstar_burst = (2*ptime2 - 40*ptime + 200)/289.0
       AccLumEff   = 1.
    else
       Lstar_burst = 0.0
       AccLumEff   = 1.
    endif
 case(2)    ! No burning, but ramp up AccLumEff after initial transient
    Lstar_burst = 0.0
    if( ptime < 0.25 ) then
       AccLumEff   = 0.0
    else if( ptime < 1.0  ) then
       AccLumEff   = (ptime*4. - 1.)/3.
    else
       AccLumEff   = 1.0
    endif
 case(3)    ! No burning, but ramp up AccLumEff half as fast as in case(2)
    Lstar_burst = 0.0
    if( ptime < 0.25 ) then
       AccLumEff   = 0.0
    else if( ptime < 1.75  ) then
       AccLumEff   = (ptime*4. - 1.)/6.
    else
       AccLumEff   = 1.0
    endif
 case(4)    ! Remove initial transient; ramp up acclum over 3/4 second; allow to settle to 10s;
    ! impose a half eddington burst at 10s
    Lstar_burst = 0.0
    if( ptime < 0.25 ) then        !remove initial transient
       AccLumEff   = 0.0
    else if( ptime < 1.00  ) then  !ramp up AccLum
       AccLumEff   = (ptime*4. - 1.)/3.
    else if( ptime < 10.00 ) then  !Settle disc
       AccLumEff   = 1.0
    else if( ptime < 10.5 ) then   !Burst rise
       AccLumEff = 1.0
       Lstar_burst = (ptime-10)
    else if( ptime < 20.5 ) then   !Burst decay
       AccLumEff = 1.0
       Lstar_burst = (4*ptime2 - 164*ptime + 1681)/800.
    else                           !Post burst
       AccLumEff   = 1.0
       Lstar_burst = 0.0
    endif
 case(5)   ! Eddington luminosity from beginning of simulation
    Lstar_burst = 1.0
    AccLumEff   = 1.0
 case(6)    ! Remove initial transient; ramp up acclum over 3/4 second; allow to settle to 10s;
    ! impose a half eddington burst at 10s
    Lstar_burst = 0.0
    if( ptime < 0.25 ) then        !remove initial transient
       AccLumEff   = 0.0
    else if( ptime < 1.00  ) then  !ramp up AccLum
       AccLumEff   = (ptime*4. - 1.)/3.
    else if( ptime < 10.00 ) then  !Settle disc
       AccLumEff   = 1.0
    else if( ptime < 10.5 ) then   !Burst rise
       AccLumEff = 1.0
       Lstar_burst = (ptime-10)*2
    else if( ptime < 20.5 ) then   !Burst decay
       AccLumEff = 1.0
       Lstar_burst = (4*ptime2 - 164*ptime + 1681)/400.
    else                           !Post burst
       AccLumEff   = 1.0
       Lstar_burst = 0.0
    endif
 case(7)    !Fit to a type-1 nonPRE burst from 1636-536. Since I'll be starting from an already settled simulation, no need to adjust time
    ptime2  = ptime-3
    Lstar_burst = 0.0
    AccLumEff = 0.0
    if( ptime > 2. .and. ptime <3. ) then
       Lstar_burst=0
       AccLumEff = ptime-2.
    endif
    if( ptime>3 ) then
       Lstar_burst = 9.4*ptime2**1.6*exp(-2.9*ptime2**0.4)
       AccLumEff = 1.
    endif
 case(8)     !0.1 L_Edd from beginning of simulation, used for prtest
    Lstar_burst = 0.10
    AccLumEff   = 0.00
 end select

 if( AccLumEff > 1.e-16 ) then     !If we are including luminosity feedback
    LumAcc = AccLumEff * get_AccLum( dmdt, Mstar )
 else
    LumAcc = 0.0
 endif
 Lstar = Lstar_burst + LumAcc

end subroutine set_Lstar

!----------------------------------------------------
!+
!   Converts accretion rate into luminosity
!+
!----------------------------------------------------

real function get_AccLum( dmdt, Mstar )       !Luminosity from accretion
 use physcon, only:gg
 use units, only:udist,umass,utime
 real, intent(in) :: dmdt, Mstar
 real :: ggcode, Rstar

 ggcode = gg / (udist**3/(utime**2*umass))
 Rstar = 1.

 get_AccLum =  (ggcode * Mstar * dmdt / Rstar) / LEdd

end function get_AccLum

!----------------------------------------------
!+
!  function computing kappa opacity from hydrogen
!  fraction X and temperature kT (in keV).
!  Returns opacity in cm^2 g^{-1}
!
!  See Lewin et al 1993, SSR, 62, 233 (p. 276 4.13b)
!+
!----------------------------------------------

real function calc_kappa( X, kT )
 real, intent(in) :: X, kT
 real :: k0, tempcorr
 k0 = 0.2*(1.0 + X)
 tempcorr = 1.0 + ( kT/39.2 )**0.86
 calc_kappa = k0/tempcorr
end function calc_kappa

!----------------------------------------------
!+
!  function computing the beta parameter
!+
!----------------------------------------------
real function beta(x,y,z)
 use physcon,              only:c, gg, fourpi, pi, roottwo, rpiontwo
 use io,                   only:fatal
 use units,   only:umass,udist
 real, intent(in) :: x,y,z
 real             :: r, theta, phi, rcyl, H, kappa, tau
 integer          :: rbin, thetabin, phibin, zbin

 beta = 0.
 rcyl = sqrt( x*x + y*y )
 phi = pi + atan2(y,x)
 r = sqrt(x**2 + y**2 + z**2)
 theta=acos(z/r)

 select case( rad_trans )
 case( 0 )
    r = sqrt(x**2 + y**2 + z**2)
    theta = acos( z/r )
    if( theta > thetamax ) theta = thetamax
    if( theta < thetamin ) theta = thetamin
    if( r>rmin.and.r<rmax ) then
       if( made_grid_points == 0 ) then
          beta = 0.
       else
          beta = beta_by_interp(r, theta, phi)
       endif
    else
       beta = 0.
    endif

 case( 1 )
    rcyl = sqrt( x*x + y*y )
    if( rcyl>rmin.and.rcyl<rmax ) then
       call get_grid_bins( rcyl, 0., rbin, thetabin, 0., phibin )

       beta = w92betagrid( rbin )
    else
       beta = 0
    endif

 case( 2 )
    call get_grid_bins( rcyl, z, rbin, zbin, phi, phibin )
    beta = beta_by_interp_cyl(rcyl, z, phi)
 case( 3 )
    r = sqrt(x**2 + y**2)
    H = calc_scaleheight(r)
    kappa = calc_kappa( frac_X, 2. ) * (umass/(udist*udist))
    tau = 1.0-erf( abs(z) / (roottwo*H) )
    tau = tau * rpiontwo * kappa * calc_sigma(r) * H
    beta = exp(-tau)

 end select

 beta = beta * Lstar

 if( beta > Lstar ) then
    beta = Lstar    !hopefully unnecessary sanity checks
 endif
 if( beta < 0. ) beta = 0.

end function beta

!----------------------------------------------
!+
!  Calculates a beta by calling bilin_interp
!+
!---------------------------------------------
real function beta_by_interp(r, theta, phi)
 real, intent(in) :: r, theta, phi
 real             :: betain, betaout
 integer          :: rbin, thetabin, phibin

 beta_by_interp = 0
 rbin=0
 thetabin=0
 phibin=0
 call get_grid_bins( r, theta, rbin, thetabin, phi, phibin )

 if( rbin >= nr-1 ) rbin = nr-1

 if( r<=rmin ) then
    betain  = 0.
    betaout = 0.
 else if( frac_X < 0. ) then
    betain  = 1.
    betaout = 1.
 else if( rbin>=nr-1 ) then
    betain  = lbetagrid(nr-1, thetabin, phibin)
    betain  = min(exp(betain), 1.)
    betaout = betain
 else
    betain  = bilin_interp( lbetagrid, ringrid(rbin), theta, phi)
    betain  = exp(betain)

    betaout = bilin_interp( lbetagrid, ringrid(rbin+1), theta, phi )
    betaout = exp(betaout)
 endif

 if( npgrid( rbin, thetabin, phibin )==0) then
    beta_by_interp = betain
 else
    beta_by_interp = (betaout + (betain-betaout)/npgrid( rbin, thetabin, phibin ))

 endif

end function beta_by_interp

!----------------------------------------------
!+
!  Calculates a beta by calling bilin_interp on a cylindrical grid
!+
!---------------------------------------------
real function beta_by_interp_cyl(r, z, phi)
 real, intent(in) :: r, z, phi
 real             :: betain, betaout
 integer          :: rbin, zbin, phibin
 real             :: znew
 znew=z
 if( z > zmax ) znew = zmax - eps
 if( z < zmin ) znew = zmin + eps
 rbin=0
 zbin=0
 phibin=0
 call get_grid_bins( r, z, rbin, zbin, phi, phibin )
 if( r<=rmin ) then
    betain  = 0.
    betaout = 0.
 else if( frac_X < 0. ) then
    betain  = 1.
    betaout = 1.
 else if( rbin>=nr-1 ) then
    betain  = bilin_interp_cyl( cyl_lbetagrid, r, z, phi )
    betain  = min(exp(betain), 1.)
    betaout = betain
 else
    betain  = bilin_interp_cyl( cyl_lbetagrid, r, z, phi )
    betain  = exp(betain)
    betaout = bilin_interp_cyl( cyl_lbetagrid, r, z, phi )
    betaout = exp(betaout)
 endif

 if( cyl_npgrid( rbin, zbin, phibin )==0) then
    beta_by_interp_cyl = betain
 else
    beta_by_interp_cyl = (betaout + (betain-betaout)/cyl_npgrid( rbin, zbin, phibin ))
 endif

end function beta_by_interp_cyl

!----------------------------------------------
!+
!  Finds a value by bilinearly interpolating on a grid
!+
!----------------------------------------------

real function bilin_interp( array, r, theta, phi )
 use physcon,    only: twopi
 real, intent(in) :: array(0:nr-1, 0:nth-1, 0:nph-1), phi, theta, r
 real :: t1, t2, p1, p2, ft1p1, ft2p1,ft1p2, ft2p2, tmid, pmid
 real :: ftp, dta, dtb, dpa, dpb, dummy
 integer :: tbin1, tbin2, pbin1, pbin2, tbin0, pbin0, rbin0

 bilin_interp = 0.
 tbin0=0
 if( r>rmin) then
    call get_grid_bins( r, theta, rbin0, tbin0, phi, pbin0 )
 else
    call get_grid_bins( rmin+eps, theta, rbin0, tbin0, phi, pbin0 )
    rbin0=0
 endif
 call get_grid_points( thmingrid, tbin0, nth, thetamax, t1, t2, tmid )

 if( tbin0 == 0.and.theta<tmid ) then
    call get_bracket_grid_points( thmingrid, 0, nth, thetamax, t1, t2 )
    tbin1 = 0
 else if( tbin0 == nth-1.and.theta>tmid ) then
    call get_bracket_grid_points( thmingrid, nth-2, nth, thetamax, t1, t2 )
    tbin1 = nth-2
 else if( theta > tmid ) then
    call get_bracket_grid_points( thmingrid, tbin0, nth, thetamax, t1, t2 )
    tbin1 = tbin0
 else
    call get_bracket_grid_points( thmingrid, tbin0-1, nth, thetamax, t1, t2 )
    tbin1 = tbin0-1
 endif

 call get_grid_points( phmingrid, pbin0, nph, phimax, p1, p2, pmid )

 if( pbin0 == 0.and.phi<pmid ) then      ! These ones "wrap around"
    call get_bracket_grid_points( phmingrid, 0, nph, phimax, dummy, p2 )
    call get_bracket_grid_points( phmingrid, nph-2, nph, phimax, p1, dummy )
    p1 = p1 - twopi
    pbin1 = nph-1
    pbin2 = 0
 else if( pbin0 == nph-1.and.phi>pmid ) then
    call get_bracket_grid_points( phmingrid, nph-2, nph, phimax, p1, dummy )
    call get_bracket_grid_points( phmingrid, 0, nph, phimax, dummy, p2 )
    p2 = p2 + twopi
    pbin1 = nph-1
    pbin2 = 0
 else if( phi>pmid ) then
    call get_bracket_grid_points( phmingrid, pbin0, nph, phimax, p1, p2 )
    pbin1 = pbin0
    pbin2 = pbin0+1
 else
    call get_bracket_grid_points( phmingrid, pbin0-1, nph, phimax, p1, p2 )
    pbin1 = pbin0-1
    pbin2 = pbin0
 endif

 tbin2=tbin1+1

 if( pbin1<0 ) pbin1 = nph-1
 if( pbin2>nph-1 ) pbin2 = 0

 ft1p1 = array( rbin0, tbin1, pbin1 )
 ft1p2 = array( rbin0, tbin1, pbin2 )
 ft2p1 = array( rbin0, tbin2, pbin1 )
 ft2p2 = array( rbin0, tbin2, pbin2 )

 dta = theta - t1
 dtb = t2 - theta
 dpa = phi - p1
 dpb = p2 - phi

 ftp =    ft1p1*dtb*dpb &
          + ft2p1*dta*dpb &
          + ft1p2*dtb*dpa &
          + ft2p2*dta*dpa

 bilin_interp = ftp/( (t2-t1)*(p2-p1) )

end function bilin_interp

!----------------------------------------------
!+
!  Finds a value by bilinearly interpolating on a cylindrical grid
!+
!----------------------------------------------

real function bilin_interp_cyl( array, r, z, phi )
 use physcon,    only: twopi
 real, intent(in) :: array(0:nr-1, 0:nz, 0:nph-1), phi, z, r
 real :: r1, r2, p1, p2, fr1p1, fr2p1,fr1p2, fr2p2, rmid, pmid
 real :: frp, dra, drb, dpa, dpb, dummy
 integer :: rbin1, rbin2, pbin1, pbin2, zbin0, pbin0, rbin0
 real    :: znew
 znew=z
 zbin0=0
 if( z>zmax ) znew = zmax - eps
 if( z<zmin ) znew = zmin + eps
 call get_grid_bins( r, znew, rbin0, zbin0, phi, pbin0 )
 call get_grid_points( cyl_ringrid, rbin0, nr, rmax, r1, r2, rmid )

 if( rbin0 == 0.and.r<rmid ) then
    call get_bracket_grid_points( cyl_ringrid, 0, nr, rmax, r1, r2 )
    rbin1 = 0
 else if( rbin0 == nr-1.and.r>rmid ) then
    call get_bracket_grid_points( cyl_ringrid, nr-2, nr, rmax, r1, r2 )
    rbin1 = nr-2
 else if( r > rmid ) then
    call get_bracket_grid_points( cyl_ringrid, rbin0, nr, rmax, r1, r2 )
    rbin1 = rbin0
 else
    call get_bracket_grid_points( cyl_ringrid, rbin0-1, nr, rmax, r1, r2 )
    rbin1 = rbin0-1
 endif

 call get_grid_points( phmingrid, pbin0, nph, phimax, p1, p2, pmid )

 if( pbin0 == 0.and.phi<pmid ) then      ! These ones "wrap around"
    call get_bracket_grid_points( phmingrid, 0, nph, phimax, dummy, p2 )
    call get_bracket_grid_points( phmingrid, nph-2, nph, phimax, p1, dummy )
    p1 = p1 - twopi
    pbin1 = nph-1
    pbin2 = 0
 else if( pbin0 == nph-1.and.phi>pmid ) then
    call get_bracket_grid_points( phmingrid, nph-2, nph, phimax, p1, dummy )
    call get_bracket_grid_points( phmingrid, 0, nph, phimax, dummy, p2 )
    p2 = p2 + twopi
    pbin1 = nph-1
    pbin2 = 0
 else if( phi>pmid ) then
    call get_bracket_grid_points( phmingrid, pbin0, nph, phimax, p1, p2 )
    pbin1 = pbin0
    pbin2 = pbin0+1
 else
    call get_bracket_grid_points( phmingrid, pbin0-1, nph, phimax, p1, p2 )
    pbin1 = pbin0-1
    pbin2 = pbin0
 endif

 rbin2=rbin1+1

 if( pbin1<0 ) pbin1 = nph-1
 if( pbin2>nph-1 ) pbin2 = 0

 fr1p1 = array( rbin1, zbin0, pbin1 )
 fr1p2 = array( rbin1, zbin0, pbin2 )
 fr2p1 = array( rbin2, zbin0, pbin1 )
 fr2p2 = array( rbin2, zbin0, pbin2 )

 dra = r - r1
 drb = r2 - r
 dpa = phi - p1
 dpb = p2 - phi

 frp =    fr1p1*drb*dpb &
          + fr2p1*dra*dpb &
          + fr1p2*drb*dpa &
          + fr2p2*dra*dpa

 bilin_interp_cyl = frp/( (r2-r1)*(p2-p1) )

end function bilin_interp_cyl

!----------------------------------------------
!+
!  Returns the natural logarithm of a number,
!  or a very large negative number if given a negative
!+
!----------------------------------------------

real function careful_log( x )
 real, intent(in) :: x
 real, parameter  :: xbig = range(x)*log(10.)
 if( x <= 0. ) then
    careful_log = -100.
 else
    careful_log = max(log(x), -xbig/2.)
 endif
end function careful_log

!----------------------------------------------
!+
!  function computing the disc scale height
!+
!----------------------------------------------
real function calc_scaleheight( r )
 use eos, only:polyk, qfacdisc
 real, intent(in) :: r
 real :: omega, cs
 if( r > 0. ) then
    omega = 1.0/r**(1.5)
    cs = sqrt(polyk) * r**(-qfacdisc)
    calc_scaleheight = (cs/omega)
 else
    calc_scaleheight = -1000.       !sentinel value for star interior
 endif
end function calc_scaleheight

!----------------------------------------------
!+
!  function computing the disc surface density
!+
!----------------------------------------------
real function calc_sigma( r )
 use units,   only:umass
 use physcon, only:solarm
 real, intent(in) :: r
 real :: R_in = 1., Mdisc
 Mdisc = 1.4d0*solarm/umass*5.e-16
 if( r>r_In ) then
    calc_sigma = sqrt(R_in) * Mdisc * r**(-3./2.)*(1-sqrt(R_in/r))
 else
    calc_sigma = 0.
 endif

end function calc_sigma

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_lumin_nsdisc(iunit)
 use infile_utils,         only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to the neutron star disc'
 call write_inopt(BurstProfile,'BurstProfile',&
                  'Burst Profile',iunit)
 ! Between 0 and 1    = constant luminosity
 ! Any negative value = burst profile as described in set_Lstar

 call write_inopt(frac_X,'frac_X',&
                  'Hydrogen fraction (-ve for zero opacity)',iunit)

 call write_inopt(rad_trans, 'rad_trans', &
                  'Radiation transport prescription', iunit)

end subroutine write_options_lumin_nsdisc

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_lumin_nsdisc(name,valstring,imatch,igotall,ierr)
 use io,                   only:fatal, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_lumin_nsdisc'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('BurstProfile')
    read(valstring,*,iostat=ierr) BurstProfile
 case('frac_X')
    read(valstring,*,iostat=ierr) frac_X
 case('rad_trans')
    read(valstring,*,iostat=ierr) rad_trans
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)

end subroutine read_options_lumin_nsdisc

end module lumin_nsdisc
