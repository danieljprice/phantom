!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine comparing time between dumps
!
!  Produces three output files:
!     radial.out      - rho, kappa, beta calculated on the radial grid
!     radialinterp.out- the same quantities interpolated onto a Cartesian grid
!     applied.out     - beta for all the actual particles
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: lumin_nsdisc, physcon, units
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'prdrag'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use physcon
 use lumin_nsdisc, only: beta, nr, nth, rmin, rmax, thetamin, thetamax, &
                         make_beta_grids, densitygrid, tauradgrid, betagrid, &
                         get_grid_points, ringrid, thmingrid, bilin_interp, careful_log, &
                         get_grid_bins, sphere_segment, Lstar, lbetagrid, ltauradgrid, ldensitygrid, &
                         get_bracket_grid_points
 use units, only: udist, umass, utime

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 integer :: ierr
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass, time
 integer :: rbin, tbin, ri, zi, ipart
 real :: r, theta, rcyl, z, rho, tau, betac, cell_rin, cell_rout, cell_thin, cell_thout, cell_rmid, cell_thmid
 real :: r1, r2

 print*,' Hello Hauke, time in file = ',time
 open( unit=106, file='radial.out', status='replace',  iostat=ierr)
 if ( ierr /= 0 ) stop 'error opening radial.out'

 call make_beta_grids( xyzh, particlemass, npart )

 write(106,*), "#r_rad rbin theta thetabin r_cyl z ln(rho) ln(tau) ln(beta)"

 call get_bracket_grid_points( ringrid, 50, nr, rmax, r1, r2 )

 do rbin=0, nr-1
    call get_grid_points( ringrid, rbin, nr, rmax, cell_rin, cell_rout, r )
    do tbin=0, nth-1
       call get_grid_points( thmingrid, tbin, nth, thetamax, cell_thin, cell_thout, theta )
       write( 106,* ), r, rbin, theta, tbin, r*sin(theta), r*cos(theta), ldensitygrid( rbin, tbin ), &
                    ltauradgrid( rbin, tbin), lbetagrid( rbin, tbin )
    enddo
 enddo

 close( 106 )
 open( unit=106, file='radialinterp.out', status='replace',  iostat=ierr)
 if ( ierr /= 0 ) stop 'error opening radialinterp.out'

 write(106,*), "#r_rad rbin theta thetabin r_cyl z rho tau beta"
 Lstar=1
 do ri=0, 600
    rcyl = rmin + ri/1.
    do zi=0, 250

       z = zi/1.
       r = sqrt( rcyl**2 + z**2 )
       theta = atan(rcyl/z)
       call get_grid_bins( r, theta, rbin, tbin )
       call get_grid_points( ringrid, rbin, nr, rmax, cell_rin, cell_rout, cell_rmid )
       call get_grid_points( thmingrid, tbin, nth, thetamax, cell_thin, cell_thout, cell_thmid )
       rho  = exp(bilin_interp( ldensitygrid, r, theta ))
       tau  = exp(bilin_interp( ltauradgrid, r, theta ))
       betac = exp(bilin_interp( lbetagrid, r, theta ))
       write( 106, * ), r, rbin, theta, tbin, rcyl, z, rho, tau, betac

    enddo

 enddo

 close(106)

 open( unit=106, file='applied.out', status='replace',  iostat=ierr)
 if ( ierr /= 0 ) stop 'error opening applied.out'

 write(106,*), "#x y z r_cyl beta r_bin th_bin"

 do ipart=1, npart


    r = sqrt(xyzh(1,ipart)**2 + xyzh(2,ipart)**2 + xyzh(3,ipart)**2)
    if ( r>1 ) then

       call get_grid_bins( r, acos(abs(xyzh(3,ipart))/r), rbin, tbin )
       write(106,*), xyzh(1,ipart), xyzh(2,ipart), xyzh(3,ipart), r, &
        beta(xyzh(1,ipart), xyzh(2,ipart), xyzh(3,ipart),0.0), rbin, tbin
    endif
 enddo

 close(106)



end subroutine do_analysis

end module
