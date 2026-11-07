!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric
!
! Generic module for a tabulated metric, e.g. from the Einstein Toolkit
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - evolve_metric   : *evolve the metric with the NR code (else read metric_file)*
!   - metric_file     : *file from which to read tabulated metric (if present)*
!   - ngrid_evolve    : *grid nodes per dimension for the evolved metric*
!   - nr_deposit      : *T_munu deposit (0=SPH, 1=SPH at nodes, 2=zero moment, 3=SPH exact)*
!   - nr_deposit_hmin : *minimum h for the deposit at the nodes (grid spacings)*
!   - nr_dust_source  : *deposit rho*, rho* u_i once per advance (pressureless matter only)*
!   - nr_refine_h     : *refine the NR grid where h < nr_refine_h dx (with amr_maxlevel > 0)*
!   - nr_solve_init   : *solve the full constraints for the initial data (else linear)*
!   - redshift_init   : *redshift at which the scale factor is 1 (cosmology, <0 = none)*
!
! :Dependencies: dump_utils, infile_utils, inverse4x4, metric_et_utils,
!   nr_api
!
 implicit none
 character(len=*), parameter :: metric_type = 'et'
 integer,          parameter :: imetric     = 7
 ! This are dummy parameters to stop the compiler complaing
 ! Not used anywhere in the code - Needs a fix!
 real, public  :: mass1 = 1.       ! mass of central object
 real, public  :: a     = 0.0       ! spin of central object
 real, public  :: charge= 0.       ! charge of central object

contains
!----------------------------------------------------------------
!+
!  Compute the metric tensor in both covariant (gcov) and
!  contravariant (gcon) form. Here we merely interpolate
!  these values from the global grid.
!+
!----------------------------------------------------------------
pure subroutine get_metric_cartesian(position,gcov,gcon,sqrtg)
 use metric_et_utils, only:gridinit,interpolate_gcov
 use inverse4x4,      only:inv4x4
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: gconi(0:3,0:3),det

 ! The subroutine that computes the metric tensor for a given position
 ! In this case it is interpolated from the global grid values
 if (.not.gridinit) then
    ! This is required for phantomsetup
    ! As no grid information has been passed to phantom from ET
    ! So interpolation cannot be performed
    gcov = 0.
    gcov(0,0) = -1.
    gcov(1,1) = 1.
    gcov(2,2) = 1.
    gcov(3,3) = 1.
    if (present(gcon)) gcon = gcov
    if (present(sqrtg)) sqrtg = 1.
 else
    call interpolate_gcov(position,gcov)
    ! inverse and determinant are computed from the interpolated metric,
    ! so that they are consistent with it
    if (present(gcon) .or. present(sqrtg)) then
       call inv4x4(gcov,gconi,det)
       if (present(gcon)) gcon = gconi
       if (present(sqrtg)) sqrtg = sqrt(abs(det))
    endif
 endif

end subroutine get_metric_cartesian

!-----------------------------------------------------------------------
!+
!  dummy routine to get the metric in spherical coordinates (not used)
!+
!-----------------------------------------------------------------------
pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: r2,sintheta

 gcov = 0.

 r2       = position(1)**2
 sintheta = sin(position(2))

 gcov(0,0) = -1.
 gcov(1,1) = 1.
 gcov(2,2) = r2
 gcov(3,3) = r2*sintheta**2

 if (present(gcon)) then
    gcon      = 0.
    gcon(0,0) = -1.
    gcon(1,1) = 1.
    gcon(2,2) = 1./r2
    gcov(3,3) = 1./gcov(3,3)
 endif

 if (present(sqrtg)) sqrtg = r2*sintheta

end subroutine get_metric_spherical

!-----------------------------------------------------------------------
!+
!  cartesian metric derivatives, interpolates the derivatives from
!  the grid
!+
!-----------------------------------------------------------------------
pure subroutine metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
 use metric_et_utils, only:gridinit,interpolate_gcov_derivs
 real, intent(in)  :: position(3)
 real, intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
 real :: dgcov(0:3,0:3,3)

 if (.not.gridinit) then
    dgcov = 0.
 else
    call interpolate_gcov_derivs(position,dgcov)
 endif
 dgcovdx = dgcov(:,:,1)
 dgcovdy = dgcov(:,:,2)
 dgcovdz = dgcov(:,:,3)

end subroutine metric_cartesian_derivatives

!-----------------------------------------------------------------------
!+
!  dummy routine for spherical metric derivatives, not used
!+
!-----------------------------------------------------------------------
pure subroutine metric_spherical_derivatives(position,dgcovdr, dgcovdtheta, dgcovdphi)
 real, intent(in)  :: position(3)
 real, intent(out) :: dgcovdr(0:3,0:3),dgcovdtheta(0:3,0:3),dgcovdphi(0:3,0:3)
 real :: r, theta

 r     = position(1)
 theta = position(2)

 dgcovdr     = 0.
 dgcovdtheta = 0.
 dgcovdphi   = 0.

 dgcovdr(2,2) = 2*r
 dgcovdr(3,3) = 2*r*sin(theta)**2

 dgcovdtheta(3,3) = 2*r**2*cos(theta)*sin(theta)

end subroutine metric_spherical_derivatives

!-----------------------------------------------------------------------
!+
!  dummy routine to convert cartesian to spherical coordinates
!+
!-----------------------------------------------------------------------
pure subroutine cartesian2spherical(xcart,xspher)
 real, intent(in)  :: xcart(3)
 real, intent(out) :: xspher(3)
 real :: x,y,z
 real :: r,theta,phi

 x  = xcart(1)
 y  = xcart(2)
 z  = xcart(3)

 r     = sqrt(x**2+y**2+z**2)
 theta = acos(z/r)
 phi   = atan2(y,x)

 xspher   = (/r,theta,phi/)

end subroutine cartesian2spherical

!-------------------------------------------------------------------------------
!+
!  Subroutine to update the metric inputs if time dependent
!+
!-------------------------------------------------------------------------------
subroutine update_metric(time)
 use metric_et_utils, only:advance_metric
 real, intent(in) :: time

 ! advance the metric if it is evolved by a numerical relativity code
 if (associated(advance_metric)) call advance_metric(time)

end subroutine update_metric

!-----------------------------------------------------------------------
!+
!  Check if a particle should be accreted by the black hole
!+
!-----------------------------------------------------------------------
subroutine accrete_particles_metric(xi,yi,zi,mi,ti,accradius,accreted)
 real,    intent(in)  :: xi,yi,zi,mi,ti,accradius
 logical, intent(out) :: accreted

 accreted = .false.

end subroutine accrete_particles_metric

!-----------------------------------------------------------------------
!+
!  writes relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_metric(hdr,time,accradius,ierr)
 use dump_utils,      only:dump_h,add_to_rheader
 use metric_et_utils, only:evolve_metric,proper_time,scale_factor,redshift_init
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time,accradius
 integer,      intent(out)   :: ierr

 ierr = 0
 ! for an evolved metric, the proper time (age) and mean scale factor of the
 ! grid, and the redshift if the run is cosmological
 if (evolve_metric .and. proper_time >= 0.) then
    call add_to_rheader(proper_time,'proper_time',hdr,ierr)
    call add_to_rheader(scale_factor,'scale_factor',hdr,ierr)
    if (redshift_init >= 0.) call add_to_rheader((1. + redshift_init)/scale_factor - 1.,'redshift',hdr,ierr)
 endif

end subroutine write_headeropts_metric

!-----------------------------------------------------------------------
!+
!  reads relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_metric(hdr,ierr)
 use dump_utils,      only:dump_h,extract
 use metric_et_utils, only:proper_time
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 integer :: ierrs

 ierr  = 0
 ! proper time of an evolved metric, so that it continues on restart
 call extract('proper_time',proper_time,hdr,ierrs)
 if (ierrs /= 0) proper_time = -1.

end subroutine read_headeropts_metric

!-----------------------------------------------------------------------
!+
!  writes metric options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_metric(iunit)
 use infile_utils,    only:write_inopt
 use metric_et_utils, only:metric_file,evolve_metric,ngrid_evolve,nr_solve_init,nr_deposit,nr_deposit_hmin,&
                            redshift_init,nr_dust_source,nr_refine_h
 use nr_api,          only:nr_linked,nr_write_options
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to the '//trim(metric_type)//' metric'
 call write_inopt(metric_file,'metric_file','file from which to read tabulated metric (if present)',iunit)
 if (nr_linked()) then
    call write_inopt(evolve_metric,'evolve_metric','evolve the metric with the NR code (else read metric_file)',iunit)
    call write_inopt(ngrid_evolve,'ngrid_evolve','grid nodes per dimension for the evolved metric',iunit)
    call write_inopt(nr_solve_init,'nr_solve_init','solve the full constraints for the initial data (else linear)',iunit)
    call write_inopt(nr_deposit,'nr_deposit','T_munu deposit (0=SPH, 1=SPH at nodes, 2=zero moment, 3=SPH exact)',iunit)
    call write_inopt(nr_deposit_hmin,'nr_deposit_hmin','minimum h for the deposit at the nodes (grid spacings)',iunit)
    call write_inopt(redshift_init,'redshift_init','redshift at which the scale factor is 1 (cosmology, <0 = none)',iunit)
    call write_inopt(nr_dust_source,'nr_dust_source','deposit rho*, rho* u_i once per advance (pressureless matter only)',iunit)
    call write_inopt(nr_refine_h,'nr_refine_h','refine the NR grid where h < nr_refine_h dx (with amr_maxlevel > 0)',iunit)
    call nr_write_options(iunit)
 endif

end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(db,nerr)
 use infile_utils,    only:inopts,read_inopt
 use metric_et_utils, only:metric_file,evolve_metric,ngrid_evolve,nr_solve_init,nr_deposit,nr_deposit_hmin,&
                            redshift_init,nr_dust_source,nr_refine_h
 use nr_api,          only:nr_linked,nr_read_options
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(metric_file,'metric_file',db,errcount=nerr)
 if (nr_linked()) then
    call read_inopt(evolve_metric,'evolve_metric',db,errcount=nerr)
    call read_inopt(ngrid_evolve,'ngrid_evolve',db,min=5,errcount=nerr)
    call read_inopt(nr_solve_init,'nr_solve_init',db,errcount=nerr)
    call read_inopt(nr_deposit,'nr_deposit',db,min=0,max=3,errcount=nerr)
    call read_inopt(nr_deposit_hmin,'nr_deposit_hmin',db,min=0.,errcount=nerr)
    call read_inopt(redshift_init,'redshift_init',db,errcount=nerr)
    call read_inopt(nr_dust_source,'nr_dust_source',db,errcount=nerr)
    call read_inopt(nr_refine_h,'nr_refine_h',db,min=0.,default=1.,errcount=nerr)
    call nr_read_options(db,nerr)
 endif

end subroutine read_options_metric

end module metric

