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
! :Owner: Spencer Magnall
!
! :Runtime parameters:
!   - metric_file : *file from which to read tabulated metric (if present)*
!
! :Dependencies: dump_utils, infile_utils, inverse4x4, metric_et_utils
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
 real, intent(in) :: time

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
 use dump_utils, only:dump_h
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time,accradius
 integer,      intent(out)   :: ierr

 ierr = 0

end subroutine write_headeropts_metric

!-----------------------------------------------------------------------
!+
!  reads relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_metric(hdr,ierr)
 use dump_utils, only:dump_h
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr

 ierr  = 0

end subroutine read_headeropts_metric

!-----------------------------------------------------------------------
!+
!  writes metric options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_metric(iunit)
 use infile_utils,    only:write_inopt
 use metric_et_utils, only:metric_file
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to the '//trim(metric_type)//' metric'
 call write_inopt(metric_file,'metric_file','file from which to read tabulated metric (if present)',iunit)

end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(db,nerr)
 use infile_utils,    only:inopts,read_inopt
 use metric_et_utils, only:metric_file
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(metric_file,'metric_file',db,errcount=nerr)

end subroutine read_options_metric

end module metric

