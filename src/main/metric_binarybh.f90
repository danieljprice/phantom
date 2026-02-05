!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric
!
! Superposed Kerr-Schild metric for a binary black hole system
!
! :References: 
!   Combi & Ressler (2024) arXiv:2403.13308
!   Original file from: https://zenodo.org/records/10841021
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, timestep
!
 use iso_c_binding, only:c_double
 implicit none
 character(len=*), parameter :: metric_type = 'binarybh'
 integer,          parameter :: imetric = 5

 ! extern void SuperposedBBH(const double *xx, double gcov[][NDIM], const double *traj_array)
 interface
    pure subroutine SuperposedBBH(xx,gcov,traj_array) bind(c,name='SuperposedBBH')
     import c_double
     real(c_double), intent(in)  :: xx(3)
     real(c_double), intent(out) :: gcov(4,4)
     real(c_double), intent(in)  :: traj_array(20)
    end subroutine SuperposedBBH
 end interface

 integer, parameter :: nparams = 20
 real, public :: bh_trajectory(nparams)

 real, private :: bh1_spinx,bh1_spiny, bh1_spinz, bh2_spinx, bh2_spiny, bh2_spinz
 real, public :: mass1 = 1.0
 real, public :: mass2 = 0.
 real, public :: a = 0.    ! black hole 1 spin
 real, public :: a2 = 0.   ! black hole 2 spin
 character(len=128), public :: trajectory_file = 'cbwaves.txt'

contains

!-------------------------------------------------------------------------------
!+
!  Subroutine to update the metric inputs if time dependent
!+
!-------------------------------------------------------------------------------
subroutine update_metric(time)
 real, intent(in) :: time
 real :: x1(3),x2(3),v1(3),v2(3)

 x1 = [0.,0.,0.]
 x2 = [10000000.,0.,0.]
 v1 = [0.,0.,0.]
 v2 = [0.,0.,0.]
 mass1 = 1.
 mass2 = 0.
 !call get_trajectory_from_file(time,x1,x2,v1,v2)

 bh_trajectory = 0.  
 bh_trajectory(1:3) = x1
 bh_trajectory(4:6) = x2
 bh_trajectory(7:9) = v1
 bh_trajectory(10:12) = v2
 bh_trajectory(13) = bh1_spinx
 bh_trajectory(14) = bh1_spiny
 bh_trajectory(15) = bh1_spinz
 bh_trajectory(16) = bh2_spinx
 bh_trajectory(17) = bh2_spiny
 bh_trajectory(18) = bh2_spinz
 bh_trajectory(19) = mass1
 bh_trajectory(20) = mass2

end subroutine update_metric

!----------------------------------------------------------------
!+
!  Read the binary black hole trajectory from file with linear
!  interpolation in time to match the specified time.
!+
!----------------------------------------------------------------
subroutine get_trajectory_from_file(time,x1,x2,v1,v2)
 use io, only:error
 real, intent(in) :: time
 real, intent(out) :: x1(3),x2(3),v1(3),v2(3)
 integer :: iu,ierr,nlines
 real :: t_prev,t_next
 real :: x1_prev(3),x2_prev(3),v1_prev(3),v2_prev(3)
 real :: x1_next(3),x2_next(3),v1_next(3),v2_next(3)
 real :: frac

 x1 = [0.,0.,0.]
 x2 = [0.,0.,0.]
 v1 = [0.,0.,0.]
 v2 = [0.,0.,0.]

 open(newunit=iu,file=trajectory_file,status='old',action='read',iostat=ierr)
 if (ierr /= 0) then
    call error('metric_binary_bh','could not open trajectory file '//trim(trajectory_file))
 endif
 read(iu,*,iostat=ierr) nlines
 if (ierr /= 0) then
    call error('metric_binary_bh','could not read nlines from '//trim(trajectory_file))
 endif
 read(iu,*,iostat=ierr) t_prev,x1_prev,x2_prev,v1_prev,v2_prev
 if (ierr /= 0) then
    call error('metric_binary_bh','could not read first trajectory line from '//trim(trajectory_file))
 endif
 if (time <= t_prev) then
    x1 = x1_prev
    x2 = x2_prev
    v1 = v1_prev
    v2 = v2_prev
    close(iu)
    return
 endif
 do
    read(iu,*,iostat=ierr) t_next,x1_next,x2_next,v1_next,v2_next
    if (ierr /= 0) then
       call error('metric_binary_bh','end of file reached before time '//trim(trajectory_file))
    endif
    if (time <= t_next) then
       frac = (time - t_prev) / (t_next - t_prev)
       x1 = x1_prev + frac * (x1_next - x1_prev)
       x2 = x2_prev + frac * (x2_next - x2_prev)
       v1 = v1_prev + frac * (v1_next - v1_prev)
       v2 = v2_prev + frac * (v2_next - v2_prev)
       close(iu)
       return
    endif
    t_prev = t_next
    x1_prev = x1_next
    x2_prev = x2_next
    v1_prev = v1_next
    v2_prev = v2_next
 end do
 print*,' time = ',time,' binary separation = ',sqrt(dot_product(x1 - x2,x1 - x2))

end subroutine get_trajectory_from_file

!----------------------------------------------------------------
!+
!  Compute the metric tensor in both covariant (gcov) and
!  contravariant (gcon) form
!+
!----------------------------------------------------------------
pure subroutine get_metric_cartesian(xx,gcov,gcon,sqrtg)
 use inverse4x4, only:inv4x4
 real, intent(in)  :: xx(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: det

 call SuperposedBBH(xx,gcov,bh_trajectory)

 if (present(gcon)) then
    gcon = 0.
    call inv4x4(gcov,gcon,det)
    if (present(sqrtg)) sqrtg = sqrt(-det)
 endif

end subroutine get_metric_cartesian

pure subroutine metric_cartesian_derivatives(xx,dgcovdx,dgcovdy,dgcovdz)
 real,    intent(in)  :: xx(3)
 real,    intent(out) :: dgcovdx(0:3,0:3),dgcovdy(0:3,0:3),dgcovdz(0:3,0:3)
 real(c_double), parameter :: eps = 1e-10_c_double
 real(c_double), parameter :: two_eps = 2._c_double * eps
 real(c_double) :: gcov1(0:3,0:3), gcov2(0:3,0:3)

 dgcovdx = 0.
 dgcovdy = 0.
 dgcovdz = 0.

 ! metric has no time dependence; centred differencing: (f(x+h)-f(x-h))/(2h)
 ! x direction
 call get_metric_cartesian(xx+[eps,0.,0.],gcov1)
 call get_metric_cartesian(xx+[-eps,0.,0.],gcov2)
 dgcovdx = (gcov1 - gcov2) / two_eps

 ! y direction
 call get_metric_cartesian(xx+[0.,eps,0.],gcov1)
 call get_metric_cartesian(xx+[0.,-eps,0.],gcov2)
 dgcovdy = (gcov1 - gcov2) / two_eps

 ! z direction
 call get_metric_cartesian(xx+[0.,0.,eps],gcov1)
 call get_metric_cartesian(xx+[0.,0.,-eps],gcov2)
 dgcovdz = (gcov1 - gcov2) / two_eps

end subroutine metric_cartesian_derivatives


!----------------------------------------------------------------
!+
!  The metric tensor in SPHERICAL-like form
!  (these are dummy routines for compatibility with other metrics)
!+
!----------------------------------------------------------------
pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg

 gcov = 0.
 if (present(gcon)) gcon = 0.
 if (present(sqrtg)) sqrtg = 0.

end subroutine get_metric_spherical

pure subroutine cartesian2spherical(xcart,xspher)
 real, intent(in)  :: xcart(3)
 real, intent(out) :: xspher(3)

 xspher = xcart

end subroutine cartesian2spherical

pure subroutine metric_spherical_derivatives(position,dgcovdr,dgcovdtheta,dgcovdphi)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi

 dgcovdr     = 0.
 dgcovdtheta = 0.
 dgcovdphi   = 0. 

end subroutine metric_spherical_derivatives

!----------------------------------------------------------------
!+
!  Check if a particle should be accreted by either black hole
!+
!----------------------------------------------------------------
subroutine accrete_particles_metric(xi,yi,zi,mi,ti,accradius1,accradius2,accreted)
 real,    intent(in)  :: xi,yi,zi,mi,ti,accradius1,accradius2
 logical, intent(out) :: accreted
 real :: r1,r2,x1(3),x2(3)

 x1 = bh_trajectory(1:3)
 x2 = bh_trajectory(4:6)

 r1 = ((xi-x1(1))**2 + (yi-x1(2))**2 + (zi-x1(3))**2)/mass1**2
 r2 = ((xi-x2(1))**2 + (yi-x2(2))**2 + (zi-x2(3))**2)/mass2**2
 if (r1 < accradius1**2 .or. r2 < accradius2**2) then
    accreted = .true.
 else
    accreted = .false.
 endif

end subroutine accrete_particles_metric

!-----------------------------------------------------------------------
!+
!  writes relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_binarybh(hdr,time,accradius1,accradius2,ierr)
 use dump_utils, only:lentag,dump_h,add_to_rheader
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time
 real,         intent(in)    :: accradius1,accradius2
 integer,      intent(out)   :: ierr
 character(len=lentag)       :: tags(16)
 real    :: rheader(16)
 integer :: i

 ierr = 0
 rheader(1:3) = bh_trajectory(1:3)
 rheader(4) = mass1
 rheader(5) = accradius1
 rheader(6:8) = bh_trajectory(4:6)
 rheader(9) = mass2
 rheader(10) = accradius2
 rheader(11:13) = bh_trajectory(7:9)
 rheader(14:16) = bh_trajectory(10:12)

 !  rheader(17) = accretedmass1
 !  rheader(18) = accretedmass2

 tags(1:16) = (/'x1 ','y1 ','z1 ','m1 ','h1 ','x2 ','y2 ','z2 ','m2 ', &
                'h2 ','vx1','vy1','vz1','vx2','vy2','vz2'/)

 do i=1,16
    call add_to_rheader(rheader(i),tags(i),hdr,ierr)
 enddo

end subroutine write_headeropts_binarybh

!-----------------------------------------------------------------------
!+
!  reads relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_binarybh(hdr,ierr)
 use dump_utils, only:dump_h,extract
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr

 ierr  = 0
! call extract('accretedmass1',accretedmass1,hdr,ierr1)
! call extract('accretedmass2',accretedmass2,hdr,ierr2)

! if (ierr1 /= 0 .or. ierr2 /= 0) then
 !   write(*,*) ' ERROR extracting accretedmass1 and accretedmass2 from file'
 !   ierr = 1
 !endif

end subroutine read_headeropts_binarybh

!-----------------------------------------------------------------------
!+
!  writes metric options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_metric(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# Options relating to the '//trim(metric_type)//' metric'
 call write_inopt(bh1_spinx,'bh1_spinx','spin parameter for black hole 1',iunit)
 call write_inopt(bh1_spiny,'bh1_spiny','spin parameter for black hole 1',iunit)
 call write_inopt(bh1_spinz,'bh1_spinz','spin parameter for black hole 1',iunit)
 call write_inopt(bh2_spinx,'bh2_spinx','spin parameter for black hole 2',iunit)
 call write_inopt(bh2_spiny,'bh2_spiny','spin parameter for black hole 2',iunit)
 call write_inopt(bh2_spinz,'bh2_spinz','spin parameter for black hole 2',iunit)
 call write_inopt(mass1,'mass1','mass of black hole 1',iunit)
 call write_inopt(mass2,'mass2','mass of black hole 2',iunit)
 call write_inopt(trajectory_file,'trajectory_file','file containing binary black hole trajectory',iunit)

end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(db,nerr)
 use infile_utils, only:inopts,read_inopt
 use io,           only:error
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(bh1_spinx,'bh1_spinx',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(bh1_spiny,'bh1_spiny',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(bh1_spinz,'bh1_spinz',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(bh2_spinx,'bh2_spinx',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(bh2_spiny,'bh2_spiny',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(bh2_spinz,'bh2_spinz',db,errcount=nerr,min=-1.,max=1.)
 call read_inopt(mass1,'mass1',db,errcount=nerr,min=0.,max=1.e12)
 call read_inopt(mass2,'mass2',db,errcount=nerr,min=0.,max=1.e12)
 call read_inopt(trajectory_file,'trajectory_file',db,errcount=nerr)

 a = sqrt(bh1_spinx**2 + bh1_spiny**2 + bh1_spinz**2)
 a2 = sqrt(bh2_spinx**2 + bh2_spiny**2 + bh2_spinz**2)
 if (a > 1.) then
    call error('metric','black hole spin: a > 1 for black hole 1')
    nerr = nerr + 1
 endif
 if (a2 > 1.) then
    call error('metric','black hole spin: a > 1 for black hole 2')
    nerr = nerr + 1
 endif

end subroutine read_options_metric

end module metric
