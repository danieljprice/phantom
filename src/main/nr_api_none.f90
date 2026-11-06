!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module nr_api
!
! Stub of the interface to the numerical relativity library (libnr), used
! when Phantom is compiled without NR=yes. The real module is supplied by
! the library (nr_api.mod in $(NRDIR)/lib).
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils
!
 use infile_utils, only:inopts
 implicit none

 integer, parameter, public :: nr_dp = kind(1.d0)

 abstract interface
  subroutine tmunu_func(t,gcov,metricderivs,tmunu)
   import nr_dp
   real(nr_dp), intent(in)  :: t,gcov(0:,0:,:,:,:),metricderivs(0:,0:,:,:,:,:)
   real(nr_dp), intent(out) :: tmunu(0:,0:,:,:,:)
  end subroutine tmunu_func
 end interface

 public :: tmunu_func
 public :: nr_linked,nr_init,nr_evolve_to,nr_get_metric,nr_time,nr_constraints
 public :: nr_formulation,nr_finalise,nr_dt_expansion
 public :: nr_write_options,nr_read_options

 private

contains

logical function nr_linked()
 nr_linked = .false.
end function nr_linked

subroutine nr_init(n,origin,dx,t0,alpha,beta,g3,kij,get_tmunu,ierr)
 integer,     intent(in)  :: n(3)
 real(nr_dp), intent(in)  :: origin(3),dx(3),t0
 real(nr_dp), intent(in)  :: alpha(:,:,:),beta(:,:,:,:),g3(:,:,:,:,:),kij(:,:,:,:,:)
 procedure(tmunu_func)    :: get_tmunu
 integer,     intent(out) :: ierr
 ierr = 1
end subroutine nr_init

subroutine nr_evolve_to(t_end,nsteps)
 real(nr_dp), intent(in)            :: t_end
 integer,     intent(out), optional :: nsteps
 if (present(nsteps)) nsteps = 0
end subroutine nr_evolve_to

subroutine nr_get_metric(gcov,gcon,sqrtg,metricderivs)
 real(nr_dp), intent(out) :: gcov(0:,0:,:,:,:),gcon(0:,0:,:,:,:),sqrtg(:,:,:)
 real(nr_dp), intent(out) :: metricderivs(0:,0:,:,:,:,:)
 gcov = 0.
 gcon = 0.
 sqrtg = 0.
 metricderivs = 0.
end subroutine nr_get_metric

real(nr_dp) function nr_time()
 nr_time = 0.
end function nr_time

subroutine nr_constraints(cnorm,hnorm)
 real(nr_dp), intent(out) :: cnorm
 real(nr_dp), intent(out), optional :: hnorm
 cnorm = 0.
 if (present(hnorm)) hnorm = 0.
end subroutine nr_constraints

function nr_formulation() result(name)
 character(len=32) :: name
 name = 'none'
end function nr_formulation

subroutine nr_finalise()
end subroutine nr_finalise

real(nr_dp) function nr_dt_expansion(host_rate)
 real(nr_dp), intent(in), optional :: host_rate
 nr_dt_expansion = huge(1.)
end function nr_dt_expansion

subroutine nr_write_options(iunit)
 integer, intent(in) :: iunit
end subroutine nr_write_options

subroutine nr_read_options(db,nerr)
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
end subroutine nr_read_options

end module nr_api
