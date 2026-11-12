!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module omputils
!
! Utility subroutines specific to openMP parallelisation
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!

!$ use dim, only:maxptmass
 implicit none
!$ integer, parameter :: nlocks = maxptmass
!$ integer(kind=8), dimension(0:nlocks) :: ipart_omp_lock

 integer :: omp_num_threads

 ! chunk size for the particle loops in substepping and cons2prim:
 ! zero means static, otherwise dynamic (set when the cores have unequal speeds)
 integer :: loop_chunk = 0
 integer :: saved_kind,saved_chunk

contains
!----------------------------------------------------------------
!+
!  info routine
!+
!----------------------------------------------------------------
subroutine info_omp
!$ integer, external :: omp_get_num_threads

!$omp parallel default(none)
!$omp single
!$ print "(a,i4,a)",' Running in openMP on',omp_get_num_threads(),' threads'
!$omp end single
!$omp end parallel
!$ if (loop_chunk > 0) print "(a,i0,a)",' Unequal cpu speeds: particle loops use schedule(dynamic,',loop_chunk,')'

!$ if (.false.) then
    print "(a)",' openMP parallelisation is OFF'
!$ endif

end subroutine info_omp

!----------------------------------------------------------------
!+
!  initialisation routine necessary if locks are used.
!  Particle loops use a dynamic schedule if the cores have
!  unequal speeds, otherwise static
!+
!----------------------------------------------------------------
subroutine init_omp(unequal_speeds)
 logical, intent(in) :: unequal_speeds
!$ integer :: i
!$ external :: omp_init_lock
!$ integer, external :: omp_get_num_threads

 omp_num_threads = 1

!$ do i = 0, nlocks
!$  call omp_init_lock(ipart_omp_lock(i))
!$ enddo

!$omp parallel default(none) shared(omp_num_threads)
!$ omp_num_threads = omp_get_num_threads()
!$omp end parallel

 if (unequal_speeds) loop_chunk = 2048

end subroutine init_omp

!----------------------------------------------------------------
!+
!  interface to omp_get_num_threads that works even if
!  openMP is off
!+
!----------------------------------------------------------------
integer function omp_thread_num()
!$ integer, external :: omp_get_thread_num

 omp_thread_num = 0
!$ omp_thread_num = omp_get_thread_num()

end function omp_thread_num

!----------------------------------------------------------------
!+
!  set the openMP runtime schedule for the particle loops
!  that use schedule(runtime). Call restore_loop_schedule
!  after the loop, so other schedule(runtime) loops still
!  follow OMP_SCHEDULE
!+
!----------------------------------------------------------------
subroutine use_particle_schedule()
!$ use omp_lib, only:omp_get_schedule,omp_set_schedule,omp_sched_kind,omp_sched_static,omp_sched_dynamic
!$ integer(kind=omp_sched_kind) :: kind

!$ call omp_get_schedule(kind,saved_chunk)
!$ saved_kind = kind
!$ if (loop_chunk > 0) then
!$    call omp_set_schedule(omp_sched_dynamic,loop_chunk)
!$ else
!$    call omp_set_schedule(omp_sched_static,0)
!$ endif

end subroutine use_particle_schedule

!----------------------------------------------------------------
!+
!  restore the runtime schedule saved by use_particle_schedule
!+
!----------------------------------------------------------------
subroutine restore_loop_schedule()
!$ use omp_lib, only:omp_set_schedule,omp_sched_kind

!$ call omp_set_schedule(int(saved_kind,kind=omp_sched_kind),saved_chunk)

end subroutine restore_loop_schedule

end module omputils
