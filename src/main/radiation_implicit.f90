module radiation_implicit
 implicit none
 integer, parameter :: ierr_failed_to_converge = 1

contains

 subroutine do_radiation_implicit(dt,dtmax,npart,rad,xyzh,vxyzu,drad,ierr)
  integer, intent(in) :: npart
  real, intent(in) :: dt,dtmax,rad(:,:),xyzh(:,:),vxyzu(:,:),drad(:,:)
  integer, intent(out) :: ierr
  integer :: nsubsteps,i,nit
  logical :: failed
  real :: dtsub,errorE,errorU
  
  ierr = 0
  nsubsteps = 1

  dtsub = dt
  do i = 1,nsubsteps
     call do_radiation_onestep(dtsub,dtmax,npart,rad,xyzh,vxyzu,drad,failed,nit,errorE,errorU)
     if (failed) ierr = ierr_failed_to_converge
  enddo

 end subroutine do_radiation_implicit


 subroutine do_radiation_onestep(dt,dtmax,npart,rad,xyzh,vxyzu,drad,failed,nit,errorE,errorU)
  use units, only:get_c_code,get_radconst_code
  integer, intent(in) :: npart
  real, intent(in) :: dt,dtmax,rad(:,:),xyzh(:,:),vxyzu(:,:),drad(:,:)
  logical, intent(out) :: failed
  integer, intent(out) :: nit
  real, intent(out) :: errorE,errorU
  real :: ccode,acode
  integer, allocatable :: ivar(:,:),ijvar(:)
  integer :: ncompact
  real, allocatable :: varij(:,:),varij2(:,:)

  failed = .false.
  nit = 0
  errorE = 0.
  errorU = 0.

  ccode = get_c_code()
  acode = get_radconst_code()
  !dtimax = dt/imaxstep

  call get_compacted_neighbour_list(ivar,ijvar,ncompact)
  call fill_pairwise_arrays(varij,varij2)
  call compute_drad(varij,varij2,drad)

 end subroutine do_radiation_onestep


 subroutine fill_pairwise_arrays(varij,varij2)
  real, allocatable, intent(out) :: varij(:,:),varij2(:,:)


 end subroutine fill_pairwise_arrays


 subroutine get_compacted_neighbour_list(ivar,ijvar,ncompact)
  integer, allocatable, intent(out) :: ivar(:,:),ijvar(:)
  integer, intent(out) :: ncompact 


 end subroutine get_compacted_neighbour_list


 subroutine compute_drad(varij,varij2,drad)
  real, intent(in) :: varij(:,:),varij2(:,:)
  real, intent(out) :: drad(:,:)

 end subroutine compute_drad

end module radiation_implicit