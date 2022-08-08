!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module radiation_implicit
!
! Implicit scheme for radiative transfer in flux limited diffusion
! approximation
!
! :References: Wurster, Bate & Monaghan (2006)
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
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
 integer :: icell,i,j,n,ip
 real :: dx,dy,dz
 logical :: iactivei,iamtypei,iamdusti,iamgasi

 !$omp parallel do schedule(runtime)
 over_cells: do icell=1,int(ncells)
    i = ifirstincell(icell)

    !--skip empty cells AND inactive cells
    if (i <= 0) cycle over_cells

    !
    !--get the neighbour list and fill the cell cache
    !
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.true.)

    over_parts: do ip = inoderange(1,icell),inoderange(2,icell)
       i = inodeparts(ip)

       if (maxphase==maxp) then
          call get_partinfo(iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
       else
          iactivei = .true.
          iamtypei = igas
          iamdusti = .false.
          iamgasi  = .true.
       endif

       if (.not.iactivei .or. .not.iamgasi) then ! skip if particle is inactive or not gas
          cycle over_parts
       endif

       loop_over_neigh: do n = 1,nneigh

          j = listneigh(n)
          !--do self contribution separately to avoid problems with 1/sqrt(0.)
          if (j==i) cycle loop_over_neigh

          if (ifilledcellcache .and. n <= isizecellcache) then
             ! positions from cache are already mod boundary
             dx = xpartveci(ixi) - xyzcache(n,1)
             dy = xpartveci(iyi) - xyzcache(n,2)
             dz = xpartveci(izi) - xyzcache(n,3)
          else
             dx = xpartveci(ixi) - xyzh(1,j)
             dy = xpartveci(iyi) - xyzh(2,j)
             dz = xpartveci(izi) - xyzh(3,j)
          endif
#ifdef PERIODIC
          if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
          if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
          if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
          rij2 = dx*dx + dy*dy + dz*dz
          q2i = rij2*hi21
       !--hj is in the cell cache but not in the neighbour cache
       !  as not accessed during the density summation
          if (ifilledcellcache .and. n <= maxcellcache) then
             hj1 = xyzcache(n,4)
          else
             hj1 = 1./xyzh(4,j)
          endif
          hj21 = hj1*hj1
          q2j  = rij2*hj21
!
!--do interaction if r/h < compact support size
!
          is_sph_neighbour: if (q2i < radkern2 .or. q2j < radkern2) then

!$omp critical(listcompact)
             ncompact = ncompact + 1
             ncompacthere = ncompact
             icompacthere = icompact
             icompact = icompact + nneighlocal
!$omp end critical (listcompact)
             if (icompacthere+nneighlocal > icompactmax) then
                call fatal('radiation-implicit','not enough memory allocated for neighbour list')
             endif
             ivar(1,ncompacthere) = nneighlocal
             ivar(2,ncompacthere) = icompacthere
             ivar(3,ncompacthere) = ipart

             do k = 1, nneighlocal
                j = neighlist(k)
                ijvar(icompacthere + k) = j
             enddo

          endif is_sph_neighbour
       enddo loop_over_neigh
    enddo over_parts
 enddo over_cells
 !$omp end parallel do

end subroutine get_compacted_neighbour_list


subroutine compute_drad(varij,varij2,drad)
 real, intent(in) :: varij(:,:),varij2(:,:)
 real, intent(out) :: drad(:,:)

end subroutine compute_drad

end module radiation_implicit
