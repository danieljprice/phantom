module sdar_group
!
! this module contains everything to identify
! and integrate regularized groups...
!
! :References: Makino et Aarseth 2002,Wang et al. 2020, Wang et al. 2021, Rantala et al. 2023
!
! :Owner: Daniel Price
!
 implicit none
 public :: group_identify
 public :: evolve_groups
 ! parameters for group identification
 real, public :: r_neigh  = 0.0
 real, public :: t_crit   = 0.0
 real, public :: C_bin    = 0.0
 real, public :: r_search = 0.0
 private
contains

!
!
! Group identification routines
!
!
subroutine group_identify(nptmass,xyzmh_ptmass,vxyz_ptmass,group_info,nmatrix)
 real,    intent(in)               :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in)               :: group_info(:,:)
 integer(kind=1), intent(inout)    :: nmatrix(:,:)
 integer, intent(in) :: nptmass

 ngroup = 0
 n_ingroup = 0
 n_sing = 0
 call matrix_construction(xyzmh_ptmass,vxyz_ptmass,nmatrix,nptmass)
 call form_group(group_info,nmatrix,nptmass)

end subroutine group_identify


subroutine form_group(group_info,nmatrix,nptmass)
 use part, only : igid,igarg,igsize,igcum
 integer(kind=1), intent(in) :: nmatrix(:,:)
 integer, intent(out):: group_info(:,:)
 integer, intent(in) :: nptmass
 integer :: i
 logical :: visited(nptmass)
 do i=1,nptmass
    if(.not.visited(i)) then
       n_ingroup = n_ingroup + 1
       call dfs(i,i,visited,group_info,nmatrix,nptmass,n_ingroup)
       if (group_info(igsize,i)>1)then
          ngroup = ngroup + 1
          group_info(igcum,ngroup+1) = group_info(igsize,i) + group_info(igcum,ngroup)
       else
          n_ingroup= n_ingroup - 1
          group_info(igsize,i) = 0
          group_info(igarg,nptmass-n_sing) = i
          group_info(igid,nptmass-n_sing) = 0
          n_sing = n_sing + 1
       endif
    endif
 enddo
end subroutine form_group

recursive subroutine dfs(inode,iroot,visited,group_info,nmatrix,npt,n_ingroup)
 use part, only : igid,igarg,igsize,igcum
 integer, intent(in) :: inode,npt,iroot
 integer(kind=1), intent(in) :: nmatrix(:,:)
 integer, intent(inout) :: group_info(:,:)
 integer, intent(inout) :: n_ingroup
 logical, intent(inout) :: visited(:)
 integer :: j
 !print*,nping,inode
 group_info(igarg,n_ingroup) = inode
 group_info(igid,n_ingroup) = iroot
 group_info(igsize,iroot) = group_info(igsize,iroot)+1
 visited(inode) = .true.
 do j=1,npt
    if (nmatrix(inode,j)==1 .and. (visited(j).eqv..false.)) then
       n_ingroup = n_ingroup + 1
       call dfs(j,iroot,visited,group_info,nmatrix,npt,n_ingroup)
    endif
 enddo
end subroutine dfs


subroutine matrix_construction(xyzmh_ptmass,vxyz_ptmass,nmatrix,nptmass)
 use utils_kepler, only: bindE,extract_a,extract_e,extract_ea
 integer(kind=1), intent(out):: nmatrix(:,:)
 real,    intent(in) :: xyzmh_ptmass(:,:)
 real,    intent(in) :: vxyz_ptmass(:,:)
 integer, intent(in) :: nptmass
 real :: xi,yi,zi,vxi,vyi,vzi,mi
 real :: dx,dy,dz,dvx,dvy,dvz,r2,r,v2,mu
 real :: aij,eij,B,rperi
 integer :: i,j

 nmatrix = 0.

 do i=1,nptmass
    xi = xyzmh_ptmass(1,i)
    yi = xyzmh_ptmass(2,i)
    zi = xyzmh_ptmass(3,i)
    mi = xyzmh_ptmass(4,i)
    vxi = vxyz_ptmass(1,i)
    vyi = vxyz_ptmass(2,i)
    vzi = vxyz_ptmass(3,i)
    do j=1,nptmass
       if(i==j) cycle
       if(j>i) then
          dx = xi - xyzmh_ptmass(1,j)
          dy = yi - xyzmh_ptmass(2,j)
          dz = zi - xyzmh_ptmass(3,j)
          r2 = dx**2+dy**2+dz**2
          r = sqrt(r2)
          if (r<r_neigh) then
             nmatrix(i,j) = 1
             cycle
          else if (r>r_search) then
             nmatrix(i,j) = 0
             cycle
          endif
          mu = mi + xyzmh_ptmass(4,j)
          dvx = vxi - vxyz_ptmass(1,j)
          dvy = vyi - vxyz_ptmass(2,j)
          dvz = vzi - vxyz_ptmass(3,j)
          v2 = dvx**2+dvy**2+dvz**2
          call bindE(v2,r,mu,B)
          call extract_a(r,mu,v2,aij)
          if (B<0) then
             if (aij<r_neigh) then
                nmatrix(i,j) = 1
             endif
          else
             call extract_e(dx,dy,dz,dvx,dvy,dvz,mu,r,eij)
             rperi = aij*(1-eij)
             if (rperi<r_neigh .and. C_bin*(r2/v2)<t_crit) then
                nmatrix(i,j) = 1
             endif
          endif
       else
          nmatrix(i,j) = nmatrix(j,i)
       endif
    enddo
 enddo
end subroutine matrix_construction

!
!
! Routines needed to integrate subgroups
!
!

subroutine evolve_groups(n_group,tnext,group_info,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad)
 use part, only: igid,igarg,igsize,igcum
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),gtgrad(:,:)
 integer, intent(in)    :: group_info(:,:)
 integer, intent(inout) :: n_group
 real,    intent(in)    :: tnext
 integer :: i,j,start_id,end_id,gsize
 real :: W,tcoord
 !$omp parallel do default(none)&
 !$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)&
 !$omp shared(tnext)&
 !$omp private(i,j,start_id,end_id,gsize,W,tcoord)&
 do i=1,n_group
    start_id = group_info(igcum,i) + 1
    end_id   = group_info(igcum,i+1)
    gsize    = end_id - start_id
    call integrate_to_time(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad)
 enddo


end subroutine evolve_groups

subroutine integrate_to_time()

end subroutine integrate_to_time



subroutine drift_TTL(tcoord,W,h,xyzmh_ptmass,vxyz_ptmass)
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(inout) :: tcoord
 real, intent(in)    :: h,W
 real :: dtd

 dtd = h/W

 tcoord = tcoord + dtd

 do i=s_id,e_id
    xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + dtd*vxyz_ptmass(1,i)
    xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + dtd*vxyz_ptmass(2,i)
    xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + dtd*vxyz_ptmass(3,i)
 enddo

end subroutine drift_TTL

subroutine kick_TTL(h,W,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad,s_id,e_id)
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),gtgrad(:,:)
 real, intent(in)    :: h
 real, intent(inout) :: W
 integer,intent(in)  :: s_id,e_id
 real :: om,dw,dt
 integer :: i

 call get_force_TTL(xyzmh_ptmass,om,fxyz_ptmass,gtgrad,s_id,e_id)


 dtk = h/om
 do i=s_id,e_id
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + (0.5*dtk)*fxyz_ptmass(1,i)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + (0.5*dtk)*fxyz_ptmass(2,i)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + (0.5*dtk)*fxyz_ptmass(3,i)
 enddo

 dw = 0.
 do i=s_id,e_id
    dw = dw + vxyz_ptmass(1,i)*gtgrad(1,i) + &
              vxyz_ptmass(2,i)*gtgrad(2,i) + &
              vxyz_ptmass(3,i)*gtgrad(3,i)
 enddo

 W = W + dw*dtk

 do i=s_id,e_id
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + (0.5*dtk)*fxyz_ptmass(1,i)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + (0.5*dtk)*fxyz_ptmass(2,i)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + (0.5*dtk)*fxyz_ptmass(3,i)
 enddo


end subroutine kick_TTL

subroutine get_force_TTL(xyzmh_ptmass,om,fxyz_ptmass,gtgrad,s_id,e_id)
 real, intent(in)    :: xyzmh_ptmass(:,:)
 real, intent(inout) :: fxyz_ptmass(:,:),gtgrad(:,:)
 real, intent(out)   :: om
 integer, intent(in) :: s_id,e_id
 real    :: mi,mj,xi,yi,zi,dx,dy,dz,r2
 real    :: gravf,gtki
 integer :: i,j
 om = 0.
 do i=s_id,e_id
    fxyz_ptmass(1,i) = 0.
    fxyz_ptmass(2,i) = 0.
    fxyz_ptmass(3,i) = 0.
    gtgrad(1,i) = 0.
    gtgrad(2,i) = 0.
    gtgrad(3,i) = 0.
 enddo

 do i=s_id,e_id
    gtki = 0.
    xi = xyzmh_ptmass(1,j)
    yi = xyzmh_ptmass(2,j)
    zi = xyzmh_ptmass(3,j)
    mi = xyzmh_ptmass(4,j)
    do j=s_id,e_id
       if (i==j) cycle
       dx = xi - xyzmh_ptmass(1,j)
       dy = yi - xyzmh_ptmass(2,j)
       dz = zi - xyzmh_ptmass(3,j)
       r2 = dx**2+dy**2+dz**3
       r  = sqrt(r)
       mj = xyzmh_ptmass(4,j)
       gravf = xyzmh_ptmass(4,j)*(1./r2*r)
       gtki  = gtki + mj*(1./r)
       fxyz_ptmass(1,i) = fxyz_ptmass(1,i) + dx*gravf
       fxyz_ptmass(2,i) = fxyz_ptmass(2,i) + dy*gravf
       fxyz_ptmass(3,i) = fxyz_ptmass(3,i) + dz*gravf
       gtgrad(1,i) = gtgrad(1,i) + dx*gravf*mi
       gtgrad(2,i) = gtgrad(2,i) + dy*gravf*mi
       gtgrad(3,i) = gtgrad(3,i) + dz*gravf*mi
    enddo
    om = om + gtki*m(i)
 enddo

 om = om*0.5

end subroutine get_force_TTL

subroutine initial_OM(xyzmh_ptmass,om,s_id,e_id)
 real, intent(in)    :: xyzmh_ptmass(:,:)
 real, intent(out)   :: om
 integer, intent(in) :: s_id,e_id
 integer :: i,j
 real    :: gtki,dx,dy,dz,xi,yi,zi,r1

 om = 0.

 do i=s_id,e_id
    gtki = 0.
    xi = xyzmh_ptmass(1,i)
    yi = xyzmh_ptmass(2,i)
    zi = xyzmh_ptmass(3,i)
    mi = xyzmh_ptmass(4,i)
    do j=s_id,e_id
       if (i == j) cycle
       dx = xi - xyzmh_ptmass(1,j)
       dy = yi - xyzmh_ptmass(2,j)
       dz = zi - xyzmh_ptmass(3,j)
       r1 = 1./sqrt(dx**2+dy**2+dz**2)
       gtki = gtki + xyzmh_ptmass(4,j)*r1
    enddo
    om = om + gtki*mi
 enddo

 om = om*0.5

end subroutine initial_OM

end module sdar_group
