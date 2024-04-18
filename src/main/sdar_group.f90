module sdar_group
!
! this module contains everything to identify
! and integrate regularized groups...
!
! :References: Makkino et Aarseth 2002,Wang et al. 2020, Wang et al. 2021, Rantala et al. 2023
!
! :Owner: Yann BERNARD
!
 use utils_sdar
 implicit none
 public :: group_identify
 public :: evolve_groups
 ! parameters for group identification
 real, parameter :: eta_pert = 20
 real, parameter :: time_error = 2.5e-14
 real, parameter :: max_step = 100000000
 real, parameter, public :: r_neigh  = 0.0001
 real,            public :: t_crit   = 1.e-9
 real,            public :: C_bin    = 0.02
 real,            public :: r_search = 100.*r_neigh
 private
contains

!-----------------------------------------------
!
! Group identification routines
!
!-----------------------------------------------
subroutine group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,vxyz_ptmass,group_info,nmatrix)
 use io ,only:id,master,iverbose,iprint
 integer, intent(in)            :: nptmass
 real,    intent(in)            :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout)         :: group_info(3,nptmass)
 integer(kind=1), intent(inout) :: nmatrix(nptmass,nptmass)
 integer, intent(inout)         :: n_group,n_ingroup,n_sing

 n_group = 0
 n_ingroup = 0
 n_sing = 0

 call matrix_construction(xyzmh_ptmass,vxyz_ptmass,nmatrix,nptmass)
 call form_group(group_info,nmatrix,nptmass,n_group,n_ingroup,n_sing)

 if (id==master .and. iverbose>1) then
    write(iprint,"(i6,a,i6,a,i6,a)") n_group," groups identified, ",n_ingroup," in a group, ",n_sing," singles..."
 endif

end subroutine group_identify


subroutine form_group(group_info,nmatrix,nptmass,n_group,n_ingroup,n_sing)
 use part, only : igarg,igcum,igid
 integer,         intent(in)    :: nptmass
 integer(kind=1), intent(inout) :: nmatrix(nptmass,nptmass)
 integer,         intent(inout) :: group_info(3,nptmass)
 integer,         intent(inout) :: n_group,n_ingroup,n_sing
 integer :: i,ncg
 logical :: visited(nptmass)
 visited = .false.
 group_info(igcum,1) = 0
 do i=1,nptmass
    if(.not.visited(i)) then
       n_ingroup = n_ingroup + 1
       call dfs(i,visited,group_info,nmatrix,nptmass,n_ingroup,ncg)
       if (ncg>1)then
          n_group = n_group + 1
          group_info(igcum,n_group+1) = (ncg) + group_info(igcum,n_group)
       else
          n_ingroup = n_ingroup - 1
          group_info(igarg,nptmass-n_sing) = i
          group_info(igid,nptmass-n_sing) = i
          n_sing = n_sing + 1
       endif
    endif
 enddo
end subroutine form_group

subroutine dfs(iroot,visited,group_info,nmatrix,nptmass,n_ingroup,ncg)
 use part, only : igarg,igid
 integer,         intent(in)    :: nptmass,iroot
 integer,         intent(out)   :: ncg
 integer(kind=1), intent(in)    :: nmatrix(nptmass,nptmass)
 integer,         intent(inout) :: group_info(3,nptmass)
 integer,         intent(inout) :: n_ingroup
 logical,         intent(inout) :: visited(nptmass)
 integer :: stack(nptmass)
 integer :: j,stack_top,inode

 stack_top = 0
 ncg = 1
 inode = iroot
 group_info(igarg,n_ingroup) = inode
 group_info(igid,n_ingroup) = iroot
 stack_top = stack_top + 1
 stack(stack_top) = inode
 visited(inode) = .true.
 do while(stack_top>0)
    inode = stack(stack_top)
    stack_top = stack_top - 1
    do j= 1,nptmass
       if (nmatrix(inode,j)==1 .and. .not.(visited(j))) then
          n_ingroup = n_ingroup + 1
          ncg = ncg + 1
          stack_top = stack_top + 1
          stack(stack_top) = j
          visited(j) = .true.
          group_info(igarg,n_ingroup) = j
          group_info(igid,n_ingroup) = iroot
       endif
    enddo
 enddo
end subroutine dfs


subroutine matrix_construction(xyzmh_ptmass,vxyz_ptmass,nmatrix,nptmass)
 use utils_kepler, only: Espec,extract_a,extract_e,extract_ea
 integer, intent(in) :: nptmass
 integer(kind=1), intent(out):: nmatrix(nptmass,nptmass)
 real,    intent(in) :: xyzmh_ptmass(:,:)
 real,    intent(in) :: vxyz_ptmass(:,:)
 real :: xi,yi,zi,vxi,vyi,vzi,mi
 real :: dx,dy,dz,dvx,dvy,dvz,r2,r,v2,mu
 real :: aij,eij,B,rperi
 integer :: i,j
!
!!TODO MPI Proof version of the matrix construction
!

 !$omp parallel do default(none) &
 !$omp shared(nptmass,C_bin,t_crit,nmatrix) &
 !$omp shared(xyzmh_ptmass,vxyz_ptmass,r_search) &
 !$omp private(xi,yi,zi,mi,vxi,vyi,vzi,i,j) &
 !$omp private(dx,dy,dz,r,r2) &
 !$omp private(dvx,dvy,dvz,v2) &
 !$omp private(mu,aij,eij,B,rperi)
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
       call Espec(v2,r,mu,B)
       call extract_a(r,mu,v2,aij)
       if (B<0) then
          if (aij<r_neigh) then
             nmatrix(i,j) = 1
          endif
       else
          call extract_e(dx,dy,dz,dvx,dvy,dvz,mu,r,eij)
          rperi = aij*(1-eij)
          if (rperi<r_neigh .and. C_bin*sqrt(r2/v2)<t_crit) then
             nmatrix(i,j) = 1
          endif
       endif
    enddo
 enddo
 !$omp end parallel do
end subroutine matrix_construction

!---------------------------------------------
!
! Routines needed to integrate subgroups
!
!---------------------------------------------

subroutine evolve_groups(n_group,nptmass,time,tnext,group_info,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad)
 use part, only: igarg,igcum
 use io, only: id,master
 use mpiutils,only:bcast_mpi
 integer, intent(in)    :: n_group,nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),gtgrad(:,:)
 integer, intent(in)    :: group_info(:,:)
 real,    intent(in)    :: tnext,time
 integer :: i,start_id,end_id,gsize
 if (n_group>0) then
    if(id==master) then
       !$omp parallel do default(none)&
       !$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)&
       !$omp shared(tnext,time,group_info,gtgrad)&
       !$omp private(i,start_id,end_id,gsize)
       do i=1,n_group
          start_id = group_info(igcum,i) + 1
          end_id   = group_info(igcum,i+1)
          gsize    = (end_id - start_id) + 1
          call integrate_to_time(start_id,end_id,gsize,time,tnext,xyzmh_ptmass,vxyz_ptmass,group_info,fxyz_ptmass,gtgrad)
       enddo
       !$omp end parallel do
    endif
 endif

 call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
 call bcast_mpi(vxyz_ptmass(:,1:nptmass))


end subroutine evolve_groups

subroutine integrate_to_time(start_id,end_id,gsize,time,tnext,xyzmh_ptmass,vxyz_ptmass,group_info,fxyz_ptmass,gtgrad)
 use part, only: igarg
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:), &
                        fxyz_ptmass(:,:),gtgrad(:,:)
 integer, intent(in) :: group_info(:,:)
 integer, intent(in) :: start_id,end_id,gsize
 real,    intent(in) :: tnext,time
 real, allocatable   :: bdata(:)
 real    :: ds(2)
 real    :: time_table(ck_size)
 integer :: switch
 integer :: step_count_int,step_count_tsyn,n_step_end
 real    :: dt,ds_init,dt_end,step_modif,t_old,W_old
 real    :: W,tcoord
 logical :: t_end_flag,backup_flag,ismultiple
 integer :: i,prim,sec


 tcoord = time

 ismultiple = gsize > 2

 call initial_int(xyzmh_ptmass,fxyz_ptmass,vxyz_ptmass,group_info,W,start_id,end_id,ismultiple,ds_init)

 allocate(bdata(gsize*6))

 step_count_int  = 0
 step_count_tsyn = 0
 n_step_end = 0
 t_end_flag = .false.
 backup_flag = .true.
 ds(:) = ds_init
 switch = 1

 !print*,ds_init, tcoord,tnext,W

 do while (.true.)

    if (backup_flag) then
       call backup_data(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,group_info,fxyz_ptmass,bdata)
    else
       call restore_state(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,group_info,tcoord,t_old,W,W_old,bdata)
    endif
    t_old = tcoord
    W_old = W
    if (gsize>2) then
       do i=1,ck_size
          call drift_TTL (tcoord,W,ds(switch)*cks(i),xyzmh_ptmass,vxyz_ptmass,group_info,start_id,end_id)
          time_table(i) = tcoord
          call kick_TTL  (ds(switch)*dks(i),W,xyzmh_ptmass,vxyz_ptmass,group_info,fxyz_ptmass,gtgrad,start_id,end_id)
       enddo
    else
       prim = group_info(igarg,start_id)
       sec  = group_info(igarg,end_id)
       call oneStep_bin(tcoord,W,ds(switch),xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad,time_table,prim,sec)
    endif
    dt = tcoord - t_old

    step_count_int = step_count_int + 1

    if(step_count_int > max_step) then
       print*,"MAX STEP NUMBER, ABORT !!!"
       call abort
    endif

    if ((.not.t_end_flag).and.(dt<0.)) then
       !print*,"neg dt !!!",tnext,dt
       call regularstepfactor((abs(tnext/dt))**(1./6.),step_modif)
       step_modif = min(max(step_modif,0.0625),0.5)
       ds(switch) = ds(switch)*step_modif
       ds(3-switch) = ds(switch)

       backup_flag = .false.
       continue
    endif

    if (tcoord < tnext - time_error) then
       if (t_end_flag .and. (ds(switch)==ds(3-switch))) then
          step_count_tsyn = step_count_tsyn + 1
          dt_end = tnext - tcoord
          if (dt<0.) then
             call regularstepfactor((abs(tnext/dt))**(1./6.),step_modif)
             step_modif = min(max(step_modif,0.0625),0.5)
             ds(switch)   = ds(switch)*step_modif
             ds(3-switch) = ds(switch)
          else if ((n_step_end > 1) .and. (dt<0.3*dt_end)) then
             ds(3-switch) = ds(switch) * dt_end/dt
          else
             n_step_end = n_step_end + 1
          endif
       endif
       ds(switch) = ds(3-switch)
       switch = 3 - switch
       if (dt>0) then
          backup_flag = .true.
       else
          backup_flag  = .false.
       endif

    else if (tcoord > tnext + time_error) then
       t_end_flag = .true.
       backup_flag = .false.
       n_step_end = 0
       step_count_tsyn = step_count_tsyn + 1

       call new_ds_sync_sup(ds,time_table,tnext,switch)
    else
       exit
    endif
 enddo

 !print*,step_count_int,tcoord,tnext,ds_init

 deallocate(bdata)

end subroutine integrate_to_time


subroutine regularstepfactor(fac_in,fac_out)
 real, intent(in) :: fac_in
 real, intent(out):: fac_out
 fac_out = 1.0
 if (fac_in<1) then
    do while (fac_out>fac_in)
       fac_out = fac_out*0.5
    enddo
 else
    do while(fac_out<=fac_in)
       fac_out = fac_out *2
    enddo
    fac_out = fac_out*0.5
 endif
end subroutine regularstepfactor

subroutine new_ds_sync_sup(ds,time_table,tnext,switch)
 real,    intent(inout) :: ds(:)
 real,    intent(in)    :: time_table(:)
 real,    intent(in)    :: tnext
 integer, intent(in)    :: switch
 integer :: i,k
 real :: tp,dtc,dstmp
 do i=1,ck_size
    k = cck_sorted_id(i)
    if(tnext<time_table(k)) exit
 enddo

 if (i==1) then
    ds(switch) = ds(switch)*cck_sorted(i)*(tnext/time_table(k))
    ds(3-switch) = ds(switch)

 else
    tp  = time_table(cck_sorted_id(i-1))
    dtc = time_table(k)-tp
    dstmp = ds(switch)
    ds(switch) = ds(switch)*cck_sorted(i-1)
    ds(3-switch) = dstmp*(cck_sorted(i)-cck_sorted(i-1))*min(1.0,(tnext-tp+time_error)/dtc)
 endif

end subroutine new_ds_sync_sup



subroutine backup_data(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,group_info,fxyz_ptmass,bdata)
 use part, only: igarg
 real, intent(in)   ::xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 integer,intent(in) :: group_info(:,:)
 real, intent(out)  ::bdata(:)
 integer,intent(in) :: start_id,end_id
 integer :: i,j,k
 j=0
 do k=start_id,end_id
    i = group_info(igarg,k)
    bdata(j*6+1) = xyzmh_ptmass(1,i)
    bdata(j*6+2) = xyzmh_ptmass(2,i)
    bdata(j*6+3) = xyzmh_ptmass(3,i)
    bdata(j*6+4) = vxyz_ptmass(1,i)
    bdata(j*6+5) = vxyz_ptmass(2,i)
    bdata(j*6+6) = vxyz_ptmass(3,i)
    j = j + 1
 enddo

end subroutine backup_data


subroutine restore_state(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,group_info,tcoord,t_old,W,W_old,bdata)
 use part, only: igarg
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 integer,intent(in)  :: group_info(:,:)
 real, intent(out)   :: tcoord,W
 real, intent(in)    :: t_old,W_old
 real, intent(in)    :: bdata(:)
 integer, intent(in) :: start_id,end_id
 integer :: k,i,j
 j = 0
 do k=start_id,end_id
    i = group_info(igarg,k)
    xyzmh_ptmass(1,i) = bdata(j*6+1)
    xyzmh_ptmass(2,i) = bdata(j*6+2)
    xyzmh_ptmass(3,i) = bdata(j*6+3)
    vxyz_ptmass(1,i)  = bdata(j*6+4)
    vxyz_ptmass(2,i)  = bdata(j*6+5)
    vxyz_ptmass(3,i)  = bdata(j*6+6)
    j = j + 1
 enddo
 tcoord   = t_old
 W = W_old

end subroutine restore_state


subroutine drift_TTL(tcoord,W,h,xyzmh_ptmass,vxyz_ptmass,group_info,s_id,e_id)
 use part, only: igarg
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer,intent(in)  :: group_info(:,:)
 real, intent(inout) :: tcoord
 real, intent(in)    :: h,W
 integer,intent(in)  :: s_id,e_id
 integer :: k,i
 real :: dtd

 dtd = h/W

 tcoord = tcoord + dtd

 do k=s_id,e_id
    i = group_info(igarg,k)
    xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + dtd*vxyz_ptmass(1,i)
    xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + dtd*vxyz_ptmass(2,i)
    xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + dtd*vxyz_ptmass(3,i)
 enddo

end subroutine drift_TTL

subroutine kick_TTL(h,W,xyzmh_ptmass,vxyz_ptmass,group_info,fxyz_ptmass,gtgrad,s_id,e_id)
 use part, only: igarg
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),gtgrad(:,:)
 integer,intent(in)  :: group_info(:,:)
 real, intent(in)    :: h
 real, intent(inout) :: W
 integer,intent(in)  :: s_id,e_id
 real :: om,dw,dtk
 integer :: i,k

 call get_force_TTL(xyzmh_ptmass,group_info,fxyz_ptmass,gtgrad,om,s_id,e_id)


 dtk = h/om
 do k=s_id,e_id
    i=group_info(igarg,k)
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + (0.5*dtk)*fxyz_ptmass(1,i)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + (0.5*dtk)*fxyz_ptmass(2,i)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + (0.5*dtk)*fxyz_ptmass(3,i)
 enddo

 dw = 0.
 do k=s_id,e_id
    i=group_info(igarg,k)
    dw = dw + vxyz_ptmass(1,i)*gtgrad(1,i) + &
              vxyz_ptmass(2,i)*gtgrad(2,i) + &
              vxyz_ptmass(3,i)*gtgrad(3,i)
 enddo

 W = W + dw*dtk

 do k=s_id,e_id
    i=group_info(igarg,k)
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + (0.5*dtk)*fxyz_ptmass(1,i)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + (0.5*dtk)*fxyz_ptmass(2,i)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + (0.5*dtk)*fxyz_ptmass(3,i)
 enddo


end subroutine kick_TTL

subroutine oneStep_bin(tcoord,W,ds,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad,time_table,i,j)
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),gtgrad(:,:),time_table(:)
 real, intent(in)    :: ds
 real, intent(inout) :: tcoord,W
 integer, intent(in) :: i,j
 integer :: k
 real :: dtd,dtk,dvel1(3),dvel2(3),dw,om

 do k = 1,ck_size
    dtd = ds*cks(k)/W
    tcoord = tcoord + dtd
    !if (i == 1) print*, fxyz_ptmass(1,i),i,j
    time_table(k) = tcoord

    xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + dtd*vxyz_ptmass(1,i)
    xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + dtd*vxyz_ptmass(2,i)
    xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + dtd*vxyz_ptmass(3,i)
    xyzmh_ptmass(1,j) = xyzmh_ptmass(1,j) + dtd*vxyz_ptmass(1,j)
    xyzmh_ptmass(2,j) = xyzmh_ptmass(2,j) + dtd*vxyz_ptmass(2,j)
    xyzmh_ptmass(3,j) = xyzmh_ptmass(3,j) + dtd*vxyz_ptmass(3,j)

    call get_force_TTL_bin(xyzmh_ptmass,fxyz_ptmass,gtgrad,om,i,j)

    dtk = ds*dks(k)/om

    dvel1(1) = 0.5*dtk*fxyz_ptmass(1,i)
    dvel1(2) = 0.5*dtk*fxyz_ptmass(2,i)
    dvel1(3) = 0.5*dtk*fxyz_ptmass(3,i)
    dvel2(1) = 0.5*dtk*fxyz_ptmass(1,j)
    dvel2(2) = 0.5*dtk*fxyz_ptmass(2,j)
    dvel2(3) = 0.5*dtk*fxyz_ptmass(3,j)

    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + dvel1(1)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + dvel1(2)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + dvel1(3)
    vxyz_ptmass(1,j) = vxyz_ptmass(1,j) + dvel2(1)
    vxyz_ptmass(2,j) = vxyz_ptmass(2,j) + dvel2(2)
    vxyz_ptmass(3,j) = vxyz_ptmass(3,j) + dvel2(3)

    dw = gtgrad(1,i)*vxyz_ptmass(1,i)+&
         gtgrad(2,i)*vxyz_ptmass(2,i)+&
         gtgrad(3,i)*vxyz_ptmass(3,i)+&
         gtgrad(1,j)*vxyz_ptmass(1,j)+&
         gtgrad(2,j)*vxyz_ptmass(2,j)+&
         gtgrad(3,j)*vxyz_ptmass(3,j)

    W = W + dw*dtk

    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + dvel1(1)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + dvel1(2)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + dvel1(3)
    vxyz_ptmass(1,j) = vxyz_ptmass(1,j) + dvel2(1)
    vxyz_ptmass(2,j) = vxyz_ptmass(2,j) + dvel2(2)
    vxyz_ptmass(3,j) = vxyz_ptmass(3,j) + dvel2(3)
 enddo


end subroutine oneStep_bin


subroutine get_force_TTL(xyzmh_ptmass,group_info,fxyz_ptmass,gtgrad,om,s_id,e_id)
 use part, only: igarg
 real, intent(in)    :: xyzmh_ptmass(:,:)
 real, intent(inout) :: fxyz_ptmass(:,:),gtgrad(:,:)
 integer,intent(in)  :: group_info(:,:)
 real, intent(out)   :: om
 integer, intent(in) :: s_id,e_id
 real    :: mi,mj,xi,yi,zi,dx,dy,dz,r2,r
 real    :: gravf,gtki
 integer :: i,j,k,l
 om = 0.


 do k=s_id,e_id
    i = group_info(igarg,k)
    fxyz_ptmass(1,i) = 0.
    fxyz_ptmass(2,i) = 0.
    fxyz_ptmass(3,i) = 0.
    gtgrad(1,i) = 0.
    gtgrad(2,i) = 0.
    gtgrad(3,i) = 0.
    gtki = 0.
    xi = xyzmh_ptmass(1,i)
    yi = xyzmh_ptmass(2,i)
    zi = xyzmh_ptmass(3,i)
    mi = xyzmh_ptmass(4,i)
    do l=s_id,e_id
       if (k==l) cycle
       j = group_info(igarg,l)
       dx = xi - xyzmh_ptmass(1,j)
       dy = yi - xyzmh_ptmass(2,j)
       dz = zi - xyzmh_ptmass(3,j)
       r2 = dx**2+dy**2+dz**2
       r  = sqrt(r2)
       mj = xyzmh_ptmass(4,j)
       gravf = mj*(1./(r2*r))
       gtki  = gtki + mj*(1./r)
       fxyz_ptmass(1,i) = fxyz_ptmass(1,i) + dx*gravf
       fxyz_ptmass(2,i) = fxyz_ptmass(2,i) + dy*gravf
       fxyz_ptmass(3,i) = fxyz_ptmass(3,i) + dz*gravf
       gtgrad(1,i) = gtgrad(1,i) + dx*gravf*mi
       gtgrad(2,i) = gtgrad(2,i) + dy*gravf*mi
       gtgrad(3,i) = gtgrad(3,i) + dz*gravf*mi
    enddo
    om = om + gtki*mi
 enddo

 om = om*0.5

end subroutine get_force_TTL

subroutine get_force_TTL_bin(xyzmh_ptmass,fxyz_ptmass,gtgrad,om,i,j)
 real,    intent(in)    :: xyzmh_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:),gtgrad(:,:)
 integer, intent(in)    :: i,j
 real,    intent(out)   :: om
 real :: dx,dy,dz,r2,r,ddr3,mi,mj
 real :: gravfi,gtki,gravfj,gtkj

 mi = xyzmh_ptmass(4,i)
 mj = xyzmh_ptmass(4,j)
 dx = xyzmh_ptmass(1,i) - xyzmh_ptmass(1,j)
 dy = xyzmh_ptmass(2,i) - xyzmh_ptmass(2,j)
 dz = xyzmh_ptmass(3,i) - xyzmh_ptmass(3,j)
 r2 = dx**2+dy**2+dz**2
 r  = sqrt(r2)
 ddr3 = (1./(r2*r))
 gravfi = mj*ddr3
 gravfj = mi*ddr3
 gtki  = mj*(1./r)
 gtkj  = mi*(1./r)

 fxyz_ptmass(1,i) = -dx*gravfi
 fxyz_ptmass(2,i) = -dy*gravfi
 fxyz_ptmass(3,i) = -dz*gravfi
 fxyz_ptmass(1,j) = dx*gravfj
 fxyz_ptmass(2,j) = dy*gravfj
 fxyz_ptmass(3,j) = dz*gravfj

 gtgrad(1,i) = -dx*gravfi*mi
 gtgrad(2,i) = -dy*gravfi*mi
 gtgrad(3,i) = -dz*gravfi*mi
 gtgrad(1,j) = dx*gravfj*mj
 gtgrad(2,j) = dy*gravfj*mj
 gtgrad(3,j) = dz*gravfj*mj

 om = gtki*mi

end subroutine get_force_TTL_bin

subroutine initial_int(xyzmh_ptmass,fxyz_ptmass,vxyz_ptmass,group_info,om,s_id,e_id,ismultiple,ds_init)
 use utils_kepler, only :extract_a_dot,extract_a,Espec
 use part, only:igarg
 real, intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(inout) :: fxyz_ptmass(:,:)
 integer,intent(in)  :: group_info(:,:)
 real, intent(out)   :: om,ds_init
 logical, intent(in) :: ismultiple
 integer, intent(in) :: s_id,e_id
 real    :: mi,mj,mu,xi,yi,zi,dx,dy,dz,r,r2
 real    :: vxi,vyi,vzi,v2,vi,dvx,dvy,dvz,v,rdotv,axi,ayi,azi,acc,gravfi
 real    :: gravf,gtki
 real    :: Edot,E,semi,semidot
 integer :: k,l,i,j

 Edot = 0.
 E = 0.
 om = 0.

 do k=s_id,e_id
    i = group_info(igarg,k)
    gtki = 0.
    gravfi = 0.
    xi = xyzmh_ptmass(1,i)
    yi = xyzmh_ptmass(2,i)
    zi = xyzmh_ptmass(3,i)
    mi = xyzmh_ptmass(4,i)
    vxi = vxyz_ptmass(1,i)
    vyi = vxyz_ptmass(2,i)
    vzi = vxyz_ptmass(3,i)
    fxyz_ptmass(1,i) = 0.
    fxyz_ptmass(2,i) = 0.
    fxyz_ptmass(3,i) = 0.
    do l=s_id,e_id
       if (k==l) cycle
       j = group_info(igarg,l)
       dx = xi - xyzmh_ptmass(1,j)
       dy = yi - xyzmh_ptmass(2,j)
       dz = zi - xyzmh_ptmass(3,j)
       dvx = vxi - vxyz_ptmass(1,j)
       dvy = vyi - vxyz_ptmass(2,j)
       dvz = vzi - vxyz_ptmass(3,j)
       r2 = dx**2+dy**2+dz**2
       r  = sqrt(r2)
       mj = xyzmh_ptmass(4,j)
       gravf = xyzmh_ptmass(4,j)*(1./r2*r)
       gtki  = gtki + mj*(1./r)
       fxyz_ptmass(1,i) = fxyz_ptmass(1,i) + dx*gravf
       fxyz_ptmass(2,i) = fxyz_ptmass(2,i) + dy*gravf
       fxyz_ptmass(3,i) = fxyz_ptmass(3,i) + dz*gravf
       if (ismultiple) then
          rdotv = dx*dvx + dy*dvy + dz*dvz
          gravfi = gravfi + gravf*rdotv
       else
          v2 = dvx**2 + dvy**2 + dvz**2
          v = sqrt(v2)
       endif

    enddo
    om = om + gtki*mi
    axi = fxyz_ptmass(1,i)
    ayi = fxyz_ptmass(2,i)
    azi = fxyz_ptmass(3,i)
    acc = sqrt(axi**2 + ayi**2 + azi**2)
    if (ismultiple) then
       vi = sqrt(vxi**2 + vyi**2 + vzi**2)
       Edot = Edot + mi*(vi*acc - gravfi)
       E = E + 0.5*mi*vi**2 - gtki
    else
       mu = mi*mj
       call extract_a_dot(r2,r,mu,v2,v,acc,semidot)
       call extract_a(r,mu,v2,semi)
    endif
 enddo

 om = om*0.5

 if (ismultiple) then
    ds_init = eta_pert * abs(E/Edot)
 else
    ds_init = eta_pert * abs(semi/semidot)
 endif

end subroutine initial_int

end module sdar_group
