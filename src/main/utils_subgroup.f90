!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_subgroup
!
! utils_subgroup
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 real,    parameter :: time_error = 2.5e-12
 real,    parameter :: max_step   = 1000000
 integer, parameter :: ck_size = 8
 integer,dimension(ck_size),parameter :: cck_sorted_id=(/6,1,3,4,5,7,2,8/)
 real,dimension(ck_size),parameter :: cks        = (/0.3922568052387800,0.5100434119184585,-0.4710533854097566,&
                                                     0.0687531682525181,0.0687531682525181,-0.4710533854097566,&
                                                     0.5100434119184585,0.3922568052387800/)
 real,dimension(ck_size),parameter :: cck_sorted = (/0.0976997828427615,0.3922568052387800,0.4312468317474820,&
                                                     0.5000000000000000,0.5687531682525181,0.6077431947612200,&
                                                     0.9023002171572385,1.0000000000000000/)
 real,dimension(ck_size),parameter :: dks        = (/0.7845136104775600,0.2355732133593570,-1.1776799841788701,&
                                                     1.3151863206839063,-1.1776799841788701,0.2355732133593570,&
                                                     0.7845136104775600,0.0000000000000000/)

contains

!----------------------------------------------
!
! Routines to synchronise time integration
!
!----------------------------------------------
subroutine converge_to_tend(ds,dt,t_end,tcoord,time_table,switch,do_store,do_sync,is_end,&
                            nstep_int,nstep_tsync,nstep_end)
 use io, only:fatal
 real,    intent(inout) :: ds(:)
 real,    intent(in)    :: dt,t_end,tcoord,time_table(:)
 integer, intent(inout) :: switch,nstep_int,nstep_tsync,nstep_end
 logical, intent(inout) :: do_store,do_sync,is_end
 real :: step_modif,dt_end

 nstep_int = nstep_int + 1

 if (nstep_int > max_step) call fatal("subgroup",'max number of step reached, abort...')

 if ((.not.do_sync).and.(dt<0.)) then !-- end is not reached but time goes backward...
    call regularstepfactor((abs(t_end/dt))**(1./6.),step_modif)
    step_modif = min(max(step_modif,0.0625),0.5)
    ds(switch) = ds(switch)*step_modif
    ds(3-switch) = ds(switch)
    do_store = .false.
    return
 endif

 if (tcoord < t_end - time_error) then
    if (do_sync .and. (ds(switch)==ds(3-switch))) then
       nstep_tsync = nstep_tsync + 1
       dt_end = t_end - tcoord
       if (dt<0.) then
          call regularstepfactor((abs(t_end/dt))**(1./6.),step_modif)
          step_modif = min(max(step_modif,0.0625),0.5)
          ds(switch)   = ds(switch)*step_modif
          ds(3-switch) = ds(switch)
       elseif ((nstep_end > 1) .and. (dt<0.3*dt_end)) then
          ds(3-switch) = ds(switch) * dt_end/dt
       else
          nstep_end = nstep_end + 1
       endif
    endif
    ds(switch) = ds(3-switch)
    switch = 3 - switch
    if (dt>0) then
       do_store = .true.
    else
       do_store  = .false.
    endif

 elseif (tcoord > t_end + time_error) then
    do_sync     = .true.
    do_store    = .false.
    nstep_end   = 0
    nstep_tsync = nstep_tsync + 1

    call new_ds_sync_sup(ds,time_table,t_end,switch)
 else
    is_end = .true.
 endif

end subroutine converge_to_tend

subroutine new_ds_sync_sup(ds,time_table,t_end,switch)
 real,    intent(inout) :: ds(:)
 real,    intent(in)    :: time_table(:)
 real,    intent(in)    :: t_end
 integer, intent(in)    :: switch
 integer :: i,k
 real :: tp,dtc,dstmp
 do i=1,ck_size
    k = cck_sorted_id(i)
    if (t_end<time_table(k)) exit
 enddo

 if (i==1) then
    ds(switch) = ds(switch)*cck_sorted(i)*(t_end/time_table(k))
    ds(3-switch) = ds(switch)

 else
    tp  = time_table(cck_sorted_id(i-1))
    dtc = time_table(k)-tp
    dstmp = ds(switch)
    ds(switch) = ds(switch)*cck_sorted(i-1)
    ds(3-switch) = dstmp*(cck_sorted(i)-cck_sorted(i-1))*min(1.0,(t_end-tp+time_error)/dtc)
 endif

end subroutine new_ds_sync_sup

subroutine subgroup_step_init(nstep_int,nstep_tsync,nstep_end,do_sync,do_store,is_end,ds,ds_init,switch)
 integer, intent(out) :: nstep_int,nstep_tsync,nstep_end,switch
 logical, intent(out) :: do_sync,do_store,is_end
 real,    intent(out) :: ds(:)
 real,    intent(in)  :: ds_init

 nstep_int   = 0
 nstep_tsync = 0
 nstep_end   = 0
 do_sync     = .false.
 do_store    = .true.
 is_end      = .false.
 ds          = ds_init
 switch      = 1

end subroutine subgroup_step_init

subroutine regularstepfactor(fac_in,fac_out)
 real, intent(in)  :: fac_in
 real, intent(out) :: fac_out

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

subroutine store_state(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,gstate)
 use part, only:igarg,ikappa
 real,    intent(in)  ::xyzmh_ptmass(:,:),vxyz_ptmass(:,:),bin_info(:,:)
 integer, intent(in)  :: group_info(:,:)
 real,    intent(out) ::gstate(:)
 integer, intent(in)  :: start_id,end_id
 integer :: i,j,k
 j=0
 do k=start_id,end_id
    i = group_info(igarg,k)
    gstate(j*7+1) = xyzmh_ptmass(1,i)
    gstate(j*7+2) = xyzmh_ptmass(2,i)
    gstate(j*7+3) = xyzmh_ptmass(3,i)
    gstate(j*7+4) = vxyz_ptmass(1,i)
    gstate(j*7+5) = vxyz_ptmass(2,i)
    gstate(j*7+6) = vxyz_ptmass(3,i)
    gstate(j*7+7) = bin_info(ikappa,i)
    j = j + 1
 enddo

end subroutine store_state

subroutine restore_state(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,tcoord,t_old,W,W_old,gstate)
 use part, only:igarg,ikappa
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),bin_info(:,:)
 integer, intent(in)    :: group_info(:,:)
 real,    intent(out)   :: tcoord,W
 real,    intent(in)    :: t_old,W_old
 real,    intent(in)    :: gstate(:)
 integer, intent(in)    :: start_id,end_id
 integer :: k,i,j
 j = 0
 do k=start_id,end_id
    i = group_info(igarg,k)
    xyzmh_ptmass(1,i)  = gstate(j*7+1)
    xyzmh_ptmass(2,i)  = gstate(j*7+2)
    xyzmh_ptmass(3,i)  = gstate(j*7+3)
    vxyz_ptmass(1,i)   = gstate(j*7+4)
    vxyz_ptmass(2,i)   = gstate(j*7+5)
    vxyz_ptmass(3,i)   = gstate(j*7+6)
    bin_info(ikappa,i) = gstate(j*7+7)
    j = j + 1
 enddo
 tcoord = t_old
 W      = W_old

end subroutine restore_state

end module utils_subgroup
