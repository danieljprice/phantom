!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module subgroup
!
! This module contains everything to identify and integrate regularized groups.
! TTL (Mikkola and Aarseth 2002) is used to regularize this subgroups. Indentification is done
! using a fixed searching radius and few arguments on bounding systems (see Rantala et al. 2023)
! Slow down method is now directly implemented in the integration. Kappa is computed using Mikkola
! Aarseth (1996) criterion...
!
! :References: Mikkola et Aarseth 2002,Wang et al. 2020, Wang et al. 2021, Rantala et al. 2023
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: io, mpiutils, part, timing, utils_kepler, utils_subgroup
!
 use utils_subgroup
 implicit none
 public :: group_identify
 public :: evolve_groups
 public :: get_pot_subsys
 public :: init_subgroup
 public :: update_kappa
 !
 !-- parameters for group identification
 !
 real, parameter :: time_error = 2.5e-12
 real, parameter :: max_step   = 1000000
 real, parameter :: C_bin      = 0.02
 real, public    :: r_neigh    = 0.001 ! default value assume udist = 1 pc
 real            :: elli_res   = 1/128.
 real            :: hyper_res  = 1/256.
 real            :: r_search

 !
 !-- parameter for Slow Down method
 !
 real, parameter :: kref = 1.e-6

 private
contains
!-----------------------------------------------
!
! Initialisation routine
!
!-----------------------------------------------
subroutine init_subgroup
 r_search = 100.*r_neigh
end subroutine init_subgroup

!-----------------------------------------------
!
! Group identification routines (Subgroups + binary orbital parameters)
!
!-----------------------------------------------
subroutine group_identify(nptmass,n_group,n_ingroup,n_sing,xyzmh_ptmass,vxyz_ptmass, &
                          group_info,bin_info,nmatrix,dtext,new_ptmass)
 use io,     only:id,master,iverbose,iprint
 use timing, only:get_timings,increment_timer,itimer_sg_id
 integer,           intent(in)    :: nptmass
 real,              intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,              intent(inout) :: bin_info(:,:)
 integer,           intent(inout) :: group_info(4,nptmass)
 integer,           intent(inout) :: n_group,n_ingroup,n_sing
 integer(kind=1),   intent(inout) :: nmatrix(nptmass,nptmass)
 logical, optional, intent(in)    :: new_ptmass
 real, optional,    intent(in)    :: dtext
 real(kind=4) :: t1,t2,tcpu1,tcpu2
 logical      :: large_search,reset_nm


 large_search = present(dtext)
 reset_nm     = present(new_ptmass)
 n_group = 0
 n_ingroup = 0
 n_sing = 0
 if (nptmass > 0) then

    call get_timings(t1,tcpu1)
    group_info(:,:) = 0

    if (reset_nm) nmatrix = 0

    if (large_search) then
       call matrix_construction(xyzmh_ptmass,vxyz_ptmass,nmatrix,nptmass,dtext)
    else
       call matrix_construction(xyzmh_ptmass,vxyz_ptmass,nmatrix,nptmass)
    endif
    call form_group(group_info,nmatrix,nptmass,n_group,n_ingroup,n_sing)

    if (n_group > 0) call find_binaries(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,n_group)

    call get_timings(t2,tcpu2)
    call increment_timer(itimer_sg_id,t2-t1,tcpu2-tcpu1)

 endif

 if (id==master .and. iverbose>1) then
    write(iprint,"(i6,a,i6,a,i6,a)") n_group," groups identified, ",n_ingroup," in a group, ",n_sing," singles..."
 endif

end subroutine group_identify

!------------------------------------------------------------------
!
! routine to find binary properties in multiple and binary systems
!
!------------------------------------------------------------------
subroutine find_binaries(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,n_group)
 use part,         only: igarg,igcum,icomp,isemi,iecc,iapo,iorb

 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: group_info(:,:)
 real,    intent(inout) :: bin_info(:,:)
 integer, intent(in)    :: n_group
 integer :: i,k,l,start_id,end_id,gsize
 real    :: akl,ekl,apokl,Tkl

 bin_info(isemi,:) = 0.
 bin_info(iecc,:) = 0.
 bin_info(iapo,:) = 0.
 bin_info(iorb,:) = 0.

!$omp parallel do default(none)&
!$omp shared(n_group,xyzmh_ptmass,vxyz_ptmass)&
!$omp shared(group_info,bin_info)&
!$omp private(start_id,end_id,gsize)&
!$omp private(akl,ekl,apokl,Tkl,k,l,i)
 do i=1,n_group
    start_id = group_info(igcum,i) + 1
    end_id   = group_info(igcum,i+1)
    gsize    = (end_id - start_id) + 1
    if (gsize > 2) then
       call binaries_in_multiples(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,&
                                  gsize,start_id,end_id)
    else
       k = group_info(igarg,start_id)
       l = group_info(igarg,end_id)
       group_info(icomp,end_id)   = k
       group_info(icomp,start_id) = l
       !
       !-- Compute and store main orbital parameters needed for SDAR method
       !
       call get_orbparams(xyzmh_ptmass,vxyz_ptmass,akl,ekl,apokl,Tkl,k,l)
       bin_info(isemi,k) = akl
       bin_info(isemi,l) = akl
       bin_info(iecc,k)  = ekl
       bin_info(iecc,l)  = ekl
       bin_info(iapo,k)  = apokl
       bin_info(iapo,l)  = apokl
       bin_info(iorb,k)  = Tkl
       bin_info(iorb,l)  = Tkl
    endif
 enddo
!$omp end parallel do

end subroutine find_binaries

!--------------------------------------------------------------------------
!
! specialized routine to find orbital parameters of binaries in multiples
!
!--------------------------------------------------------------------------
subroutine binaries_in_multiples(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,gsize,start_id,end_id)
 use part,         only: igarg,icomp,isemi,iecc,iapo,iorb
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: group_info(:,:)
 real,    intent(inout) :: bin_info(:,:)
 integer, intent(in)    :: start_id, end_id,gsize
 integer, allocatable :: r2min_id(:)
 real    :: akl,ekl,apokl,Tkl
 integer :: np,ns,j,k,l

 group_info(icomp,start_id:end_id) = -1
 allocate(r2min_id(gsize))
 call get_r2min(xyzmh_ptmass,group_info,r2min_id,start_id,end_id)
 do j=start_id,end_id
    np = (j-start_id) + 1
    k = group_info(igarg,j)
    if (group_info(icomp,j) < 0) then
       ns = r2min_id(np)
       if (r2min_id(ns) == np) then ! We found a binary into a subgroup : tag as binary component and compute parameters
          l = group_info(igarg,ns+(start_id-1))
          group_info(icomp,j)               = l
          group_info(icomp,ns+(start_id-1)) = k
          !
          !-- Compute and store main orbital parameters needed for SDAR method
          !
          call get_orbparams(xyzmh_ptmass,vxyz_ptmass,akl,ekl,apokl,Tkl,k,l)
          bin_info(isemi,k) = akl
          bin_info(isemi,l) = akl
          bin_info(iecc,k)  = ekl
          bin_info(iecc,l)  = ekl
          bin_info(iapo,k)  = apokl
          bin_info(iapo,l)  = apokl
          bin_info(iorb,k)  = Tkl
          bin_info(iorb,l)  = Tkl
       else  ! No matches... Only a single
          group_info(icomp,j) = k
          bin_info(isemi,k) = 0.
          bin_info(iecc,k)  = 0.
          bin_info(iapo,k)  = 0.
          bin_info(iorb,k)  = 0.
       endif
    endif
 enddo
 deallocate(r2min_id)


end subroutine binaries_in_multiples

!--------------------------------------------
!
! routine to find common nearest neighbours
!
!--------------------------------------------
subroutine get_r2min(xyzmh_ptmass,group_info,r2min_id,start_id,end_id)
 use part, only : igarg,igcum
 real   , intent(in)    :: xyzmh_ptmass(:,:)
 integer, intent(in)    :: group_info(:,:)
 integer, intent(out)   :: r2min_id(:)
 integer, intent(in)    :: start_id,end_id
 integer :: i,j,k,l,n
 real :: dr(3),r2,r2min

 r2min_id = 0

 do i=start_id,end_id
    n = (i-start_id)+1
    j = group_info(igarg,i)
    r2 = 0.
    r2min = huge(r2)
    do k=start_id,end_id
       l = group_info(igarg,k)
       if (j == l) cycle
       dr(1) = xyzmh_ptmass(1,j) - xyzmh_ptmass(1,l)
       dr(2) = xyzmh_ptmass(2,j) - xyzmh_ptmass(2,l)
       dr(3) = xyzmh_ptmass(3,j) - xyzmh_ptmass(3,l)
       r2 = dr(1)**2+dr(2)**2+dr(3)**2
       if (r2 < r2min) then
          r2min       = r2
          r2min_id(n) = (k-start_id)+1
       endif
    enddo
 enddo

end subroutine get_r2min

!----------------------------------------------------------
!
! routine to extract main orbital parameters needed for SD
!
!----------------------------------------------------------
subroutine get_orbparams(xyzmh_ptmass,vxyz_ptmass,aij,eij,apoij,Tij,i,j)
 use utils_kepler, only: extract_e,extract_a,extract_T
 real, intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(out)   :: aij,eij,apoij,Tij
 integer, intent(in) :: i,j
 real :: dv(3),dr(3),mu,r,v2

 dv(1) = vxyz_ptmass(1,j)  - vxyz_ptmass(1,i)
 dv(2) = vxyz_ptmass(2,j)  - vxyz_ptmass(2,i)
 dv(3) = vxyz_ptmass(3,j)  - vxyz_ptmass(3,i)
 dr(1) = xyzmh_ptmass(1,j) - xyzmh_ptmass(1,i)
 dr(2) = xyzmh_ptmass(2,j) - xyzmh_ptmass(2,i)
 dr(3) = xyzmh_ptmass(3,j) - xyzmh_ptmass(3,i)
 mu    = xyzmh_ptmass(4,i) + xyzmh_ptmass(4,j)

 r = sqrt(dr(1)**2+dr(2)**2+dr(3)**2)
 v2 = dv(1)**2+dv(2)**2+dv(3)**2

 call extract_a(r,mu,v2,aij)

 call extract_e(dr(1),dr(2),dr(3),dv(1),dv(2),dv(3),mu,r,eij)

 call extract_T(mu,aij,Tij)

 apoij = aij*(1+eij)

end subroutine get_orbparams

!--------------------------------------------------------------------------
!
! interface routine to read the adjacency matrix to identify groups member
!
!--------------------------------------------------------------------------
subroutine form_group(group_info,nmatrix,nptmass,n_group,n_ingroup,n_sing)
 use part, only : igarg,igcum,igid,icomp
 integer,         intent(in)    :: nptmass
 integer(kind=1), intent(inout) :: nmatrix(nptmass,nptmass)
 integer,         intent(inout) :: group_info(4,nptmass)
 integer,         intent(inout) :: n_group,n_ingroup,n_sing
 logical, allocatable :: visited(:)
 integer :: i,ncg
 allocate(visited(nptmass))
 visited = .false.
 group_info(igcum,1) = 0
 do i=1,nptmass
    if (.not.visited(i)) then
       n_ingroup = n_ingroup + 1
       call dfs(i,visited,group_info,nmatrix,nptmass,n_ingroup,ncg)
       if (ncg>1) then
          n_group = n_group + 1
          group_info(igcum,n_group+1) = (ncg) + group_info(igcum,n_group)
       else
          n_ingroup = n_ingroup - 1
          group_info(igarg,nptmass-n_sing) = i
          group_info(igid,nptmass-n_sing) = i
          group_info(icomp,nptmass-n_sing) = i
          n_sing = n_sing + 1
       endif
    endif
 enddo
 deallocate(visited)

end subroutine form_group

!--------------------------------------------------------------------------
!
! simple deep first search algorithm to form subgroups of point masses
!
!--------------------------------------------------------------------------
subroutine dfs(iroot,visited,group_info,nmatrix,nptmass,n_ingroup,ncg)
 use part, only : igarg,igid,icomp
 integer,         intent(in)    :: nptmass,iroot
 integer,         intent(out)   :: ncg
 integer(kind=1), intent(in)    :: nmatrix(nptmass,nptmass)
 integer,         intent(inout) :: group_info(4,nptmass)
 integer,         intent(inout) :: n_ingroup
 logical,         intent(inout) :: visited(nptmass)
 integer, allocatable :: stack(:)
 integer              :: j,stack_top,inode

 allocate(stack(nptmass))

 stack_top = 0
 ncg = 1
 inode = iroot
 group_info(igarg,n_ingroup) = inode
 group_info(igid,n_ingroup)  = iroot
 group_info(icomp,n_ingroup) = -1   ! icomp to -1 -> need to be identified later
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
          group_info(icomp,n_ingroup) = -1 ! icomp to -1 -> need to be identified later
       endif
    enddo
 enddo
 deallocate(stack)
end subroutine dfs

!------------------------------------------------------------------------------------------
!
! Adjacency matrix construction routine using fixed searching radius (Rantala et al. 2023)
!
!------------------------------------------------------------------------------------------
subroutine matrix_construction(xyzmh_ptmass,vxyz_ptmass,nmatrix,nptmass,dtext)
 use utils_kepler, only: extract_a,extract_e,extract_ea
 integer,         intent(in) :: nptmass
 real,            intent(in) :: xyzmh_ptmass(:,:)
 real,            intent(in) :: vxyz_ptmass(:,:)
 integer(kind=1), intent(out):: nmatrix(nptmass,nptmass)
 real, optional,  intent(in) :: dtext
 real :: xi,yi,zi,vxi,vyi,vzi,mi,mj
 real :: dx,dy,dz,dvx,dvy,dvz,r2,r,v2,mu
 real :: aij,eij,rperi,dtexti
 integer :: i,j
 if (present(dtext)) then
    dtexti = dtext
 else
    dtexti = 0.
 endif
!
!!TODO MPI Proof version of the matrix construction
!

 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nptmass,dtexti,nmatrix,r_neigh) &
 !$omp shared(xyzmh_ptmass,vxyz_ptmass,r_search) &
 !$omp private(xi,yi,zi,mi,vxi,vyi,vzi,i,j) &
 !$omp private(dx,dy,dz,r,r2,mj) &
 !$omp private(dvx,dvy,dvz,v2) &
 !$omp private(mu,aij,eij,rperi)
 do i=1,nptmass
    xi = xyzmh_ptmass(1,i)
    yi = xyzmh_ptmass(2,i)
    zi = xyzmh_ptmass(3,i)
    mi = xyzmh_ptmass(4,i)
    vxi = vxyz_ptmass(1,i)
    vyi = vxyz_ptmass(2,i)
    vzi = vxyz_ptmass(3,i)
    if (mi <= 0. ) then
       nmatrix(i,:) = 0 ! killed point masses can't be in a group
       cycle
    endif
    do j=1,nptmass
       if (i==j) cycle
       mj = xyzmh_ptmass(4,j)
       if (mj <= 0. ) then
          nmatrix(i,j) = 0 ! killed point masses can't be in a group
          cycle
       endif
       dx = xi - xyzmh_ptmass(1,j)
       dy = yi - xyzmh_ptmass(2,j)
       dz = zi - xyzmh_ptmass(3,j)
       r2 = dx**2+dy**2+dz**2
       r = sqrt(r2)
       !
       !-- searching routine
       !
       if (r<r_neigh) then ! if inside neighbour radius set to 1
          nmatrix(i,j) = 1
       elseif (r<r_search) then ! if inside searching radius need to check
          mu = mi + xyzmh_ptmass(4,j)
          dvx = vxi - vxyz_ptmass(1,j)
          dvy = vyi - vxyz_ptmass(2,j)
          dvz = vzi - vxyz_ptmass(3,j)
          v2 = dvx**2+dvy**2+dvz**2
          call extract_a(r,mu,v2,aij)
          if (aij>0) then ! check if the system is bounded
             if (aij<r_neigh) then ! if yes then check on the semi major axis
                nmatrix(i,j) = 1
             else
                nmatrix(i,j) = 0
             endif
          else ! if hyperbolic encounter, then check the periapsis
             call extract_e(dx,dy,dz,dvx,dvy,dvz,mu,r,eij)
             rperi = aij*(1-eij)
             if ((rperi < r_neigh) .and. (C_bin*sqrt(r2/v2) < dtexti)) then
                nmatrix(i,j) = 1
             else
                nmatrix(i,j) = 0
             endif
          endif
       else
          nmatrix(i,j) = 0
       endif
    enddo
 enddo
 !$omp end parallel do
end subroutine matrix_construction

!---------------------------------------------
!
! Interface routine to integrate subgroups
!
!---------------------------------------------
subroutine evolve_groups(n_group,nptmass,time,tnext,group_info,bin_info, &
                         xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad)
 use part,     only:igarg,igcum
 use io,       only:id,master
 use mpiutils, only:bcast_mpi
 use timing,   only:get_timings,increment_timer,itimer_sg_evol
 integer, intent(in)    :: n_group,nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),gtgrad(:,:)
 real,    intent(inout) :: bin_info(:,:)
 integer, intent(inout) :: group_info(:,:)
 real,    intent(in)    :: tnext,time
 integer      :: i,start_id,end_id,gsize
 real(kind=4) :: t1,t2,tcpu1,tcpu2




 if (n_group>0) then

    call get_timings(t1,tcpu1)

    call find_binaries(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,n_group)

    if (id==master) then

       !$omp parallel do default(none)&
       !$omp shared(xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass)&
       !$omp shared(tnext,time,group_info,bin_info,gtgrad,n_group)&
       !$omp private(i,start_id,end_id,gsize)
       do i=1,n_group
          start_id = group_info(igcum,i) + 1
          end_id   = group_info(igcum,i+1)
          gsize    = (end_id - start_id) + 1
          call integrate_to_time(start_id,end_id,gsize,time,tnext,xyzmh_ptmass,vxyz_ptmass,&
                                 bin_info,group_info,fxyz_ptmass,gtgrad)
       enddo
       !$omp end parallel do
    endif

    call get_timings(t2,tcpu2)
    call increment_timer(itimer_sg_evol,t2-t1,tcpu2-tcpu1)
 endif

 call bcast_mpi(xyzmh_ptmass(:,1:nptmass))
 call bcast_mpi(vxyz_ptmass(:,1:nptmass))



end subroutine evolve_groups

!------------------------------------------------------------------------------------
!
! Main integration routine to evolve subgroups, containing Kick and Drift routines
! and time synchronisation algorithm. cf : Wang et al. (2020)
!
!------------------------------------------------------------------------------------
subroutine integrate_to_time(start_id,end_id,gsize,time,tnext,xyzmh_ptmass,vxyz_ptmass,&
                            bin_info,group_info,fxyz_ptmass,gtgrad)
 use part, only: igarg,ikap,isemi
 use io,   only: fatal
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:), &
                           fxyz_ptmass(:,:),gtgrad(:,:),bin_info(:,:)
 integer, intent(inout) :: group_info(:,:)
 integer, intent(in)    :: start_id,end_id,gsize
 real,    intent(in)    :: tnext,time
 real, allocatable      :: bdata(:)
 real    :: ds(2)
 real    :: time_table(ck_size)
 integer :: switch
 integer :: step_count_int,step_count_tsyn,n_step_end
 real    :: dt,ds_init,dt_end,step_modif,t_old,W_old
 real    :: W,tcoord,kappa1
 logical :: t_end_flag,backup_flag,ismultiple
 integer :: i,prim,sec


 tcoord = time

 ismultiple = gsize > 2

 if (ismultiple) then
    call get_kappa(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,gsize,start_id,end_id)
    call get_force_TTL(xyzmh_ptmass,group_info,bin_info,fxyz_ptmass,gtgrad,W,start_id,end_id,ds_init=ds_init)
 else
    prim = group_info(igarg,start_id)
    sec  = group_info(igarg,end_id)
    !
    !-- We need to compute the force a the beginning of the step ( and kappa if slow down)
    !
    call get_kappa_bin(xyzmh_ptmass,bin_info,prim,sec)
    if (bin_info(ikap,prim) >= 1.) then
       kappa1 = 1./bin_info(ikap,prim)
    else
       kappa1 = 1.
       call fatal('subgroup','kappa value bellow 1... something wrong here!',var='kappa',val=bin_info(ikap,prim))
    endif
    call get_force_TTL_bin(xyzmh_ptmass,fxyz_ptmass,gtgrad,W,kappa1,prim,sec,&
                           ds_init=ds_init,semiij=bin_info(isemi,prim))
 endif


 allocate(bdata(gsize*7))

 step_count_int  = 0
 step_count_tsyn = 0
 n_step_end = 0
 t_end_flag = .false.
 backup_flag = .true.
 ds(:) = ds_init
 switch = 1



 do while (.true.)
    if (backup_flag) then
       call backup_data(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,bdata)
    else
       call restore_state(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,tcoord,t_old,W,W_old,bdata)
    endif
    t_old = tcoord
    W_old = W
    if (ismultiple) then
       do i=1,ck_size
          call drift_TTL (tcoord,W,ds(switch)*cks(i),xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,gsize,start_id,end_id)
          time_table(i) = tcoord
          call kick_TTL  (ds(switch)*dks(i),W,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,fxyz_ptmass,gtgrad,start_id,end_id)
       enddo
    else
       prim = group_info(igarg,start_id)
       sec  = group_info(igarg,end_id)
       call oneStep_bin(tcoord,W,ds(switch),kappa1,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad,time_table,prim,sec)
    endif
    dt = tcoord - t_old

    step_count_int = step_count_int + 1

    if (step_count_int > max_step) then
       print*,"MAX STEP NUMBER, ABORT !!!"
       print*,step_count_int,step_count_tsyn,tcoord,tnext,ds_init,ds(switch)
       call abort()
    endif

    if ((.not.t_end_flag).and.(dt<0.)) then
       !print*,"neg dt !!!",tnext,dt,step_count_int
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
          elseif ((n_step_end > 1) .and. (dt<0.3*dt_end)) then
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

    elseif (tcoord > tnext + time_error) then
       t_end_flag = .true.
       backup_flag = .false.
       n_step_end = 0
       step_count_tsyn = step_count_tsyn + 1

       call new_ds_sync_sup(ds,time_table,tnext,switch)
    else
       exit
    endif
 enddo

 !print*,"integrate : ",step_count_int,step_count_tsyn,tcoord,tnext,ds_init

 deallocate(bdata)

end subroutine integrate_to_time


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

subroutine new_ds_sync_sup(ds,time_table,tnext,switch)
 real,    intent(inout) :: ds(:)
 real,    intent(in)    :: time_table(:)
 real,    intent(in)    :: tnext
 integer, intent(in)    :: switch
 integer :: i,k
 real :: tp,dtc,dstmp
 do i=1,ck_size
    k = cck_sorted_id(i)
    if (tnext<time_table(k)) exit
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



subroutine backup_data(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,bdata)
 use part, only: igarg,ikappa
 real,    intent(in)  ::xyzmh_ptmass(:,:),vxyz_ptmass(:,:),bin_info(:,:)
 integer, intent(in)  :: group_info(:,:)
 real,    intent(out) ::bdata(:)
 integer, intent(in)  :: start_id,end_id
 integer :: i,j,k
 j=0
 do k=start_id,end_id
    i = group_info(igarg,k)
    bdata(j*7+1) = xyzmh_ptmass(1,i)
    bdata(j*7+2) = xyzmh_ptmass(2,i)
    bdata(j*7+3) = xyzmh_ptmass(3,i)
    bdata(j*7+4) = vxyz_ptmass(1,i)
    bdata(j*7+5) = vxyz_ptmass(2,i)
    bdata(j*7+6) = vxyz_ptmass(3,i)
    bdata(j*7+7) = bin_info(ikappa,i)
    j = j + 1
 enddo

end subroutine backup_data


subroutine restore_state(start_id,end_id,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,tcoord,t_old,W,W_old,bdata)
 use part, only: igarg,ikappa
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),bin_info(:,:)
 integer, intent(in)    :: group_info(:,:)
 real,    intent(out)   :: tcoord,W
 real,    intent(in)    :: t_old,W_old
 real,    intent(in)    :: bdata(:)
 integer, intent(in)    :: start_id,end_id
 integer :: k,i,j
 j = 0
 do k=start_id,end_id
    i = group_info(igarg,k)
    xyzmh_ptmass(1,i)  = bdata(j*7+1)
    xyzmh_ptmass(2,i)  = bdata(j*7+2)
    xyzmh_ptmass(3,i)  = bdata(j*7+3)
    vxyz_ptmass(1,i)   = bdata(j*7+4)
    vxyz_ptmass(2,i)   = bdata(j*7+5)
    vxyz_ptmass(3,i)   = bdata(j*7+6)
    bin_info(ikappa,i) = bdata(j*7+7)
    j = j + 1
 enddo
 tcoord   = t_old
 W = W_old

end subroutine restore_state


!---------------------------------------
!
! TTL Drift routine for multiples only.
!
!---------------------------------------
subroutine drift_TTL(tcoord,W,h,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,gsize,s_id,e_id)
 use part, only: igarg,icomp,ikap
 use io,   only: fatal
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),bin_info(:,:)
 integer, intent(in)    :: group_info(:,:)
 real,    intent(inout) :: tcoord
 real,    intent(in)    :: h,W
 integer, intent(in)    :: s_id,e_id,gsize
 integer, allocatable :: binstack(:)
 integer :: k,i,compi,n
 real    :: dtd,vcom(3),kappai,kappa1i

 allocate(binstack((gsize/4)+1))
 binstack = 0
 n = 0

 dtd = h/W

 tcoord = tcoord + dtd

 do k=s_id,e_id
    i = group_info(igarg,k)
    compi = group_info(icomp,k)
    if (compi/=i) then ! It's a binary. identify companion and drift binary.
       kappai = bin_info(ikap,i)
       if (kappai >= 1.) then
          kappa1i = 1./kappai
       else
          kappa1i = 1.
          call fatal('subgroup','kappa value bellow 1... something wrong here!',var='kappai',val=kappai)
       endif
       if (any(binstack == i)) cycle! If already treated i will be in binstack
       call get_bin_com(i,compi,xyzmh_ptmass,vxyz_ptmass,vcom)
       n = n + 1 ! stack level
       binstack(n) = compi

       call correct_com_drift(xyzmh_ptmass,vxyz_ptmass,vcom,kappa1i,dtd,i,compi)

    else
       xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + dtd*vxyz_ptmass(1,i)
       xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + dtd*vxyz_ptmass(2,i)
       xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + dtd*vxyz_ptmass(3,i)
    endif
 enddo

 deallocate(binstack)

end subroutine drift_TTL

!---------------------------------------
!
! TTL Kick routine for multiples only.
!
!---------------------------------------
subroutine kick_TTL(h,W,xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,fxyz_ptmass,gtgrad,s_id,e_id)
 use part, only: igarg,ikap,icomp
 use io,   only: fatal
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:)
 real,    intent(inout) :: gtgrad(:,:),bin_info(:,:)
 integer, intent(inout) :: group_info(:,:)
 real,    intent(in)    :: h
 real,    intent(inout) :: W
 integer, intent(in)    :: s_id,e_id
 integer, allocatable :: binstack(:)
 real :: om,dw,dtk,kappa1i,kappai,om_old,vcom(3)
 integer :: i,k,n,gsize,compi


 gsize = (e_id-s_id+1)
 allocate(binstack((gsize/4)+1))
 binstack = 0
 n = 0

 if (h==0.) then
    call binaries_in_multiples(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,&
                               gsize,s_id,e_id)
    call get_force_TTL(xyzmh_ptmass,group_info,bin_info,fxyz_ptmass,gtgrad,om_old,s_id,e_id,.true.)
    call get_kappa(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,gsize,s_id,e_id)
    call get_force_TTL(xyzmh_ptmass,group_info,bin_info,fxyz_ptmass,gtgrad,om,s_id,e_id)
    W = W + (om-om_old) ! correct W after updating kappa...
 else

    call get_force_TTL(xyzmh_ptmass,group_info,bin_info,fxyz_ptmass,gtgrad,om,s_id,e_id)

    dtk = h/om
    do k=s_id,e_id
       i=group_info(igarg,k)
       vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + (0.5*dtk)*fxyz_ptmass(1,i)
       vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + (0.5*dtk)*fxyz_ptmass(2,i)
       vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + (0.5*dtk)*fxyz_ptmass(3,i)
    enddo

    dw = 0.
    do k=s_id,e_id
       i      = group_info(igarg,k)
       compi  = group_info(icomp,k)
       if (i/=compi) then

          kappai = bin_info(ikap,i)
          if (kappai >= 1.) then
             kappa1i = 1./kappai
          else
             kappa1i = 1.
             call fatal('subgroup','kappa value bellow 1... something wrong here!',var='kappai',val=kappai)
          endif
          if (any(binstack == i)) cycle! If already treated i will be in binstack
          call get_bin_com(i,compi,xyzmh_ptmass,vxyz_ptmass,vcom)
          n = n+1 ! stack level
          binstack(n) = compi

          call correct_W_SD(dw,vxyz_ptmass,gtgrad,vcom,kappa1i,i,compi)

       else
          dw = dw + (vxyz_ptmass(1,i)*gtgrad(1,i) + &
                    vxyz_ptmass(2,i)*gtgrad(2,i) + &
                    vxyz_ptmass(3,i)*gtgrad(3,i))
       endif
    enddo

    W = W + dw*dtk

    do k=s_id,e_id
       i=group_info(igarg,k)
       vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + (0.5*dtk)*fxyz_ptmass(1,i)
       vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + (0.5*dtk)*fxyz_ptmass(2,i)
       vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + (0.5*dtk)*fxyz_ptmass(3,i)
    enddo
 endif

 deallocate(binstack)

end subroutine kick_TTL

!--------------------------------------------------------------------------
!
! Compressed and optimized Drift-Kick routine for binaries (group size = 2)
!
!--------------------------------------------------------------------------
subroutine oneStep_bin(tcoord,W,ds,kappa1,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,gtgrad,time_table,i,j)
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:),fxyz_ptmass(:,:),gtgrad(:,:),time_table(:)
 real, intent(in)    :: ds,kappa1
 real, intent(inout) :: tcoord,W
 integer, intent(in) :: i,j
 integer :: k
 real :: dtd,dtk,dvel1(3),dvel2(3),dw,om
 real :: vcom(3)



 do k = 1,ck_size
    dtd = ds*cks(k)/W
    tcoord = tcoord + dtd
    time_table(k) = tcoord

    call get_bin_com(i,j,xyzmh_ptmass,vxyz_ptmass,vcom)

    if (kappa1 < 1.0) then
       call correct_com_drift(xyzmh_ptmass,vxyz_ptmass,vcom,kappa1,dtd,i,j)
    else
       xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + dtd*vxyz_ptmass(1,i)
       xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + dtd*vxyz_ptmass(2,i)
       xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + dtd*vxyz_ptmass(3,i)
       xyzmh_ptmass(1,j) = xyzmh_ptmass(1,j) + dtd*vxyz_ptmass(1,j)
       xyzmh_ptmass(2,j) = xyzmh_ptmass(2,j) + dtd*vxyz_ptmass(2,j)
       xyzmh_ptmass(3,j) = xyzmh_ptmass(3,j) + dtd*vxyz_ptmass(3,j)

    endif

    call get_force_TTL_bin(xyzmh_ptmass,fxyz_ptmass,gtgrad,om,kappa1,i,j)

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

    W = W + dw*dtk*kappa1

    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + dvel1(1)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + dvel1(2)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + dvel1(3)
    vxyz_ptmass(1,j) = vxyz_ptmass(1,j) + dvel2(1)
    vxyz_ptmass(2,j) = vxyz_ptmass(2,j) + dvel2(2)
    vxyz_ptmass(3,j) = vxyz_ptmass(3,j) + dvel2(3)

 enddo


end subroutine oneStep_bin

!------------------------------------------------------------------
!
! SD method alters binary intrinsec motion but conserve CoM motion.
! this routine will update the position of binaries by slowing down
! the internal motion but correct the CoM drift operation...
!
!------------------------------------------------------------------
subroutine correct_com_drift(xyzmh_ptmass,vxyz_ptmass,vcom,kappa1,dtd,i,j)
 real, intent(inout) :: xyzmh_ptmass(:,:)
 real, intent(in)    :: vxyz_ptmass(:,:),vcom(3)
 real, intent(in)    :: kappa1,dtd
 integer, intent(in) :: i,j
 real :: vrel(3)


 vrel(1) = vxyz_ptmass(1,i) - vcom(1)
 vrel(2) = vxyz_ptmass(2,i) - vcom(2)
 vrel(3) = vxyz_ptmass(3,i) - vcom(3)

 xyzmh_ptmass(1,i) = xyzmh_ptmass(1,i) + dtd*vrel(1)*kappa1 + vcom(1)*dtd
 xyzmh_ptmass(2,i) = xyzmh_ptmass(2,i) + dtd*vrel(2)*kappa1 + vcom(2)*dtd
 xyzmh_ptmass(3,i) = xyzmh_ptmass(3,i) + dtd*vrel(3)*kappa1 + vcom(3)*dtd

 vrel(1) = vxyz_ptmass(1,j) - vcom(1)
 vrel(2) = vxyz_ptmass(2,j) - vcom(2)
 vrel(3) = vxyz_ptmass(3,j) - vcom(3)

 xyzmh_ptmass(1,j) = xyzmh_ptmass(1,j) + dtd*vrel(1)*kappa1 + vcom(1)*dtd
 xyzmh_ptmass(2,j) = xyzmh_ptmass(2,j) + dtd*vrel(2)*kappa1 + vcom(2)*dtd
 xyzmh_ptmass(3,j) = xyzmh_ptmass(3,j) + dtd*vrel(3)*kappa1 + vcom(3)*dtd



end subroutine correct_com_drift


!------------------------------------------------------------------
!
! Correction method to compute the new value of W when SD on
!
!------------------------------------------------------------------
subroutine correct_W_SD(dW,vxyz_ptmass,gtgrad,vcom,kappa1,i,j)
 real, intent(inout) :: dW
 real, intent(in)    :: vxyz_ptmass(:,:),gtgrad(:,:),vcom(3)
 real, intent(in)    :: kappa1
 integer, intent(in) :: i,j
 real :: vrel(3)


 vrel(1) = vxyz_ptmass(1,i) - vcom(1)
 vrel(2) = vxyz_ptmass(2,i) - vcom(2)
 vrel(3) = vxyz_ptmass(3,i) - vcom(3)

 dW = dW + (vrel(1)*kappa1 + vcom(1))*gtgrad(1,i)
 dW = dW + (vrel(2)*kappa1 + vcom(2))*gtgrad(2,i)
 dW = dW + (vrel(3)*kappa1 + vcom(3))*gtgrad(3,i)

 vrel(1) = vxyz_ptmass(1,j) - vcom(1)
 vrel(2) = vxyz_ptmass(2,j) - vcom(2)
 vrel(3) = vxyz_ptmass(3,j) - vcom(3)

 dW = dW + (vrel(1)*kappa1 + vcom(1))*gtgrad(1,j)
 dW = dW + (vrel(2)*kappa1 + vcom(2))*gtgrad(2,j)
 dW = dW + (vrel(3)*kappa1 + vcom(3))*gtgrad(3,j)



end subroutine correct_W_SD

!---------------------------------------
!
! TTL Force routine for multiples only.
! Potential and initial time step are
! computed here as well.
!
!---------------------------------------
subroutine get_force_TTL(xyzmh_ptmass,group_info,bin_info,fxyz_ptmass,gtgrad,om,s_id,e_id,potonly,energ,ds_init)
 use part, only: igarg,ikap,icomp,isemi
 use io,   only: fatal
 real,              intent(in)    :: xyzmh_ptmass(:,:)
 real,              intent(inout) :: fxyz_ptmass(:,:),gtgrad(:,:),bin_info(:,:)
 integer,           intent(in)    :: group_info(:,:)
 real,              intent(out)   :: om
 integer,           intent(in)    :: s_id,e_id
 logical, optional, intent(in)    :: potonly
 logical, optional, intent(in)    :: energ
 real,    optional, intent(out)   :: ds_init
 real    :: mi,mj,xi,yi,zi,dx,dy,dz,r2,ddr,ddr3,dsi,mcomp,semii
 real    :: gravf,gtk,gtki,gravfi(3),gtgradi(3),kappa1i,kappai
 integer :: i,j,k,l,compi
 logical :: init
 om = 0.


 if (present(ds_init)) then
    init = .true.
    ds_init = huge(om)
 else
    init = .false.
 endif


 do k=s_id,e_id
    i         = group_info(igarg,k)
    compi     = group_info(icomp,k)
    kappai = bin_info(ikap,i)
    if (kappai >= 1.) then
       kappa1i = 1./kappai
    else
       kappa1i = 1.
       call fatal('subgroup','kappa value bellow 1... something wrong here!',var='kappai',val=kappai)
    endif
    gravfi(1) = 0.
    gravfi(2) = 0.
    gravfi(3) = 0.
    gtgradi(1) = 0.
    gtgradi(2) = 0.
    gtgradi(3) = 0.
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
       ddr  = 1./sqrt(r2)
       mj = xyzmh_ptmass(4,j)
       if (j == compi) then
          if (present(potonly) .and. present(energ)) then
             gtk = mj*ddr
          else
             gtk = mj*ddr*kappa1i
          endif
       else
          gtk = mj*ddr
       endif
       gtki  = gtki + gtk
       if (.not.present(potonly)) then
          ddr3 = ddr*ddr*ddr
          if (j == compi) then
             gravf = -mj*ddr3*kappa1i
          else
             gravf = -mj*ddr3
          endif
          gravfi(1) = gravfi(1) + dx*gravf
          gravfi(2) = gravfi(2) + dy*gravf
          gravfi(3) = gravfi(3) + dz*gravf
          gtgradi(1) = gtgradi(1) + dx*gravf*mi
          gtgradi(2) = gtgradi(2) + dy*gravf*mi
          gtgradi(3) = gtgradi(3) + dz*gravf*mi
       endif
    enddo
    if (.not.present(potonly)) then
       fxyz_ptmass(4,i) = -gtki
       fxyz_ptmass(1,i) = gravfi(1)
       fxyz_ptmass(2,i) = gravfi(2)
       fxyz_ptmass(3,i) = gravfi(3)
       gtgrad(1,i) = gtgradi(1)
       gtgrad(2,i) = gtgradi(2)
       gtgrad(3,i) = gtgradi(3)
    endif

    if (init) then
       if (compi /=i) then
          semii = bin_info(isemi,i)
          mcomp = xyzmh_ptmass(4,compi)
          if (semii >= 0) then
             dsi = mi*mcomp*sqrt(semii/(mi+mcomp))*elli_res
          else
             dsi = mi*mcomp*sqrt(-semii/(mi+mcomp))*hyper_res
          endif
          ds_init = min(ds_init,dsi)
       endif
    endif
    om = om + gtki*mi
 enddo

 om = om*0.5

end subroutine get_force_TTL

!--------------------------------------------------------
!
! routine that compute the slowing down factor depending
! on outside pertubartions for multiples only
!
!--------------------------------------------------------
subroutine get_kappa(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,gsize,s_id,e_id)
 use part,         only:igarg,icomp,ipert,ikap,iapo,iecc,iorb,isemi,ipertg
 use utils_kepler, only:extract_a,extract_e
 use dim ,         only:use_sinktree
 real   , intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real   , intent(inout) :: bin_info(:,:)
 integer, intent(in)    :: group_info(:,:)
 integer, intent(in)    :: s_id,e_id,gsize
 integer, allocatable :: binstack(:)
 integer :: k,l,i,j,compi,n
 real :: pouti,r2,v2,drdv,dr(3),dv(3),ddr,ddr3,r,m1,m2,mj,mui,muij
 real :: kappa,kappa_max,rapo,rapo3,Ti,semij
 real :: vcom(3),xcom(3),eij,e2,rvrm,vr_max2,ronvr,timescale
 allocate(binstack(gsize))
 binstack = 0
 n = 0
 timescale = huge(timescale)

 do k=s_id,e_id
    i     = group_info(igarg,k)
    compi = group_info(icomp,k)
    if (compi == i .or. bin_info(isemi,i)<0.) then
       bin_info(ikap,i) = 1.
    else
       if (any(binstack == i)) cycle
       call get_bin_com(i,compi,xyzmh_ptmass,vxyz_ptmass,vcom,xcom)
       n = n+1 ! level of the stack
       binstack(n) = compi
       pouti = bin_info(ipert,i)
       if (use_sinktree) pouti = pouti + bin_info(ipertg,i)
       Ti    = bin_info(iorb,i)
       m1 = xyzmh_ptmass(4,i)
       m2 = xyzmh_ptmass(4,compi)

       do l=s_id,e_id
          if (k == l) cycle
          j = group_info(igarg,l)
          if (j == compi) cycle
          mj   = xyzmh_ptmass(4,j)

          dr(1)   = xcom(1) - xyzmh_ptmass(1,j)
          dr(2)   = xcom(2) - xyzmh_ptmass(2,j)
          dr(3)   = xcom(3) - xyzmh_ptmass(3,j)
          dv(1)   = vcom(1) - vxyz_ptmass(1,j)
          dv(2)   = vcom(2) - vxyz_ptmass(2,j)
          dv(3)   = vcom(3) - vxyz_ptmass(3,j)
          r2   = dr(1)**2+dr(2)**2+dr(3)**2
          v2   = dv(1)**2+dv(2)**2+dv(3)**2
          drdv = dr(1)*dv(1)+dr(2)*dv(2)+dr(3)*dv(3)
          ddr  = 1./sqrt(r2)
          r    = 1/ddr
          ddr3 = ddr*ddr*ddr
          pouti = pouti + mj*ddr3
          muij = m1 + m2 + mj

          call extract_a(r,v2,muij,semij)
          if (semij<0.) then
             timescale = min(timescale, r2/v2)

          else
             call extract_e(dr(1),dr(2),dr(3),dv(1),dv(2),dv(3),muij,r,eij)
             e2   = eij**2
             rvrm = semij*(1-e2)
             if (r < rvrm) then
                vr_max2 = e2*muij/rvrm
                timescale = min(timescale, semij*semij/vr_max2)
             endif
             ronvr = r2/abs(drdv)
             timescale = min(timescale, ronvr**2)
          endif

       enddo

       mui    = (m1*m2)/(m1+m2)
       rapo  = bin_info(iapo,i)
       rapo3 = rapo*rapo*rapo
       timescale = sqrt(timescale)

       kappa_max = max(0.001*timescale/Ti,1.0)
       kappa     = kref/((rapo3/mui)*pouti)
       kappa     = min(kappa_max,kappa)
       kappa     = max(1.0,kappa)
       bin_info(ikap,i)     = kappa
       bin_info(ikap,compi) = kappa

    endif
 enddo

 deallocate(binstack)

end subroutine get_kappa


!---------------------------------------
!
! TTL Force routine for binaries only.
! Potential and initial time step are
! computed here as well.
!
!---------------------------------------
subroutine get_force_TTL_bin(xyzmh_ptmass,fxyz_ptmass,gtgrad,om,kappa1,i,j,potonly,ds_init,semiij)
 real,    intent(in)    :: xyzmh_ptmass(:,:)
 real,    intent(inout) :: fxyz_ptmass(:,:),gtgrad(:,:)
 integer, intent(in)    :: i,j
 real,    intent(in)    :: kappa1
 real,    intent(out)   :: om
 logical, optional, intent(in)    :: potonly
 real,    optional, intent(out)   :: ds_init
 real,    optional, intent(in)    :: semiij
 real :: dx,dy,dz,r2,ddr,ddr3,mi,mj
 real :: gravfi,gravfj,gtki,gtkj,fxi,fyi,fzi,fxj,fyj,fzj

 mi = xyzmh_ptmass(4,i)
 mj = xyzmh_ptmass(4,j)
 dx = xyzmh_ptmass(1,i) - xyzmh_ptmass(1,j)
 dy = xyzmh_ptmass(2,i) - xyzmh_ptmass(2,j)
 dz = xyzmh_ptmass(3,i) - xyzmh_ptmass(3,j)
 r2 = dx**2+dy**2+dz**2
 ddr  = 1./sqrt(r2)
 ddr3 = ddr*ddr*ddr

 if (kappa1<1.0 .and. .not.present(potonly)) then
    gravfi = kappa1*mj*ddr3
    gravfj = kappa1*mi*ddr3
    gtki   = kappa1*mj*ddr
    gtkj   = kappa1*mi*ddr
 else
    gravfi = mj*ddr3
    gravfj = mi*ddr3
    gtki   = mj*ddr
    gtkj   = mi*ddr
 endif



 if (.not.present(potonly)) then
    fxi = -dx*gravfi
    fyi = -dy*gravfi
    fzi = -dz*gravfi
    fxj =  dx*gravfj
    fyj =  dy*gravfj
    fzj =  dz*gravfj
    fxyz_ptmass(4,i) = -gtki
    fxyz_ptmass(4,j) = -gtkj
    fxyz_ptmass(1,i) = fxi
    fxyz_ptmass(2,i) = fyi
    fxyz_ptmass(3,i) = fzi
    fxyz_ptmass(1,j) = fxj
    fxyz_ptmass(2,j) = fyj
    fxyz_ptmass(3,j) = fzj
    gtgrad(1,i) = -dx*gravfi*mi
    gtgrad(2,i) = -dy*gravfi*mi
    gtgrad(3,i) = -dz*gravfi*mi
    gtgrad(1,j) = dx*gravfj*mj
    gtgrad(2,j) = dy*gravfj*mj
    gtgrad(3,j) = dz*gravfj*mj
 endif

 om = gtki*mi

 if (present(ds_init) .and. .not.present(potonly)) then
    if (semiij >= 0) then
       ds_init = mi*mj*sqrt(semiij/(mi+mj))*elli_res
    else
       ds_init = mi*mj*sqrt(-semiij/(mi+mj))*hyper_res
    endif
 endif


end subroutine get_force_TTL_bin


!--------------------------------------------------------
!
! routine that compute the slowing down factor depending
! on outside pertubartions for binaries only
!
!--------------------------------------------------------
subroutine get_kappa_bin(xyzmh_ptmass,bin_info,i,j)
 use part, only:ipert,iapo,ikap,isemi,iecc,ipertg
 use dim , only:use_sinktree
 real, intent(inout) :: bin_info(:,:)
 real, intent(in)    :: xyzmh_ptmass(:,:)
 integer, intent(in) :: i,j
 real    :: kappa,m1,m2,pert,mu,rapo,rapo3
 logical :: isellip


 isellip = bin_info(isemi,i) > 0.
 m1 = xyzmh_ptmass(4,i)
 m2 = xyzmh_ptmass(4,j)
 mu = (m1*m2)/(m1+m2)
 pert = bin_info(ipert,i)
 if (use_sinktree) pert = pert + bin_info(ipertg,i)
 rapo = bin_info(iapo,i)
 rapo3 = rapo*rapo*rapo
 kappa = kref/((rapo3/mu)*pert)
 !print*,xyzmh_ptmass(2,i),pert,kappa,rapo,bin_info(isemi,i),bin_info(iecc,i)
 if (kappa > 1. .and. isellip) then
    bin_info(ikap,i) = kappa
    bin_info(ikap,j) = kappa
 else
    bin_info(ikap,i) = 1.
    bin_info(ikap,j) = 1.
 endif

end subroutine get_kappa_bin

!--------------------------------------------------------
!
! interface routine that update the slowing down factor
! for each subgroups in the simulation
!
!--------------------------------------------------------
subroutine update_kappa(xyzmh_ptmass,vxyz_ptmass,bin_info,group_info,n_group)
 use part, only:igcum,igarg
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(inout) :: bin_info(:,:)
 integer, intent(in)    :: group_info(:,:)
 integer, intent(in)    :: n_group
 integer :: i,start_id,end_id,prim,sec,gsize

 do i=1,n_group
    start_id = group_info(igcum,i) + 1
    end_id   = group_info(igcum,i+1)
    gsize    = (end_id - start_id) + 1
    if (gsize>2) then
       call get_kappa(xyzmh_ptmass,vxyz_ptmass,group_info,bin_info,gsize,start_id,end_id)
    else
       prim = group_info(igarg,start_id)
       sec = group_info(igarg,end_id)
       call get_kappa_bin(xyzmh_ptmass,bin_info,prim,sec)
    endif
 enddo
end subroutine update_kappa

subroutine get_bin_com(i,j,xyzmh_ptmass,vxyz_ptmass,vcom,xcom)
 real,    intent(in)        :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(out)       :: vcom(3)
 integer, intent(in)        :: i,j
 real, intent(out),optional :: xcom(3)
 real :: mtot,m1,m2

 m1 = xyzmh_ptmass(4,i)
 m2 = xyzmh_ptmass(4,j)
 mtot = m1 + m2

 vcom(1) = (m1*vxyz_ptmass(1,i)+m2*vxyz_ptmass(1,j))/mtot
 vcom(2) = (m1*vxyz_ptmass(2,i)+m2*vxyz_ptmass(2,j))/mtot
 vcom(3) = (m1*vxyz_ptmass(3,i)+m2*vxyz_ptmass(3,j))/mtot

 if (present(xcom)) then
    xcom(1) = (m1*xyzmh_ptmass(1,i)+m2*xyzmh_ptmass(1,j))/mtot
    xcom(2) = (m1*xyzmh_ptmass(2,i)+m2*xyzmh_ptmass(2,j))/mtot
    xcom(3) = (m1*xyzmh_ptmass(3,i)+m2*xyzmh_ptmass(3,j))/mtot
 endif

end subroutine get_bin_com

!--------------------------------------------------------
!
! Routine to compute potential energy in subgroups
!
!--------------------------------------------------------
subroutine get_pot_subsys(n_group,group_info,bin_info,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,&
                          gtgrad,epot_sinksink)
 use part,     only:igarg,igcum,ikap
 use io,       only:id,master,fatal
 use mpiutils, only:bcast_mpi
 integer, intent(in)    :: n_group
 real,    intent(inout) :: xyzmh_ptmass(:,:),fxyz_ptmass(:,:),gtgrad(:,:)
 real,    intent(inout) :: bin_info(:,:),vxyz_ptmass(:,:)
 integer, intent(in)    :: group_info(:,:)
 real,    intent(inout) :: epot_sinksink
 integer :: i,start_id,end_id,gsize,prim,sec
 real :: phitot,phigroup,kappa1

 phitot = 0.


 if (n_group>0) then
    call update_kappa(xyzmh_ptmass,vxyz_ptmass,bin_info,group_info,n_group)
    if (id==master) then

       !$omp parallel do default(none)&
       !$omp shared(xyzmh_ptmass,fxyz_ptmass)&
       !$omp shared(group_info,gtgrad,n_group,bin_info)&
       !$omp private(i,start_id,end_id,gsize,prim,sec,phigroup,kappa1)&
       !$omp reduction(+:phitot)
       do i=1,n_group
          start_id = group_info(igcum,i) + 1
          end_id   = group_info(igcum,i+1)
          gsize    = (end_id - start_id) + 1
          if (gsize>2) then
             call get_force_TTL(xyzmh_ptmass,group_info,bin_info,fxyz_ptmass,gtgrad,phigroup,start_id,end_id,.true.,.true.)
          else
             prim = group_info(igarg,start_id)
             sec = group_info(igarg,end_id)
             if (bin_info(ikap,prim) >= 1.) then
                kappa1 = 1./bin_info(ikap,prim)
             else
                kappa1 = 1.
                call fatal('subgroup','kappa value bellow 1... something wrong here!(energy)',var='kappa',val=bin_info(ikap,prim))
             endif
             call get_force_TTL_bin(xyzmh_ptmass,fxyz_ptmass,gtgrad,phigroup,kappa1,prim,sec,.true.)
          endif
          phitot = phitot + phigroup
       enddo
       !$omp end parallel do
    endif
 endif

 epot_sinksink = epot_sinksink - phitot
 call bcast_mpi(epot_sinksink) ! broadcast to other MPI threads

end subroutine get_pot_subsys


end module subgroup
