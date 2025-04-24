!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpitree
!
! This module handles the MPI exchange of information during the kdtree routines
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, dim, dtypekdtree, io, mpi
!

#ifdef MPI
 use mpi
#endif
 use io,             only:id,nprocs

 implicit none

 interface reduce_group
  module procedure reduce_group_real,reduce_group_int
 end interface reduce_group

 public :: get_group_cofm
 public :: reduce_group

 public :: tree_sync
 public :: tree_bcast
 public :: finish_tree_comms

#ifdef MPI
 integer :: globallevel
#endif

 integer,allocatable :: comm_cofm(:)  ! only comms up to globallevel are used
 integer,allocatable :: comm_owner(:) ! only comms up to globallevel are used

contains

subroutine allocate_tree_comms_arrays
 use allocutils, only:allocate_array
 use dim,        only:mpi
 if (mpi) then
    call allocate_array('comm_cofm',     comm_cofm,     nprocs)
    call allocate_array('comm_owner',    comm_owner,    nprocs)
    call init_tree_comms
 endif
end subroutine allocate_tree_comms_arrays

subroutine deallocate_tree_comms_arrays
 if (allocated(comm_cofm    )) deallocate(comm_cofm    )
 if (allocated(comm_owner   )) deallocate(comm_owner   )
end subroutine deallocate_tree_comms_arrays

!----------------------------------------------------------------
!+
!  initialise communicators for tree construction
!+
!----------------------------------------------------------------
subroutine init_tree_comms()
#ifdef MPI
 integer :: level,groupsize,color
 integer                            :: mpierr

 globallevel = int(ceiling(log(real(nprocs)) / log(2.0)))

 ! cofm group is all procs at level 0
 call MPI_COMM_DUP(MPI_COMM_WORLD, comm_cofm(1), mpierr)
 ! owner comm is never needed at level 0, left uninitialised

 do level = 1, globallevel
    groupsize = 2**(globallevel - level)

    ! cofm group
    color = id / groupsize
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, id, comm_cofm(level+1), mpierr)

    ! owner group
    if (groupsize == 1) then
       call MPI_COMM_DUP(MPI_COMM_WORLD, comm_owner(level+1), mpierr)
    else
       if (mod(id, groupsize) == 0) then
          color = 0
       else
          color = 1
       endif
       call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, id, comm_owner(level+1), mpierr)
    endif
 enddo
#endif
end subroutine init_tree_comms

subroutine finish_tree_comms()
#ifdef MPI
 integer :: level
 integer :: mpierr

 do level = 0, globallevel
    call MPI_COMM_FREE(comm_cofm(level+1), mpierr)
 enddo
 do level = 1, globallevel
    call MPI_COMM_FREE(comm_owner(level+1), mpierr)
 enddo
#endif
end subroutine finish_tree_comms

!----------------------------------------------------------------
!+
!  get the COFM of the group
!+
!----------------------------------------------------------------
subroutine get_group_cofm(xyzcofm,totmass_node,level,cofmsum,totmassg)
 real,      intent(in)        :: xyzcofm(3)
 real,      intent(in)        :: totmass_node
 integer,   intent(in)        :: level

 real,      intent(out)       :: cofmsum(3)
 real,      intent(out)       :: totmassg

#ifdef MPI
 real                         :: cofmpart(3)
 integer                      :: mpierr

 cofmpart = xyzcofm * totmass_node
 call MPI_ALLREDUCE(totmass_node,totmassg,1,MPI_REAL8,MPI_SUM,comm_cofm(level+1),mpierr)
 call MPI_ALLREDUCE(cofmpart,cofmsum,3,MPI_REAL8,MPI_SUM,comm_cofm(level+1),mpierr)
 cofmsum = cofmsum / totmassg
#else
 cofmsum = xyzcofm*totmass_node
 totmassg = totmass_node
#endif

end subroutine get_group_cofm

!----------------------------------------------------------------
!+
!  tree group reductions
!+
!----------------------------------------------------------------

function reduce_group_real(x,string,level) result(xg)
 use io, only:fatal
 real,               intent(in)        :: x
 character(len=*),   intent(in)        :: string
 integer,            intent(in)        :: level
 real                                  :: xg

#ifdef MPI
 real                                  :: isend, ired
 integer                               :: mpierr

 isend = x

 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,1,MPI_REAL8,MPI_SUM,comm_cofm(level+1),mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,1,MPI_REAL8,MPI_MAX,comm_cofm(level+1),mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,1,MPI_REAL8,MPI_MIN,comm_cofm(level+1),mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select

 xg = ired
#else
 xg = x
#endif

end function reduce_group_real

function reduce_group_int(x,string,level) result(xg)
 use io, only:fatal
 integer,            intent(in)        :: x
 character(len=*),   intent(in)        :: string
 integer,            intent(in)        :: level
 integer                               :: xg

#ifdef MPI
 integer                               :: isend, ired
 integer                               :: mpierr

 isend = x

 select case(trim(string))
 case('+')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_SUM,comm_cofm(level+1),mpierr)
 case('max')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_MAX,comm_cofm(level+1),mpierr)
 case('min')
    call MPI_ALLREDUCE(isend,ired,1,MPI_INTEGER,MPI_MIN,comm_cofm(level+1),mpierr)
 case default
    call fatal('reduceall (mpi)','unknown reduction operation')
 end select

 xg = ired
#else
 xg = x
#endif
end function reduce_group_int

!----------------------------------------------------------------
!+
!  synchronize the global tree, placing nodes in the correct position
!+
!----------------------------------------------------------------
subroutine tree_sync(node_in,n_in,node_synced,n_synced,ifirstingroup,level)
 use dtypekdtree, only:get_mpitype_of_kdnode,kdnode

 integer, intent(in)         :: ifirstingroup,level
 integer, intent(in)         :: n_in      ! nodes sent per proc
 integer, intent(in)         :: n_synced  ! nodes in the synchronised array
 type(kdnode), intent(in)    :: node_in(n_in)
 type(kdnode), intent(inout) :: node_synced(n_synced)

#ifdef MPI
 integer                     :: dtype_kdnode
 integer                     :: mpierr

!  If there is only 1 owner, do a direct copy
 if (n_in == n_synced) then
    node_synced(:) = node_in(:)
 else
    call get_mpitype_of_kdnode(dtype_kdnode)
    ! skip if we are not an owner
    if (id == ifirstingroup) then
       ! perform node exchange
       call MPI_ALLGATHER(node_in,n_in,dtype_kdnode, &
                          node_synced,n_in,dtype_kdnode, &
                          comm_owner(level+1),mpierr)
    endif
 endif
#endif

end subroutine tree_sync

!----------------------------------------------------------------
!+
!  broadcast the tree from owners to non-owners
!+
!----------------------------------------------------------------
subroutine tree_bcast(node, nnode, level)
 use dtypekdtree, only:get_mpitype_of_kdnode,kdnode

 integer,      intent(in)        :: nnode
 type(kdnode), intent(inout)     :: node(nnode)
 integer,      intent(in)        :: level

#ifdef MPI
 integer                         :: dtype_kdnode
 integer                         :: mpierr

 call get_mpitype_of_kdnode(dtype_kdnode)
 ! the bcast root is relative to the communicator (i.e. it needs to be 0, not ifirstingroup)
 call MPI_BCAST(node, nnode, dtype_kdnode, 0, comm_cofm(level+1), mpierr)
#endif

end subroutine tree_bcast

end module mpitree
