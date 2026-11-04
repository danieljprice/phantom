!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpitree
!
! This module handles the MPI exchange of tree nodes between tasks
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dtypekdtree, io, mpi
!

#ifdef MPI
 use mpi
#endif
 use io, only:nprocs

 implicit none

 public :: allgatherv_nodes

 private

contains

!----------------------------------------------------------------
!+
!  gather blocks of tree nodes of different sizes from all tasks
!  (as allgatherv_mpi in mpiutils, for the kdnode type)
!+
!----------------------------------------------------------------
subroutine allgatherv_nodes(node_in,node_out,counts,displs)
 use dtypekdtree, only:get_mpitype_of_kdnode,kdnode
 type(kdnode), intent(in)  :: node_in(:)
 type(kdnode), intent(out) :: node_out(:)
 integer,      intent(in)  :: counts(nprocs),displs(nprocs)
#ifdef MPI
 integer :: dtype_kdnode,mpierr

 call get_mpitype_of_kdnode(dtype_kdnode)
 call MPI_ALLGATHERV(node_in,size(node_in),dtype_kdnode,node_out,counts,displs,dtype_kdnode,&
                     MPI_COMM_WORLD,mpierr)
#else
 node_out(1:size(node_in)) = node_in
#endif

end subroutine allgatherv_nodes

end module mpitree
