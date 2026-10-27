!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpiradiation
!
! Remote neighbours for the implicit radiation solver under MPI.
!
! Before each implicit solve, every active leaf with neighbours on other
! tasks is sent to those tasks. Each task returns its particles that are
! neighbours of the leaf's solver particles (either particle inside the
! kernel of the other), with the data that stay fixed during the solve.
! These become remote particles, numbered npart+1 ... npart+nremote in the
! solver's arrays. After each flux calculation, the owners send the values
! that change during the solve (E, opacity and flux limiter) for the same
! particles, in the same order.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, io, kdtree, mpi, mpiutils, neighcompact,
!   neighkdtree, part
!
#ifdef MPI
 use mpi
#endif
 implicit none

 !--number of remote particles, and their positions, masses, densities and
 !  kernel gradient corrections (zeta*omega)
 integer, public :: nremote = 0
 real, allocatable, public :: xyzh_remote(:,:),pmass_remote(:),rho_remote(:),zetaomega_remote(:)
 !--remote particles that may be neighbours of each local leaf (compressed rows):
 !  those of leaf icell are icand(icand_start(icell):icand_start(icell+1)-1)
 integer, allocatable, public :: icand_start(:),icand(:)
 !--local particles sent to other tasks as remote particles, in the order sent
 integer, public :: nsendlist = 0
 integer, allocatable, public :: isendlist(:)

 public :: setup_remote_neighbours,exchange_remote_values

 private

#ifdef MPI
 integer, parameter :: nvals    = 3     ! values exchanged each sweep: E, opacity, flux limiter
 integer, parameter :: nstatic  = 7     ! x, y, z, h, mass, density, zeta*omega
 integer, parameter :: nleafhdr = 6     ! particle count, centre, size and search radius of a leaf
 integer, parameter :: maxcellcache = 10000
 integer, parameter :: itag_remote = 7301

 !--particles sent to and received from each task: isendlist(isendoff(iproc)+1:
 !  isendoff(iproc)+nsendto(iproc)) go to task iproc-1, and the remote particles
 !  irecvoff(iproc)+1 ... irecvoff(iproc)+nrecvfrom(iproc) come from it
 integer, allocatable :: nsendto(:),isendoff(:),nrecvfrom(:),irecvoff(:)
 real,    allocatable :: sendvals(:),recvvals(:)

 !--list of integers, used for the neighbours found for each requested leaf
 type intlist
    integer, allocatable :: i(:)
 end type intlist
#endif

contains

!-----------------------------------------------------------------------
!+
!  find the remote neighbours of the local solver particles for the
!  implicit solve, and the local particles to send to other tasks
!  (must be called by every task, after the tree has been built)
!+
!-----------------------------------------------------------------------
subroutine setup_remote_neighbours(npart,xyzh)
 use io,          only:nprocs
 use neighkdtree, only:ncells
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
#ifdef MPI
 integer, allocatable :: iexport(:,:),nleafto(:),ileafoff(:)
 real,    allocatable :: reqsend(:),reqrecv(:),repsend(:),reprecv(:)
 integer :: nreq(nprocs),nreqin(nprocs),ireqinoff(nprocs),nrep(nprocs),nrepin(nprocs),irepinoff(nprocs)
#endif

 if (allocated(icand_start)) deallocate(icand_start)
 allocate(icand_start(int(ncells)+1))
 icand_start = 1
 nremote = 0
 nsendlist = 0
 if (nprocs == 1) return
#ifdef MPI
 if (.not.allocated(nsendto)) allocate(nsendto(nprocs),isendoff(nprocs),nrecvfrom(nprocs),irecvoff(nprocs))
 call find_exported_leaves(xyzh,iexport,nleafto,ileafoff)
 call build_requests(xyzh,iexport,nleafto,ileafoff,reqsend,nreq)
 call alltoallv_reals(reqsend,nreq,reqrecv,nreqin,ireqinoff)
 call answer_requests(npart,xyzh,reqrecv,nreqin,ireqinoff,repsend,nrep)
 call alltoallv_reals(repsend,nrep,reprecv,nrepin,irepinoff)
 call read_replies(reprecv,irepinoff,iexport,nleafto,ileafoff)
#endif

end subroutine setup_remote_neighbours

!-----------------------------------------------------------------------
!+
!  send the current E, opacity and flux limiter of the particles in the
!  send lists, and receive those of the remote particles
!  (must be called by every task, from one thread)
!+
!-----------------------------------------------------------------------
subroutine exchange_remote_values(EU0,npart)
#ifdef MPI
 use io,       only:nprocs
 use mpiutils, only:MPI_DEFAULT_REAL
#endif
 integer, intent(in)    :: npart
 real,    intent(inout) :: EU0(:,:)
#ifdef MPI
 integer :: iproc,k,j,nreq,ierr
 integer :: ireq(2*nprocs)

 if (nprocs == 1) return
 !--post the receives
 nreq = 0
 do iproc=1,nprocs
    if (nrecvfrom(iproc) > 0) then
       nreq = nreq + 1
       call MPI_IRECV(recvvals(nvals*irecvoff(iproc)+1),nvals*nrecvfrom(iproc),MPI_DEFAULT_REAL, &
                      iproc-1,itag_remote,MPI_COMM_WORLD,ireq(nreq),ierr)
    endif
 enddo
 !--pack and send
 do k=1,nsendlist
    j = isendlist(k)
    sendvals(nvals*(k-1)+1) = EU0(1,j)
    sendvals(nvals*(k-1)+2) = EU0(4,j)
    sendvals(nvals*(k-1)+3) = EU0(5,j)
 enddo
 do iproc=1,nprocs
    if (nsendto(iproc) > 0) then
       nreq = nreq + 1
       call MPI_ISEND(sendvals(nvals*isendoff(iproc)+1),nvals*nsendto(iproc),MPI_DEFAULT_REAL, &
                      iproc-1,itag_remote,MPI_COMM_WORLD,ireq(nreq),ierr)
    endif
 enddo
 call MPI_WAITALL(nreq,ireq,MPI_STATUSES_IGNORE,ierr)
 !--unpack into the remote entries
 do k=1,nremote
    EU0(1,npart+k) = recvvals(nvals*(k-1)+1)
    EU0(4,npart+k) = recvvals(nvals*(k-1)+2)
    EU0(5,npart+k) = recvvals(nvals*(k-1)+3)
 enddo
#endif

end subroutine exchange_remote_values

#ifdef MPI
!-----------------------------------------------------------------------
!+
!  number of solver particles in a leaf
!+
!-----------------------------------------------------------------------
integer function nsolver_particles(icell,xyzh)
 use kdtree,       only:inodeparts,inoderange
 use neighcompact, only:is_solver_particle
 integer, intent(in) :: icell
 real,    intent(in) :: xyzh(:,:)
 integer :: ip

 nsolver_particles = 0
 do ip=inoderange(1,icell),inoderange(2,icell)
    if (is_solver_particle(inodeparts(ip),xyzh)) nsolver_particles = nsolver_particles + 1
 enddo

end function nsolver_particles

!-----------------------------------------------------------------------
!+
!  whether two particles (x,y,z,h) are neighbours: either particle
!  lies inside the kernel of the other
!+
!-----------------------------------------------------------------------
logical function is_pair(xyzhi,xyzhj)
 use dim,      only:periodic
 use boundary, only:mod_boundary
 use kernel,   only:radkern2
 real, intent(in) :: xyzhi(4),xyzhj(4)
 real :: dx,dy,dz,rij2

 dx = xyzhi(1) - xyzhj(1)
 dy = xyzhi(2) - xyzhj(2)
 dz = xyzhi(3) - xyzhj(3)
 if (periodic) call mod_boundary(dx,dy,dz)
 rij2 = dx*dx + dy*dy + dz*dz
 is_pair = (rij2 < radkern2*xyzhi(4)**2 .or. rij2 < radkern2*xyzhj(4)**2)

end function is_pair

!-----------------------------------------------------------------------
!+
!  start offsets (from zero) from counts
!+
!-----------------------------------------------------------------------
subroutine set_offsets(ncount,ioff)
 integer, intent(in)  :: ncount(:)
 integer, intent(out) :: ioff(:)
 integer :: i

 ioff(1) = 0
 do i=2,size(ncount)
    ioff(i) = ioff(i-1) + ncount(i-1)
 enddo

end subroutine set_offsets

!-----------------------------------------------------------------------
!+
!  the active leaves with solver particles and neighbours on other tasks,
!  as (leaf, task) pairs grouped by task: those for task iproc-1 are
!  entries ileafoff(iproc)+1 ... ileafoff(iproc)+nleafto(iproc) of iexport
!+
!-----------------------------------------------------------------------
subroutine find_exported_leaves(xyzh,iexport,nleafto,ileafoff)
 use io,          only:id,nprocs
 use neighkdtree, only:ncells,leaf_is_active,get_neighbour_list,listneigh
 real,    intent(in)  :: xyzh(:,:)
 integer, allocatable, intent(out) :: iexport(:,:),nleafto(:),ileafoff(:)
 integer, allocatable :: ifound(:,:)
 real,    allocatable, save :: xyzcache(:,:)
 !$omp threadprivate(xyzcache)
 logical :: remote_export(nprocs)
 integer :: icell,iproc,nneigh,n,nfound

 nfound = 0
 allocate(ifound(2,1024))
 !$omp parallel default(none) shared(ncells,leaf_is_active,xyzh,nprocs,id,nfound,ifound) &
 !$omp private(icell,iproc,nneigh,remote_export)
 if (.not.allocated(xyzcache)) allocate(xyzcache(4,maxcellcache))
 !$omp do schedule(runtime)
 do icell=1,int(ncells)
    if (leaf_is_active(icell) <= 0) cycle
    if (nsolver_particles(icell,xyzh) == 0) cycle
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true., &
                            remote_export=remote_export)
    remote_export(id+1) = .false.
    do iproc=1,nprocs
       if (.not.remote_export(iproc)) cycle
       !$omp critical (radremote_export)
       if (nfound >= size(ifound,2)) call grow_pairs(ifound)
       nfound = nfound + 1
       ifound(1,nfound) = icell
       ifound(2,nfound) = iproc
       !$omp end critical (radremote_export)
    enddo
 enddo
 !$omp enddo
 !$omp end parallel

 !--group by task
 allocate(nleafto(nprocs),ileafoff(nprocs),iexport(2,nfound))
 nleafto = 0
 do n=1,nfound
    nleafto(ifound(2,n)) = nleafto(ifound(2,n)) + 1
 enddo
 call set_offsets(nleafto,ileafoff)
 nleafto = 0
 do n=1,nfound
    iproc = ifound(2,n)
    nleafto(iproc) = nleafto(iproc) + 1
    iexport(:,ileafoff(iproc)+nleafto(iproc)) = ifound(:,n)
 enddo

end subroutine find_exported_leaves

!-----------------------------------------------------------------------
!+
!  double the size of a list of (leaf, task) pairs
!+
!-----------------------------------------------------------------------
subroutine grow_pairs(ipairs)
 integer, allocatable, intent(inout) :: ipairs(:,:)
 integer, allocatable :: itmp(:,:)

 allocate(itmp(2,2*size(ipairs,2)))
 itmp(:,1:size(ipairs,2)) = ipairs
 call move_alloc(itmp,ipairs)

end subroutine grow_pairs

!-----------------------------------------------------------------------
!+
!  pack the requests: for each exported leaf, the number of solver
!  particles, the leaf centre, size and search radius, then x, y, z, h
!  of its solver particles
!+
!-----------------------------------------------------------------------
subroutine build_requests(xyzh,iexport,nleafto,ileafoff,reqsend,nreq)
 use io,           only:nprocs
 use neighkdtree,  only:get_cell_location
 use kdtree,       only:inodeparts,inoderange
 use neighcompact, only:is_solver_particle
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(in)  :: iexport(:,:),nleafto(:),ileafoff(:)
 real, allocatable, intent(out) :: reqsend(:)
 integer, intent(out) :: nreq(:)
 integer :: iproc,l,icell,ip,i,ioff,ihdr,np
 real    :: xpos(3),xsizei,rcuti

 do iproc=1,nprocs
    nreq(iproc) = 0
    do l=ileafoff(iproc)+1,ileafoff(iproc)+nleafto(iproc)
       nreq(iproc) = nreq(iproc) + nleafhdr + 4*nsolver_particles(iexport(1,l),xyzh)
    enddo
 enddo
 allocate(reqsend(max(sum(nreq),1)))

 ioff = 0
 do iproc=1,nprocs
    do l=ileafoff(iproc)+1,ileafoff(iproc)+nleafto(iproc)
       icell = iexport(1,l)
       ihdr  = ioff
       ioff  = ioff + nleafhdr
       np = 0
       do ip=inoderange(1,icell),inoderange(2,icell)
          i = inodeparts(ip)
          if (.not.is_solver_particle(i,xyzh)) cycle
          np = np + 1
          reqsend(ioff+1:ioff+4) = xyzh(1:4,i)
          ioff = ioff + 4
       enddo
       call get_cell_location(icell,xpos,xsizei,rcuti)
       reqsend(ihdr+1)        = real(np)
       reqsend(ihdr+2:ihdr+4) = xpos
       reqsend(ihdr+5)        = xsizei
       reqsend(ihdr+6)        = rcuti
    enddo
 enddo

end subroutine build_requests

!-----------------------------------------------------------------------
!+
!  answer the requests from other tasks: find the local particles that
!  are neighbours of the solver particles of each leaf. The reply to each
!  task is the number of particles sent, their static data (each particle
!  once), then for each leaf the number of its neighbours and their
!  positions in that list
!+
!-----------------------------------------------------------------------
subroutine answer_requests(npart,xyzh,reqrecv,nreqin,ireqinoff,repsend,nrep)
 use io,   only:nprocs
 use part, only:rho,gradh,get_pmass
 use dim,  only:igradzeta,igradomega
 integer, intent(in)  :: npart,nreqin(:),ireqinoff(:)
 real,    intent(in)  :: xyzh(:,:),reqrecv(:)
 real, allocatable, intent(out) :: repsend(:)
 integer, intent(out) :: nrep(:)
 type(intlist), allocatable :: found(:)
 integer, allocatable :: ileafstart(:),ileafproc(:),imarkproc(:),imarkpos(:)
 integer :: nleafin,iproc,ioff,l,k,j,nfound,nlist,ipos,ipos0,ntotfound

 !--locate the leaves in the requests
 nleafin = 0
 do iproc=1,nprocs
    ioff = ireqinoff(iproc)
    do while (ioff < ireqinoff(iproc)+nreqin(iproc))
       nleafin = nleafin + 1
       ioff = ioff + nleafhdr + 4*nint(reqrecv(ioff+1))
    enddo
 enddo
 allocate(ileafstart(nleafin),ileafproc(nleafin),found(nleafin))
 l = 0
 do iproc=1,nprocs
    ioff = ireqinoff(iproc)
    do while (ioff < ireqinoff(iproc)+nreqin(iproc))
       l = l + 1
       ileafstart(l) = ioff
       ileafproc(l)  = iproc
       ioff = ioff + nleafhdr + 4*nint(reqrecv(ioff+1))
    enddo
 enddo

 call find_leaf_neighbours(nleafin,ileafstart,reqrecv,xyzh,found)

 !--the list of particles sent to each task, each particle once per task;
 !  the neighbours of each leaf are replaced by their positions in that list
 ntotfound = 0
 do l=1,nleafin
    ntotfound = ntotfound + size(found(l)%i)
 enddo
 allocate(imarkproc(npart),imarkpos(npart))
 imarkproc = 0
 imarkpos  = 0
 if (allocated(isendlist)) deallocate(isendlist)
 allocate(isendlist(max(ntotfound,1)))
 nsendlist = 0
 do iproc=1,nprocs
    nlist = 0
    do l=1,nleafin
       if (ileafproc(l) /= iproc) cycle
       do k=1,size(found(l)%i)
          j = found(l)%i(k)
          if (imarkproc(j) /= iproc) then
             nlist = nlist + 1
             imarkproc(j) = iproc
             imarkpos(j)  = nlist
             isendlist(nsendlist+nlist) = j
          endif
          found(l)%i(k) = imarkpos(j)
       enddo
    enddo
    nsendto(iproc) = nlist
    nsendlist = nsendlist + nlist
 enddo
 call set_offsets(nsendto,isendoff)

 !--pack the replies
 allocate(repsend(nprocs + nstatic*nsendlist + nleafin + ntotfound))
 ipos = 0
 do iproc=1,nprocs
    ipos0 = ipos
    repsend(ipos+1) = real(nsendto(iproc))
    ipos = ipos + 1
    do k=isendoff(iproc)+1,isendoff(iproc)+nsendto(iproc)
       j = isendlist(k)
       repsend(ipos+1:ipos+4) = xyzh(1:4,j)
       repsend(ipos+5) = get_pmass(j,.true.)
       repsend(ipos+6) = rho(j)
       repsend(ipos+7) = real(gradh(igradzeta,j))*real(gradh(igradomega,j))
       ipos = ipos + nstatic
    enddo
    do l=1,nleafin
       if (ileafproc(l) /= iproc) cycle
       nfound = size(found(l)%i)
       repsend(ipos+1) = real(nfound)
       repsend(ipos+2:ipos+1+nfound) = real(found(l)%i(1:nfound))
       ipos = ipos + 1 + nfound
    enddo
    nrep(iproc) = ipos - ipos0
 enddo
 if (allocated(sendvals)) deallocate(sendvals)
 allocate(sendvals(max(nvals*nsendlist,1)))

end subroutine answer_requests

!-----------------------------------------------------------------------
!+
!  for each requested leaf, the local particles that are neighbours of
!  at least one of its solver particles
!+
!-----------------------------------------------------------------------
subroutine find_leaf_neighbours(nleafin,ileafstart,reqrecv,xyzh,found)
 use neighkdtree, only:get_neighbour_list,listneigh
 use part,        only:isdead_or_accreted
 use dim,         only:maxpsph
 integer,       intent(in)    :: nleafin,ileafstart(:)
 real,          intent(in)    :: reqrecv(:),xyzh(:,:)
 type(intlist), intent(inout) :: found(:)
 real,    allocatable, save :: xyzcache(:,:)
 integer, allocatable, save :: ilist(:)
 !$omp threadprivate(xyzcache,ilist)
 integer :: l,ioff,np,nneigh,n,j,ip,ipi,nkeep

 !$omp parallel default(none) shared(nleafin,ileafstart,reqrecv,xyzh,found,maxpsph) &
 !$omp private(l,ioff,np,nneigh,n,j,ip,ipi,nkeep)
 if (.not.allocated(xyzcache)) allocate(xyzcache(4,maxcellcache))
 if (allocated(ilist)) then
    if (size(ilist) < size(xyzh,2)) deallocate(ilist)
 endif
 if (.not.allocated(ilist)) allocate(ilist(size(xyzh,2)))
 !$omp do schedule(runtime)
 do l=1,nleafin
    ioff = ileafstart(l)
    np = nint(reqrecv(ioff+1))
    call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true., &
                            cell_xpos=reqrecv(ioff+2:ioff+4),cell_xsizei=reqrecv(ioff+5), &
                            cell_rcuti=reqrecv(ioff+6))
    nkeep = 0
    over_neigh: do n=1,nneigh
       j = listneigh(n)
       if (j > maxpsph) cycle over_neigh
       if (isdead_or_accreted(xyzh(4,j))) cycle over_neigh
       do ip=1,np
          ipi = ioff + nleafhdr + 4*(ip-1)
          if (is_pair(reqrecv(ipi+1:ipi+4),xyzh(1:4,j))) then
             nkeep = nkeep + 1
             ilist(nkeep) = j
             cycle over_neigh
          endif
       enddo
    enddo over_neigh
    allocate(found(l)%i(nkeep))
    found(l)%i = ilist(1:nkeep)
 enddo
 !$omp enddo
 !$omp end parallel

end subroutine find_leaf_neighbours

!-----------------------------------------------------------------------
!+
!  read the replies: store the remote particles, and the remote candidates
!  of each exported leaf in compressed rows
!+
!-----------------------------------------------------------------------
subroutine read_replies(reprecv,irepinoff,iexport,nleafto,ileafoff)
 use io,          only:nprocs
 use neighkdtree, only:ncells
 real,    intent(in) :: reprecv(:)
 integer, intent(in) :: irepinoff(:),iexport(:,:),nleafto(:),ileafoff(:)
 integer, allocatable :: ifill(:)
 integer :: iproc,ipos,np,k,l,icell,nfound,ncand,nc

 !--count the remote particles and the candidates of each leaf
 nc = int(ncells)
 icand_start = 0
 do iproc=1,nprocs
    ipos = irepinoff(iproc)
    np = nint(reprecv(ipos+1))
    nrecvfrom(iproc) = np
    ipos = ipos + 1 + nstatic*np
    do l=ileafoff(iproc)+1,ileafoff(iproc)+nleafto(iproc)
       nfound = nint(reprecv(ipos+1))
       icell  = iexport(1,l)
       icand_start(icell+1) = icand_start(icell+1) + nfound
       ipos = ipos + 1 + nfound
    enddo
 enddo
 call set_offsets(nrecvfrom,irecvoff)
 nremote = sum(nrecvfrom)
 icand_start(1) = 1
 do icell=1,nc
    icand_start(icell+1) = icand_start(icell+1) + icand_start(icell)
 enddo
 ncand = icand_start(nc+1) - 1

 if (allocated(xyzh_remote)) deallocate(xyzh_remote,pmass_remote,rho_remote,zetaomega_remote,icand,recvvals)
 allocate(xyzh_remote(4,max(nremote,1)),pmass_remote(max(nremote,1)),rho_remote(max(nremote,1)), &
          zetaomega_remote(max(nremote,1)),icand(max(ncand,1)),recvvals(max(nvals*nremote,1)))

 !--store them
 allocate(ifill(nc))
 ifill = icand_start(1:nc) - 1
 do iproc=1,nprocs
    ipos = irepinoff(iproc) + 1
    do k=irecvoff(iproc)+1,irecvoff(iproc)+nrecvfrom(iproc)
       xyzh_remote(1:4,k)  = reprecv(ipos+1:ipos+4)
       pmass_remote(k)     = reprecv(ipos+5)
       rho_remote(k)       = reprecv(ipos+6)
       zetaomega_remote(k) = reprecv(ipos+7)
       ipos = ipos + nstatic
    enddo
    do l=ileafoff(iproc)+1,ileafoff(iproc)+nleafto(iproc)
       icell  = iexport(1,l)
       nfound = nint(reprecv(ipos+1))
       do k=1,nfound
          ifill(icell) = ifill(icell) + 1
          icand(ifill(icell)) = irecvoff(iproc) + nint(reprecv(ipos+1+k))
       enddo
       ipos = ipos + 1 + nfound
    enddo
 enddo

end subroutine read_replies

!-----------------------------------------------------------------------
!+
!  send nsend(iproc) reals of sendbuf to each task iproc-1 (in task
!  order), receiving nrecv(iproc) reals from each, starting at
!  recvbuf(irecvoff(iproc)+1)
!+
!-----------------------------------------------------------------------
subroutine alltoallv_reals(sendbuf,nsend,recvbuf,nrecv,irecvoff_out)
 use mpiutils, only:MPI_DEFAULT_REAL
 real,    intent(in)  :: sendbuf(:)
 integer, intent(in)  :: nsend(:)
 real, allocatable, intent(out) :: recvbuf(:)
 integer, intent(out) :: nrecv(:),irecvoff_out(:)
 integer :: isendoff_l(size(nsend)),ierr

 call MPI_ALLTOALL(nsend,1,MPI_INTEGER,nrecv,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 call set_offsets(nsend,isendoff_l)
 call set_offsets(nrecv,irecvoff_out)
 allocate(recvbuf(max(sum(nrecv),1)))
 call MPI_ALLTOALLV(sendbuf,nsend,isendoff_l,MPI_DEFAULT_REAL,recvbuf,nrecv,irecvoff_out, &
                    MPI_DEFAULT_REAL,MPI_COMM_WORLD,ierr)

end subroutine alltoallv_reals
#endif

end module mpiradiation
