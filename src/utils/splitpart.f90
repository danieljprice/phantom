module splitpart
 implicit none
 integer, parameter :: ires = 1 ! use 12 particles per sphere

contains

subroutine split_particles(npart,npartoftype,xyzh,massoftype,ierr)
 use injectutils, only:get_parts_per_sphere
 use icosahedron, only:pixel2vector,compute_corners,compute_matrices
 use part,        only:copy_particle
 use io,          only:error
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: massoftype(:)
 integer, intent(out)   :: ierr
 integer :: npart_per_sphere,iold,inew,j
 real    :: dhfac,dx(3),sep
 real    :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)

 ierr = 0
 npart_per_sphere = get_parts_per_sphere(ires)
 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 !--check there is enough memory
 if (size(xyzh(1,:)) < npart*(npart_per_sphere+1)) then
    call error('split_particles','not enough memory, increase MAXP and recompile')
    ierr = 1
    return
 endif

 !--divide mass of all particles by 13
 massoftype(:) = massoftype(:) / (npart_per_sphere + 1)

 !--update npartoftype
 npartoftype(:) = npartoftype*(npart_per_sphere+1)

 !--amend smoothing length of original particles
 dhfac = 1./(npart_per_sphere + 1)**(1./3.)
 do iold=1,npart
    xyzh(4,iold) = xyzh(4,iold)*dhfac
 enddo

 !--find positions of the new particles
 inew = npart ! put new particles after original parts
 do iold=1,npart
    sep = xyzh(4,iold)
    ! add 12 new child particles
    do j=0,npart_per_sphere-1
       inew = inew + 1
       ! copy all properties from original particle
       call copy_particle(iold,inew)
       ! find positions of original particles
       call pixel2vector(j,ires,geodesic_R,geodesic_v,dx)
       xyzh(1:3,inew) = xyzh(1:3,iold) + sep*dx(:)
    enddo
    ! amend smoothing length of original particle
    npart = npart + npart_per_sphere
 enddo

end subroutine split_particles

end module splitpart
