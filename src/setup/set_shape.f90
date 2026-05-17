!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module shape
!
! Populate a lattice then crop it to a user-specified shape.
! The shape is read from a simple text file:
!   sphere radius
!   ellipsoid ax by cz
!   box hx hy hz
!   cylinder radius height [x|y|z]
!   mesh path/to/shape.obj [scale]
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon, random, stretchmap, unifdis
   use unifdis, only:set_unifdis
 implicit none

 public :: set_shape

 private

contains

subroutine set_shape(lattice,id,master,np_requested,x0,rmax,hfact,np,xyzh,nptot,objfile)
 character(len=*), intent(in)    :: lattice
 integer,          intent(in)    :: id,master
 integer,          intent(in)    :: np_requested
 real,             intent(in)    :: x0(3),rmax,hfact
 integer,          intent(out)   :: np
 real,             intent(out)   :: xyzh(:,:)
 integer(kind=8),  intent(inout) :: nptot
 character(len=*), intent(in)    :: objfile
 integer :: i,iter,np_try,np_keep,ncube
 real    :: xmin,xmax,ymin,ymax,zmin,zmax,delta,ratio
 character(len=32) :: shape_kind
 character(len=1)  :: axis
 real :: p1,p2,p3
 character(len=256) :: meshfile
 real, allocatable :: vertices(:,:)
 integer, allocatable :: faces(:,:)
 real :: bmin(3),bmax(3)
 logical :: mesh_ok

 xmin = -rmax; xmax = rmax
 ymin = -rmax; ymax = rmax
 zmin = -rmax; zmax = rmax

 call read_shape_file(objfile,rmax,shape_kind,p1,p2,p3,axis,meshfile,id,master)
 mesh_ok = .false.
 if (trim(shape_kind) == 'mesh') then
    call load_obj_mesh(trim(meshfile),p1,vertices,faces,bmin,bmax,mesh_ok,id,master)
    if (.not.mesh_ok) then
      shape_kind = 'sphere'
      p1 = rmax
      p2 = rmax
      p3 = rmax
      axis = 'z'
    endif
 endif

 np = 0
 nptot = 0_8
 ncube = max(6,nint(np_requested**(1.0/3.0)))
 delta = 2.0*rmax/real(ncube)

 do iter=1,8
    np = 0
    nptot = 0_8
    call set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,delta,hfact,np,xyzh,.false., &
                     nptot=nptot,verbose=(id==master),centre=.true.)

    np_try = np
    if (trim(shape_kind) == 'mesh' .and. mesh_ok) then
       np_keep = crop_particles_mesh(np,xyzh,vertices,faces,bmin,bmax)
    else
       np_keep = crop_particles(np,xyzh,shape_kind,p1,p2,p3,axis)
    endif
    np = np_keep
    nptot = int(np_keep,kind=8)

    if (np_try <= 0) exit
    if (np_requested <= 0) exit
    if (abs(np_keep-np_requested) <= max(8,np_requested/50)) exit

    ratio = real(np_keep)/real(np_requested)
    if (ratio <= 0.) then
       delta = 0.9*delta
    else
       delta = delta*ratio**(1.0/3.0)
    endif
 enddo

 if (id==master) then
    write(*,"(1x,a,1x,a)") 'shape type:',trim(shape_kind)
    write(*,"(1x,a,1x,i9,a,i9)") 'particles kept:',np,' (target ',np_requested,')'
    write(*,"(1x,a,3(es10.3,1x))") 'shifting origin to ',x0(:)
 endif

 do i=1,np
    xyzh(1,i) = xyzh(1,i) + x0(1)
    xyzh(2,i) = xyzh(2,i) + x0(2)
    xyzh(3,i) = xyzh(3,i) + x0(3)
 enddo

 if (allocated(vertices)) deallocate(vertices)
 if (allocated(faces)) deallocate(faces)

end subroutine set_shape

integer function crop_particles(np,xyzh,shape_kind,p1,p2,p3,axis) result(np_keep)
 integer,          intent(in)    :: np
 real,             intent(inout) :: xyzh(:,:)
 character(len=*), intent(in)    :: shape_kind
 character(len=1), intent(in)    :: axis
 real,             intent(in)    :: p1,p2,p3
 integer :: i,j

 j = 0
 do i=1,np
    if (inside_shape(xyzh(1,i),xyzh(2,i),xyzh(3,i),shape_kind,p1,p2,p3,axis)) then
       j = j + 1
       if (j /= i) xyzh(:,j) = xyzh(:,i)
    endif
 enddo
 np_keep = j
end function crop_particles

integer function crop_particles_mesh(np,xyzh,vertices,faces,bmin,bmax) result(np_keep)
 integer, intent(in) :: np
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in) :: vertices(:,:)
 integer, intent(in) :: faces(:,:)
 real,    intent(in) :: bmin(3),bmax(3)
 integer :: i,j
 real :: x,y,z

 j = 0
 do i=1,np
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    if (x < bmin(1) .or. x > bmax(1)) cycle
    if (y < bmin(2) .or. y > bmax(2)) cycle
    if (z < bmin(3) .or. z > bmax(3)) cycle
    if (point_inside_mesh((/x,y,z/),vertices,faces)) then
       j = j + 1
       if (j /= i) xyzh(:,j) = xyzh(:,i)
    endif
 enddo
 np_keep = j
end function crop_particles_mesh

logical function inside_shape(x,y,z,shape_kind,p1,p2,p3,axis) result(ok)
 real,             intent(in) :: x,y,z,p1,p2,p3
 character(len=*), intent(in) :: shape_kind
 character(len=1), intent(in) :: axis
 real :: r2,rr

 ok = .false.
 select case(trim(shape_kind))
 case('sphere')
    ok = (x*x + y*y + z*z <= p1*p1)
 case('ellipsoid')
    if (p1 > 0. .and. p2 > 0. .and. p3 > 0.) then
       rr = (x/p1)**2 + (y/p2)**2 + (z/p3)**2
       ok = (rr <= 1.0)
    endif
 case('box')
    ok = (abs(x) <= p1 .and. abs(y) <= p2 .and. abs(z) <= p3)
 case('cylinder')
    select case(axis)
    case('x')
       r2 = y*y + z*z
       ok = (r2 <= p1*p1 .and. abs(x) <= 0.5*p2)
    case('y')
       r2 = x*x + z*z
       ok = (r2 <= p1*p1 .and. abs(y) <= 0.5*p2)
    case default
       r2 = x*x + y*y
       ok = (r2 <= p1*p1 .and. abs(z) <= 0.5*p2)
    end select
 end select
end function inside_shape

subroutine read_shape_file(objfile,rmax,shape_kind,p1,p2,p3,axis,meshfile,id,master)
 character(len=*), intent(in)  :: objfile
 real,             intent(in)  :: rmax
 character(len=32),intent(out) :: shape_kind
 real,             intent(out) :: p1,p2,p3
 character(len=1), intent(out) :: axis
 character(len=256),intent(out) :: meshfile
 integer,          intent(in)  :: id,master
 integer :: iu,ios
 character(len=256) :: line
 character(len=32)  :: key
 logical :: loaded

 shape_kind = 'sphere'
 p1 = rmax
 p2 = rmax
 p3 = rmax
 axis = 'z'
 meshfile = ''
 loaded = .false.

 if (is_obj_path(trim(objfile))) then
    shape_kind = 'mesh'
    meshfile = trim(objfile)
    p1 = 1.0
    p2 = 0.0
    p3 = 0.0
    return
 endif

 open(newunit=iu,file=trim(objfile),status='old',action='read',iostat=ios)
 if (ios /= 0) then
    if (id==master) write(*,"(1x,a,1x,a)") 'shape file not found, defaulting to sphere:',trim(objfile)
    return
 endif

 do
    read(iu,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#') cycle

    key = ''
    read(line,*,iostat=ios) key
    if (ios /= 0) cycle

    select case(trim(key))
    case('sphere')
       read(line,*,iostat=ios) key,p1
       if (ios == 0) then
          shape_kind = 'sphere'
          loaded = .true.
       endif
    case('ellipsoid')
       read(line,*,iostat=ios) key,p1,p2,p3
       if (ios == 0) then
          shape_kind = 'ellipsoid'
          loaded = .true.
       endif
    case('box')
       read(line,*,iostat=ios) key,p1,p2,p3
       if (ios == 0) then
          shape_kind = 'box'
          loaded = .true.
       endif
    case('cylinder')
       axis = 'z'
       read(line,*,iostat=ios) key,p1,p2,axis
       if (ios /= 0) then
          read(line,*,iostat=ios) key,p1,p2
       endif
       if (ios == 0) then
          shape_kind = 'cylinder'
          loaded = .true.
       endif
    case('mesh')
       p1 = 1.0
       p2 = 0.0
       p3 = 0.0
       read(line,*,iostat=ios) key,meshfile,p1
       if (ios /= 0) then
          read(line,*,iostat=ios) key,meshfile
       endif
       if (ios == 0) then
          shape_kind = 'mesh'
          loaded = .true.
       endif
    end select

    if (loaded) exit
 enddo

 close(iu)

 if (.not.loaded) then
    shape_kind = 'sphere'
    p1 = rmax
    p2 = rmax
    p3 = rmax
    axis = 'z'
    meshfile = ''
    if (id==master) write(*,"(1x,a,1x,a)") 'could not parse shape file, using sphere:',trim(objfile)
 endif

end subroutine read_shape_file

subroutine load_obj_mesh(filename,scale,vertices,faces,bmin,bmax,ok,id,master)
 character(len=*), intent(in) :: filename
 real,             intent(in) :: scale
 real, allocatable, intent(out) :: vertices(:,:)
 integer, allocatable, intent(out) :: faces(:,:)
 real, intent(out) :: bmin(3),bmax(3)
 logical, intent(out) :: ok
 integer, intent(in) :: id,master
 integer :: iu,ios,nv,nftri,iv,ifc
 character(len=512) :: line
 character(len=1) :: tag
 integer :: i1,i2,i3,i4
 real :: x,y,z

 ok = .false.
 nv = 0
 nftri = 0

 open(newunit=iu,file=trim(filename),status='old',action='read',iostat=ios)
 if (ios /= 0) then
    if (id==master) write(*,"(1x,a,1x,a)") 'mesh file not found:',trim(filename)
    return
 endif

 do
    read(iu,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#') cycle
    read(line,*,iostat=ios) tag
    if (ios /= 0) cycle
    if (trim(tag) == 'v') then
       nv = nv + 1
    elseif (trim(tag) == 'f') then
       i4 = 0
       call parse_face4(line,i1,i2,i3,i4,ios)
       if (ios == 0) then
          nftri = nftri + 1
          if (i4 /= 0) nftri = nftri + 1
       endif
    endif
 enddo

 if (nv < 4 .or. nftri < 4) then
    close(iu)
    if (id==master) write(*,"(1x,a)") 'mesh file has too few vertices/faces'
    return
 endif

 allocate(vertices(3,nv))
 allocate(faces(3,nftri))
 rewind(iu)

 iv = 0
 ifc = 0
 do
    read(iu,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#') cycle
    read(line,*,iostat=ios) tag
    if (ios /= 0) cycle
    if (trim(tag) == 'v') then
       read(line,*,iostat=ios) tag,x,y,z
       if (ios /= 0) cycle
       iv = iv + 1
       vertices(1,iv) = x
       vertices(2,iv) = y
       vertices(3,iv) = z
    elseif (trim(tag) == 'f') then
       call parse_face4(line,i1,i2,i3,i4,ios)
       if (ios /= 0) cycle
       call normalize_face_indices(i1,i2,i3,i4,nv)
       if (.not.valid_face(i1,i2,i3,nv)) cycle
       ifc = ifc + 1
       faces(:,ifc) = (/i1,i2,i3/)
       if (i4 /= 0) then
          if (.not.valid_face(i1,i3,i4,nv)) cycle
          ifc = ifc + 1
          faces(:,ifc) = (/i1,i3,i4/)
       endif
    endif
 enddo
 close(iu)

 if (ifc < 4 .or. iv < 4) then
    if (allocated(vertices)) deallocate(vertices)
    if (allocated(faces)) deallocate(faces)
    if (id==master) write(*,"(1x,a)") 'mesh parse failed after reading file'
    return
 endif

 bmin = minval(vertices,dim=2)
 bmax = maxval(vertices,dim=2)
 call center_and_scale_vertices(vertices,bmin,bmax,scale)
 bmin = minval(vertices,dim=2)
 bmax = maxval(vertices,dim=2)
 ok = .true.

 if (id==master) then
    write(*,"(1x,a,1x,a)") 'loaded mesh:',trim(filename)
    write(*,"(1x,a,i8,a,i8)") 'mesh vertices=',iv,' triangles=',ifc
 endif
end subroutine load_obj_mesh

subroutine parse_face4(line,i1,i2,i3,i4,ios)
 character(len=*), intent(in) :: line
 integer, intent(out) :: i1,i2,i3,i4,ios
 character(len=64) :: tag,t1,t2,t3,t4

 i1 = 0; i2 = 0; i3 = 0; i4 = 0
 t1 = ''; t2 = ''; t3 = ''; t4 = ''
 read(line,*,iostat=ios) tag,t1,t2,t3,t4
 if (ios /= 0) then
    ios = 0
    read(line,*,iostat=ios) tag,t1,t2,t3
    if (ios /= 0) return
 endif
 i1 = parse_face_token(t1)
 i2 = parse_face_token(t2)
 i3 = parse_face_token(t3)
 if (len_trim(t4) > 0) i4 = parse_face_token(t4)
 if (i1 == 0 .or. i2 == 0 .or. i3 == 0) ios = 1
end subroutine parse_face4

integer function parse_face_token(tok) result(idx)
 character(len=*), intent(in) :: tok
 integer :: p,ios
 character(len=64) :: first
 idx = 0
 p = index(tok,'/')
 if (p > 1) then
    first = tok(:p-1)
 else
    first = tok
 endif
 read(first,*,iostat=ios) idx
 if (ios /= 0) idx = 0
end function parse_face_token

subroutine normalize_face_indices(i1,i2,i3,i4,nv)
 integer, intent(inout) :: i1,i2,i3,i4
 integer, intent(in) :: nv
 if (i1 < 0) i1 = nv + i1 + 1
 if (i2 < 0) i2 = nv + i2 + 1
 if (i3 < 0) i3 = nv + i3 + 1
 if (i4 < 0) i4 = nv + i4 + 1
end subroutine normalize_face_indices

logical function valid_face(i1,i2,i3,nv) result(ok)
 integer, intent(in) :: i1,i2,i3,nv
 ok = (i1 >= 1 .and. i1 <= nv .and. i2 >= 1 .and. i2 <= nv .and. i3 >= 1 .and. i3 <= nv)
end function valid_face

subroutine center_and_scale_vertices(vertices,bmin,bmax,scale)
 real, intent(inout) :: vertices(:,:)
 real, intent(in) :: bmin(3),bmax(3),scale
 real :: c(3),extent,maxextent,s
 integer :: i

 c = 0.5*(bmin + bmax)
 extent = maxval(bmax - bmin)
 if (extent <= 0.) then
    s = scale
 else
    s = scale * 2.0 / extent
 endif
 do i=1,size(vertices,2)
    vertices(:,i) = (vertices(:,i) - c)*s
 enddo
end subroutine center_and_scale_vertices

logical function point_inside_mesh(p,vertices,faces) result(is_inside)
 real, intent(in) :: p(3)
 real, intent(in) :: vertices(:,:)
 integer, intent(in) :: faces(:,:)
 integer :: i,nhit
 real, parameter :: dir(3) = (/1.0,0.39271,0.18391/)
 real, parameter :: eps = 1.0e-7

 nhit = 0
 do i=1,size(faces,2)
    if (ray_hits_triangle(p,dir,vertices(:,faces(1,i)),vertices(:,faces(2,i)),vertices(:,faces(3,i)),eps)) nhit = nhit + 1
 enddo
 is_inside = mod(nhit,2) == 1
end function point_inside_mesh

logical function ray_hits_triangle(orig,dir,v0,v1,v2,eps) result(hit)
 real, intent(in) :: orig(3),dir(3),v0(3),v1(3),v2(3),eps
 real :: edge1(3),edge2(3),h(3),s(3),q(3)
 real :: a,f,u,v,t

 edge1 = v1 - v0
 edge2 = v2 - v0
 h = cross(dir,edge2)
 a = dot_product(edge1,h)
 if (abs(a) < eps) then
    hit = .false.
    return
 endif
 f = 1.0/a
 s = orig - v0
 u = f*dot_product(s,h)
 if (u < 0.0 .or. u > 1.0) then
    hit = .false.
    return
 endif
 q = cross(s,edge1)
 v = f*dot_product(dir,q)
 if (v < 0.0 .or. u + v > 1.0) then
    hit = .false.
    return
 endif
 t = f*dot_product(edge2,q)
 hit = (t > eps)
end function ray_hits_triangle

pure function cross(a,b) result(c)
 real, intent(in) :: a(3),b(3)
 real :: c(3)
 c(1) = a(2)*b(3) - a(3)*b(2)
 c(2) = a(3)*b(1) - a(1)*b(3)
 c(3) = a(1)*b(2) - a(2)*b(1)
end function cross

logical function is_obj_path(path) result(ok)
 character(len=*), intent(in) :: path
 integer :: n
 ok = .false.
 n = len_trim(path)
 if (n < 4) return
 if (path(n-3:n) == '.obj' .or. path(n-3:n) == '.OBJ') ok = .true.
end function is_obj_path

end module shape
