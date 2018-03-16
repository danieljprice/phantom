!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION:
!  Handles Roche Lobe injection
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    Mdot -- mass injection rate at L1, in Msun/yr
!
!  DEPENDENCIES: infile_utils, io, part, partinject, physcon, setbinary,
!    units
!+
!--------------------------------------------------------------------------
module inject
 implicit none
 character(len=*), parameter, public :: inject_type = 'rochelobe'

 public :: inject_particles, write_options_inject, read_options_inject

 real, private :: Mdot = 0.
 real, private :: Mdotcode = 0.

contains

!-----------------------------------------------------------------------
!+
!  Main routine handling injection at the L1 point.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass, &
           npart,npartoftype)
 use io,        only:fatal,iverbose
 use part,      only:nptmass,massoftype,igas,hfact
 use partinject,only:add_or_update_particle
 use setbinary, only:L1_point
 use physcon,   only:pi
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real :: m1,m2,q,radL1
 real :: xyzi(3),xyzL1(3),vxyz(3),dr(3),dir(3),x1(3),x2(3),x0(3),v1(3),v2(3),dv(3)
 real :: delta,h,u,Minject,vinject,time_between_walls,local_time,ymin,zmin,rcyl,rcyl2
 integer :: N,i,iy,iz,i_part,part_type,handled_walls
 integer :: outer_wall, inner_wall, inner_handled_wall, particles_per_wall
!
!--find the L1 point
!
 if (nptmass < 2) call fatal('inject_rochelobe','not enough point masses for roche lobe injection')
 if (nptmass > 2) call fatal('inject_rochelobe','too many point masses for roche lobe injection')
 x1 = xyzmh_ptmass(1:3,1)
 x2 = xyzmh_ptmass(1:3,2)
 dr = x2 - x1

 m1 = xyzmh_ptmass(4,1)
 m2 = xyzmh_ptmass(4,2)
 q  = m2/m1
 radL1      = L1_point(m1/m2)                     ! find L1 point given binary mass ratio
 xyzL1(1:3) = xyzmh_ptmass(1:3,1) + radL1*dr(:)   ! set as vector position

 ! get centre of mass of binary
 x0 = (m1*x1 + m2*x2)/(m1 + m2)

 ! get direction vector in which to point the stream
 dir = x2 - xyzL1

 ! get angular velocity of binary
 !v1 = vxyz_ptmass(1:3,1)
 !v2 = vxyz_ptmass(1:3,2)
 !dv = v2 - v1
 !call get_v_spherical(dir,dv,vr,vphi,vtheta)
 !r = sqrt(dot_product(r,r))
 !omega0 = vphi/r

 N = 4              ! number of particles in x axis of cylinder
 rcyl = 0.01*radL1  ! radius of injection cylinder, as fraction of L1 distance
 rcyl2 = rcyl*rcyl
 delta = 2.*rcyl/N  ! particle separation in cylinder
 ymin = -rcyl   ! to ensure flow is centred around L1 point
 zmin = -rcyl
 h = hfact*delta

 handled_walls = 4
 vinject = 0.05

 print*,'injecting at ',xyzi(1:3)
!
!--inject material at the L1 point, with the orbital motion of the secondary
!
!  deltat = time - dtlast
!  Minject = deltat*Mdotcode
!  ninject = int(Minject/massoftype(igas)) - 1
!  deltatp = deltat/real(ninject - 1)
!  vnew(:) = vxyz_ptmass(1:3,2)
!  vinject = sqrt(dot_product(vnew,vnew))

 time_between_walls = delta/vinject
 outer_wall = ceiling((time-dtlast)/time_between_walls)
 inner_wall = ceiling(time/time_between_walls)-1
 inner_handled_wall = inner_wall+handled_walls
 particles_per_wall = int(0.25*pi*N**2)  ! cross section of cylinder

 print *, "t = ", time
 print *, "dt last = ", dtlast
 print *, "delta t = ", time_between_walls
 print *, "Injecting wall ", inner_wall, " to ", outer_wall
 print *, "Handling wall ", inner_handled_wall, " to ", inner_wall-1
 print *, ' v = ', vinject
 print *, '*** ', time, dtlast, time_between_walls, inner_wall, outer_wall

 do i=inner_handled_wall,outer_wall,-1
    local_time = time - i*time_between_walls
    if (i  >  inner_wall) then
       ! Handled wall
       i_part = (inner_handled_wall-i)*particles_per_wall
       part_type = igas
    else
       ! Outer wall
       i_part = npart
       part_type = igas
    endif
    xyzi(1) = local_time * vinject
    print *, '==== ', i, xyzi(1)
    do iy = 1,N
       do iz = 1,N
          xyzi(2) = ymin + (iy-.5)*delta
          xyzi(3) = zmin + (iz-.5)*delta
          ! crop to cylinder
          if (xyzi(2)**2 + xyzi(3)**2 < rcyl2) then
             ! rotate direction so that x axis is pointing towards secondary
             vxyz = (/ vinject, 0., 0. /)
             call rotate_into_plane(xyzi,vxyz,dir)
             ! add position offset
             xyzi = xyzi + xyzL1
             ! add rotation velocity of primary
             !vxyz = vxyz + vxyz_ptmass(1:3,1)
             i_part = i_part + 1
             ! Another brick in the wall
             call add_or_update_particle(part_type, xyzi, vxyz, h, u, i_part, npart, npartoftype, xyzh, vxyzu)
          endif
       enddo
    enddo
    !read*
 enddo

end subroutine inject_particles

subroutine get_v_spherical(r1,v1,vr,vphi,vtheta)
 real, intent(in)  :: r1(3),v1(3)
 real, intent(out) :: vr,vphi,vtheta
 real :: rcyl,r,rcyl2

 rcyl2  = r1(1)**2 + r1(2)**2
 rcyl   = sqrt(rcyl2)
 r      = sqrt(rcyl2 + r1(3)**2)
 vr     = dot_product(v1,r1/r)
 vphi   = v1(1)*(-r1(2)/rcyl) + v1(2)*(r1(1)/rcyl)
 vtheta = dot_product(v1,(/r1(1)*r1(3),r1(2)*r1(3),-rcyl*rcyl/)/(r*rcyl))

end subroutine get_v_spherical

!
! Routine to rotate one vector to point in the direction of another
!
subroutine rotate_into_plane(r1,v1,ref)
 real, intent(inout) :: r1(3),v1(3)
 real, intent(in)    :: ref(3)
 real :: dr,dphi,dtheta,r,phi,theta
 real :: cosphi,sinphi,costheta,sintheta
 real :: vr,vphi,vtheta

 dr     = sqrt(dot_product(ref,ref))
 dphi   = atan2(ref(2),ref(1))
 dtheta = 0. !acos(ref(3)/dr)

 r     = sqrt(dot_product(r1,r1))
 phi   = atan2(r1(2),r1(1))
 theta = acos(r1(3)/r)

 call get_v_spherical(r1,v1,vr,vphi,vtheta)

 phi   = phi + dphi
 theta = theta + dtheta

 !print*,'dtheta = ',acos(ref(3)/dr)
 cosphi = cos(phi)
 sinphi = sin(phi)
 sintheta = sin(theta)
 costheta = cos(theta)

 r1(1) = r*cosphi*sintheta
 r1(2) = r*sinphi*sintheta
 r1(3) = r*costheta

 v1(1) = vr*cosphi*sintheta - vphi*sinphi + vtheta*cosphi*costheta
 v1(2) = vr*sinphi*sintheta + vphi*cosphi + vtheta*sinphi*costheta
 v1(3) = vr*costheta - vtheta*sintheta

end subroutine rotate_into_plane

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(Mdot,'Mdot','mass injection rate at L1, in Msun/yr',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal,error
 use physcon, only:solarm,years
 use units,   only:umass,utime
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 select case(trim(name))
 case('Mdot')
    read(valstring,*,iostat=ierr) Mdot
    ngot = ngot + 1
    if (Mdot  <  0.) call fatal(label,'Mdot < 0 in input options')
!
!--convert mass injection rate to code units
!
    Mdotcode = Mdot*(umass/solarm)/(utime/years)
    print*,' DEBUG: Mdot is ',Mdot,' Msun/yr, which is ',Mdotcode,' in code units'
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_inject

end module inject
