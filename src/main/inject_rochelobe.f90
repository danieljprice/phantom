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
 real, private :: spd_inject = 0.
 real, private :: chi = 0.001
 real, private :: dNdt = 100.

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
 use physcon,   only:pi,twopi,solarm,years,gg
 use random,    only:ran2, rayleigh_deviate
 use units,     only:udist, umass, utime
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real :: m1,m2,q,radL1,h,u,theta_s,A,mu,theta_rand,r_rand,dNdt_code,Porb,r12
 real :: xyzL1(3),xyzi(3),vxyz(3),dr(3),x1(3),x2(3),x0(3),dxyz(3),vxyzL1(3),v1(3),v2(3)
 integer :: i_part,part_type,s1,wall_i,particles_to_place

 
!
!--find the L1 point
!
 if (nptmass < 2) call fatal('inject_rochelobe','not enough point masses for roche lobe injection')
 if (nptmass > 2) call fatal('inject_rochelobe','too many point masses for roche lobe injection')
 x1 = xyzmh_ptmass(1:3,1)
 x2 = xyzmh_ptmass(1:3,2)
 v1 = vxyz_ptmass(1:3,1)
 v2 = vxyz_ptmass(1:3,1)
 dr = x2 - x1
 r12 = sqrt(dot_product(dr,dr))
 m1 = xyzmh_ptmass(4,1)
 m2 = xyzmh_ptmass(4,2)
 q  = m2/m1
 mu = 1./(1 + q)
 radL1      = L1_point(m1/m2)                     ! find L1 point given binary mass ratio
 A  = mu / abs(radL1 - 1. + mu)**3 + (1. - mu)/abs(radL1 + mu)  ! See Lubow & Shu 1975
 theta_s = -acos( -4./(3.*A)+(1-8./(9.*A))**0.5)/2.             ! See Lubow & Shu 1975
 xyzL1(1:3) = xyzmh_ptmass(1:3,1) + radL1*dr(:)   ! set as vector position
 vxyzL1 = v1*sqrt( dot_product(xyzL1-x0, xyzL1-x0)/dot_product(x1-x0, x1-x0) ) ! orbital motion of L1 point
    !unclear if this is OK with eccentric orbits, but if you have Roche Lobe overflow, orbits should be
    !circularised anyway

 ! mass of gas particles is set by mass accretion rate and particle injection rate
 Mdotcode  = Mdot*(solarm/years)/(umass/utime)
 Porb      = twopi * sqrt( (r12*udist)**3 / (gg*(m1+m2)*umass) )
 dNdt_code = dNdt*utime / Porb
 massoftype(igas) = Mdotcode/dNdt_code

 ! get centre of mass of binary
 x0 = (m1*x1 + m2*x2)/(m1 + m2)

 ! how many particles do we need to place?
 if(npartoftype(igas)<8) then
    particles_to_place = 8-npartoftype(igas)      ! Seems to need at least eight gas particles to not crash
 else
    particles_to_place = max(0, int(0.5 + (time*Mdotcode/massoftype(igas)) - npartoftype(igas) ))
 endif
 do wall_i=1,particles_to_place

    ! calculate particle offset
    theta_rand = ran2(s1)*twopi
    r_rand = rayleigh_deviate(s1)*chi
    dxyz=(/0.0, cos(theta_rand), sin(theta_rand)/)*r_rand   ! Stream is placed randomly in a cylinder
                                                            ! with a Gaussian density distribution
 
    ! prepare to add a particle
    part_type = igas
    vxyz = (/ cos(theta_s), sin(theta_s), 0.0 /)*spd_inject
    h = hfact*chi
    u = 0.0   !get back to this
    i_part = npart + 1
    call rotate_into_plane(dxyz,vxyz,x2-xyzL1)
    vxyz = vxyz + vxyzL1
    xyzi = xyzL1 + dxyz

    !add the particle

    call add_or_update_particle(part_type, xyzi, vxyz, h, u, i_part, npart, npartoftype, xyzh, vxyzu)
    !print*, "Adding particle at ", xyzi
    !print*, "With velocity ", vxyz
    !print*, "L1 point calculated to be at ", xyzL1

 enddo

!print*, "White Dwarf is at ", x1
!print*, "Companion is at   ", x2
!print*, "Theta_s was ", theta_s*180/pi
!print*, "A was", A  
!print*, "Mu was", mu
!print*, "Mdot is ",Mdot
!print*, "Mdotcode is ", Mdotcode
!print*, "But Mdotcode should be ", Mdot*(solarm/years)/(umass/utime),"= ",Mdot*(solarm/years)/(umass/utime)*umass/utime," g/s"
!print*, "Because, now, umass= ",umass," and utime= ",utime

!print*, "Particlemass is ", massoftype(igas), " or ", massoftype(igas)*umass, " g"


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
 call write_inopt(dNdt,'dNdt','particle injection rate in particles/binary orbit',iunit)
 call write_inopt(spd_inject,'spd_inject','speed of injection at L1 in cm/s',iunit)
 call write_inopt(chi,'chi','width of injection stream in cm',iunit)

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
    print*,' DEBUG: Mdot is ',Mdot,' Msun/yr, which is ',Mdotcode,' in code units'
    
 case('dNdt')
    read(valstring,*,iostat=ierr) dNdt
    ngot = ngot + 1
    if (dNdt < 0.) call fatal(label,'dNdt < 0 in input options')
 case('spd_inject')
    read(valstring,*,iostat=ierr) spd_inject
    ngot = ngot + 1
 case('chi')
    read(valstring,*,iostat=ierr) chi
    ngot = ngot + 1
    if (dNdt <= 0.) call fatal(label,'chi <= 0 in input options')
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_inject

end module inject
