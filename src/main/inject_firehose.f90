!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION:
!  Injection module for "firehose" simulations, used for numerical experiments
!  in Liptai et al. paper on tidal disruption events
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
 character(len=*), parameter, public :: inject_type = 'firehose'

 public :: inject_particles, write_options_inject, read_options_inject
 public :: init_inject

 real, private :: Mdot = 0.
 real, private :: Mdotcode = 0.
 real, private :: mach = 2.

contains

subroutine init_inject(ierr)
 integer, intent(out) :: ierr

 ierr = 0

end subroutine init_inject
!-----------------------------------------------------------------------
!+
!  Main routine handling injection at the L1 point.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass, &
           npart,npartoftype,dtinject)
  use io,        only:fatal
  use part,      only:igas,hfact
  use partinject,only:add_or_update_particle
  use physcon,   only:pi,solarr,au
  use units,     only:udist
  use eos,       only:gamma
  real,    intent(in)    :: time, dtlast
  real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
  integer, intent(inout) :: npart
  integer, intent(inout) :: npartoftype(:)
  real,    intent(out)   :: dtinject
  real :: Rp,Rtidal,Rstar,beta
  real :: xyzi(3),xyzL1(3),vxyz(3),cs
  real :: delta,h,u,vinject,time_between_walls,local_time,ymin,zmin,rcyl,rcyl2
  integer :: N,i,iy,iz,i_part,part_type,handled_walls
  integer :: outer_wall, inner_wall, inner_handled_wall, particles_per_wall

  ! get the location of the injection point
  ! the following are just to base this on parameters used in our TDE simulation
  ! but Rp in the end is just a location along the x axis
  beta = 1.    ! penetration factor
  Rstar = solarr/udist
  Rtidal = (solarr/udist)*(1e6)**(1./3.)
  Rp = Rtidal/beta
  !print*,' Rtidal = ',Rtidal,' = ',Rtidal*udist/au,' au'
  !print*,' Rp = ',Rp
  !print*,' solarr = ',Rstar
  !print*,'  1 au = ',au/udist
  !print*,' Rschwarzschild = ',2.*udist/au,' au, or ',2*udist/solarr,' Rsun'

  ! specify the injection point
  xyzL1(1:3) = (/Rp,0.,0./)

  N = 32              ! number of particles in x axis of cylinder
  rcyl = Rstar !0.01*Rp  ! radius of injection cylinder, as fraction of L1 distance
  rcyl2 = rcyl*rcyl
  delta = 2.*rcyl/N  ! particle separation in cylinder
  ymin = -rcyl       ! to ensure flow is centred around injection point
  zmin = -rcyl
  h = hfact*delta
  handled_walls = 4

  ! give the injection velocity in terms of the velocity of a parabolic orbit at Rp
  ! then use this to specify the sound speed from the desired Mach number
  vinject = sqrt(2./Rp)
  cs = vinject/mach
  u = cs**2/(gamma*(gamma-1))

  print*,'injecting at R=',Rp,' with  v = ',sqrt(2./Rp),' Mach # = ',vinject/cs
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
      xyzi(2) = local_time * vinject
      print *, '==== ', i, xyzi(1)
      do iy = 1,N
        do iz = 1,N
           xyzi(1) = ymin + (iy-.5)*delta
           xyzi(3) = zmin + (iz-.5)*delta
           ! crop to cylinder
           if (xyzi(1)**2 + xyzi(3)**2 < rcyl2) then
              ! rotate direction so that x axis is pointing towards secondary
              vxyz = (/ 0., vinject, 0. /)
              ! add position offset
              xyzi = xyzi + xyzL1
              i_part = i_part + 1
              ! Another brick in the wall
              call add_or_update_particle(part_type, xyzi, vxyz, h, u, i_part, npart, npartoftype, xyzh, vxyzu)
           endif
        enddo
      enddo
      !read*
  enddo
  dtinject = huge(dtinject)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(Mdot,'Mdot','mass injection rate at L1, in Msun/yr',iunit)
 call write_inopt(mach,'mach','Mach number of injected stream',iunit)

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
    print*,' Mdot is ',Mdot,' Msun/yr, which is ',Mdotcode,' in code units'
 case('mach')
   read(valstring,*,iostat=ierr) mach
   ngot = ngot + 1
   if (mach <=  0.) call fatal(label,'mach number <= 0 in input options')

 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_inject

end module inject
