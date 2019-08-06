!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION:
!  Handles uniform distribution injection
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    wind_resolution  -- resolution of the wind -- DO NOT CHANGE AFTER RUNNING SETUP --
!    wind_temperature -- temperature of the wind (Kelvin)
!
!  DEPENDENCIES: boundary, eos, infile_utils, io, part, partinject,
!    physcon, units
!+
!--------------------------------------------------------------------------
module inject
 implicit none
 character(len=*), parameter, public :: inject_type = 'unifwind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 real, public :: wind_density = 7.2d-16
 real, public :: wind_velocity = 29.
 real, public :: wind_resolution = 64.
 real, public :: wind_temperature = 1700.
 private

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use part,      only:hfact,igas,iboundary,massoftype
 use partinject,only:add_or_update_particle
 use io,        only:iprint
 use units,     only:umass,udist,utime
 use physcon,   only:Rg
 use eos,       only:gamma
 use boundary,  only:xmin,xmax,ymin,ymax,zmin,zmax
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 integer, parameter :: handled_walls = 3
 real, parameter :: mu = 1.26 ! Used in Bowen (1988)
 real :: rho, v, energy_to_temperature_ratio, u, h, delta, time_between_walls
 integer :: N, outer_wall, inner_wall, inner_handled_wall, particles_per_wall
 integer :: i, iy, iz, i_part, injected_particle, part_type
 real :: local_time, vxyz(3), pxyz(3)

 rho = wind_density / (umass/udist**3)
 v = wind_velocity * 1.d5 / (udist/utime)
 N = int(wind_resolution)
 energy_to_temperature_ratio = Rg/(mu*(gamma-1.))/(udist/utime)**2
 u = wind_temperature * energy_to_temperature_ratio
 delta = (ymax-ymin)/N
 time_between_walls = delta/v
 h = hfact * delta / 2.

 outer_wall = ceiling((time-dtlast)/time_between_walls)
 inner_wall = ceiling(time/time_between_walls)-1
 inner_handled_wall = inner_wall+handled_walls
 particles_per_wall = N**2

 !print *, "t = ", time
 !print *, "dt last = ", dtlast
 !print *, "delta t = ", time_between_walls
 !print *, "Injecting wall ", inner_wall, " to ", outer_wall
 !print *, "Handling wall ", inner_handled_wall, " to ", inner_wall-1
 print *, ' v = ', v
 print *, '*** ', time, dtlast, time_between_walls, inner_wall, outer_wall

 vxyz = (/ v, 0., 0. /)
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
    pxyz(1) = local_time * v
    print *, '==== ', i, pxyz(1)
    do iy = 1,N
       do iz = 1,N
          pxyz(2) = ymin + (iy-.5)*delta
          pxyz(3) = zmin + (iz-.5)*delta
          i_part = i_part + 1
          call add_or_update_particle(part_type, pxyz, vxyz, h, u, i_part, npart, npartoftype, xyzh, vxyzu) ! Another brick in the wall
       enddo
    enddo
 enddo
 !
 !-- timestep constraint
 !
 dtinject = time_between_walls

end subroutine

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(wind_velocity,'wind_velocity', &
      'velocity at which wind is injected (km/s) -- DO NOT CHANGE AFTER RUNNING SETUP --',iunit)
 call write_inopt(wind_density,'wind_density', &
      'wind density (g/cm³) -- DO NOT CHANGE AFTER RUNNING SETUP --',iunit)
 call write_inopt(wind_temperature,'wind_temperature','temperature of the wind (Kelvin)',iunit)
 call write_inopt(wind_resolution,'wind_resolution','resolution of the wind -- DO NOT CHANGE AFTER RUNNING SETUP --',iunit)
end subroutine

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only: fatal, error, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('wind_velocity')
    read(valstring,*,iostat=ierr) wind_velocity
    ngot = ngot + 1
    if (wind_velocity < 0.)    call fatal(label,'invalid setting for wind_velocity (<0)')
    if (wind_velocity > 1.e10) call error(label,'wind_velocity is huge!!!')
 case('wind_density')
    read(valstring,*,iostat=ierr) wind_density
    ngot = ngot + 1
    if (wind_density < 0.)    call fatal(label,'invalid setting for wind_density (<0)')
 case('wind_temperature')
    read(valstring,*,iostat=ierr) wind_temperature
    ngot = ngot + 1
    if (wind_temperature < 0.)    call fatal(label,'invalid setting for wind_temperature (<0)')
 case('wind_resolution')
    read(valstring,*,iostat=ierr) wind_resolution
    ngot = ngot + 1
    if (int(wind_resolution) < 1) call fatal(label,'wind_resolution must be bigger than zero')
 end select

 igotall = (ngot >= 4)
end subroutine

end module inject
