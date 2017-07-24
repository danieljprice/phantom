module extern_gr
 implicit none

 public :: get_gr_force, update_grforce_leapfrog

 private

contains

!---------------------------------------------------------------
!+
!  subroutine to compute the force due to spacetime curvature
!+
!---------------------------------------------------------------
subroutine get_gr_force(xyzi,veli,densi,ui,pi,fexti)
 use force_gr, only:get_forcegr
 real, intent(in)  :: xyzi(3),veli(3),densi,ui,pi
 real, intent(out) :: fexti(3)
 logical :: its_a_testparticle
 its_a_testparticle = .false.

 if (its_a_testparticle) then
    call get_sourceterms(xyzi,veli,fexti)
 else
    call get_forcegr(xyzi,veli,densi,ui,pi,fexti)
 endif

end subroutine get_gr_force

subroutine update_grforce_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dt,xi,yi,zi,densi,ui,pi)
 use io,             only:fatal
 real, intent(in)    :: dt,xi,yi,zi
 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(inout) :: fexti(3)
 real, intent(in)    :: densi,ui,pi
 real                :: fextv(3)
 real                :: v1x, v1y, v1z, v1xold, v1yold, v1zold, vhalf2, erri, dton2
 logical             :: converged
 integer             :: its, itsmax
 integer, parameter  :: maxitsext = 50 ! maximum number of iterations on external force
 real, parameter :: tolv = 1.e-2
 real, parameter :: tolv2 = tolv*tolv
 real,dimension(3) :: pos,vel

 itsmax = maxitsext
 its = 0
 converged = .false.
 dton2 = 0.5*dt

 v1x = vhalfx
 v1y = vhalfy
 v1z = vhalfz
 vhalf2 = vhalfx*vhalfx + vhalfy*vhalfy + vhalfz*vhalfz
 fextv = 0. ! to avoid compiler warning

 iterations : do while (its < itsmax .and. .not.converged)
    its = its + 1
    erri = 0.
    v1xold = v1x
    v1yold = v1y
    v1zold = v1z
    pos = (/xi,yi,zi/)
    vel = (/v1x,v1y,v1z/)
    call get_gr_force(pos,vel,densi,ui,pi,fextv)
!    xi = pos(1)
!    yi = pos(2)
!    zi = pos(3)
    v1x = vel(1)
    v1y = vel(2)
    v1z = vel(3)

    v1x = vhalfx + dton2*(fxi + fextv(1))
    v1y = vhalfy + dton2*(fyi + fextv(2))
    v1z = vhalfz + dton2*(fzi + fextv(3))

    erri = (v1x - v1xold)**2 + (v1y - v1yold)**2 + (v1z - v1zold)**2
    erri = erri / vhalf2
    converged = (erri < tolv2)

 enddo iterations

 if (its >= maxitsext) call fatal('update_grforce_leapfrog','VELOCITY ITERATIONS ON EXTERNAL FORCE NOT CONVERGED!!')

 fexti(1) = fextv(1)
 fexti(2) = fextv(2)
 fexti(3) = fextv(3)

 fxi = fxi + fexti(1)
 fyi = fyi + fexti(2)
 fzi = fzi + fexti(3)

end subroutine update_grforce_leapfrog


! !-----------------------------------------------------------------------
! !+
! !  writes input options to the input file
! !+
! !-----------------------------------------------------------------------
! subroutine write_options_ltforce(iunit)
!  use infile_utils, only:write_inopt
!  use physcon, only:pi
!  integer, intent(in) :: iunit
!
!  blackhole_spin_angle = blackhole_spin_angle*(180.0/pi)
!  write(iunit,"(/,a)") '# options relating to Lense-Thirring precession'
!  call write_inopt(blackhole_spin,'blackhole_spin','spin of central black hole (-1 to 1)',iunit)
!  call write_inopt(blackhole_spin_angle, &
!                  'blackhole_spin_angle','black hole spin angle w.r.t. x-y plane (0 to 180)',iunit)
!  blackhole_spin_angle = blackhole_spin_angle*(pi/180.0)
!
! end subroutine write_options_ltforce

! !-----------------------------------------------------------------------
! !+
! !  reads input options from the input file
! !+
! !-----------------------------------------------------------------------
! subroutine read_options_ltforce(name,valstring,imatch,igotall,ierr)
!  use io,      only:fatal
!  use physcon, only:pi
!  character(len=*), intent(in)  :: name,valstring
!  logical,          intent(out) :: imatch,igotall
!  integer,          intent(out) :: ierr
!  integer, save :: ngot = 0
!  character(len=30), parameter :: label = 'read_options_ltforce'
!
!  imatch  = .true.
!  igotall = .false.
!
!  select case(trim(name))
!  case('blackhole_spin')
!     read(valstring,*,iostat=ierr) blackhole_spin
!     if (blackhole_spin > 1 .or. blackhole_spin < -1.) then
!        call fatal(label,'invalid spin parameter for black hole')
!     endif
!     ngot = ngot + 1
!  case('blackhole_spin_angle')
!     read(valstring,*,iostat=ierr) blackhole_spin_angle
!     if (blackhole_spin_angle > 180. .or. blackhole_spin_angle < 0.) then
!        call fatal(label,'invalid spin angle for black hole (should be between 0 and 180 degrees)')
!     else
!        blackhole_spin_angle = blackhole_spin_angle*(pi/180.0)
!        sin_spinangle = sin(blackhole_spin_angle)
!        cos_spinangle = cos(blackhole_spin_angle)
!     endif
!  case default
!     imatch = .false.
!  end select
!
!  igotall = (ngot >= 1)
!
! end subroutine read_options_ltforce

end module extern_gr
