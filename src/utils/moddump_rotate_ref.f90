!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Solid body rotation for all particles to a given mean inclination and mean position angle
!
! :References: None
!
! :Owner: Antoine Alaguero
!
! :Runtime parameters:
!   - Omega     : *current disc position angle [deg]*
!   - incl      : *current disc inclination [deg]*
!   - ref_Omega : *reference position angle [deg]*
!   - ref_incl  : *reference inclination [deg]*
!
! :Dependencies: infile_utils, moddump_utils, part, physcon, prompting,
!   vectorutils
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters (angles in degrees)
 real :: incl      = 50.864   ! current disc inclination
 real :: Omega     = 34.36    ! current disc position angle
 real :: ref_incl  = 54.6     ! reference inclination (disc rotated into this plane)
 real :: ref_Omega = 53.0     ! reference position angle

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use physcon,              only:pi
 use part,                 only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use vectorutils,          only:rotatevec

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: alpha,gamma
 real                   :: a,b,c,d,e,f,g,h,i
 real                   :: temp_x,temp_y,temp_z
 real                   :: temp_vx,temp_vy,temp_vz
 real                   :: temp(3),temp_v(3)
 integer                :: j

 if (prompt_for_params) call read_interactive_moddumpfile()

 ! Rotation angles & coeffs
 alpha = (ref_incl - incl) *pi/180    !about x
 gamma = (ref_Omega - Omega) *pi/180     !about z

 !
 a = cos(alpha)
 b = -sin(alpha)*cos(gamma)
 c = sin(alpha)*sin(gamma)
 d = sin(alpha)
 e = cos(alpha)*cos(gamma)
 f = -cos(alpha)*sin(gamma)
 g = 0
 h = sin(gamma)
 i = cos(gamma)

 ! New positions & velocities
 do j = 1,npart
    temp(1) = xyzh(1,j)
    temp(2) = xyzh(2,j)
    temp(3) = xyzh(3,j)
    call rotatevec(temp, (/1.0,0.,0./), alpha)
    call rotatevec(temp, (/0.,0.,1.0/), gamma)

    temp_v(1) = vxyzu(1,j)
    temp_v(2) = vxyzu(2,j)
    temp_v(3) = vxyzu(3,j)
    call rotatevec(temp_v, (/1.0,0.,0./), alpha)
    call rotatevec(temp_v, (/0.,0.,1.0/), gamma)

!    temp_x = a*xyzh(1,j) + b*xyzh(2,j) + c*xyzh(3,j)
!    temp_y = d*xyzh(1,j) + e*xyzh(2,j) + f*xyzh(3,j)
!    temp_z = g*xyzh(1,j) + h*xyzh(2,j) + i*xyzh(3,j)
!
!    temp_vx = a*vxyzu(1,j) + b*vxyzu(2,j) + c*vxyzu(3,j)
!    temp_vy = d*vxyzu(1,j) + e*vxyzu(2,j) + f*vxyzu(3,j)
!    temp_vz = g*vxyzu(1,j) + h*vxyzu(2,j) + i*vxyzu(3,j)

    temp_x = temp(1)
    temp_y = temp(2)
    temp_z = temp(3)
    xyzh(1,j) = temp_x
    xyzh(2,j) = temp_y
    xyzh(3,j) = temp_z

    temp_vx = temp_v(1)
    temp_vy = temp_v(2)
    temp_vz = temp_v(3)
    vxyzu(1,j) = temp_vx
    vxyzu(2,j) = temp_vy
    vxyzu(3,j) = temp_vz

 enddo

 do j = 1,nptmass

    temp(1) = xyzmh_ptmass(1,j)
    temp(2) = xyzmh_ptmass(2,j)
    temp(3) = xyzmh_ptmass(3,j)
    call rotatevec(temp, (/1.0,0.,0./), alpha)
    call rotatevec(temp, (/0.,0.,1.0/), gamma)

    temp_v(1) = vxyz_ptmass(1,j)
    temp_v(2) = vxyz_ptmass(2,j)
    temp_v(3) = vxyz_ptmass(3,j)
    call rotatevec(temp_v, (/1.0,0.,0./), alpha)
    call rotatevec(temp_v, (/0.,0.,1.0/), gamma)

    temp_x = temp(1)
    temp_y = temp(2)
    temp_z = temp(3)
    xyzmh_ptmass(1,j) = temp_x
    xyzmh_ptmass(2,j) = temp_y
    xyzmh_ptmass(3,j) = temp_z

    temp_vx = temp_v(1)
    temp_vy = temp_v(2)
    temp_vz = temp_v(3)
    vxyz_ptmass(1,j) = temp_vx
    vxyz_ptmass(2,j) = temp_vy
    vxyz_ptmass(3,j) = temp_vz
!    temp_x = a*xyzmh_ptmass(1,j) + b*xyzmh_ptmass(2,j) + c*xyzmh_ptmass(3,j)
!    temp_y = d*xyzmh_ptmass(1,j) + e*xyzmh_ptmass(2,j) + f*xyzmh_ptmass(3,j)
!    temp_z = g*xyzmh_ptmass(1,j) + h*xyzmh_ptmass(2,j) + i*xyzmh_ptmass(3,j)
!
!    temp_vx = a*vxyz_ptmass(1,j) + b*vxyz_ptmass(2,j) + c*vxyz_ptmass(3,j)
!    temp_vy = d*vxyz_ptmass(1,j) + e*vxyz_ptmass(2,j) + f*vxyz_ptmass(3,j)
!    temp_vz = g*vxyz_ptmass(1,j) + h*vxyz_ptmass(2,j) + i*vxyz_ptmass(3,j)
!
!
!    xyzmh_ptmass(1,j) = temp_x
!    xyzmh_ptmass(2,j) = temp_y
!    xyzmh_ptmass(3,j) = temp_z
!
!    vxyz_ptmass(1,j) = temp_vx
!    vxyz_ptmass(2,j) = temp_vy
!    vxyz_ptmass(3,j) = temp_vz
!
 enddo

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter current disc inclination (deg)',incl)
 call prompt('Enter current disc position angle Omega (deg)',Omega)
 call prompt('Enter reference inclination (deg)',ref_incl)
 call prompt('Enter reference position angle Omega (deg)',ref_Omega)

end subroutine read_interactive_moddumpfile

subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 23
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(incl,'incl',db,errcount=nerr)
 call read_inopt(Omega,'Omega',db,errcount=nerr)
 call read_inopt(ref_incl,'ref_incl',db,errcount=nerr)
 call read_inopt(ref_Omega,'ref_Omega',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# reference-frame rotation angles (degrees)'
 call write_inopt(incl,'incl','current disc inclination [deg]',iunit)
 call write_inopt(Omega,'Omega','current disc position angle [deg]',iunit)
 call write_inopt(ref_incl,'ref_incl','reference inclination [deg]',iunit)
 call write_inopt(ref_Omega,'ref_Omega','reference position angle [deg]',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
