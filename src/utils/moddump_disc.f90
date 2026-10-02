!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! rndisc moddump routine: adds a warp in the disc/adds magnetic field
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - H_warp  : *warp width [code units]*
!   - HonR    : *disc aspect ratio H/R (must match the setup)*
!   - R_warp  : *warp radius [code units]*
!   - beta    : *plasma beta for the added toroidal field*
!   - incl    : *sine of inclination angle (0->1)*
!   - posangl : *position angle [deg]*
!
! :Dependencies: infile_utils, moddump_utils, part, prompting, setdisc
!
 use moddump_utils, only:prompt_for_params,write_moddump_header, &
                         init_moddump=>init_moddump_empty
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters (written to / read from the prefix.moddump file)
 real :: HonR    = 0.02   ! disc aspect ratio H/R (must match the setup)
 real :: R_warp  = 2.321  ! warp radius [code units]
 real :: H_warp  = 0.0    ! warp width [code units]
 real :: incl    = 0.5    ! sine of inclination angle (0->1)
 real :: posangl = 0.0    ! position angle [deg]
 real :: beta    = 10.0   ! plasma beta for the added toroidal field (MHD only)

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use setdisc,       only:set_incline_or_warp
 use part,          only:Bxyz,mhd,rho,igas

 integer, intent(in)    :: npartoftype(:)
 real,    intent(in)    :: massoftype(:)
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: npart_start_count,npart_tot,ii,i
 real    :: Bzero,pmassii,phi
 real    :: r2,r,omega,cs,pressure,psimax
 real    :: vphiold2,vphiold,vadd,vphicorr2

 if (prompt_for_params) call read_interactive_moddumpfile()

! Similar to that in set_disc
 npart_start_count=1
 npart_tot=npart
!
!---------------------------------------------
! Call setwarp to actually calculate the warp
 call set_incline_or_warp(xyzh,vxyzu,npart_tot,npart_start_count,posangl,incl,&
                          R_warp,H_warp,psimax)
!---------------------------------------------
 do i=npart_start_count,npart_tot
    xyzh(1,i)=xyzh(1,i)
    xyzh(2,i)=xyzh(2,i)
    xyzh(3,i)=xyzh(3,i)
    xyzh(4,i)=xyzh(4,i)

    vxyzu(1,i)=vxyzu(1,i)
    vxyzu(2,i)=vxyzu(2,i)
    vxyzu(3,i)=vxyzu(3,i)
 enddo
 print*,' Disc is now warped '

! Add magnetic field
 if (mhd) then

! Set up a magnetic field just in Bphi
    do ii = 1,npart
       r2 = xyzh(1,ii)**2 + xyzh(2,ii)**2 + xyzh(3,ii)**2
       r = sqrt(r2)
       phi = atan2(xyzh(2,ii),xyzh(1,ii))
       omega = r**(-1.5)
       cs = HonR*r*omega
       pmassii = massoftype(igas)
       pressure = cs**2*rho(ii)
       Bzero = sqrt(2.*pressure/beta)
       Bxyz(1,ii) = -Bzero*sin(phi)
       Bxyz(2,ii) = Bzero*cos(phi)

       ! Calculate correction in v_phi due to B
       vphiold = (-xyzh(2,ii)*vxyzu(1,ii) + xyzh(1,ii)*vxyzu(2,ii))/r
       vphiold2 = vphiold**2
       vphicorr2 = -2.*cs**2
!    if (vphicorr2 > vphi
       vadd = sqrt(vphiold2 + vphicorr2)
       vxyzu(1,ii) = vxyzu(1,ii) + sin(phi)*(vphiold - vadd)
       vxyzu(2,ii) = vxyzu(2,ii) - cos(phi)*(vphiold - vadd)

    enddo

    Bxyz(3,:) = 0.0

    print*,'Magnetic field added.'

! Set up poloidal magnetic field throughout (part of) the disc
!!  Bzero = sqrt(2.*polyk*rhosum/beta)
!!    do ii=1,npart
!!      if (abs(xyzh(3,ii)) < HonR) then  ! to only set the field up in a section of the disc
!!       theta=atan2(xyzh(2,ii),xyzh(1,ii))
!!       Bxyz(1,ii) = 0. !real(Bzero*sin(theta),kind=4)
!!       Bxyz(2,ii) = 0. !real(-Bzero*cos(theta),kind=4)
!!       Bxyz(3,ii) = 0.
!!      endif
!!    enddo
 endif

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt
 use part,      only:mhd

 call prompt('Enter disc aspect ratio H/R (must match the setup)',HonR,0.)
 call prompt('Enter warp radius R_warp (code units)',R_warp,0.)
 call prompt('Enter warp width H_warp (code units)',H_warp,0.)
 call prompt('Enter sine of inclination angle (0->1)',incl,0.,1.)
 call prompt('Enter position angle (deg)',posangl)
 if (mhd) call prompt('Enter plasma beta for added toroidal field',beta,0.)

end subroutine read_interactive_moddumpfile

subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use part,         only:mhd
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 23
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(HonR,'HonR',db,errcount=nerr,min=0.)
 call read_inopt(R_warp,'R_warp',db,errcount=nerr,min=0.)
 call read_inopt(H_warp,'H_warp',db,errcount=nerr,min=0.)
 call read_inopt(incl,'incl',db,errcount=nerr,min=0.,max=1.)
 call read_inopt(posangl,'posangl',db,errcount=nerr)
 if (mhd) call read_inopt(beta,'beta',db,errcount=nerr,min=0.)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 use part,         only:mhd
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# disc warp parameters'
 call write_inopt(HonR,'HonR','disc aspect ratio H/R (must match the setup)',iunit)
 call write_inopt(R_warp,'R_warp','warp radius [code units]',iunit)
 call write_inopt(H_warp,'H_warp','warp width [code units]',iunit)
 call write_inopt(incl,'incl','sine of inclination angle (0->1)',iunit)
 call write_inopt(posangl,'posangl','position angle [deg]',iunit)
 if (mhd) call write_inopt(beta,'beta','plasma beta for the added toroidal field',iunit)
 close(iunit)

end subroutine write_moddump

end module moddump
