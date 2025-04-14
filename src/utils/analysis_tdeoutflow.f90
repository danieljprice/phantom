!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Computes the outflow profile in a TDE simulation
!
! :References: None
!
! :Owner: Fitz) Hu
!
! :Runtime parameters:
!   - phi_max   : *max phi (in deg) (-ve = ignore)*
!   - phi_min   : *min phi (in deg) (-ve = ignore)*
!   - r_in      : *radius to count outflow (in cm)*
!   - theta_max : *max theta (in deg) (-ve = ignore)*
!   - theta_min : *min theta (in deg) (-ve = ignore)*
!
! :Dependencies: infile_utils, io, part, physcon, readwrite_dumps, units
!
 implicit none
 character(len=10), parameter, public :: analysistype = 'tdeoutflow'
 public :: do_analysis

 private

 character(len=7) :: ana
 real, dimension(:), allocatable    :: rad_all,vr_all,v_all
 real, dimension(:), allocatable    :: theta,plot_theta,phi,vr,vtheta,vphi
 logical, dimension(:), allocatable :: cap
 real    :: m_accum, m_cap
 real    :: vr_accum_mean, vr_accum_max, vr_cap_mean, vr_cap_max
 real    :: r_accum_maxv, r_cap_maxv
 real    :: v_accum_mean, v_cap_mean
 real    :: e_accum, e_cap
 integer :: n_accum, n_cap
 real    :: shock_v, rad_min, rad_max, shock_e, shock_m!, shock_rho
 real    :: shock_v_tde, rad_min_tde, rad_max_tde, shock_e_tde, shock_m_tde!, shock_rho
 real    :: shock_v_cnm, rad_min_cnm, rad_max_cnm, shock_e_cnm, shock_m_cnm!, shock_rho

 !---- These can be changed in the params file
 real    :: rad_cap = 1.e16 ! radius where the outflow in captured (in cm)
 real    :: drad_cap = 4.7267e14 ! thickness of the shell to capture outflow (in cm)
 real    :: v_min = 0.
 real    :: v_max = 1.
 real    :: theta_min = -180.
 real    :: theta_max = 180.
 real    :: phi_min = -90.
 real    :: phi_max = 90.

 !--- shock detection global var
 integer           :: npart_cnm = -1, npart_tde, npart_tde_reserve=-1
 real, allocatable :: ent_bg(:)
 logical, allocatable :: counted(:),accreted(:)
 real :: told,r_in=1.e14
 logical :: first = .true.

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use readwrite_dumps, only: opened_full_dump
 use units,           only: udist,utime,unit_energ,umass!,unit_density
 use physcon,         only: solarm,days,c
 use part,            only: pxyzu
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 character(len=120) :: output
 character(len=30)  :: filename,outfile
 integer            :: i,ierr,npart_new,npart_tde_old
 logical            :: iexist
 real               :: toMsun,todays,dt
 real :: mout,vrout,vout,macc

 toMsun = umass/solarm
 todays = utime/days

 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif

! Print the analysis being done
 write(*,'(" Performing analysis type ",A)') analysistype
 write(*,'(" Input file name is ",A)') dumpfile

 ! Read black hole mass from params file
 filename = 'analysis_'//trim(analysistype)//'.params'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_tdeparams(filename,ierr)
 if (.not.iexist.or.ierr/=0) then
    call write_tdeparams(filename)
    print*,' Edit '//trim(filename)//' and rerun phantomanalysis'
    stop
 endif

 ! input to code unit
 r_in = r_in / udist

 ! allocate memory
 if (allocated(rad_all)) deallocate(rad_all(npart),vr_all(npart),v_all(npart))
 allocate(rad_all(npart),vr_all(npart),v_all(npart))
 call to_rad(npart,xyzh,vxyzu,rad_all,vr_all,v_all)

 write(*,'(a)') ' Analysing the outflow ...'

 print*, 'Counting outflow from', r_in

 if (first) then
    allocate(counted(npart),accreted(npart))
    counted = .false.
    accreted = .false.
    mout = 0.
    vrout = 0.
    vout = 0.
    macc = 0.
    dt = 1.
 else
    call outflow_analysis(npart,pmass,xyzh,vxyzu,rad_all,vr_all,v_all,mout,vrout,vout,macc)
    dt = time - told
 endif
 told = time

 outfile='outflow'
 inquire(file=outfile,exist=iexist)
 if (iexist .and. .not. first) then
    open(iunit,file=outfile,status='old',access='append')
 elseif (iexist) then
    open(iunit,file=outfile,status='replace')
 else
    open(iunit,file=outfile,status='new')
 endif

 if (first) then
    write(iunit,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
           1,'time [s]',    &
           2,'mout [g/s]',      &
           3,'vrout [cm/s]',  &
           4,'vout [cm/s]',  &
           5,'macc [g/s]'
 endif

 write(iunit,'(5(es18.10,1X))') &
           time*utime, &
           mout/dt*umass/utime,   &
           vrout, &
           vout, &
           macc/dt*umass/utime
 close(iunit)
 first = .false.

end subroutine do_analysis

subroutine to_rad(npart,xyzh,vxyzu,rad,vr,v)
 integer, intent(in) :: npart
 real, intent(in) :: xyzh(:,:),vxyzu(:,:)
 real, intent(out) :: rad(:),vr(:),v(:)
 integer :: i
 real :: xyz(1:3),vxyz(1:3)

 do i = 1,npart
    xyz = xyzh(1:3,i)
    vxyz = vxyzu(1:3,i)
    rad(i) = sqrt(dot_product(xyz,xyz))
    vr(i) = dot_product(xyz,vxyz)/rad(i)
    v(i) = sqrt(dot_product(vxyz,vxyz))
 enddo

end subroutine to_rad
!--------------------------------------------------------------------------------------------------------------------
!
!-- Actual subroutine where the analysis is done!
!
!--------------------------------------------------------------------------------------------------------------------
subroutine outflow_analysis(npart,pmass,xyzh,vxyzu,rad_all,vr_all,v_all,mout,vrout,vout,macc)
 use io, only: fatal
 use part, only:isdead_or_accreted
 integer, intent(in) :: npart
 real, intent(in)    :: pmass,xyzh(:,:),vxyzu(:,:),rad_all(:),vr_all(:),v_all(:)
 real, intent(out)   :: mout,vrout,vout,macc
 integer :: i,nout,nacc
 real    :: ri,vi,x,y,z
 real    :: thetai,phii,vri
 real    :: vrsum,vsum

 nout = 0
 nacc = 0
 vrsum = 0.
 vsum = 0.

 do i = 1,npart
    ri = rad_all(i)
    vi = v_all(i)
    vri = vr_all(i)
    if (isdead_or_accreted(xyzh(4,i))) then
       nacc = nacc + 1
       accreted(i) = .true.
    elseif (ri > r_in) then
       if (.not. counted(i)) then
          if (theta_min < -180. .or. theta_min > 180.) theta_min = -180.
          if (theta_max < theta_min .or. theta_max > 180.) theta_max = 180.
          if (phi_min < -90. .or. phi_min > 90.) phi_min = -90.
          if (phi_max < phi_min .or. phi_max > 90.) phi_max = 90.

          x = xyzh(1,i)
          y = xyzh(2,i)
          z = xyzh(3,i)
          thetai = atan2d(y,x)
          phii = atan2d(z,sqrt(x**2+y**2))

          if ((thetai >= theta_min .and. thetai <= theta_max) .and. (phii >= phi_min .and. phii <= phi_max)) then
             nout = nout + 1
             vrsum = vrsum + vri
             vsum = vsum + vi
          endif
          counted(i) = .true.
       endif
    else
       counted(i) = .false.
    endif
 enddo
 mout = nout * pmass
 vrout = vrsum / nout
 vout = vsum / nout
 macc = nacc * pmass

end subroutine outflow_analysis

!----------------------------------------------------------------
!+
!  Read/write tde information from/to params file
!+
!----------------------------------------------------------------
subroutine write_tdeparams(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing analysis options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a,/)") '# options when performing TDE outflow analysis'

 call write_inopt(r_in,'r_in','radius to count outflow (in cm)',iunit)

 call write_inopt(theta_min,'theta_min','min theta (in deg) (-ve = ignore)',iunit)
 call write_inopt(theta_max,'theta_max','max theta (in deg) (-ve = ignore)',iunit)

 call write_inopt(phi_min,'phi_min','min phi (in deg) (-ve = ignore)',iunit)
 call write_inopt(phi_max,'phi_max','max phi (in deg) (-ve = ignore)',iunit)

 close(iunit)

end subroutine write_tdeparams

subroutine read_tdeparams(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter        :: iunit = 21
 integer                   :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading analysis options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(r_in,'r_in',db,min=0.,errcount=nerr)

 call read_inopt(theta_min,'theta_min',db,max=360.,errcount=nerr)
 call read_inopt(theta_max,'theta_max',db,max=360.,errcount=nerr)

 call read_inopt(phi_min,'phi_min',db,max=180.,errcount=nerr)
 call read_inopt(phi_max,'phi_max',db,max=180.,errcount=nerr)

 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of params file: re-writing...'
    ierr = nerr
 endif

end subroutine read_tdeparams

end module analysis

