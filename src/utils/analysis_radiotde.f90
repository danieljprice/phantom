!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
!   - drad_cap  : *capture thickness (in cm) (-ve for all particles at outer radius)*
!   - phi_max   : *max phi (in deg)*
!   - phi_min   : *min phi (in deg)*
!   - rad_cap   : *capture inner radius (in cm)*
!   - theta_max : *max theta (in deg)*
!   - theta_min : *min theta (in deg)*
!   - v_max     : *max velocity (in c)*
!   - v_min     : *min velocity (in c)*
!
! :Dependencies: infile_utils, io, physcon, readwrite_dumps, units
!
 implicit none
 character(len=8), parameter, public :: analysistype = 'radiotde'
 public :: do_analysis

 private

 real, dimension(:), allocatable    :: theta,plot_theta,phi,vr,vtheta,vphi
 logical, dimension(:), allocatable :: cap

 !---- These can be changed in the params file
 real    :: rad_cap = 1.e16 ! radius where the outflow in captured (in cm)
 real    :: drad_cap = 4.7267e14 ! thickness of the shell to capture outflow (in cm)
 real    :: v_min = 0.
 real    :: v_max = 1.
 real    :: theta_min = -180.
 real    :: theta_max = 180.
 real    :: phi_min = -90.
 real    :: phi_max = 90.
 real    :: m_accum, m_cap, vr_accum_mean, vr_cap_mean, v_accum_mean, v_cap_mean, e_accum, e_cap
 integer :: n_accum, n_cap

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use readwrite_dumps, only: opened_full_dump
 use units,           only: utime,udist,unit_energ,umass
 use physcon,         only: solarm,days
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 character(len=120) :: output
 character(len=30)  :: filename
 integer            :: i,ierr
 logical            :: iexist
 real               :: toMsun,todays

 m_accum = 0.
 n_accum = 0
 m_cap = 0.
 n_cap = 0
 e_accum = 0.
 e_cap = 0.

 toMsun = umass/solarm
 todays = utime/days

 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a8,i5.5)") 'outflow_',numfile
 write(*,'("Output file name is ",A)') output

 ! Read black hole mass from params file
 filename = 'analysis_'//trim(analysistype)//'.params'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_tdeparams(filename,ierr)
 if (.not.iexist.or.ierr/=0) then
    call write_tdeparams(filename)
    print*,' Edit '//trim(filename)//' and rerun phantomanalysis'
    stop
 endif

 rad_cap = rad_cap/udist
 if (drad_cap < 0.) then
    drad_cap = huge(0.)
 else
    drad_cap = drad_cap/udist
 endif
 print*, 'Capture particles from', rad_cap, 'to', rad_cap+drad_cap

 allocate(theta(npart),plot_theta(npart),phi(npart),vr(npart),vtheta(npart),vphi(npart),cap(npart))
 cap = .false.

 call tde_analysis(npart,pmass,xyzh,vxyzu)

 if (n_cap > 0) then
    open(iunit,file=output)
    write(iunit,'("# ",es20.12,"   # TIME")') time
    write(iunit,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'theta',      &
       2,'thetap',  &
       3,'phi', &
       4,'vr',     &
       5,'vtheta',   &
       6,'vphi'

    do i = 1,npart
       if (cap(i)) then
          write(iunit,'(6(es18.10,1X))') &
          theta(i),   &
          plot_theta(i), &
          phi(i),    &
          vr(i),   &
          vtheta(i),    &
          vphi(i)
       endif
    enddo
    close(iunit)
 endif

 deallocate(theta,plot_theta,phi,vr,vtheta,vphi,cap)

 inquire(file='outflows',exist=iexist)
 if (iexist) then
    open(iunit,file='outflows',status='old',position='append')
 else
    open(iunit,file='outflows',status='new')
    write(iunit,'(9(A15,1X))') '# time', 'm_cap[msun]', 'm_accum[msun]', 'vr_accum_mean[c]', 'vr_cap_mean[c]', &
                               'v_accum_mean[c]', 'v_cap_mean[c]', 'e_accum[erg]', 'e_cap[erg]'
 endif
 write(iunit,'(9(es18.10,1x))') &
    time*todays, &
    m_cap*toMsun, &
    m_accum*toMsun, &
    vr_accum_mean, &
    vr_cap_mean, &
    v_accum_mean, &
    v_cap_mean, &
    e_accum*unit_energ, &
    e_cap*unit_energ
 close(iunit)

 write(*,'(I8,1X,A2,1X,I8,1X,A34)') n_cap, 'of', npart, 'particles are in the capture shell'
 write(*,'(I8,1X,A2,1X,I8,1X,A40)') n_accum, 'of', npart, 'particles are outside the capture radius'

end subroutine do_analysis

!--------------------------------------------------------------------------------------------------------------------
!
!-- Actual subroutine where the analysis is done!
!
!--------------------------------------------------------------------------------------------------------------------
subroutine tde_analysis(npart,pmass,xyzh,vxyzu)
 integer, intent(in) :: npart
 real, intent(in)    :: pmass,xyzh(:,:),vxyzu(:,:)
 integer :: i
 real    :: r,v,x,y,z,xyz(1:3),vx,vy,vz,vxyz(1:3)
 real    :: thetai,phii,vri
 real    :: vr_accum_add,vr_cap_add,v_accum_add,v_cap_add

 vr_accum_add = 0.
 vr_cap_add = 0.
 v_accum_add = 0.
 v_cap_add = 0.

 do i = 1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    xyz = (/x,y,z/)
    vx = vxyzu(1,i)
    vy = vxyzu(2,i)
    vz = vxyzu(3,i)
    vxyz = (/vx,vy,vz/)
    r = sqrt(dot_product(xyz,xyz))
    v = sqrt(dot_product(vxyz,vxyz))
    if (r > rad_cap) then
       m_accum = m_accum + pmass
       n_accum = n_accum + 1
       e_accum = e_accum + 0.5*pmass*v**2
       vri = dot_product(vxyz,xyz)/r
       vr_accum_add = vr_accum_add + vri
       v_accum_add = v_accum_add + v
       if (r-rad_cap < drad_cap .and. (v  >=  v_min .and. v  <=  v_max)) then
          thetai = atan2d(y,x)
          phii = atan2d(z,sqrt(x**2+y**2))
          if ((thetai  >=  theta_min .and. thetai  <=  theta_max) .and. (phii  >=  phi_min .and. phii  <=  phi_max)) then
             m_cap = m_cap + pmass
             n_cap = n_cap + 1
             cap(i) = .true.
             theta(i) = thetai
             phi(i) = phii
             plot_theta(i) = theta(i) * sqrt(cosd(phi(i)))
             vr(i) = vri
             vtheta(i) = -sind(phii)*vx + cosd(phii)*vy
             vphi(i) = cosd(thetai)*cosd(phii)*vx + cosd(thetai)*sind(phii)*vy - sind(thetai)*vz
             e_cap = e_cap + 0.5*pmass*v**2
             vr_cap_add = vr_cap_add + vri
             v_cap_add = v_cap_add + v
          endif
       endif
    endif
 enddo
 vr_accum_mean = vr_accum_add/n_accum
 v_accum_mean = v_accum_add/n_accum
 vr_cap_mean = vr_cap_add/n_cap
 v_cap_mean = v_cap_add/n_cap

end subroutine tde_analysis

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
 write(iunit,"(a,/)") '# options when performing radio TDE analysis'

 call write_inopt(rad_cap,'rad_cap','capture inner radius (in cm)',iunit)
 call write_inopt(drad_cap,'drad_cap','capture thickness (in cm) (-ve for all particles at outer radius)',iunit)

 call write_inopt(v_min,'v_min','min velocity (in c)',iunit)
 call write_inopt(v_max,'v_max','max velocity (in c)',iunit)

 call write_inopt(theta_min,'theta_min','min theta (in deg)',iunit)
 call write_inopt(theta_max,'theta_max','max theta (in deg)',iunit)

 call write_inopt(phi_min,'phi_min','min phi (in deg)',iunit)
 call write_inopt(phi_max,'phi_max','max phi (in deg)',iunit)

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

 print "(a)",'reading analysis options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(rad_cap,'rad_cap',db,min=0.,errcount=nerr)
 call read_inopt(drad_cap,'drad_cap',db,errcount=nerr)

 call read_inopt(v_min,'v_min',db,min=0.,max=1.,errcount=nerr)
 call read_inopt(v_max,'v_max',db,min=0.,max=1.,errcount=nerr)

 call read_inopt(theta_min,'theta_min',db,min=-180.,max=180.,errcount=nerr)
 call read_inopt(theta_max,'theta_max',db,min=-180.,max=180.,errcount=nerr)

 call read_inopt(phi_min,'phi_min',db,min=-90.,max=90.,errcount=nerr)
 call read_inopt(phi_max,'phi_max',db,min=-90.,max=90.,errcount=nerr)

 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of params file: re-writing...'
    ierr = nerr
 endif

end subroutine read_tdeparams

end module analysis

