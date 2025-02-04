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
!   - drad_cap  : *capture thickness (in cm) (-ve for all particles at outer radius)*
!   - npart_tde : *npart in tde sims (-ve=10*npart of cnm)*
!   - phi_max   : *max phi (in deg)*
!   - phi_min   : *min phi (in deg)*
!   - rad_cap   : *capture inner radius (in cm)*
!   - theta_max : *max theta (in deg)*
!   - theta_min : *min theta (in deg)*
!   - v_max     : *max velocity (in c)*
!   - v_min     : *min velocity (in c)*
!
! :Dependencies: infile_utils, io, part, physcon, readwrite_dumps, units
!
 implicit none
 character(len=8), parameter, public :: analysistype = 'radiotde'
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

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use readwrite_dumps, only: opened_full_dump
 use units,           only: utime,udist,unit_energ,umass!,unit_density
 use physcon,         only: solarm,days
 use part,            only: pxyzu
 character(len=*),   intent(in) :: dumpfile
 integer,            intent(in) :: numfile,npart,iunit
 real,               intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,               intent(in) :: pmass,time
 character(len=120) :: output
 character(len=30)  :: filename
 integer            :: i,ierr,npart_new,npart_tde_old
 logical            :: iexist
 real               :: toMsun,todays

 m_accum = 0.
 n_accum = 0
 m_cap = 0.
 n_cap = 0
 e_accum = 0.
 e_cap = 0.
 ana = 'shock'

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

 ! read background entropy
 if (npart_cnm < 0) then
    if (npart_tde_reserve < 0) npart_tde_reserve = 10*npart
    allocate(ent_bg(npart_tde_reserve+npart)) ! save more memory for later injection
    npart_cnm = npart
    call record_background(pxyzu(4,:),0,npart,ent_bg)
    write(*,'(I9,1x,a16)') npart_cnm, 'particles in CNM'
 endif
! not meaningful and will not do anything if cut-and-put
 npart_tde_old = npart_tde
 npart_tde = npart - npart_cnm
 npart_new = npart_tde - npart_tde_old
 if (npart_new > 0) call record_background(pxyzu(4,:),npart_tde_old+npart_cnm,npart_new,ent_bg)

 ! allocate memory
 allocate(rad_all(npart),vr_all(npart),v_all(npart))
 call to_rad(npart,xyzh,vxyzu,rad_all,vr_all,v_all)

 select case (trim(ana))
 case ('outflow')
    write(*,'(a)') ' Analysing the outflow ...'
    write(output,"(a8,i5.5)") 'outflow_',numfile
    write(*,'(" Output file name is ",A)') output

    rad_cap = rad_cap/udist
    if (drad_cap < 0.) then
       drad_cap = huge(0.)
    else
       drad_cap = drad_cap/udist
    endif
    print*, 'Capture particles from', rad_cap, 'to', rad_cap+drad_cap

    allocate(theta(npart),plot_theta(npart),phi(npart),vr(npart),vtheta(npart), &
             vphi(npart),cap(npart))
    cap = .false.

    call outflow_analysis(npart,pmass,xyzh,vxyzu,rad_all,vr_all,v_all)

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

    deallocate(theta,plot_theta,phi,vr,vtheta,vphi,cap,rad_all,vr_all,v_all)

    inquire(file='outflows',exist=iexist)
    if (iexist) then
       open(iunit,file='outflows',status='old',position='append')
    else
       open(iunit,file='outflows',status='new')
       write(iunit,'(14(A,1X))') '#', 'time', 'm_cap[msun]', 'm_accum[msun]', 'vr_accum_mean[c]', 'vr_accum_max[c]', &
                                 'r_accum_maxv[cm]', 'vr_cap_mean[c]', 'vr_cap_max[c]', 'r_cap_maxv[cm]',  &
                                 'v_accum_mean[c]', 'v_cap_mean[c]', 'e_accum[erg]', 'e_cap[erg]'
    endif
    write(iunit,'(13(es18.10,1x))') &
        time*todays, &
        m_cap*toMsun, &
        m_accum*toMsun, &
        vr_accum_mean, &
        vr_accum_max, &
        r_accum_maxv*udist, &
        vr_cap_mean, &
        vr_cap_max, &
        r_cap_maxv*udist, &
        v_accum_mean, &
        v_cap_mean, &
        e_accum*unit_energ, &
        e_cap*unit_energ
    close(iunit)

    write(*,'(I8,1X,A2,1X,I8,1X,A34)') n_cap, 'of', npart, 'particles are in the capture shell'
    write(*,'(I8,1X,A2,1X,I8,1X,A40)') n_accum, 'of', npart, 'particles are outside the capture radius'

 case ('shock')
    write(*,'(a)') ' Analysing the shock ...'

    call shock_analysis(npart,pmass,rad_all,vr_all,pxyzu(4,:))

    deallocate(rad_all,vr_all,v_all)

    inquire(file='shock',exist=iexist)
    if (iexist) then
       open(iunit,file='shock',status='old',position='append')
    else
       open(iunit,file='shock',status='new')
       write(iunit,'(17(A,1x))') '#', 'time', 'rad_min[cm]', 'rad_max[cm]', 'velocity[c]', 'mass[Msun]', 'energy[erg]', & !'density[g/cm-3]'
                                             'rad_min_tde[cm]', 'rad_max_tde[cm]', 'vel_tde[c]', 'mass_tde[Msun]', 'ene_tde[erg]', &
                                             'rad_min_cnm[cm]', 'rad_max_cnm[cm]', 'vel_cnm[c]', 'mass_cnm[Msun]', 'ene_cnm[erg]'
    endif
    if (rad_max > 0.) then
       write(iunit,'(16(es18.10,1x))') &
       time*todays, &
       rad_min*udist, rad_max*udist, shock_v, shock_m*umass/solarm, shock_e*unit_energ, &
       rad_min_tde*udist, rad_max_tde*udist, shock_v_tde, shock_m_tde*umass/solarm, shock_e_tde*unit_energ, &
       rad_min_cnm*udist, rad_max_cnm*udist, shock_v_cnm, shock_m_cnm*umass/solarm, shock_e_cnm*unit_energ !shock_rho*unit_density
    endif
    close(iunit)

 case default
    write(*,'(a)') " Unknown analysis type. Do 'outflow' or 'shock'"
    stop
 end select

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
subroutine outflow_analysis(npart,pmass,xyzh,vxyzu,rad_all,vr_all,v_all)
 integer, intent(in) :: npart
 real, intent(in)    :: pmass,xyzh(:,:),vxyzu(:,:),rad_all(:),vr_all(:),v_all(:)
 integer :: i
 real    :: r,v,x,y,z,vx,vy,vz
 real    :: thetai,phii,vri
 real    :: vr_accum_add,vr_cap_add,v_accum_add,v_cap_add

 vr_accum_add = 0.
 vr_cap_add = 0.
 v_accum_add = 0.
 v_cap_add = 0.
 vr_accum_max = 0.
 vr_cap_max = 0.

 do i = 1,npart
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    vx = vxyzu(1,i)
    vy = vxyzu(2,i)
    vz = vxyzu(3,i)
    r = rad_all(i)
    v = v_all(i)
    if (r > rad_cap) then
       m_accum = m_accum + pmass
       n_accum = n_accum + 1
       e_accum = e_accum + 0.5*pmass*v**2
       vri = vr_all(i)
       vr_accum_add = vr_accum_add + vri
       v_accum_add = v_accum_add + v
       if (vri > vr_accum_max) then
          vr_accum_max = vri
          r_accum_maxv = r
       endif
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
             if (vri > vr_cap_max) then
                vr_cap_max = vri
                r_cap_maxv = r
             endif
          endif
       endif
    endif
 enddo
 vr_accum_mean = vr_accum_add/n_accum
 v_accum_mean = v_accum_add/n_accum
 vr_cap_mean = vr_cap_add/n_cap
 v_cap_mean = v_cap_add/n_cap

end subroutine outflow_analysis

subroutine record_background(ent,npart_old,npart_new,ent_bg)
 real, intent(in)    :: ent(:)
 integer, intent(in) :: npart_old,npart_new
 real, intent(inout) :: ent_bg(:)
 integer, parameter  :: iunit=235
 integer             :: i

 print*, 'Record background entropy of ', npart_new, ' particles'

 do i=1,npart_new
    ent_bg(npart_old+i) = ent(npart_old+i)*1.1 ! give some range for self evolution
    !(is there a reasonable choice instead of arbitrary?)
 enddo

end subroutine record_background

subroutine shock_analysis(npart,pmass,rad_all,vr_all,ent)
 use units,   only: udist
 use physcon, only: au,pi
 integer, intent(in) :: npart
 real, intent(in) :: pmass,rad_all(:),vr_all(:),ent(:)
 integer :: i,n,n_cnm,n_tde
 real    :: ri,half_m,ei,vi
 !
 !------Determine the shock
 !
 n = 0
 n_cnm = 0.
 n_tde = 0.
 shock_e = 0.
 shock_e_cnm = 0.
 shock_e_tde = 0.
 shock_v = 0. ! take max vel
 shock_v_cnm = 0.
 shock_v_tde = 0.
 rad_max = 0.
 rad_max_cnm = 0.
 rad_max_tde = 0.
 rad_min = huge(0.)
 rad_min_cnm = huge(0.)
 rad_min_tde = huge(0.)
 half_m = pmass*0.5

 do i = 1,npart
    if (ent(i) > ent_bg(i)) then
       ri = rad_all(i)
       vi = vr_all(i)
       ei = half_m*vi**2
       n = n + 1
       if (vi > shock_v) shock_v = vi
       if (ri < rad_min) rad_min = ri
       if (ri > rad_max) rad_max = ri
       shock_e = shock_e + ei

       if (i > npart_cnm) then
          ! tde outflow
          n_tde = n_tde + 1
          if (vi > shock_v_tde) shock_v_tde = vi
          if (ri < rad_min_tde) rad_min_tde = ri
          if (ri > rad_max_tde) rad_max_tde = ri
          shock_e_tde = shock_e_tde + ei
       else
          ! cnm
          n_cnm = n_cnm + 1
          if (vi > shock_v_cnm) shock_v_cnm = vi
          if (ri < rad_min_cnm) rad_min_cnm = ri
          if (ri > rad_max_cnm) rad_max_cnm = ri
          shock_e_cnm = shock_e_cnm + ei
       endif
    endif
 enddo

 write(*,'(a14,1x,es8.1,1x,a5,1x,es8.1,1x,a2)') ' Shock is from', rad_min*udist/au, 'au to', rad_max*udist/au, 'au'

 shock_m = shock_e*2./shock_v**2  !pmass*n
 shock_m_cnm = shock_e_cnm*2./shock_v_cnm**2 !pmass*n_cnm
 shock_m_tde = shock_e_tde*2./shock_v_tde**2 !pmass*n_tde
 !shock_rho = shock_m*4./3.*pi*(rad_max**3-rad_min**3)

end subroutine shock_analysis

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
 call write_inopt(ana,'analysis',"analysis type: 'outflow' or 'shock'",iunit)

 select case (trim(ana))
 case ('outflow')
    call write_inopt(rad_cap,'rad_cap','capture inner radius (in cm)',iunit)
    call write_inopt(drad_cap,'drad_cap','capture thickness (in cm) (-ve for all particles at outer radius)',iunit)

    call write_inopt(v_min,'v_min','min velocity (in c)',iunit)
    call write_inopt(v_max,'v_max','max velocity (in c)',iunit)

    call write_inopt(theta_min,'theta_min','min theta (in deg)',iunit)
    call write_inopt(theta_max,'theta_max','max theta (in deg)',iunit)

    call write_inopt(phi_min,'phi_min','min phi (in deg)',iunit)
    call write_inopt(phi_max,'phi_max','max phi (in deg)',iunit)
 case ('shock')
    call write_inopt(npart_tde_reserve,'npart_tde','npart in tde sims (-ve=10*npart of cnm)',iunit)
 case default
 end select

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

 call read_inopt(ana,'analysis',db,errcount=nerr)

 select case (trim(ana))
 case ('outflow')
    call read_inopt(rad_cap,'rad_cap',db,min=0.,errcount=nerr)
    call read_inopt(drad_cap,'drad_cap',db,errcount=nerr)

    call read_inopt(v_min,'v_min',db,min=0.,max=1.,errcount=nerr)
    call read_inopt(v_max,'v_max',db,min=0.,max=1.,errcount=nerr)

    call read_inopt(theta_min,'theta_min',db,min=-180.,max=180.,errcount=nerr)
    call read_inopt(theta_max,'theta_max',db,min=-180.,max=180.,errcount=nerr)

    call read_inopt(phi_min,'phi_min',db,min=-90.,max=90.,errcount=nerr)
    call read_inopt(phi_max,'phi_max',db,min=-90.,max=90.,errcount=nerr)
 case ('shock')
    call read_inopt(npart_tde_reserve,'npart_tde',db,errcount=nerr)
 case default
 end select

 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of params file: re-writing...'
    ierr = nerr
 endif

end subroutine read_tdeparams

end module analysis

