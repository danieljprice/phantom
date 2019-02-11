!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine calculating the mass and radius of a disc
!  As per Price & Bate (2007), Wurster, Price & Bate (2016), disc particles
!  are those with rho > 1e-13 g cm^-3 & the radius contains 99% of this mass.
!  Generalised for gas, dust and stellar discs.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, dim, eos, infile_utils, io, kernel, nicil,
!    options, part, physcon, sortutils, units
!+
!--------------------------------------------------------------------------
module analysis
 use dim,         only: maxp,maxvxyzu,mhd_nonideal
 use options,     only: alphaB
 use part,        only: maxptmass,n_R,n_electronT
 use part,        only: isdead_or_accreted,iamtype,iphase,igas,massoftype,maxphase,rhoh
 use eos,         only: ieos,init_eos,equationofstate,init_eos, &
                        get_spsound,get_temperature,get_temperature_from_ponrho
 use nicil,       only: nicil_initialise,nicil_get_ion_n,nicil_get_eta,unit_eta,n_data_out
 use physcon,     only: pi
 implicit none
 character(len=20), parameter, public :: analysistype = 'discRM'
 public :: do_analysis
 !--Free parameters
 real,    private, parameter :: dthreshcgs  = 1.0e-13 ! density (g cm^-3) above which gas particles are considered 'in' the disc
 real,    private, parameter :: massfrac    = 0.99    ! disc radius which contains massfrac percent of the mass
 real,    private, parameter :: angle       = 20.0    ! Use all gas within angle of the midplane; if <=0, then use dthreshcgs
 real,    private, parameter :: dringthresh = 0.1     ! density threshold relative to the max density in the disc to define its radius
 real,    private            :: rmax_au     = 128.0   ! Maximum radial distance of interest
 real,    private            :: rmin_au     = 0.01    ! Minimum radial distance of interest (only for log r
 real,    private, parameter :: rmergeAU    =  12.0   ! if sink separation is smaller, sinks are treated as one
 integer, private, parameter :: rthreshAU   = 128     ! radius (AU) of sphere in which to gather statistics (is integer since it used in filenames)
 integer,          parameter :: nbins       = 256     ! number of bins use to calculate disc properties
 integer, private, parameter :: ninring     = 116     ! number of particles in a ring when determine radial cutoff of disc (relatively insensitive to results)
 real,    private            :: rcom_au     = 64.0    ! Radius within which will be used to reset the CoM
 logical, private            :: map_all_R   = .true.  ! when radially binning, use all particles above the density threshhold (true)
 ! or just those in the disc (false)
 logical, private            :: discRM_only = .false. ! only calculate the disc mass and radius
 logical, private            :: ignor_lowrho= .true.  ! only calculate radial profiles if rho > dthresh
 logical, private            :: printrad    = .true.  ! Print radial profile at each dump
 logical, private            :: printvol    = .false. ! Print time evolution of the volume characteristics
 logical, private            :: ignore_sinks= .false. ! Assume one object only, and never centre on sink particles
 logical, private            :: sinks_only  = .false. ! Only perform analysis if sinks are present
 logical, private            :: sinks_two   = .true.  ! Perform analysis on only first two sinks (if ignore_sinks=.false.)
 logical, private            :: calc_mu     = .true.  ! Calculate the evolution of the mass to flux ratio (independent of discRM_only)
 logical, private            :: calc_mu_verb= .true.  ! If the above is true, output M(R) & <B> data to file
 logical, private            :: calc_eta    = .true.  ! Calculate global eta values on rho > rhocrit (independent of discRM_only)
 logical, private            :: use_etaart_old = .true. ! Use the original artificial resistivity as used in WPB2016
 logical, private            :: log_rbin    = .true.  ! use logarithmic radial bins
 ! parameters to control the mu-evolution calculation; set the radii below
 integer, private, parameter :: nmu_global  = 5       ! number of gloabl mu-values to calculate (i.e. at which radii)
 integer, private, parameter :: nmu_sink    = 6       ! number of mu-values to calculate around first two sinks (i.e. at which radii)
 !--NOT TRUE FREE PARAMETERS
 integer, private, parameter :: junit       = 47      ! unit number for *discRMnx.dat files
 integer, private, parameter :: kunit       = 62      ! unit number for *_vol*RM.dat files
 integer, private, parameter :: punit       = 23      ! unit number for rhosurf_*.dat files
 integer, private, parameter :: bunit       = 49      ! unit number for rhosurfB_*.dat file
 integer, private, parameter :: qunit       = 29      ! unit number for *in file
 integer, private, parameter :: eunit       = 39      ! unit number for global_eta* file
 integer, private            :: npt_prev    =  0      ! number of sink particles on the previous step
 integer, private            :: num_old     =  0      ! this is a magical term: do not modify
 logical, private            :: firstcall   = .true.  ! this is a magical term: do not modify
 logical, private            :: firstlog    = .true.  ! this is a magical term: do not modify
 logical, private            :: adjustcom   = .false. ! this will reset the origin to the com within Rcom of the sink
 real,    private            :: rthresh,rmerge2,h_acc2,dr,rcom2,rbin2max,rbin2maxish
 real,    private            :: rbins2(nbins)
 real,    private, allocatable :: etaart(:)
 real,    private            :: rmu_global(nmu_global),rmu_sink(nmu_sink+1)
 logical, private            :: no_file(0:maxptmass+1)
 !
 private
 !
contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,          only: mhd
 use sortutils,    only: indexx
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use centreofmass, only: reset_centreofmass
 use part,         only: igas,idust,istar,xyzmh_ptmass,vxyz_ptmass,nptmass,Bxyz
 use units,        only: udist,umass,unit_density,unit_velocity,unit_Bfield
 use physcon,      only: au,solarm
#ifdef NONIDEALMHD
 use io,           only: fatal
 use units,        only: utime
 use nicil,        only: use_ohm,use_hall,use_ambi, &
                          ion_rays,ion_thermal,use_massfrac,massfrac_X,massfrac_Y, &
                          g_cnst,fdg,rho_bulk,a0_grain,an_grain,ax_grain, &
                          zeta_cgs,mass_MionR_mp
#endif
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,iunit
 integer,          intent(inout) :: npart
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in) :: particlemass,time
 integer                      :: i,j,isink,isink0,isinkN,itype,ibin,ierr,ndens
 real                         :: dthresh,mdisc,rdisc,msink, &
                                 xi,yi,hi,rhoi,rmax,rmin,rtmp2,m_low_dens,logr_min,dlogr
 integer                      :: indx(maxp),imerged(nptmass)
 real                         :: rad2(maxp),xyz_tmp(3),vxyz_tmp(3)
 real                         :: eta(6),eta1,eta_1,eta_2,eta_3,eta_4,eta_5,eta_6,rsepmin2
 real                         :: etaohm,etahall,etaambi
 real                         :: mu_global(nmu_global),mass_mu_global(nmu_global),B_mu_global(nmu_global)
 real                         :: mu_sink(2,nmu_sink),  mass_mu_sink(2,nmu_sink),  B_mu_sink(2,nmu_sink)
 logical                      :: iexist,calc_rad_prof
 character(len=  3)           :: rval,csink
 character(len=200)           :: fileout1,fileout2,fileout3,fileout6,fileout7,fileout8,fileout9,filelog
#ifdef NONIDEALMHD
 integer                                 :: iloc
 character(len=200)                      :: infile,fileprefix
 type(inopts), dimension(:), allocatable :: db
#endif
 !
 ! Yell if contradictory commands
 if ( sinks_only .and. ignore_sinks ) then
    print*, "You cannot both ignore sinks and only centre on sinks.  Aborting task."
    return
 endif
 ! If no sinks present, abort analysis if requested
 if ( sinks_only ) then
    if ( nptmass==0 ) return
    isink0 = 1
 else
    isink0 = 0
 endif
 ! No special treatment of sinks, if requested
 if ( ignore_sinks ) then
    isinkN = 0
 else
    isinkN = nptmass
    if (sinks_two) isinkN = min(isinkN,2)
 endif
 !
 ! Initialise parameters
 !
 dthresh       = dthreshcgs/unit_density    ! gas  density threshold in code units
 rthresh       = real(rthreshAU)*au/udist
 rmerge2       = ((rmergeAU)*au/udist)**2
 m_low_dens    =   0.0
 itype         = igas
 eta           =   0.0
 imerged       =     0
 rmu_global(1) =   4.0                      ! code units = 0.013pc = 2680au = initial radius of cloud
 rmu_global(2) =   1.8                      ! code units           = 1206au ~ size of sphere at 1.51tff
 rmu_global(3) = 500.0*au/udist
 rmu_global(4) = 200.0*au/udist
 rmu_global(5) =   0.0                      ! dynamically calculated as R_sep+0.5*(r_disc1+r_disc_2)
 rmu_sink(1)   =  30.0*au/udist
 rmu_sink(2)   =  60.0*au/udist
 rmu_sink(3)   =  90.0*au/udist
 rmu_sink(4)   = 120.0*au/udist
 rmu_sink(5)   = 200.0*au/udist
 rmu_sink(6)   =   0.0                      ! dynamically calculated at R_disc
 mu_sink       =   0.0
 !
 ! If num <= num_old, then treat as firstcall again
 if (num <= num_old) firstcall = .true.
 num_old = num
 !
 ! Initalise eos, read in important nonideal MHD details from the .in file, and initialise
 if ( firstcall ) then
    firstcall = .false.
    ! Re-initialise no_file (will be unnecessary if only analysing one file simulation)
    no_file   = .true.
    ! Initialise the equation of state
    call init_eos(ieos,ierr)
    rcom2 = (rcom_au*au/udist)**2
    ! Initialise the radial bins for calculation of disc mass
    rmax = rmax_au*au/udist ! Maximum radial distance of interest (code units)
    rmin = rmin_au*au/udist ! Minimum radial distance of interest (code units)
    if (log_rbin) then
       logr_min = log10(rmin)
       dlogr    = (log10(rmax) - logr_min )/nbins
       do i = 1,nbins
          rbins2(i) = (10**(logr_min +(i-1)*dlogr ))**2
       enddo
    else
       dr   = rmax/float(nbins)
       do i = 1,nbins
          rbins2(i) = (float(i)*dr)**2
       enddo
    endif
    print*, rbins2
    rbin2max = 2.0*rbins2(nbins)
    rbin2maxish = 0.9*rbin2max
#ifdef NONIDEALMHD
    iloc = index(dumpfile,'_0')
    if (iloc > 1) then
       fileprefix = trim(dumpfile(1:iloc-1))
    else
       fileprefix = trim(dumpfile)
    endif
    infile = trim(fileprefix)//'.in'
    inquire(file=trim(infile),exist=iexist)
    if (iexist) then
       call open_db_from_file(db,infile,qunit,ierr)
       call read_inopt(use_ohm,      'use_ohm',      db,ierr)
       call read_inopt(use_hall,     'use_hall',     db,ierr)
       call read_inopt(use_ambi,     'use_ambi',     db,ierr)
       call read_inopt(ion_rays,     'ion_rays',     db,ierr)
       call read_inopt(ion_thermal,  'ion_thermal',  db,ierr)
       call read_inopt(use_massfrac, 'use_massfrac', db,ierr)
       call read_inopt(massfrac_X,   'massfrac_X',   db,ierr)
       call read_inopt(massfrac_Y,   'massfrac_Y',   db,ierr)
       call read_inopt(g_cnst,       'g_cnst',       db,ierr)
       call read_inopt(fdg,          'fdg',          db,ierr)
       call read_inopt(rho_bulk,     'rho_bulk',     db,ierr)
       call read_inopt(a0_grain,     'a0_grain',     db,ierr)
       call read_inopt(an_grain,     'an_grain',     db,ierr)
       call read_inopt(ax_grain,     'ax_grain',     db,ierr)
       call read_inopt(zeta_cgs,     'zeta_cgs',     db,ierr)
       call read_inopt(mass_MionR_mp,'mass_MionR_mp',db,ierr)
       call close_db(db)
       close(qunit)
       n_R         = 0.0
       n_electronT = 0.0
    endif
    ! Initialise nicil
    call nicil_initialise(utime,udist,umass,unit_Bfield,ierr)
    if (ierr/=0) call fatal('initial','error initialising nicil (the non-ideal MHD library)')
#else
    calc_eta = .false.
#endif
 endif
 !
 ! Calculate artificial resistivity
 !
 allocate(etaart(maxp))
 etaart = 0.
 if (use_etaart_old) then
    print*, "THIS IS NOT A TRUE REPRESENTATION OF ETA_art since it uses a different vsig!"
!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(npart,xyzh,alphaB,iphase,massoftype,etaart,Bxyz,dthresh) &
!$omp private(i,hi,itype,rhoi)
!$omp do
    do i = 1,npart
       hi = xyzh(4,i)
       if (.not.isdead_or_accreted(hi)) then
          if (maxphase==maxp) itype = iamtype(iphase(i))
          rhoi = rhoh(hi,massoftype(itype))
          if (rhoi > dthresh .and. itype==igas) then ! to save time since we never care about low density material
             etaart(i) = etaart_old(hi,rhoi,alphaB,Bxyz(1:3,i))
          endif
       endif
    enddo
!$omp enddo
!$omp end parallel
 else
    print*, "THIS IS LIKELY NOT A TRUE REPRESENTATION OF ETA-art since the algorithm is only the same in spirit!"
    print*, "starting to calculate etaart"
!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(npart,xyzh,vxyzu,iphase,massoftype,etaart,Bxyz,dthresh) &
!$omp private(i,hi,itype,rhoi)
!$omp do
    do i = 1,npart
       hi = xyzh(4,i)
       if (.not.isdead_or_accreted(hi)) then
          if (maxphase==maxp) itype = iamtype(iphase(i))
          rhoi = rhoh(hi,massoftype(itype))
          if (rhoi > dthresh .and. itype==igas) then ! to save time since we never care about low density material
             etaart(i) = etaart_new(i,npart,massoftype(itype),xyzh,vxyzu)
          endif
       endif
    enddo
!$omp enddo
!$omp end parallel
    print*, "completed calculating etaart"
 endif
 !
 ! Perform the calculations for each sink; if no sink, reset to the centre of mass
 !
 over_sinks: do isink = isink0,isinkN
    calc_rad_prof = .not. ignor_lowrho
    if (isink > 0) calc_rad_prof = .true.
    rsepmin2 = huge(rsepmin2)
    ! Skip 'merged' sinks
    if (isink > 0) then
       do j = 1,isinkN
          if (j/=isink) then
             rtmp2    = (xyzmh_ptmass(1,isink)-xyzmh_ptmass(1,j))**2 &
                      + (xyzmh_ptmass(2,isink)-xyzmh_ptmass(2,j))**2
             rsepmin2 = min(rsepmin2,rtmp2)
          endif
       enddo
       if (imerged(isink) > 0) cycle over_sinks
       ! Determine if two sinks should be merged
       do j = isink+1,isinkN
          rtmp2 = (xyzmh_ptmass(1,isink)-xyzmh_ptmass(1,j))**2 &
                + (xyzmh_ptmass(2,isink)-xyzmh_ptmass(2,j))**2 &
                + (xyzmh_ptmass(3,isink)-xyzmh_ptmass(3,j))**2
          if (rtmp2 < rmerge2) then
             write(filelog,'(3a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'.log'
             inquire(file=filelog,exist=iexist)
             if ( firstlog .or. .not.iexist ) then
                firstlog = .false.
                open(unit=9969,file=filelog,status='replace')
             else
                open(unit=9969,file=filelog,position='append')
             endif
             write(9969,'(a,f6.2,3(a,I5),a)') "At time = ",time," (num=",num,") sinks ",isink," & ",j," were treated at one"
             close(9969)
             xyzmh_ptmass(1:3,isink) = 0.5*(xyzmh_ptmass(1:3,isink)+xyzmh_ptmass(1:3,j) )
             xyzmh_ptmass(4,  isink) = xyzmh_ptmass(4, isink) + xyzmh_ptmass(4,j)
             imerged(j) = 1
          endif
       enddo
    endif
    !
    ! Open files: note, give no sink file id as if there was a sink.  Maybe we should add nptmasses to the list
    write(rval, '(I3.3)') rthreshAU
    write(csink,'(I3.3)') isink
    if (isink==0) then
       write(fileout1,'(3a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_discRM.dat'
       write(fileout2,'(3a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_discRMnx.dat'
       write(fileout3,'(5a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_vol',rval,'RM.dat'
    else
       write(fileout1,'(5a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_S',csink,'discRM.dat'
       write(fileout2,'(5a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_S',csink,'discRMnx.dat'
       write(fileout3,'(7a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_S',csink,'vol',rval,'RM.dat'
    endif
    if ( no_file(isink) ) then
       open(iunit,file=fileout1,status='replace')
       call write_header_file1(iunit)
       if (.not. discRM_only) then
          open(junit,file=fileout2,status='replace')
          call write_header_file2(junit)
       endif
       if ( printvol ) then
          open(kunit,file=fileout3,status='replace')
          call write_header_file3(kunit)
       endif
    else
       open(iunit,file=fileout1,position='append')
       if (.not. discRM_only) open(junit,file=fileout2,position='append')
       if (printvol) open(kunit,file=fileout3,position='append')
    endif
    !
    if (isink==0) then
       ! Reset the centre of mass to the origin
       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
       h_acc2 = 0.0
    else
       ! Reset the centre of mass & velocity to be on sink particle isink
       xyz_tmp  = xyzmh_ptmass(1:3,isink)
       vxyz_tmp = vxyz_ptmass (1:3,isink)
       call reset_origin(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,xyz_tmp,vxyz_tmp)
       msink  = xyzmh_ptmass(4,isink)
       h_acc2 = xyzmh_ptmass(5,isink)*xyzmh_ptmass(5,isink)
       !--Adjust the CoM, if requested
       if ( adjustcom ) then
          call adjust_origin(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,xyz_tmp,vxyz_tmp,dthresh,particlemass)
          call reset_origin(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,xyz_tmp,vxyz_tmp)
       endif
    endif
    !
    ! Initialise parameters
    itype = igas
    mdisc = 0.0
    rad2  = rbin2max
    ndens = 0
    eta_1 = 0.0
    eta_2 = 0.0
    eta_3 = 0.0
    eta_4 = 0.0
    eta_5 = 0.0
    eta_6 = 0.0
    !  give warning, if necessary
    if (isink > 0 .and. rsepmin2 < rbins2(nbins)) then
       write(*,'(a,Es12.5)')'min sep of this BH and its non-merged neighbour: ',sqrt(rsepmin2)*udist/au
    endif
    ! Determines disc mass and properties
!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(npart,isink,isink0,isinkN,xyzh,Bxyz,n_R,n_electronT,etaart,iphase) &
!$omp shared(calc_eta,particlemass,dthresh,rsepmin2,rad2,dr,calc_rad_prof,rbins2,log_rbin) &
!$omp private(i,xi,yi,hi,rhoi,rtmp2,ibin,etaohm,etahall,etaambi) &
!$omp firstprivate(itype) &
!$omp reduction(+:ndens,mdisc,m_low_dens,eta_1,eta_2,eta_3,eta_4,eta_5,eta_6)
!$omp do
    do i = 1,npart
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       hi = xyzh(4,i)
       if (.not.isdead_or_accreted(hi)) then
          if (maxphase==maxp) itype = iamtype(iphase(i))
          if (itype==igas) then
             rhoi = rhoh(hi,particlemass)
             if (rhoi > dthresh) then
                calc_rad_prof = .true.
                rtmp2  = xi*xi + yi*yi
                if (log_rbin) then
                   ibin = 1
                   do while(rbins2(ibin) < rtmp2 .and. ibin < nbins)
                      ibin = ibin + 1
                   enddo
                else
                   ibin = min(int(sqrt(rtmp2)/dr)+1,nbins)
                endif
                if (rtmp2 < rsepmin2*0.25 .and. ibin < nbins) then
                   ndens   = ndens + 1
                   mdisc   = mdisc + particlemass
                   rad2(i) = rtmp2
                endif
                if (isink==isink0 .and. calc_eta) then
                   call get_eta_global(etaohm,etahall,etaambi,rhoi,n_R(:,i),n_electronT(i),Bxyz(1:3,i))
                   eta_1 = eta_1 + etaart(i)
                   eta_2 = eta_2 + etaohm
                   eta_3 = eta_3 + etahall
                   eta_4 = eta_4 + abs(etahall)
                   eta_5 = eta_5 + etaambi
                   eta_6 = eta_6 + 1.0
                endif
             else
                m_low_dens = m_low_dens + particlemass
             endif
          else
             xyzh(4,i) = -abs(xyzh(4,i)) ! just to kill all non-gas particle for the duration of the analysis
          endif
       endif
    enddo
!$omp enddo
!$omp end parallel
    eta(1) = eta_1
    eta(2) = eta_2
    eta(3) = eta_3
    eta(4) = eta_4
    eta(5) = eta_5
    eta(6) = eta_6
    !
    ! Sort the rad2 array (this occurs in increasing order of distance)
    call indexx(npart,rad2,indx)
    !
    ! Determine the radii which contains massfrac of the mass (total & component-wise)
    if (isink > 0) then
       call get_mass_and_radius(npart,ndens,rad2,xyzh(3,:),xyzh(4,:),indx,particlemass,mdisc,rdisc)
    else
       call get_radius(npart,rdisc,msink,(msink+mdisc)*massfrac,rad2,massoftype(igas),indx)
    endif
    !
    ! Call analysis to get the (r,phi,z) components of the B & V fields;  this is for gas only!
    if ( (.not. discRM_only .and. calc_rad_prof) .or. angle > 0.0 ) then
       call doanalysisRPZ(csink,dumpfile,num,npart,xyzh,vxyzu,Bxyz,particlemass,dthresh &
                         ,au,udist,umass,solarm,unit_velocity,unit_Bfield,rdisc**2,time,mhd,rthresh**2,msink)
    endif
    !
    ! Calculate the evolution of mu
    if (calc_mu) then
       if (isink==1 .or. isink==2) then
          if (isink==1) rmu_sink(nmu_sink+1) = rdisc
          rmu_sink(nmu_sink) = rdisc
          call get_mu(npart,nptmass,nmu_sink,rmu_sink(1:nmu_sink),mu_sink(isink,:),mass_mu_sink(isink,:), &
                         B_mu_sink(isink,:),xyzh,xyzmh_ptmass,Bxyz,particlemass)
       endif
    endif
    !
    ! Convert mass & radius to M_sun & AU, respectively & Print results to file
    m_low_dens = m_low_dens*umass/solarm
    mdisc      = mdisc*umass/solarm
    msink      = msink*umass/solarm
    rdisc      = rdisc*udist/au
    write(iunit, '(I18,1x,7(1pe18.10,1x))') num,time,m_low_dens,msink,mdisc,mdisc+msink,rdisc
    !
    close(iunit)
    close(junit)
    close(kunit)
 enddo over_sinks
 !
 ! Print the globally averaged eta-values, for particles with rho > rho_crit
 if (calc_eta) then
    write(fileout6,'(3a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_eta.dat'
    if ( no_file(maxptmass+1) ) then
       open(eunit,file=fileout6,status='replace')
       call write_header_file6(eunit)
    else
       open(eunit,file=fileout6,position='append')
    endif
    if (eta(6) > 0.0) then
       eta1 = 1.0/eta(6)
    else
       eta1 = 0.0
    endif
    write(eunit,'(I18,1x,6(1pe18.10,1x),2I18)') num,time,eta(1:5)*eta1*unit_eta,nptmass,npt_prev
    npt_prev = nptmass
    close(eunit)
 endif
 !
 ! Calculate the global evolution and print the evolution of mu to file
 if (calc_mu) then
    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
    if (nptmass > 1) then
       rmu_global(nmu_global) = sqrt( (xyzmh_ptmass(1,1) - xyzmh_ptmass(1,2))**2   &
                              +       (xyzmh_ptmass(2,1) - xyzmh_ptmass(2,2))**2   &
                              +       (xyzmh_ptmass(3,1) - xyzmh_ptmass(3,2))**2 ) &
                              + 0.5*(rmu_sink(nmu_sink)+rmu_sink(nmu_sink+1))
    else
       rmu_global(nmu_global) = 0.
    endif
    call get_mu(npart,nptmass,nmu_global,rmu_global,mu_global,mass_mu_global,B_mu_global, &
                xyzh,xyzmh_ptmass,Bxyz,particlemass)
    !
    write(fileout7,'(3a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_mu.dat'
    write(fileout8,'(3a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_mu_mass.dat'
    write(fileout9,'(3a)') 'analysisout_',trim(dumpfile(1:INDEX(dumpfile,'_')-1)),'_mu_B.dat'
    if ( no_file(maxptmass+1) ) then
       open(eunit,file=fileout7,status='replace')
       call write_header_file7(eunit)
       if (calc_mu_verb) then
          open(eunit+1,file=fileout8,status='replace')
          open(eunit+2,file=fileout9,status='replace')
          call write_header_file7(eunit+1)
          call write_header_file7(eunit+2)
       endif
    else
       open(eunit,file=fileout7,position='append')
       if (calc_mu_verb) then
          open(eunit+1,file=fileout8,position='append')
          open(eunit+2,file=fileout9,position='append')
       endif
    endif
    write(eunit,'(I18,1x,21(1pe18.10,1x))') num,time,mu_global,rmu_global(nmu_global)/au*udist,&
       mu_sink(1,:),rmu_sink(nmu_sink+1)/au*udist,mu_sink(2,:),rmu_sink(nmu_sink)/au*udist
    close(eunit)
    if (calc_mu_verb) then
       write(eunit+1,'(I18,1x,21(1pe18.10,1x))') num,time,mass_mu_global,rmu_global(nmu_global)/au*udist,&
          mass_mu_sink(1,:),rmu_sink(nmu_sink+1)/au*udist,mass_mu_sink(2,:),rmu_sink(nmu_sink)/au*udist
       write(eunit+2,'(I18,1x,21(1pe18.10,1x))') num,time,B_mu_global,rmu_global(nmu_global)/au*udist,&
          B_mu_sink(1,:),rmu_sink(nmu_sink+1)/au*udist,B_mu_sink(2,:),rmu_sink(nmu_sink)/au*udist
       close(eunit+1)
       close(eunit+2)
    endif
 endif
 do isink=isink0,isinkN
    no_file(isink) = .false.
 enddo
 no_file(maxptmass+1) = .false.
!
end subroutine do_analysis
!
!----------------------------------------------------------------
!+
!  Determine the radius of the disc, assuming a single disc
!+
!----------------------------------------------------------------
subroutine get_radius(npart,rdisc,msink,discmasslim,rad2,pmass,indx)
 integer, intent(in)  :: npart,indx(:)
 real,    intent(out) :: rdisc
 real,    intent(in)  :: msink,discmasslim,rad2(:),pmass
 integer              :: i,j
 real                 :: discmass
 !
 i        = 1
 discmass = msink
 do while (i <= npart .and. discmass < discmasslim)
    j = indx(i)
    discmass = discmass + pmass
    i = i + 1
 enddo
 if (i==1) then
    rdisc = 0.0
 else if (i==npart+1) then
    rdisc = huge(rdisc)
 else
    rdisc = sqrt(rad2(j))
 endif
 !
end subroutine get_radius
!
!----------------------------------------------------------------
!+
!  Determine the mass and radius of the disc, assuming a binary
!  system disc
!  The edge of the disc is defined as where the density falls to
!  1/10th the max density, as calculated in rings
!  Other tested criteria include forcing mdisc < 0.99mdisc_max
!  (no change), and to force each ring to have rho > rho_crit, which
!  reduces the radius slightly, most noticeably when they are already
!  large.  Thus, 1/10th (where 1/10 is a free parameter) is our
!  only criteria
!+
!----------------------------------------------------------------
subroutine get_mass_and_radius(npart,ndens,rad2,zdir,hdir,indx,pmass,mdisc,rdisc)
 real,    intent(in)    :: rad2(:),zdir(:),hdir(:),pmass
 integer, intent(in)    :: npart,ndens,indx(:)
 real,    intent(out)   :: rdisc,mdisc
 integer                :: i,j,k,ninring_loc
 real                   :: zmin,zmax,totmass,notdisc,dring,rmin2,rmax2,dmax,dr0,rad2old,rad2now,rad2next
 logical                :: indisc
 !
 ! Initialise parameters, including several to simply avoid compiler warnings
 indisc  = .true.
 rmin2   = 0.0
 rmax2   = 0.0
 mdisc   = 0.0
 notdisc = 0.0
 totmass = 0.0
 dring   = 0.0
 dmax    = 0.0
 rad2old = 0.0
 j       = 0
 ! Modify the number of particles in ring if there are too few high density particles
 if (ndens > 4*ninring) then
    ninring_loc = ninring
 else if (ndens > ninring) then
    ninring_loc = int(0.25*ninring)
 else if (ndens > 0.5*ninring) then
    ninring_loc = int(0.1*ninring)
 else
    write(*,'(a,I6,a)') "There are ",ndens," particles above the density threhhold.  Assuming no disc"
    mdisc = 0.0
    rdisc = 0.0
    return
 endif
 !
 do while (j < npart)
    j       = j + 1
    i       = indx(j)
    k       = 1
    zmin    = zdir(i)
    zmax    = zdir(i)
    totmass = pmass
    do while (k <= ninring_loc .and. j < npart)
       j       = j + 1
       i       = indx(j)
       rad2now = rad2(i)
       if (rad2now > rbin2maxish) then
          k = ninring_loc ! to exit loop since we are now into dummy radii
       else
          k       = k + 1
          dr0      = sqrt(rad2now) - sqrt(rad2old)
          if (dr0 < 10.0*hdir(i)) then
             totmass  = totmass + pmass
             zmin     = min(zmin,zdir(i))
             zmax     = max(zmax,zdir(i))
             rad2next = rad2(indx(j+1))
             if (rad2next > rad2now .and. rad2next < rbin2maxish ) then
                rmax2 = ( 0.5*(sqrt(rad2now) + sqrt(rad2next)) )**2
             else
                rmax2 = rad2now
             endif
          else
             k = ninring_loc ! to exit loop
          endif
       endif
       rad2old = rad2now
    enddo
    dring = totmass/(pi*(rmax2-rmin2)*(zmax-zmin))
    dmax  = max(dmax,dring)
    if (dring > dringthresh*dmax) then
       ! this shell is in the disc, so temporarally record it as the outer values
       rdisc = sqrt(rmax2)
       mdisc = mdisc + totmass
       if (.not. indisc) then
          ! there was a small gap, but that material should be included in the total mass
          mdisc   = mdisc + notdisc
          notdisc = 0.0
          indisc  = .true.
       endif
    else
       notdisc = notdisc + totmass
       indisc  = .false.
    endif
    rmin2 = rmax2
 enddo
 !
end subroutine get_mass_and_radius
!----------------------------------------------------------------
!+
!  Calculate the global eta values for rho > rho_crit
!+
!----------------------------------------------------------------
subroutine get_eta_global(etaohm,etahall,etaambi,rhoi,n_Ri,n_electronTi,B)
 real,   intent(inout) :: n_Ri(:),n_electronTi
 real,   intent(in)    :: rhoi
 real,   intent(in)    :: B(3)
 real,   intent(out)   :: etaohm,etahall,etaambi
 integer               :: ierr
 real                  :: temperature,B2i,vdummy(maxvxyzu),xdummy(3)
 !
 ! Calculate temperature
 xdummy      = 0.0
 vdummy      = 0.0
 temperature = get_temperature(ieos,xdummy,rhoi,vdummy)
 !
 ! Calculate the artificial coefficient
 B2i         = sqrt( dot_product(B,B) )
 ! Calculate physical coefficients
 call nicil_get_ion_n(real(rhoi),temperature,n_Ri(:),n_electronTi,ierr)
 call nicil_get_eta(etaohm,etahall,etaambi,sqrt(B2i),real(rhoi),temperature, &
                           n_Ri(:),n_electronTi,ierr)
 !
end subroutine get_eta_global
!----------------------------------------------------------------
!+
!  Calculate artificial resistivity assuming old resistivity
!+
!----------------------------------------------------------------
real function etaart_old(hi,rhoi,alphaB_in,B)
 real, intent(in)    :: hi,rhoi,alphaB_in
 real, intent(in)    :: B(3)
 real                :: spsoundi,B2i,valfven2i,vsigi,vdummy(maxvxyzu),xdummy(3)
 !
 xdummy   = 0.0
 vdummy   = 0.0
 spsoundi = get_spsound(ieos,xdummy,rhoi,vdummy)
 !
 ! Calculate the artificial coefficient
 B2i       = sqrt( dot_product(B,B) )
 valfven2i = B2i/rhoi
 vsigi     = sqrt(valfven2i + spsoundi*spsoundi)
 etaart_old = 0.5*hi*vsigi*alphaB_in
 !
end function etaart_old
!----------------------------------------------------------------
!+
!  Calculate artificial resistivity assuming new resistivity
!  This method defines alphaB == 1
!  NOTE! This method is still likely incorrect compared to Phantom
!        since this used v and dr, whereas Phantom more accurately
!        uses dv and dr, and has merged the div.V into this section
!+
!----------------------------------------------------------------
real function etaart_new(ipart,npart,pmass,xyzh,vxyzu)
 use kernel, only: get_kernel,radkern2,cnormk
 integer,      intent(in)    :: ipart,npart
 real,         intent(in)    :: pmass,xyzh(:,:),vxyzu(:,:)
 integer                     :: j
 real                        :: xi,yi,zi,hi,hi2,vxi,vyi,vzi,xj,yj,zj,hj,hj2,dx2,dy2,dz2,rad2
 real                        :: runix,runiy,runiz,rhoj
 real                        :: qi,q2i,radkern2i,wkern,grkern
 real                        :: vsigx,vsigy,vsigz,vsigB
 !
 xi  = xyzh(1,ipart)
 yi  = xyzh(2,ipart)
 zi  = xyzh(3,ipart)
 hi  = xyzh(4,ipart)
 hi2 = hi*hi
 vxi = vxyzu(1,ipart)
 vyi = vxyzu(2,ipart)
 vzi = vxyzu(3,ipart)
 etaart_new = 0.0
 do j = 1,npart
    xj  = xyzh(1,j)
    yj  = xyzh(2,j)
    zj  = xyzh(3,j)
    hj  = xyzh(4,j)
    hj2 = hj*hj
    radkern2i = radkern2*hi2
    if (.not.isdead_or_accreted(hj) .and. j/=ipart) then
       dx2 = (xi-xj)**2
       if (dx2 < radkern2i) then
          dy2 = (yi-yj)**2
          if (dy2 < radkern2i) then
             dz2 = (zi-zj)**2
             if (dz2 < radkern2i) then
                rad2 = dx2 + dy2 + dz2
                if (rad2 < radkern2i) then
                   rhoj       = rhoh(hj,pmass)
                   runix      = (xi-xj)/sqrt(rad2)
                   runiy      = (yi-yj)/sqrt(rad2)
                   runiz      = (zi-zj)/sqrt(rad2)
                   vsigx      = vyi*runiz - vzi*runiy
                   vsigy      = vzi*runix - vxi*runiz
                   vsigz      = vxi*runiy - vyi*runix
                   vsigB      = sqrt(vsigx*vsigx + vsigy*vsigy + vsigz*vsigz)
                   q2i        = rad2/hi2
                   qi         = sqrt(q2i)
                   call get_kernel(q2i,qi,wkern,grkern)
                   etaart_new = etaart_new + vsigB*wkern/rhoj
                endif
             endif
          endif
       endif
    endif
 enddo
 etaart_new = 0.5*etaart_new*hi* (pmass*cnormk/hi**3) ! the terms in brackets is from the kernel weighting
 !
end function etaart_new
!----------------------------------------------------------------
!+
!  Calculates the evolution of the mass-to-flux ratio
!+
!----------------------------------------------------------------
subroutine get_mu(npart,nptmass,nrad,rad_mu,mu,mass,B,xyzh,xyzmh_ptmass,Bxyz,pmassi)
 integer, intent(in)  :: npart,nptmass,nrad
 real,    intent(in)  :: pmassi,rad_mu(:),xyzh(:,:),xyzmh_ptmass(:,:)
 real,    intent(in)  :: Bxyz(:,:)
 real,    intent(out) :: mu(nrad),mass(nrad),B(nrad)
 integer              :: i,j
 real                 :: rmasstoflux_crit
 real                 :: xi,yi,zi,hi,hi3,mi,Bxi,Byi,Bzi,Bi,rad2,rad2_mu(nrad)
 real                 :: mass_thread(nrad),vol(nrad),vol_thread(nrad),B_thread(nrad)
 !
 rad2_mu = rad_mu*rad_mu
 mass    = 0.0
 vol     = 0.0
 B       = 0.0
 mu      = 0.0
!$omp parallel default(none) &
!$omp shared(npart,nrad,xyzh,Bxyz,iphase,pmassi,rad2_mu) &
!$omp shared(mass,vol,B) &
!$omp private(i,j,xi,yi,zi,hi,hi3,Bxi,Byi,Bzi,Bi,rad2) &
!$omp private(mass_thread,vol_thread,B_thread)
 mass_thread = 0.0
 vol_thread  = 0.0
 B_thread    = 0.0
!$omp do
 do i = 1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (hi > 0.0) then
       hi3 = hi*hi*hi
       Bxi = Bxyz(1,i)
       Byi = Bxyz(2,i)
       Bzi = Bxyz(3,i)
       rad2 = xi*xi + yi*yi + zi*zi  ! Note that we are already centred on the point of interest
       Bi   = sqrt(Bxi*Bxi + Byi*Byi + Bzi*Bzi)
       do j = 1,nrad
          if (rad2 < rad2_mu(j)) then
             mass_thread(j) = mass_thread(j) + pmassi
             vol_thread(j)  = vol_thread(j)  + hi3
             B_thread(j)    = B_thread(j)    + hi3*Bi
          endif
       enddo
    endif
 enddo
!$omp enddo
!$omp critical(collatedata)
 do j = 1,nrad
    mass(j) = mass(j) + mass_thread(j)
    vol(j)  = vol(j)  + vol_thread(j)
    B(j)    = B(j)    + B_thread(j)
 enddo
!$omp end critical(collatedata)
!$omp end parallel
 !
 ! Add mass contribution of point masses
 do i = 1,nptmass
    xi   = xyzmh_ptmass(1,i)
    yi   = xyzmh_ptmass(2,i)
    zi   = xyzmh_ptmass(3,i)
    mi   = xyzmh_ptmass(4,i)
    rad2 = xi*xi + yi*yi + zi*zi  ! Note that we are already centred on the point of interest
    do j = 1,nrad
       if (rad2 < rad2_mu(j)) mass(j) = mass(j) + mi
    enddo
 enddo
 !
 ! Calculate the mass to flux ratio
 do j = 1,nrad
    if (vol(j) > 0.0) B(j) = B(j)/vol(j)
    if (rad2_mu(j) > 0.0 .and. B(j) > 0.0) mu(j) = mass(j)/(pi*rad2_mu(j)*B(j))
 enddo
 rmasstoflux_crit = 2./3.*0.53*sqrt(5./pi)
 mu               = mu/rmasstoflux_crit
 !
end subroutine get_mu
!----------------------------------------------------------------
!+
!  Calculate radial profiles at the given time
!+
!----------------------------------------------------------------
subroutine doanalysisRPZ(csink,dumpfile,num,npart,xyzh,vxyzu,Bxyz,particlemass,dthreshg &
                        ,au,udist,umass,solarm,unit_velocity,unit_Bfield &
                        ,rdisc2,time,mhd,rthresh2,dmassp)
 character(len=*), intent(in)    :: dumpfile,csink
 integer,          intent(in)    :: npart,num
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(inout) :: vxyzu(:,:)
 real,             intent(in)    :: rdisc2,time,rthresh2,dmassp
 real,             intent(in)    :: Bxyz(:,:)
 real,             intent(in)    :: particlemass,dthreshg,au,udist,umass,solarm,unit_velocity,unit_Bfield
 logical,          intent(in)    :: mhd
 !
 ! Indicies of the disc bins (Dbins): TIME AVERAGED
 !  Note: Cbins is similar, but different radial cutoff; currently left with hardcoded indicies
 integer, parameter           :: iDvr      =  1 ! vr
 integer, parameter           :: iDvphi    =  2 ! vphi
 integer, parameter           :: iDvz      =  3 ! vz
 integer, parameter           :: iDbr      =  4 ! Br
 integer, parameter           :: iDbphi    =  5 ! Bphi
 integer, parameter           :: iDbp      =  6 ! Bp = sqrt(Bz**2+Br**2)
 integer, parameter           :: iDbz      =  7 ! Bz
 integer, parameter           :: iDby      =  8 ! By
 integer, parameter           :: iDbx      =  9 ! Bx
 integer, parameter           :: iDb       = 10 ! B
 integer, parameter           :: iDbrat    = 11 ! Bphi/Bp
 integer, parameter           :: iDbeta    = 12 ! plasma beta
 integer, parameter           :: iDLx      = 13 ! x-angular momentum
 integer, parameter           :: iDLy      = 14 ! y-angular momentum
 integer, parameter           :: iDLz      = 15 ! z-angular momentum
 integer, parameter           :: iDetaF    = 16 ! eta_art
 integer, parameter           :: iDetaO    = 17 ! eta_ohm
 integer, parameter           :: iDetaH    = 18 ! eta_hall
 integer, parameter           :: iDetaHb   = 19 ! abs(eta_hall)
 integer, parameter           :: iDetaA    = 20 ! eta_ambi
 integer, parameter           :: iDetaHp   = 21 ! eta_hall > 0
 integer, parameter           :: iDetaHn   = 22 ! eta_hall < 0
 integer, parameter           :: iDetaOA   = 23 ! eta_ohm/eta_art
 integer, parameter           :: iDetaHA   = 24 ! eta_hall/eta_art
 integer, parameter           :: iDetaAA   = 25 ! eta_ambi/eta_art
 integer, parameter           :: iDtemA    = 26 ! temperature (mean)
 integer, parameter           :: iDtemX    = 27 ! temperature (maximum)
 integer, parameter           :: iDnn      = 28 ! n_neutral
 integer, parameter           :: iDngm     = 29 ! n_grain(Z=-1)
 integer, parameter           :: iDngz     = 30 ! n_grain(Z= 0)
 integer, parameter           :: iDngp     = 31 ! n_grain(Z=+1)
 integer, parameter           :: iDnhion   = 32 ! n_{hydrogen-helium ion species}
 integer, parameter           :: iDnmion   = 33 ! n_{metallic ion species}
 integer, parameter           :: iD        = 33 ! the number of bins in this array
 ! Indicies of values within the disc and with rho > rhocrit: RADIAL PROFILE
 integer, parameter           :: iHn       =  1 ! total number of particles
 integer, parameter           :: iHv       =  2 ! velocity
 integer, parameter           :: iHvr      =  3 ! radial velocity
 integer, parameter           :: iHvphi    =  4 ! angular velocity
 integer, parameter           :: iHb       =  5 ! magnetic field
 integer, parameter           :: iHbr      =  6 ! radial magnetic field
 integer, parameter           :: iHbphi    =  7 ! angular magnetic field
 integer, parameter           :: iHbz      =  8 ! Bz
 integer, parameter           :: iHby      =  9 ! By
 integer, parameter           :: iHbx      = 10 ! Bx
 integer, parameter           :: iHbeta    = 11 ! plasma beta
 integer, parameter           :: iHeart    = 12 ! artificial resistivity
 integer, parameter           :: iHeohm    = 13 ! ohmic resistivity
 integer, parameter           :: iHehall   = 14 ! Hall effect
 integer, parameter           :: iHehalla  = 15 ! abs(Hall effect)
 integer, parameter           :: iHeambi   = 16 ! Ambipolar diffusion
 integer, parameter           :: iHehallp  = 17 ! Hall effect > 0
 integer, parameter           :: iHehalln  = 18 ! Hall effect < 0
 integer, parameter           :: iHnn      = 19 ! n_neutral
 integer, parameter           :: iHngm     = 20 ! n_grain(Z=-1)
 integer, parameter           :: iHngz     = 21 ! n_grain(Z= 0)
 integer, parameter           :: iHngp     = 22 ! n_grain(Z=+1)
 integer, parameter           :: iHnhion   = 23 ! n_{hydrogen-helium ion species}
 integer, parameter           :: iHnmion   = 24 ! n_{metallic ion species}
 integer, parameter           :: iH        = 24 ! The number of bins in this array
 ! Indicies of the Volume arrays
 integer, parameter           :: iVN       =  1 ! total number of particles
 integer, parameter           :: iVNvphip  =  2 ! number of particles with v_phi > 0
 integer, parameter           :: iVNvphin  =  3 ! number of particles with v_phi < 0
 integer, parameter           :: iVNrhoh   =  4 ! number of high density particles
 integer, parameter           :: iVNrhol   =  5 ! number of low density particles
 integer, parameter           :: iVb       =  6 ! magnitude of magnetic fiels
 integer, parameter           :: iVbeta    =  7 ! plasma beta
 integer, parameter           :: iVLx      =  8 ! x-angular momentum
 integer, parameter           :: iVLy      =  9 ! y-angular momentum
 integer, parameter           :: iVLz      = 10 ! z-angular momentum
 integer, parameter           :: iVvphi    = 11 ! vphi
 integer, parameter           :: iVvphip   = 12 ! vphi > 0
 integer, parameter           :: iVvphin   = 13 ! vphi < 0
 ! Indicies of the rotational fractions
 integer, parameter           :: iFN       =  1 ! total number of particles
 integer, parameter           :: iFNrhoh   =  2 ! number of high density particles
 integer, parameter           :: iFNrhol   =  3 ! number of low density particles
 integer, parameter           :: iFNvphip  =  4 ! number of particles with v_phi > 0              (disc only)
 integer, parameter           :: iFNvphiph =  5 ! number of high density particles with v_phi > 0 (disc only)
 integer, parameter           :: iFNvphipl =  6 ! number of low density particles with v_phi > 0  (disc only)
 integer, parameter           :: iFNvpvrn  =  4 ! min(vphi/vr) (volume only)
 integer, parameter           :: iFNvpvra  =  5 ! ave(vphi/vr) (volume only)
 integer, parameter           :: iFNvpvrx  =  6 ! max(vphi/vr) (volume only)
 ! Local variables
 integer                      :: i,j,k,p,nobj,ierr
 integer                      :: volN(5),ibins(4,nbins)
 real                         :: hi,rhoi,rtmp2,rtmp3d2,B2i,Bi,anglei
 real                         :: vr,vphi,Br,Bphi,plasmab,temperature,rr1
 real                         :: spsoundi,ponrhoi,angx,angy,angz,ang,massint
 real                         :: xi,yi,zi,vxi,vyi,vzi,Bxi,Byi,Bzi
 real                         :: Bi1,rho1i,etaohm,etahall,etaambi,etaart1,volL
 real                         :: volP(13),Dbins(iD,nbins),Cbins(7,nbins),Hbins(2,iH)
 real                         :: fracrotVol(6),vrat,fracrotDisc(6,nbins),ReyK(7),Rey(7)
 real                         :: data_out(n_data_out)
 character(len=200)           :: fileout4,fileout5
 !
 ! Initialise variables
 nobj        = 0
 ibins       = 0
 volN        = 0
 Dbins       = 0.0
 Cbins       = 0.0
 Hbins       = 0.0
 volP        = 0.0
 fracrotVol  = 0.0
 fracrotDisc = 0.0
 etaohm      = 0.0
 etahall     = 0.0
 etaambi     = 0.0
 data_out    = 0.0
 !
 !--Bin the data
 parts: do i = 1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (hi < tiny (hi)) cycle parts      ! dead (or gas) particle
    rtmp2   = xi*xi + yi*yi
    rtmp3d2 = rtmp2 + zi*zi

    anglei = 0.0
    if (angle > 0.0) then
       anglei = abs(180./pi*asin(zi/sqrt(rtmp3d2)))
       if (anglei > angle) anglei = 0.0 ! i.e. define as outside of our range
    endif
    !--Ensure particle is in at least one radius of interest
    if (rtmp3d2 < rthresh2 .or. rtmp2 < rdisc2 .or. map_all_R .or. anglei > 0.0) then
       !--Properties of the particle
       vxi = vxyzu(1,i)
       vzi = vxyzu(2,i)
       vyi = vxyzu(3,i)
       if (mhd) then
          Bxi = Bxyz(1,i)
          Byi = Bxyz(2,i)
          Bzi = Bxyz(3,i)
       else
          Bxi = 0.0
          Byi = 0.0
          Bzi = 0.0
       endif
       if (rtmp2 > 0.0) then
          rr1  = 1.0/sqrt(rtmp2)
       else
          rr1  = 0.0
       endif
       if (log_rbin) then
          j = 1
          do while (rbins2(j) < rtmp2 .and. j < nbins)
             j = j + 1
          enddo
       else
          j = int(sqrt(rtmp2)/dr)+1
       endif
       B2i  = Bxi*Bxi + Byi*Byi + Bzi*Bzi
       Bi   = sqrt(B2i)
       vr   = ( xi*vxi + yi*vyi )*rr1
       vphi = ( xi*vyi - yi*vxi )*rr1
       angx =   yi*vzi - zi*vyi
       angy =   zi*vxi - xi*vzi
       angz =   xi*vyi - yi*vxi
       rhoi = rhoh(hi,particlemass)
       if (maxvxyzu >= 4) then
          call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,vxyzu(4,i))
       else
          call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi)
       endif
       if (B2i > 0.0) then
          plasmab  = 2.0*ponrhoi*rhoi/B2i
          Bi1      = 1.0/Bi
       else
          plasmab  = 0.0
          Bi1      = 0.0
       endif
       rho1i     = 1./rhoi
       temperature = get_temperature_from_ponrho(ponrhoi)
       if ( mhd_nonideal .and. mhd ) then
          call nicil_get_ion_n(real(rhoi),temperature,n_R(:,i),n_electronT(i),ierr)
          call nicil_get_eta(etaohm,etahall,etaambi,sqrt(B2i),rhoi,temperature, &
                           n_R(:,i),n_electronT(i),ierr,data_out)
       endif
       if (etaart(i) > 0.0) then
          etaart1 = 1.0/etaart(i)
       else
          etaart1 = 0.0
       endif
       !
       if (rtmp3d2 < rthresh2) then
          volN(iVN)      = volN(iVN)    + 1
          volP(iVb)      = volP(iVb)    + Bi
          volP(iVbeta)   = volP(iVbeta) + plasmab
          volP(iVLx)     = volP(iVLx)   + angx
          volP(iVLy)     = volP(iVLy)   + angy
          volP(iVLz)     = volP(iVLz)   + angz
          volP(iVvphi)   = volP(iVvphi) + vphi
          if (vphi > 0.0) then
             volN(iVNvphip) = volN(iVNvphip) + 1
             volP(iVvphip)  = volP(iVvphip)  + vphi
          else
             volN(iVNvphin) = volN(iVNvphin) + 1
             volP(iVvphin)  = volP(iVvphin)  + vphi
          endif
          if (j < nbins) then
             ibins(2,j) = ibins(2,j) + 1
             Cbins(2,j) = Cbins(2,j) + Bi
             Cbins(3,j) = Cbins(3,j) + plasmab
             Cbins(4,j) = Cbins(4,j) + angx
             Cbins(5,j) = Cbins(5,j) + angy
             Cbins(6,j) = Cbins(6,j) + angz
             fracrotDisc(iFN,j) = fracrotDisc(iFN,j) + 1.0
             if (rhoi > dthreshg) then
                fracrotDisc(iFNrhoh,j) = fracrotDisc(iFNrhoh,j) + 1.0
             else
                fracrotDisc(iFNrhol,j) = fracrotDisc(iFNrhol,j) + 1.0
             endif
             if (vphi > 0.0) then
                fracrotDisc(iFNvphip,j) = fracrotDisc(iFNvphip,j) + 1.0
                if (rhoi > dthreshg) then
                   fracrotDisc(iFNvphiph,j) = fracrotDisc(iFNvphiph,j) + 1.0
                else
                   fracrotDisc(iFNvphipl,j) = fracrotDisc(iFNvphipl,j) + 1.0
                endif
             endif
          endif
          if (rhoi > dthreshg) then
             volN(iVNrhoh) = volN(iVNrhoh) + 1
          else
             volN(iVNrhol) = volN(iVNrhol) + 1
          endif
          if (vphi > 0.0) then
             fracrotVol(iFN) = fracrotVol(iFN) + 1
             if (rhoi > dthreshg) then
                fracrotVol(iFNrhoh) = fracrotVol(iFNrhoh) + 1
             else
                fracrotVol(iFNrhol) = fracrotVol(iFNrhol) + 1
             endif
          endif
          vrat       = abs(vphi/vr)
          if (vrat > 0.0) then
             fracrotVol(iFNvpvrn) = min(fracrotVol(iFNvpvrn),vrat)
             fracrotVol(iFNvpvra) =     fracrotVol(iFNvpvra)+vrat
             fracrotVol(iFNvpvrx) = max(fracrotVol(iFNvpvrx),vrat)
          endif
       endif
       !
       !--If the particle is in the disc (or at least has a high enough density this could be)
       if ( (rhoi > dthreshg .and. (rtmp2 < rdisc2 .or. map_all_R)) .or. anglei > 0.0 ) then
          if (j < nbins) then
             ibins(1,j) = ibins(1,j) + 1
             Br         = ( xi*Bxi + yi*Byi )*rr1
             Bphi       = ( xi*Byi - yi*Bxi )*rr1
             !
             Dbins(iDvr,   j) = Dbins(iDvr,   j) + vr
             Dbins(iDvphi, j) = Dbins(iDvphi, j) + vphi
             Dbins(iDvz,   j) = Dbins(iDvz,   j) + vzi
             Dbins(iDbr,   j) = Dbins(iDbr,   j) + Br
             Dbins(iDbphi, j) = Dbins(iDbphi, j) + Bphi
             Dbins(iDbp  , j) = Dbins(iDbp,   j) + sqrt(Br**2 + Bzi**2)
             Dbins(iDbz,   j) = Dbins(iDbz,   j) + Bzi
             Dbins(iDby,   j) = Dbins(iDby,   j) + Byi
             Dbins(iDbx,   j) = Dbins(iDbx,   j) + Bxi
             Dbins(iDb,    j) = Dbins(iDb,    j) + Bi
             Dbins(iDbrat, j) = Dbins(iDbrat, j) + Bphi/sqrt(Br**2 + Bzi**2)
             Dbins(iDbeta, j) = Dbins(iDbeta, j) + plasmab
             Dbins(iDLx,   j) = Dbins(iDLx,   j) + angx
             Dbins(iDLy,   j) = Dbins(iDLy,   j) + angy
             Dbins(iDLz,   j) = Dbins(iDLz,   j) + angz
             Dbins(iDetaF, j) = Dbins(iDetaF, j) + etaart(i)
             Dbins(iDetaO, j) = Dbins(iDetaO, j) + etaohm
             Dbins(iDetaH, j) = Dbins(iDetaH, j) + etahall
             Dbins(iDetaHb,j) = Dbins(iDetaHb,j) + abs(etahall)
             Dbins(iDetaA, j) = Dbins(iDetaA, j) + etaambi
             Dbins(iDetaOA,j) = Dbins(iDetaOA,j) + etaohm*etaart1
             Dbins(iDetaHA,j) = Dbins(iDetaHA,j) + etahall*etaart1
             Dbins(iDetaAA,j) = Dbins(iDetaAA,j) + etaambi*etaart1
             Dbins(iDtemA, j) = Dbins(iDtemA, j) + temperature
             Dbins(iDtemX, j) = max(Dbins(iDtemX,j), temperature)
             Dbins(iDnn,   j) = Dbins(iDnn,   j) + data_out( 7)
             Dbins(iDngm,  j) = Dbins(iDngm,  j) + data_out(12)
             Dbins(iDngz,  j) = Dbins(iDngz,  j) + data_out(13)
             Dbins(iDngp,  j) = Dbins(iDngp,  j) + data_out(14)
             Dbins(iDnhion,j) = Dbins(iDnhion,j) + data_out( 8)
             Dbins(iDnmion,j) = Dbins(iDnmion,j) + data_out( 9)
             if (etahall > 0.0) then
                Dbins(iDetaHp,j) = Dbins(iDetaHp,j) + etahall
                ibins(3,      j) = ibins(3,      j) + 1
             else if (etahall < 0.0) then
                Dbins(iDetaHn,j) = Dbins(iDetaHn,j) + etahall
                ibins(4,      j) = ibins(4,      j) + 1
             endif
             if (rtmp2 < rdisc2) then
                p = 2                   ! care only about gas in the disc
             else
                p = 1                   ! care about all gas with rho > rhocrit
             endif
             do k = 1,p
                Hbins(k,iHn)      = Hbins(k,iHn)      + 1.0
                Hbins(k,iHv)      = Hbins(k,iHv)      + sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
                Hbins(k,iHvr)     = Hbins(k,iHvr)     + vr
                Hbins(k,iHvphi)   = Hbins(k,iHvphi)   + vphi
                Hbins(k,iHb)      = Hbins(k,iHb)      + Bi
                Hbins(k,iHbr)     = Hbins(k,iHbr)     + Br
                Hbins(k,iHbphi)   = Hbins(k,iHbphi)   + Bphi
                Hbins(k,iHbz)     = Hbins(k,iHbz)     + Bzi
                Hbins(k,iHby)     = Hbins(k,iHby)     + Byi
                Hbins(k,iHbx)     = Hbins(k,iHbx)     + Bxi
                Hbins(k,iHbeta)   = Hbins(k,iHbeta)   + plasmab
                Hbins(k,iHeart)   = Hbins(k,iHeart)   + etaart(i)
                Hbins(k,iHeohm)   = Hbins(k,iHeohm)   + etaohm
                Hbins(k,iHehall)  = Hbins(k,iHehall)  + etahall
                Hbins(k,iHehalla) = Hbins(k,iHehalla) + abs(etahall)
                Hbins(k,iHeambi)  = Hbins(k,iHeambi)  + etaambi
                Hbins(k,iHnn)     = Hbins(k,iHnn)     + data_out( 7)
                Hbins(k,iHngm)    = Hbins(k,iHngm)    + data_out(12)
                Hbins(k,iHngz)    = Hbins(k,iHngz)    + data_out(13)
                Hbins(k,iHngp)    = Hbins(k,iHngp)    + data_out(14)
                Hbins(k,iHnhion)  = Hbins(k,iHnhion)  + data_out( 8)
                Hbins(k,iHnmion)  = Hbins(k,iHnmion)  + data_out( 9)
                if (etahall > 0.0) then
                   Hbins(k,iHehallp)  = Hbins(k,iHehallp)  + etahall
                else if (etahall < 0.0) then
                   Hbins(k,iHehalln)  = Hbins(k,iHehalln)  + etahall
                endif
             enddo
          endif
       endif
    endif
 enddo parts
 deallocate(etaart)
 !
 angx = 0.0
 angy = 0.0
 angz = 0.0
 do i = 1,nbins-1
    if (ibins(1,i) > 0) then
       Dbins(iDvr:  iDvz,  i) = Dbins(iDvr:  iDvz,  i)/float(ibins(1,i))
       Dbins(iDbr:  iDb,   i) = Dbins(iDbr:  iDb,   i)/float(ibins(1,i))
       Dbins(iDbrat:iDbeta,i) = Dbins(iDbrat:iDbeta,i)/float(ibins(1,i))
       Dbins(iDetaF:iDtemA,i) = Dbins(iDetaF:iDtemA,i)/float(ibins(1,i))
       Dbins(iDnn  :iD,    i) = Dbins(iDnn  :iD,    i)/float(ibins(1,i))
       if (ibins(3,i) > 0)      Dbins(iDetaHp,i) = Dbins(iDetaHp,i)/float(ibins(3,i))
       if (ibins(4,i) > 0)      Dbins(iDetaHn,i) = Dbins(iDetaHn,i)/float(ibins(4,i))
       if (rbins2(i) < rdisc2) then
          angx     = angx + Dbins(iDLx,i)
          angy     = angy + Dbins(iDLy,i)
          angz     = angz + Dbins(iDLz,i)
       endif
       Dbins(iDLx:iDLz,i) = Dbins(iDLx:iDLz,i)/float(ibins(1,i))
    endif
    if (ibins(2,i) > 0) then
       Cbins(2,i) = Cbins(2,i)/float(ibins(2,i))
       Cbins(3,i) = Cbins(3,i)/float(ibins(2,i))
       Cbins(7,i) = sqrt( Cbins(4,i)*Cbins(4,i) + Cbins(5,i)*Cbins(5,i) + Cbins(6,i)*Cbins(6,i))
       Cbins(7,i) = Cbins(7,i)/float(ibins(2,i))
    endif
    if (fracrotDisc(1,i) > 0) then
       fracrotDisc(4,i) = fracrotDisc(4,i)/fracrotDisc(1,i)
       if (fracrotDisc(2,i) > 0) fracrotDisc(5,i) = fracrotDisc(5,i)/fracrotDisc(2,i)
       if (fracrotDisc(3,i) > 0) fracrotDisc(6,i) = fracrotDisc(6,i)/fracrotDisc(3,i)
       fracrotDisc(2,i) = fracrotDisc(2,i)/fracrotDisc(1,i)
    endif
 enddo
 ang  = sqrt(angx*angx + angy*angy + angz*angz)
 volL = sqrt(volP(iVLx)*volP(iVLx) + volP(iVLy)*volP(iVLy) + volP(iVLy)*volP(iVLy) )
 if (Hbins(2,iHn) > 0.0) then
    ang        = ang       /Hbins(2,iHn)
    do i = iHvr,iH
       Hbins(:,i) = Hbins(:,i)/Hbins(:,iHn)
    enddo
    nobj = 0
    do i = 1,nbins-1
       nobj = nobj + ibins(3,i)
    enddo
    if (nobj > 0) Hbins(2,iHehallp) = Hbins(2,iHehallp)/nobj
    nobj = 0
    do i = 1,nbins-1
       nobj = nobj + ibins(4,i)
    enddo
    if (nobj > 0) Hbins(2,iHehallp) = Hbins(2,iHehalln)/nobj
    Hbins(1,iHehallp) = 0.0 ! since we do not have a proper number count
    Hbins(1,iHehalln) = 0.0 ! since we do not have a proper number count
 endif
 if (volN(iVN) > 0) then
    volP(iVb)    = volP(iVb)    / volN(iVN)
    volP(iVbeta) = volP(iVbeta) / volN(iVN)
    volP(iVvphi) = volP(iVvphi) / volN(iVN)
    fracrotVol(iFN) = fracrotVol(iFN)/volN(iVN)
    if (volN(iVNrhoh) > 0) fracrotVol(iFNrhoh) = fracrotVol(iFNrhoh)/volN(iVNrhoh)
    if (volN(iVNrhol) > 0) fracrotVol(iFNrhol) = fracrotVol(iFNrhol)/volN(iVNrhol)
    fracrotVol(iFNvpvra) = fracrotVol(iFNvpvra)/volN(iVN)
 endif
 if (volN(iVNvphip) > 0) volP(iVvphip) = volP(iVvphip)/volN(iVNvphip)
 if (volN(iVNvphin) > 0) volP(iVvphin) = volP(iVvphin)/volN(iVNvphin)
 !
 ! Write results to file
 if ( printrad ) then
    if (csink=="000") then
       write(fileout4,'(3a)') 'rhosurf_',trim(dumpfile),'.dat'
       write(fileout5,'(3a)') 'rhosurfB_',trim(dumpfile),'.dat'
    else
       write(fileout4,'(5a)') 'rhosurf_',trim(dumpfile),'_S',csink,'.dat'
       write(fileout5,'(5a)') 'rhosurfB_',trim(dumpfile),'_S',csink,'.dat'
    endif
    open(punit,file=fileout4)
    open(bunit,file=fileout5)
    call write_header_file4(punit)
    call write_header_file5(bunit)
    Rey  = 0.
    ReyK = dmassp
    massint = 0.0
    do i = 1,nbins-1
       massint = massint + ibins(1,i)*particlemass*umass/solarm
       ReyK(1:6) = 0.
       do j = iDetaF,iDetaHn
          if (abs(Dbins(j,i)) > 0.0) then
             ReyK(j) = sqrt( ReyK(7)*sqrt(rbins2(i)) )/abs(Dbins(j,i))
             if (rbins2(i) < rdisc2) Rey(j)  = Rey(j) + ibins(1,i)*ReyK(j)
          endif
       enddo
       ReyK(7) = ReyK(7) + ibins(1,i)*particlemass
       if (rbins2(i) < rdisc2) Rey(7) = Rey(7) + ibins(1,i)
       ! Radial data at a given time
       write(punit,'(47(1pe18.10,1x))') sqrt(rbins2(i))*udist/au,Dbins(iDvr:iDvz,i)*unit_velocity &
                                 , Dbins(iDbr:iDbz,i)*unit_Bfield,Dbins(iDbrat:iDbeta,i) &
                                 , Dbins(iDLx:iDLz,i)*udist*unit_velocity,ibins(1,i)*particlemass*umass/solarm &
                                 , ibins(2,i)*particlemass*umass/solarm,Cbins(2,i)*unit_Bfield &
                                 , Cbins(3,i),Cbins(7,i)*udist*unit_velocity &
                                 , fracrotDisc(2,i),fracrotDisc(4:6,i),Dbins(iDetaF:iDetaHn,i)*unit_eta,ReyK(1:6) &
                                 , Dbins(iDtemA:iDtemX,i),massint,Dbins(iDetaOA:iDetaAA,i) &
                                 , Dbins(iDnn:iDnmion,i)/udist**3
       write(bunit,'(9(1pe18.10,1x))') sqrt(rbins2(i))*udist/au,Dbins(iDbr:iDb,i)*unit_Bfield,Dbins(iDbrat,i)*unit_Bfield
    enddo
    close(punit)
    close(bunit)
 endif
 !
 ! Write time averaged quantities: disc properties & all gas with rho > rho_crit
 if (Rey(7) > 0) Rey(1:6) = Rey(1:6)/Rey(7)
 write(junit,'(I18,1x,44(1pe18.10,1x))') num, time, Hbins(2,iHv:iHvphi)*unit_velocity,Hbins(2,iHb:iHbx)*unit_Bfield,&
                                     Hbins(2,iHbeta),ang*udist*unit_velocity,Hbins(2,iHeart:iHehalln)*unit_eta,&
                                     Rey(1:3),Rey(6),Hbins(2,iHnn:iHnmion)/udist**3,&
                                     Hbins(1,iHv:iHvphi)*unit_velocity,Hbins(1,iHb:iHbx)*unit_Bfield,Hbins(1,iHbeta),&
                                     Hbins(1,iHeart:iHeambi)*unit_eta
 if ( printvol ) then
    ! Write time averaged quantities: volume properties
    write(kunit,'(I18,1x,16(1pe18.10,1x))') num, time, float(volN(iVN))*particlemass*umass/solarm, &
                                     dmassp*umass/solarm,volP(iVb)*unit_Bfield, &
                                     volP(iVvphi)*unit_velocity,volP(iVvphip)*unit_velocity,&
                                     volP(iVvphin)*unit_velocity,volL*udist*unit_velocity,volP(iVbeta),&
                                     fracrotVol,float(volN(iVNrhoh))/float(volN(iVN))
 endif
 !
end subroutine
!----------------------------------------------------------------
!+
!  Resets the particle origin and velocities to that of a given location
!+
!----------------------------------------------------------------
subroutine reset_origin(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,xyz_tmp,vxyz_tmp)
 integer, intent(in)    :: npart,nptmass
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(in)    :: xyz_tmp(:),vxyz_tmp(:)
 integer                :: i

 do i = 1,npart
    xyzh (1:3,i) = xyzh (1:3,i) - xyz_tmp
    vxyzu(1:3,i) = vxyzu(1:3,i) - vxyz_tmp
 enddo
 do i = 1,nptmass
    xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - xyz_tmp
    vxyz_ptmass (1:3,i) = vxyz_ptmass (1:3,i) - vxyz_tmp
 enddo

end subroutine reset_origin
!----------------------------------------------------------------
!+
!  Adjusts the origin to the CoM of the a sphere of radius rcom
!  about the current origin; only using sink particles and particles
!  with rho > 10rhothresh
!+
!----------------------------------------------------------------
subroutine adjust_origin(npart,nptmass,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,xyz0,vxyz0,dthresh,pmassi)
 integer, intent(in)    :: npart,nptmass
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(out)   :: xyz0(3),vxyz0(3)
 real,    intent(in)    :: dthresh,pmassi
 integer                :: i
 real                   :: mcom,rhoi,hi,pmass,rtmp2

 mcom  = 0.0
 xyz0  = 0.0
 vxyz0 = 0.0
 do i = 1,nptmass
    rtmp2 = dot_product(xyzmh_ptmass(1:3,i),xyzmh_ptmass(1:3,i))
    if (rtmp2 < rcom2) then
       pmass = xyzmh_ptmass(4,i)
       mcom  =  mcom + pmass
       xyz0  =  xyz0 + xyzmh_ptmass(1:3,i)*pmass
       vxyz0 = vxyz0 +  vxyz_ptmass(1:3,i)*pmass
    endif
 enddo
 do i = 1,npart
    hi = xyzh(4,i)
    if (hi > 0.0) then
       rhoi = rhoh(hi,pmassi)
       if (rhoi > 10.*dthresh) then
          rtmp2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
          if (rtmp2 < rcom2) then
             mcom  =  mcom + pmassi
             xyz0  =  xyz0 + xyzh(1:3,i)*pmassi
             vxyz0 = vxyz0 + vxyzu(1:3,i)*pmassi
          endif
       endif
    endif
 enddo
 if (mcom > 0.0) then
    xyz0  =  xyz0/mcom
    vxyz0 = vxyz0/mcom
    write(*,'(a,3Es18.6)') "shifting the CoM by ",xyz0
 endif

end subroutine adjust_origin
!----------------------------------------------------------------
!+
!  Writes headers
!+
!----------------------------------------------------------------
! Write header for fileout1 = analysisout_*_discRM.dat
subroutine write_header_file1(inunit)
 integer, intent(in)  :: inunit
 write(inunit,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'num',      &
        2,'time',     &
        3,'m_low',    &
        4,'m_sink',   &
        5,'m_disc',   &
        6,'m_system', &
        7,'r_disc'
end subroutine write_header_file1
!
! Write header for fileout2 = analysisout_*_discRMnx.dat
subroutine write_header_file2(inunit)
 integer, intent(in)  :: inunit
 write(inunit,"('#',45(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'num',       &
        2,'time',      &
        3,'D v_ave',   & ! D == these values are computed only using gas in the disc
        4,'D v_r',     &
        5,'D v_phi',   &
        6,'D B_ave',   &
        7,'D B_r',     &
        8,'D B_phi',   &
        9,'D B_z',     &
       10,'D B_y',     &
       11,'D B_x',     &
       12,'D beta_P',  &
       13,'D ang mom', &
       14,'D eta_art', &
       15,'D eta_ohm', &
       16,'D eta_hall',&
       17,'D |e_hall|',&
       18,'D eta_ambi',&
       19,'D eta_h>0', &
       20,'D eta_h<0', &
       21,'D R_art',   &
       22,'D R_ohm',   &
       23,'D R_hall',  &
       24,'D R_ambi',  &
       25,'D n_n',     &
       26,'D n_g(Z=-1)',&
       27,'D n_g(Z=0)', &
       28,'D n_g(Z=+1)',&
       29,'D n_h-ion', &
       30,'D n_m-ion', &
       31,'G v_ave',   & ! G == these values are computed using all local gas with rho > rho_crit
       32,'G v_r',     &
       33,'G v_phi',   &
       34,'G B_ave',   &
       35,'G B_r',     &
       36,'G B_phi',   &
       37,'G B_z',     &
       38,'G B_y',     &
       39,'G B_x',     &
       40,'G beta_P',  &
       41,'G eta_art', &
       42,'G eta_ohm', &
       43,'G eta_hall',&
       44,'G |e_hall|',&
       45,'G eta_ambi'
end subroutine write_header_file2
!
! Write header for fileout3 = analysisout_*_vol*RM.dat
subroutine write_header_file3(inunit)
 integer, intent(in)  :: inunit
 write(inunit,"('#',12(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'num',      &
        2,'time',     &
        3,'m_low',    &
        4,'m_sink',   &
        5,'m_disc',   &
        6,'m_d,gas',  &
        7,'m_d,dust', &
        8,'m_d,star', &
        9,'r_disc',   &
       10,'r_d,gas',  &
       11,'r_d,dust', &
       12,'r_d,star'
end subroutine write_header_file3
!
! Write header for fileout4 = rhosurf_*.dat
subroutine write_header_file4(inunit)
 integer, intent(in)  :: inunit
 write(inunit,"('#',47(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',          &
        2,'v_r',        &
        3,'v_phi',      &
        4,'v_z',        &
        5,'B_r',        &
        6,'B_phi',      &
        7,'B_p',        &
        8,'B_z',        &
        9,'B_phi/B_p',  &
       10,'beta_P',     &
       11,'angmomX',    &
       12,'angmomY',    &
       13,'angmomZ',    &
       14,'mass',       &
       15,'mass (r<X)', &
       16,'B (r<X)',    &
       17,'b_P(r<X)',   &
       18,'L (r<X)',    &
       19,'f high rho', &
       20,'f v_phi > 0',&
       21,'f_hd v_p>0', &
       22,'f_ld v_p>0', &
       23,'eta_art',    &
       24,'eta_ohm',    &
       25,'eta_hall',   &
       26,'|eta_hall|', &
       27,'eta_ambi',   &
       28,'eta_h>0',    &
       29,'eta_h<0',    &
       30,'R_art',      &
       31,'R_ohm',      &
       32,'R_hall',     &
       33,'R_ambi',     &
       34,'R_hall>0',   &
       35,'R_hall<0',   &
       36,'T (max)',    &
       37,'T (ave)',    &
       38,'M(r)',       &
       39,'e_O/e_a',    &
       40,'|e_H|/e_a',  &
       41,'e_A/e_a',    &
       42,'n_n',        &
       43,'n_g(Z=-1)',  &
       44,'n_g(Z=0)',   &
       45,'n_g(Z=+1)',  &
       46,'n_{h-ion}',  &
       47,'n_{m-ionm}'
end subroutine write_header_file4
!
! Write header for fileout5 = rhosurfM_*.dat
subroutine write_header_file5(inunit)
 integer, intent(in)  :: inunit
 write(inunit,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',    &
        2,'B_r',  &
        3,'B_phi',&
        4,'B_p',  &
        5,'B_z',  &
        6,'B_y',  &
        7,'B_x',  &
        8,'B',    &
        9,'Bphi/Bp'
end subroutine write_header_file5
!
! Write header for fileout6 = analysisout_*_eta.dat
subroutine write_header_file6(inunit)
 integer, intent(in)  :: inunit
 write(inunit,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'num',        &
        2,'time',       &
        3,'eta_art',    &
        4,'eta_ohm',    &
        5,'eta_hall',   &
        6,'|eta_hall|', &
        7,'eta_ambi',   &
        8,'nptmass',    &
        9,'nptmassPrev'
end subroutine write_header_file6
!
! Write header for fileout7 = analysisout_*_mu.dat
subroutine write_header_file7(inunit)
 integer, intent(in)  :: inunit
 write(inunit,"('#',22(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'num',        &
        2,'time',       &
        3,'mu_g(2680au',&
        4,'mu_g(1206au',&
        5,'mu_g(500au)',&
        6,'mu_g(200au)',&
        7,'mu_g(Dg)',   &
        8,'Dg (au)',    &
        9,'mu_1( 30au)',&
       10,'mu_1( 60au)',&
       11,'mu_1( 90au)',&
       12,'mu_1(120au)',&
       13,'mu_1(200au)',&
       14,'mu_1(D1)',   &
       15,'D1 (au)',    &
       16,'mu_2( 30au)',&
       17,'mu_2( 60au)',&
       18,'mu_2( 90au)',&
       19,'mu_2(120au)',&
       20,'mu_2(200au)',&
       21,'mu_2(D2)',   &
       22,'D2 (au)'
end subroutine write_header_file7
!-----------------------------------------------------------------------
!
end module
