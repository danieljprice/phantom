!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! This will calculate relevant sink & gas stas for our converging flow study
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: dim, part, physcon, sortutils, units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'ConvergStats'
 integer, parameter :: nbins_max = 256
 integer, parameter :: nres      =   5
 integer            :: num0,iloop,nbins(nres)
 real               :: vmin,vmax,mmin,mmax,dmin,dmax,dtmax0,mass01,tnext,dtsh
 real               :: dv(nres),dm(nres),dd(nres)
 real               :: vmd(6,nbins_max,nres)
 logical, private   :: firstcall = .true.

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,          only: mhd
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,rhoh,isdead_or_accreted
 use physcon,      only: solarm,km
 use units,        only: unit_velocity,unit_density,umass
 use sortutils,    only: indexx
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i,j,k,imass0,nptmass_alive
 integer                      :: mvbins(5,nbins_max,nres),llsink(nptmass)
 real                         :: vmin_kms,vmax_kms,mmin_msun,mmax_msun,dmin_cgs,dmax_cgs,dtsh_cgs
 real                         :: vel,msink,rhoi,fcheck,mtot_sink,mtotl_sink,mave_sink,mavel_sink,std,stdl
 real                         :: msink34,msink50,msink68,mtot
 real                         :: vbinsr(4,nbins_max,nres),vbinsrNorm(4,nres)
 real                         :: mvbinsNorm(5,nres)
 logical                      :: iexist
 character(len=200)           :: fileoutSV,fileoutSH,fileoutGH,fileoutSA,wfmt

 vmin_kms  = -35.
 vmax_kms  = 35.
 mmin_msun = 1.
 mmax_msun = 1.0d8
 dmin_cgs  = 1.d-30
 dmax_cgs  = 1.d-10
 dtsh_cgs  = 1.d-23          ! will only bin gas with density higher than this threshold
 dtmax0    = 0.0335296079237 ! get the initial dump separations; this is simulation dependent
 !
 ! Initialise values & Open file
 !
 fileoutSV = trim(dumpfile(1:index(dumpfile,'_')-1))//'_SinkVel.dat'
 fileoutSA = trim(dumpfile(1:index(dumpfile,'_')-1))//'_SinkAverages.dat'
 inquire(file=fileoutSV,exist=iexist)
 if ( .not.iexist .or. firstcall ) then
    firstcall = .false.
    num0      = num
    tnext     = dtmax0
    iloop     = 0
    vmin = vmin_kms*km/unit_velocity
    vmax = vmax_kms*km/unit_velocity
    mmin = mmin_msun*solarm/umass
    mmax = mmax_msun*solarm/umass
    dmin = dmin_cgs/unit_density
    dmax = dmax_cgs/unit_density
    dtsh = dtsh_cgs/unit_density
    vmd  = 0.
    do k = 1,nres
       nbins(k) = nbins_max/2**(k-1)
       dv(k) = (vmax-vmin)/nbins(k)
       dm(k) = (log10(mmax)-log10(mmin))/nbins(k)
       dd(k) = (log10(dmax)-log10(dmin))/nbins(k)
       ! bin edges
       do i = 1,nbins(k)
          vmd(1,i,k) = vmin + (i-1)*dv(k)
          vmd(2,i,k) = 10**(log10(mmin) +(i-1)*dm(k) )
          vmd(3,i,k) = 10**(log10(dmin) +(i-1)*dd(k) )
       enddo
       ! bin centres for plotting
       do i = 1,nbins(k)-1
          vmd(4,i,k) = 0.5*(vmd(1,i,k) + vmd(1,i+1,k)) * unit_velocity/km
          vmd(5,i,k) = 10**(0.5*(log10(vmd(2,i,k)) + log10(vmd(2,i+1,k)) )) * umass/solarm
          vmd(6,i,k) = 10**(0.5*(log10(vmd(3,i,k)) + log10(vmd(3,i+1,k)) )) * unit_density
       enddo
    enddo
    imass0 = 0
    do i = 1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (rhoh(xyzh(4,i),particlemass) > (1.1d-23)/unit_density) imass0 = imass0+1
       endif
    enddo
    print*, "The total initial mass of the clouds is ",imass0*particlemass," M_sun"
    mass01 = 1./(imass0*particlemass)
    ! open files
    open(iunit,file=fileoutSV,status='replace')
    write(iunit,  "('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'idump',   &
          2,'time/dt', &
          3,'time',    &
          4,'ID',      &
          5,'mass',    &
          6,'vx',      &
          7,'vy',      &
          8,'vz'
    do k = 1,nres
       write(fileoutSH,'(2a,I3.3,a)') trim(dumpfile(1:index(dumpfile,'_')-1)),'_SinkHisto',nbins(k),'.dat'
       open(iunit+(2*k+1),file=fileoutSH,status='replace')
       write(iunit+(2*k+1),"('#',10(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'idump',   &
          2,'time/dt', &
          3,'time',    &
          4,'mass bin',&
          5,'N mass',  &
          6,'vel bin', &
          7,'N vel' ,  &
          8,'N velx',  &
          9,'N vely',  &
         10,'N velz'
       write(fileoutGH,'(2a,I3.3,a)') trim(dumpfile(1:index(dumpfile,'_')-1)),'_GasHisto',nbins(k),'.dat'
       open(iunit+(2*k+2),file=fileoutGH,status='replace')
       write(iunit+(2*k+2),"('#',10(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'idump',      &
          2,'time/dt',    &
          3,'time',       &
          4,'density bin',&
          5,'N mass',     &
          6,'vel bin',    &
          7,'N vel' ,     &
          8,'N velx',     &
          9,'N vely',     &
         10,'N velz'
    enddo
    open(iunit+20,file=fileoutSA,status='replace')
    write(iunit+20,  "('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'idump',     &
          2,'time/dt',   &
          3,'time',      &
          4,'Nsink',     &
          5,'Msink_ave', &
          6,'Msink_std', &
          7,'Msink_avel',&
          8,'Msink_stdl',&
          9,'Msink_mean'

 else
    iloop = iloop + 1
    tnext = time + dtmax0
    open(iunit   ,file=fileoutSV,position='append')
    open(iunit+20,file=fileoutSA,position='append')
    do k = 1,nres
       write(fileoutSH,'(2a,I3.3,a)') trim(dumpfile(1:index(dumpfile,'_')-1)),'_SinkHisto',nbins(k),'.dat'
       write(fileoutGH,'(2a,I3.3,a)') trim(dumpfile(1:index(dumpfile,'_')-1)),'_GasHisto',nbins(k),'.dat'
       open(iunit+(2*k+1),file=fileoutSH,position='append')
       open(iunit+(2*k+2),file=fileoutGH,position='append')
    enddo
 endif

 ! print sink stats
 print*, 'printing sink stats'
 do i = 1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       write(iunit,'(2(I18,1x),es18.10,I18,1x,4(es18.10,1x))') &
       num,iloop,time,i,xyzmh_ptmass(4,i),vxyz_ptmass(1:3,i)*unit_velocity/km
    endif
 enddo

 ! determine sink mass & velocity distribution
 print*, 'printing sink distribuitions'
 nptmass_alive = 0
 mvbins     = 0
 vbinsr     = 0.
 mtot_sink  = 0
 mtotl_sink = 0
 do i = 1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       nptmass_alive = nptmass_alive + 1
       msink      = xyzmh_ptmass(4,i)
       mtot_sink  = mtot_sink  + msink
       mtotl_sink = mtotl_sink + log10(msink)
       vel        = sqrt(dot_product(vxyz_ptmass(1:3,i),vxyz_ptmass(1:3,i)))
       do k = 1,nres
          do j = 1,nbins(k)-1
             if (vmd(1,j,k) < vel              .and. vel              < vmd(1,j+1,k) ) vbinsr(4,j,k) = vbinsr(4,j,k) + msink
             if (vmd(1,j,k) < vxyz_ptmass(3,i) .and. vxyz_ptmass(3,i) < vmd(1,j+1,k) ) vbinsr(3,j,k) = vbinsr(3,j,k) + msink
             if (vmd(1,j,k) < vxyz_ptmass(2,i) .and. vxyz_ptmass(2,i) < vmd(1,j+1,k) ) vbinsr(2,j,k) = vbinsr(2,j,k) + msink
             if (vmd(1,j,k) < vxyz_ptmass(1,i) .and. vxyz_ptmass(1,i) < vmd(1,j+1,k) ) vbinsr(1,j,k) = vbinsr(1,j,k) + msink
             if (vmd(1,j,k) < vxyz_ptmass(3,i) .and. vxyz_ptmass(3,i) < vmd(1,j+1,k) ) mvbins(5,j,k) = mvbins(5,j,k) + 1
             if (vmd(1,j,k) < vxyz_ptmass(2,i) .and. vxyz_ptmass(2,i) < vmd(1,j+1,k) ) mvbins(4,j,k) = mvbins(4,j,k) + 1
             if (vmd(1,j,k) < vxyz_ptmass(1,i) .and. vxyz_ptmass(1,i) < vmd(1,j+1,k) ) mvbins(3,j,k) = mvbins(3,j,k) + 1
             if (vmd(1,j,k) < vel              .and. vel              < vmd(1,j+1,k) ) mvbins(2,j,k) = mvbins(2,j,k) + 1
             if (vmd(2,j,k) < msink            .and. msink            < vmd(2,j+1,k) ) mvbins(1,j,k) = mvbins(1,j,k) + 1
          enddo
       enddo
    endif
 enddo
 ! properly normalise the bins
 mvbinsNorm = 0.
 vbinsrNorm = 0.
 do k = 1,nres
    do i = 1,nbins(k)-1
       mvbinsNorm(1,  k) = mvbinsNorm(1,  k) + mvbins(1,  i,k)*dm(k)
       mvbinsNorm(2:5,k) = mvbinsNorm(2:5,k) + mvbins(2:5,i,k)*dv(k)
       vbinsrNorm(1:4,k) = vbinsrNorm(1:4,k) + vbinsr(1:4,i,k)*dv(k)
    enddo
 enddo
 mvbinsNorm        = mvbinsNorm/nptmass_alive
 vbinsrNorm        = vbinsrNorm/mtot_sink
 mvbinsNorm(2:5,:) = mvbinsNorm(2:5,:)*unit_velocity/km
 vbinsrNorm(1:4,:) = vbinsrNorm(1:4,:)*unit_velocity/km

 ! print to file
 do k = 1,nres
    do i = 1,nbins(k)
       if (nptmass_alive > 0) then
          write(iunit+(2*k+1),'(2(I18,1x),8(es18.6,1x))') &
          num,iloop,time,vmd(5,i,k),mvbins(1,i,k)/mvbinsNorm(1,k), &
                         vmd(4,i,k),mvbins(2:5,i,k)/mvbinsNorm(2:5,k)
       endif
    enddo
 enddo

 ! determine the standard deviation
 std  = 0.
 stdl = 0.
 mave_sink  = mtot_sink/nptmass_alive
 mavel_sink = mtotl_sink/nptmass_alive
 do i = 1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       msink = xyzmh_ptmass(4,i)
       std   = std  + (msink-mave_sink)**2
       stdl  = stdl + (log10(msink)-mavel_sink)**2
    else
       xyzmh_ptmass(4,i) = 50000.
    endif
 enddo
 std  = sqrt(std /nptmass_alive)
 stdl = sqrt(stdl/nptmass_alive)

 ! manually count the mean & standard deviations to permit asymertic results.
 call indexx(nptmass,xyzmh_ptmass(4,1:nptmass),llsink)
 msink34 = xyzmh_ptmass(4,llsink(int(0.341*nptmass_alive)))
 msink50 = xyzmh_ptmass(4,llsink(int(0.500*nptmass_alive)))
 msink68 = xyzmh_ptmass(4,llsink(int(0.682*nptmass_alive)))

 ! print results
 write(wfmt,'(a)') '(2(I18,1x),(es18.6,1x),(I18,1x),5(es18.6,1x))'
 write(iunit+20,wfmt) num,iloop,time,nptmass_alive,mave_sink-std,std,mavel_sink-stdl,stdl,msink34
 write(iunit+20,wfmt) num,iloop,time,nptmass_alive,mave_sink    ,std,mavel_sink     ,stdl,msink50
 write(iunit+20,wfmt) num,iloop,time,nptmass_alive,mave_sink+std,std,mavel_sink+stdl,stdl,msink68
 write(iunit+20,'(a)') ' '

 ! determine gas density mass & velocity distribution
 print*, 'calculating gas distributions'
 mvbins = 0
 mtot   = 0.
!$omp parallel do default(none) &
!$omp shared(npart,xyzh,vxyzu,particlemass,nbins,vmd,dtsh) &
!$omp private(i,j,k,vel,rhoi) &
!$omp reduction(+:mvbins,mtot)
 do i = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       rhoi = rhoh(xyzh(4,i),particlemass)
       if (rhoi > dtsh) then
          vel  = sqrt(dot_product(vxyzu(1:3,i),vxyzu(1:3,i)))
          mtot = mtot + 1.0
          do k = 1,nres
             do j = 1,nbins(k)-1
                if (vmd(1,j,k) < vxyzu(3,i) .and. vxyzu(3,i) < vmd(1,j+1,k) ) mvbins(5,j,k) = mvbins(5,j,k) + 1
                if (vmd(1,j,k) < vxyzu(2,i) .and. vxyzu(2,i) < vmd(1,j+1,k) ) mvbins(4,j,k) = mvbins(4,j,k) + 1
                if (vmd(1,j,k) < vxyzu(1,i) .and. vxyzu(1,i) < vmd(1,j+1,k) ) mvbins(3,j,k) = mvbins(3,j,k) + 1
                if (vmd(1,j,k) < vel        .and. vel        < vmd(1,j+1,k) ) mvbins(2,j,k) = mvbins(2,j,k) + 1
                if (vmd(3,j,k) < rhoi       .and. rhoi       < vmd(3,j+1,k) ) mvbins(1,j,k) = mvbins(1,j,k) + 1
             enddo
          enddo
       endif
    endif
 enddo
!omp end parallel do

 ! properly normalise the bins
 mvbinsNorm = 0.
 do k = 1,nres
    do i = 1,nbins(k)
       mvbinsNorm(1,  k) = mvbinsNorm(1,  k) + mvbins(1  ,i,k)*dd(k)
       mvbinsNorm(2:5,k) = mvbinsNorm(2:5,k) + mvbins(2:5,i,k)*dv(k)
    enddo
 enddo
 mtot = mtot*particlemass
 mvbinsNorm = mvbinsNorm!*mass01
 mvbinsNorm = mvbinsNorm/mtot
 mvbinsNorm(2:5,:) = mvbinsNorm(2:5,:)*unit_velocity/km

 print*, 'printing gas distributions'
 do k = 1,nres
    fcheck = 0.
    do i = 1,nbins(k)-1
       write(iunit+(2*k+2),'(2(I18,1x),8(es18.10,1x))') num,iloop,time, &
       vmd(6,i,k),mvbins(1,i,k)/mvbinsNorm(1,k), &
       vmd(4,i,k),mvbins(2:5,i,k)/mvbinsNorm(2:5,k)
       fcheck = fcheck + mvbins(1,i,k)/mvbinsNorm(1,k)
    enddo
    print*, 'The total fraction of gas accounted for is ',fcheck
 enddo
 close(iunit)
 close(iunit+20)
 do k = 1,nres
    close(iunit+(2*k+1))
    close(iunit+(2*k+2))
 enddo

end subroutine do_analysis

end module analysis
