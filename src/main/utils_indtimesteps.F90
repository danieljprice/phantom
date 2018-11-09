!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: timestep_ind
!
!  DESCRIPTION:
!  Parameters and routines related to
!  individual particle timesteps
!  (routines should ONLY be called if -DIND_TIMESTEPS set)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, mpiutils, part
!+
!--------------------------------------------------------------------------
module timestep_ind
 use dim, only:maxp
 implicit none
 integer(kind=1)    :: nbinmax
 integer(kind=1)    :: ibinnow = 0  ! required to prevent compiler warning
 integer            :: istepfrac
 integer            :: nactive
 integer(kind=8)    :: nactivetot
 integer, parameter :: maxbins = 30 ! limit of 2**(n+1) with reals
 integer, parameter :: itdt    = 1
 integer, parameter :: ithdt   = 2
 integer, parameter :: itdt1   = 3
 integer, parameter :: ittwas  = 4

contains

!----------------------------------------------------------------
!+
!  get timestep from ibin setting
!+
!----------------------------------------------------------------
pure real function get_dt(dtmax,ibini)
 real,            intent(in) :: dtmax
 integer(kind=1), intent(in) :: ibini

 get_dt = dtmax/2**ibini

end function get_dt

#ifdef IND_TIMESTEPS
!----------------------------------------------------------------
!+
!  If dt is read in from a dump file, then initialise ibin & ibin_old.
!  Although always called if dt is read in, it is only necessary
!  for restarts.  This is necessary to prevent nearby particles from
!  being initialised with an ibin that is very different than that of
!  their neighbours.
!+
!----------------------------------------------------------------
subroutine init_ibin(npart,dtmax)
 use part, only: ibin,ibin_old,dt_in
 integer, intent(in) :: npart
 real,    intent(in) :: dtmax
 real(kind=4)        :: twoepsilon,dt_ini
 integer             :: i,j
 logical             :: find_bin

 twoepsilon = 2.0*epsilon(dt_in(1))
!$omp parallel default(none) &
!$omp shared(npart,dtmax,twoepsilon,dt_in,ibin,ibin_old) &
!$omp private(i,j,dt_ini,find_bin)
!$omp do
 do i = 1,npart
    j        = 0
    dt_ini = dt_in(i)
    if (dt_ini > twoepsilon) then
       find_bin = .true.
    else
       find_bin = .false.
    endif
    do while (find_bin .and. j < maxbins)
       if ( dtmax/2**j + twoepsilon > dt_ini) then
          ibin(i)     = int(j,kind=1)
          ibin_old(i) = int(j,kind=1)
       else
          find_bin   = .false.
       endif
       j = j + 1
    enddo
 enddo
!$omp enddo
!$omp end parallel

end subroutine init_ibin
!----------------------------------------------------------------
!+
!  routine to set iactive flag determining whether a particle's
!  forces are to be evaluated on the current timestep.
!+
!----------------------------------------------------------------
subroutine set_active_particles(npart,nactive,nalive,iphase,ibin,xyzh)
 use io,   only:iprint,fatal
 use part, only:isdead_or_accreted,iamtype,isetphase,maxp,all_active,iboundary
 integer,         intent(in)    :: npart
 integer,         intent(out)   :: nactive,nalive
 integer(kind=1), intent(inout) :: iphase(maxp),ibin(maxp)
 real,            intent(in)    :: xyzh(4,maxp)
 integer                        :: i,itype
 integer(kind=1)                :: ibini
 logical                        :: iactivei

 nactive = 0
 nalive  = 0
!$omp parallel default(none) &
!$omp shared(npart,nbinmax,ibin,iprint,istepfrac,iphase,xyzh) &
!$omp private(i,itype,iactivei,ibini) &
!$omp reduction(+:nactive,nalive)
!$omp do schedule(guided,10)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       nalive = nalive + 1
       itype  = iamtype(iphase(i))
       ! boundary particles are never active
       if (itype==iboundary) ibin(i) = 0
       ibini = ibin(i)
       !--sanity check
       if (ibini > nbinmax) then
          write(iprint,*) 'ERROR in timesteps: ibin = ',ibini,' nbinmax = ',nbinmax
          call fatal('set_active_particles','timestep bin exceeds max bins')
       endif
       if (mod(istepfrac,2**(nbinmax-ibini))==0) then
          iactivei = .true.
          nactive  = nactive + 1
       else
          iactivei = .false.
       endif
       iphase(i) = isetphase(itype,iactivei)
    endif
 enddo
!$omp enddo
!$omp end parallel
 !
 !--Determine the current maximum ibin that is active
 !
 i       = 0
 ibinnow = nbinmax
 do while (ibinnow == nbinmax .and. i < nbinmax)
    if (mod(istepfrac,2**(nbinmax-i))==0) ibinnow = int(i,kind=1)
    i = i + 1
 enddo
 !
 !--Determine activity to determine if stressmax needs to be reset
 !
 if (nactive==nalive) then
    all_active = .true.
 else
    all_active = .false.
 endif

end subroutine set_active_particles
#endif

!----------------------------------------------------------------
!+
!  routine to do bookwork needed when total number of bins changes
!  (dt changes and istepfrac is correspondingly adjusted)
!+
!----------------------------------------------------------------
subroutine change_nbinmax(nbinmax,nbinmaxprev,istepfrac,dtmax,dt)
 use io, only:iprint,id,master,iverbose,die,fatal
 integer(1), intent(in)    :: nbinmax,nbinmaxprev
 integer,    intent(inout) :: istepfrac
 real,       intent(in)    :: dtmax
 real,       intent(out)   :: dt

 if (id==master .and. iverbose >= 2) then
    write(iprint,"(1x,45('*'))")
    write(iprint,*) '* current step fraction = ',istepfrac,' / ',2**nbinmaxprev
 endif
 !
 !--the number of bins should only decrease when bins are synchronised
 !
 if (nbinmax < nbinmaxprev .and. mod(istepfrac,2**(nbinmaxprev-nbinmax)) /= 0) then
    write(iprint,*) 'error: istepfrac not multiple of ',2**(nbinmaxprev-nbinmax),' when changing nbinmax from ', &
                    nbinmaxprev,'->',nbinmax
    call die
 endif
 if (nbinmax > nbinmaxprev) then
    istepfrac = istepfrac*2**(nbinmax-nbinmaxprev)
    !print*,' adjust factor = ',2**(nbinmax-nbinmaxprev)
 else
    istepfrac = istepfrac/2**(nbinmaxprev - nbinmax)
    ! print*,' adjust factor = ',2**(nbinmaxprev-nbinmax)
 endif
 if (nbinmax > maxbins) call fatal('change_nbinmax','step too small: maximum number of timestep bins (30) exceeded')
 dt = dtmax/2**nbinmax

 if (id==master .and. iverbose >= 2) then
    write(iprint,*) '* number of timestep bins old = ',nbinmaxprev,', new = ',nbinmax
    write(iprint,*) '* new istepfrac = ',istepfrac,' / ',2**nbinmax
    write(iprint,"(1x,45('*'))")
 endif

 return
end subroutine change_nbinmax

!----------------------------------------------------------------
!+
!  routine to place particle in new timestep bin based on dt
!+
!----------------------------------------------------------------
subroutine get_newbin(dti,dtmax,ibini,allow_decrease,limit_maxbin,dtchar)
 use io, only:fatal
 real,            intent(in)    :: dti,dtmax
 integer(kind=1), intent(inout) :: ibini
 logical,         intent(in), optional :: allow_decrease,limit_maxbin
 character(len=*),intent(in), optional :: dtchar
 integer(kind=1) :: ibin_oldi
 integer         :: ibin_newi
 logical         :: iallow_decrease,ilimit_maxbin
 real, parameter :: vsmall = epsilon(vsmall)
 real, parameter :: dlog2 = 1.4426950408889634d0 ! dlog2 = 1./log(2.)

 if (present(allow_decrease)) then
    iallow_decrease = allow_decrease
 else
    iallow_decrease = .true.
 endif

 ! for super-timestepping only
 if (present(limit_maxbin)) then
    ilimit_maxbin = limit_maxbin
 else
    ilimit_maxbin = .true.
 endif

 ibin_oldi = ibini
 if (dti > dtmax) then
    ibin_newi = 0
 elseif (dti < tiny(dti)) then
    ibin_newi = maxbins
 else
    ibin_newi = max(int(log(2.*dtmax/dti)*dlog2-vsmall),0)
 endif
 if (ibin_newi > maxbins .and. ilimit_maxbin) then
    if (present(dtchar)) then
       write(*,'(a,Es16.7)') 'get_newbin: dt_ibin(0)   = ', dtmax
       write(*,'(a,Es16.7)') 'get_newbin: dt_ibin(max) = ', dtmax/2**(maxbins-1)
       write(*,'(2a)'      ) 'get_newbin: dt = ', dtchar
    endif
    call fatal('get_newbin','step too small: bin would exceed maximum',var='dt',val=dti)
 endif

 if (ibin_newi > ibin_oldi) then
    !--timestep can go down at any time
    ibini = int(ibin_newi,kind=1)
 elseif (ibin_newi < ibin_oldi .and. ibin_oldi <= nbinmax .and. iallow_decrease) then
    !--timestep can only go up if bins are synchronised (ie. new bin is active)
    !  move up only one bin at a time
    if (mod(istepfrac,2**(nbinmax-(ibin_oldi-1)))==0) then
       ibini = ibini - 1_1
    endif
    !print*,i,' changing from ',ibin_oldi,' to ',ibini, ' istepfrac = ',istepfrac,' time = ',time
 endif
 !print*,'dti = ',dtmax/2**ibini,dti,'dtmax = ',dtmax,' istep = ',2**ibini

 return
end subroutine get_newbin

!----------------------------------------------------------------
!+
!  Decreasing dtmax
!  This can only occur at full dumps, thus particles are synchronised
!+
!----------------------------------------------------------------
subroutine decrease_dtmax(npart,nbins,time,dtmax_ifactor,dtmax,ibin,ibin_wake,ibin_sts,&
                          ibin_dts)
 integer,         intent(in)    :: npart,nbins,dtmax_ifactor
 integer(kind=1), intent(inout) :: ibin(:),ibin_wake(:),ibin_sts(:)
 real,            intent(in)    :: time
 real,            intent(inout) :: dtmax,ibin_dts(4,0:nbins)
 integer                        :: i
 integer(kind=1)                :: ibin_rat

 !--Determine the new dtmax and update ibin_dts for particle waking
 if (dtmax_ifactor > 0) then
    dtmax    =  dtmax/dtmax_ifactor
    ibin_rat = int((log10(real(dtmax_ifactor)-1.0)/log10(2.0))+1,kind=1)
    do i = 0,nbins-ibin_rat
       ibin_dts(:,i) = ibin_dts(:,i+ibin_rat)
    enddo
    do i = nbins-ibin_rat+1,nbins
       ibin_dts(itdt,  i) = get_dt(dtmax,int(i,kind=1))
       ibin_dts(ithdt, i) = 0.5*ibin_dts(itdt,i)
       ibin_dts(itdt1, i) = 1.0/ibin_dts(itdt,i)
       ibin_dts(ittwas,i) = time + 0.5*get_dt(dtmax,int(i,kind=1))
    enddo
 else if (dtmax_ifactor < 0) then
    dtmax    = -dtmax*dtmax_ifactor
    ibin_rat = -int((log10(real(-dtmax_ifactor)-1.0)/log10(2.0))+1,kind=1)
    do i = nbins,abs(ibin_rat),-1
       ibin_dts(:,i) = ibin_dts(:,i+ibin_rat)
    enddo
    do i = 0,abs(ibin_rat)-1
       ibin_dts(itdt,  i) = get_dt(dtmax,int(i,kind=1))
       ibin_dts(ithdt, i) = 0.5*ibin_dts(itdt,i)
       ibin_dts(itdt1, i) = 1.0/ibin_dts(itdt,i)
       ibin_dts(ittwas,i) = time + 0.5*get_dt(dtmax,int(i,kind=1))
    enddo
 else
    return
 endif
 !
 !--Modify the ibins
!$omp parallel default(none) &
!$omp shared(npart,ibin,ibin_wake,ibin_rat) &
#ifdef STS_TIMESTEPS
!$omp shared(ibin_sts) &
#endif
!$omp private(i)
!$omp do
 do i=1,npart
    ibin(i)      = max(0_1,ibin(i)     -ibin_rat)
    ibin_wake(i) = max(0_1,ibin_wake(i)-ibin_rat)
#ifdef STS_TIMESTEPS
    ibin_sts(i)  = max(0_1,ibin_sts(i) -ibin_rat)
#endif
 enddo
!$omp enddo
!$omp end parallel
 nbinmax = max(0_1,nbinmax-ibin_rat)

end subroutine decrease_dtmax
!----------------------------------------------------------------
!+
!  this subroutine prints a summary of timestep bins
!  for individual particle timesteps
!+
!----------------------------------------------------------------
subroutine write_binsummary(npart,nbinmax,dtmax,timeperbin,iphase,ibin,xyzh)
 use io,       only:iprint,id,master,error
 use part,     only:isdead_or_accreted,iamtype,maxphase,labeltype,maxtypes
 use mpiutils, only:reduce_mpi
 integer,         intent(in) :: npart
 integer(kind=1), intent(in) :: nbinmax
 real,            intent(in) :: dtmax
 real(kind=4),    intent(in) :: timeperbin(0:maxbins)
 integer(kind=1), intent(in) :: iphase(maxp), ibin(maxp)
 real,            intent(in) :: xyzh(4,maxp)
 real(kind=4)       :: timeperbintot(0:maxbins)
 real(kind=4)       :: dtimetot
 integer            :: ninbin(0:nbinmax),itypelist(maxtypes)
 integer            :: noftypeinbin(0:nbinmax,maxtypes)
 integer            :: i,ibini,np,itype,ntypesprint
 integer(kind=8)    :: ntot
 character(len=120) :: fmtstring
 character(len=20)  :: fmtstring2

 if (npart <= 0) return
 ninbin = 0
 noftypeinbin = 0
 np = 0
 over_part: do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       np = np + 1
       ibini = ibin(i)
       if (ibini > nbinmax .or. ibini < 0) then
          call error('write_bin','timestep bin exceeds maximum',var='ibin',ival=ibini)
          cycle over_part
       endif
       ninbin(ibini) = ninbin(ibini) + 1
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
          if (itype > 0 .and. itype <= maxtypes) then
             noftypeinbin(ibini,itype) = noftypeinbin(ibini,itype) + 1
          endif
       endif
    endif
 enddo over_part

 !ninbin(:) = int(reduce_mpi('+',ninbin(:)))
 ntypesprint = 0
 itypelist(:) = 0
 do itype=1,maxtypes
    !noftypeinbin(:,itype) = int(reduce_mpi('+',noftypeinbin(:,itype)))
    if (any(noftypeinbin(:,itype) > 0)) then
       ntypesprint = ntypesprint + 1
       itypelist(ntypesprint) = itype
    endif
 enddo
 ntot = np !reduce_mpi('+',np)
 timeperbintot(:) = reduce_mpi('+',timeperbin(:))

 if (id==master) then
    dtimetot = sum(timeperbintot(0:nbinmax))
    if (dtimetot > 0.) dtimetot = 1._4/dtimetot

    if (maxphase==maxp .and. any(noftypeinbin(:,2:) > 0)) then
       !--multiphase printout
       write(fmtstring2,"(a,i2,a)") '(',17 + 32 + 11*ntypesprint,'(''-''))'
       write(iprint,fmtstring2)
       write(iprint,"(a)",ADVANCE='no') '| bin |    dt    '   ! 17 chars
       write(fmtstring,"(a,i2,a)") '(',ntypesprint,'(a11))'
       write(iprint,fmtstring,ADVANCE='no') ('| n('//trim(labeltype(itypelist(itype)))//')     ',itype=1,ntypesprint)
       write(iprint,"(a)") '|    npart |   frac  | cpufrac |' ! 32 chars
       write(iprint,fmtstring2)
       write(fmtstring,"(a,i2,a)") "('|',i3,3x,es10.3,",ntypesprint+1,"(1x,i10),2x,2pf6.2,'%',3x,2pf6.2,'% |')"
       do i=0,nbinmax
          write(iprint,fmtstring) i,dtmax/2**i,(noftypeinbin(i,itypelist(itype)),itype=1,ntypesprint),&
                                  ninbin(i),ninbin(i)/real(ntot),timeperbintot(i)*dtimetot
       enddo
       write(iprint,fmtstring2)
    else
       write(iprint,"(a)") '-------------------------------------------------'
       write(iprint,"(a)") '| bin |    dt    |   npart  |   frac  | cpufrac |'
       write(iprint,"(a)") '-------------------------------------------------'
       do i=0,nbinmax
          write(iprint,10) i,dtmax/2**i,ninbin(i),100.*ninbin(i)/real(ntot),100.*timeperbintot(i)*dtimetot
       enddo
       write(iprint,"(a)") '-------------------------------------------------'
    endif
 endif
10 format('|',i3,3x,es10.3,1x,i10,2x,f6.2,'%',2x,f6.2,'% |')

end subroutine write_binsummary

!----------------------------------------------------------------
!+
!  routine to keep track of cpu time spent in each timestep bin
!+
!----------------------------------------------------------------
subroutine update_time_per_bin(dtcpu,istepfrac,nbinmax,timeperbin,inbin)
 real(kind=4),    intent(in)    :: dtcpu
 integer,         intent(in)    :: istepfrac
 integer(kind=1), intent(in)    :: nbinmax
 real(kind=4),    intent(inout) :: timeperbin(0:maxbins)
 integer,         intent(out)   :: inbin
 integer :: i

 !--work out which bin we are updating
 inbin = nbinmax + 1 ! to prevent compiler warnings
 do i=nbinmax,0,-1
    if (mod(istepfrac,2**(nbinmax-i))==0) inbin = i
 enddo

 !--update time counter for this bin
 timeperbin(inbin) = timeperbin(inbin) + dtcpu

end subroutine update_time_per_bin

!-----------------------------------------------------------------
!+
!  routine to print out the timestep information to the log file
!  this version handles individual timesteps
!+
!-----------------------------------------------------------------
subroutine print_dtlog_ind(iprint,ifrac,nfrac,time,dt,nactive,tcpu,np)
 use io, only:formatreal,formatint
 integer,         intent(in) :: iprint,ifrac,nfrac
 real,            intent(in) :: time,dt
 integer(kind=8), intent(in) :: nactive
 real(kind=4),    intent(in) :: tcpu
 integer,         intent(in) :: np
 character(len=120) :: string
 character(len=14) :: tmp
 integer, save :: nplast = 0

 call formatint(ifrac,tmp)
 string = '> step '//tmp
 call formatint(nfrac,tmp)
 string = trim(string)//' / '//tmp
 write(tmp,"(g14.7)") time
 string = trim(string)//' t = '//trim(adjustl(tmp))
 call formatreal(dt,tmp)
 string = trim(string)//' dt = '//tmp
 call formatint(nactive,tmp)
 string = trim(string)//' moved '//tmp
 call formatreal(real(tcpu),tmp)
 string = trim(string)//' in '//trim(tmp)//' cpu-s <'
 !
 ! only print particle number if it differs from
 ! the last time we printed it out
 !
 if (np /= nplast) then
    nplast = np
    write(tmp,"(i12)") np
    string = trim(string)//' | np = '//trim(adjustl(tmp))//' |'
 endif

 write(iprint,"(a)") trim(string)
! write(iprint,5) ifrac,2**nbinmaxprev,time,dt,nactivetot,tcpu2-tcpu1
!5   format('> step ',i6,' /',i6,2x,'t = ',es14.7,1x,'dt = ',es10.3,' moved ',i10,' in ',f8.2,' cpu-s <')

end subroutine print_dtlog_ind

end module timestep_ind
