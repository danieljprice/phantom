!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: io_summary
!
!  DESCRIPTION:
!  Print a summary given parameters every X steps
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module io_summary
 implicit none
 integer, parameter :: maxiosum = 39&         ! Number of values to summarise
#ifdef RADIATION
 +1&
#endif
 +0

 integer, parameter :: maxrhomx = 16         ! Number of maximum possible rhomax' per set
 integer, parameter :: maxisink =  5         ! Maximum number of sink particles's accretion details to track
 !--Array indicies for various parameters
 !  Timesteps
 integer, parameter :: iosumdtf   =  1       ! dtforce (gas particles)
 integer, parameter :: iosumdtfng =  2       ! dtforce (non-gas particles)
 integer, parameter :: iosumdtd   =  3       ! dtdrag (gas)
 integer, parameter :: iosumdtdd  =  4       ! dtdrag (dust)
 integer, parameter :: iosumdtv   =  5       ! dtviscous
 integer, parameter :: iosumdtc   =  6       ! dtcool
 integer, parameter :: iosumdto   =  7       ! dtohmic
 integer, parameter :: iosumdth   =  8       ! dthall
 integer, parameter :: iosumdta   =  9       ! dtambipolar
 integer, parameter :: iosumdte   = 10       ! dtdust

 !  Dust Parameters
 integer, parameter :: iosumdge   = 11       ! Stokes drag regime
 integer, parameter :: iosumdgs   = 12       ! supersonic Epstein regime
 integer, parameter :: iosumdgr   = 13       ! ensuring h < t_s*c_s
 !  Super-timetsepping
 integer, parameter :: iosumstse  = 14       ! STS enabled with N -> Nsts
 integer, parameter :: iosumsts   = 15       ! Nreal for STS enabled with N -> Nsts
 integer, parameter :: iosumstsi  = 16       ! ***Nviaibin for STS enabled with N -> Nsts
 integer, parameter :: iosumstsm  = 17       ! STS enabled with N -> Nsts*Nmega
 integer, parameter :: iosumstsr  = 18       ! Nreal for STS enabled with N -> Nsts*Nmega
 integer, parameter :: iosumstsri = 19       ! ***Nviaibin for STS enabled with N -> Nsts*Nmega
 integer, parameter :: iosumstsd  = 20       ! STS disabled since Nsupersteps=Nreal for small Nreal
 integer, parameter :: iosumstsdi = 21       ! ***Nviaibin for STS disabled since Nsupersteps=Nreal for small Nreal
 integer, parameter :: iosumstso  = 22       ! STS disabled since Nsupersteps=Nreal for large Nreal
 integer, parameter :: iosumstsoi = 23       ! ***Nviaibin for STS disabled since Nsupersteps=Nreal for large Nreal
 integer, parameter :: iosumstsnn = 24       ! number of particles using Nsts
 integer, parameter :: iosumstsnm = 25       ! number of particles using Nsts*Nmega
 integer, parameter :: iosumstsns = 26       ! number of particles using Nsupersteps=Nreal for small Nreal
 integer, parameter :: iosumstsnl = 27       ! number of particles using Nsupersteps=Nreal for large Nreal
 !  Substeps for dtextf < dthydro
 integer, parameter :: iosumextr  = 28       ! ratio due to external force .or. sink-sink
 integer, parameter :: iosumextt  = 29       ! dtmin due to external force .or. sink-sink
 integer, parameter :: iosumexter = 30       ! ratio due to external force
 integer, parameter :: iosumextet = 31       ! dtmin due to external force
 integer, parameter :: iosumextsr = 32       ! ratio due to sink-sink
 integer, parameter :: iosumextst = 33       ! dtmin due to sink-sink
 !  restricted h jump
 integer, parameter :: iosumhup   = 34       ! jump up
 integer, parameter :: iosumhdn   = 35       ! jump down
 !  velocity-dependent force iterations
 integer, parameter :: iosumtvi   = 36       ! number of iterations
 integer, parameter :: iosumtve   = 37       ! errmax
 integer, parameter :: iosumtvv   = 38       ! vmean
 !  particle waking
 integer, parameter :: iowake     = 39       ! number of woken particles
 !  Number of steps
#ifdef RADIATION
 integer, parameter :: iosumdtr   = maxiosum
#endif

 integer, parameter :: iosum_nreal = maxiosum+1  ! number of 'real' steps taken
 integer, parameter :: iosum_nsts  = maxiosum+2  ! number of 'actual' steps (including STS) taken
 !
 !  Frequency of output based number of steps; 0 to turn off; if < 0 then every 2**{-iosum_nprint} steps
 integer,         parameter :: iosum_nprint  = 0
 !  Frequency of output based upon wall time (in seconds); <=0 to turn off
 real(kind=4),    parameter :: dt_wall_print = 0.
 !
 !--Local values and arrays
 integer,           private :: nrhomax, iosum_print, iosum_isink
 integer(kind=8),   private :: iosum_npart(maxiosum)
 integer,           private :: iosum_nstep(maxiosum+2)
 real,              private :: iosum_ave  (maxiosum  ), iosum_max  (maxiosum)
 integer,           private :: iosum_rxi  (maxrhomx  ), iosum_rxp  (maxrhomx), iosum_rxf(7,maxrhomx)
 real,              private :: iosum_rxa  (maxrhomx  ), iosum_rxx  (maxrhomx)
 integer,           private :: accretefail(3)
 logical,           private :: print_dt,print_sts,print_ext,print_dust,print_tolv,print_h,print_wake
 logical,           private :: print_afail,print_early
 real(kind=4),      private :: dtsum_wall
 character(len=19), private :: freason(8)
 !
 !--Public values and arrays
 integer, public  :: iosum_ptmass(5,maxisink)
 logical, public  :: print_acc
 !
contains
!
!----------------------------------------------------------------
!+
!  Initialises values for summary
!  To set iosum_print, let n = iosum_nprint:
!    if n=0 then iosum_print=0:      prints every dump
!    if n>0 then iosum_print=n:      prints every n (real) steps
!    if n<0 then iosum_print=2^{-n}: prints every 2^{-n} (real) steps
!+
!----------------------------------------------------------------
subroutine summary_initialise
 !
 call summary_reset
 !
 if (iosum_nprint < 0) then
    iosum_print = 2**(-iosum_nprint)
 else
    iosum_print = iosum_nprint
 endif
 freason(1) = 'not all active:    '
 freason(2) = 'E_therm/E_grav>0.5:'
 freason(3) = 'a_grav+b_grav > 1: '
 freason(4) = 'E_tot > 0:         '
 freason(5) = 'Not gas:           '
 freason(6) = 'div v > 0:         '
 freason(7) = '2h > h_acc:        '
 freason(8) = '                   '
 !
end subroutine summary_initialise
!----------------------------------------------------------------
!+
!  Resets all values
!+
!----------------------------------------------------------------
subroutine summary_reset
 !
 nrhomax      = 0
 iosum_nstep  = 0
 iosum_npart  = 0
 iosum_ave    = 0.0
 iosum_max    = 0.0
 iosum_rxi    = 0
 iosum_rxp    = 0
 iosum_rxa    = 0.0
 iosum_rxx    = 0.0
 iosum_rxf    = 0
 iosum_ptmass = 0
 accretefail  = 0
 dtsum_wall   = 0.0
 print_dt     = .false.
 print_sts    = .false.
 print_ext    = .false.
 print_dust   = .false.
 print_tolv   = .false.
 print_h      = .false.
 print_acc    = .false.
 print_afail  = .false.
 print_early  = .false.
 print_wake   = .false.
 !
end subroutine summary_reset
!----------------------------------------------------------------
!+
!  Counts number of steps
!+
!----------------------------------------------------------------
subroutine summary_counter(ival,dtstep_wall)
 integer,                intent(in) :: ival
 real(kind=4), optional, intent(in) :: dtstep_wall
 !
 iosum_nstep(ival) = iosum_nstep(ival) + 1
 if (present(dtstep_wall)) dtsum_wall = dtsum_wall + dtstep_wall
 !
end subroutine summary_counter
!----------------------------------------------------------------
!+
!  Adds details to variable tracker
!+
!----------------------------------------------------------------
subroutine summary_variable(cval,ival,nval,meanvalue,maxvalue,addnval)
 integer,           intent(in) :: ival,nval
 real,              intent(in) :: meanvalue
 character(len=*),  intent(in) :: cval
 real,    optional, intent(in) :: maxvalue
 logical, optional, intent(in) :: addnval
 !
 iosum_nstep(ival) = iosum_nstep(ival) + 1
 iosum_npart(ival) = iosum_npart(ival) + nval
 iosum_ave  (ival) = iosum_ave  (ival) + meanvalue
 if (present(maxvalue)) then
    iosum_max  (ival) = max( iosum_max(ival), maxvalue  )
 else
    iosum_max  (ival) = max( iosum_max(ival), meanvalue )
 endif
 ! addval is ONLY used to track limiting timestep for global timestepping
 ! The limiting dt will be called twice per step, hence the subtraction
 ! to undo the second addition
 if (present(addnval)) then
    if (addnval) then
       iosum_npart(ival) = iosum_npart(ival) + 1
       iosum_nstep(ival) = iosum_nstep(ival) - 1
    endif
 endif
 !
 if (trim(cval)=='dt'   ) print_dt   = .true.
 if (trim(cval)=='sts'  ) print_sts  = .true.
 if (trim(cval)=='ext'  ) print_ext  = .true.
 if (trim(cval)=='dust' ) print_dust = .true.
 if (trim(cval)=='tolv' ) print_tolv = .true.
 if (trim(cval)=='hupdn') print_h    = .true.
 if (trim(cval)=='wake' ) print_wake = .true.
 !
end subroutine summary_variable
!----------------------------------------------------------------
!+
!  Adds details to rhomax variable tracker
!+
!----------------------------------------------------------------
subroutine summary_variable_rhomax(ipart,inrho,iprint,nptmass)
 integer,         intent(in)    :: ipart,iprint,nptmass
 real,            intent(in)    :: inrho
 integer                        :: i,j
 !
 !--Determine the index where this particle belongs
 j = -1
 if (nrhomax==0) then
    nrhomax = 1
    j       = nrhomax
 else
    i = 1
    do while (j==-1 .and. i <= nrhomax)
       if (iosum_rxi(i)==ipart) j = i
       i = i + 1
    enddo
    if (j == -1) then
       nrhomax = nrhomax + 1
       j       = nrhomax
    endif
 endif
 !--Reset summary if we have too many particle
 if (j > maxrhomx) then
    nrhomax     = maxrhomx
    print_early = .true.
    call summary_printout(iprint,nptmass)
    nrhomax = 1
    j       = 1
 endif
 iosum_isink = j
 !--Fill in the array
 iosum_rxi(j) = ipart
 iosum_rxp(j) = iosum_rxp(j) + 1
 iosum_rxa(j) = iosum_rxa(j) + inrho
 iosum_rxx(j) = max(iosum_rxx(j),inrho)
 !
end subroutine summary_variable_rhomax
!----------------------------------------------------------------
!+
!  Summaries why a particle failed to turn into a sink particle
!+
!----------------------------------------------------------------
subroutine summary_ptmass_fail(ifail)
 integer,           intent(in) :: ifail
 !
 iosum_rxf(ifail,iosum_isink) = iosum_rxf(ifail,iosum_isink) + 1
 !
end subroutine summary_ptmass_fail
!----------------------------------------------------------------
!+
!  Modifies the iosum_ptmass array after all npart particles are
!  checked for accretion
!+
!----------------------------------------------------------------
subroutine summary_accrete(nptmass)
 integer,           intent(in) :: nptmass
 integer                       :: i, imax
 !
 if (nptmass > maxisink) then
    imax = 1
 else
    imax = nptmass
 endif
 do i = 1,imax
    if (iosum_ptmass(1,i) > 0) then
       iosum_ptmass(3,i) = iosum_ptmass(3,i) + iosum_ptmass(1,i)
       iosum_ptmass(4,i) = max(iosum_ptmass(1,i),iosum_ptmass(4,i))
       iosum_ptmass(5,i) = iosum_ptmass(5,i) + 1
       iosum_ptmass(1,i) = 0
    endif
 enddo
 !
end subroutine summary_accrete
!----------------------------------------------------------------
!+
!  Summaries how many particles failed to accrete
!+
!----------------------------------------------------------------
subroutine summary_accrete_fail(nfail)
 integer,           intent(in) :: nfail
 !
 if (nfail > 0) then
    accretefail(1) = accretefail(1) + 1
    accretefail(2) = accretefail(2) + nfail
    accretefail(3) = max(accretefail(3), nfail)
    print_afail    = .true.
 endif
 !
end subroutine summary_accrete_fail
!----------------------------------------------------------------
!+
!  Determine if it is the time to print the summary,
!  based upon 1) number of steps since printing of last summary
!             2) wall time since printing of last summary
!+
!----------------------------------------------------------------
logical function summary_printnow()
 !
 summary_printnow = .false.
 !
 if ( ( iosum_print /= 0 .and. mod(iosum_nstep(iosum_nreal),iosum_print)==0 ) .or. &
      ( dt_wall_print >  0.0 .and. dtsum_wall > dt_wall_print) )                   &
   summary_printnow = .true.
 !
end function summary_printnow
!----------------------------------------------------------------
!+
!  Print details to screen IF values are non-zero
!+
!----------------------------------------------------------------
subroutine summary_printout(iprint,nptmass)
 integer, intent(in) :: iprint,nptmass
 integer             :: i,j,k,nrealsteps
 real                :: iosum_rpart(maxiosum)
 logical             :: get_averages,print_summary,findfail
 !
 !--do not print if zero steps since last summary
 !  required in case called in a parallel loop prior to a fatal call
 if (iosum_nstep(iosum_nreal)==0) return
 nrealsteps = iosum_nstep(iosum_nreal)
 iosum_nstep(iosum_nreal) = 0
 !
 !--summarise logicals for cleanliness
 !
 if (print_dt .or. print_dust .or. print_sts .or. print_ext .or. print_tolv .or. print_h .or. print_wake) then
    get_averages = .true.
 else
    get_averages = .false.
 endif
 if (get_averages .or. print_acc .or. print_afail .or. nrhomax > 0) then
    print_summary = .true.
 else
    print_summary = .false.
 endif
 !
 !--Get averages
 !
 if ( get_averages ) then
    do i = 1,maxiosum
       if (iosum_nstep(i)/=0) then
          iosum_ave  (i) = iosum_ave(i)/real(iosum_nstep(i))
          iosum_rpart(i) = real(iosum_npart(i))/real(iosum_nstep(i))
       endif
    enddo
 endif
 !
 !--Print Header
 !
 if (print_summary) then
    write(iprint,'(a)') '------------------------------------------------------------------------------'
    if (print_early) write(iprint,'(a)') '|** Printing EARLY since nrhomax > maxrhomx in ptmass                      **|'
#ifdef STS_TIMESTEPS
    write(iprint,10)    '|** Number of (real) steps since last summary:      ',nrealsteps,              '**|'
    write(iprint,10)    '|** Number of steps (incl. STS) since last summary: ',iosum_nstep(iosum_nsts) ,'**|'
#else
    write(iprint,10)    '|** Number of steps since last summary:             ',nrealsteps,              '**|'
#endif
    if ( dtsum_wall > real(3600.0,kind=4) ) then
       write(iprint,20) '|** Wall time since last summary: ',dtsum_wall/real(3600.0,kind=4),' hours       **|'
    else if ( dtsum_wall > real(60.0,kind=4) ) then
       write(iprint,20) '|** Wall time since last summary: ',dtsum_wall/real(60.0,kind=4),  ' minutes     **|'
    else
       write(iprint,20) '|** Wall time since last summary: ',dtsum_wall,                    ' seconds     **|'
    endif
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
10 format(a,i10,13x,a)
20 format(a,18x,F10.3,a)
 !
 !-- Summary of Timesteps
 if ( print_dt ) then
    write(iprint,'(a)') '|* force: timesteps constraining dt compared to Courant by factor of:       *|'
#ifdef IND_TIMESTEPS
    write(iprint,'(a)') '| dt                 |  #times     | mean # part | mean fac    | max fac     |'
    if (iosum_nstep(iosumdtf  )/=0) write(iprint,30) '| dtforce (gas)      |' &
      ,iosum_nstep(iosumdtf  ),'|',iosum_rpart(iosumdtf  ),'|',iosum_ave(iosumdtf  ),'|',iosum_max(iosumdtf  ),'|'
    if (iosum_nstep(iosumdtfng)/=0) write(iprint,30) '| dtforce (non-gas)  |' &
      ,iosum_nstep(iosumdtfng),'|',iosum_rpart(iosumdtfng),'|',iosum_ave(iosumdtfng),'|',iosum_max(iosumdtfng),'|'
    if (iosum_nstep(iosumdtd  )/=0) write(iprint,30) '| dtdrag (gas)       |' &
      ,iosum_nstep(iosumdtd  ),'|',iosum_rpart(iosumdtd  ),'|',iosum_ave(iosumdtd  ),'|',iosum_max(iosumdtd  ),'|'
    if (iosum_nstep(iosumdtdd )/=0) write(iprint,30) '| dtdrag (dust)      |' &
      ,iosum_nstep(iosumdtdd ),'|',iosum_rpart(iosumdtdd ),'|',iosum_ave(iosumdtdd ),'|',iosum_max(iosumdtdd ),'|'
    if (iosum_nstep(iosumdtv  )/=0) write(iprint,30) '| dtvisc             |' &
      ,iosum_nstep(iosumdtv  ),'|',iosum_rpart(iosumdtv  ),'|',iosum_ave(iosumdtv  ),'|',iosum_max(iosumdtv  ),'|'
    if (iosum_nstep(iosumdtc  )/=0) write(iprint,30) '| dtcool             |' &
      ,iosum_nstep(iosumdtc  ),'|',iosum_rpart(iosumdtc  ),'|',iosum_ave(iosumdtc  ),'|',iosum_max(iosumdtc  ),'|'
    if (iosum_nstep(iosumdto  )/=0) write(iprint,30) '| dtohmic            |' &
      ,iosum_nstep(iosumdto  ),'|',iosum_rpart(iosumdto  ),'|',iosum_ave(iosumdto  ),'|',iosum_max(iosumdto  ),'|'
    if (iosum_nstep(iosumdth  )/=0) write(iprint,30) '| dthall             |' &
      ,iosum_nstep(iosumdth  ),'|',iosum_rpart(iosumdth  ),'|',iosum_ave(iosumdth  ),'|',iosum_max(iosumdth  ),'|'
    if (iosum_nstep(iosumdta  )/=0) write(iprint,30) '| dtambi             |' &
      ,iosum_nstep(iosumdta  ),'|',iosum_rpart(iosumdta  ),'|',iosum_ave(iosumdta  ),'|',iosum_max(iosumdta  ),'|'
    if (iosum_nstep(iosumdte  )/=0) write(iprint,30) '| dtdust             |' &
      ,iosum_nstep(iosumdte  ),'|',iosum_rpart(iosumdte  ),'|',iosum_ave(iosumdte  ),'|',iosum_max(iosumdte  ),'|'
#ifdef RADIATION
    if (iosum_nstep(iosumdtr  )/=0) write(iprint,30) '| dtradiation         |' &
      ,iosum_nstep(iosumdtr  ),'|',iosum_rpart(iosumdtr  ),'|',iosum_ave(iosumdtr  ),'|',iosum_max(iosumdtr  ),'|'
#endif
30  format(a,i13,a,3(f13.2,a))
#else
    write(iprint,'(a)') '| dt                    |       < dt_Courant      |     #times smallest      |'
    if (iosum_nstep(iosumdtf)/=0) &
     write(iprint,40)   '| dtforce/dtcool/dtdrag | ',iosum_nstep(iosumdtf),' |',iosum_npart(iosumdtf),'   |'
    if (iosum_nstep(iosumdtv)/=0) &
     write(iprint,40)   '| dtvisc                | ',iosum_nstep(iosumdtv),' |',iosum_npart(iosumdtv),'   |'
    if (iosum_nstep(iosumdto)/=0) &
     write(iprint,40)   '| dtohmic               | ',iosum_nstep(iosumdto),' |',iosum_npart(iosumdto),'   |'
    if (iosum_nstep(iosumdth)/=0) &
     write(iprint,40)   '| dthall                | ',iosum_nstep(iosumdth),' |',iosum_npart(iosumdth),'   |'
    if (iosum_nstep(iosumdta)/=0) &
     write(iprint,40)   '| dtambi                | ',iosum_nstep(iosumdta),' |',iosum_npart(iosumdta),'   |'
#ifdef RADIATION
    if (iosum_nstep(iosumdtr)/=0) &
     write(iprint,40)   '| dtradiation           | ',iosum_nstep(iosumdtr),' |',iosum_npart(iosumdtr),'   |'
#endif
40  format(a,2(i23,a))
#endif
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
 !
 !--Summary of super-timestepping
 if ( print_sts ) then
    write(iprint,'(a)') '|* Super-timestepping                                                       *|'
#ifdef IND_TIMESTEPS
    write(iprint,'(a)') '|         |#     |      npart      |    N_used   |   N_real    | N_via ibin  |'
    write(iprint,'(a)') '|         |times |  mean  |   max  |  mean | max |  mean | max |  mean | max |'
    if (iosum_nstep(iosumstse)/=0) write(iprint,50) &
     '|  Nsts   |',iosum_nstep(iosumstse),'|',iosum_ave(iosumstsnn),'|',int(iosum_max(iosumstsnn)),'|' &
                                             ,iosum_ave(iosumstse ),'|',int(iosum_max(iosumstse )),'|' &
                                             ,iosum_ave(iosumsts  ),'|',int(iosum_max(iosumsts  )),'|' &
                                             ,iosum_ave(iosumstsi ),'|',int(iosum_max(iosumstsi )),'|'
    if (iosum_nstep(iosumstsm)/=0) write(iprint,50) &
     '|NstsNmega|',iosum_nstep(iosumstsm),'|',iosum_ave(iosumstsnm),'|',int(iosum_max(iosumstsnm)),'|' &
                                             ,iosum_ave(iosumstsm ),'|',int(iosum_max(iosumstsm )),'|' &
                                             ,iosum_ave(iosumstsr ),'|',int(iosum_max(iosumstsr )),'|' &
                                             ,iosum_ave(iosumstsri),'|',int(iosum_max(iosumstsri)),'|'
    if (iosum_nstep(iosumstsd)/=0) write(iprint,60) &
     '|Off N=N_s|',iosum_nstep(iosumstsd),'|',iosum_ave(iosumstsns),'|',int(iosum_max(iosumstsns)),'|','|' &
                                             ,iosum_ave(iosumstsd ),'|',int(iosum_max(iosumstsd )),'|' &
                                             ,iosum_ave(iosumstsdi),'|',int(iosum_max(iosumstsdi)),'|'
    if (iosum_nstep(iosumstso)/=0) write(iprint,60) &
     '|Off N=N_b|',iosum_nstep(iosumstso),'|',iosum_ave(iosumstsnl),'|',int(iosum_max(iosumstsnl)),'|','|' &
                                             ,iosum_ave(iosumstso ),'|',int(iosum_max(iosumstso )),'|' &
                                             ,iosum_ave(iosumstsoi),'|',int(iosum_max(iosumstsoi)),'|'
50  format(a,i6,a,f8.1,a,i8,a,      3(f7.1,a,i5,a))
60  format(a,i6,a,f8.1,a,i8,a,13x,a,2(f7.1,a,i5,a))
#else
    write(iprint,'(a)') '|         |#     |       npart       |      N_used       |       N_real      |'
    write(iprint,'(a)') '|         |times |   mean  |   max   |   mean  |   max   |   mean  |   max   |'
    if (iosum_nstep(iosumstse)/=0) write(iprint,70) &
     '|  Nsts   |',iosum_nstep(iosumstse),'|',iosum_ave(iosumstsnn),'|',int(iosum_max(iosumstsnn)),'|' &
                                             ,iosum_ave(iosumstse ),'|',int(iosum_max(iosumstse )),'|' &
                                             ,iosum_ave(iosumsts  ),'|',int(iosum_max(iosumsts  )),'|'
    if (iosum_nstep(iosumstsm)/=0) write(iprint,70) &
     '|NstsNmega|',iosum_nstep(iosumstsm),'|',iosum_ave(iosumstsnm),'|',int(iosum_max(iosumstsnm)),'|' &
                                             ,iosum_ave(iosumstsm ),'|',int(iosum_max(iosumstsm )),'|' &
                                             ,iosum_ave(iosumstsr ),'|',int(iosum_max(iosumstsr )),'|'
    if (iosum_nstep(iosumstsd)/=0) write(iprint,80) &
     '|Off N=N_s|',iosum_nstep(iosumstsd),'|',iosum_ave(iosumstsns),'|',int(iosum_max(iosumstsns)),'|','|' &
                                             ,iosum_ave(iosumstsd ),'|',int(iosum_max(iosumstsd )),'|'
    if (iosum_nstep(iosumstso)/=0) write(iprint,80) &
     '|Off N=N_b|',iosum_nstep(iosumstso),'|',iosum_ave(iosumstsnl),'|',int(iosum_max(iosumstsnl)),'|','|' &
                                             ,iosum_ave(iosumstso ),'|',int(iosum_max(iosumstso )),'|'
70  format(a,i6,a,3(f8.1,1x,a,i8,1x,a))
80  format(a,i6,a,f8.1,1x,a,i8,1x,a,19x,a,f8.1,1x,a,i8,1x,a)
#endif
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
 !
 !--Summary of substepping since dtextf < dthydro
 if ( print_ext ) then
    write(iprint,'(a)') '|* Sub-Steps since dtextf < dthydro                                         *|'
    write(iprint,'(a)') '|         |#     |N_substep|  dthydro/dtextf   |             dt              |'
    write(iprint,'(a)') '|         |times |  mean   |   mean  |   max   |     mean     |     min      |'
    if (iosum_nstep(iosumextr)/=0) write(iprint,90) &
     '|         |',iosum_nstep(iosumextr ),'|',iosum_rpart(iosumextr ),'|' &
                                              ,iosum_ave(iosumextr   ),'|',    iosum_max(iosumextr ) ,'|' &
                                              ,iosum_ave(iosumextt   ),'|',1.0/iosum_max(iosumextt ) ,'|'
    if (iosum_nstep(iosumexter)/=0) write(iprint,90) &
     '|external |',iosum_nstep(iosumexter),'|',iosum_rpart(iosumexter),'|' &
                                              ,iosum_ave(iosumexter),'|',    iosum_max(iosumexter) ,'|' &
                                              ,iosum_ave(iosumextet),'|',1.0/iosum_max(iosumextet) ,'|'
    if (iosum_nstep(iosumextsr)/=0) write(iprint,90) &
     '|sink     |',iosum_nstep(iosumextsr),'|',iosum_rpart(iosumextsr),'|' &
                                              ,iosum_ave(iosumextsr),'|',    iosum_max(iosumextsr) ,'|' &
                                              ,iosum_ave(iosumextst),'|',1.0/iosum_max(iosumextst) ,'|'
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
90 format(a,i6,a,f8.2,1x,a,f8.2,1x,a,f8.2,1x,a,Es13.3,1x,a,Es13.3,1x,a)
 !
 !--Summary of Dust terms
 if ( print_dust ) then
    write(iprint,'(a)') '|* force: Dust warnings                                                     *|'
    if (iosum_nstep(iosumdgs)/=0 .or. iosum_nstep(iosumdge)/=0) &
     write(iprint,'(a)') '| regime  |  #times  | mean # part | mean %super |   max %super              |'
    if (iosum_nstep(iosumdgs)/=0) write(iprint,100) '| Stokes  | '&
      ,iosum_nstep(iosumdgs),'|',iosum_rpart(iosumdgs),'|',iosum_ave(iosumdgs),'|',iosum_max(iosumdgs),'              |'
    if (iosum_nstep(iosumdge)/=0) write(iprint,100) '| Epstein | '&
      ,iosum_nstep(iosumdge),'|',iosum_rpart(iosumdge),'|',iosum_ave(iosumdge),'|',iosum_max(iosumdge),'              |'
    if (iosum_nstep(iosumdgr)/=0) &
     write(iprint,'(a)') '|         |  #times  | mean # part |  mean fac   |   max fac                 |'
    if (iosum_nstep(iosumdgr)/=0) write(iprint,100) &
    '| h>cs*ts | ',iosum_nstep(iosumdgr),'|',iosum_rpart(iosumdgr),'|',iosum_ave(iosumdgr),'|',iosum_max(iosumdgr),'              |'
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
100 format(a,i9,a,3(f13.2,a))
 !
 !--Summary of velocity-dependent force interations (i.e. tolv in step)
 !  (Note: iterations includes the final, converged iteration)
 if (print_tolv) then
    write(iprint,'(a)') '|* velocity-dependent force; errmax & vmean for non-converged iteration     *|'
    write(iprint,'(a)') '|#     |   iterations    |         errmax          |          vmean          |'
    write(iprint,'(a)') '|times |  mean  |   max  |   mean     |    max     |    mean    |    max     |'
    if (iosum_nstep(iosumtvi)/=0) write(iprint,110)    '|' &
      ,iosum_nstep(iosumtvi),' |',iosum_ave(iosumtvi),' |',int(iosum_max(iosumtvi)),' |' &
                                 ,iosum_ave(iosumtve),' |',    iosum_max(iosumtve) ,' |' &
                                 ,iosum_ave(iosumtvv),' |',    iosum_max(iosumtvv) ,' |'
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
110 format(a,i5,a,f7.2,a,i7,a,4(Es11.3,a))
 !
 !--Summary woken particles
 if ( print_wake ) then
    write(iprint,'(a)') '|* particles woken                                                          *|'
    write(iprint,'(a)') '| #steps  | mean # part/step |  max # part/step                              |'
    write(iprint,115) '|',iosum_nstep(iowake),'|',iosum_ave(iowake),'|',int(iosum_max(iowake)),'|'
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
115 format(a,i9,a,f18.2,a,i17,30x,a)
 !
 !--Summary of Restricted h jumps
 if ( print_h ) then
    write(iprint,'(a)') '|* dens: WARNING: restricted h jump                                         *|'
    write(iprint,'(a)') '|direction|  #steps  | mean # part |  max # part                             |'
    if (iosum_nstep(iosumhup)/=0) write(iprint,120) '| up      | ' &
      ,iosum_nstep(iosumhup),'|',iosum_ave(iosumhup),'|',int(iosum_max(iosumhup)),'|'
    if (iosum_nstep(iosumhdn)/=0) write(iprint,120) '| down    | ' &
      ,iosum_nstep(iosumhdn),'|',iosum_ave(iosumhdn),'|',int(iosum_max(iosumhdn)),'|'
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
120 format(a,i9,a,f13.2,a,i10,31x,a)
 !
 !--Summary of rhomax particles
 if (nrhomax > 0) then
    write(iprint,'(a)') '|* ptmass: rhomax values & reasons for not becoming a sink particle         *|'
    write(iprint,'(a)') '| part no |  #times  | # failures |max rho (cgs)|Reason & times failed       |'
    do i = 1,nrhomax
       j = 1
       findfail = .true.
       do while( findfail )
          if ( iosum_rxf(j,i)==0 ) then
             j = j + 1
             if (j==8) findfail = .false.
          else
             findfail = .false.
          endif
       enddo
       if (j == 8) then
          ! only one particle tested and it is accreted
          write(iprint,130) '|',iosum_rxi(i),'|',iosum_rxp(i),'|',sum(iosum_rxf(1:7,i)),'|',iosum_rxx(i),'|','|'
       else
          ! list particle and its initial reason to not create a sink
          write(iprint,140) &
           '|',iosum_rxi(i),'|',iosum_rxp(i),'|',sum(iosum_rxf(1:7,i)),'|',iosum_rxx(i),'|',freason(j),iosum_rxf(j,i),'|'
       endif
       ! for the above particle, list the remaining reasons it failed to form a sink
       if (j < 7) then
          do k = j+1,7
             if (iosum_rxf(k,i)/=0) write(iprint,150) '|','|','|','|','|',freason(k),iosum_rxf(k,i),'|'
          enddo
       endif
    enddo
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
130 format(a,i9,a,i10,a,i12,a,Es13.6,a,28x,a)
140 format(a,i9,a,i10,a,i12,a,Es13.6,2a,i9,a)
150 format(a,9x,a,10x,a,12x,a,13x,2a,i9,a)
 !
 !--Summary of accretion events
 if ( print_acc ) then
    write(iprint,'(a)') '|* ptmass: number of accreted particles                                     *|'
    write(iprint,'(a)') '| sink no |  # steps | max # acc  | total # | # unconditionally              |'
    write(iprint,'(a)') '|         |  w/ acc  | per step   | accreted| accreted                       |'
    if ( nptmass > maxisink ) then
       if (iosum_ptmass(3,1)/=0) write(iprint,160) &
        '| ALL ',nptmass,'|',iosum_ptmass(5,1),'|',iosum_ptmass(4,1),'|',iosum_ptmass(3,1),'|',iosum_ptmass(2,1),'|'
    else
       do i = 1,nptmass
          if (iosum_ptmass(3,i)/=0) write(iprint,170) &
           '|',i,         '|',iosum_ptmass(5,i),'|',iosum_ptmass(4,i),'|',iosum_ptmass(3,i),'|',iosum_ptmass(2,i),'|'
       enddo
    endif
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
 if ( print_afail ) then
    write(iprint,'(a)') '|* ptmass: particles within hacc that failed to accrete                     *|'
    write(iprint,'(a)') '| # sinks | # steps  | max # fails| total #                                  |'
    write(iprint,'(a)') '|         | w/ fails | per step   | fails                                    |'
    write(iprint,180)   '|',nptmass,'|', accretefail(1),'|',accretefail(3),'|',accretefail(2),'|'
    write(iprint,'(a)') '------------------------------------------------------------------------------'
 endif
160 format(a,i4,a,i10,a,i12,a,i9,a,i9,23x,a)
170 format(a,i9,a,i10,a,i12,a,i9,a,i9,23x,a)
180 format(a,I9,a,I10,a,I12,a,i9,33x,a)
 !
 !--Reset all values
 call summary_reset
 !
end subroutine summary_printout
!----------------------------------------------------------------
end module
