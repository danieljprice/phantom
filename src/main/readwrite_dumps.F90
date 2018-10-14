!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: readwrite_dumps
!
!  DESCRIPTION:
!  This module contains all routines related
!  to the data format.
!
!  For Phantom, the format is identical to sphNG
!  (although with fewer arrays dumped)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dump_utils, eos, externalforces, fileutils,
!    gitinfo, initial_params, io, lumin_nsdisc, mpi, mpiutils, options,
!    part, setup_params, sphNGutils, timestep, units
!+
!--------------------------------------------------------------------------
module readwrite_dumps
 use dump_utils, only:lenid,ndatatypes,i_int,i_int1,i_int2,i_int4,i_int8,&
                      i_real,i_real4,i_real8,int1,int2,int1o,int2o,dump_h,lentag
 implicit none
 character(len=80), parameter, public :: &    ! module version
    modid="$Id$"

 public :: write_smalldump,write_fulldump,read_smalldump,read_dump,write_gadgetdump
 public :: get_blocklimits
 logical, public    :: opened_full_dump       ! for use in analysis files if user wishes to skip small dumps
 logical, public    :: dt_read_in             ! to determine if dt has been read in so that ibin & ibinold can be set on restarts
 integer, parameter :: maxphead = 256         ! max items in header
 integer, parameter, public :: is_small_dump = 1978
 integer, parameter, public :: is_not_mhd = 1979

 private

contains
!--------------------------------------------------------------------
!+
!  utility to determine whether to read a particular block
!  in the dump file, in whole or in part.
!  Allows limited changes to number of threads.
!+
!--------------------------------------------------------------------
subroutine get_blocklimits(npartblock,nblocks,nthreads,id,iblock,noffset,npartread)
 use io, only:die,fatal
 integer(kind=8), intent(in)  :: npartblock
 integer,         intent(in)  :: nblocks,nthreads,id,iblock
 integer,         intent(out) :: noffset,npartread
 integer                      :: nblocksperthread,nthreadsperblock
 character(len=15), parameter :: tag = 'get_blocklimits'
!
!--check for errors in input
!
 if (npartblock < 0) call fatal(tag,'block in dump file has npartinblock < 0')
 if (npartblock > huge(npartread)) call fatal(tag,'npart in block exceeds 32 bit limit')
!
!--usual situation: nblocks = nprocessors
!  read whole block if id = iblock
!
 if (nblocks==nthreads) then
    if (id==iblock-1) then
       !--read whole block
       npartread = int(npartblock)
       noffset   = 0
    else
       !--do not read block
       npartread = 0
       noffset   = 0
    endif

 elseif (nblocks > nthreads .and. mod(nblocks,nthreads)==0) then
!
!--if more blocks than processes and nblocks exactly divisible by nthreads,
!  then just read more than one block per thread
!
    nblocksperthread = nblocks/nthreads
    if (id==(iblock-1)/nblocksperthread) then
       npartread = int(npartblock)
       noffset   = 0
    else
       npartread = 0
       noffset   = 0
    endif

 elseif (nthreads > nblocks .and. mod(nthreads,nblocks)==0) then
!
!--if more threads than blocks, and exactly divisible, read fractions of blocks only
!
    nthreadsperblock = nthreads/nblocks
    if (id/nthreadsperblock==iblock-1) then
       npartread = int((npartblock-1)/nthreadsperblock) + 1
       noffset   = mod(id,nthreadsperblock)*npartread

       if (mod(id,nthreadsperblock)==nthreadsperblock-1) then
          !--last thread has remainder for non-exactly divisible numbers of particles
          npartread = int(npartblock) - (nthreadsperblock-1)*npartread
          !--die if we would need to load balance between more than the last processor.
          if (npartread < 0) then
             print*,' npart to read from last block =',npartread
             call fatal(tag,'error assigning npart to last thread')
          endif
       endif
    else
       npartread = 0
       noffset   = 0
    endif
 else
    noffset = 0
    npartread = 0
    print*,' ERROR: rearrangement of ',nblocks,' blocks to ',nthreads,' threads not implemented'
    call die
 endif

end subroutine get_blocklimits

!--------------------------------------------------------------------
!+
!  utility for initialising each thread
!+
!--------------------------------------------------------------------
subroutine start_threadwrite(id,iunit,filename)
#ifdef MPI
 use mpi
 use mpiutils, only:status,mpierr
#endif
 use io, only:fatal,iverbose
 implicit none
 integer, intent(in) :: id, iunit
 character(len=*), intent(in) :: filename
 integer :: nowgo,ierr

 if (iverbose >= 3) print *,id,' : starting write...'
 nowgo = 0
 if (id  >  0) then
#ifdef MPI
    call MPI_RECV(nowgo,1,MPI_INTEGER,id-1,99,MPI_COMM_WORLD,status,mpierr)
#endif
    open(unit=iunit,file=filename,status='old',form='unformatted',position='append',iostat=ierr)
    if (ierr /= 0) then
       call fatal('start_threadwrite','can''t append to dumpfile '//trim(filename))
    else
       if (iverbose >= 3) print*,'thread ',id,': opened file '//trim(filename)
    endif
 endif

 return
end subroutine start_threadwrite

!--------------------------------------------------------------------
!+
!  utility for finalising each thread
!+
!--------------------------------------------------------------------
subroutine end_threadwrite(id)
 use io, only:iverbose
#ifdef MPI
 use mpi
 use mpiutils, only:mpierr
 use io, only:nprocs
#endif
 implicit none
 integer, intent(in) :: id
#ifdef MPI
 integer :: nowgo
#endif

 if (iverbose >= 3) print *,' thread ',id,' : finished write.'
#ifdef MPI
 if (id  <  nprocs-1) then
    nowgo = 1
    call MPI_SEND(nowgo,1,MPI_INTEGER,id+1,99,MPI_COMM_WORLD,mpierr)
 endif
#endif

 return
end subroutine end_threadwrite

!--------------------------------------------------------------------
!+
!  contruct header string based on compile-time options
!  these are for information only (ie. not important for restarting)
!+
!--------------------------------------------------------------------
character(len=lenid) function fileident(firstchar,codestring)
 use part,    only:h2chemistry,mhd,maxBevol,npartoftype,idust,gravity,lightcurve
 use options, only:use_dustfrac
 use dim,     only:use_dustgrowth,phantom_version_string
 use gitinfo, only:gitsha
 character(len=2), intent(in) :: firstchar
 character(len=*), intent(in), optional :: codestring
 character(len=10) :: datestring, timestring
 character(len=30) :: string
!
!--print date and time stamp in file header
!
 call date_and_time(datestring,timestring)
 datestring = datestring(7:8)//'/'//datestring(5:6)//'/'//datestring(1:4)
 timestring = timestring(1:2)//':'//timestring(3:4)//':'//timestring(5:)

 string = ' '
 if (gravity) string = trim(string)//'+grav'
 if (npartoftype(idust) > 0) string = trim(string)//'+dust'
 if (use_dustfrac) string = trim(string)//'+1dust'
 if (h2chemistry) string = trim(string)//'+H2chem'
 if (lightcurve) string = trim(string)//'+lightcurve'
 if (use_dustgrowth) string = trim(string)//'+dustgrowth'

 if (present(codestring)) then
    fileident = firstchar//':'//trim(codestring)//':'//trim(phantom_version_string)//':'//gitsha
 else
    fileident = firstchar//':Phantom'//':'//trim(phantom_version_string)//':'//gitsha
 endif

 if (mhd) then
    if (maxBevol==4) then
       fileident = trim(fileident)//' (mhd+clean'//trim(string)//')  : '//trim(datestring)//' '//trim(timestring)
    else
       fileident = trim(fileident)//' (mhd'//trim(string)//')  : '//trim(datestring)//' '//trim(timestring)
    endif
 else
    fileident = trim(fileident)//' (hydro'//trim(string)//'): '//trim(datestring)//' '//trim(timestring)
 endif

end function fileident

!--------------------------------------------------------------------
!+
!  extract various options used in Phantom from the fileid string
!+
!--------------------------------------------------------------------
subroutine get_options_from_fileid(fileid,tagged,phantomdump,smalldump,&
                                   use_onefluiddust,ierr)
 character(len=lenid), intent(in)  :: fileid
 logical,              intent(out) :: tagged,phantomdump,smalldump,use_onefluiddust
 integer,              intent(out) :: ierr
!
!--if file is a small dump, return an error code but still read what
!  can be read from a small dump
!
 ierr = 0
 tagged      = .false.
 smalldump   = .false.
 phantomdump = .false.
 if (fileid(2:2)=='T') tagged = .true.
 if (fileid(1:1) /= 'F') then
    !write(*,*) 'ERROR! file header indicates file is not a full dump'
    ierr = 1
    if (fileid(1:1)=='S') smalldump = .true.
 endif
 if (index(fileid,'Phantom') /= 0) then
    phantomdump = .true.
 elseif (index(fileid,'sphNG') /= 0) then
    phantomdump = .false.
    write(*,*) 'reading dump in sphNG format'
 else
    write(*,*) 'WARNING: could not determine Phantom/sphNG from fileident'
    write(*,*) '(assuming sphNG...)'
    phantomdump = .false.
 endif
 if (index(fileid,'+1dust') /= 0) then
    use_onefluiddust = .true.
 else
    use_onefluiddust = .false.
 endif

end subroutine get_options_from_fileid

!--------------------------------------------------------------------
!+
!  extract dump size used in Phantom from the fileid string
!+
!--------------------------------------------------------------------
subroutine get_dump_size(fileid,smalldump)
 character(len=lenid), intent(in)  :: fileid
 logical,              intent(out) :: smalldump
 !
 if (fileid(1:1)=='S') then
    smalldump = .true.
 else
    smalldump = .false.
 endif

end subroutine get_dump_size


!--------------------------------------------------------------------
!+
!  subroutine to write output to full dump file
!  (this is everything needed to restart a run)
!+
!-------------------------------------------------------------------
subroutine write_fulldump(t,dumpfile,ntotal,iorder,sphNG)
 use dim,   only:maxp,maxvxyzu,maxalpha,ndivcurlv,ndivcurlB,maxgrav,gravity,use_dust,&
                 lightcurve,maxlum,store_temperature,use_dustgrowth
 use eos,   only:utherm,ieos,equationofstate,done_init_eos,init_eos
 use io,    only:idump,iprint,real4,id,master,error,warning,nprocs
 use part,  only:xyzh,xyzh_label,vxyzu,vxyzu_label,Bevol,Bxyz,Bxyz_label,npart,npartoftype,maxtypes, &
                 alphaind,rhoh,divBsymm,maxphase,iphase,iamtype_int1,iamtype_int11, &
                 nptmass,nsinkproperties,xyzmh_ptmass,xyzmh_ptmass_label,vxyz_ptmass,vxyz_ptmass_label,&
                 maxptmass,get_pmass,h2chemistry,nabundances,abundance,abundance_label,mhd,maxBevol,&
                 divcurlv,divcurlv_label,divcurlB,divcurlB_label,poten,dustfrac,deltav,deltav_label,tstop,&
                 dustfrac_label,tstop_label,dustprop,dustprop_label,temperature,St,ndusttypes,ndustsmall
 use options,    only:use_dustfrac
 use dump_utils, only:tag,open_dumpfile_w,allocate_header,&
                 free_header,write_header,write_array,write_block_header
 use mpiutils,   only:reduce_mpi,reduceall_mpi
#ifdef IND_TIMESTEPS
 use timestep,   only:dtmax
 use part,       only:ibin
#endif
#ifdef PRDRAG
 use lumin_nsdisc, only:beta
#endif
#ifdef LIGHTCURVE
 use part,  only:luminosity
#endif
#ifdef NONIDEALMHD
 use dim,  only:mhd_nonideal
 use part, only:eta_nimhd,eta_nimhd_label
#endif
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in), optional :: iorder(:)
 logical,          intent(in), optional :: sphNG
 integer(kind=8),  intent(in), optional :: ntotal

 integer, parameter :: isteps_sphNG = 0, iphase0 = 0
 integer(kind=8)    :: ilen(4)
 integer            :: nums(ndatatypes,4)
 integer            :: i,ipass,k,l
 integer            :: ierr,ierrs(20)
 integer            :: nblocks,nblockarrays,narraylengths
 integer(kind=8)    :: nparttot,npartoftypetot(maxtypes)
 logical            :: sphNGdump, write_itype, use_gas
 character(len=lenid)  :: fileid
 type(dump_h)          :: hdr
 real, allocatable :: temparr(:)
 real :: ponrhoi,rhoi,spsoundi

!
!--collect global information from MPI threads
!
!--allow non-MPI calls to create MPI dump files
#ifdef MPI
 nparttot = reduceall_mpi('+',npart)
 npartoftypetot = reduceall_mpi('+',npartoftype)
#else
 if (present(ntotal)) then
    nparttot = ntotal
    npartoftypetot = npartoftype
    if (all(npartoftypetot==0)) then
       npartoftypetot(1) = ntotal
    endif
 else
    nparttot = npart
    npartoftypetot = npartoftype
 endif
#endif
 nblocks = nprocs

 sphNGdump = .false.
 if (present(sphNG)) then
    sphNGdump = sphNG
    if (sphNG) write(iprint,*) 'ERROR: sphNG output no longer supported'
 endif

 fileid = fileident('FT','Phantom')

 if (maxphase==maxp) then
    use_gas = .false.
 else
    use_gas = .true.
 endif

!--number of blocks per thread : for hydro these are hydro + point masses
!  block 3 is blank (rt in sphNG), block 4 is for mhd.
!
 narraylengths = 2
 if (mhd) narraylengths = 4
!
!--open dumpfile
!
 masterthread: if (id==master) then

    write(iprint,"(/,/,'-------->   TIME = ',g12.4,"// &
              "': full dump written to file ',a,'   <--------',/)")  t,trim(dumpfile)

    call open_dumpfile_w(idump,dumpfile,fileid,ierr)
    if (ierr /= 0) then
       call error('write_fulldump','error creating new dumpfile '//trim(dumpfile))
       return
    endif
!
!--single values
!
    hdr = allocate_header(nint=maxphead,nreal=maxphead,err=ierr)

    call fill_header(sphNGdump,t,nparttot,npartoftypetot,nblocks,nptmass,hdr,ierr)
    if (ierr /= 0) call warning('write_fulldump','error filling header arrays')

    call write_header(idump,hdr,ierr)
    if (ierr /= 0) call error('write_fulldump','error writing header to dumpfile')

    call free_header(hdr,ierr)
    if (ierr /= 0) call error('write_fulldump','error deallocating header')
!
!--arrays
!
!--total number of blocks
!  each thread has up to 4 blocks (hydro variables, sink particles, radiative transfer and MHD)
!  repeated nblocks times (once for each MPI process)
!
    nblockarrays = narraylengths*nblocks
    write (idump, iostat=ierr) nblockarrays

 endif masterthread

 call start_threadwrite(id,idump,dumpfile)

 ierrs = 0
 nums = 0
 ilen = 0_8
 if (sphNGdump) then
    write_itype = .true.
 else
    write_itype = any(npartoftypetot(2:) > 0)
 endif
 do ipass=1,2
    do k=1,ndatatypes
       !
       ! Block 1 arrays (hydrodynamics)
       !
       ilen(1) = int(npart,kind=8)
       if (write_itype) call write_array(1,iphase,'itype',npart,k,ipass,idump,nums,ierrs(1),func=iamtype_int11)
       call write_array(1,xyzh,xyzh_label,3,npart,k,ipass,idump,nums,ierrs(2))
       if (use_dustgrowth) then
          call write_array(1,dustprop,dustprop_label,3,npart,k,ipass,idump,nums,ierrs(3))
          call write_array(1,St,'St',npart,k,ipass,idump,nums,ierrs(3))
       endif
       call write_array(1,vxyzu,vxyzu_label,maxvxyzu,npart,k,ipass,idump,nums,ierrs(4))
       if (h2chemistry)  call write_array(1,abundance,abundance_label,nabundances,npart,k,ipass,idump,nums,ierrs(5))
       if (use_dust) &
          call write_array(1,dustfrac,dustfrac_label,ndusttypes,npart,k,ipass,idump,nums,ierrs(7))
       if (use_dust) call write_array(1,tstop,tstop_label,ndustsmall,npart,k,ipass,idump,nums,ierrs(8))
       if (use_dustfrac) then
          do l=1,ndustsmall
             call write_array(1,deltav(:,l,:),deltav_label,3,npart,k,ipass,idump,nums,ierrs(10))
          enddo
       endif
       if (store_temperature) call write_array(1,temperature,'T',npart,k,ipass,idump,nums,ierrs(12))

       ! write pressure to file
       if ((ieos==8 .or. ieos==9 .or. ieos==10 .or. ieos==15) .and. k==i_real) then
          if (.not. allocated(temparr)) allocate(temparr(npart))
          if (.not.done_init_eos) call init_eos(ieos,ierr)
          do i=1,npart
             rhoi = rhoh(xyzh(4,i),get_pmass(i,use_gas))
             if (maxvxyzu >=4 ) then
                if (store_temperature) then
                   ! cases where the eos stores temperature (ie Helmholtz)
                   call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i),temperature(i))
                else
                   call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xyzh(1,i),xyzh(2,i),xyzh(3,i),vxyzu(4,i))
                endif
             else
                call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xyzh(1,i),xyzh(2,i),xyzh(3,i))
             endif
             temparr(i) = ponrhoi*rhoi
          enddo
          call write_array(1,temparr,'pressure',npart,k,ipass,idump,nums,ierrs(13))
       endif

       ! smoothing length written as real*4 to save disk space
       call write_array(1,xyzh,xyzh_label,1,npart,k,ipass,idump,nums,ierrs(14),use_kind=4,index=4)
       if (maxalpha==maxp) call write_array(1,alphaind,(/'alpha'/),1,npart,k,ipass,idump,nums,ierrs(15))
       !if (maxalpha==maxp) then ! (uncomment this to write alphaloc to the full dumps)
       !   call write_array(1,alphaind,(/'alpha ','alphaloc'/),2,npart,k,ipass,idump,nums,ierrs(10))
       !endif
       if (ndivcurlv >= 1) call write_array(1,divcurlv,divcurlv_label,ndivcurlv,npart,k,ipass,idump,nums,ierrs(16))
       if (gravity .and. maxgrav==maxp) then
          call write_array(1,poten,'poten',npart,k,ipass,idump,nums,ierrs(17))
       endif
#ifdef IND_TIMESTEPS
       call write_array(1,dtmax/2**ibin(1:npart),'dt',npart,k,ipass,idump,nums,ierrs(18),use_kind=4)
#endif
#ifdef PRDRAG
       if (k==i_real) then
          if (.not.allocated(temparr)) allocate(temparr(npart))
          do i=1,npart
             temparr(i) = real4(beta(xyzh(1,i), xyzh(2,i), xyzh(3,i)))
          enddo
          call write_array(1,temparr,'beta_pr',npart,k,ipass,idump,nums,ierrs(19))
       endif
#endif
#ifdef LIGHTCURVE
       if (lightcurve) then
          call write_array(1,luminosity,'luminosity',npart,k,ipass,idump,nums,ierrs(20))
       endif
#endif
       if (any(ierrs(1:20) /= 0)) call error('write_dump','error writing hydro arrays')
    enddo

    do k=1,ndatatypes
       !
       ! Block 2 arrays (sink particles)
       !
       if (.not. sphNGdump .and. nptmass > 0 .and. nptmass <= maxptmass) then
          ilen(2) = int(nptmass,kind=8)
          call write_array(2,xyzmh_ptmass,xyzmh_ptmass_label,nsinkproperties,nptmass,k,ipass,idump,nums,ierrs(1))
          call write_array(2,vxyz_ptmass,vxyz_ptmass_label,3,nptmass,k,ipass,idump,nums,ierrs(2))
          if (any(ierrs(1:2) /= 0)) call error('write_dump','error writing sink particle arrays')
       endif
    enddo

    do k=1,ndatatypes
       !
       ! Block 4 arrays (MHD)
       !
       if (mhd) then
          ilen(4) = int(npart,kind=8)
          call write_array(4,Bxyz,Bxyz_label,3,npart,k,ipass,idump,nums,ierrs(1))
          if (maxBevol >= 4) then
             call write_array(4,Bevol(4,:),'psi',npart,k,ipass,idump,nums,ierrs(1))
          endif
          if (ndivcurlB >= 1) then
             call write_array(4,divcurlB,divcurlB_label,ndivcurlB,npart,k,ipass,idump,nums,ierrs(2))
          else
             call write_array(4,divBsymm,'divBsymm',npart,k,ipass,idump,nums,ierrs(2))
          endif
          if (any(ierrs(1:2) /= 0)) call error('write_dump','error writing MHD arrays')
#ifdef NONIDEALMHD
          if (mhd_nonideal) then
             call write_array(4,eta_nimhd,eta_nimhd_label,4,npart,k,ipass,idump,nums,ierrs(1))
             if (ierrs(1) /= 0) call error('write_dump','error writing non-ideal MHD arrays')
          endif
#endif
       endif
    enddo

    if (ipass==1) call write_block_header(narraylengths,ilen,nums,idump,ierr)
 enddo
 if (allocated(temparr)) deallocate(temparr)

 if (ierr /= 0) write(iprint,*) 'error whilst writing dumpfile '//trim(dumpfile)

 close(unit=idump)
 call end_threadwrite(id)

end subroutine write_fulldump

!--------------------------------------------------------------------
!+
!  subroutine to write output to small dump file
!  (ie. minimal output...)
!
!  note that small dumps are always SINGLE PRECISION
!  (faked to look like the default real is real*4)
!+
!-------------------------------------------------------------------

subroutine write_smalldump(t,dumpfile)
 use dim,        only:maxp,maxtypes,use_dust,lightcurve,use_dustgrowth
 use io,         only:idump,iprint,real4,id,master,error,warning,nprocs
 use part,       only:xyzh,xyzh_label,npart,npartoftype,Bxyz,Bxyz_label,&
                      maxphase,iphase,h2chemistry,nabundances,&
                      nptmass,nsinkproperties,xyzmh_ptmass,xyzmh_ptmass_label,&
                      abundance,abundance_label,mhd,dustfrac,iamtype_int11,&
                      dustprop,dustprop_label,dustfrac_label,St,ndusttypes
 use dump_utils, only:open_dumpfile_w,dump_h,allocate_header,free_header,&
                      write_header,write_array,write_block_header
 use mpiutils,   only:reduceall_mpi
#ifdef LIGHTCURVE
 use part,       only:luminosity
#endif
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer(kind=8) :: ilen(4)
 integer         :: nums(ndatatypes,4)
 integer         :: ierr,ipass,k
 integer         :: nblocks,nblockarrays,narraylengths
 integer(kind=8) :: nparttot,npartoftypetot(maxtypes)
 logical         :: write_itype
 type(dump_h)    :: hdr

!
!--collect global information from MPI threads
!
 nparttot = reduceall_mpi('+',npart)
 npartoftypetot = reduceall_mpi('+',npartoftype)
 nblocks = nprocs

 narraylengths = 2
 if (mhd) narraylengths = 4

 masterthread: if (id==master) then
!
!--open dumpfile
!
    write(iprint,"(/,/,'-------->   TIME = ',g12.4,"// &
              "': small dump written to file ',a,'   <--------',/)")  t,trim(dumpfile)

    call open_dumpfile_w(idump,dumpfile,fileident('ST'),ierr,singleprec=.true.)
    if (ierr /= 0) then
       call error('write_smalldump','can''t create new dumpfile '//trim(dumpfile))
       return
    endif
!
!--single values
!
    hdr = allocate_header(nint=maxphead,nreal=maxphead,err=ierr)

    call fill_header(.false.,t,nparttot,npartoftypetot,nblocks,nptmass,hdr,ierr)
    if (ierr /= 0) call warning('write_smalldump','error filling header arrays')

    call write_header(idump,hdr,ierr,singleprec=.true.)
    if (ierr /= 0) call error('write_smalldump','error writing header to dumpfile')

    call free_header(hdr,ierr)
    if (ierr /= 0) call error('write_smalldump','error deallocating header')
!
!--arrays: number of array lengths
!
    nblockarrays = narraylengths*nblocks
    write (idump, iostat=ierr) nblockarrays
    if (ierr /= 0) call error('write_smalldump','error writing nblockarrays')

 endif masterthread

 call start_threadwrite(id,idump,dumpfile)

 nums = 0
 ilen = 0_8
 write_itype = (maxphase==maxp .and. any(npartoftypetot(2:) > 0))
 do ipass=1,2
    do k=1,ndatatypes
       !
       !--Block 1 (hydrodynamics)
       !
       ilen(1) = npart
       if (write_itype) call write_array(1,iphase,'itype',npart,k,ipass,idump,nums,ierr,func=iamtype_int11)
       call write_array(1,xyzh,xyzh_label,3,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       if (use_dustgrowth) then
          call write_array(1,dustprop,dustprop_label,3,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
          call write_array(1,St,'St',npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       endif
       if (h2chemistry .and. nabundances >= 1) &
          call write_array(1,abundance,abundance_label,1,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       if (use_dust) &
          call write_array(1,dustfrac,dustfrac_label,ndusttypes,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       call write_array(1,xyzh,xyzh_label,4,npart,k,ipass,idump,nums,ierr,index=4,use_kind=4)
#ifdef LIGHTCURVE
       if (lightcurve) call write_array(1,luminosity,'luminosity',npart,k,ipass,idump,nums,ierr,singleprec=.true.)
#endif
    enddo
    !
    !--Block 2 (sinks)
    !
    if (nptmass > 0) then
       ilen(2) = nptmass
       call write_array(2,xyzmh_ptmass,xyzmh_ptmass_label,nsinkproperties,nptmass,&
                        i_real,ipass,idump,nums,ierr,singleprec=.true.)
    endif
    !
    !--Block 4 (MHD)
    !
    if (mhd) then
       ilen(4) = npart
       do k=1,ndatatypes
          call write_array(4,Bxyz,Bxyz_label,3,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       enddo
    endif

    if (ipass==1) call write_block_header(narraylengths,ilen,nums,idump,ierr)
 enddo

 close(unit=idump)
 call end_threadwrite(id)
 return

end subroutine write_smalldump

!--------------------------------------------------------------------
!+
!  subroutine to read dump from file
!  needs to be able to read Phantom dumps as in write_fulldump
!  and also from standard sphNG dump files
!+
!-------------------------------------------------------------------

subroutine read_dump(dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc)
 use dim,      only:maxp,maxvxyzu,maxalpha,maxgrav,gravity,lightcurve,maxlum,mhd
 use io,       only:real4,master,iverbose,error,warning ! do not allow calls to fatal in this routine
 use part,     only:xyzh,vxyzu,massoftype,npart,npartoftype,maxtypes,iphase, &
                    maxphase,isetphase,nptmass,nsinkproperties,maxptmass,get_pmass, &
                    xyzmh_ptmass,vxyz_ptmass
 use dump_utils, only:skipblock,skip_arrays,check_tag,lenid,ndatatypes,read_header,read_array, &
                      open_dumpfile_r,get_error_text,ierr_realsize,free_header,read_block_header
 use mpiutils,   only:reduce_mpi,reduceall_mpi
 use sphNGutils, only:convert_sinks_sphNG
 use options,    only:use_dustfrac
 character(len=*),  intent(in)  :: dumpfile
 real,              intent(out) :: tfile,hfactfile
 integer,           intent(in)  :: idisk1,iprint,id,nprocs
 integer,           intent(out) :: ierr
 logical, optional, intent(in)  :: headeronly
 logical, optional, intent(in)  :: dustydisc

 integer               :: number
 integer               :: iblock,nblocks,i1,i2,noffset,npartread,narraylengths
 integer(kind=8)       :: ilen(4)
 integer               :: nums(ndatatypes,4)
 integer(kind=8)       :: nparttot,nhydrothisblock,npartoftypetot(maxtypes),npartoftypetotact(maxtypes)
 logical               :: tagged,phantomdump,smalldump,extradust
 real                  :: dumr,alphafile
 character(len=lenid)  :: fileidentr
 type(dump_h)          :: hdr
 integer               :: i,ierrh

 if (id==master) write(iprint,"(/,1x,a,i3)") '>>> reading setup from file: '//trim(dumpfile)//' on unit ',idisk1
 opened_full_dump = .true.
 dt_read_in       = .false.
 !
 ! open dump file
 !
 call open_dumpfile_r(idisk1,dumpfile,fileidentr,ierr)

 !
 ! exit with error if file not readable by current routine
 !
 if (ierr /= 0) then
    call get_dump_size(fileidentr,smalldump)
    if (smalldump) then
       call error('read_dump','file is not a Phantom full dump')
       ierr = is_small_dump
       close(idisk1)
       return
    else
       call error('read_dump',get_error_text(ierr))
    endif
    !
    ! give helpful hint if precision of file is wrong
    !
    if (ierr == ierr_realsize) then
       select case(kind(dumr))
       case(4)
          write (*,"(a,/)") '   re-compile with DOUBLEPRECISION=yes'
       case(8)
          write (*,"(a,/)") '   re-compile with DOUBLEPRECISION=no'
       end select
    endif
    close(idisk1)
    return
 endif
 if (id==master) write(iprint,*) trim(fileidentr)

 ! extract file type from the fileid string
 call get_options_from_fileid(fileidentr,tagged,phantomdump,smalldump,use_dustfrac,ierr)

 !
 ! read header from the dump file
 !
 call read_header(idisk1,hdr,tagged,ierr)
 if (ierr /= 0) then
    call error('read_dump','error reading header from file')
    return
 endif
 !
 ! for backwards compatibility with old phantom files
 ! fake the tags as if they had been read from the file
 !
 if (.not.tagged) call fake_header_tags(hdr,phantomdump,mhd,maxtypes)

 call unfill_header(hdr,phantomdump,tagged,nparttot, &
                    nblocks,npart,npartoftype, &
                    tfile,hfactfile,alphafile,iprint,id,nprocs,ierr)
 if (ierr /= 0) then
    call error('read_dump','error extracting necessary information from file header')
    call free_header(hdr,ierrh)
    return
 endif

 call free_header(hdr,ierr)
!
!--arrays
!
!--number of array lengths
!
 read (idisk1, end=100) number
 narraylengths = number/nblocks
 if (iverbose >= 2 .and. id==master) then
    write(iprint,"(a,i3)") ' number of array sizes = ',narraylengths
    write(iprint,"(a,i3)") ' number of blocks      = ',nblocks
 endif
!
!--check this
!
 if (mhd .and. narraylengths < 4) then
    write (*,*) 'WARNING! readdump: MHD data not present in dumpfile'
    !ierr = 7
    !return
 elseif (narraylengths < 2 .or. narraylengths > 4) then
    write (*,*) 'error 7 in readdump, narraylengths=',narraylengths
    ierr = 7
    return
 endif

 npart = 0
 i2 = 0

 overblocks: do iblock=1,nblocks
! print*,' thread ',id,' block ',iblock
    nums = 0
    call read_block_header(narraylengths,ilen,nums,idisk1,ierr)
!
!--check block header for errors
!
    call check_block_header(narraylengths,nblocks,ilen,nums,nparttot,nhydrothisblock,nptmass,ierr)
    if (ierr /= 0) then
       call error('read_dump','error in array headers')
       return
    endif
!
!--exit after reading the file header if the optional argument
!  "headeronly" is present and set to true
!
    if (present(headeronly)) then
       if (headeronly) return
    endif
!
!--determine if extra dust quantites should be read
!
    if (present(dustydisc)) then
       extradust = .true.
    else
       extradust = .false.
    endif
!
!--determine whether or not to read this particular block
!  onto this particular thread, either in whole or in part
!  Also handles MPI -> non-MPI dump conversion and vice-versa.
!  Can be used by non-MPI codes to read isolated blocks only.
!
    call get_blocklimits(nhydrothisblock,nblocks,nprocs,id,iblock,noffset,npartread)
    i1 = i2 + 1
    i2 = i1 + (npartread - 1)
    npart = npart + npartread
#ifdef MPI
    if (npart > maxp) then
       write(*,*) 'npart > maxp in readwrite_dumps'
       ierr = 1
       return
    endif
#endif
    if (npartread <= 0 .and. nptmass <= 0) then
       print*,' SKIPPING BLOCK npartread = ',npartread
       call skipblock(idisk1,nums(:,1),nums(:,2),nums(:,3),nums(:,4),tagged,ierr)
       if (ierr /= 0) then
          print*,' error skipping block'
          return
       endif
       cycle overblocks
    elseif (npartread > 0) then
#ifdef MPI
       write(*,"(a,i5,2(a,i10),a,i5,a,i10,'-',i10)") &
     'thread ',id,' reading particles ',noffset+1,':',noffset+npartread,', from block ',iblock,' lims=',i1,i2
#else
       write(*,"(2(a,i10),a,i5,a,i10,'-',i10)") &
     ' reading particles ',noffset+1,':',noffset+npartread,', from block ',iblock,' lims=',i1,i2
#endif
    else
       write(*,"(a,i10,a)") ' WARNING! block contains no SPH particles, reading ',nptmass,' point mass particles only'
    endif

    call read_phantom_arrays(i1,i2,noffset,narraylengths,nums,npartread,npartoftype,&
                          massoftype,nptmass,nsinkproperties,phantomdump,tagged,.false.,&
                          extradust,tfile,alphafile,idisk1,iprint,ierr)

    if (ierr /= 0) call warning('read_dump','error reading arrays from file')

 enddo overblocks

 !
 ! determine npartoftype
 !
 npartoftypetot = npartoftype
 call count_particle_types(npartoftype)
 npartoftypetotact = reduceall_mpi('+',npartoftype)
 do i = 1,maxtypes
    if (npartoftypetotact(i) /= npartoftypetot(i)) then
       write(*,*) 'npartoftypetot    =',npartoftypetot
       write(*,*) 'npartoftypetotact =',npartoftypetotact
       call error('read_dump','particle type counts do not match header')
       ierr = 8
    endif
 enddo

 !
 ! convert sinks from sphNG -> Phantom
 !
 if (.not.phantomdump .and. nptmass > 0 .and. maxphase==maxp) then
    call convert_sinks_sphNG(npart,nptmass,iphase,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,ierr)
 endif

 if (sum(npartoftype)==0) npartoftype(1) = npart
 if (narraylengths >= 4) then
    if (id==master) write(iprint,"(a,/)") ' <<< finished reading (MHD) file '
 else
    if (id==master) write(iprint,"(a,/)") ' <<< finished reading (hydro) file '
 endif
 close(idisk1)
 return

100 close (idisk1)
 if (sum(npartoftype)==0) npartoftype(1) = npart
 write(iprint,"(a,/)") ' <<< ERROR! end of file reached in data read'
 ierr = 666
 return

end subroutine read_dump

!--------------------------------------------------------------------
!+
!  subroutine to read a small dump from file, as written
!  in write_smalldump
!+
!-------------------------------------------------------------------

subroutine read_smalldump(dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc)
 use dim,      only:maxp,maxvxyzu,mhd,maxBevol
 use io,       only:real4,master,iverbose,error,warning ! do not allow calls to fatal in this routine
 use part,     only:npart,npartoftype,maxtypes,nptmass,nsinkproperties,maxptmass, &
                    massoftype
 use dump_utils, only:skipblock,skip_arrays,check_tag,open_dumpfile_r,get_error_text,&
                      ierr_realsize,read_header,extract,free_header,read_block_header
 use mpiutils,   only:reduce_mpi,reduceall_mpi
 use options,    only:use_dustfrac
 character(len=*),  intent(in)  :: dumpfile
 real,              intent(out) :: tfile,hfactfile
 integer,           intent(in)  :: idisk1,iprint,id,nprocs
 integer,           intent(out) :: ierr
 logical, optional, intent(in)  :: headeronly
 logical, optional, intent(in)  :: dustydisc

 integer               :: number
 integer               :: iblock,nblocks,i1,i2,noffset,npartread,narraylengths
 integer(kind=8)       :: ilen(4)
 integer               :: nums(ndatatypes,4)
 integer(kind=8)       :: nparttot,nhydrothisblock,npartoftypetot(maxtypes),npartoftypetotact(maxtypes)
 logical               :: tagged,phantomdump,smalldump,extradust
 real                  :: alphafile
 character(len=lenid)  :: fileidentr
 type(dump_h)          :: hdr
 integer               :: i

 if (id==master) write(iprint,"(/,1x,a,i3)") '>>> reading small dump file: '//trim(dumpfile)//' on unit ',idisk1
 opened_full_dump = .false.
 !
 ! open dump file
 !
 call open_dumpfile_r(idisk1,dumpfile,fileidentr,ierr,singleprec=.true.)
 if (ierr /= 0) then
    call error('read_smalldump',get_error_text(ierr))
    if (ierr == ierr_realsize) then
       if (id==master) write(*,*) ' *** this file appears to be a small dump written in double precision ***'
    endif
    ierr = 1
    return
 endif
 if (id==master) write(iprint,*) trim(fileidentr)

 ! extract file type from the fileid string
 call get_options_from_fileid(fileidentr,tagged,phantomdump,smalldump,use_dustfrac,ierr)

 if (.not.smalldump) then
    if (id==master) call error('read_smalldump','this routine only works for small dump files, aborting...')
    ierr = 2
    return
 else
    if (id==master) call warning('read_smalldump','*** VELOCITY WILL BE MISSING FROM SMALL DUMP FILES ***')
 endif
 if (.not.phantomdump) then
    if (id==master) call error('read_smalldump','this routine only works for phantom small dump files, aborting...')
    ierr = 3
    return
 endif

!
!--single values
!
 call read_header(idisk1,hdr,tagged,ierr,singleprec=.true.)
 if (ierr /= 0) then
    call error('read_smalldump','error reading header from file')
    return
 endif
 !
 ! for backwards compatibility with old phantom files
 ! fake the tags as if they had been read from the file
 !
 if (.not.tagged) call fake_header_tags(hdr,phantomdump,mhd,maxtypes)

 call unfill_header(hdr,phantomdump,tagged,nparttot, &
                    nblocks,npart,npartoftype, &
                    tfile,hfactfile,alphafile,iprint,id,nprocs,ierr)
 if (ierr /= 0) then
    call error('read_smalldump','error extracting header information')
    call free_header(hdr,ierr)
    return
 endif

 call free_header(hdr,ierr)
!
!--arrays
!
!--number of array lengths
!
 read (idisk1, end=100) number
 narraylengths = number/nblocks
 if (iverbose >= 2 .and. id==master) then
    write(iprint,"(a,i3)") ' number of array sizes = ',narraylengths
    write(iprint,"(a,i3)") ' number of blocks      = ',nblocks
 endif
!
!--check this
!
 if (mhd .and. narraylengths < 4) then
    if (id==master) write (*,*) 'WARNING! readdump: MHD data not present in dumpfile'
    !ierr = 7
    !return
 elseif (narraylengths < 2 .or. narraylengths > 4) then
    if (id==master) write (*,*) 'error 7 in readdump, narraylengths=',narraylengths
    ierr = 7
    return
 endif

 npart = 0
 i2 = 0
 nums = 0

 overblocks: do iblock=1,nblocks
    nums = 0
    call read_block_header(narraylengths,ilen,nums,idisk1,ierr)
!
!--check block header for errors
!
    call check_block_header(narraylengths,nblocks,ilen,nums,nparttot,nhydrothisblock,nptmass,ierr)
    if (ierr /= 0) then
       call error('read_dump','error in array headers')
       return
    endif
!
!--exit after reading the file header if the optional argument
!  "headeronly" is present and set to true
!
    if (present(headeronly)) then
       if (headeronly) return
    endif
!
!--determine if extra dust quantites should be read
!
    if (present(dustydisc)) then
       extradust = .true.
    else
       extradust = .false.
    endif
!
!--determine whether or not to read this particular block
!  onto this particular thread, either in whole or in part
!  Also handles MPI -> non-MPI dump conversion and vice-versa.
!  Can be used by non-MPI codes to read isolated blocks only.
!
    call get_blocklimits(nhydrothisblock,nblocks,nprocs,id,iblock,noffset,npartread)
    i1 = i2 + 1
    i2 = i1 + (npartread - 1)
    npart = npart + npartread
#ifdef MPI
    if (npart > maxp) then
       write(*,*) 'npart > maxp in readwrite_dumps'
       ierr = 1
       return
    endif
#endif
    if (npartread <= 0 .and. nptmass <= 0) then
       call skipblock(idisk1,nums(:,1),nums(:,2),nums(:,3),nums(:,4),tagged,ierr)
       if (ierr /= 0) then
          print*,' error skipping block'
          return
       endif
       cycle overblocks
    elseif (npartread > 0) then
#ifdef MPI
       write(*,"(a,i5,2(a,i10),a,i5,a,i10,'-',i10)") &
     'thread ',id,' reading particles ',noffset+1,':',noffset+npartread,', from block ',iblock,' lims=',i1,i2
#else
       write(*,"(2(a,i10),a,i5,a,i10,'-',i10)") &
     ' reading particles ',noffset+1,':',noffset+npartread,', from block ',iblock,' lims=',i1,i2
#endif
    else
       write(*,"(a,i10,a)") ' WARNING! block contains no SPH particles, reading ',nptmass,' point mass particles only'
    endif

    call read_phantom_arrays(i1,i2,noffset,narraylengths,nums,npartread,npartoftype,&
                          massoftype,nptmass,nsinkproperties,phantomdump,tagged,smalldump,&
                          extradust,tfile,alphafile,idisk1,iprint,ierr)

    if (ierr /= 0) call warning('read_dump','error reading arrays from file')

 enddo overblocks

 !
 ! determine npartoftype
 !
 npartoftypetot = npartoftype
 call count_particle_types(npartoftype)
 npartoftypetotact = reduceall_mpi('+',npartoftype)
 do i = 1,maxtypes
    if (npartoftypetotact(i) /= npartoftypetot(i)) then
       write(*,*) 'npartoftypetot    =',npartoftypetot
       write(*,*) 'npartoftypetotact =',npartoftypetotact
       call error('read_dump','particle type counts do not match header')
    endif
 enddo

 if (sum(npartoftype)==0) npartoftype(1) = npart
 if (narraylengths >= 4) then
    if (id==master) write(iprint,"(a,/)") ' <<< finished reading (MHD) file '
 else
    if (id==master) write(iprint,"(a,/)") ' <<< finished reading (hydro) file '
 endif
 close(idisk1)
 return

100 close (idisk1)
 if (sum(npartoftype)==0) npartoftype(1) = npart
 write(iprint,"(a,/)") ' <<< ERROR! end of file reached in data read'
 ierr = 666
 return

end subroutine read_smalldump

!--------------------------------------------------------------------
!+
!  read arrays from the main block in the file into the relevant
!  phantom modules
!+
!-------------------------------------------------------------------
subroutine read_phantom_arrays(i1,i2,noffset,narraylengths,nums,npartread,npartoftype,&
                               massoftype,nptmass,nsinkproperties,phantomdump,tagged,singleprec,&
                               extradust,tfile,alphafile,idisk1,iprint,ierr)
 use dump_utils, only:read_array,match_tag
 use dim,        only:use_dust,h2chemistry,maxalpha,maxp,gravity,maxgrav,maxvxyzu,maxBevol, &
                      store_temperature,use_dustgrowth,maxdusttypes
 use part,       only:xyzh,xyzh_label,vxyzu,vxyzu_label,dustfrac,abundance,abundance_label, &
                      alphaind,poten,xyzmh_ptmass,xyzmh_ptmass_label,vxyz_ptmass,vxyz_ptmass_label, &
                      Bevol,Bxyz,Bxyz_label,nabundances,iphase,idust,tstop,deltav,dustfrac_label, &
                      tstop_label,deltav_label,temperature,dustprop,dustprop_label,St
#ifdef IND_TIMESTEPS
 use part,       only:dt_in
#endif
 integer, intent(in)   :: i1,i2,noffset,narraylengths,nums(:,:),npartread,npartoftype(:),idisk1,iprint
 real,    intent(in)   :: massoftype(:)
 integer, intent(in)   :: nptmass,nsinkproperties
 logical, intent(in)   :: phantomdump,singleprec,tagged
 logical, intent(in)   :: extradust
 real,    intent(in)   :: tfile,alphafile
 integer, intent(out)  :: ierr
 logical               :: got_dustfrac(maxdusttypes),got_tstop(maxdusttypes),got_deltav(3,maxdusttypes)
 logical               :: match
 logical               :: got_iphase,got_xyzh(4),got_vxyzu(4),got_abund(nabundances),got_alpha,got_poten
 logical               :: got_sink_data(nsinkproperties),got_sink_vels(3),got_Bxyz(3)
 logical               :: got_psi,got_temp,got_dustprop(3),got_St
 character(len=lentag) :: tag,tagarr(64)
 integer :: k,i,iarr,ik,ndustfraci,ntstopi,ndustveli

!
!--read array type 1 arrays
!
 got_iphase    = .false.
 got_xyzh      = .false.
 got_vxyzu     = .false.
 got_dustfrac  = .false.
 got_tstop     = .false.
 got_deltav    = .false.
 got_abund     = .false.
 got_alpha     = .false.
 got_poten     = .false.
 got_sink_data = .false.
 got_sink_vels = .false.
 got_Bxyz      = .false.
 got_psi       = .false.
 got_temp      = .false.
 got_dustprop  = .false.
 got_St        = .false.

 ndustfraci = 0
 ntstopi    = 0
 ndustveli  = 0

 over_arraylengths: do iarr=1,narraylengths

    do k=1,ndatatypes
       ik = k
       if (singleprec .and. k==i_real) ik = i_real4
       if (.not.tagged) call fake_array_tags(iarr,k,tagarr,phantomdump)
       do i=1,nums(k,iarr)
          match = .false.
          if (tagged) then
             read(idisk1,end=100) tag
          else
             tag = tagarr(i)
          endif
          !write(*,*) 'CHECKING '//trim(tag)
          select case(iarr)
          case(1)
             call read_array(iphase,'itype',got_iphase,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             call read_array(xyzh, xyzh_label, got_xyzh, ik,i1,i2,noffset,idisk1,tag,match,ierr)
             call read_array(vxyzu,vxyzu_label,got_vxyzu,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             if (use_dust) then
                if (any(tag == dustfrac_label)) then
                   ndustfraci = ndustfraci + 1
                   call read_array(dustfrac(ndustfraci,:),dustfrac_label(ndustfraci),got_dustfrac(ndustfraci), &
                                   ik,i1,i2,noffset,idisk1,tag,match,ierr)
                endif
             endif
             if (extradust) then
                if (any(tag == tstop_label)) then
                   ntstopi = ntstopi + 1
                   call read_array(tstop(ntstopi,:),tstop_label(ntstopi),got_tstop(ntstopi),&
                                   ik,i1,i2,noffset,idisk1,tag,match,ierr)
                endif
                if (any(tag == deltav_label)) then
                   !--use deltavx to identify each new dust species
                   if (tag == deltav_label(1)) ndustveli = ndustveli + 1
                   call read_array(deltav(:,ndustveli,:),deltav_label,got_deltav(:,ndustveli), &
                                   ik,i1,i2,noffset,idisk1,tag,match,ierr)
                endif
             endif
             if (use_dustgrowth) then
                call read_array(dustprop,dustprop_label,got_dustprop,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                call read_array(St,'St',got_St,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (h2chemistry) then
                call read_array(abundance,abundance_label,got_abund,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (store_temperature) then
                call read_array(temperature,'T',got_temp,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (maxalpha==maxp) call read_array(alphaind(1,:),'alpha',got_alpha,ik,i1,i2,noffset,idisk1,tag,match,ierr)

             !
             ! read gravitational potential if it is in the file
             !
             if (gravity .and. maxgrav==maxp) call read_array(poten,'poten',got_poten,ik,i1,i2,noffset,idisk1,tag,match,ierr)
#ifdef IND_TIMESTEPS
             !
             ! read dt if it is in the file
             !
             call read_array(dt_in,'dt',dt_read_in,ik,i1,i2,noffset,idisk1,tag,match,ierr)
#endif
          case(2)
             call read_array(xyzmh_ptmass,xyzmh_ptmass_label,got_sink_data,ik,1,nptmass,0,idisk1,tag,match,ierr)
             call read_array(vxyz_ptmass, vxyz_ptmass_label, got_sink_vels,ik,1,nptmass,0,idisk1,tag,match,ierr)
          case(4)
             call read_array(Bxyz,Bxyz_label,got_Bxyz,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             call read_array(Bevol(4,:),'psi',got_psi,ik,i1,i2,noffset,idisk1,tag,match,ierr)
          end select
          if (.not.match) then
             !write(*,*) 'skipping '//trim(tag)
             read(idisk1,end=100) ! skip unknown array
          endif
       enddo
    enddo

 enddo over_arraylengths
 !
 ! check for errors
 !
 call check_arrays(i1,i2,npartoftype,npartread,nptmass,nsinkproperties,massoftype,&
                   alphafile,tfile,phantomdump,got_iphase,got_xyzh,got_vxyzu,got_alpha, &
                   got_abund,got_dustfrac,got_sink_data,got_sink_vels,got_Bxyz,got_psi,got_dustprop,got_St, &
                   got_temp,iphase,xyzh,vxyzu,alphaind,xyzmh_ptmass,Bevol,iprint,ierr)

 return
100 continue
 write(iprint,"(a,/)") ' <<< ERROR! end of file reached in data read'

end subroutine read_phantom_arrays

!--------------------------------------------------------------------
!+
!  small utility to see if a parameter is different between the
!  code and the dump file
!+
!-------------------------------------------------------------------
subroutine checkparam(valfile,valcode,string)
 use io, only:iprint,id,master
 real,             intent(in) :: valfile,valcode
 character(len=*), intent(in) :: string

 if (id==master) then
    if (abs(valfile-valcode) > tiny(valcode)) then
       write(iprint,*) 'comment: '//trim(string)//' was ',valfile,' now ',valcode
    endif
 endif

 return
end subroutine checkparam

!------------------------------------------------------------
!+
!  subroutine to perform sanity checks on the array headers
!+
!------------------------------------------------------------
subroutine check_block_header(narraylengths,nblocks,ilen,nums,nparttot,nhydrothisblock,nptmass,ierr)
 use dim, only:maxptmass
 use io,  only:warning
 integer,         intent(in)  :: narraylengths,nblocks
 integer(kind=8), intent(in)  :: ilen(narraylengths),nparttot
 integer,         intent(in)  :: nums(:,:)
 integer(kind=8), intent(out) :: nhydrothisblock
 integer,         intent(out) :: nptmass,ierr

 nhydrothisblock = ilen(1)
 if (nblocks==1 .and. nhydrothisblock < nparttot) then
    ierr = 8
    write (*,*) 'ERROR in read_dump: npart wrong',nhydrothisblock,nparttot
    return
 endif
 if (narraylengths >= 2) then
    if (ilen(2) > 0) then
       nptmass = int(ilen(2),kind=4)
       if (nptmass > maxptmass) then
          write (*,*) 'error in readdump: nptmass = ',nptmass,' > maxptmass (',maxptmass,&
                      '): recompile using make MAXPTMASS=',nptmass
          ierr = 9
          return
       endif
    else
       if (ilen(2) < 0) then
          write(*,*) 'error in readdump: nptmass < 0 in dump header'
          ierr = 10
          return
       endif
       nptmass = 0
    endif
 else
    nptmass = 0
 endif

 if (nptmass == 0 .and. nums(i_real,2) > 0) then
    call warning('read_dump','got nptmass = 0 from header, but file has sink info (skipping...)')
 endif

 if (narraylengths >= 4) then
    if (ilen(4) /= nhydrothisblock) then
       write(*,*) 'ERROR: MHD dimensions differ from hydro ',ilen(4),nhydrothisblock
       ierr = 9
       return
    endif
 endif

end subroutine check_block_header

!---------------------------------------------------------------
!+
!  make sure required arrays have been read from Phantom file
!  and perform basic sanity checks
!+
!---------------------------------------------------------------
subroutine check_arrays(i1,i2,npartoftype,npartread,nptmass,nsinkproperties,massoftype,&
                        alphafile,tfile,phantomdump,got_iphase,got_xyzh,got_vxyzu,got_alpha, &
                        got_abund,got_dustfrac,got_sink_data,got_sink_vels,got_Bxyz,got_psi,got_dustprop,got_St, &
                        got_temp,iphase,xyzh,vxyzu,alphaind,xyzmh_ptmass,Bevol,iprint,ierr)
 use dim,  only:maxp,maxvxyzu,maxalpha,maxBevol,mhd,h2chemistry,store_temperature,use_dustgrowth
 use eos,  only:polyk,gamma
 use part, only:maxphase,isetphase,set_particle_type,igas,ihacc,ihsoft,imacc,&
                xyzmh_ptmass_label,vxyz_ptmass_label,get_pmass,rhoh,dustfrac
 use io,   only:warning,id,master
 use options,    only:alpha,use_dustfrac
 use sphNGutils, only:itype_from_sphNG_iphase,isphNG_accreted
 integer,         intent(in)    :: i1,i2,npartoftype(:),npartread,nptmass,nsinkproperties
 real,            intent(in)    :: massoftype(:),alphafile,tfile
 logical,         intent(in)    :: phantomdump,got_iphase,got_xyzh(:),got_vxyzu(:),got_alpha,got_dustprop(:),got_St
 logical,         intent(in)    :: got_abund(:),got_dustfrac(:),got_sink_data(:),got_sink_vels(:),got_Bxyz(:)
 logical,         intent(in)    :: got_psi, got_temp
 integer(kind=1), intent(inout) :: iphase(:)
 real,            intent(inout) :: vxyzu(:,:), Bevol(:,:)
 real(kind=4),    intent(inout) :: alphaind(:,:)
 real,            intent(inout) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer,         intent(in)    :: iprint
 integer,         intent(out)   :: ierr
 logical :: use_gas
 integer :: i,itype,nread
 !
 ! particle type information
 !
 if (maxphase==maxp) then
    if (got_iphase) then
       if (phantomdump) then
          do i=i1,i2
             itype = int(iphase(i))
             iphase(i) = isetphase(itype,iactive=.true.)
          enddo
       else
          ! convert from sphNG
          do i=i1,i2
             itype = itype_from_sphNG_iphase(iphase(i))
             iphase(i) = isetphase(itype,iactive=.true.)
             if (itype==isphNG_accreted) then
                !  mark accreted/unknown particle types as dead according
                !  to Phantom (i.e., give negative smoothing length)
                xyzh(4,i) = -abs(xyzh(4,i))
                call set_particle_type(i,igas) ! to give an allowed particle type
             endif
          enddo
       endif
    elseif (any(npartoftype(2:) > 0)) then
       !
       !--iphase is required if there is more than one particle type
       !
       write(*,*) 'error in rdump: need type information but iamtype not present in dump file'
       ierr = 8
       return
    else
       !
       !--iphase does not need to be read if there is only one particle type
       !  but to start the code it should be set such that all particles are active
       !
       do i=i1,i2
          iphase(i) = isetphase(igas,iactive=.true.)
       enddo
    endif
 endif
 if (maxphase==maxp) then
    use_gas = .false.
 else
    use_gas = .true.
 endif

 !
 ! hydrodynamics arrays
 !
 if (any(.not.got_xyzh)) then
    if (id==master .and. i1==1) write(*,*) 'ERROR: x, y, z or h not found in file'
    ierr = 9
    return
 endif
 if (any(.not.got_vxyzu(1:3))) then
    if (id==master .and. i1==1) write(*,*) 'ERROR: missing velocity information from file'
 endif
 if (maxvxyzu==4 .and. .not.got_vxyzu(4)) then
    if (gamma < 1.01) then
       do i=i1,i2
          vxyzu(4,i) = 1.5*polyk
          !print*,'u = ',vxyzu(4,i)
       enddo
       if (id==master .and. i1==1) write(*,*) 'WARNING: u not in file but setting u = 3/2 * cs^2'
    else
       do i=i1,i2
          vxyzu(4,i) = (1.0/(gamma-1.0))*polyk*rhoh(xyzh(4,i),get_pmass(i,use_gas))**(gamma - 1.)
          !print*,'u = ',vxyzu(4,i)
       enddo
       if (id==master .and. i1==1) write(*,*) 'WARNING: u not in file but setting u = (K*rho**(gamma-1))/(gamma-1)'
    endif
 endif
 if (h2chemistry .and. .not.all(got_abund)) then
    if (id==master) write(*,*) 'error in rdump: using H2 chemistry, but abundances not found in dump file'
    ierr = 9
    return
 endif
 if (store_temperature .and. .not.got_temp) then
    if (id==master .and. i1==1) write(*,*) 'WARNING: missing temperature information from file'
 endif
 if (maxalpha==maxp) then
    if (got_alpha) then
       if (alphafile < 0.99 .and. tfile > 0.) then
          if (any(alphaind(1,i1:i2) > 1.0 .or. alphaind(1,i1:i2) < 0.)) then
             if (id==master) write(iprint,*) 'ERROR! AV alpha < 0 or alpha > 1 in dump file: using alpha'
             alphaind(1,i1:i2) = real(alpha,kind=4)
          endif
       endif
    else
       if (id==master .and. i1==1) write(*,*) 'WARNING: alpha not found in file'
       alphaind(1,i1:i2) = real(alpha,kind=4)
    endif
 endif
 if (any(massoftype <= 0. .and. npartoftype /= 0) .and. npartread > 0) then
    if (id==master .and. i1==1) write(*,*) 'ERROR! mass not set in read_dump (Phantom)'
    ierr = 12
    return
 endif
 if (use_dustfrac .and. .not. all(got_dustfrac)) then
    if (id==master .and. i1==1) write(*,*) 'ERROR! using one-fluid dust, but no dust fraction found in dump file'
    if (id==master .and. i1==1) write(*,*) ' Setting dustfrac = 0'
    dustfrac = 0.
    !ierr = 13
    return
 endif
 if (use_dustgrowth .and. .not.got_dustprop(1)) then
    write(*,*) 'ERROR! using dustgrowth, but no grain size found in dump file'
    return
 endif
 if (use_dustgrowth .and. .not.got_dustprop(2)) then
    write(*,*) 'ERROR! using dustgrowth, but no grain density found in dump file'
    return
 endif
 if (use_dustgrowth .and. .not.got_dustprop(3)) then
    write(*,*) 'ERROR! using dustgrowth, but no ratio vrel/vfrag found in dump file'
    return
 endif
 if (use_dustgrowth .and. .not.got_St) then
    write(*,*) 'ERROR! using dustgrowth, but no Stokes number found in dump file'
    return
 endif
 !
 ! sink particle arrays
 !
 nread = 0
 if (nptmass > 0) then
    do i=1,nsinkproperties
       if (.not.got_sink_data(i)) then
          if (i <= 5) then
             if (id==master) write(*,*) 'ERROR! sink particle '//trim(xyzmh_ptmass_label(i))//' not found'
             ierr = 10
             return
          else
             if (id==master) write(*,*) 'WARNING! sink particle '//trim(xyzmh_ptmass_label(i))//' not found'
          endif
       endif
    enddo
    if (.not.all(got_sink_vels(1:3))) then
       if (id==master .and. i1==1) write(*,*) 'WARNING! sink particle velocities not found'
    endif
    if (id==master .and. i1==1) then
       print "(2(a,i2),a)",' got ',nsinkproperties,' sink properties from ',nptmass,' sink particles'
       if (nptmass > 0) print "(1x,47('-'),/,1x,a,'|',4(a9,1x,'|'),/,1x,47('-'))",&
                              'ID',' Mass    ',' Racc    ',' Macc    ',' hsoft   '
       do i=1,min(nptmass,999)
          print "(i3,'|',4(1pg9.2,1x,'|'))",i,xyzmh_ptmass(4,i),xyzmh_ptmass(ihacc,i),xyzmh_ptmass(imacc,i),xyzmh_ptmass(ihsoft,i)
       enddo
       if (nptmass > 0) print "(1x,47('-'))"
    endif
 endif

 !
 ! MHD arrays
 !
 if (mhd) then
    if (.not.all(got_Bxyz(1:3))) then
       if (id==master .and. i1==1) write(*,*) 'WARNING: MHD but magnetic field arrays not found in Phantom dump file'
    endif
    if (maxBevol==4 .and. .not.got_psi) then
       if (id==master .and. i1==1) write(*,*) 'WARNING! div B cleaning field (Psi) not found in Phantom dump file: assuming psi=0'
       Bevol(maxBevol,i1:i2) = 0.
    endif
 endif

end subroutine check_arrays

!--------------------------------------------------------------------
!+
!  utility to extract header variables to phantom
!+
!-------------------------------------------------------------------
subroutine unfill_header(hdr,phantomdump,got_tags,nparttot, &
                         nblocks,npart,npartoftype, &
                         tfile,hfactfile,alphafile,iprint,id,nprocs,ierr)
 use dim,        only:maxp,maxdustlarge,use_dust
 use io,         only:master ! check this
 use eos,        only:isink
 use part,       only:maxtypes,igas,idust,ndustsmall,ndustlarge,ndusttypes
 use units,      only:udist,umass,utime,set_units_extra,set_units
 use dump_utils, only:extract,dump_h
 use fileutils,  only:make_tags_unique
 type(dump_h),    intent(in)  :: hdr
 logical,         intent(in)  :: phantomdump,got_tags
 integer(kind=8), intent(out) :: nparttot
 integer,         intent(out) :: nblocks,npart,npartoftype(maxtypes)
 real,            intent(out) :: tfile,hfactfile,alphafile
 integer,         intent(in)  :: iprint,id,nprocs
 integer,         intent(out) :: ierr
 integer         :: nparttoti,npartoftypetoti(maxtypes),ntypesinfile
 integer         :: ierr1,ierrs(3),i,counter
 integer(kind=8) :: npartoftypetot(maxtypes),ntypesinfile8
 character(len=10) :: dust_label(maxdustlarge)

 ierr = 0
 nparttot = 0
 npartoftypetot(:) = 0
 npart = 0
 npartoftype(:) = 0
 isink = 0
 call extract('ntypes',ntypesinfile,hdr,ierr1)
 if (ierr1 /= 0 .or. ntypesinfile < 1) then
    if (phantomdump .and. got_tags) then
       ierr = 4
       return
    else
       ntypesinfile = 5
    endif
 endif

 ! extract quantities from integer header
 call extract('nparttot',nparttoti,hdr,ierr1)
 if (ierr1 /= 0) then
    ierr = 5
    return
 endif
 if (ntypesinfile > maxtypes) then
    write(*,*) 'WARNING: number of particle types in file exceeds array limits'
    ntypesinfile = maxtypes
 endif
 call extract('npartoftype',npartoftypetoti(1:ntypesinfile),hdr,ierr1)
 if (ierr1 /= 0) then
    npartoftype(1) = nparttoti  ! assume only gas particles
 endif
 call extract('nblocks',nblocks,hdr,ierr1,default=1)
 if (ierr1 /= 0) write(*,*) 'number of MPI blocks not read: assuming 1'

 nparttot = int(nparttoti,kind=8)
 npartoftypetot = int(npartoftypetoti,kind=8)
 if (nblocks==1) then
    npartoftype(1:ntypesinfile) = int(npartoftypetot(1:ntypesinfile))
    if (npartoftype(idust) > 0) write(*,*) 'n(gas) = ',npartoftype(igas)
    counter = 0
    do i=1,maxdustlarge
       if (npartoftype(idust+i-1) > 0) then
          counter = counter + 1
       endif
    enddo
    dust_label = 'dust'
    call make_tags_unique(counter,dust_label)
    do i=1,counter
       write(*,*) 'n('//trim(dust_label(i))//') = ',npartoftype(idust+i-1)
    enddo
 endif
 call extract('isink',isink,hdr,ierr1)

!--non-MPI dumps
 if (nprocs==1) then
    if (nparttoti > maxp) then
       write (*,*) 'ERROR in readdump: number of particles exceeds MAXP: recompile with MAXP=',nparttoti
       ierr = 4
       return
    elseif (nparttoti > huge(npart)) then
       write (*,*) 'ERROR in readdump: number of particles exceeds 32 bit limit, must use int(kind=8)''s ',nparttoti
       ierr = 4
       return
    endif
 endif
 if (nblocks==1) then
    npart = int(nparttoti)
    nparttot = npart
    if (id==master) write (iprint,*) 'npart = ',npart
 endif
 if (got_tags) then
    call extract('ntypes',ntypesinfile8,hdr,ierr1)
    ntypesinfile = int(ntypesinfile8)
 elseif (ntypesinfile < 1) then
    ntypesinfile = 5
 endif
 call extract('nparttot',nparttot,hdr,ierr1)
 if (nblocks > 1) then
    call extract('npartoftype',npartoftype(1:ntypesinfile),hdr,ierr1)
 endif
 if (id==master) write(*,*) 'npart(total) = ',nparttot

 call extract('udist',udist,hdr,ierrs(1))
 call extract('umass',umass,hdr,ierrs(2))
 call extract('utime',utime,hdr,ierrs(3))
 if (all(ierrs(1:3)==0)) then
    call set_units_extra()
 else
    write(iprint,*) 'ERROR reading units from dump file, assuming default'
    call set_units()  ! use default units
 endif

!--default real
 call unfill_rheader(hdr,phantomdump,ntypesinfile,&
                     tfile,hfactfile,alphafile,iprint,ierr)

 if (use_dust) then
    call extract('ndustsmall',ndustsmall,hdr,ierrs(1))
    call extract('ndustlarge',ndustlarge,hdr,ierrs(2))
    ndusttypes = ndustsmall + ndustlarge
    if (any(ierrs(1:2) /= 0)) then
       write(*,*) 'ERROR reading number of small/large grain types from file header'
    endif
 endif

 if (ierr /= 0) return

 if (id==master) write(iprint,*) 'time = ',tfile

end subroutine unfill_header

!--------------------------------------------------------------------
!+
!  subroutine to fill the real header with various things
!+
!-------------------------------------------------------------------
subroutine fill_header(sphNGdump,t,nparttot,npartoftypetot,nblocks,nptmass,hdr,ierr)
 use eos,            only:polyk,gamma,polyk2,qfacdisc,isink
 use options,        only:tolh,alpha,alphau,alphaB,iexternalforce,ieos
 use part,           only:massoftype,hfact,Bextx,Bexty,Bextz,ndustsmall,ndustlarge,&
                          idust,grainsize,graindens
 use initial_params, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use setup_params,   only:rhozero
 use timestep,       only:dtmax,C_cour,C_force
 use externalforces, only:write_headeropts_extern
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax
 use dump_utils,     only:reset_header,add_to_rheader,add_to_header,add_to_iheader,num_in_header
 use dim,            only:use_dust,maxtypes,use_dustgrowth, &
                          phantom_version_major,phantom_version_minor,phantom_version_micro
 use units,          only:udist,umass,utime,unit_Bfield
 logical,         intent(in)    :: sphNGdump
 real,            intent(in)    :: t
 integer(kind=8), intent(in)    :: nparttot,npartoftypetot(:)
 integer,         intent(in)    :: nblocks,nptmass
 type(dump_h),    intent(inout) :: hdr
 integer,         intent(out)   :: ierr
 integer :: number

 ierr = 0
 ! default int
 call add_to_iheader(int(nparttot),'nparttot',hdr,ierr)
 call add_to_iheader(maxtypes,'ntypes',hdr,ierr)
 call add_to_iheader(int(npartoftypetot(1:maxtypes)),'npartoftype',hdr,ierr)
 call add_to_iheader(nblocks,'nblocks',hdr,ierr)
 call add_to_iheader(isink,'isink',hdr,ierr)
 call add_to_iheader(nptmass,'nptmass',hdr,ierr)
 call add_to_iheader(ndustlarge,'ndustlarge',hdr,ierr)
 call add_to_iheader(ndustsmall,'ndustsmall',hdr,ierr)
 call add_to_iheader(idust,'idust',hdr,ierr)
 call add_to_iheader(phantom_version_major,'majorv',hdr,ierr)
 call add_to_iheader(phantom_version_minor,'minorv',hdr,ierr)
 call add_to_iheader(phantom_version_micro,'microv',hdr,ierr)

 ! int*8
 call add_to_header(nparttot,'nparttot',hdr,ierr)
 call add_to_header(int(maxtypes,kind=8),'ntypes',hdr,ierr)
 call add_to_header(npartoftypetot(1:maxtypes),'npartoftype',hdr,ierr)

 ! int*4
 call add_to_header(iexternalforce,'iexternalforce',hdr,ierr)
 call add_to_header(ieos,'ieos',hdr,ierr)

 ! default real variables
 call add_to_rheader(t,'time',hdr,ierr)
 call add_to_rheader(dtmax,'dtmax',hdr,ierr)
 call add_to_rheader(gamma,'gamma',hdr,ierr)
 call add_to_rheader(rhozero,'rhozero',hdr,ierr)
 call add_to_rheader(1.5*polyk,'RK2',hdr,ierr)
 if (sphNGdump) then ! number = 23
    call add_to_rheader(0.,'escaptot',hdr,ierr)
    call add_to_rheader(0.,'tkin',hdr,ierr)
    call add_to_rheader(0.,'tgrav',hdr,ierr)
    call add_to_rheader(0.,'tterm',hdr,ierr)
    call add_to_rheader(0.,'anglostx',hdr,ierr)
    call add_to_rheader(0.,'anglosty',hdr,ierr)
    call add_to_rheader(0.,'anglostz',hdr,ierr)
    call add_to_rheader(0.,'specang',hdr,ierr)
    call add_to_rheader(0.,'ptmassin',hdr,ierr)
    call add_to_rheader(0.,'tmag',hdr,ierr)
    call add_to_rheader(Bextx,'Bextx',hdr,ierr)
    call add_to_rheader(Bexty,'Bexty',hdr,ierr)
    call add_to_rheader(Bextz,'Bextz',hdr,ierr)
    call add_to_rheader(0.,'hzero',hdr,ierr)
    call add_to_rheader(1.5*polyk2,'uzero_n2',hdr,ierr)
    call add_to_rheader(0.,'hmass',hdr,ierr)
    call add_to_rheader(0.,'gapfac',hdr,ierr)
    call add_to_rheader(0.,'pmassinitial',hdr,ierr)
 else ! number = 49
    call add_to_rheader(hfact,'hfact',hdr,ierr)
    call add_to_rheader(tolh,'tolh',hdr,ierr)
    call add_to_rheader(C_cour,'C_cour',hdr,ierr)
    call add_to_rheader(C_force,'C_force',hdr,ierr)
    call add_to_rheader(alpha,'alpha',hdr,ierr)
    call add_to_rheader(alphau,'alphau',hdr,ierr)
    call add_to_rheader(alphaB,'alphaB',hdr,ierr)
    call add_to_rheader(polyk2,'polyk2',hdr,ierr)
    call add_to_rheader(qfacdisc,'qfacdisc',hdr,ierr)
    call add_to_rheader(massoftype,'massoftype',hdr,ierr) ! array
    call add_to_rheader(Bextx,'Bextx',hdr,ierr)
    call add_to_rheader(Bexty,'Bexty',hdr,ierr)
    call add_to_rheader(Bextz,'Bextz',hdr,ierr)
    call add_to_rheader(0.,'dum',hdr,ierr)
    if (iexternalforce /= 0) call write_headeropts_extern(iexternalforce,hdr,t,ierr)
    call add_to_rheader(xmin,'xmin',hdr,ierr)
    call add_to_rheader(xmax,'xmax',hdr,ierr)
    call add_to_rheader(ymin,'ymin',hdr,ierr)
    call add_to_rheader(ymax,'ymax',hdr,ierr)
    call add_to_rheader(zmin,'zmin',hdr,ierr)
    call add_to_rheader(zmax,'zmax',hdr,ierr)
    call add_to_rheader(get_conserv,'get_conserv',hdr,ierr)
    call add_to_rheader(etot_in,'etot_in',hdr,ierr)
    call add_to_rheader(angtot_in,'angtot_in',hdr,ierr)
    call add_to_rheader(totmom_in,'totmom_in',hdr,ierr)
    call add_to_rheader(mdust_in,'mdust_in',hdr,ierr)
    if (use_dust) then
       call add_to_rheader(grainsize,'grainsize',hdr,ierr)
       call add_to_rheader(graindens,'graindens',hdr,ierr)
    endif
 endif

 ! real*8
 call add_to_header(udist,'udist',hdr,ierr)
 call add_to_header(umass,'umass',hdr,ierr)
 call add_to_header(utime,'utime',hdr,ierr)
 call add_to_header(unit_Bfield,'umagfd',hdr,ierr)

 if (ierr /= 0) write(*,*) ' ERROR: arrays too small writing rheader'

 number = num_in_header(hdr%realtags)
 if (number >= maxphead) then
    write(*,*) 'error: header arrays too small for number of items in header: will be truncated'
 endif

end subroutine fill_header

!--------------------------------------------------------------------
!+
!  subroutine to set runtime parameters having read the real header
!+
!-------------------------------------------------------------------
subroutine unfill_rheader(hdr,phantomdump,ntypesinfile,&
                          tfile,hfactfile,alphafile,iprint,ierr)
 use io,             only:id,master
 use dim,            only:maxp,maxvxyzu,use_dust
 use eos,            only:polyk,gamma,polyk2,qfacdisc,extract_eos_from_hdr
 use options,        only:ieos,tolh,alpha,alphau,alphaB,iexternalforce
 use part,           only:massoftype,hfact,Bextx,Bexty,Bextz,mhd,periodic,&
                          maxtypes,grainsize,graindens
 use initial_params, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use setup_params,   only:rhozero
 use timestep,       only:dtmax,C_cour,C_force
 use externalforces, only:read_headeropts_extern
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax,set_boundary
 use dump_utils,     only:extract
 type(dump_h), intent(in)  :: hdr
 logical,      intent(in)  :: phantomdump
 integer,      intent(in)  :: iprint,ntypesinfile
 real,         intent(out) :: tfile,hfactfile,alphafile
 integer,      intent(out) :: ierr

 integer, parameter :: lu = 173
 integer            :: ierrs(10),iextern_in_file
 real               :: rk2,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,dtmaxi
 real               :: alphaufile,alphaBfile,C_courfile,C_forcefile,tolhfile
 logical            :: iexist

 ierr  = 0
 call extract('time',tfile,hdr,ierr)
 if (ierr/=0)  call extract('gt',tfile,hdr,ierr)  ! this is sphNG's label for time
 call extract('dtmax',dtmaxi,hdr,ierr)
 call checkparam(dtmaxi,dtmax,'dtmax')

 call extract('gamma',gamma,hdr,ierr)
 call extract('rhozero',rhozero,hdr,ierr)
 call extract('RK2',rk2,hdr,ierr)
 polyk = 2./3.*rk2
 if (maxvxyzu >= 4) then
    if (id==master) write(iprint,*) 'adiabatic eos: gamma = ',gamma
 else
    write(iprint,*) 'setting isothermal sound speed^2 (polyk) = ',polyk,' gamma = ',gamma
    if (polyk <= tiny(polyk)) then
       write(iprint,*) 'WARNING! sound speed zero in dump!, polyk = ',polyk
    endif
 endif

 Bextx = 0.
 Bexty = 0.
 Bextz = 0.
 if (phantomdump) then
    call extract('hfact',hfactfile,hdr,ierr)
    call extract('tolh',tolhfile,hdr,ierr)
    call extract('C_cour',C_courfile,hdr,ierr)
    call extract('C_force',C_forcefile,hdr,ierr)
    call checkparam(hfactfile,hfact,'hfact')
    call checkparam(tolhfile,tolh,'tolh')
    call checkparam(C_courfile,C_cour,'C_cour')
    call checkparam(C_forcefile,C_force,'C_force')

    call extract('alpha',alphafile,hdr,ierr)
    call checkparam(alphafile,alpha,'alpha')
    if (maxvxyzu >= 4) then
       call extract('alphau',alphaufile,hdr,ierr)
       call checkparam(alphaufile,alphau,'alphau')
    else
       alphaufile = 0.
    endif
    if (mhd) then
       call extract('alphaB',alphaBfile,hdr,ierr)
       call checkparam(alphaBfile,alphaB,'alphaB')
    endif
    call extract('polyk2',polyk2,hdr,ierr)
    call extract('qfacdisc',qfacdisc,hdr,ierr)
    if (extract_eos_from_hdr) call extract('ieos',ieos,hdr,ierr)
    if (ieos==3) then
       if (qfacdisc <= tiny(qfacdisc)) then
          write(iprint,*) 'ERROR: qfacdisc <= 0'
          ierr = 2
       else
          write(iprint,*) 'qfacdisc = ',qfacdisc
       endif
    endif
    call extract('massoftype',massoftype(1:ntypesinfile),hdr,ierr)
    if (ierr /= 0) then
       write(*,*) '*** ERROR reading massoftype from dump header ***'
       ierr = 4
    endif
    call extract('iexternalforce',iextern_in_file,hdr,ierrs(1))
    if (iexternalforce /= 0) then
       call read_headeropts_extern(iexternalforce,hdr,ierrs(1))
       if (ierrs(1) /= 0) ierr = 5
    elseif (iextern_in_file /= 0) then
       call read_headeropts_extern(iextern_in_file,hdr,ierrs(1))
       if (ierrs(1) /= 0) ierr = 5
    endif
 else
    massoftype(1) = 0.
    polyk2 = 0.
    hfactfile = 0.
 endif

 if (periodic) then
    call extract('xmin',xmini,hdr,ierrs(1))
    call extract('xmax',xmaxi,hdr,ierrs(2))
    call extract('ymin',ymini,hdr,ierrs(3))
    call extract('ymax',ymaxi,hdr,ierrs(4))
    call extract('zmin',zmini,hdr,ierrs(5))
    call extract('zmax',zmaxi,hdr,ierrs(6))
    if (any(ierrs(1:6) /= 0)) then
       write(*,"(2(/,a))") ' ERROR: dump does not contain boundary positions', &
                           '        but we are using periodic boundaries'
       inquire(file='bound.tmp',exist=iexist)
       if (iexist) then
          open(unit=lu,file='bound.tmp')
          read(lu,*) xmini,xmaxi,ymini,ymaxi,zmini,zmaxi
          close(lu)
          call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
          write(*,"(a,6(es10.3,1x))") ' READ from bound.tmp ',xmin,xmax,ymin,ymax,zmin,zmax
       else
          write(*,"(3(/,a),/,/,a)") ' To silence this error and restart from an older dump file ', &
                           ' create an ascii file called "bound.tmp" in the current directory', &
                           ' with xmin,xmax,ymin,ymax,zmin & zmax in it, e.g.: ', &
                           ' 0. 1. 0. 1. 0. 1.'
          ierr = 5  ! spit fatal error
       endif
    else
       call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
    endif
 endif

 if (mhd) then
    call extract('Bextx',Bextx,hdr,ierrs(1))
    call extract('Bexty',Bexty,hdr,ierrs(2))
    call extract('Bextz',Bextz,hdr,ierrs(3))
    if (any(ierrs(1:3) /= 0)) then
       write(*,*) 'ERROR reading external field (setting to zero)'
    else
       write(*,*) 'External field found, Bext = ',Bextx,Bexty,Bextz
    endif
 endif

 ! values to track that conserved values remain conserved
 call extract('get_conserv',get_conserv,hdr,ierrs(1))
 call extract('etot_in',    etot_in,    hdr,ierrs(2))
 call extract('angtot_in',  angtot_in,  hdr,ierrs(3))
 call extract('totmom_in',  totmom_in,  hdr,ierrs(4))
 call extract('mdust_in',   mdust_in,   hdr,ierrs(5))
 if (any(ierrs(1:5) /= 0)) then
    write(*,*) 'ERROR reading values to verify conservation laws.  Resetting initial values.'
    get_conserv = 1.0
 endif

 if (abs(gamma-1.) > tiny(gamma) .and. maxvxyzu < 4) then
    write(*,*) 'WARNING! compiled for isothermal equation of state but gamma /= 1, gamma=',gamma
 endif

 !--pull grain size and density arrays
 if (use_dust) then
    call extract('grainsize',grainsize,hdr,ierrs(1))
    call extract('graindens',graindens,hdr,ierrs(2))
    if (any(ierrs(1:2) /= 0)) then
       write(*,*) 'ERROR reading grain size/density from file header'
    endif
 endif

 return
end subroutine unfill_rheader


!-----------------------------------------------------------------
!+
!  if tags not read, give expected order of variables in header
!  this is for backwards compatibility with old (untagged) format
!+
!-----------------------------------------------------------------
subroutine fake_header_tags(hdr,phantomdump,mhd,maxtypes)
 type(dump_h), intent(inout) :: hdr
 logical,      intent(in)    :: phantomdump,mhd
 integer,      intent(in)    :: maxtypes
 character(len=lentag) :: tagarr(49)
 integer :: nread,n

 ! default int
 tagarr(1) = 'nparttot'
 tagarr(2:6) =  'npartoftype'
 tagarr(7) = 'nblocks'
 tagarr(8) = 'isink'
 if (allocated(hdr%inttags)) then
    n = min(size(hdr%inttags),8)
    hdr%inttags(1:n) = tagarr(1:n)
 endif

 ! int*8
 tagarr(1) = 'nparttot'
 tagarr(2:1+maxtypes) = 'npartoftype'
 if (allocated(hdr%int8tags)) then
    n = min(size(hdr%int8tags),1+maxtypes)
    hdr%int8tags(1:n) = tagarr(1:n)
 endif

 ! default real
 nread = hdr%nums(i_real)
 tagarr = ''
 if (nread > 5) then
    tagarr(1:5) = (/'time   ','dtmax  ','gamma  ','rhozero','RK2    '/)
 endif
 if (phantomdump) then
    if (nread >= 14) then
       tagarr(6:14) = (/'hfact   ','tolh    ','C_cour  ','C_force ', &
                       'alpha   ','alphau  ','alphaB  ','polyk2  ','qfacdisc'/)
    endif
    if (nread >= 19) tagarr(15:19) = 'massoftype'
    if (mhd .and. nread >= 22) tagarr(20:22) = (/'Bextx','Bexty','Bextz'/)

    ! 20 quantities related to external binary potential
    if (nread >= 24) then
       tagarr(24:40) = (/'x1 ','y1 ','z1 ','m1 ','h1 ','x2 ','y2 ','z2 ','m2 ', &
                        'h2 ','vx1','vy1','vz1','vx2','vy2','vz2','a0 '/)
       tagarr(41:43) = (/'direction    ','accretedmass1','accretedmass2'/)
    endif

    if (nread >= 49) tagarr(44:49) = (/'xmin','xmax','ymin','ymax','zmin','zmax'/)
 else
    if (mhd .and. nread >= 18) tagarr(16:18) = (/'Bextx','Bexty','Bextz'/)
 endif
 if (allocated(hdr%realtags)) then
    n = min(size(hdr%realtags),nread)
    hdr%realtags(1:n) = tagarr(1:n)
 endif

 ! real*8
 tagarr(1:3) = (/'udist','umass','utime'/)
 if (allocated(hdr%real8tags)) then
    n = min(size(hdr%real8tags),3)
    hdr%real8tags(1:n) = tagarr(1:n)
 endif

end subroutine fake_header_tags

!-----------------------------------------------------------------
!+
!  if tags not read, give expected order of variables in header
!  this is for backwards compatibility with old (untagged) format
!+
!-----------------------------------------------------------------
subroutine fake_array_tags(iblock,ikind,tags,phantomdump)
 use dim, only:maxalpha,maxp,maxvxyzu,h2chemistry
 integer, intent(in) :: iblock,ikind
 logical, intent(in) :: phantomdump
 character(len=lentag), intent(out) :: tags(:)
 integer :: ilen

 tags = ''
 select case(iblock)
 case(1) ! hydro arrays
    select case(ikind)
    case(i_int1)
       tags = (/'itype'/)
    case(i_real)
       if (phantomdump) then
          ilen = 6
          tags(1:ilen) = (/'x ','y ','z ','vx','vy','vz'/)
       else
          ilen = 8
          tags(1:ilen) = (/'x ','y ','z ','m ','h ','vx','vy','vz'/)
       endif
       if (maxvxyzu >= 4) then
          ilen = ilen + 1
          tags(ilen) = 'u'
       endif
       if (h2chemistry) then
          tags(ilen+1:ilen+5) = (/'h2ratio','abHIq  ','abhpq  ','abeq   ','abco   '/)
          ilen = ilen + 5
       endif
    case(i_real4)
       ilen = 1
       tags(ilen) = 'h'
       ilen = ilen + 1
       tags(ilen) = 'alpha'
    end select
 case(2) ! sink particle arrays
    if (ikind==i_real) then
       tags(1:4) = (/'x','y','z','m'/)
       tags(5:9) = (/'hsoft    ','maccreted','spinx    ','spiny    ','spinz    '/)
       tags(10:12) = (/'vx','vy','vz'/)
    endif
 case(4) ! MHD arrays
    if (ikind==i_real4) then
       tags = (/'Bx','By','Bz'/)
    endif
 end select

end subroutine fake_array_tags

!--------------------------------------------------------------------
!+
!  subroutine to write output to full dump file
!  in GADGET format
!+
!-------------------------------------------------------------------
subroutine write_gadgetdump(dumpfile,t,xyzh,particlemass,vxyzu,rho,utherm,npart)
 use io,       only:iprint,idump,real4
#ifdef PERIODIC
 use boundary, only:dxbound
#endif
 real,             intent(in) :: t,particlemass,utherm
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: npart
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: rho(:)

 integer(kind=4) :: particleid(size(rho))
 integer :: npartoftype(6),nall(6),ncrap(6)
 real(kind=8) :: massoftype(6)
 real(kind=8)                          :: time,boxsize
 real(kind=8), parameter               :: dumz = 0.d0
 real(kind=4) :: unused(15)
 integer, parameter :: iflagsfr = 0, iflagfeedback = 0, iflagcool = 0
 integer, parameter :: nfiles = 1
 integer            :: ierr,i,j
!
!--open dumpfile
!
 write(iprint,"(/,/,'-------->   TIME = ',f12.4,"// &
              "': full dump written to file ',a,'   <--------',/)")  t,trim(dumpfile)

 write(iprint,*) 'writing to unit ',idump
 open(unit=idump,file=dumpfile,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'error: can''t create new dumpfile ',trim(dumpfile)
    stop
 endif

 npartoftype(:) = 0
 npartoftype(1) = npart
 nall(:)  = npartoftype(:)
 ncrap(:) = 0
 time     = t
#ifdef PERIODIC
 boxsize = dxbound
#else
 boxsize = 0.
#endif

 massoftype(:) = 0.
 massoftype(1) = particlemass
 unused(:) = 0

 do i=1,npart
    particleid(i) = i
 enddo
 write(idump,iostat=ierr) npartoftype(1:6),massoftype(1:6),time,dumz, &
                          iflagsfr,iflagfeedback,nall(1:6),iflagcool,nfiles,boxsize, &
                          dumz,dumz,dumz,iflagsfr,iflagsfr,ncrap(1:6),iflagsfr,unused(:)

 write(idump,iostat=ierr) ((real4(xyzh(j,i)),j=1,3),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing positions'
    return
 endif
 write(idump,iostat=ierr) ((real4(vxyzu(j,i)),j=1,3),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing velocities'
    return
 endif
 write(idump,iostat=ierr) (particleid(i),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing particle ID'
    return
 endif
 if (size(vxyzu(:,1)) >= 4) then
    write(idump,iostat=ierr) (real4(vxyzu(4,i)),i=1,npart)
 else
    write(idump,iostat=ierr) (real4(utherm),i=1,npart)
 endif
 if (ierr /= 0) then
    print*,' error writing utherm'
    return
 endif
 write(idump,iostat=ierr) (real4(rho(i)),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing rho'
    return
 endif
 write(idump,iostat=ierr) (real4(xyzh(4,i)),i=1,npart)
 if (ierr /= 0) then
    print*,' error writing h'
    return
 endif
 print*,' finished writing file -- OK'

 return
end subroutine write_gadgetdump

subroutine count_particle_types(npartoftype)
 use part, only:iphase,iamtype,npart
 integer, intent(out) :: npartoftype(:)
 integer :: i, itype

 npartoftype(:) = 0
 do i = 1, npart
    itype = iamtype(iphase(i))
    npartoftype(itype) = npartoftype(itype) + 1
 enddo
end subroutine count_particle_types

end module readwrite_dumps
