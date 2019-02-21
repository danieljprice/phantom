!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!    gitinfo, initial_params, io, lumin_nsdisc, memory, mpi, mpiutils,
!    options, part, setup_params, sphNGutils, timestep, units
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
subroutine write_fulldump(t,dumpfile,ntotal)
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer(kind=8),  intent(in), optional :: ntotal
 if (present(ntotal)) then
    call write_dump(t,dumpfile,fulldump=.true.,ntotal=ntotal)
 else
    call write_dump(t,dumpfile,fulldump=.true.)
 endif
end subroutine

!--------------------------------------------------------------------
!+
!  subroutine to write output to small dump file
!  (ie. minimal output...)
!
!  note that small dumps are always SINGLE PRECISION
!+
!-------------------------------------------------------------------
subroutine write_smalldump(t,dumpfile)
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 call write_dump(t,dumpfile,fulldump=.false.)
end subroutine

! Generic subroutine for writing a dumpfile
subroutine write_dump(t,dumpfile,fulldump,ntotal)
 use output_hdf5,    only:open_hdf5file,close_hdf5file,outputfile_id
 use output_hdf5,    only:write_hdf5_header,write_hdf5_arrays,write_hdf5_arrays_small
 use dim,            only:maxp,maxvxyzu,gravity,maxalpha,mhd,mhd_nonideal,use_dust,use_dustgrowth
 use dim,            only:phantom_version_major,phantom_version_minor,phantom_version_micro,store_temperature
 use dim,            only:phantom_version_string
 use gitinfo,        only:gitsha
 use eos,            only:ieos,equationofstate,done_init_eos,init_eos,polyk,gamma,polyk2,qfacdisc,isink
 use io,             only:nprocs,fatal,id,master,iprint
 use options,        only:tolh,alpha,alphau,alphaB,iexternalforce,use_dustfrac
 use part,           only:xyzh,vxyzu,Bevol,Bxyz,npart,npartoftype,maxtypes, &
                          alphaind,rhoh,divBsymm,maxphase,iphase, &
                          nptmass,xyzmh_ptmass,vxyz_ptmass,&
                          get_pmass,abundance,&
                          divcurlv,divcurlB,poten,dustfrac,deltav,tstop,&
                          dustprop,temperature,St,ndustsmall,luminosity,&
                          eta_nimhd,massoftype,hfact,Bextx,Bexty,Bextz,&
                          ndustlarge,idust,grainsize,graindens,&
                          h2chemistry,lightcurve,maxBevol,&
                          ndivcurlB,ndivcurlv,ndusttypes
#ifdef IND_TIMESTEPS
 use part,           only:ibin
#endif
 use mpiutils,       only:reduce_mpi,reduceall_mpi
 use lumin_nsdisc,   only:beta
 use initial_params, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use setup_params,   only:rhozero
 use timestep,       only:dtmax,C_cour,C_force
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax
 use units,          only:udist,umass,utime,unit_Bfield
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 logical,          intent(in) :: fulldump
 integer(kind=8),  intent(in), optional :: ntotal
 integer           :: i
 integer           :: ierr,nblocks
 integer(kind=8)   :: nparttot,npartoftypetot(maxtypes)
 logical           :: use_gas,ind_timesteps,const_av,prdrag,isothermal
 real              :: ponrhoi,rhoi,spsoundi
 real, allocatable, dimension(:) :: pressure,dtind,beta_pr
 character(len=100):: fileident
 character(len=10) :: datestring, timestring
 character(len=30) :: string
 character(len=9)  :: dumptype
 integer :: error

 if (id==master) then
    if (fulldump) then
       write(iprint,"(/,/,'-------->   TIME = ',g12.4,"//"': full dump written to file ',a,'   <--------',/)")  t,trim(dumpfile)
    else
       write(iprint,"(/,/,'-------->   TIME = ',g12.4,"//"': small dump written to file ',a,'   <--------',/)")  t,trim(dumpfile)
    endif
 endif

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

#ifdef IND_TIMESTEPS
 ind_timesteps = .true.
#else
 ind_timesteps = .false.
#endif
#ifdef PRDRAG
 prdrag = .true.
#else
 prdrag = .false.
#endif
#ifdef ISOTHERMAL
 isothermal = .true.
#else
 isothermal = .false.
#endif

 if (maxphase==maxp) then
    use_gas = .false.
 else
    use_gas = .true.
 endif

 if (fulldump) then
    allocate(pressure(nparttot),beta_pr(nparttot),dtind(nparttot))

    ! Compute pressure and beta_pr array
    if (.not.done_init_eos) call init_eos(ieos,ierr)
    !$omp parallel do default(none) &
    !$omp shared(xyzh,vxyzu,ieos,nparttot,pressure,beta_pr,temperature,use_gas,prdrag) &
    !$omp private(i,ponrhoi,spsoundi,rhoi)
    do i=1,int(nparttot)
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
       pressure(i) = ponrhoi*rhoi
       if (prdrag) beta_pr(i)  = beta(xyzh(1,i), xyzh(2,i), xyzh(3,i))
    enddo
    !$omp end parallel do

    ! Compute dtind array
    if (ind_timesteps) dtind = dtmax/2**ibin(1:npart)
 endif

! Check if constant AV
 if (maxp==maxalpha) then
    const_av = .false.
 else
    const_av = .true.
 endif

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

 if (fulldump) then
    dumptype = 'fulldump '
 else
    dumptype = 'smalldump'
 endif
 fileident = trim(dumptype)//': '//'Phantom'//' '//trim(phantom_version_string)//' '//gitsha

 if (mhd) then
    if (maxBevol==4) then
       fileident = trim(fileident)//' (mhd+clean'//trim(string)//')  : '//trim(datestring)//' '//trim(timestring)
    else
       fileident = trim(fileident)//' (mhd'//trim(string)//')  : '//trim(datestring)//' '//trim(timestring)
    endif
 else
    fileident = trim(fileident)//' (hydro'//trim(string)//'): '//trim(datestring)//' '//trim(timestring)
 endif

 call open_hdf5file(trim(dumpfile)//'.h5',outputfile_id,error)
 if (error/=0) call fatal('write_fulldump_hdf5','could not open file')
 call write_hdf5_header(outputfile_id,error,trim(fileident),maxtypes,nblocks,isink,nptmass,ndustlarge,ndustsmall,idust,&
                        phantom_version_major,phantom_version_minor,phantom_version_micro,                             &
                        int(nparttot),int(npartoftypetot),iexternalforce,ieos,t,dtmax,gamma,rhozero,                   &
                        polyk,hfact,tolh,C_cour,C_force,alpha,alphau,alphaB,polyk2,qfacdisc,                           &
                        massoftype,Bextx,Bexty,Bextz,xmin,xmax,ymin,ymax,zmin,zmax,get_conserv,                        &
                        etot_in,angtot_in,totmom_in,mdust_in,grainsize,graindens,udist,umass,utime,unit_Bfield         )
 if (error/=0) call fatal('write_fulldump_hdf5','could not write header')
 if (fulldump) then
    call write_hdf5_arrays(outputfile_id,error,xyzh,vxyzu,int(iphase),pressure,alphaind,dtind,poten,xyzmh_ptmass,            &
                           vxyz_ptmass,Bxyz,Bevol,divcurlB,divBsymm,eta_nimhd,                                               &
                           dustfrac(1:ndusttypes,:),tstop(1:ndustsmall,:),deltav(:,1:ndustsmall,:),dustprop,st,              &
                           abundance,temperature,divcurlv,luminosity,beta_pr,                                                &
                           const_av,ind_timesteps,gravity,nptmass,mhd,maxBevol,ndivcurlB,mhd_nonideal,use_dust,              &
                           use_dustfrac,use_dustgrowth,h2chemistry,store_temperature,ndivcurlv,lightcurve,prdrag,isothermal  )
 else
    call write_hdf5_arrays_small(outputfile_id,error,xyzh,int(iphase),xyzmh_ptmass,Bxyz,dustfrac,dustprop,st,   &
                                 abundance,luminosity,nptmass,mhd,use_dust,use_dustgrowth,h2chemistry,lightcurve)
 endif
 if (error/=0) call fatal('write_fulldump_hdf5','could not write arrays')
 call close_hdf5file(outputfile_id,error)
 if (error/=0) call fatal('write_fulldump_hdf5','could not close file')

 if (fulldump) deallocate(pressure,beta_pr,dtind)

end subroutine write_dump

!--------------------------------------------------------------------
!+
!  subroutine to read dump from file
!  needs to be able to read Phantom dumps as in write_fulldump
!  and also from standard sphNG dump files
!+
!-------------------------------------------------------------------

subroutine read_dump(dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc)
 use memory,   only:allocate_memory
 use dim,      only:maxp,maxvxyzu,gravity,lightcurve,mhd
#ifdef INJECT_PARTICLES
 use dim,      only:maxp_hard
#endif
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
 !--Allocate main arrays
 !
#ifdef INJECT_PARTICLES
 call allocate_memory(maxp_hard)
#else
 call allocate_memory(int(nparttot / nprocs))
#endif

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
 use memory,   only:allocate_memory
 use dim,      only:maxvxyzu,mhd,maxBevol
#ifdef MPI
 use dim,      only:maxp
#endif
#ifdef INJECT_PARTICLES
 use dim,      only:maxp_hard
#endif
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
 !--Allocate main arrays
 !
#ifdef INJECT_PARTICLES
 call allocate_memory(maxp_hard)
#else
 call allocate_memory(int(nparttot / nprocs))
#endif

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
                xyzmh_ptmass_label,vxyz_ptmass_label,get_pmass,rhoh,dustfrac,ndusttypes
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
 if (npartread > 0) then
    do i = 1, size(massoftype)
       if (npartoftype(i) > 0) then
          if (massoftype(i) <= 0.0) then
             if (id==master .and. i1==1) write(*,*) 'ERROR! mass not set in read_dump (Phantom)'
             ierr = 12
             return
          endif
       endif
    enddo
 endif
 if (use_dustfrac .and. .not. all(got_dustfrac(1:ndusttypes))) then
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
 use dim,        only:maxp_hard,maxdustlarge,use_dust
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
    if (nparttoti > maxp_hard) then
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
!  subroutine to set runtime parameters having read the real header
!+
!-------------------------------------------------------------------
subroutine unfill_rheader(hdr,phantomdump,ntypesinfile,&
                          tfile,hfactfile,alphafile,iprint,ierr)
 use io,             only:id,master
 use dim,            only:maxvxyzu,use_dust
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
 use dim, only:maxvxyzu,h2chemistry
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
