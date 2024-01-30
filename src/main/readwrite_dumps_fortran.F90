!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module readwrite_dumps_fortran
!
! This module contains all routines related
!  to the data format.
!
!  For Phantom, the format is identical to sphNG
!  (although with fewer arrays dumped)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, boundary_dyn, checkconserved, dim, dump_utils,
!   dust, dust_formation, eos, externalforces, fileutils, io, lumin_nsdisc,
!   memory, metric_tools, mpi, mpiutils, options, part,
!   readwrite_dumps_common, setup_params, sphNGutils, timestep, units
!
 use dump_utils, only:lenid,ndatatypes,i_int,i_int1,i_int2,i_int4,i_int8,&
                      i_real,i_real4,i_real8,int1,int2,int1o,int2o,dump_h,lentag
 use readwrite_dumps_common, only:check_arrays,fileident,get_options_from_fileid
 implicit none

 public :: write_smalldump_fortran,write_fulldump_fortran,read_smalldump_fortran,read_dump_fortran

 logical, target, public    :: opened_full_dump_fortran       ! for use in analysis files if user wishes to skip small dumps
 logical, target, public    :: dt_read_in_fortran             ! to determine if dt has been read in so that ibin & ibinold can be set on restarts
 integer, parameter :: maxphead = 256         ! max items in header
 integer, parameter :: is_small_dump = 1978
 integer, parameter :: is_not_mhd = 1979

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

end subroutine end_threadwrite

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
subroutine write_fulldump_fortran(t,dumpfile,ntotal,iorder,sphNG)
 use dim,   only:maxp,maxvxyzu,maxalpha,ndivcurlv,ndivcurlB,maxgrav,gravity,use_dust,&
                 lightcurve,use_dustgrowth,store_dust_temperature,gr,do_nucleation,&
                 ind_timesteps,mhd_nonideal,use_krome,h2chemistry,update_muGamma
 use eos,   only:ieos,eos_is_non_ideal,eos_outputs_mu,eos_outputs_gasP
 use io,    only:idump,iprint,real4,id,master,error,warning,nprocs
 use part,  only:xyzh,xyzh_label,vxyzu,vxyzu_label,Bevol,Bevol_label,Bxyz,Bxyz_label,npart,maxtypes, &
                 npartoftypetot,update_npartoftypetot, &
                 alphaind,rhoh,divBsymm,maxphase,iphase,iamtype_int1,iamtype_int11, &
                 nptmass,nsinkproperties,xyzmh_ptmass,xyzmh_ptmass_label,vxyz_ptmass,vxyz_ptmass_label,&
                 maxptmass,get_pmass,nabundances,abundance,abundance_label,mhd,&
                 divcurlv,divcurlv_label,divcurlB,divcurlB_label,poten,dustfrac,deltav,deltav_label,tstop,&
                 dustfrac_label,tstop_label,dustprop,dustprop_label,eos_vars,eos_vars_label,ndusttypes,ndustsmall,VrelVf,&
                 VrelVf_label,dustgasprop,dustgasprop_label,dust_temp,pxyzu,pxyzu_label,dens,& !,dvdx,dvdx_label
                 rad,rad_label,radprop,radprop_label,do_radiation,maxirad,maxradprop,itemp,igasP,igamma,&
                 iorig,iX,iZ,imu,nucleation,nucleation_label,n_nucleation,tau,itau_alloc,tau_lucy,itauL_alloc,&
                 luminosity,eta_nimhd,eta_nimhd_label
 use part,  only:metrics,metricderivs,tmunus
 use options,    only:use_dustfrac,use_var_comp,icooling
 use dump_utils, only:tag,open_dumpfile_w,allocate_header,&
                 free_header,write_header,write_array,write_block_header
 use mpiutils,   only:reduce_mpi,reduceall_mpi
 use timestep,   only:dtmax,idtmax_n,idtmax_frac
 use part,       only:ibin,krome_nmols,T_gas_cool
#ifdef PRDRAG
 use lumin_nsdisc, only:beta
#endif
 use metric_tools, only:imetric, imet_et
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in), optional :: iorder(:)
 logical,          intent(in), optional :: sphNG
 integer(kind=8),  intent(in), optional :: ntotal

 integer, parameter :: isteps_sphNG = 0, iphase0 = 0
 integer(kind=8)    :: ilen(4)
 integer            :: nums(ndatatypes,4)
 integer            :: ipass,k,l,ioffset
 integer            :: ierr,ierrs(30)
 integer            :: nblocks,nblockarrays,narraylengths
 integer(kind=8)    :: nparttot
 logical            :: sphNGdump,write_itype,use_gas
 character(len=lenid)  :: fileid
 character(len=120)    :: blankarray
 type(dump_h)          :: hdr
 real, allocatable :: temparr(:)
!
!--collect global information from MPI threads
!
!--allow non-MPI calls to create MPI dump files
#ifdef MPI
 nparttot = reduceall_mpi('+',npart)
 call update_npartoftypetot
#else
 if (present(ntotal)) then
    nparttot = ntotal
    call update_npartoftypetot
    if (all(npartoftypetot==0)) then
       npartoftypetot(1) = ntotal
    endif
 else
    nparttot = npart
    call update_npartoftypetot
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

    if (idtmax_frac==0) then
       write(iprint,"(/,/,'-------->   TIME = ',g12.4,': full dump written to file ',a,'   <--------',/)")  t,trim(dumpfile)
    else
       ioffset = max(0,len(trim(dumpfile))-1)
       write(blankarray,'(a)') ' '
       write(iprint,"(/,/,'-------->   TIME = ',g12.4,': full dump written to file ',a,'   <--------')")  &
       t,trim(dumpfile)
       write(iprint,"('-------->                        Writing sub-dumps: ',I4,' of',I4,a,'<--------',/)")  &
       idtmax_frac,idtmax_n,blankarray(1:ioffset)
    endif
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
          call write_array(1,dustprop,dustprop_label,2,npart,k,ipass,idump,nums,ierrs(3))
          call write_array(1,VrelVf,VrelVf_label,npart,k,ipass,idump,nums,ierrs(3))
          call write_array(1,dustgasprop,dustgasprop_label,4,npart,k,ipass,idump,nums,ierrs(3))
       endif
       if (h2chemistry)  call write_array(1,abundance,abundance_label,nabundances,npart,k,ipass,idump,nums,ierrs(5))
       if (use_dust) call write_array(1,dustfrac,dustfrac_label,ndusttypes,npart,k,ipass,idump,nums,ierrs(7))
       if (use_dust) call write_array(1,tstop,tstop_label,ndustsmall,npart,k,ipass,idump,nums,ierrs(8))
       if (use_dustfrac) then
          do l=1,ndustsmall
             call write_array(1,deltav(:,l,:),deltav_label,3,npart,k,ipass,idump,nums,ierrs(10))
          enddo
       endif
       if (gr) then
          call write_array(1,pxyzu,pxyzu_label,maxvxyzu,npart,k,ipass,idump,nums,ierrs(8))
          call write_array(1,dens,'dens prim',npart,k,ipass,idump,nums,ierrs(8))
          if (imetric==imet_et) then
             ! Output metric if imetric=iet
             call write_array(1,metrics(1,1,1,:), 'gtt (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             !  call write_array(1,metrics(1,2,1,:), 'gtx (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             !  call write_array(1,metrics(1,3,1,:), 'gty (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             !  call write_array(1,metrics(1,2,1,:), 'gtz (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             !  call write_array(1,metrics(1,2,1,:), 'gtx (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             call write_array(1,metrics(2,2,1,:), 'gxx (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             call write_array(1,metrics(3,3,1,:), 'gyy (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             call write_array(1,metrics(4,4,1,:), 'gzz (covariant)',npart,k,ipass,idump,nums,ierrs(8))

             call write_array(1,metricderivs(1,1,1,:), 'dxgtt (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             call write_array(1,metricderivs(2,2,1,:), 'dxgxx (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             call write_array(1,metricderivs(3,3,1,:), 'dxgyy (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             call write_array(1,metricderivs(4,4,1,:), 'dxgzz (covariant)',npart,k,ipass,idump,nums,ierrs(8))

             call write_array(1,tmunus(1,1,:),  'tmunutt (covariant)',npart,k,ipass,idump,nums,ierrs(8))
          endif
       endif
       if (eos_is_non_ideal(ieos) .or. (.not.store_dust_temperature .and. icooling > 0)) then
          call write_array(1,eos_vars(itemp,:),eos_vars_label(itemp),npart,k,ipass,idump,nums,ierrs(12))
       endif
       if (eos_is_non_ideal(ieos)) call write_array(1,eos_vars(igamma,:),eos_vars_label(igamma),npart,k,ipass,idump,nums,ierrs(12))

       call write_array(1,vxyzu,vxyzu_label,maxvxyzu,npart,k,ipass,idump,nums,ierrs(4))
       ! write pressure to file
       if ((eos_outputs_gasP(ieos) .or. eos_is_non_ideal(ieos)) .and. k==i_real) then
          call write_array(1,eos_vars,eos_vars_label,1,npart,k,ipass,idump,nums,ierrs(13),index=igasP)
       endif
       ! write X, Z, mu to file
       if (eos_outputs_mu(ieos)) then
          call write_array(1,eos_vars,eos_vars_label,1,npart,k,ipass,idump,nums,ierrs(13),index=imu)
          if (use_var_comp) then
             call write_array(1,eos_vars,eos_vars_label,1,npart,k,ipass,idump,nums,ierrs(13),index=iX)
             call write_array(1,eos_vars,eos_vars_label,1,npart,k,ipass,idump,nums,ierrs(13),index=iZ)
          endif
       endif

       ! smoothing length written as real*4 to save disk space
       call write_array(1,xyzh,xyzh_label,1,npart,k,ipass,idump,nums,ierrs(14),use_kind=4,index=4)
       if (maxalpha==maxp) call write_array(1,alphaind,(/'alpha'/),1,npart,k,ipass,idump,nums,ierrs(15))
       !if (maxalpha==maxp) then ! (uncomment this to write alphaloc to the full dumps)
       !   call write_array(1,alphaind,(/'alpha ','alphaloc'/),2,npart,k,ipass,idump,nums,ierrs(10))
       !endif
       if (ndivcurlv >= 1) call write_array(1,divcurlv,divcurlv_label,ndivcurlv,npart,k,ipass,idump,nums,ierrs(16))
       !if (maxdvdx==maxp) call write_array(1,dvdx,dvdx_label,9,npart,k,ipass,idump,nums,ierrs(17))
       if (gravity .and. maxgrav==maxp) then
          call write_array(1,poten,'poten',npart,k,ipass,idump,nums,ierrs(17))
       endif
       if (ind_timesteps) then
          if (.not.allocated(temparr)) allocate(temparr(npart))
          temparr(1:npart) = dtmax/2.**ibin(1:npart)
          call write_array(1,temparr,'dt',npart,k,ipass,idump,nums,ierrs(18),use_kind=4)
       endif
       call write_array(1,iorig,'iorig',npart,k,ipass,idump,nums,ierrs(29))

#ifdef PRDRAG
       if (k==i_real) then
          if (.not.allocated(temparr)) allocate(temparr(npart))
          do l=1,npart
             temparr(l) = real4(beta(xyzh(1,l), xyzh(2,l), xyzh(3,l)))
          enddo
          call write_array(1,temparr,'beta_pr',npart,k,ipass,idump,nums,ierrs(19))
       endif
#endif
       if (lightcurve) then
          call write_array(1,luminosity,'luminosity',npart,k,ipass,idump,nums,ierrs(20))
       endif

       if (use_krome) then
          call write_array(1,abundance,abundance_label,krome_nmols,npart,k,ipass,idump,nums,ierrs(21))
          call write_array(1,T_gas_cool,'temp',npart,k,ipass,idump,nums,ierrs(24))
       endif
       if (update_muGamma .or. use_krome) then
          call write_array(1,eos_vars(imu,:),eos_vars_label(imu),npart,k,ipass,idump,nums,ierrs(12))
          call write_array(1,eos_vars(igamma,:),eos_vars_label(igamma),npart,k,ipass,idump,nums,ierrs(12))
       endif
       if (do_nucleation) then
          call write_array(1,nucleation,nucleation_label,n_nucleation,npart,k,ipass,idump,nums,ierrs(25))
       endif
       If (itau_alloc == 1) then
          call write_array(1,tau,'tau',npart,k,ipass,idump,nums,ierrs(30))
       endif
       If (itauL_alloc == 1) then
          call write_array(1,tau_lucy,'tau_lucy',npart,k,ipass,idump,nums,ierrs(30))
       endif
       if (store_dust_temperature) then
          call write_array(1,dust_temp,'Tdust',npart,k,ipass,idump,nums,ierrs(26))
       endif
       if (do_radiation) then
          call write_array(1,rad,rad_label,maxirad,npart,k,ipass,idump,nums,ierrs(27))
          call write_array(1,radprop,radprop_label,maxradprop,npart,k,ipass,idump,nums,ierrs(28))
       endif
       if (any(ierrs(1:28) /= 0)) call error('write_dump','error writing hydro arrays')
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
          call write_array(4,Bxyz,Bxyz_label,3,npart,k,ipass,idump,nums,ierrs(1)) ! Bx,By,Bz
          call write_array(4,Bevol,Bevol_label,1,npart,k,ipass,idump,nums,ierrs(1),index=4) ! psi
          if (ndivcurlB >= 1) then
             call write_array(4,divcurlB,divcurlB_label,ndivcurlB,npart,k,ipass,idump,nums,ierrs(2))
          else
             call write_array(4,divBsymm,'divBsymm',npart,k,ipass,idump,nums,ierrs(2))
          endif
          if (any(ierrs(1:2) /= 0)) call error('write_dump','error writing MHD arrays')
          if (mhd_nonideal) then
             call write_array(4,eta_nimhd,eta_nimhd_label,4,npart,k,ipass,idump,nums,ierrs(1))
             if (ierrs(1) /= 0) call error('write_dump','error writing non-ideal MHD arrays')
          endif
       endif
    enddo
    if (ipass==1) call write_block_header(narraylengths,ilen,nums,idump,ierr)
 enddo
 if (allocated(temparr)) deallocate(temparr)

 if (ierr /= 0) write(iprint,*) 'error whilst writing dumpfile '//trim(dumpfile)

 close(unit=idump)
 call end_threadwrite(id)

end subroutine write_fulldump_fortran

!--------------------------------------------------------------------
!+
!  subroutine to write output to small dump file
!  (ie. minimal output...)
!
!  note that small dumps are always SINGLE PRECISION
!  (faked to look like the default real is real*4)
!+
!-------------------------------------------------------------------

subroutine write_smalldump_fortran(t,dumpfile)
 use dim,        only:maxp,maxtypes,use_dust,lightcurve,use_dustgrowth,h2chemistry
 use io,         only:idump,iprint,real4,id,master,error,warning,nprocs
 use part,       only:xyzh,xyzh_label,npart,Bxyz,Bxyz_label,&
                      npartoftypetot,update_npartoftypetot,&
                      maxphase,iphase,nabundances,&
                      nptmass,nsinkproperties,xyzmh_ptmass,xyzmh_ptmass_label,&
                      abundance,abundance_label,mhd,dustfrac,iamtype_int11,&
                      dustprop,dustprop_label,dustfrac_label,ndusttypes,&
                      rad,rad_label,do_radiation,maxirad,luminosity
 use dump_utils, only:open_dumpfile_w,dump_h,allocate_header,free_header,&
                      write_header,write_array,write_block_header
 use mpiutils,   only:reduceall_mpi
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer(kind=8) :: ilen(4)
 integer         :: nums(ndatatypes,4)
 integer         :: ierr,ipass,k
 integer         :: nblocks,nblockarrays,narraylengths
 integer(kind=8) :: nparttot
 logical         :: write_itype
 type(dump_h)    :: hdr
!
!--collect global information from MPI threads
!
 nparttot = reduceall_mpi('+',npart)
 call update_npartoftypetot
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
          call write_array(1,dustprop,dustprop_label,2,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       endif
       if (h2chemistry .and. nabundances >= 1) &
          call write_array(1,abundance,abundance_label,1,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       if (use_dust) &
          call write_array(1,dustfrac,dustfrac_label,ndusttypes,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       call write_array(1,xyzh,xyzh_label,4,npart,k,ipass,idump,nums,ierr,index=4,use_kind=4)

       if (lightcurve) call write_array(1,luminosity,'luminosity',npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       if (do_radiation) call write_array(1,rad,rad_label,maxirad,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
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

end subroutine write_smalldump_fortran

!--------------------------------------------------------------------
!+
!  subroutine to read dump from file
!  needs to be able to read Phantom dumps as in write_fulldump
!  and also from standard sphNG dump files
!+
!-------------------------------------------------------------------

subroutine read_dump_fortran(dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc)
 use memory,   only:allocate_memory
 use dim,      only:maxp,maxvxyzu,gravity,lightcurve,mhd,maxp_hard,inject_parts,mpi
 use io,       only:real4,master,iverbose,error,warning ! do not allow calls to fatal in this routine
 use part,     only:xyzh,vxyzu,massoftype,npart,npartoftype,maxtypes,iphase, &
                    maxphase,isetphase,nptmass,nsinkproperties,maxptmass,get_pmass, &
                    xyzmh_ptmass,vxyz_ptmass
 use dump_utils,   only:skipblock,skip_arrays,check_tag,lenid,ndatatypes,read_header, &
                        open_dumpfile_r,get_error_text,ierr_realsize,free_header,read_block_header
 use mpiutils,     only:reduce_mpi,reduceall_mpi
 use sphNGutils,   only:convert_sinks_sphNG,mass_sphng
 use options,      only:use_dustfrac
 use boundary_dyn, only:dynamic_bdy
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
 logical               :: tagged,phantomdump,smalldump
 real                  :: dumr,alphafile
 character(len=lenid)  :: fileidentr
 type(dump_h)          :: hdr
 integer               :: i,ierrh

 if (id==master) write(iprint,"(/,1x,a,i3)") '>>> reading setup from file: '//trim(dumpfile)//' on unit ',idisk1
 opened_full_dump_fortran = .true.
 dt_read_in_fortran       = .false.
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
 call read_header(idisk1,hdr,ierr,tagged=tagged)
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
!--allocate main arrays
!
    if (iblock==1) then
       if (dynamic_bdy .or. inject_parts) then
          if (mpi) then
             call allocate_memory(max(nparttot,int(maxp_hard/nprocs,kind=8)))
          else
             call allocate_memory(max(nparttot,int(maxp_hard,kind=8)))
          endif
       else
          call allocate_memory(nparttot)
       endif
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

    if (.not. phantomdump) then
       print *, "allocating arrays for nptmass=", nptmass
       allocate(mass_sphng(maxp_hard))
    endif

    call read_phantom_arrays(i1,i2,noffset,narraylengths,nums,npartread,npartoftype,&
                          massoftype,nptmass,nsinkproperties,phantomdump,tagged,.false.,&
                          tfile,alphafile,idisk1,iprint,ierr)

    if (ierr /= 0) call warning('read_dump','error reading arrays from file')

 enddo overblocks

 !
 ! check npartoftype- breaks for sphng "type19" sinks
 ! post-itype conversion
 !
 if (maxphase==maxp .and. phantomdump) then
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
 endif

 !
 ! convert sinks from sphNG -> Phantom
 !
 if (.not.phantomdump .and. nptmass > 0 .and. maxphase==maxp) then
    call convert_sinks_sphNG(npart,nptmass,iphase,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,ierr)
 endif

 call check_npartoftype(npartoftype,npart)
 if (.not. phantomdump) then
    deallocate(mass_sphng)
 endif

 if (narraylengths >= 4) then
    if (id==master) write(iprint,"(a,/)") ' <<< finished reading (MHD) file '
 else
    if (id==master) write(iprint,"(a,/)") ' <<< finished reading (hydro) file '
 endif
 close(idisk1)
 return

100 close (idisk1)
 call check_npartoftype(npartoftype,npart)
 write(iprint,"(a,/)") ' <<< ERROR! end of file reached in data read'
 ierr = 666

end subroutine read_dump_fortran

!--------------------------------------------------------------------
!+
!  sanity check on npartoftype
!+
!-------------------------------------------------------------------
subroutine check_npartoftype(npartoftype,npart)
 integer, intent(inout) :: npartoftype(:)
 integer, intent(in)    :: npart

 if (sum(npartoftype)==0) then
    print*,'WARNING: npartoftype not set in file, ASSUMING ALL PARTICLES ARE GAS'
    npartoftype(1) = npart
 endif

end subroutine check_npartoftype

!--------------------------------------------------------------------
!+
!  subroutine to read a small dump from file, as written
!  in write_smalldump
!+
!-------------------------------------------------------------------
subroutine read_smalldump_fortran(dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc)
 use memory,   only:allocate_memory
 use dim,      only:maxvxyzu,mhd,maxphase,maxp
 use io,       only:real4,master,iverbose,error,warning ! do not allow calls to fatal in this routine
 use part,     only:npart,npartoftype,maxtypes,nptmass,nsinkproperties,maxptmass, &
                    massoftype
 use dump_utils,   only:skipblock,skip_arrays,check_tag,open_dumpfile_r,get_error_text,&
                        ierr_realsize,read_header,extract,free_header,read_block_header
 use mpiutils,     only:reduce_mpi,reduceall_mpi
 use options,      only:use_dustfrac
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
 logical               :: tagged,phantomdump,smalldump
 real                  :: alphafile
 character(len=lenid)  :: fileidentr
 type(dump_h)          :: hdr
 integer               :: i

 if (id==master) write(iprint,"(/,1x,a,i3)") '>>> reading small dump file: '//trim(dumpfile)//' on unit ',idisk1
 opened_full_dump_fortran = .false.
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
 call read_header(idisk1,hdr,ierr,singleprec=.true.,tagged=tagged)
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
 !--Allocate main arrays (no need for extra space here re: particle injection
 !  as small dumps are only read for visualisation/analysis purposes)
 !
 call allocate_memory(nparttot)
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
                          tfile,alphafile,idisk1,iprint,ierr)

    if (ierr /= 0) call warning('read_dump','error reading arrays from file')

 enddo overblocks

 !
 ! determine npartoftype
 !
 npartoftypetot = npartoftype
 if (maxphase==maxp) then
    call count_particle_types(npartoftype)
    npartoftypetotact = reduceall_mpi('+',npartoftype)
    do i = 1,maxtypes
       if (npartoftypetotact(i) /= npartoftypetot(i)) then
          write(*,*) 'npartoftypetot    =',npartoftypetot
          write(*,*) 'npartoftypetotact =',npartoftypetotact
          call error('read_dump','particle type counts do not match header')
       endif
    enddo
 endif

 call check_npartoftype(npartoftype,npart)
 if (narraylengths >= 4) then
    if (id==master) write(iprint,"(a,/)") ' <<< finished reading (MHD) file '
 else
    if (id==master) write(iprint,"(a,/)") ' <<< finished reading (hydro) file '
 endif
 close(idisk1)
 return

100 close (idisk1)
 call check_npartoftype(npartoftype,npart)
 write(iprint,"(a,/)") ' <<< ERROR! end of file reached in data read'
 ierr = 666
 return

end subroutine read_smalldump_fortran

!--------------------------------------------------------------------
!+
!  read arrays from the main block in the file into the relevant
!  phantom modules
!+
!-------------------------------------------------------------------
subroutine read_phantom_arrays(i1,i2,noffset,narraylengths,nums,npartread,npartoftype,&
                               massoftype,nptmass,nsinkproperties,phantomdump,tagged,singleprec,&
                               tfile,alphafile,idisk1,iprint,ierr)
 use dump_utils, only:read_array,match_tag
 use dim,        only:use_dust,h2chemistry,maxalpha,maxp,gravity,maxgrav,maxvxyzu,do_nucleation, &
                      use_dustgrowth,maxdusttypes,ndivcurlv,maxphase,gr,store_dust_temperature,&
                      ind_timesteps,use_krome
 use part,       only:xyzh,xyzh_label,vxyzu,vxyzu_label,dustfrac,dustfrac_label,abundance,abundance_label, &
                      alphaind,poten,xyzmh_ptmass,xyzmh_ptmass_label,vxyz_ptmass,vxyz_ptmass_label, &
                      Bevol,Bxyz,Bxyz_label,nabundances,iphase,idust, &
                      eos_vars,eos_vars_label,maxeosvars,dustprop,dustprop_label,divcurlv,divcurlv_label,iX,iZ,imu, &
                      VrelVf,VrelVf_label,dustgasprop,dustgasprop_label,pxyzu,pxyzu_label,dust_temp, &
                      rad,rad_label,radprop,radprop_label,do_radiation,maxirad,maxradprop,ifluxx,ifluxy,ifluxz, &
                      nucleation,nucleation_label,n_nucleation,ikappa,tau,itau_alloc,tau_lucy,itauL_alloc,&
                      ithick,ilambda,iorig,dt_in,krome_nmols,T_gas_cool
 use sphNGutils, only:mass_sphng,got_mass,set_gas_particle_mass
 integer, intent(in)   :: i1,i2,noffset,narraylengths,nums(:,:),npartread,npartoftype(:),idisk1,iprint
 real,    intent(in)   :: massoftype(:)
 integer, intent(in)   :: nptmass,nsinkproperties
 logical, intent(in)   :: phantomdump,singleprec,tagged
 real,    intent(in)   :: tfile,alphafile
 integer, intent(out)  :: ierr
 logical               :: match
 logical               :: got_dustfrac(maxdusttypes)
 logical               :: got_iphase,got_xyzh(4),got_vxyzu(4),got_abund(nabundances),got_alpha(1),got_poten
 logical               :: got_sink_data(nsinkproperties),got_sink_vels(3),got_Bxyz(3)
 logical               :: got_krome_mols(krome_nmols),got_krome_T,got_krome_gamma,got_krome_mu
 logical               :: got_eosvars(maxeosvars),got_nucleation(n_nucleation),got_ray_tracer
 logical               :: got_psi,got_Tdust,got_dustprop(2),got_VrelVf,got_dustgasprop(4)
 logical               :: got_divcurlv(4),got_rad(maxirad),got_radprop(maxradprop),got_pxyzu(4),got_iorig
 character(len=lentag) :: tag,tagarr(64)
 integer :: k,i,iarr,ik,ndustfraci

!
!--read array type 1 arrays
!
 got_iphase      = .false.
 got_xyzh        = .false.
 got_vxyzu       = .false.
 got_dustfrac    = .false.
 got_abund       = .false.
 got_alpha       = .false.
 got_poten       = .false.
 got_sink_data   = .false.
 got_sink_vels   = .false.
 got_Bxyz        = .false.
 got_psi         = .false.
 got_eosvars     = .false.
 got_dustprop    = .false.
 got_VrelVf      = .false.
 got_dustgasprop = .false.
 got_divcurlv    = .false.
 got_Tdust       = .false.
 got_krome_mols  = .false.
 got_krome_gamma = .false.
 got_krome_mu    = .false.
 got_krome_T     = .false.
 got_nucleation  = .false.
 got_ray_tracer  = .false.
 got_rad         = .false.
 got_radprop     = .false.
 got_pxyzu       = .false.
 got_iorig       = .false.

 ndustfraci = 0
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
             if (maxphase==maxp) call read_array(iphase,'itype',got_iphase,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             call read_array(xyzh, xyzh_label, got_xyzh, ik,i1,i2,noffset,idisk1,tag,match,ierr)
             call read_array(vxyzu,vxyzu_label,got_vxyzu,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             if (.not. phantomdump) then
                call read_array(iphase,'iphase',got_iphase,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                call read_array(mass_sphng,'m',got_mass,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (gr) call read_array(pxyzu,pxyzu_label,got_pxyzu,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             if (use_dustgrowth) then
                call read_array(dustprop,dustprop_label,got_dustprop,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                call read_array(VrelVf,VrelVf_label,got_VrelVf,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                call read_array(dustgasprop,dustgasprop_label,got_dustgasprop,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (use_dust) then
                if (any(tag == dustfrac_label)) then
                   ndustfraci = ndustfraci + 1
                   call read_array(dustfrac(ndustfraci,:),dustfrac_label(ndustfraci),got_dustfrac(ndustfraci), &
                                   ik,i1,i2,noffset,idisk1,tag,match,ierr)
                endif
             endif
             if (h2chemistry) then
                call read_array(abundance,abundance_label,got_abund,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (use_krome) then
                call read_array(abundance,abundance_label,got_krome_mols,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                call read_array(T_gas_cool,'temp',got_krome_T,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (do_nucleation) then
                call read_array(nucleation,nucleation_label,got_nucleation,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (itau_alloc == 1) then
                call read_array(tau,'tau',got_ray_tracer,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (itauL_alloc == 1) then
                call read_array(tau_lucy,'tau_lucy',got_ray_tracer,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (store_dust_temperature) then
                call read_array(dust_temp,'Tdust',got_Tdust,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             call read_array(eos_vars,eos_vars_label,got_eosvars,ik,i1,i2,noffset,idisk1,tag,match,ierr)

             if (maxalpha==maxp) call read_array(alphaind,(/'alpha'/),got_alpha,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             !
             ! read divcurlv if it is in the file
             !
             if (ndivcurlv >= 1) call read_array(divcurlv,divcurlv_label,got_divcurlv,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             !
             ! read gravitational potential if it is in the file
             !
             if (gravity .and. maxgrav==maxp) call read_array(poten,'poten',got_poten,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             !
             ! read dt if it is in the file
             !
             if (ind_timesteps) call read_array(dt_in,'dt',dt_read_in_fortran,ik,i1,i2,noffset,idisk1,tag,match,ierr)

             ! read particle ID's
             call read_array(iorig,'iorig',got_iorig,ik,i1,i2,noffset,idisk1,tag,match,ierr)

             if (do_radiation) then
                call read_array(rad,rad_label,got_rad,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                call read_array(radprop,radprop_label,got_radprop,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
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
 call check_arrays(i1,i2,noffset,npartoftype,npartread,nptmass,nsinkproperties,massoftype,&
                   alphafile,tfile,phantomdump,got_iphase,got_xyzh,got_vxyzu,got_alpha, &
                   got_krome_mols,got_krome_gamma,got_krome_mu,got_krome_T, &
                   got_abund,got_dustfrac,got_sink_data,got_sink_vels,got_Bxyz,got_psi,got_dustprop,got_pxyzu,got_VrelVf, &
                   got_dustgasprop,got_rad,got_radprop,got_Tdust,got_eosvars,got_nucleation,got_iorig,iphase,&
                   xyzh,vxyzu,pxyzu,alphaind,xyzmh_ptmass,Bevol,iorig,iprint,ierr)
 if (.not. phantomdump) then
    print *, "Calling set_gas_particle_mass"
    call set_gas_particle_mass(mass_sphng)
 endif
 return
100 continue
 write(iprint,"(a,/)") ' <<< ERROR! end of file reached in data read'

end subroutine read_phantom_arrays

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

!--------------------------------------------------------------------
!+
!  utility to extract header variables to phantom
!+
!-------------------------------------------------------------------
subroutine unfill_header(hdr,phantomdump,got_tags,nparttot, &
                         nblocks,npart,npartoftype, &
                         tfile,hfactfile,alphafile,iprint,id,nprocs,ierr)
 use dim,        only:maxdustlarge,use_dust
 use io,         only:master ! check this
 use eos,        only:isink
 use part,       only:maxtypes,igas,idust,ndustsmall,ndustlarge,ndusttypes,&
                      npartoftypetot
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
 integer         :: nparttoti,npartoftypetoti(maxtypes),ntypesinfile,nptinfile
 integer         :: ierr1,ierrs(3),i,counter
 integer(kind=8) :: ntypesinfile8
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
    write(*,*) 'READING ONLY FIRST ',maxtypes,' OF ',ntypesinfile,' particle types'
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
    if (nparttoti > huge(npart)) then
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
 endif
 if (ntypesinfile > maxtypes) then
    write(*,*) 'WARNING: number of particle types in file exceeds array limits'
    write(*,*) 'READING ONLY FIRST ',maxtypes,' OF ',ntypesinfile,' particle types'
    ntypesinfile = maxtypes
 endif
 call extract('nparttot',nparttot,hdr,ierr1)
 if (nblocks > 1) then
    call extract('npartoftype',npartoftype(1:ntypesinfile),hdr,ierr1)
 endif
 if (id==master) write(*,*) 'npart(total) = ',nparttot
!
!--number of dust species
!
 if (use_dust) then
    call extract('ndustsmall',ndustsmall,hdr,ierrs(1))
    call extract('ndustlarge',ndustlarge,hdr,ierrs(2))
    if (any(ierrs(1:2) /= 0)) then
       call extract('ndustfluids',ndustsmall,hdr,ierrs(1)) ! for backwards compatibility
       if (ierrs(1) /= 0) write(*,*) 'ERROR reading number of small/large grain types from file header'
    endif
    ndusttypes = ndustsmall + ndustlarge
 endif
!
!--units
!
 call extract('udist',udist,hdr,ierrs(1))
 call extract('umass',umass,hdr,ierrs(2))
 call extract('utime',utime,hdr,ierrs(3))
 if (all(ierrs(1:3)==0)) then
    call set_units_extra()
 else
    write(iprint,*) 'ERROR reading units from dump file, assuming default'
    call set_units()  ! use default units
 endif
 ! get nptmass from header, needed to figure out if gwinspiral info is sensible
 call extract('nptmass',nptinfile,hdr,ierrs(1))
!--default real
 call unfill_rheader(hdr,phantomdump,ntypesinfile,nptinfile,&
                     tfile,hfactfile,alphafile,iprint,ierr)
 if (ierr /= 0) return

 if (id==master) write(iprint,*) 'time = ',tfile

end subroutine unfill_header

!--------------------------------------------------------------------
!+
!  subroutine to fill the real header with various things
!+
!-------------------------------------------------------------------
subroutine fill_header(sphNGdump,t,nparttot,npartoftypetot,nblocks,nptmass,hdr,ierr)
 use eos,            only:write_headeropts_eos,polyk2
 use options,        only:tolh,alpha,alphau,alphaB,iexternalforce,ieos
 use part,           only:massoftype,hfact,Bextx,Bexty,Bextz,ndustsmall,ndustlarge,&
                          idust,grainsize,graindens,ndusttypes
 use checkconserved, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use setup_params,   only:rhozero
 use timestep,       only:dtmax_user,idtmax_n_next,idtmax_frac_next,C_cour,C_force
 use externalforces, only:write_headeropts_extern
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax
 use boundary_dyn,   only:dynamic_bdy,dxyz,rho_bkg_ini,irho_bkg_ini
 use dump_utils,     only:reset_header,add_to_rheader,add_to_header,add_to_iheader,num_in_header
 use dim,            only:use_dust,maxtypes,use_dustgrowth,do_nucleation, &
                          phantom_version_major,phantom_version_minor,phantom_version_micro,periodic,idumpfile
 use units,          only:udist,umass,utime,unit_Bfield
 use dust_formation, only:write_headeropts_dust_formation

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
 call add_to_iheader(nptmass,'nptmass',hdr,ierr)
 call add_to_iheader(ndustlarge,'ndustlarge',hdr,ierr)
 call add_to_iheader(ndustsmall,'ndustsmall',hdr,ierr)
 call add_to_iheader(idust,'idust',hdr,ierr)
 call add_to_iheader(idtmax_n_next,'idtmax_n',hdr,ierr)
 call add_to_iheader(idtmax_frac_next,'idtmax_frac',hdr,ierr)
 call add_to_iheader(idumpfile,'idumpfile',hdr,ierr)
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
 call write_headeropts_eos(ieos,hdr,ierr)

 ! default real variables
 call add_to_rheader(t,'time',hdr,ierr)
 call add_to_rheader(dtmax_user,'dtmax',hdr,ierr)
 call add_to_rheader(rhozero,'rhozero',hdr,ierr)
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
    call add_to_rheader(massoftype,'massoftype',hdr,ierr) ! array
    if (do_nucleation) call write_headeropts_dust_formation(hdr,ierr)
    call add_to_rheader(Bextx,'Bextx',hdr,ierr)
    call add_to_rheader(Bexty,'Bexty',hdr,ierr)
    call add_to_rheader(Bextz,'Bextz',hdr,ierr)
    call add_to_rheader(0.,'dum',hdr,ierr)
    if (iexternalforce /= 0) call write_headeropts_extern(iexternalforce,hdr,t,ierr)
    if (periodic) then
       call add_to_rheader(xmin,'xmin',hdr,ierr)
       call add_to_rheader(xmax,'xmax',hdr,ierr)
       call add_to_rheader(ymin,'ymin',hdr,ierr)
       call add_to_rheader(ymax,'ymax',hdr,ierr)
       call add_to_rheader(zmin,'zmin',hdr,ierr)
       call add_to_rheader(zmax,'zmax',hdr,ierr)
    endif
    if (dynamic_bdy) then
       call add_to_rheader(dxyz,'dxyz',hdr,ierr)
       call add_to_iheader(irho_bkg_ini,'irho_bkg_ini',hdr,ierr)
       call add_to_rheader(rho_bkg_ini,'rho_bkg_ini',hdr,ierr)
    endif
    call add_to_rheader(get_conserv,'get_conserv',hdr,ierr)
    call add_to_rheader(etot_in,'etot_in',hdr,ierr)
    call add_to_rheader(angtot_in,'angtot_in',hdr,ierr)
    call add_to_rheader(totmom_in,'totmom_in',hdr,ierr)
    call add_to_rheader(mdust_in(1:ndusttypes),'mdust_in',hdr,ierr)
    if (use_dust) then
       call add_to_rheader(grainsize(1:ndusttypes),'grainsize',hdr,ierr)
       call add_to_rheader(graindens(1:ndusttypes),'graindens',hdr,ierr)
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
subroutine unfill_rheader(hdr,phantomdump,ntypesinfile,nptmass,&
                          tfile,hfactfile,alphafile,iprint,ierr)
 use io,             only:id,master
 use dim,            only:maxvxyzu,nElements,use_dust,use_dustgrowth,use_krome,do_nucleation,idumpfile
 use eos,            only:extract_eos_from_hdr, read_headeropts_eos
 use options,        only:ieos,iexternalforce
 use part,           only:massoftype,Bextx,Bexty,Bextz,mhd,periodic,&
                          maxtypes,grainsize,graindens,ndusttypes
 use checkconserved, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use setup_params,   only:rhozero
 use externalforces, only:read_headeropts_extern,extract_iextern_from_hdr
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax,set_boundary
 use boundary_dyn,   only:dynamic_bdy,dxyz,irho_bkg_ini,rho_bkg_ini,rho_bkg_ini1
 use dump_utils,     only:extract
 use dust,           only:grainsizecgs,graindenscgs
 use units,          only:unit_density,udist
 use timestep,       only:idtmax_n,idtmax_frac
 use dust_formation, only:read_headeropts_dust_formation
 type(dump_h), intent(in)  :: hdr
 logical,      intent(in)  :: phantomdump
 integer,      intent(in)  :: iprint,ntypesinfile,nptmass
 real,         intent(out) :: tfile,hfactfile,alphafile
 integer,      intent(out) :: ierr

 integer, parameter :: lu = 173
 integer            :: ierrs(10),iextern_in_file
 real               :: xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,dtmaxi
 real               :: alphaufile,alphaBfile,C_courfile,C_forcefile,tolhfile
 logical            :: iexist

 ierr  = 0
 call extract('time',tfile,hdr,ierr)
 if (ierr/=0)  call extract('gt',tfile,hdr,ierr)  ! this is sphNG's label for time
 call extract('dtmax',dtmaxi,hdr,ierr)
 call extract('rhozero',rhozero,hdr,ierr)
 Bextx = 0.
 Bexty = 0.
 Bextz = 0.
 if (phantomdump) then
    call extract('hfact',hfactfile,hdr,ierr)
    call extract('tolh',tolhfile,hdr,ierr)
    call extract('C_cour',C_courfile,hdr,ierr)
    call extract('C_force',C_forcefile,hdr,ierr)
    call extract('alpha',alphafile,hdr,ierr)
    if (maxvxyzu >= 4) then
       call extract('alphau',alphaufile,hdr,ierr)
    else
       alphaufile = 0.
    endif
    if (mhd) then
       call extract('alphaB',alphaBfile,hdr,ierr)
    endif

    if (extract_eos_from_hdr) call extract('ieos',ieos,hdr,ierr)

    call extract('massoftype',massoftype(1:ntypesinfile),hdr,ierr)
    if (ierr /= 0) then
       write(*,*) '*** ERROR reading massoftype from dump header ***'
       ierr = 4
    endif
    if (do_nucleation) then
       call read_headeropts_dust_formation(hdr,ierr)
       if (ierr /= 0) ierr = 6
    endif

    call extract('iexternalforce',iextern_in_file,hdr,ierrs(1))
    if (extract_iextern_from_hdr) iexternalforce = iextern_in_file
    if (iexternalforce /= 0) then
       call read_headeropts_extern(iexternalforce,hdr,nptmass,ierrs(1))
       if (ierrs(1) /= 0) ierr = 5
    elseif (iextern_in_file /= 0) then
       call read_headeropts_extern(iextern_in_file,hdr,nptmass,ierrs(1))
       if (ierrs(1) /= 0) ierr = 5
    endif

    call extract('idtmax_n',idtmax_n,hdr,ierr,default=1)
    call extract('idtmax_frac',idtmax_frac,hdr,ierr)
    call extract('idumpfile',idumpfile,hdr,ierr)
 else
    massoftype(1) = 0.
    hfactfile = 0.
 endif

 call read_headeropts_eos(ieos,hdr,ierr)

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

 if (dynamic_bdy) then
    call extract('irho_bkg_ini',irho_bkg_ini,hdr,ierrs(1))
    call extract('rho_bkg_ini',rho_bkg_ini,hdr,ierrs(1))
    call extract('dxyz',dxyz,hdr,ierrs(2))
    if (rho_bkg_ini > 0.) then
       rho_bkg_ini1 = 1.0/rho_bkg_ini
    else
       rho_bkg_ini1 = 0.
    endif
 endif

 if (mhd) then
    call extract('Bextx',Bextx,hdr,ierrs(1))
    call extract('Bexty',Bexty,hdr,ierrs(2))
    call extract('Bextz',Bextz,hdr,ierrs(3))
    if (id==master) then
       if (any(ierrs(1:3) /= 0)) then
          write(*,*) 'ERROR reading external field (setting to zero)'
       else
          write(*,*) 'External field found, Bext = ',Bextx,Bexty,Bextz
       endif
    endif
 endif

 ! values to track that conserved values remain conserved
 call extract('get_conserv',get_conserv,hdr,ierrs(1))
 call extract('etot_in',    etot_in,    hdr,ierrs(2))
 call extract('angtot_in',  angtot_in,  hdr,ierrs(3))
 call extract('totmom_in',  totmom_in,  hdr,ierrs(4))
 call extract('mdust_in',   mdust_in(1:ndusttypes), hdr,ierrs(5))
 if (any(ierrs(1:4) /= 0)) then
    write(*,*) 'ERROR reading values to verify conservation laws.  Resetting initial values.'
    get_conserv = 1.0
 endif


 !--pull grain size and density arrays if they are in the header
 !-- i.e. if dustgrowth is not ON
 if (use_dust .and. .not.use_dustgrowth) then
    call extract('grainsize',grainsize(1:ndusttypes),hdr,ierrs(1))
    call extract('graindens',graindens(1:ndusttypes),hdr,ierrs(2))
    if (any(ierrs(1:2) /= 0)) then
       write(*,*) 'ERROR reading grain size/density from file header'
       grainsize(1) = real(grainsizecgs/udist)
       graindens(1) = real(graindenscgs/unit_density)
    endif
 endif

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

end module readwrite_dumps_fortran
