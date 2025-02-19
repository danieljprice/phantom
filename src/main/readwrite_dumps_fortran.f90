!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
! :Dependencies: boundary_dyn, dim, dump_utils, eos, eos_stamatellos, io,
!   memory, metric_tools, mpiutils, options, part, readwrite_dumps_common,
!   sphNGutils, timestep
!

! :Dependencies: boundary_dyn, dim, dump_utils, eos, eos_stamatellos, io,
!   memory, metric_tools, mpiutils, options, part, readwrite_dumps_common,
!   sphNGutils, timestep
!
 use dump_utils, only:lenid,ndatatypes,i_int,i_int1,i_int2,i_int4,i_int8,&
                      i_real,i_real4,i_real8,int1,int2,int1o,int2o,dump_h,lentag
 use readwrite_dumps_common, only:check_arrays,fileident,get_options_from_fileid,fill_header,unfill_header
 implicit none

 public :: write_smalldump_fortran,write_fulldump_fortran,read_smalldump_fortran,read_dump_fortran,unfill_header

 logical, target, public    :: opened_full_dump_fortran       ! for use in analysis files if user wishes to skip small dumps
 logical, target, public    :: dt_read_in_fortran             ! to determine if dt has been read in so that ibin & ibinold can be set on restarts
 integer, parameter :: maxphead = 256         ! max items in header
 integer, parameter :: is_small_dump = 1978
 integer, parameter :: is_not_mhd = 1979

 private

contains

!--------------------------------------------------------------------
!+
!  subroutine to write output to full dump file
!  (this is everything needed to restart a run)
!+
!-------------------------------------------------------------------
subroutine write_fulldump_fortran(t,dumpfile,ntotal,iorder,sphNG)
 use dim,   only:maxp,maxvxyzu,maxalpha,ndivcurlv,ndivcurlB,maxgrav,gravity,use_dust,&
                 lightcurve,use_dustgrowth,store_dust_temperature,gr,do_nucleation,&
                 ind_timesteps,mhd_nonideal,use_krome,h2chemistry,update_muGamma,mpi,use_apr,&
                 store_sf_ptmass
 use eos,   only:ieos,eos_is_non_ideal,eos_outputs_mu,eos_outputs_gasP
 use io,    only:idump,iprint,real4,id,master,error,warning,nprocs
 use part,  only:xyzh,xyzh_label,vxyzu,vxyzu_label,Bevol,Bevol_label,Bxyz,Bxyz_label,npart,maxtypes, &
                 npartoftypetot,update_npartoftypetot, &
                 alphaind,rhoh,divBsymm,maxphase,iphase,iamtype_int1,iamtype_int11, &
                 nptmass,nsinkproperties,xyzmh_ptmass,xyzmh_ptmass_label,vxyz_ptmass,vxyz_ptmass_label, sf_ptmass, &
                 sf_ptmass_label,maxptmass,get_pmass,nabundances,abundance,abundance_label,mhd,&
                 divcurlv,divcurlv_label,divcurlB,divcurlB_label,poten,dustfrac,deltav,deltav_label,tstop,&
                 dustfrac_label,tstop_label,dustprop,dustprop_label,eos_vars,eos_vars_label,ndusttypes,ndustsmall,VrelVf,&
                 VrelVf_label,dustgasprop,dustgasprop_label,filfac,filfac_label,dust_temp,pxyzu,pxyzu_label,dens,& !,dvdx,dvdx_label
                 rad,rad_label,radprop,radprop_label,do_radiation,maxirad,maxradprop,itemp,igasP,igamma,&
                 iorig,iX,iZ,imu,nucleation,nucleation_label,n_nucleation,tau,itau_alloc,tau_lucy,itauL_alloc,&
                 luminosity,eta_nimhd,eta_nimhd_label,apr_level
 use part,  only:metrics,metricderivs,tmunus
 use options,    only:use_dustfrac,use_porosity,use_var_comp,icooling
 use dump_utils, only:tag,open_dumpfile_w,allocate_header,&
                 free_header,write_header,write_array,write_block_header
 use mpiutils,   only:reduce_mpi,reduceall_mpi,start_threadwrite,end_threadwrite
 use timestep,   only:dtmax,idtmax_n,idtmax_frac
 use part,       only:ibin,krome_nmols,T_gas_cool
 use metric_tools, only:imetric, imet_et
 use eos_stamatellos, only:ttherm_store,ueqi_store,opac_store
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in), optional :: iorder(:)
 logical,          intent(in), optional :: sphNG
 integer(kind=8),  intent(in), optional :: ntotal

 integer, parameter :: isteps_sphNG = 0, iphase0 = 0
 integer(kind=8)    :: ilen(4)
 integer            :: nums(ndatatypes,4)
 integer            :: ipass,k,l,ioffset
 integer            :: ierr,nerr
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
 if (mpi) then
    nparttot = reduceall_mpi('+',npart)
    call update_npartoftypetot
 else
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
 endif
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
    write (idump,iostat=ierr) nblockarrays

 endif masterthread

 call start_threadwrite(id,idump,dumpfile)

 nerr = 0
 nums = 0
 ilen = 0_8
 if (sphNGdump) then
    write_itype = .true.
 else
    write_itype = any(npartoftypetot(2:) > 0)
 endif
 do ipass=1,2
    do k=1,ndatatypes
       nerr = 0
       !
       ! Block 1 arrays (hydrodynamics)
       !
       ilen(1) = int(npart,kind=8)
       if (write_itype) call write_array(1,iphase,'itype',npart,k,ipass,idump,nums,nerr,func=iamtype_int11)
       call write_array(1,xyzh,xyzh_label,3,npart,k,ipass,idump,nums,nerr)
       if (use_dustgrowth) then
          call write_array(1,dustprop,dustprop_label,2,npart,k,ipass,idump,nums,nerr)
          call write_array(1,VrelVf,VrelVf_label,npart,k,ipass,idump,nums,nerr)
          call write_array(1,dustgasprop,dustgasprop_label,4,npart,k,ipass,idump,nums,nerr)
          if (use_porosity) call write_array(1,filfac,filfac_label,npart,k,ipass,idump,nums,nerr)
       endif
       if (h2chemistry)  call write_array(1,abundance,abundance_label,nabundances,npart,k,ipass,idump,nums,nerr)
       if (use_dust) call write_array(1,dustfrac,dustfrac_label,ndusttypes,npart,k,ipass,idump,nums,nerr)
       if (use_dust) call write_array(1,tstop,tstop_label,ndustsmall,npart,k,ipass,idump,nums,nerr)
       if (use_dustfrac) then
          do l=1,ndustsmall
             call write_array(1,deltav(:,l,:),deltav_label,3,npart,k,ipass,idump,nums,nerr)
          enddo
       endif
       if (gr) then
          call write_array(1,pxyzu,pxyzu_label,maxvxyzu,npart,k,ipass,idump,nums,nerr)
          call write_array(1,dens,'dens prim',npart,k,ipass,idump,nums,nerr)
          if (imetric==imet_et) then
             ! Output metric if imetric=iet
             call write_array(1,metrics(1,1,1,:), 'gtt (covariant)',npart,k,ipass,idump,nums,nerr)
             !  call write_array(1,metrics(1,2,1,:), 'gtx (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             !  call write_array(1,metrics(1,3,1,:), 'gty (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             !  call write_array(1,metrics(1,2,1,:), 'gtz (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             !  call write_array(1,metrics(1,2,1,:), 'gtx (covariant)',npart,k,ipass,idump,nums,ierrs(8))
             call write_array(1,metrics(2,2,1,:), 'gxx (covariant)',npart,k,ipass,idump,nums,nerr)
             call write_array(1,metrics(3,3,1,:), 'gyy (covariant)',npart,k,ipass,idump,nums,nerr)
             call write_array(1,metrics(4,4,1,:), 'gzz (covariant)',npart,k,ipass,idump,nums,nerr)

             call write_array(1,metricderivs(1,1,1,:), 'dxgtt (covariant)',npart,k,ipass,idump,nums,nerr)
             call write_array(1,metricderivs(2,2,1,:), 'dxgxx (covariant)',npart,k,ipass,idump,nums,nerr)
             call write_array(1,metricderivs(3,3,1,:), 'dxgyy (covariant)',npart,k,ipass,idump,nums,nerr)
             call write_array(1,metricderivs(4,4,1,:), 'dxgzz (covariant)',npart,k,ipass,idump,nums,nerr)

             call write_array(1,tmunus(1,1,:),  'tmunutt (covariant)',npart,k,ipass,idump,nums,nerr)
          endif
       endif
       if (eos_is_non_ideal(ieos) .or. (.not.store_dust_temperature .and. icooling > 0)) then
          call write_array(1,eos_vars(itemp,:),eos_vars_label(itemp),npart,k,ipass,idump,nums,nerr)
       endif
       if (eos_is_non_ideal(ieos)) call write_array(1,eos_vars(igamma,:),eos_vars_label(igamma),npart,k,ipass,idump,nums,nerr)

       call write_array(1,vxyzu,vxyzu_label,maxvxyzu,npart,k,ipass,idump,nums,nerr)
       ! write pressure to file
       if ((eos_outputs_gasP(ieos) .or. eos_is_non_ideal(ieos)) .and. k==i_real) then
          call write_array(1,eos_vars,eos_vars_label,1,npart,k,ipass,idump,nums,nerr,index=igasP)
       endif
       ! write X, Z, mu to file
       if (eos_outputs_mu(ieos)) then
          call write_array(1,eos_vars,eos_vars_label,1,npart,k,ipass,idump,nums,nerr,index=imu)
          if (use_var_comp) then
             call write_array(1,eos_vars,eos_vars_label,1,npart,k,ipass,idump,nums,nerr,index=iX)
             call write_array(1,eos_vars,eos_vars_label,1,npart,k,ipass,idump,nums,nerr,index=iZ)
          endif
       endif
       ! write stamatellos cooling values
       if (icooling == 9) then
          call write_array(1,ueqi_store,'ueqi',npart,k,ipass,idump,nums,nerr)
          call write_array(1,ttherm_store,'ttherm',npart,k,ipass,idump,nums,nerr)
          call write_array(1,opac_store,'opacity',npart,k,ipass,idump,nums,nerr)
       endif
       ! smoothing length written as real*4 to save disk space
       call write_array(1,xyzh,xyzh_label,1,npart,k,ipass,idump,nums,nerr,use_kind=4,index=4)
       if (maxalpha==maxp) call write_array(1,alphaind,(/'alpha'/),1,npart,k,ipass,idump,nums,nerr)
       !if (maxalpha==maxp) then ! (uncomment this to write alphaloc to the full dumps)
       !   call write_array(1,alphaind,(/'alpha ','alphaloc'/),2,npart,k,ipass,idump,nums,ierrs(10))
       !endif
       if (ndivcurlv >= 1) call write_array(1,divcurlv,divcurlv_label,ndivcurlv,npart,k,ipass,idump,nums,nerr)
       !if (maxdvdx==maxp) call write_array(1,dvdx,dvdx_label,9,npart,k,ipass,idump,nums,ierrs(17))
       if (gravity .and. maxgrav==maxp) then
          call write_array(1,poten,'poten',npart,k,ipass,idump,nums,nerr)
       endif
       if (ind_timesteps) then
          if (.not.allocated(temparr)) allocate(temparr(npart))
          temparr(1:npart) = dtmax/2.**ibin(1:npart)
          call write_array(1,temparr,'dt',npart,k,ipass,idump,nums,nerr,use_kind=4)
       endif
       call write_array(1,iorig,'iorig',npart,k,ipass,idump,nums,nerr)

       if (lightcurve) then
          call write_array(1,luminosity,'luminosity',npart,k,ipass,idump,nums,nerr)
       endif

       if (use_apr) then
          call write_array(1,apr_level,'apr_level',npart,k,ipass,idump,nums,nerr)
       endif

       if (use_krome) then
          call write_array(1,abundance,abundance_label,krome_nmols,npart,k,ipass,idump,nums,nerr)
          call write_array(1,T_gas_cool,'temp',npart,k,ipass,idump,nums,nerr)
       endif
       if (update_muGamma .or. use_krome) then
          call write_array(1,eos_vars(imu,:),eos_vars_label(imu),npart,k,ipass,idump,nums,nerr)
          call write_array(1,eos_vars(igamma,:),eos_vars_label(igamma),npart,k,ipass,idump,nums,nerr)
       endif
       if (do_nucleation) then
          call write_array(1,nucleation,nucleation_label,n_nucleation,npart,k,ipass,idump,nums,nerr)
       endif
       If (itau_alloc == 1) then
          call write_array(1,tau,'tau',npart,k,ipass,idump,nums,nerr)
       endif
       If (itauL_alloc == 1) then
          call write_array(1,tau_lucy,'tau_lucy',npart,k,ipass,idump,nums,nerr)
       endif
       if (store_dust_temperature) then
          call write_array(1,dust_temp,'Tdust',npart,k,ipass,idump,nums,nerr)
       endif
       if (do_radiation) then
          call write_array(1,rad,rad_label,maxirad,npart,k,ipass,idump,nums,nerr)
          call write_array(1,radprop,radprop_label,maxradprop,npart,k,ipass,idump,nums,nerr)
       endif
       if (nerr > 0) call error('write_dump','error writing hydro arrays')
    enddo

    nerr = 0
    do k=1,ndatatypes
       !
       ! Block 2 arrays (sink particles)
       !
       if (.not. sphNGdump .and. nptmass > 0 .and. nptmass <= maxptmass) then
          ilen(2) = int(nptmass,kind=8)
          call write_array(2,xyzmh_ptmass,xyzmh_ptmass_label,nsinkproperties,nptmass,k,ipass,idump,nums,nerr)
          call write_array(2,vxyz_ptmass,vxyz_ptmass_label,3,nptmass,k,ipass,idump,nums,nerr)
          if (store_sf_ptmass) then
             call write_array(2,sf_ptmass,sf_ptmass_label,2,nptmass,k,ipass,idump,nums,nerr)
          endif
          if (nerr > 0) call error('write_dump','error writing sink particle arrays')
       endif
    enddo

    do k=1,ndatatypes
       !
       ! Block 4 arrays (MHD)
       !
       if (mhd) then
          nerr = 0
          ilen(4) = int(npart,kind=8)
          call write_array(4,Bxyz,Bxyz_label,3,npart,k,ipass,idump,nums,nerr) ! Bx,By,Bz
          call write_array(4,Bevol,Bevol_label,1,npart,k,ipass,idump,nums,nerr,index=4) ! psi
          if (ndivcurlB >= 1) then
             call write_array(4,divcurlB,divcurlB_label,ndivcurlB,npart,k,ipass,idump,nums,nerr)
          else
             call write_array(4,divBsymm,'divBsymm',npart,k,ipass,idump,nums,nerr)
          endif
          if (nerr > 0) call error('write_dump','error writing MHD arrays')
          if (mhd_nonideal) then
             nerr = 0
             call write_array(4,eta_nimhd,eta_nimhd_label,4,npart,k,ipass,idump,nums,nerr)
             if (nerr > 0) call error('write_dump','error writing non-ideal MHD arrays')
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
 use dim,        only:maxp,maxtypes,use_dust,lightcurve,use_dustgrowth,&
                      h2chemistry,use_apr
 use options,    only:use_porosity
 use io,         only:idump,iprint,real4,id,master,error,warning,nprocs
 use part,       only:xyzh,xyzh_label,npart,Bxyz,Bxyz_label,&
                      npartoftypetot,update_npartoftypetot,&
                      maxphase,iphase,nabundances,&
                      nptmass,nsinkproperties,xyzmh_ptmass,xyzmh_ptmass_label,&
                      abundance,abundance_label,mhd,dustfrac,iamtype_int11,&
                      dustprop,dustprop_label,dustfrac_label,&
                      filfac,filfac_label,ndusttypes,&
                      rad,rad_label,do_radiation,maxirad,luminosity,apr_level
 use dump_utils, only:open_dumpfile_w,dump_h,allocate_header,free_header,&
                      write_header,write_array,write_block_header
 use mpiutils,   only:reduceall_mpi,start_threadwrite,end_threadwrite
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
       call error('write_smalldump','could not write new dumpfile '//trim(dumpfile))
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
    write (idump,iostat=ierr) nblockarrays
    if (ierr /= 0) call error('write_smalldump','error writing nblockarrays')

 endif masterthread

 call start_threadwrite(id,idump,dumpfile)

 nums = 0
 ierr = 0
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
          if (use_porosity) call write_array(1,filfac,filfac_label,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       endif
       if (h2chemistry .and. nabundances >= 1) &
          call write_array(1,abundance,abundance_label,1,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       if (use_dust) &
          call write_array(1,dustfrac,dustfrac_label,ndusttypes,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       call write_array(1,xyzh,xyzh_label,4,npart,k,ipass,idump,nums,ierr,index=4,use_kind=4)

       if (lightcurve) call write_array(1,luminosity,'luminosity',npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       if (do_radiation) call write_array(1,rad,rad_label,maxirad,npart,k,ipass,idump,nums,ierr,singleprec=.true.)
       if (use_apr) then
          call write_array(1,apr_level,'apr_level',npart,k,ipass,idump,nums,ierr,func=iamtype_int11)
       endif
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
 use dump_utils,   only:get_dump_size,skipblock,skip_arrays,check_tag,lenid,ndatatypes,read_header, &
                        open_dumpfile_r,get_error_text,ierr_realsize,free_header,read_block_header,&
                        get_blocklimits
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
 character(len=12)     :: string
 type(dump_h)          :: hdr
 integer               :: i,ierrh

 if (id==master .and. iverbose >= 0) write(iprint,"(/,1x,a,i3)") '>>> reading setup from file: '//trim(dumpfile)//' on unit ',idisk1
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
 if (id==master .and. iverbose >= 0) write(iprint,*) trim(fileidentr)

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
 elseif (narraylengths < 1 .or. narraylengths > 4) then
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
    call get_blocklimits(nhydrothisblock,nblocks,nprocs,id,iblock,noffset,npartread,ierr)
    if (ierr /= 0) then
       call error('read_dump','could not map blocks in dump to number of threads')
       return
    endif
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
       string = ''
       if (nprocs > 1) write(string,'(a,i5)') 'thread',iblock
       if (iverbose >= 0) write(*,"(2(a,i10),a,i5,a,i10,'-',i10)") trim(string)//' reading particles ',noffset+1,&
           ':',noffset+npartread,', from block ',iblock,' lims=',i1,i2
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
    if (id==master .and. iverbose >= 0) write(iprint,"(a,/)") ' <<< finished reading (MHD) file '
 else
    if (id==master .and. iverbose >= 0) write(iprint,"(a,/)") ' <<< finished reading (hydro) file '
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
                        ierr_realsize,read_header,extract,free_header,read_block_header,get_blocklimits
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
 character(len=12)     :: string
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
 if (id==master .and. iverbose >= 0) write(iprint,*) trim(fileidentr)

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
    call get_blocklimits(nhydrothisblock,nblocks,nprocs,id,iblock,noffset,npartread,ierr)
    if (ierr /= 0) then
       call error('read_dump','could not map blocks in dump to number of threads')
       return
    endif
    i1 = i2 + 1
    i2 = i1 + (npartread - 1)
    npart = npart + npartread
    if (npart > maxp) then
       write(*,*) 'npart > maxp in readwrite_dumps'
       ierr = 1
       return
    endif
    if (npartread <= 0 .and. nptmass <= 0) then
       call skipblock(idisk1,nums(:,1),nums(:,2),nums(:,3),nums(:,4),tagged,ierr)
       if (ierr /= 0) then
          print*,' error skipping block'
          return
       endif
       cycle overblocks
    elseif (npartread > 0) then
       string = ''
       if (nprocs > 1) write(string,'(a,i5)') 'thread',iblock
       write(*,"(2(a,i10),a,i5,a,i10,'-',i10)") trim(string)//' reading particles ',noffset+1,&
           ':',noffset+npartread,', from block ',iblock,' lims=',i1,i2
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
                      ind_timesteps,use_krome,use_apr,store_sf_ptmass,mhd
 use part,       only:xyzh,xyzh_label,vxyzu,vxyzu_label,dustfrac,dustfrac_label,abundance,abundance_label, &
                      alphaind,poten,xyzmh_ptmass,xyzmh_ptmass_label,vxyz_ptmass,vxyz_ptmass_label,sf_ptmass, &
                      sf_ptmass_label,Bevol,Bxyz,Bxyz_label,nabundances,iphase,idust, &
                      eos_vars,eos_vars_label,maxeosvars,dustprop,dustprop_label,divcurlv,divcurlv_label,iX,iZ,imu, &
                      VrelVf,VrelVf_label,dustgasprop,dustgasprop_label,filfac,filfac_label,pxyzu,pxyzu_label,dust_temp, &
                      rad,rad_label,radprop,radprop_label,do_radiation,maxirad,maxradprop,ifluxx,ifluxy,ifluxz, &
                      nucleation,nucleation_label,n_nucleation,ikappa,tau,itau_alloc,tau_lucy,itauL_alloc,&
                      ithick,ilambda,iorig,dt_in,krome_nmols,T_gas_cool,apr_level
 use sphNGutils, only:mass_sphng,got_mass,set_gas_particle_mass
 use options,    only:use_porosity
 integer, intent(in)   :: i1,i2,noffset,narraylengths,nums(:,:),npartread,npartoftype(:),idisk1,iprint
 real,    intent(in)   :: massoftype(:)
 integer, intent(in)   :: nptmass,nsinkproperties
 logical, intent(in)   :: phantomdump,singleprec,tagged
 real,    intent(in)   :: tfile,alphafile
 integer, intent(out)  :: ierr
 logical               :: match
 logical               :: got_dustfrac(maxdusttypes)
 logical               :: got_iphase,got_xyzh(4),got_vxyzu(4),got_abund(nabundances),got_alpha(1),got_poten
 logical               :: got_sink_data(nsinkproperties),got_sink_vels(3),got_sink_sfprop(2),got_Bxyz(3)
 logical               :: got_krome_mols(krome_nmols),got_krome_T,got_krome_gamma,got_krome_mu
 logical               :: got_eosvars(maxeosvars),got_nucleation(n_nucleation),got_ray_tracer
 logical               :: got_psi,got_Tdust,got_dustprop(2),got_VrelVf,got_dustgasprop(4)
 logical               :: got_filfac,got_divcurlv(4),got_rad(maxirad),got_radprop(maxradprop),got_pxyzu(4),&
                          got_iorig,got_apr_level
 character(len=lentag) :: tag,tagarr(64)
 integer :: k,i,iarr,ik,ndustfraci
 real, allocatable :: tmparray(:)

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
 got_sink_sfprop  = .false.
 got_Bxyz        = .false.
 got_psi         = .false.
 got_eosvars     = .false.
 got_dustprop    = .false.
 got_VrelVf      = .false.
 got_filfac      = .false.
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
 got_apr_level   = .false.

 ndustfraci = 0
 if (use_dust .or. mhd) allocate(tmparray(max(size(dustfrac,2),size(Bevol,2))))
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
                if (use_porosity) call read_array(filfac,filfac_label,got_filfac,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (use_dust) then
                if (any(tag == dustfrac_label)) then
                   ndustfraci = ndustfraci + 1
                   call read_array(tmparray,dustfrac_label(ndustfraci),got_dustfrac(ndustfraci), &
                                   ik,i1,i2,noffset,idisk1,tag,match,ierr)
                   dustfrac(ndustfraci,i1:i2) = tmparray(i1:i2)
                endif
             endif
             if (h2chemistry) then
                call read_array(abundance,abundance_label,got_abund,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (use_krome) then
                call read_array(abundance,abundance_label,got_krome_mols,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                call read_array(T_gas_cool,'temp',got_krome_T,ik,i1,i2,noffset,idisk1,tag,match,ierr)
             endif
             if (use_apr) then
                call read_array(apr_level,'apr_level',got_apr_level,ik,i1,i2,noffset,idisk1,tag,match,ierr)
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
             if (store_sf_ptmass) then
                call read_array(sf_ptmass,sf_ptmass_label,got_sink_sfprop,ik,1,nptmass,0,idisk1,tag,match,ierr)
             endif
          end select
          select case(iarr)   ! MHD arrays can either be in block 1 or block 4
          case(1,4)
             if (mhd) then
                call read_array(Bxyz,Bxyz_label,got_Bxyz,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                call read_array(tmparray,'psi',got_psi,ik,i1,i2,noffset,idisk1,tag,match,ierr)
                Bevol(4,i1:i2) = tmparray(i1:i2)
             endif
          end select
          if (.not.match) then
             !write(*,*) 'skipping '//trim(tag)
             read(idisk1,end=100) ! skip unknown array
          endif
       enddo
    enddo

 enddo over_arraylengths
 if (allocated(tmparray)) deallocate(tmparray)
 !
 ! check for errors
 !
 call check_arrays(i1,i2,noffset,npartoftype,npartread,nptmass,nsinkproperties,massoftype,&
                   alphafile,tfile,phantomdump,got_iphase,got_xyzh,got_vxyzu,got_alpha, &
                   got_krome_mols,got_krome_gamma,got_krome_mu,got_krome_T, &
                   got_abund,got_dustfrac,got_sink_data,got_sink_vels,got_sink_sfprop,got_Bxyz, &
                   got_psi,got_dustprop,got_pxyzu,got_VrelVf,got_dustgasprop,got_rad, &
                   got_radprop,got_Tdust,got_eosvars,got_nucleation,got_iorig,  &
                   got_apr_level,iphase,xyzh,vxyzu,pxyzu,alphaind,xyzmh_ptmass,Bevol,iorig,iprint,ierr)
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
