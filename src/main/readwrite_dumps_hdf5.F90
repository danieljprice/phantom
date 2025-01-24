!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module readwrite_dumps_hdf5
!
! This module contains all routines related to the HDF5 data format.
!
! :References: None
!
! :Owner: Daniel Mentiplay
!
! :Runtime parameters: None
!
! :Dependencies: boundary, checkconserved, dim, eos, extern_binary,
!   extern_gwinspiral, externalforces, io, lumin_nsdisc, memory, mpiutils,
!   options, part, readwrite_dumps_common, setup_params, timestep, units,
!   utils_dumpfiles_hdf5
!
 use readwrite_dumps_common, only:check_arrays,fileident,get_options_from_fileid
 use utils_dumpfiles_hdf5,   only:create_hdf5file,         &
                                  open_hdf5file,           &
                                  close_hdf5file,          &
                                  hdf5_file_id,            &
                                  write_hdf5_header,       &
                                  write_hdf5_arrays,       &
                                  write_hdf5_arrays_small, &
                                  read_hdf5_header,        &
                                  read_hdf5_arrays,        &
                                  header_hdf5,             &
                                  got_arrays_hdf5,         &
                                  arrays_options_hdf5,     &
                                  externalforce_hdf5

 implicit none

 public :: read_dump_hdf5,read_smalldump_hdf5
 public :: write_smalldump_hdf5,write_fulldump_hdf5,write_dump_hdf5

 ! for use in analysis files if user wishes to skip small dumps
 logical, target, public :: opened_full_dump_hdf5
 ! to determine if dt has been read in so that ibin & ibinold can be set on restarts
 logical, target, public :: dt_read_in_hdf5

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
subroutine write_fulldump_hdf5(t,dumpfile,ntotal,iorder,sphNG)
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer(kind=8),  intent(in), optional :: ntotal
 integer,          intent(in), optional :: iorder(:)
 logical,          intent(in), optional :: sphNG
 if (present(ntotal)) then
    call write_dump_hdf5(t,dumpfile,fulldump=.true.,ntotal=ntotal)
 else
    call write_dump_hdf5(t,dumpfile,fulldump=.true.)
 endif
end subroutine write_fulldump_hdf5

!--------------------------------------------------------------------
!+
!  subroutine to write output to small dump file
!  (ie. minimal output...)
!
!  note that small dumps are always SINGLE PRECISION
!+
!-------------------------------------------------------------------
subroutine write_smalldump_hdf5(t,dumpfile)
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 call write_dump_hdf5(t,dumpfile,fulldump=.false.)
end subroutine write_smalldump_hdf5

!--------------------------------------------------------------------
!+
!  Generic subroutine for writing a dumpfile
!+
!-------------------------------------------------------------------
subroutine write_dump_hdf5(t,dumpfile,fulldump,ntotal,dtind)
 use dim,            only:maxp,maxvxyzu,gravity,maxalpha,mhd,mhd_nonideal,      &
                          use_dust,use_dustgrowth,phantom_version_major,        &
                          phantom_version_minor,phantom_version_micro,          &
                          phantom_version_string,use_krome,   &
                          store_dust_temperature,do_radiation,gr,do_nucleation, &
                          mpi,idumpfile
 use eos,            only:ieos,polyk,gamma,polyk2,qfacdisc,isink
 use io,             only:fatal,id,master,iprint
 use options,        only:tolh,alpha,alphau,alphaB,iexternalforce,use_dustfrac
 use part,           only:xyzh,vxyzu,Bevol,Bxyz,npart,npartoftype,maxtypes,    &
                          alphaind,rhoh,divBsymm,iphase,nptmass,               &
                          xyzmh_ptmass,vxyz_ptmass,get_pmass,abundance,        &
                          divcurlv,divcurlB,poten,dustfrac,deltav,tstop,       &
                          dustprop,VrelVf,dustgasprop,ndustsmall,              &
                          luminosity,eta_nimhd,massoftype,hfact,Bextx,Bexty,   &
                          Bextz,ndustlarge,idust,idustbound,grainsize,         &
                          graindens,h2chemistry,lightcurve,ndivcurlB,          &
                          ndivcurlv,pxyzu,dens,T_gas_cool,                     &
                          dust_temp,rad,radprop,itemp,igasP,eos_vars,iorig,    &
                          npartoftypetot,update_npartoftypetot
 use part,           only:nucleation
#ifdef IND_TIMESTEPS
 use part,           only:ibin
#endif
 use mpiutils,       only:reduce_mpi,reduceall_mpi
 use checkconserved, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use setup_params,   only:rhozero
 use timestep,       only:dtmax,C_cour,C_force,dtmax_user,idtmax_n_next,idtmax_frac_next
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax
 use units,          only:udist,umass,utime,unit_Bfield
 use externalforces, only:iext_gwinspiral,iext_binary,iext_corot_binary
 use extern_binary,  only:binary_posvel,a0,direction,accretedmass1,accretedmass2
 use extern_gwinspiral, only:Nstar
 use lumin_nsdisc,   only:beta
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 logical,          intent(in) :: fulldump
 integer(kind=8),  intent(in), optional :: ntotal
 real(kind=4),     intent(in), optional :: dtind(:)

 integer            :: i
 integer            :: ierr
 integer(kind=8)    :: nparttot
 logical            :: ind_timesteps,const_av,prdrag,isothermal
 real, allocatable  :: dtin(:),beta_pr(:)
 character(len=200) :: fileid,fstr,sstr
 real :: posmh(10)
 real :: vels(6)

 type (header_hdf5) :: hdr
 type (arrays_options_hdf5) :: array_options
 type (externalforce_hdf5) :: extern

 fstr = "(/,/,'-------->   TIME = ',g12.4,"//"': full dump written to file ',a,'   <--------',/)"
 sstr = "(/,/,'-------->   TIME = ',g12.4,"//"': small dump written to file ',a,'   <--------',/)"
 if (id==master) then
    if (fulldump) then
       write(iprint,fstr) t, trim(dumpfile)
    else
       write(iprint,sstr) t, trim(dumpfile)
    endif
 endif

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
       npartoftypetot = npartoftype
       if (all(npartoftypetot==0)) then
          npartoftypetot(1) = ntotal
       endif
    else
       nparttot = npart
       npartoftypetot = npartoftype
    endif
 endif

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

 if (fulldump) then
    allocate(beta_pr(npart),dtin(npart))

    ! Compute beta_pr array
    if (prdrag) then
       do i=1,int(npart)
          beta_pr(i) = beta(xyzh(1,i), xyzh(2,i), xyzh(3,i))
       enddo
    endif

    ! Compute dtind array
#ifdef IND_TIMESTEPS
    if (present(dtind)) then
       dtin = dtind
    else
       dtin = dtmax/2**ibin(1:npart)
    endif
#endif

 endif

! Check if constant AV
 if (maxp==maxalpha) then
    const_av = .false.
 else
    const_av = .true.
 endif

 ! generate fileid
 if (fulldump) then
    fileid = fileident('FT')
 else
    fileid = fileident('ST')
 endif

 ! create the HDF file
 call create_hdf5file(trim(dumpfile)//'.h5',hdf5_file_id,ierr)
 if (ierr/=0) call fatal('write_fulldump_hdf5','could not open file')

 ! construct header derived type
 hdr%fileident = trim(fileid)
 hdr%ntypes = maxtypes
 hdr%isink = isink
 hdr%nptmass = nptmass
 hdr%ndustlarge = ndustlarge
 hdr%ndustsmall = ndustsmall
 hdr%idust = idust
 hdr%idustbound = idustbound
 hdr%phantom_version_major = phantom_version_major
 hdr%phantom_version_minor = phantom_version_minor
 hdr%phantom_version_micro = phantom_version_micro
 hdr%nparttot = int(nparttot)
 hdr%npartoftypetot = reshape(int(npartoftypetot),shape(hdr%npartoftypetot),pad=[0])
 hdr%iexternalforce = iexternalforce
 hdr%ieos = ieos
 hdr%time = t
 hdr%dtmax = dtmax
 hdr%dtmax_user = dtmax_user
 hdr%idumpfile = idumpfile
 hdr%idtmax_n_next = idtmax_n_next
 hdr%idtmax_frac_next = idtmax_frac_next
 hdr%gamma = gamma
 hdr%rhozero = rhozero
 hdr%polyk = polyk
 hdr%hfact = hfact
 hdr%tolh = tolh
 hdr%C_cour = C_cour
 hdr%C_force = C_force
 hdr%alpha = alpha
 hdr%alphau = alphau
 hdr%alphaB = alphaB
 hdr%polyk2 = polyk2
 hdr%qfacdisc = qfacdisc
 hdr%massoftype = reshape(massoftype,shape(hdr%massoftype),pad=[0.0])
 hdr%Bextx = Bextx
 hdr%Bexty = Bexty
 hdr%Bextz = Bextz
 hdr%xmin = xmin
 hdr%xmax = xmax
 hdr%ymin = ymin
 hdr%ymax = ymax
 hdr%zmin = zmin
 hdr%zmax = zmax
 hdr%get_conserv = get_conserv
 hdr%etot_in = etot_in
 hdr%angtot_in = angtot_in
 hdr%totmom_in = totmom_in
 hdr%mdust_in = reshape(mdust_in,shape(hdr%mdust_in),pad=[0.0])
 if (.not. use_dustgrowth) then
    hdr%grainsize = reshape(grainsize,shape(hdr%grainsize),pad=[0.0])
    hdr%graindens = reshape(graindens,shape(hdr%graindens),pad=[0.0])
 endif
 hdr%udist = udist
 hdr%umass = umass
 hdr%utime = utime
 hdr%unit_Bfield = unit_Bfield

 ! contstruct external force derived type
 select case(hdr%iexternalforce)
 case(iext_gwinspiral)
    extern%Nstar = Nstar
 case(iext_binary,iext_corot_binary)
    call binary_posvel(t,posmh,vels)
    extern%xyzmh1 = posmh(1:5)
    extern%xyzmh2 = posmh(6:10)
    extern%vxyz1 = vels(1:3)
    extern%vxyz2 = vels(4:6)
    extern%accretedmass1 = accretedmass1
    extern%accretedmass2 = accretedmass2
    extern%a0 = a0
    extern%direction = direction
 end select

 ! write the header to the HDF file
 call write_hdf5_header(hdf5_file_id,hdr,extern,ierr)
 if (ierr/=0) call fatal('write_fulldump_hdf5','could not write header')

 ! create options derived type for writing arrays
 array_options%ieos = ieos
 array_options%isothermal = isothermal
 array_options%const_av = const_av
 array_options%ind_timesteps = ind_timesteps
 array_options%gravity = gravity
 array_options%mhd = mhd
 array_options%mhd_nonideal = mhd_nonideal
 array_options%use_dust = use_dust
 array_options%use_dustfrac = use_dustfrac
 array_options%use_dustgrowth = use_dustgrowth
 array_options%h2chemistry = h2chemistry
 array_options%lightcurve = lightcurve
 array_options%prdrag = prdrag
 array_options%ndivcurlB = ndivcurlB
 array_options%ndivcurlv = ndivcurlv
 array_options%ndustsmall = ndustsmall
 array_options%ndustlarge = ndustlarge
 array_options%store_dust_temperature = store_dust_temperature
 array_options%radiation = do_radiation
 array_options%krome = use_krome
 array_options%gr = gr
 array_options%nucleation = do_nucleation

 ! write the arrays to file
 if (fulldump) then
    call write_hdf5_arrays(hdf5_file_id, & ! File ID
                           ierr,         & ! Error code
                           npart,        & ! Num particles
                           nptmass,      & ! Num sinks
                           xyzh,         & !---------
                           vxyzu,        & ! Arrays
                           iphase,       & !
                           eos_vars,     & !
                           alphaind,     & !
                           dtin,         & !
                           poten,        & !
                           xyzmh_ptmass, & !
                           vxyz_ptmass,  & !
                           Bxyz,         & !
                           Bevol,        & !
                           divcurlB,     & !
                           divBsymm,     & !
                           eta_nimhd,    & !
                           iorig,        & !
                           dustfrac,     & !
                           tstop,        & !
                           deltav,       & !
                           dustprop,     & !
                           VrelVf,       & !
                           dustgasprop,  & !
                           abundance,    & !
                           divcurlv,     & !
                           luminosity,   & !
                           beta_pr,      & !
                           pxyzu,        & !
                           dens,         & !
                           T_gas_cool,   & !
                           nucleation,   & !
                           dust_temp,    & !
                           rad,          & !
                           radprop,      & !---------
                           array_options)  ! Options
 else
    call write_hdf5_arrays_small(hdf5_file_id, & ! File ID
                                 ierr,         & ! Error code
                                 npart,        & ! Num particles
                                 nptmass,      & ! Num sinks
                                 xyzh,         & !--------
                                 iphase,       & !
                                 xyzmh_ptmass, & ! Arrays
                                 Bxyz,         & !
                                 dustfrac,     & !
                                 dustprop,     & !
                                 abundance,    & !
                                 luminosity,   & !
                                 rad,          & !--------
                                 array_options)  ! Options
 endif
 if (ierr/=0) call fatal('write_fulldump_hdf5','could not write arrays')

 call close_hdf5file(hdf5_file_id,ierr)
 if (ierr/=0) call fatal('write_fulldump_hdf5','could not close file')

 if (fulldump) deallocate(beta_pr,dtin)

end subroutine write_dump_hdf5

!--------------------------------------------------------------------
!+
!  subroutine to read a full dump from file
!+
!-------------------------------------------------------------------
subroutine read_dump_hdf5(                                                    &
   dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc &
)
 character(len=*),  intent(in)  :: dumpfile
 real,              intent(out) :: tfile,hfactfile
 integer,           intent(in)  :: idisk1,iprint,id,nprocs
 integer,           intent(out) :: ierr
 logical, optional, intent(in)  :: headeronly
 logical, optional, intent(in)  :: dustydisc

 call read_any_dump_hdf5( &
    dumpfile,             &
    tfile,                &
    hfactfile,            &
    idisk1,               &
    iprint,               &
    id,                   &
    nprocs,               &
    ierr,                 &
    headeronly,           &
    dustydisc,            &
    acceptsmall=.false.   &
 )

end subroutine read_dump_hdf5

!--------------------------------------------------------------------
!+
!  subroutine to read a small dump from file
!+
!-------------------------------------------------------------------
subroutine read_smalldump_hdf5(                                               &
   dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc &
)
 character(len=*),  intent(in)  :: dumpfile
 real,              intent(out) :: tfile,hfactfile
 integer,           intent(in)  :: idisk1,iprint,id,nprocs
 integer,           intent(out) :: ierr
 logical, optional, intent(in)  :: headeronly
 logical, optional, intent(in)  :: dustydisc

 call read_any_dump_hdf5( &
    dumpfile,             &
    tfile,                &
    hfactfile,            &
    idisk1,               &
    iprint,               &
    id,                   &
    nprocs,               &
    ierr,                 &
    headeronly,           &
    dustydisc,            &
    acceptsmall=.true.    &
 )

end subroutine read_smalldump_hdf5

!--------------------------------------------------------------------
!+
!  subroutine to read full/small dump from file
!+
!-------------------------------------------------------------------
subroutine read_any_dump_hdf5(                                                            &
   dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc,acceptsmall &
)
 use boundary,       only:set_boundary
 use dim,            only:maxp,gravity,maxalpha,mhd,use_dust,use_dustgrowth, &
                          h2chemistry,nsinkproperties,     &
                          maxp_hard,use_krome,store_dust_temperature,        &
                          do_radiation,do_nucleation,gr,idumpfile
 use eos,            only:ieos,polyk,gamma,polyk2,qfacdisc,isink
 use checkconserved, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use io,             only:fatal,error
 use memory,         only:allocate_memory
 use options,        only:iexternalforce,use_dustfrac
 use part,           only:iphase,xyzh,vxyzu,npart,npartoftype,massoftype,      &
                          nptmass,xyzmh_ptmass,vxyz_ptmass,ndustlarge,         &
                          ndustsmall,grainsize,graindens,Bextx,Bexty,Bextz,    &
                          alphaind,poten,Bxyz,Bevol,dustfrac,deltav,dustprop,  &
                          dustgasprop,VrelVf,eos_vars,abundance,               &
                          periodic,ndusttypes,pxyzu,T_gas_cool,dust_temp,      &
                          nucleation,rad,radprop,igasP,itemp,iorig
#ifdef IND_TIMESTEPS
 use part,           only:dt_in
#endif
 use setup_params,   only:rhozero
 use timestep,       only:dtmax_user,idtmax_n,idtmax_frac
 use units,          only:udist,umass,utime,unit_Bfield,set_units_extra
 use externalforces, only:iext_gwinspiral,iext_binary,iext_corot_binary
 use extern_gwinspiral, only:Nstar
 use extern_binary,  only:a0,direction,accretedmass1,accretedmass2
 character(len=*),  intent(in)  :: dumpfile
 real,              intent(out) :: tfile,hfactfile
 integer,           intent(in)  :: idisk1,iprint,id,nprocs
 integer,           intent(out) :: ierr
 integer(kind=8)                :: nparttot
 logical, optional, intent(in)  :: headeronly
 logical, optional, intent(in)  :: dustydisc
 logical, optional, intent(in)  :: acceptsmall

 real(kind=4), allocatable :: dtind(:)

 real :: xmin,xmax,ymin,ymax,zmin,zmax
 real :: dtmaxi,tolhfile,C_courfile,C_forcefile,alphafile,alphaufile,alphaBfile
 character(len=200) :: fileid
 logical :: smalldump,isothermal,ind_timesteps,const_av,small_ok,tagged,phantomdump

 type (header_hdf5) :: hdr
 type (arrays_options_hdf5) :: array_options
 type (got_arrays_hdf5) :: got_arrays
 type (externalforce_hdf5) :: extern

 small_ok = .false.
 if (present(acceptsmall)) then
    if (acceptsmall) small_ok = .true.
 endif

 call open_hdf5file(trim(dumpfile)//'.h5',hdf5_file_id,ierr)
 if (ierr /= 0) then
    ierr = 1
    call fatal('read_dump_hdf5',trim(dumpfile)//'.h5 does not exist')
 endif

 call read_hdf5_header(hdf5_file_id,hdr,extern,ierr)
 if (ierr /= 0) then
    ierr = 2
    call error('read_dump_hdf5','cannot read header')
 endif

 !
 !--Set values from header
 !
 fileid = hdr%fileident
 isink = hdr%isink
 nptmass = hdr%nptmass
 ndustlarge = hdr%ndustlarge
 ndustsmall = hdr%ndustsmall
 npart = hdr%nparttot
 npartoftype = reshape(hdr%npartoftypetot,shape(npartoftype),pad=[0])
 iexternalforce = hdr%iexternalforce
 ieos = hdr%ieos
 tfile = hdr%time
 dtmaxi = hdr%dtmax
 dtmax_user = hdr%dtmax_user
 idumpfile = hdr%idumpfile
 idtmax_n = hdr%idtmax_n_next
 idtmax_frac = hdr%idtmax_frac_next
 gamma = hdr%gamma
 rhozero = hdr%rhozero
 polyk = hdr%polyk
 hfactfile = hdr%hfact
 tolhfile = hdr%tolh
 C_courfile = hdr%C_cour
 C_forcefile = hdr%C_force
 alphafile = hdr%alpha
 alphaufile = hdr%alphau
 alphaBfile = hdr%alphaB
 polyk2 = hdr%polyk2
 qfacdisc = hdr%qfacdisc
 massoftype = reshape(hdr%massoftype,shape(massoftype),pad=[0.0])
 Bextx = hdr%Bextx
 Bexty = hdr%Bexty
 Bextz = hdr%Bextz
 xmin = hdr%xmin
 xmax = hdr%xmax
 ymin = hdr%ymin
 ymax = hdr%ymax
 zmin = hdr%zmin
 zmax = hdr%zmax
 get_conserv = hdr%get_conserv
 etot_in = hdr%etot_in
 angtot_in = hdr%angtot_in
 totmom_in = hdr%totmom_in
 mdust_in = reshape(hdr%mdust_in,shape(mdust_in),pad=[0.0])
 if (.not. use_dustgrowth) then
    grainsize = reshape(hdr%grainsize,shape(grainsize),pad=[0.0])
    graindens = reshape(hdr%graindens,shape(graindens),pad=[0.0])
 endif
 udist = hdr%udist
 umass = hdr%umass
 utime = hdr%utime
 unit_Bfield = hdr%unit_Bfield

 call set_units_extra()
 ndusttypes = ndustsmall + ndustlarge

 call get_options_from_fileid(fileid,tagged,phantomdump,smalldump,use_dustfrac,ierr)

 !
 !--Set values from external forces
 !
 select case(hdr%iexternalforce)
 case(iext_gwinspiral)
    Nstar = extern%Nstar
 case(iext_binary,iext_corot_binary)
    a0 = extern%a0
    direction = extern%direction
    accretedmass1 = extern%accretedmass1
    accretedmass2 = extern%accretedmass2
 end select

 !
 !--Allocate main arrays
 !
#ifdef INJECT_PARTICLES
 call allocate_memory(maxp_hard)
#else
 call allocate_memory(int(min(nprocs,4)*nparttot/nprocs,8))
#endif

 if (periodic) then
    call set_boundary(xmin,xmax,ymin,ymax,zmin,zmax)
 endif

#ifdef ISOTHERMAL
 isothermal = .true.
#else
 isothermal = .false.
#endif

#ifdef IND_TIMESTEPS
 ind_timesteps = .true.
#else
 ind_timesteps = .false.
#endif

! Check if constant AV
 if (maxp==maxalpha) then
    const_av = .false.
 else
    const_av = .true.
 endif

 array_options%isothermal = isothermal
 array_options%const_av = const_av
 array_options%ind_timesteps = ind_timesteps
 array_options%gravity = gravity
 array_options%mhd = mhd
 array_options%use_dust = use_dust
 array_options%use_dustfrac = use_dustfrac
 array_options%use_dustgrowth = use_dustgrowth
 array_options%h2chemistry = h2chemistry
 array_options%store_dust_temperature = store_dust_temperature
 array_options%radiation = do_radiation
 array_options%krome = use_krome
 array_options%gr = gr
 array_options%nucleation = do_nucleation

 allocate(dtind(npart))
 call read_hdf5_arrays(hdf5_file_id,  &
                       ierr,          &
                       npart,         &
                       nptmass,       &
                       ndusttypes,    &
                       iphase,        &
                       xyzh,          &
                       vxyzu,         &
                       xyzmh_ptmass,  &
                       vxyz_ptmass,   &
                       eos_vars,      &
                       dtind,         &
                       alphaind,      &
                       poten,         &
                       Bxyz,          &
                       Bevol,         &
                       iorig,         &
                       dustfrac,      &
                       deltav,        &
                       dustprop,      &
                       VrelVf,        &
                       dustgasprop,   &
                       abundance,     &
                       pxyzu,         &
                       T_gas_cool,    &
                       nucleation,    &
                       dust_temp,     &
                       rad,           &
                       radprop,       &
                       array_options, &
                       got_arrays)
 if (ierr /= 0) then
    ierr = 3
    call error('read_dump_hdf5','cannot read arrays')
 endif

 if (.not.smalldump) then
    call check_arrays(1,                          &
                      npart,                      &
                      0,                          &
                      npartoftype,                &
                      npart,                      &
                      nptmass,                    &
                      nsinkproperties,            &
                      massoftype,                 &
                      alphafile,                  &
                      tfile,                      &
                      .true.,                     &
                      got_arrays%got_iphase,      &
                      got_arrays%got_xyzh,        &
                      got_arrays%got_vxyzu,       &
                      got_arrays%got_alpha,       &
                      got_arrays%got_krome_mols,  &
                      got_arrays%got_krome_gamma, &
                      got_arrays%got_krome_mu,    &
                      got_arrays%got_krome_T,     &
                      got_arrays%got_x,           &
                      got_arrays%got_z,           &
                      got_arrays%got_mu,          &
                      got_arrays%got_abund,       &
                      got_arrays%got_dustfrac,    &
                      got_arrays%got_sink_data,   &
                      got_arrays%got_sink_vels,   &
                      got_arrays%got_Bxyz,        &
                      got_arrays%got_psi,         &
                      got_arrays%got_dustprop,    &
                      got_arrays%got_pxyzu,       &
                      got_arrays%got_VrelVf,      &
                      got_arrays%got_dustgasprop, &
                      got_arrays%got_temp,        &
                      got_arrays%got_raden,       &
                      got_arrays%got_kappa,       &
                      got_arrays%got_Tdust,       &
                      got_arrays%got_nucleation,  &
                      got_arrays%got_iorig,       &
                      iphase,                     &
                      xyzh,                       &
                      vxyzu,                      &
                      pxyzu,                      &
                      alphaind,                   &
                      xyzmh_ptmass,               &
                      Bevol,                      &
                      iorig,                      &
                      iprint,                     &
                      ierr)
 endif
 if (ierr /= 0) then
    ierr = 4
    call error('read_dump_hdf5','error in checking arrays')
 endif

#ifdef IND_TIMESTEPS
 if (size(dt_in)/=size(dtind)) then
    call error('read_dump_hdf5','problem reading individual timesteps')
 endif
 dt_in = dtind
#endif
 deallocate(dtind)

 call close_hdf5file(hdf5_file_id,ierr)
 if (ierr /= 0) then
    ierr = 5
    call error('read_dump_hdf5','cannot close file')
 endif

 if (smalldump) then
    if (.not.small_ok) then
       call error('read_dump_hdf5',trim(dumpfile)//'.h5 is not a full dump')
       ierr = is_small_dump
    endif
 endif

end subroutine read_any_dump_hdf5

end module readwrite_dumps_hdf5
