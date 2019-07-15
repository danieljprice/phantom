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
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, eos, extern_binary, extern_gwinspiral,
!    externalforces, gitinfo, initial_params, io, lumin_nsdisc, memory,
!    mpiutils, options, part, setup_params, timestep, units,
!    utils_dumpfiles_hdf5
!+
!--------------------------------------------------------------------------
#ifdef PHANTOM2HDF5
module readwrite_dumps_hdf5
#else
module readwrite_dumps
#endif
 use utils_dumpfiles_hdf5, only:create_hdf5file,         &
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
 character(len=80), parameter, public :: &    ! module version
    modid="$Id$"

 public :: read_dump,read_smalldump,write_smalldump,write_fulldump,write_gadgetdump
 logical, public    :: opened_full_dump       ! for use in analysis files if user wishes to skip small dumps
 logical, public    :: dt_read_in             ! to determine if dt has been read in so that ibin & ibinold can be set on restarts
 integer, parameter, public :: is_small_dump = 1978
 integer, parameter, public :: is_not_mhd = 1979

 private

contains

!--------------------------------------------------------------------
!+
!  extract various options used in Phantom from the fileid string
!+
!--------------------------------------------------------------------
subroutine get_options_from_fileid(fileid,smalldump,use_onefluiddust,ierr)
 character(len=*), intent(in)  :: fileid
 logical,          intent(out) :: smalldump,use_onefluiddust
 integer,          intent(out) :: ierr
!
!--if file is a small dump, return an error code but still read what
!  can be read from a small dump
!
 ierr = 0
 smalldump = .false.
 if (fileid(1:4) /= 'full') then
    if (fileid(1:5)=='small') then
       smalldump = .true.
    else
       ierr = 1
    endif
 endif
 if (index(fileid,'+1dust') /= 0) then
    use_onefluiddust = .true.
 else
    use_onefluiddust = .false.
 endif

end subroutine get_options_from_fileid

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

!--------------------------------------------------------------------
!+
!  Generic subroutine for writing a dumpfile
!+
!-------------------------------------------------------------------
subroutine write_dump(t,dumpfile,fulldump,ntotal)
 use dim,            only:maxp,maxvxyzu,gravity,maxalpha,mhd,mhd_nonideal,   &
                          use_dust,use_dustgrowth,phantom_version_major,     &
                          phantom_version_minor,phantom_version_micro,       &
                          store_temperature,phantom_version_string,use_krome
 use eos,            only:ieos,equationofstate,done_init_eos,init_eos,polyk, &
                          gamma,polyk2,qfacdisc,isink
 use gitinfo,        only:gitsha
 use io,             only:nprocs,fatal,id,master,iprint
 use options,        only:tolh,alpha,alphau,alphaB,iexternalforce,use_dustfrac
 use part,           only:xyzh,vxyzu,Bevol,Bxyz,npart,npartoftype,maxtypes,  &
                          alphaind,rhoh,divBsymm,maxphase,iphase,nptmass,    &
                          xyzmh_ptmass,vxyz_ptmass,get_pmass,abundance,      &
                          divcurlv,divcurlB,poten,dustfrac,deltav,tstop,     &
                          dustprop,temperature,St,ndustsmall,luminosity,     &
                          eta_nimhd,massoftype,hfact,Bextx,Bexty,Bextz,      &
                          ndustlarge,idust,grainsize,graindens,              &
                          h2chemistry,lightcurve,maxBevol,                   &
                          ndivcurlB,ndivcurlv
#ifdef IND_TIMESTEPS
 use part,           only:ibin
#ifdef PHANTOM2HDF5
 use part,           only:dt_in
#endif
#endif
#ifdef KROME
 use part,           only:gamma_chem
#endif
 use mpiutils,       only:reduce_mpi,reduceall_mpi
 use lumin_nsdisc,   only:beta
 use initial_params, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use setup_params,   only:rhozero
 use timestep,       only:dtmax,C_cour,C_force
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax
 use units,          only:udist,umass,utime,unit_Bfield
 use externalforces, only:iext_gwinspiral,iext_binary,iext_corot_binary
 use extern_binary,  only:binary_posvel,a0,direction,accretedmass1,accretedmass2
 use extern_gwinspiral, only:Nstar
 real,             intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 logical,          intent(in) :: fulldump
 integer(kind=8),  intent(in), optional :: ntotal

 integer           :: i
 integer           :: ierr,nblocks
 integer(kind=8)   :: nparttot,npartoftypetot(maxtypes)
 logical           :: use_gas,ind_timesteps,const_av,prdrag,isothermal
 real              :: ponrhoi,rhoi,spsoundi
 real, allocatable :: pressure(:),dtind(:),beta_pr(:)
 character(len=100):: fileident
 character(len=10) :: datestring, timestring
 character(len=30) :: string
 character(len=9)  :: dumptype
 integer :: error
 real :: posmh(10)
 real :: vels(6)

 type (header_hdf5) :: hdr
 type (arrays_options_hdf5) :: array_options
 type (externalforce_hdf5) :: extern

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
    allocate(pressure(npart),beta_pr(npart),dtind(npart))

    ! Compute pressure and beta_pr array
    if (.not.done_init_eos) call init_eos(ieos,ierr)
    !$omp parallel do default(none) &
    !$omp shared(xyzh,vxyzu,ieos,npart,pressure,beta_pr,temperature,use_gas,prdrag) &
    !$omp private(i,ponrhoi,spsoundi,rhoi)
    do i=1,int(npart)
       rhoi = rhoh(xyzh(4,i),get_pmass(i,use_gas))
       if (maxvxyzu >=4 ) then
          if (use_krome) then
             call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xyzh(1,i),xyzh(2,i),xyzh(3,i),eni=vxyzu(4,i), &
                                  gamma_local=gamma_chem(i))
          else if (store_temperature) then
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
#ifdef IND_TIMESTEPS
#ifdef PHANTOM2HDF5
    if (ind_timesteps) dtind = dt_in
#else
    if (ind_timesteps) dtind = dtmax/2**ibin(1:npart)
#endif
#endif
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

 ! create the HDF file
 call create_hdf5file(trim(dumpfile)//'.h5',hdf5_file_id,error)
 if (error/=0) call fatal('write_fulldump_hdf5','could not open file')

 ! construct header derived type
 hdr%fileident = trim(fileident)
 hdr%ntypes = maxtypes
 hdr%isink = isink
 hdr%nptmass = nptmass
 hdr%ndustlarge = ndustlarge
 hdr%ndustsmall = ndustsmall
 hdr%idust = idust
 hdr%phantom_version_major = phantom_version_major
 hdr%phantom_version_minor = phantom_version_minor
 hdr%phantom_version_micro = phantom_version_micro
 hdr%nparttot = int(nparttot)
 hdr%npartoftypetot = reshape(int(npartoftypetot),shape(hdr%npartoftypetot),pad=[0])
 hdr%iexternalforce = iexternalforce
 hdr%ieos = ieos
 hdr%time = t
 hdr%dtmax = dtmax
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
 hdr%grainsize = reshape(grainsize,shape(hdr%grainsize),pad=[0.0])
 hdr%graindens = reshape(graindens,shape(hdr%graindens),pad=[0.0])
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
 call write_hdf5_header(hdf5_file_id,hdr,extern,error)
 if (error/=0) call fatal('write_fulldump_hdf5','could not write header')

 ! create options derived type for writing arrays
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
 array_options%store_temperature = store_temperature
 array_options%maxBevol = maxBevol
 array_options%ndivcurlB = ndivcurlB
 array_options%ndivcurlv = ndivcurlv
 array_options%ndustsmall = ndustsmall
 array_options%ndustlarge = ndustlarge

 ! write the arrays to file
 if (fulldump) then
    call write_hdf5_arrays(hdf5_file_id, & ! File ID
                           error,        & ! Error code
                           npart,        & ! # particles
                           nptmass,      & ! # sinks
                           xyzh,         & !---------
                           vxyzu,        & !
                           iphase,       & !
                           pressure,     & !
                           alphaind,     & !
                           dtind,        & !
                           poten,        & !
                           xyzmh_ptmass, & !
                           vxyz_ptmass,  & !
                           Bxyz,         & !
                           Bevol,        & ! Arrays
                           divcurlB,     & !
                           divBsymm,     & !
                           eta_nimhd,    & !
                           dustfrac,     & !
                           tstop,        & !
                           deltav,       & !
                           dustprop,     & !
                           st,           & !
                           abundance,    & !
                           temperature,  & !
                           divcurlv,     & !
                           luminosity,   & !
                           beta_pr,      & !---------
                           array_options)                         ! Options
 else
    call write_hdf5_arrays_small(hdf5_file_id, & ! File ID
                                 error,        & ! Error code
                                 npart,        & ! # particles
                                 nptmass,      & ! # sinks
                                 xyzh,         & !--------
                                 iphase,       & !
                                 xyzmh_ptmass, & !
                                 Bxyz,         & !
                                 dustfrac,     & ! Arrays
                                 dustprop,     & !
                                 st,           & !
                                 abundance,    & !
                                 luminosity,   & !--------
                                 array_options)  ! Options
 endif
 if (error/=0) call fatal('write_fulldump_hdf5','could not write arrays')

 call close_hdf5file(hdf5_file_id,error)
 if (error/=0) call fatal('write_fulldump_hdf5','could not close file')

 if (fulldump) deallocate(pressure,beta_pr,dtind)

end subroutine write_dump

!--------------------------------------------------------------------
!+
!  subroutine to read full dump from file
!+
!-------------------------------------------------------------------
subroutine read_dump(dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc)
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax
 use dim,            only:maxp,gravity,maxalpha,mhd,use_dust,use_dustgrowth, &
                          h2chemistry,store_temperature,nsinkproperties,maxp_hard
 use eos,            only:ieos,polyk,gamma,polyk2,qfacdisc,isink
 use initial_params, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use io,             only:fatal,error
 use memory,         only:allocate_memory
 use options,        only:tolh,alpha,alphau,alphaB,iexternalforce,use_dustfrac,use_moddump
 use part,           only:iphase,xyzh,vxyzu,npart,npartoftype,massoftype,     &
                          nptmass,xyzmh_ptmass,vxyz_ptmass,ndustlarge,        &
                          ndustsmall,grainsize,graindens,Bextx,Bexty,Bextz,   &
                          alphaind,poten,Bxyz,Bevol,dustfrac,deltav,dustprop, &
                          tstop,St,temperature,abundance,ndusttypes
#ifdef IND_TIMESTEPS
 use part,           only:dt_in
#endif
 use setup_params,   only:rhozero
 use timestep,       only:dtmax,C_cour,C_force
 use units,          only:udist,umass,utime,unit_Bfield,set_units_extra
 use externalforces, only:iext_gwinspiral,iext_binary,iext_corot_binary
 use extern_gwinspiral, only:Nstar
 use extern_binary,  only:a0,direction,accretedmass1,accretedmass2
 character(len=*),  intent(in)  :: dumpfile
 real,              intent(out) :: tfile,hfactfile
 integer,           intent(in)  :: idisk1,iprint,id,nprocs
 integer,           intent(out) :: ierr
 logical, optional, intent(in)  :: headeronly
 logical, optional, intent(in)  :: dustydisc

 real(kind=4), allocatable :: dtind(:)

 character(len=200) :: fileident
 integer :: errors(5)
 logical :: smalldump,isothermal,ind_timesteps,const_av

 type (header_hdf5) :: hdr
 type (arrays_options_hdf5) :: array_options
 type (got_arrays_hdf5) :: got_arrays
 type (externalforce_hdf5) :: extern

 errors(:) = 0

 call open_hdf5file(trim(dumpfile)//'.h5',hdf5_file_id,errors(1))
 if (errors(1) /= 0) then
    ierr = 1
    call fatal('read_dump',trim(dumpfile)//'.h5 does not exist')
 endif

 call read_hdf5_header(hdf5_file_id,hdr,extern,errors(2))

 !
 !--Set values from header
 !
 fileident = hdr%fileident
 isink = hdr%isink
 nptmass = hdr%nptmass
 ndustlarge = hdr%ndustlarge
 ndustsmall = hdr%ndustsmall
 npart = hdr%nparttot
 npartoftype = reshape(hdr%npartoftypetot,shape(npartoftype),pad=[0])
 iexternalforce = hdr%iexternalforce
 ieos = hdr%ieos
 tfile = hdr%time
 dtmax = hdr%dtmax
 gamma = hdr%gamma
 rhozero = hdr%rhozero
 polyk = hdr%polyk
 hfactfile = hdr%hfact
 tolh = hdr%tolh
 C_cour = hdr%C_cour
 C_force = hdr%C_force
 alpha = hdr%alpha
 alphau = hdr%alphau
 alphaB = hdr%alphaB
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
 grainsize = reshape(hdr%grainsize,shape(grainsize),pad=[0.0])
 graindens = reshape(hdr%graindens,shape(graindens),pad=[0.0])
 udist = hdr%udist
 umass = hdr%umass
 utime = hdr%utime
 unit_Bfield = hdr%unit_Bfield

 call set_units_extra()
 ndusttypes = ndustsmall + ndustlarge

 call get_options_from_fileid(fileident,smalldump,use_dustfrac,errors(3))

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
 if (.not. use_moddump) then
    call allocate_memory(int(npart / nprocs))
 else
    ! This is required for the cases when particles will be added during moddump
    call allocate_memory(maxp_hard)
 endif
#endif

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
 array_options%store_temperature = store_temperature

 if (.not.smalldump) then
    allocate(dtind(npart))
    call read_hdf5_arrays(hdf5_file_id,errors(4),npart,nptmass,iphase,xyzh,    &
                          vxyzu,xyzmh_ptmass,vxyz_ptmass,dtind,alphaind,poten, &
                          Bxyz,Bevol,dustfrac,deltav,dustprop,tstop,St,        &
                          temperature,abundance,array_options,got_arrays)

    call check_arrays(1,npart,npartoftype,nptmass,nsinkproperties,massoftype, &
                      alpha,tfile,got_arrays%got_iphase,got_arrays%got_xyzh,  &
                      got_arrays%got_vxyzu,got_arrays%got_alpha,              &
                      got_arrays%got_abund,got_arrays%got_dustfrac,           &
                      got_arrays%got_sink_data,got_arrays%got_sink_vels,      &
                      got_arrays%got_Bxyz,got_arrays%got_psi,                 &
                      got_arrays%got_dustprop,got_arrays%got_St,              &
                      got_arrays%got_temp,iphase,xyzh,vxyzu,alphaind,         &
                      xyzmh_ptmass,Bevol,iprint,ierr)

#ifdef IND_TIMESTEPS
    if (size(dt_in)/=size(dtind)) call error('read_smalldump','problem reading individual timesteps')
    dt_in = dtind
#endif
    deallocate(dtind)

 else
    call error('read_dump',trim(dumpfile)//'.h5 is not a full dump')
 endif

 call close_hdf5file(hdf5_file_id,errors(5))

 ierr = maxval(abs(errors))

 if (smalldump) ierr = is_small_dump

end subroutine read_dump

!--------------------------------------------------------------------
!+
!  subroutine to read full dump from file
!+
!-------------------------------------------------------------------
subroutine read_smalldump(dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,headeronly,dustydisc)
 use boundary,       only:xmin,xmax,ymin,ymax,zmin,zmax
 use dim,            only:maxp,gravity,maxalpha,mhd,use_dust,use_dustgrowth, &
                          h2chemistry,store_temperature,nsinkproperties
 use eos,            only:ieos,polyk,gamma,polyk2,qfacdisc,isink
 use initial_params, only:get_conserv,etot_in,angtot_in,totmom_in,mdust_in
 use io,             only:error,fatal
 use memory,         only:allocate_memory
 use options,        only:tolh,alpha,alphau,alphaB,iexternalforce,use_dustfrac
 use part,           only:iphase,xyzh,vxyzu,npart,npartoftype,massoftype,   &
                          nptmass,xyzmh_ptmass,vxyz_ptmass,ndustlarge,      &
                          ndustsmall,grainsize,graindens,Bextx,Bexty,Bextz, &
                          alphaind,poten,Bxyz,Bevol,dustfrac,deltav,        &
                          dustprop,tstop,St,temperature,abundance,ndusttypes
#ifdef IND_TIMESTEPS
 use part,           only:dt_in
#endif
#ifdef INJECT_PARTICLES
 use dim,            only:maxp_hard
#endif
 use setup_params,   only:rhozero
 use timestep,       only:dtmax,C_cour,C_force
 use units,          only:udist,umass,utime,unit_Bfield,set_units_extra
 use externalforces, only:iext_gwinspiral,iext_binary,iext_corot_binary
 use extern_gwinspiral, only:Nstar
 use extern_binary,  only:a0,direction,accretedmass1,accretedmass2
 character(len=*),  intent(in)  :: dumpfile
 real,              intent(out) :: tfile,hfactfile
 integer,           intent(in)  :: idisk1,iprint,id,nprocs
 integer,           intent(out) :: ierr
 logical, optional, intent(in)  :: headeronly
 logical, optional, intent(in)  :: dustydisc

 real(kind=4) :: dtind(npart)

 character(len=200) :: fileident
 integer :: errors(5)
 logical :: smalldump,isothermal,ind_timesteps,const_av

 type (header_hdf5) :: hdr
 type (arrays_options_hdf5) :: array_options
 type (got_arrays_hdf5) :: got_arrays
 type (externalforce_hdf5) :: extern

 call open_hdf5file(trim(dumpfile)//'.h5',hdf5_file_id,errors(1))

 call read_hdf5_header(hdf5_file_id,hdr,extern,errors(2))

 !
 !--Set values from header
 !
 fileident = hdr%fileident
 isink = hdr%isink
 nptmass = hdr%nptmass
 ndustlarge = hdr%ndustlarge
 ndustsmall = hdr%ndustsmall
 npart = hdr%nparttot
 npartoftype = reshape(hdr%npartoftypetot,shape(npartoftype),pad=[0])
 iexternalforce = hdr%iexternalforce
 ieos = hdr%ieos
 tfile = hdr%time
 dtmax = hdr%dtmax
 gamma = hdr%gamma
 rhozero = hdr%rhozero
 polyk = hdr%polyk
 hfactfile = hdr%hfact
 tolh = hdr%tolh
 C_cour = hdr%C_cour
 C_force = hdr%C_force
 alpha = hdr%alpha
 alphau = hdr%alphau
 alphaB = hdr%alphaB
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
 grainsize = reshape(hdr%grainsize,shape(grainsize),pad=[0.0])
 graindens = reshape(hdr%graindens,shape(graindens),pad=[0.0])
 udist = hdr%udist
 umass = hdr%umass
 utime = hdr%utime
 unit_Bfield = hdr%unit_Bfield

 call set_units_extra()
 ndusttypes = ndustsmall + ndustlarge

 call get_options_from_fileid(fileident,smalldump,use_dustfrac,errors(3))

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
 call allocate_memory(int(npart / nprocs))
#endif

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
 array_options%store_temperature = store_temperature

 if (smalldump) then
    call read_hdf5_arrays(hdf5_file_id,errors(4),npart,nptmass,iphase,xyzh,    &
                          vxyzu,xyzmh_ptmass,vxyz_ptmass,dtind,alphaind,poten, &
                          Bxyz,Bevol,dustfrac,deltav,dustprop,tstop,St,        &
                          temperature,abundance,array_options,got_arrays)

#ifdef IND_TIMESTEPS
    if (size(dt_in)/=size(dtind)) call error('read_smalldump','problem reading individual timesteps')
    dt_in = dtind
#endif

 else
    call error('read_smalldump',trim(dumpfile)//'.h5 is not a small dump')
 endif

 call close_hdf5file(hdf5_file_id,errors(5))

 ierr = maxval(abs(errors))

end subroutine read_smalldump

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

!---------------------------------------------------------------
!+
!  make sure required arrays have been read from Phantom file
!  and perform basic sanity checks
!+
!---------------------------------------------------------------
subroutine check_arrays(i1,i2,npartoftype,nptmass,nsinkproperties,massoftype, &
                        alphafile,tfile,got_iphase,got_xyzh,got_vxyzu,        &
                        got_alpha,got_abund,got_dustfrac,got_sink_data,       &
                        got_sink_vels,got_Bxyz,got_psi,got_dustprop,got_St,   &
                        got_temp,iphase,xyzh,vxyzu,alphaind,xyzmh_ptmass,     &
                        Bevol,iprint,ierr)

 use dim,        only:maxp,maxvxyzu,maxalpha,maxBevol,mhd,h2chemistry, &
                      store_temperature,use_dustgrowth
 use eos,        only:polyk,gamma
 use part,       only:maxphase,isetphase,igas,ihacc,ihsoft,imacc, &
                      xyzmh_ptmass_label,get_pmass,rhoh,dustfrac
 use io,         only:warning,id,master
 use options,    only:alpha,use_dustfrac

 integer,         intent(in)    :: i1,i2,npartoftype(:),nptmass,nsinkproperties
 real,            intent(in)    :: massoftype(:),alphafile,tfile
 logical,         intent(in)    :: got_iphase,got_xyzh,got_vxyzu,  &
                                   got_alpha,got_dustprop(:),      &
                                   got_St,got_abund,got_dustfrac,  &
                                   got_sink_data(:),got_sink_vels, &
                                   got_Bxyz,got_psi,got_temp
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
       do i=i1,i2
          itype = int(iphase(i))
          iphase(i) = isetphase(itype,iactive=.true.)
       enddo
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
 if (.not.got_xyzh) then
    if (id==master .and. i1==1) write(*,*) 'ERROR: x, y, z or h not found in file'
    ierr = 9
    return
 endif
 if (.not.got_vxyzu) then
    if (id==master .and. i1==1) write(*,*) 'ERROR: missing velocity information from file'
 endif
 if (maxvxyzu==4 .and. .not.got_vxyzu) then
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
 if (h2chemistry .and. .not.got_abund) then
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
 do i = 1, size(massoftype)
    if (npartoftype(i) > 0) then
       if (massoftype(i) <= 0.0) then
          if (id==master .and. i1==1) write(*,*) 'ERROR! mass not set in read_dump (Phantom)'
          ierr = 12
          return
       endif
    endif
 enddo
 if (use_dustfrac .and. .not. got_dustfrac) then
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
    if (.not.got_sink_vels) then
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
    if (.not.got_Bxyz) then
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

#ifdef PHANTOM2HDF5
end module readwrite_dumps_hdf5
#else
end module readwrite_dumps
#endif
