!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: utils_dumpfiles_hdf5
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: utils_hdf5
!+
!--------------------------------------------------------------------------
module utils_dumpfiles_hdf5
 use utils_hdf5, only:write_to_hdf5,    &
                      read_from_hdf5,   &
                      create_hdf5file,  &
                      open_hdf5file,    &
                      close_hdf5file,   &
                      open_hdf5group,   &
                      create_hdf5group, &
                      close_hdf5group,  &
                      HID_T

 implicit none

 public :: write_hdf5_header,write_hdf5_arrays,write_hdf5_arrays_small
 public :: read_hdf5_header,read_hdf5_arrays
 public :: create_hdf5file,open_hdf5file,close_hdf5file
 public :: header_hdf5,got_arrays_hdf5,arrays_options_hdf5,externalforce_hdf5

 integer(HID_T), public :: hdf5_file_id

 ! Ideally this section should come from Phantom modules.
 ! However, this module aims to be a library with no non-HDF5 dependencies.
 integer, parameter :: maxdustlarge_hdf5 = 50
 integer, parameter :: maxdustsmall_hdf5 = 50
 integer, parameter :: maxdusttypes_hdf5 = maxdustsmall_hdf5 + maxdustlarge_hdf5
 integer, parameter :: maxtypes_hdf5 = 7 + maxdustlarge_hdf5 - 1
 integer, parameter :: nsinkproperties_hdf5 = 11
 integer, parameter :: iext_binary_hdf5 = 3
 integer, parameter :: iext_gwinspiral_hdf5 = 14
 integer, parameter :: iext_corot_binary_hdf5 = 16

 type header_hdf5
    character(len=100) :: fileident
    integer            :: ntypes,                        &
                          isink,                         &
                          nptmass,                       &
                          ndustlarge,                    &
                          ndustsmall,                    &
                          idust,                         &
                          phantom_version_major,         &
                          phantom_version_minor,         &
                          phantom_version_micro,         &
                          nparttot,                      &
                          npartoftypetot(maxtypes_hdf5), &
                          iexternalforce,                &
                          ieos
    real               :: time,                          &
                          dtmax,                         &
                          gamma,                         &
                          rhozero,                       &
                          polyk,                         &
                          hfact,                         &
                          tolh,                          &
                          C_cour,                        &
                          C_force,                       &
                          alpha,                         &
                          alphau,                        &
                          alphaB,                        &
                          polyk2,                        &
                          qfacdisc,                      &
                          massoftype(maxtypes_hdf5),     &
                          Bextx,                         &
                          Bexty,                         &
                          Bextz,                         &
                          xmin,                          &
                          xmax,                          &
                          ymin,                          &
                          ymax,                          &
                          zmin,                          &
                          zmax,                          &
                          get_conserv,                   &
                          etot_in,                       &
                          angtot_in,                     &
                          totmom_in,                     &
                          mdust_in(maxdusttypes_hdf5),   &
                          grainsize(maxdusttypes_hdf5),  &
                          graindens(maxdusttypes_hdf5),  &
                          udist,                         &
                          umass,                         &
                          utime,                         &
                          unit_Bfield
 end type

 type externalforce_hdf5
    ! extern_binary
    real    :: xyzmh1(5),     &
               xyzmh2(5),     &
               vxyz1(3),      &
               vxyz2(3),      &
               accretedmass1, &
               accretedmass2, &
               a0,            &
               direction
    ! extern_gwinspiral
    integer :: Nstar(2)
 end type

 type got_arrays_hdf5
    logical :: got_iphase,                          &
               got_xyzh,                            &
               got_vxyzu,                           &
               got_dustfrac,                        &
               got_tstop,                           &
               got_deltav,                          &
               got_abund,                           &
               got_dt_in,                           &
               got_alpha,                           &
               got_poten,                           &
               got_sink_data(nsinkproperties_hdf5), &
               got_sink_vels,                       &
               got_Bxyz,                            &
               got_psi,                             &
               got_temp,                            &
               got_dustprop(2),                     &
               got_St,                              &
               got_VrelVf
 end type

 type arrays_options_hdf5
    logical :: isothermal,        &
               const_av,          &
               ind_timesteps,     &
               gravity,           &
               mhd,               &
               mhd_nonideal,      &
               use_dust,          &
               use_dustfrac,      &
               use_dustgrowth,    &
               h2chemistry,       &
               lightcurve,        &
               prdrag,            &
               store_temperature
    integer :: maxBevol,          &
               ndivcurlB,         &
               ndivcurlv,         &
               ndustsmall,        &
               ndustlarge
 end type

 private

contains


!--------------------------------------------------------------------
!+
!  write header
!+
!-------------------------------------------------------------------
subroutine write_hdf5_header(file_id,hdr,extern,error)

 integer(HID_T),            intent(in)  :: file_id
 type (header_hdf5),        intent(in)  :: hdr
 type (externalforce_hdf5), intent(in)  :: extern
 integer,                   intent(out) :: error

 integer(HID_T) :: group_id
 integer :: errors(70)

 errors(:) = 0

 ! Create header group
 call create_hdf5group(file_id,'header',group_id,errors(1))

 ! Write things to header group
 call write_to_hdf5(hdr%fileident,'fileident',group_id,errors(2))
 call write_to_hdf5(hdr%ntypes,'ntypes',group_id,errors(3))
 call write_to_hdf5(hdr%isink,'isink',group_id,errors(5))
 call write_to_hdf5(hdr%nptmass,'nptmass',group_id,errors(6))
 call write_to_hdf5(hdr%ndustlarge,'ndustlarge',group_id,errors(7))
 call write_to_hdf5(hdr%ndustsmall,'ndustsmall',group_id,errors(8))
 call write_to_hdf5(hdr%idust,'idust',group_id,errors(9))
 call write_to_hdf5(hdr%phantom_version_major,'majorv',group_id,errors(10))
 call write_to_hdf5(hdr%phantom_version_minor,'minorv',group_id,errors(11))
 call write_to_hdf5(hdr%phantom_version_micro,'microv',group_id,errors(12))
 call write_to_hdf5(hdr%nparttot,'nparttot',group_id,errors(13))
 call write_to_hdf5(hdr%npartoftypetot(1:maxtypes_hdf5),'npartoftype',group_id,errors(14))
 call write_to_hdf5(hdr%iexternalforce,'iexternalforce',group_id,errors(15))
 call write_to_hdf5(hdr%ieos,'ieos',group_id,errors(16))
 call write_to_hdf5(hdr%time,'time',group_id,errors(17))
 call write_to_hdf5(hdr%dtmax,'dtmax',group_id,errors(18))
 call write_to_hdf5(hdr%gamma,'gamma',group_id,errors(19))
 call write_to_hdf5(hdr%rhozero,'rhozero',group_id,errors(20))
 call write_to_hdf5(1.5*hdr%polyk,'RK2',group_id,errors(21))
 call write_to_hdf5(hdr%hfact,'hfact',group_id,errors(22))
 call write_to_hdf5(hdr%tolh,'tolh',group_id,errors(23))
 call write_to_hdf5(hdr%C_cour,'C_cour',group_id,errors(24))
 call write_to_hdf5(hdr%C_force,'C_force',group_id,errors(25))
 call write_to_hdf5(hdr%alpha,'alpha',group_id,errors(26))
 call write_to_hdf5(hdr%alphau,'alphau',group_id,errors(27))
 call write_to_hdf5(hdr%alphaB,'alphaB',group_id,errors(28))
 call write_to_hdf5(hdr%polyk2,'polyk2',group_id,errors(29))
 call write_to_hdf5(hdr%qfacdisc,'qfacdisc',group_id,errors(30))
 call write_to_hdf5(hdr%massoftype,'massoftype',group_id,errors(31))
 call write_to_hdf5(hdr%Bextx,'Bextx',group_id,errors(32))
 call write_to_hdf5(hdr%Bexty,'Bexty',group_id,errors(33))
 call write_to_hdf5(hdr%Bextz,'Bextz',group_id,errors(34))
 call write_to_hdf5(0.,'dum',group_id,errors(35))
 call write_to_hdf5(hdr%xmin,'xmin',group_id,errors(36))
 call write_to_hdf5(hdr%xmax,'xmax',group_id,errors(37))
 call write_to_hdf5(hdr%ymin,'ymin',group_id,errors(38))
 call write_to_hdf5(hdr%ymax,'ymax',group_id,errors(39))
 call write_to_hdf5(hdr%zmin,'zmin',group_id,errors(40))
 call write_to_hdf5(hdr%zmax,'zmax',group_id,errors(41))
 call write_to_hdf5(hdr%get_conserv,'get_conserv',group_id,errors(42))
 call write_to_hdf5(hdr%etot_in,'etot_in',group_id,errors(43))
 call write_to_hdf5(hdr%angtot_in,'angtot_in',group_id,errors(44))
 call write_to_hdf5(hdr%totmom_in,'totmom_in',group_id,errors(45))
 call write_to_hdf5(hdr%mdust_in,'mdust_in',group_id,errors(46))
 call write_to_hdf5(hdr%grainsize,'grainsize',group_id,errors(47))
 call write_to_hdf5(hdr%graindens,'graindens',group_id,errors(48))
 call write_to_hdf5(hdr%udist,'udist',group_id,errors(49))
 call write_to_hdf5(hdr%umass,'umass',group_id,errors(50))
 call write_to_hdf5(hdr%utime,'utime',group_id,errors(51))
 call write_to_hdf5(hdr%unit_Bfield,'umagfd',group_id,errors(52))

 ! Close the header group
 call close_hdf5group(group_id, errors(53))

 ! Write external force information
 select case(hdr%iexternalforce)
 case(iext_gwinspiral_hdf5)

    ! Create externalforce group
    call create_hdf5group(file_id,'header/externalforce',group_id,errors(54))

    ! Write values
    call write_to_hdf5(extern%Nstar,'Nstar',group_id,errors(55))

    ! Close the externalforce group
    call close_hdf5group(group_id,errors(56))

 case(iext_binary_hdf5,iext_corot_binary_hdf5)

    ! Create externalforce group
    call create_hdf5group(file_id,'header/externalforce',group_id,errors(54))

    ! Write values
    call write_to_hdf5(extern%xyzmh1,'xyzmh1',group_id,errors(55))
    call write_to_hdf5(extern%xyzmh2,'xyzmh2',group_id,errors(56))
    call write_to_hdf5(extern%vxyz1,'vxyz1',group_id,errors(57))
    call write_to_hdf5(extern%vxyz2,'vxyz2',group_id,errors(58))
    call write_to_hdf5(extern%accretedmass1,'accretedmass1',group_id,errors(59))
    call write_to_hdf5(extern%accretedmass2,'accretedmass2',group_id,errors(60))
    call write_to_hdf5(extern%a0,'a0',group_id,errors(61))
    call write_to_hdf5(extern%direction,'direction',group_id,errors(62))

    ! Close the externalforce group
    call close_hdf5group(group_id,errors(63))

 end select

 error = maxval(abs(errors))

end subroutine write_hdf5_header

!--------------------------------------------------------------------
!+
!  write arrays for full dump
!+
!-------------------------------------------------------------------
subroutine write_hdf5_arrays(file_id,error,npart,nptmass,xyzh,vxyzu,iphase,   &
                             pressure,alphaind,dtind,poten,xyzmh_ptmass,      &
                             vxyz_ptmass,Bxyz,Bevol,divcurlB,divBsymm,        &
                             eta_nimhd,dustfrac,tstop,deltav,dustprop,VrelVf, &
                             dustgasprop,abundance,temperature,divcurlv,     &
                             luminosity,beta_pr,array_options)

 integer(HID_T),  intent(in) :: file_id
 integer,         intent(out):: error
 integer,         intent(in) :: npart,nptmass
 real,            intent(in) :: pressure(:),dtind(:),beta_pr(:),VrelVf(:),     &
                                temperature(:),xyzh(:,:),vxyzu(:,:),Bxyz(:,:), &
                                Bevol(:,:), eta_nimhd(:,:),xyzmh_ptmass(:,:),  &
                                vxyz_ptmass(:,:),dustfrac(:,:),tstop(:,:),     &
                                dustprop(:,:),dustgasprop(:,:),                &
                                abundance(:,:),deltav(:,:,:)
 real(kind=4),    intent(in) :: poten(:),divBsymm(:),luminosity(:),            &
                                alphaind(:,:),divcurlv(:,:),divcurlB(:,:)
 integer(kind=1), intent(in) :: iphase(:)
 type (arrays_options_hdf5), intent(in)  :: array_options

 integer(HID_T) :: group_id
 integer :: errors(44),ndusttypes

 errors(:) = 0

 ! Create particles group
 call create_hdf5group(file_id,'particles',group_id,errors(1))

 ! Main arrays
 call write_to_hdf5(xyzh(1:3,1:npart),'xyz',group_id,errors(2))
 ! Write smoothing length in single precision to save disc space
 call write_to_hdf5(real(xyzh(4,1:npart),kind=4),'h',group_id,errors(3))
 call write_to_hdf5(vxyzu(1:3,1:npart),'vxyz',group_id,errors(4))
 if (.not.array_options%isothermal) call write_to_hdf5(vxyzu(4,1:npart),'u',group_id,errors(5))
 call write_to_hdf5(iphase(1:npart),'itype',group_id,errors(6))
 ! call write_to_hdf5(pressure(1:npart),'pressure',group_id,errors(7))

 if (.not.array_options%const_av)  call write_to_hdf5(alphaind(1,1:npart),'alpha',group_id,errors(8))   ! Viscosity (only ever write 'first' alpha)
 if (array_options%ind_timesteps)  call write_to_hdf5(real(dtind(1:npart),kind=4),'dt',group_id,errors(9)) ! Individual timesteps
 if (array_options%gravity)        call write_to_hdf5(poten(1:npart),'poten',group_id,errors(10))

 ! MHD arrays
 if (array_options%mhd) then
    call write_to_hdf5(Bxyz(:,1:npart),'Bxyz',group_id,errors(11))
    if (array_options%maxBevol >= 4) then
       call write_to_hdf5(Bevol(4,1:npart),'psi',group_id,errors(12))
    endif
    if (array_options%ndivcurlB >= 1) then
       call write_to_hdf5(divcurlB(1,1:npart),'divB',group_id,errors(13))
       call write_to_hdf5(divcurlB(2:4,1:npart),'curlBxyz',group_id,errors(14))
    else
       call write_to_hdf5(divBsymm(1:npart),'divBsymm',group_id,errors(16))
    endif
    if (array_options%mhd_nonideal) then
       call write_to_hdf5(eta_nimhd(1,1:npart),'eta_OR', group_id,errors(17))
       call write_to_hdf5(eta_nimhd(2,1:npart),'eta_HE', group_id,errors(18))
       call write_to_hdf5(eta_nimhd(3,1:npart),'eta_AD', group_id,errors(19))
       call write_to_hdf5(eta_nimhd(4,1:npart),'ne_on_n',group_id,errors(20))
    endif
 endif

 ! Dust arrays
 ndusttypes = array_options%ndustsmall + array_options%ndustlarge
 if (array_options%use_dust .and. ndusttypes > 0) then
    call write_to_hdf5(dustfrac(1:ndusttypes,1:npart),'dustfrac',group_id,errors(21))
    call write_to_hdf5(tstop(1:ndusttypes,1:npart),'tstop',group_id,errors(22))
    if (array_options%use_dustfrac .and. array_options%ndustsmall > 0) then
       call write_to_hdf5(deltav(:,1:array_options%ndustsmall,1:npart),'deltavxyz',group_id,errors(23))
    endif
 endif
 if (array_options%use_dustgrowth) then
    call write_to_hdf5(dustprop(1,1:npart),'grainsize',    group_id,errors(24))
    call write_to_hdf5(dustprop(2,1:npart),'graindens',    group_id,errors(25))
    call write_to_hdf5(VrelVf(1:npart),'vrel_on_vfrag',group_id,errors(26))
    ! call write_to_hdf5(dustprop(4,:),'dv_dust',group_id,errors())
    call write_to_hdf5(dustgasprop(3,1:npart),'St',group_id,errors(27))
 endif

 ! Other Arrays
 if (array_options%h2chemistry)       call write_to_hdf5(abundance(:,1:npart),'abundance',group_id,errors(28))
 if (array_options%store_temperature) call write_to_hdf5(temperature(1:npart),'T',group_id,errors(29))
 if (array_options%ndivcurlv >= 1) then
    call write_to_hdf5(divcurlv(1,1:npart),'divv',group_id,errors(30))
    if (array_options%ndivcurlv>=4) call write_to_hdf5(divcurlv(2:4,1:npart),'curlvxyz',group_id,errors(31))
 endif
 if (array_options%lightcurve) call write_to_hdf5(luminosity(1:npart),'luminosity',group_id,errors(32))
 if (array_options%prdrag)     call write_to_hdf5(real(beta_pr(1:npart),kind=4),'beta_pr',group_id,errors(33))

 ! Close the particles group
 call close_hdf5group(group_id, errors(34))

 ! Create sink group
 call create_hdf5group(file_id,'sinks',group_id,errors(35))
 if (nptmass > 0) then
    call write_to_hdf5(xyzmh_ptmass(1:3,1:nptmass),'xyz',group_id,errors(36))
    call write_to_hdf5(xyzmh_ptmass(4,1:nptmass),'m',group_id,errors(37))
    call write_to_hdf5(xyzmh_ptmass(5,1:nptmass),'h',group_id,errors(38))
    call write_to_hdf5(xyzmh_ptmass(6,1:nptmass),'hsoft',group_id,errors(39))
    call write_to_hdf5(xyzmh_ptmass(7,1:nptmass),'maccreted',group_id,errors(40))
    call write_to_hdf5(xyzmh_ptmass(8:10,1:nptmass),'spinxyz',group_id,errors(41))
    call write_to_hdf5(xyzmh_ptmass(11,1:nptmass),'tlast',group_id,errors(42))
    call write_to_hdf5(vxyz_ptmass(:,1:nptmass),'vxyz',group_id,errors(43))
 endif
 ! Close the sink group
 call close_hdf5group(group_id, errors(44))

 error = maxval(abs(errors))

end subroutine write_hdf5_arrays

!--------------------------------------------------------------------
!+
!  write arrays for small dump
!+
!-------------------------------------------------------------------
subroutine write_hdf5_arrays_small(file_id,error,npart,nptmass,xyzh,iphase,       &
                                   xyzmh_ptmass,Bxyz,dustfrac,dustprop,VrelVf,    &
                                   dustgasprop,abundance,luminosity,array_options)

 integer(HID_T),  intent(in) :: file_id
 integer,         intent(out):: error
 integer,         intent(in) :: npart,nptmass
 real,            intent(in) :: VrelVf(:),xyzh(:,:),Bxyz(:,:),xyzmh_ptmass(:,:), &
                                dustfrac(:,:),dustprop(:,:),dustgasprop(:,:),abundance(:,:)
 real(kind=4),    intent(in) :: luminosity(:)
 integer(kind=1), intent(in) :: iphase(:)
 type (arrays_options_hdf5), intent(in)  :: array_options

 integer(HID_T) :: group_id
 integer :: errors(22),ndusttypes

 errors(:) = 0

 ! Create particles group
 call create_hdf5group(file_id,'particles',group_id,errors(1))

 ! Main arrays
 call write_to_hdf5(real(xyzh(1:3,1:npart),kind=4),'xyz',group_id,errors(2))
 call write_to_hdf5(real(xyzh(4,1:npart),kind=4),'h',group_id,errors(3))
 call write_to_hdf5(iphase(1:npart),'itype',group_id,errors(4))

 ! MHD arrays
 if (array_options%mhd) then
    call write_to_hdf5(real(Bxyz(:,1:npart),kind=4),'Bxyz',group_id,errors(5))
 endif

 ! Dust arrays
 ndusttypes = array_options%ndustsmall + array_options%ndustlarge
 if (array_options%use_dust .and. ndusttypes > 0) then
    call write_to_hdf5(real(dustfrac(:,1:npart),kind=4),'dustfrac',group_id,errors(6))
 endif
 if (array_options%use_dustgrowth) then
    call write_to_hdf5(real(dustprop(1,1:npart),kind=4),'grainsize',    group_id,errors(7))
    call write_to_hdf5(real(dustprop(2,1:npart),kind=4),'graindens',    group_id,errors(8))
    call write_to_hdf5(real(VrelVf(1:npart),kind=4),'vrel_on_vfrag',group_id,errors(9))
    ! call write_to_hdf5(real(dustprop(4,1:npart),kind=4),'dv_dust',group_id,errors())
    call write_to_hdf5(real(dustgasprop(3,1:npart),kind=4),'St',group_id,errors(10))
 endif

 ! Other Arrays
 if (array_options%h2chemistry) call write_to_hdf5(real(abundance(:,1:npart),kind=4),'abundance',group_id,errors(11))
 if (array_options%lightcurve)  call write_to_hdf5(luminosity(1:npart),'luminosity',group_id,errors(12))

 ! Close the particles group
 call close_hdf5group(group_id, errors(13))

 ! Create sink group
 call create_hdf5group(file_id,'sinks',group_id,errors(14))
 if (nptmass > 0) then
    call write_to_hdf5(real(xyzmh_ptmass(1:3,1:nptmass),kind=4),'xyz',group_id,errors(15))
    call write_to_hdf5(real(xyzmh_ptmass(4,1:nptmass),kind=4),'m',group_id,errors(16))
    call write_to_hdf5(real(xyzmh_ptmass(5,1:nptmass),kind=4),'h',group_id,errors(17))
    call write_to_hdf5(real(xyzmh_ptmass(6,1:nptmass),kind=4),'hsoft',group_id,errors(18))
    call write_to_hdf5(real(xyzmh_ptmass(7,1:nptmass),kind=4),'maccreted',group_id,errors(19))
    call write_to_hdf5(real(xyzmh_ptmass(8:10,1:nptmass),kind=4),'spinxyz',group_id,errors(20))
    call write_to_hdf5(real(xyzmh_ptmass(11,1:nptmass),kind=4),'tlast',group_id,errors(21))
 endif
 ! Close the sink group
 call close_hdf5group(group_id, errors(22))

 error = maxval(abs(errors))

end subroutine write_hdf5_arrays_small

!--------------------------------------------------------------------
!+
!  read header
!+
!-------------------------------------------------------------------
subroutine read_hdf5_header(file_id,hdr,extern,error)

 integer(HID_T),            intent(in)  :: file_id
 type (header_hdf5),        intent(out) :: hdr
 type (externalforce_hdf5), intent(out) :: extern
 integer,                   intent(out) :: error

 integer(HID_T) :: group_id
 integer        :: errors(70)
 real           :: rval
 logical        :: got_val

 errors(:) = 0

 ! Open header group
 call open_hdf5group(file_id,'header',group_id,errors(1))

 ! Write things to header group
 call read_from_hdf5(hdr%fileident,'fileident',group_id,got_val,errors(2))
 call read_from_hdf5(hdr%ntypes,'ntypes',group_id,got_val,errors(3))
 call read_from_hdf5(hdr%isink,'isink',group_id,got_val,errors(5))
 call read_from_hdf5(hdr%nptmass,'nptmass',group_id,got_val,errors(6))
 call read_from_hdf5(hdr%ndustlarge,'ndustlarge',group_id,got_val,errors(7))
 call read_from_hdf5(hdr%ndustsmall,'ndustsmall',group_id,got_val,errors(8))
 call read_from_hdf5(hdr%nparttot,'nparttot',group_id,got_val,errors(13))
 call read_from_hdf5(hdr%npartoftypetot(1:hdr%ntypes),'npartoftype',group_id,got_val,errors(14))
 call read_from_hdf5(hdr%iexternalforce,'iexternalforce',group_id,got_val,errors(15))
 call read_from_hdf5(hdr%ieos,'ieos',group_id,got_val,errors(16))
 call read_from_hdf5(hdr%time,'time',group_id,got_val,errors(17))
 call read_from_hdf5(hdr%dtmax,'dtmax',group_id,got_val,errors(18))
 call read_from_hdf5(hdr%gamma,'gamma',group_id,got_val,errors(19))
 call read_from_hdf5(hdr%rhozero,'rhozero',group_id,got_val,errors(20))
 call read_from_hdf5(rval,'RK2',group_id,got_val,errors(21))
 hdr%polyk = rval/1.5
 call read_from_hdf5(hdr%hfact,'hfact',group_id,got_val,errors(22))
 call read_from_hdf5(hdr%tolh,'tolh',group_id,got_val,errors(23))
 call read_from_hdf5(hdr%C_cour,'C_cour',group_id,got_val,errors(24))
 call read_from_hdf5(hdr%C_force,'C_force',group_id,got_val,errors(25))
 call read_from_hdf5(hdr%alpha,'alpha',group_id,got_val,errors(26))
 call read_from_hdf5(hdr%alphau,'alphau',group_id,got_val,errors(27))
 call read_from_hdf5(hdr%alphaB,'alphaB',group_id,got_val,errors(28))
 call read_from_hdf5(hdr%polyk2,'polyk2',group_id,got_val,errors(29))
 call read_from_hdf5(hdr%qfacdisc,'qfacdisc',group_id,got_val,errors(30))
 call read_from_hdf5(hdr%massoftype,'massoftype',group_id,got_val,errors(31))
 call read_from_hdf5(hdr%Bextx,'Bextx',group_id,got_val,errors(32))
 call read_from_hdf5(hdr%Bexty,'Bexty',group_id,got_val,errors(33))
 call read_from_hdf5(hdr%Bextz,'Bextz',group_id,got_val,errors(34))
 call read_from_hdf5(hdr%xmin,'xmin',group_id,got_val,errors(36))
 call read_from_hdf5(hdr%xmax,'xmax',group_id,got_val,errors(37))
 call read_from_hdf5(hdr%ymin,'ymin',group_id,got_val,errors(38))
 call read_from_hdf5(hdr%ymax,'ymax',group_id,got_val,errors(39))
 call read_from_hdf5(hdr%zmin,'zmin',group_id,got_val,errors(40))
 call read_from_hdf5(hdr%zmax,'zmax',group_id,got_val,errors(41))
 call read_from_hdf5(hdr%get_conserv,'get_conserv',group_id,got_val,errors(42))
 call read_from_hdf5(hdr%etot_in,'etot_in',group_id,got_val,errors(43))
 call read_from_hdf5(hdr%angtot_in,'angtot_in',group_id,got_val,errors(44))
 call read_from_hdf5(hdr%totmom_in,'totmom_in',group_id,got_val,errors(45))
 call read_from_hdf5(hdr%mdust_in,'mdust_in',group_id,got_val,errors(46))
 call read_from_hdf5(hdr%grainsize,'grainsize',group_id,got_val,errors(47))
 call read_from_hdf5(hdr%graindens,'graindens',group_id,got_val,errors(48))
 call read_from_hdf5(hdr%udist,'udist',group_id,got_val,errors(49))
 call read_from_hdf5(hdr%umass,'umass',group_id,got_val,errors(50))
 call read_from_hdf5(hdr%utime,'utime',group_id,got_val,errors(51))
 call read_from_hdf5(hdr%unit_Bfield,'umagfd',group_id,got_val,errors(52))

 ! Close the header group
 call close_hdf5group(group_id,errors(53))

 ! Read external force information
 select case(hdr%iexternalforce)
 case(iext_gwinspiral_hdf5)

    ! Open externalforce group
    call open_hdf5group(file_id,'header/externalforce',group_id,errors(54))

    ! Read values
    call read_from_hdf5(extern%Nstar,'Nstar',group_id,got_val,errors(55))

    ! Close the externalforce group
    call close_hdf5group(group_id,errors(56))

 case(iext_binary_hdf5,iext_corot_binary_hdf5)

    ! Create externalforce group
    call open_hdf5group(file_id,'header/externalforce',group_id,errors(54))

    ! Write values
    call read_from_hdf5(extern%xyzmh1,'xyzmh1',group_id,got_val,errors(55))
    call read_from_hdf5(extern%xyzmh2,'xyzmh2',group_id,got_val,errors(56))
    call read_from_hdf5(extern%vxyz1,'vxyz1',group_id,got_val,errors(57))
    call read_from_hdf5(extern%vxyz2,'vxyz2',group_id,got_val,errors(58))
    call read_from_hdf5(extern%accretedmass1,'accretedmass1',group_id,got_val,errors(59))
    call read_from_hdf5(extern%accretedmass2,'accretedmass2',group_id,got_val,errors(60))
    call read_from_hdf5(extern%a0,'a0',group_id,got_val,errors(61))
    call read_from_hdf5(extern%direction,'direction',group_id,got_val,errors(62))

    ! Close the externalforce group
    call close_hdf5group(group_id,errors(63))

 end select

 error = maxval(abs(errors))

end subroutine read_hdf5_header

!--------------------------------------------------------------------
!+
!  read arrays for full dump
!+
!-------------------------------------------------------------------
subroutine read_hdf5_arrays(file_id,error,npart,nptmass,iphase,xyzh,vxyzu,   &
                            xyzmh_ptmass,vxyz_ptmass,dt_in,alphaind,poten,   &
                            Bxyz,Bevol,dustfrac,deltav,dustprop,tstop,       &
                            VrelVf,dustgasprop,temperature,abundance,        &
                            array_options,got_arrays)

 integer(HID_T),  intent(in)  :: file_id
 integer,         intent(in)  :: npart,nptmass
 type (arrays_options_hdf5), intent(in)  :: array_options
 type (got_arrays_hdf5),     intent(out) :: got_arrays
 integer(kind=1), intent(out) :: iphase(:)
 real,            intent(out) :: xyzh(:,:),                      &
                                 vxyzu(:,:),                     &
                                 xyzmh_ptmass(:,:),              &
                                 vxyz_ptmass(:,:),               &
                                 Bxyz(:,:),                      &
                                 Bevol(:,:),                     &
                                 dustfrac(:,:),                  &
                                 deltav(:,:,:),                  &
                                 dustprop(:,:),                  &
                                 dustgasprop(:,:),               &
                                 tstop(:,:),                     &
                                 VrelVf(:),                      &
                                 temperature(:),                 &
                                 abundance(:,:)
 real(kind=4),    intent(out) :: dt_in(:),                       &
                                 alphaind(:,:),                  &
                                 poten(:)
 integer,         intent(out) :: error

 integer(HID_T) :: group_id
 integer :: errors(44)
 logical :: got

 real(kind=4) :: rtmp(npart)

 errors(:) = 0

 got_arrays%got_iphase    = .false.
 got_arrays%got_xyzh      = .false.
 got_arrays%got_vxyzu     = .false.
 got_arrays%got_dustfrac  = .false.
 got_arrays%got_tstop     = .false.
 got_arrays%got_deltav    = .false.
 got_arrays%got_abund     = .false.
 got_arrays%got_dt_in     = .false.
 got_arrays%got_alpha     = .false.
 got_arrays%got_poten     = .false.
 got_arrays%got_sink_data = .false.
 got_arrays%got_sink_vels = .false.
 got_arrays%got_Bxyz      = .false.
 got_arrays%got_psi       = .false.
 got_arrays%got_temp      = .false.
 got_arrays%got_dustprop  = .false.
 got_arrays%got_St        = .false.
 got_arrays%got_VrelVf    = .false.

 ! Open particles group
 call open_hdf5group(file_id,'particles',group_id,errors(1))

 ! Main arrays
 call read_from_hdf5(iphase,'itype',group_id,got_arrays%got_iphase,errors(2))
 call read_from_hdf5(xyzh(1:3,:),'xyz',group_id,got,errors(3))
 if (got) got_arrays%got_xyzh = .true.
 call read_from_hdf5(rtmp,'h',group_id,got,errors(4))
 if (got) then
    xyzh(4,:) = real(rtmp)
 else
    got_arrays%got_xyzh = .false.
 endif
 call read_from_hdf5(vxyzu(1:3,:),'vxyz',group_id,got,errors(5))
 if (got) got_arrays%got_vxyzu = .true.
 if (.not.array_options%isothermal) then
    call read_from_hdf5(vxyzu(4,:),'u',group_id,got,errors(6))
    if (.not.got) got_arrays%got_vxyzu = .false.
 endif
 if (array_options%ind_timesteps)  call read_from_hdf5(dt_in,'dt',group_id,got_arrays%got_dt_in,errors(7))
 if (.not. array_options%const_av) call read_from_hdf5(alphaind(1,:),'alpha',group_id,got_arrays%got_alpha,errors(8))
 if (array_options%gravity)        call read_from_hdf5(poten,'poten',group_id,got_arrays%got_poten,errors(9))

 ! MHD arrays
 if (array_options%mhd) then
    call read_from_hdf5(Bxyz,'Bxyz',group_id,got_arrays%got_Bxyz,errors(10))
    call read_from_hdf5(Bevol(4,:),'psi',group_id,got_arrays%got_psi,errors(11))
 endif

 ! Dust arrays
 if (array_options%use_dust) then
    call read_from_hdf5(dustfrac,'dustfrac',group_id,got_arrays%got_dustfrac,errors(12))
    call read_from_hdf5(tstop,'tstop',group_id,got_arrays%got_tstop,errors(13))
 endif
 if (array_options%use_dustfrac) call read_from_hdf5(deltav,'deltavxyz',group_id,got_arrays%got_deltav,errors(14))
 if (array_options%use_dustgrowth) then
    call read_from_hdf5(dustprop(1,:),'grainsize',    group_id,got_arrays%got_dustprop(1),errors(15))
    call read_from_hdf5(dustprop(2,:),'graindens',    group_id,got_arrays%got_dustprop(2),errors(16))
    call read_from_hdf5(VrelVf(:),'vrel_on_vfrag',group_id,got_arrays%got_VrelVf,errors(17))
    ! call read_from_hdf5(dustprop(4,:),'dv_dust',group_id,got_dv_dust,errors())
    call read_from_hdf5(dustgasprop(3,:),'St',group_id,got_arrays%got_St,errors(18))
 endif

 ! Other Arrays
 if (array_options%h2chemistry) call read_from_hdf5(abundance,'abundance',group_id,got_arrays%got_abund,errors(19))
 if (array_options%store_temperature) call read_from_hdf5(temperature,'T',group_id,got_arrays%got_temp,errors(20))

 ! Close the particles group
 call close_hdf5group(group_id, errors(21))

 ! Open sinks group
 call open_hdf5group(file_id,'sinks',group_id,errors(22))

 ! Sink arrays
 if (nptmass > 0) then
    call read_from_hdf5(xyzmh_ptmass(1:3,1:nptmass),'xyz',group_id,got_arrays%got_sink_data(1),errors(23))
    got_arrays%got_sink_data(1:3) = got_arrays%got_sink_data(1)
    call read_from_hdf5(xyzmh_ptmass(4,1:nptmass),'m',group_id,got_arrays%got_sink_data(4),errors(24))
    call read_from_hdf5(xyzmh_ptmass(5,1:nptmass),'h',group_id,got_arrays%got_sink_data(5),errors(25))
    call read_from_hdf5(xyzmh_ptmass(6,1:nptmass),'hsoft',group_id,got_arrays%got_sink_data(6),errors(26))
    call read_from_hdf5(xyzmh_ptmass(7,1:nptmass),'maccreted',group_id,got_arrays%got_sink_data(7),errors(27))
    call read_from_hdf5(xyzmh_ptmass(8:10,1:nptmass),'spinxyz',group_id,got_arrays%got_sink_data(8),errors(28))
    got_arrays%got_sink_data(8:10) = got_arrays%got_sink_data(8)
    call read_from_hdf5(xyzmh_ptmass(11,1:nptmass),'tlast',group_id,got_arrays%got_sink_data(11),errors(29))
    call read_from_hdf5(vxyz_ptmass(:,1:nptmass),'vxyz',group_id,got_arrays%got_sink_vels,errors(30))
 endif

 ! Close the sinks group
 call close_hdf5group(group_id, errors(31))

 error = maxval(abs(errors))

end subroutine read_hdf5_arrays

end module utils_dumpfiles_hdf5
