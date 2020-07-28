!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module utils_dumpfiles_hdf5
!
! None
!
! :References: None
!
! :Owner: Daniel Mentiplay
!
! :Runtime parameters: None
!
! :Dependencies: utils_hdf5
!
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

 public :: write_hdf5_header,       &
           write_hdf5_arrays,       &
           write_hdf5_arrays_small, &
           read_hdf5_header,        &
           read_hdf5_arrays,        &
           create_hdf5file,         &
           open_hdf5file,           &
           close_hdf5file,          &
           header_hdf5,             &
           got_arrays_hdf5,         &
           arrays_options_hdf5,     &
           externalforce_hdf5

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
                          idustbound,                    &
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
    integer :: ieos,              &
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
subroutine write_hdf5_header(file_id, hdr, extern, error)

 integer(HID_T),            intent(in)  :: file_id
 type (header_hdf5),        intent(in)  :: hdr
 type (externalforce_hdf5), intent(in)  :: extern
 integer,                   intent(out) :: error

 integer(HID_T) :: group_id

 error = 0

 ! Create header group
 call create_hdf5group(file_id, 'header', group_id, error)

 ! Write things to header group
 call write_to_hdf5(hdr%fileident, 'fileident', group_id, error)
 call write_to_hdf5(hdr%ntypes, 'ntypes', group_id, error)
 call write_to_hdf5(hdr%isink, 'isink', group_id, error)
 call write_to_hdf5(hdr%nptmass, 'nptmass', group_id, error)
 call write_to_hdf5(hdr%ndustlarge, 'ndustlarge', group_id, error)
 call write_to_hdf5(hdr%ndustsmall, 'ndustsmall', group_id, error)
 call write_to_hdf5(hdr%idust, 'idust', group_id, error)
 call write_to_hdf5(hdr%idustbound, 'idustbound', group_id, error)
 call write_to_hdf5(hdr%phantom_version_major, 'majorv', group_id, error)
 call write_to_hdf5(hdr%phantom_version_minor, 'minorv', group_id, error)
 call write_to_hdf5(hdr%phantom_version_micro, 'microv', group_id, error)
 call write_to_hdf5(hdr%nparttot, 'nparttot', group_id, error)
 call write_to_hdf5(hdr%npartoftypetot(1:maxtypes_hdf5), 'npartoftype', group_id, error)
 call write_to_hdf5(hdr%iexternalforce, 'iexternalforce', group_id, error)
 call write_to_hdf5(hdr%ieos, 'ieos', group_id, error)
 call write_to_hdf5(hdr%time, 'time', group_id, error)
 call write_to_hdf5(hdr%dtmax, 'dtmax', group_id, error)
 call write_to_hdf5(hdr%gamma, 'gamma', group_id, error)
 call write_to_hdf5(hdr%rhozero, 'rhozero', group_id, error)
 call write_to_hdf5(1.5*hdr%polyk, 'RK2', group_id, error)
 call write_to_hdf5(hdr%hfact, 'hfact', group_id, error)
 call write_to_hdf5(hdr%tolh, 'tolh', group_id, error)
 call write_to_hdf5(hdr%C_cour, 'C_cour', group_id, error)
 call write_to_hdf5(hdr%C_force, 'C_force', group_id, error)
 call write_to_hdf5(hdr%alpha, 'alpha', group_id, error)
 call write_to_hdf5(hdr%alphau, 'alphau', group_id, error)
 call write_to_hdf5(hdr%alphaB, 'alphaB', group_id, error)
 call write_to_hdf5(hdr%polyk2, 'polyk2', group_id, error)
 call write_to_hdf5(hdr%qfacdisc, 'qfacdisc', group_id, error)
 call write_to_hdf5(hdr%massoftype, 'massoftype', group_id, error)
 call write_to_hdf5(hdr%Bextx, 'Bextx', group_id, error)
 call write_to_hdf5(hdr%Bexty, 'Bexty', group_id, error)
 call write_to_hdf5(hdr%Bextz, 'Bextz', group_id, error)
 call write_to_hdf5(0., 'dum', group_id, error)
 call write_to_hdf5(hdr%xmin, 'xmin', group_id, error)
 call write_to_hdf5(hdr%xmax, 'xmax', group_id, error)
 call write_to_hdf5(hdr%ymin, 'ymin', group_id, error)
 call write_to_hdf5(hdr%ymax, 'ymax', group_id, error)
 call write_to_hdf5(hdr%zmin, 'zmin', group_id, error)
 call write_to_hdf5(hdr%zmax, 'zmax', group_id, error)
 call write_to_hdf5(hdr%get_conserv, 'get_conserv', group_id, error)
 call write_to_hdf5(hdr%etot_in, 'etot_in', group_id, error)
 call write_to_hdf5(hdr%angtot_in, 'angtot_in', group_id, error)
 call write_to_hdf5(hdr%totmom_in, 'totmom_in', group_id, error)
 call write_to_hdf5(hdr%mdust_in, 'mdust_in', group_id, error)
 call write_to_hdf5(hdr%grainsize, 'grainsize', group_id, error)
 call write_to_hdf5(hdr%graindens, 'graindens', group_id, error)
 call write_to_hdf5(hdr%udist, 'udist', group_id, error)
 call write_to_hdf5(hdr%umass, 'umass', group_id, error)
 call write_to_hdf5(hdr%utime, 'utime', group_id, error)
 call write_to_hdf5(hdr%unit_Bfield, 'umagfd', group_id, error)

 ! Close the header group
 call close_hdf5group(group_id, error)

 ! Write external force information
 select case(hdr%iexternalforce)
 case(iext_gwinspiral_hdf5)

    ! Create externalforce group
    call create_hdf5group(file_id, 'header/externalforce', group_id, error)

    ! Write values
    call write_to_hdf5(extern%Nstar, 'Nstar', group_id, error)

    ! Close the externalforce group
    call close_hdf5group(group_id, error)

 case(iext_binary_hdf5, iext_corot_binary_hdf5)

    ! Create externalforce group
    call create_hdf5group(file_id, 'header/externalforce', group_id, error)

    ! Write values
    call write_to_hdf5(extern%xyzmh1, 'xyzmh1', group_id, error)
    call write_to_hdf5(extern%xyzmh2, 'xyzmh2', group_id, error)
    call write_to_hdf5(extern%vxyz1, 'vxyz1', group_id, error)
    call write_to_hdf5(extern%vxyz2, 'vxyz2', group_id, error)
    call write_to_hdf5(extern%accretedmass1, 'accretedmass1', group_id, error)
    call write_to_hdf5(extern%accretedmass2, 'accretedmass2', group_id, error)
    call write_to_hdf5(extern%a0, 'a0', group_id, error)
    call write_to_hdf5(extern%direction, 'direction', group_id, error)

    ! Close the externalforce group
    call close_hdf5group(group_id, error)

 end select

end subroutine write_hdf5_header

!--------------------------------------------------------------------
!+
!  write arrays for full dump
!+
!-------------------------------------------------------------------
subroutine write_hdf5_arrays( &
   file_id,                   &
   error,                     &
   npart,                     &
   nptmass,                   &
   xyzh,                      &
   vxyzu,                     &
   iphase,                    &
   pressure,                  &
   alphaind,                  &
   dtind,                     &
   poten,                     &
   xyzmh_ptmass,              &
   vxyz_ptmass,               &
   Bxyz,                      &
   Bevol,                     &
   divcurlB,                  &
   divBsymm,                  &
   eta_nimhd,                 &
   dustfrac,                  &
   tstop,                     &
   deltav,                    &
   dustprop,                  &
   VrelVf,                    &
   dustgasprop,               &
   abundance,                 &
   temperature,               &
   divcurlv,                  &
   luminosity,                &
   beta_pr,                   &
   array_options              &
)

 integer(HID_T),  intent(in) :: file_id
 integer,         intent(out):: error
 integer,         intent(in) :: npart, nptmass
 real,            intent(in) :: dtind(:),          &
                                beta_pr(:),        &
                                VrelVf(:),         &
                                temperature(:),    &
                                xyzh(:,:),         &
                                vxyzu(:,:),        &
                                Bxyz(:,:),         &
                                Bevol(:,:),        &
                                pressure(:),       &
                                eta_nimhd(:,:),    &
                                xyzmh_ptmass(:,:), &
                                vxyz_ptmass(:,:),  &
                                dustfrac(:,:),     &
                                tstop(:,:),        &
                                dustprop(:,:),     &
                                dustgasprop(:,:),  &
                                abundance(:,:),    &
                                deltav(:,:,:)
 real(kind=4),    intent(in) :: poten(:),          &
                                divBsymm(:),       &
                                luminosity(:),     &
                                alphaind(:,:),     &
                                divcurlv(:,:),     &
                                divcurlB(:,:)
 integer(kind=1), intent(in) :: iphase(:)
 type (arrays_options_hdf5), intent(in) :: array_options

 integer(HID_T) :: group_id
 integer :: ndusttypes,ieos

 error = 0
 ieos = array_options%ieos

 ! Create particles group
 call create_hdf5group(file_id, 'particles', group_id, error)

 ! Main arrays
 call write_to_hdf5(xyzh(1:3,1:npart), 'xyz', group_id, error)
 ! Write smoothing length in single precision to save disc space
 call write_to_hdf5(real(xyzh(4,1:npart), kind=4), 'h', group_id, error)
 call write_to_hdf5(vxyzu(1:3,1:npart), 'vxyz', group_id, error)
 if (.not.array_options%isothermal) call write_to_hdf5(vxyzu(4,1:npart), 'u', group_id, error)
 call write_to_hdf5(iphase(1:npart), 'itype', group_id, error)
 if (ieos==8 .or. ieos==9 .or. ieos==10 .or. ieos==15) then
    call write_to_hdf5(pressure(1:npart), 'pressure', group_id, error)
 endif

 ! Viscosity (only ever write 'first' alpha)
 if (.not.array_options%const_av) call write_to_hdf5(alphaind(1,1:npart), 'alpha', group_id, error)
 if (array_options%ind_timesteps) call write_to_hdf5(real(dtind(1:npart), kind=4), 'dt', group_id, error)
 if (array_options%gravity) call write_to_hdf5(poten(1:npart), 'poten', group_id, error)

 ! MHD arrays
 if (array_options%mhd) then
    call write_to_hdf5(Bxyz(:,1:npart), 'Bxyz', group_id, error)
    call write_to_hdf5(Bevol(4,1:npart), 'psi', group_id, error)
    if (array_options%ndivcurlB >= 1) then
       call write_to_hdf5(divcurlB(1,1:npart), 'divB', group_id, error)
       call write_to_hdf5(divcurlB(2:4,1:npart), 'curlBxyz', group_id, error)
    else
       call write_to_hdf5(divBsymm(1:npart), 'divBsymm', group_id, error)
    endif
    if (array_options%mhd_nonideal) then
       call write_to_hdf5(eta_nimhd(1,1:npart), 'eta_OR', group_id, error)
       call write_to_hdf5(eta_nimhd(2,1:npart), 'eta_HE', group_id, error)
       call write_to_hdf5(eta_nimhd(3,1:npart), 'eta_AD', group_id, error)
       call write_to_hdf5(eta_nimhd(4,1:npart), 'ne_on_n', group_id, error)
    endif
 endif

 ! Dust arrays
 ndusttypes = array_options%ndustsmall + array_options%ndustlarge
 if (array_options%use_dust .and. ndusttypes > 0) then
    call write_to_hdf5(dustfrac(1:ndusttypes,1:npart), 'dustfrac', group_id, error)
    call write_to_hdf5(tstop(1:ndusttypes,1:npart), 'tstop', group_id, error)
    if (array_options%use_dustfrac .and. array_options%ndustsmall > 0) then
       call write_to_hdf5(deltav(:,1:array_options%ndustsmall,1:npart), 'deltavxyz', group_id, error)
    endif
 endif
 if (array_options%use_dustgrowth) then
    call write_to_hdf5(dustprop(1,1:npart), 'grainsize', group_id, error)
    call write_to_hdf5(dustprop(2,1:npart), 'graindens', group_id, error)
    call write_to_hdf5(VrelVf(1:npart), 'vrel_on_vfrag', group_id, error)
    call write_to_hdf5(dustgasprop(3,1:npart), 'St', group_id, error)
 endif

 ! Other Arrays
 if (array_options%h2chemistry) call write_to_hdf5(abundance(:,1:npart), 'abundance', group_id, error)
 if (array_options%store_temperature) call write_to_hdf5(temperature(1:npart), 'T', group_id, error)
 if (array_options%ndivcurlv >= 1) then
    call write_to_hdf5(divcurlv(1,1:npart), 'divv', group_id, error)
    if (array_options%ndivcurlv>=4) call write_to_hdf5(divcurlv(2:4,1:npart), 'curlvxyz', group_id, error)
 endif
 if (array_options%lightcurve) call write_to_hdf5(luminosity(1:npart), 'luminosity', group_id, error)
 if (array_options%prdrag) call write_to_hdf5(real(beta_pr(1:npart), kind=4), 'beta_pr', group_id, error)

 ! Close the particles group
 call close_hdf5group(group_id, error)

 ! Create sink group
 call create_hdf5group(file_id, 'sinks', group_id, error)
 if (nptmass > 0) then
    call write_to_hdf5(xyzmh_ptmass(1:3,1:nptmass),'xyz',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(4,1:nptmass),'m',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(5,1:nptmass),'h',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(6,1:nptmass),'hsoft',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(7,1:nptmass),'maccreted',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(8:10,1:nptmass),'spinxyz',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(11,1:nptmass),'tlast',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(12,1:nptmass),'lum',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(13,1:nptmass),'Teff',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(14,1:nptmass),'Reff',group_id,error)
    call write_to_hdf5(xyzmh_ptmass(15,1:nptmass),'mdot',group_id,error)
    call write_to_hdf5(vxyz_ptmass(:,1:nptmass),'vxyz',group_id,error)
 endif

 ! Close the sink group
 call close_hdf5group(group_id, error)

end subroutine write_hdf5_arrays

!--------------------------------------------------------------------
!+
!  write arrays for small dump
!+
!-------------------------------------------------------------------
subroutine write_hdf5_arrays_small( &
   file_id,                         &
   error,                           &
   npart,                           &
   nptmass,                         &
   xyzh,                            &
   iphase,                          &
   xyzmh_ptmass,                    &
   Bxyz,                            &
   dustfrac,                        &
   dustprop,                        &
   VrelVf,                          &
   dustgasprop,                     &
   abundance,                       &
   luminosity,                      &
   array_options                    &
)

 integer(HID_T),  intent(in) :: file_id
 integer,         intent(out):: error
 integer,         intent(in) :: npart, nptmass
 real,            intent(in) :: xyzh(:,:),         &
                                Bxyz(:,:),         &
                                xyzmh_ptmass(:,:), &
                                dustfrac(:,:),     &
                                dustprop(:,:),     &
                                dustgasprop(:,:),  &
                                abundance(:,:),    &
                                VrelVf(:)
 real(kind=4),    intent(in) :: luminosity(:)
 integer(kind=1), intent(in) :: iphase(:)
 type (arrays_options_hdf5), intent(in)  :: array_options

 integer(HID_T) :: group_id
 integer :: ndusttypes

 error = 0

 ! Create particles group
 call create_hdf5group(file_id, 'particles', group_id, error)

 ! Main arrays
 call write_to_hdf5(real(xyzh(1:3,1:npart), kind=4), 'xyz', group_id, error)
 call write_to_hdf5(real(xyzh(4,1:npart), kind=4), 'h', group_id, error)
 call write_to_hdf5(iphase(1:npart), 'itype', group_id, error)

 ! MHD arrays
 if (array_options%mhd) then
    call write_to_hdf5(real(Bxyz(:,1:npart), kind=4), 'Bxyz', group_id, error)
 endif

 ! Dust arrays
 ndusttypes = array_options%ndustsmall + array_options%ndustlarge
 if (array_options%use_dust .and. ndusttypes > 0) then
    call write_to_hdf5(real(dustfrac(:,1:npart), kind=4), 'dustfrac', group_id, error)
 endif
 if (array_options%use_dustgrowth) then
    call write_to_hdf5(real(dustprop(1,1:npart), kind=4), 'grainsize', group_id, error)
    call write_to_hdf5(real(dustprop(2,1:npart), kind=4), 'graindens', group_id, error)
    call write_to_hdf5(real(VrelVf(1:npart), kind=4), 'vrel_on_vfrag', group_id, error)
    call write_to_hdf5(real(dustgasprop(3,1:npart), kind=4), 'St', group_id, error)
 endif

 ! Other Arrays
 if (array_options%h2chemistry) call write_to_hdf5(real(abundance(:,1:npart), kind=4), 'abundance', group_id, error)
 if (array_options%lightcurve) call write_to_hdf5(luminosity(1:npart), 'luminosity', group_id, error)

 ! Close the particles group
 call close_hdf5group(group_id, error)

 ! Create sink group
 call create_hdf5group(file_id, 'sinks', group_id, error)
 if (nptmass > 0) then
    call write_to_hdf5(real(xyzmh_ptmass(1:3,1:nptmass), kind=4), 'xyz', group_id, error)
    call write_to_hdf5(real(xyzmh_ptmass(4,1:nptmass), kind=4), 'm', group_id, error)
    call write_to_hdf5(real(xyzmh_ptmass(5,1:nptmass), kind=4), 'h', group_id, error)
    call write_to_hdf5(real(xyzmh_ptmass(6,1:nptmass), kind=4), 'hsoft', group_id, error)
    call write_to_hdf5(real(xyzmh_ptmass(7,1:nptmass), kind=4), 'maccreted', group_id, error)
    call write_to_hdf5(real(xyzmh_ptmass(8:10,1:nptmass), kind=4), 'spinxyz', group_id, error)
    call write_to_hdf5(real(xyzmh_ptmass(11,1:nptmass), kind=4), 'tlast', group_id, error)
 endif
 ! Close the sink group
 call close_hdf5group(group_id, error)

end subroutine write_hdf5_arrays_small

!--------------------------------------------------------------------
!+
!  read header
!+
!-------------------------------------------------------------------
subroutine read_hdf5_header(file_id, hdr, extern, error)

 integer(HID_T),            intent(in)  :: file_id
 type (header_hdf5),        intent(out) :: hdr
 type (externalforce_hdf5), intent(out) :: extern
 integer,                   intent(out) :: error

 integer(HID_T) :: group_id
 real           :: rval
 logical        :: got_val

 error = 0

 ! Open header group
 call open_hdf5group(file_id, 'header', group_id, error)

 ! Write things to header group
 call read_from_hdf5(hdr%fileident, 'fileident', group_id, got_val, error)
 call read_from_hdf5(hdr%ntypes, 'ntypes', group_id, got_val, error)
 call read_from_hdf5(hdr%isink, 'isink', group_id, got_val, error)
 call read_from_hdf5(hdr%nptmass, 'nptmass', group_id, got_val, error)
 call read_from_hdf5(hdr%ndustlarge, 'ndustlarge', group_id, got_val, error)
 call read_from_hdf5(hdr%ndustsmall, 'ndustsmall', group_id, got_val, error)
 call read_from_hdf5(hdr%nparttot, 'nparttot', group_id, got_val, error)
 call read_from_hdf5(hdr%npartoftypetot(1:hdr%ntypes), 'npartoftype', group_id, got_val, error)
 call read_from_hdf5(hdr%iexternalforce, 'iexternalforce', group_id, got_val, error)
 call read_from_hdf5(hdr%ieos, 'ieos', group_id, got_val, error)
 call read_from_hdf5(hdr%time, 'time', group_id, got_val, error)
 call read_from_hdf5(hdr%dtmax, 'dtmax', group_id, got_val, error)
 call read_from_hdf5(hdr%gamma, 'gamma', group_id, got_val, error)
 call read_from_hdf5(hdr%rhozero, 'rhozero', group_id, got_val, error)
 call read_from_hdf5(rval, 'RK2', group_id, got_val, error)
 hdr%polyk = rval/1.5
 call read_from_hdf5(hdr%hfact, 'hfact', group_id, got_val, error)
 call read_from_hdf5(hdr%tolh, 'tolh', group_id, got_val, error)
 call read_from_hdf5(hdr%C_cour, 'C_cour', group_id, got_val, error)
 call read_from_hdf5(hdr%C_force, 'C_force', group_id, got_val, error)
 call read_from_hdf5(hdr%alpha, 'alpha', group_id, got_val, error)
 call read_from_hdf5(hdr%alphau, 'alphau', group_id, got_val, error)
 call read_from_hdf5(hdr%alphaB, 'alphaB', group_id, got_val, error)
 call read_from_hdf5(hdr%polyk2, 'polyk2', group_id, got_val, error)
 call read_from_hdf5(hdr%qfacdisc, 'qfacdisc', group_id, got_val, error)
 call read_from_hdf5(hdr%massoftype, 'massoftype', group_id, got_val, error)
 call read_from_hdf5(hdr%Bextx, 'Bextx', group_id, got_val, error)
 call read_from_hdf5(hdr%Bexty, 'Bexty', group_id, got_val, error)
 call read_from_hdf5(hdr%Bextz, 'Bextz', group_id, got_val, error)
 call read_from_hdf5(hdr%xmin, 'xmin', group_id, got_val, error)
 call read_from_hdf5(hdr%xmax, 'xmax', group_id, got_val, error)
 call read_from_hdf5(hdr%ymin, 'ymin', group_id, got_val, error)
 call read_from_hdf5(hdr%ymax, 'ymax', group_id, got_val, error)
 call read_from_hdf5(hdr%zmin, 'zmin', group_id, got_val, error)
 call read_from_hdf5(hdr%zmax, 'zmax', group_id, got_val, error)
 call read_from_hdf5(hdr%get_conserv, 'get_conserv', group_id, got_val, error)
 call read_from_hdf5(hdr%etot_in, 'etot_in', group_id, got_val, error)
 call read_from_hdf5(hdr%angtot_in, 'angtot_in', group_id, got_val, error)
 call read_from_hdf5(hdr%totmom_in, 'totmom_in', group_id, got_val, error)
 call read_from_hdf5(hdr%mdust_in, 'mdust_in', group_id, got_val, error)
 call read_from_hdf5(hdr%grainsize, 'grainsize', group_id, got_val, error)
 call read_from_hdf5(hdr%graindens, 'graindens', group_id, got_val, error)
 call read_from_hdf5(hdr%udist, 'udist', group_id, got_val, error)
 call read_from_hdf5(hdr%umass, 'umass', group_id, got_val, error)
 call read_from_hdf5(hdr%utime, 'utime', group_id, got_val, error)
 call read_from_hdf5(hdr%unit_Bfield, 'umagfd', group_id, got_val, error)

 ! Close the header group
 call close_hdf5group(group_id, error)

 ! Read external force information
 select case(hdr%iexternalforce)
 case(iext_gwinspiral_hdf5)

    ! Open externalforce group
    call open_hdf5group(file_id, 'header/externalforce', group_id, error)

    ! Read values
    call read_from_hdf5(extern%Nstar, 'Nstar', group_id, got_val, error)

    ! Close the externalforce group
    call close_hdf5group(group_id, error)

 case(iext_binary_hdf5, iext_corot_binary_hdf5)

    ! Create externalforce group
    call open_hdf5group(file_id, 'header/externalforce', group_id, error)

    ! Write values
    call read_from_hdf5(extern%xyzmh1, 'xyzmh1', group_id, got_val, error)
    call read_from_hdf5(extern%xyzmh2, 'xyzmh2', group_id, got_val, error)
    call read_from_hdf5(extern%vxyz1, 'vxyz1', group_id, got_val, error)
    call read_from_hdf5(extern%vxyz2, 'vxyz2', group_id, got_val, error)
    call read_from_hdf5(extern%accretedmass1, 'accretedmass1', group_id, got_val, error)
    call read_from_hdf5(extern%accretedmass2, 'accretedmass2', group_id, got_val, error)
    call read_from_hdf5(extern%a0, 'a0', group_id, got_val, error)
    call read_from_hdf5(extern%direction, 'direction', group_id, got_val, error)

    ! Close the externalforce group
    call close_hdf5group(group_id, error)

 end select

end subroutine read_hdf5_header

!--------------------------------------------------------------------
!+
!  read arrays for full dump
!+
!-------------------------------------------------------------------
subroutine read_hdf5_arrays( &
   file_id,                  &
   error,                    &
   npart,                    &
   nptmass,                  &
   iphase,                   &
   xyzh,                     &
   vxyzu,                    &
   xyzmh_ptmass,             &
   vxyz_ptmass,              &
   dt_in,                    &
   alphaind,                 &
   poten,                    &
   Bxyz,                     &
   Bevol,                    &
   dustfrac,                 &
   deltav,                   &
   dustprop,                 &
   tstop,                    &
   VrelVf,                   &
   dustgasprop,              &
   temperature,              &
   abundance,                &
   array_options,            &
   got_arrays)

 integer(HID_T),  intent(in)  :: file_id
 integer,         intent(in)  :: npart, nptmass
 type (arrays_options_hdf5), intent(in)  :: array_options
 type (got_arrays_hdf5),     intent(out) :: got_arrays
 integer(kind=1), intent(out) :: iphase(:)
 real,            intent(out) :: xyzh(:,:),         &
                                 vxyzu(:,:),        &
                                 xyzmh_ptmass(:,:), &
                                 vxyz_ptmass(:,:),  &
                                 Bxyz(:,:),         &
                                 Bevol(:,:),        &
                                 dustfrac(:,:),     &
                                 deltav(:,:,:),     &
                                 dustprop(:,:),     &
                                 dustgasprop(:,:),  &
                                 tstop(:,:),        &
                                 VrelVf(:),         &
                                 temperature(:),    &
                                 abundance(:,:)
 real(kind=4),    intent(out) :: dt_in(:),          &
                                 alphaind(:,:),     &
                                 poten(:)
 integer,         intent(out) :: error

 integer(HID_T) :: group_id
 logical :: got

 real(kind=4) :: rtmp(npart)

 error = 0

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
 call open_hdf5group(file_id, 'particles', group_id, error)

 ! Main arrays
 call read_from_hdf5(iphase, 'itype', group_id, got_arrays%got_iphase, error)
 call read_from_hdf5(xyzh(1:3,:), 'xyz', group_id, got, error)
 if (got) got_arrays%got_xyzh = .true.
 call read_from_hdf5(rtmp, 'h', group_id, got, error)
 if (got) then
    xyzh(4,:) = real(rtmp)
 else
    got_arrays%got_xyzh = .false.
 endif
 call read_from_hdf5(vxyzu(1:3,:), 'vxyz', group_id, got, error)
 if (got) got_arrays%got_vxyzu = .true.
 if (.not.array_options%isothermal) then
    call read_from_hdf5(vxyzu(4,:), 'u', group_id, got, error)
    if (.not.got) got_arrays%got_vxyzu = .false.
 endif
 if (array_options%ind_timesteps) call read_from_hdf5(dt_in, 'dt', group_id, got_arrays%got_dt_in, error)
 if (.not. array_options%const_av) call read_from_hdf5(alphaind(1,:), 'alpha', group_id, got_arrays%got_alpha, error)
 if (array_options%gravity) call read_from_hdf5(poten, 'poten', group_id, got_arrays%got_poten, error)

 ! MHD arrays
 if (array_options%mhd) then
    call read_from_hdf5(Bxyz, 'Bxyz', group_id, got_arrays%got_Bxyz, error)
    call read_from_hdf5(Bevol(4,:), 'psi', group_id, got_arrays%got_psi, error)
 endif

 ! Dust arrays
 if (array_options%use_dust) then
    call read_from_hdf5(dustfrac, 'dustfrac', group_id, got_arrays%got_dustfrac, error)
    call read_from_hdf5(tstop, 'tstop', group_id, got_arrays%got_tstop, error)
 endif
 if (array_options%use_dustfrac) call read_from_hdf5(deltav, 'deltavxyz', group_id, got_arrays%got_deltav, error)
 if (array_options%use_dustgrowth) then
    call read_from_hdf5(dustprop(1,:), 'grainsize', group_id, got_arrays%got_dustprop(1), error)
    call read_from_hdf5(dustprop(2,:), 'graindens', group_id, got_arrays%got_dustprop(2), error)
    call read_from_hdf5(VrelVf(:), 'vrel_on_vfrag', group_id, got_arrays%got_VrelVf, error)
    call read_from_hdf5(dustgasprop(3,:), 'St', group_id, got_arrays%got_St, error)
 endif

 ! Other Arrays
 if (array_options%h2chemistry) call read_from_hdf5(abundance, 'abundance', group_id, got_arrays%got_abund, error)
 if (array_options%store_temperature) call read_from_hdf5(temperature, 'T', group_id, got_arrays%got_temp, error)

 ! Close the particles group
 call close_hdf5group(group_id, error)

 ! Open sinks group
 call open_hdf5group(file_id, 'sinks', group_id, error)

 ! Sink arrays
 if (nptmass > 0) then
    call read_from_hdf5(xyzmh_ptmass(1:3,1:nptmass), 'xyz', group_id, got_arrays%got_sink_data(1), error)
    got_arrays%got_sink_data(1:3) = got_arrays%got_sink_data(1)
    call read_from_hdf5(xyzmh_ptmass(4,1:nptmass), 'm', group_id, got_arrays%got_sink_data(4), error)
    call read_from_hdf5(xyzmh_ptmass(5,1:nptmass), 'h', group_id, got_arrays%got_sink_data(5), error)
    call read_from_hdf5(xyzmh_ptmass(6,1:nptmass), 'hsoft', group_id, got_arrays%got_sink_data(6), error)
    call read_from_hdf5(xyzmh_ptmass(7,1:nptmass), 'maccreted', group_id, got_arrays%got_sink_data(7), error)
    call read_from_hdf5(xyzmh_ptmass(8:10,1:nptmass), 'spinxyz', group_id, got_arrays%got_sink_data(8), error)
    got_arrays%got_sink_data(8:10) = got_arrays%got_sink_data(8)
    call read_from_hdf5(xyzmh_ptmass(11,1:nptmass), 'tlast', group_id, got_arrays%got_sink_data(11), error)
    call read_from_hdf5(vxyz_ptmass(:,1:nptmass), 'vxyz', group_id, got_arrays%got_sink_vels, error)
 endif

 ! Close the sinks group
 call close_hdf5group(group_id, error)

end subroutine read_hdf5_arrays

end module utils_dumpfiles_hdf5
