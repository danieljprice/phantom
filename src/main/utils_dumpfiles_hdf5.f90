!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
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
! :Dependencies: dim, eos, part, utils_hdf5
!
 use dim,        only:maxtypes,maxdustsmall,maxdustlarge,nabundances,nsinkproperties
 use part,       only:eos_vars_label,igasP,itemp,iX,iZ,imu,maxirad,n_nucleation
 use eos,        only:ieos,eos_is_non_ideal,eos_outputs_gasP
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

 integer, parameter :: maxdusttypes = maxdustsmall + maxdustlarge

 ! Ideally this section should come from Phantom modules.
 ! However, this module aims to be a library with minimal non-HDF5 dependencies.
 integer, parameter :: max_krome_nmols_hdf5 = 100
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
                          npartoftypetot(maxtypes),      &
                          iexternalforce,                &
                          ieos,                          &
                          idumpfile,                     &
                          idtmax_n_next,                 &
                          idtmax_frac_next
    real               :: time,                          &
                          dtmax,                         &
                          dtmax_user,                    &
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
                          massoftype(maxtypes),          &
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
                          mdust_in(maxdusttypes),        &
                          grainsize(maxdusttypes),       &
                          graindens(maxdusttypes),       &
                          udist,                         &
                          umass,                         &
                          utime,                         &
                          unit_Bfield
 end type header_hdf5

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
 end type externalforce_hdf5

 type got_arrays_hdf5
    logical :: got_iphase,                           &
               got_xyzh(4),                          &
               got_vxyzu(4),                         &
               got_dustfrac(maxdusttypes),           &
               got_deltav,                           &
               got_abund(nabundances),               &
               got_dt_in,                            &
               got_alpha,                            &
               got_poten,                            &
               got_sink_data(nsinkproperties),       &
               got_sink_vels(3),                     &
               got_Bxyz(3),                          &
               got_psi,                              &
               got_temp,                             &
               got_dustgasprop(3),                   &
               got_dustprop(2),                      &
               got_St,                               &
               got_VrelVf,                           &
               got_pxyzu(4),                         &
               got_raden(maxirad),                   &
               got_kappa,                            &
               got_Tdust,                            &
               got_iorig,                            &
               got_krome_mols(max_krome_nmols_hdf5), &
               got_krome_gamma,                      &
               got_krome_mu,                         &
               got_krome_T,                          &
               got_x,                                &
               got_z,                                &
               got_mu,                               &
               got_nucleation(n_nucleation),         &
               got_orig
 end type got_arrays_hdf5

 type arrays_options_hdf5
    logical :: isothermal,             &
               const_av,               &
               ind_timesteps,          &
               gravity,                &
               mhd,                    &
               mhd_nonideal,           &
               use_dust,               &
               use_dustfrac,           &
               use_dustgrowth,         &
               h2chemistry,            &
               lightcurve,             &
               prdrag,                 &
               store_dust_temperature, &
               radiation,              &
               krome,                  &
               nucleation,             &
               gr
    integer :: ieos,              &
               ndivcurlB,         &
               ndivcurlv,         &
               ndustsmall,        &
               ndustlarge
 end type arrays_options_hdf5

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
 call write_to_hdf5(hdr%npartoftypetot(1:maxtypes), 'npartoftype', group_id, error)
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
   eos_vars,                  &
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
   iorig,                     &
   dustfrac,                  &
   tstop,                     &
   deltav,                    &
   dustprop,                  &
   VrelVf,                    &
   dustgasprop,               &
   abundance,                 &
   divcurlv,                  &
   luminosity,                &
   beta_pr,                   &
   pxyzu,                     &
   dens,                      &
   T_gas_cool,                &
   nucleation,                &
   dust_temp,                 &
   rad,                       &
   radprop,                   &
   array_options              &
)

 integer(HID_T),  intent(in) :: file_id
 integer,         intent(out):: error
 integer,         intent(in) :: npart, nptmass
 real,            intent(in) :: dtind(:),          &
                                beta_pr(:),        &
                                VrelVf(:),         &
                                xyzh(:,:),         &
                                vxyzu(:,:),        &
                                Bxyz(:,:),         &
                                Bevol(:,:),        &
                                eos_vars(:,:),     &
                                eta_nimhd(:,:),    &
                                xyzmh_ptmass(:,:), &
                                vxyz_ptmass(:,:),  &
                                dustfrac(:,:),     &
                                tstop(:,:),        &
                                dustprop(:,:),     &
                                dustgasprop(:,:),  &
                                abundance(:,:),    &
                                deltav(:,:,:),     &
                                pxyzu(:,:),        &
                                dens(:),           &
                                T_gas_cool(:),     &
                                nucleation(:,:),   &
                                dust_temp(:),      &
                                rad(:,:),          &
                                radprop(:,:)
 real(kind=4),    intent(in) :: poten(:),          &
                                divBsymm(:),       &
                                luminosity(:),     &
                                alphaind(:,:),     &
                                divcurlv(:,:),     &
                                divcurlB(:,:)
 integer(kind=1), intent(in) :: iphase(:)
 integer(kind=8), intent(in) :: iorig(:)
 type (arrays_options_hdf5), intent(in) :: array_options

 integer(HID_T) :: group_id
 integer :: ndusttypes,ndustsmall,ieos

 error = 0
 ieos = array_options%ieos

 ! Create particles group
 call create_hdf5group(file_id, 'particles', group_id, error)

 ! Type, position, smoothing length, velocity
 call write_to_hdf5(iphase(1:npart), 'itype', group_id, error)
 call write_to_hdf5(xyzh(1:3,1:npart), 'xyz', group_id, error)
 call write_to_hdf5(real(xyzh(4,1:npart), kind=4), 'h', group_id, error)
 call write_to_hdf5(vxyzu(1:3,1:npart), 'vxyz', group_id, error)

 ! Equation of state
 if (.not.array_options%isothermal) then
    call write_to_hdf5(vxyzu(4,1:npart), 'u', group_id, error)
 endif
 if (eos_outputs_gasP(ieos)) then
    call write_to_hdf5(eos_vars(igasP,1:npart), eos_vars_label(igasP), group_id, error)
 endif
 if (eos_is_non_ideal(ieos)) then
    call write_to_hdf5(eos_vars(itemp,1:npart), eos_vars_label(itemp), group_id, error)
 endif
 if (array_options%lightcurve) then
    call write_to_hdf5(luminosity(1:npart), 'luminosity', group_id, error)
 endif
 if (array_options%prdrag) then
    call write_to_hdf5(real(beta_pr(1:npart), kind=4), 'beta_pr', group_id, error)
 endif

 ! General relativity
 if (array_options%gr) then
    call write_to_hdf5(pxyzu(1:3,1:npart), 'gr_momentum', group_id, error)
    call write_to_hdf5(pxyzu(4,1:npart), 'gr_entropy', group_id, error)
    call write_to_hdf5(dens(1:npart), 'gr_density', group_id, error)
 endif

 ! Viscosity (only ever write 'first' alpha)
 if (.not.array_options%const_av) then
    call write_to_hdf5(alphaind(1,1:npart), 'alpha', group_id, error)
 endif

 ! Individual timesteps
 if (array_options%ind_timesteps) then
    call write_to_hdf5(real(dtind(1:npart), kind=4), 'dt', group_id, error)
 endif

 ! Self-gravity
 if (array_options%gravity) then
    call write_to_hdf5(poten(1:npart), 'poten', group_id, error)
 endif

 ! MHD
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

 ! Dust
 ndustsmall = array_options%ndustsmall
 ndusttypes = array_options%ndustsmall + array_options%ndustlarge
 if (array_options%use_dust .and. ndusttypes > 0) then
    call write_to_hdf5(dustfrac(1:ndusttypes,1:npart), 'dustfrac', group_id, error)
    call write_to_hdf5(tstop(1:ndusttypes,1:npart), 'tstop', group_id, error)
    if (array_options%use_dustfrac .and. ndustsmall > 0) then
       call write_to_hdf5(deltav(:,1:ndustsmall,1:npart), 'deltavxyz', group_id, error)
    endif
 endif

 ! Dust growth
 if (array_options%use_dustgrowth) then
    call write_to_hdf5(dustprop(1,1:npart), 'grainsize', group_id, error)
    call write_to_hdf5(dustprop(2,1:npart), 'graindens', group_id, error)
    call write_to_hdf5(VrelVf(1:npart), 'vrel_on_vfrag', group_id, error)
    call write_to_hdf5(dustgasprop(3,1:npart), 'St', group_id, error)
 endif

 ! Chemistry
 if (array_options%h2chemistry) then
    call write_to_hdf5(abundance(:,1:npart), 'abundance', group_id, error)
 endif

 ! Chemistry (Krome)
 if (array_options%krome) then
    call write_to_hdf5(abundance(:,1:npart), 'abundance', group_id, error)
    call write_to_hdf5(T_gas_cool(1:npart), 'T_gas_cool', group_id, error)
 endif

 ! Nucleation
 if (array_options%nucleation) then
    call write_to_hdf5(nucleation(1,1:npart), 'nucleation_Jstar', group_id, error)
    call write_to_hdf5(nucleation(2,1:npart), 'nucleation_K0', group_id, error)
    call write_to_hdf5(nucleation(3,1:npart), 'nucleation_K1', group_id, error)
    call write_to_hdf5(nucleation(4,1:npart), 'nucleation_K2', group_id, error)
    call write_to_hdf5(nucleation(5,1:npart), 'nucleation_K3', group_id, error)
    call write_to_hdf5(nucleation(6,1:npart), 'nucleation_mu', group_id, error)
    call write_to_hdf5(nucleation(7,1:npart), 'nucleation_gamma', group_id, error)
    call write_to_hdf5(nucleation(8,1:npart), 'nucleation_S'    , group_id, error)
    call write_to_hdf5(nucleation(9,1:npart), 'nucleation_kappa', group_id, error)
 endif

 ! Radiation
 if (array_options%store_dust_temperature) then
    call write_to_hdf5(dust_temp(1:npart), 'temperature_dust', group_id, error)
 endif
 if (array_options%radiation) then
    call write_to_hdf5(rad(1,1:npart), 'radiation_xi', group_id, error)
    call write_to_hdf5(radprop(1:3,1:npart), 'radition_F', group_id, error)
    call write_to_hdf5(radprop(4,1:npart), 'radiation_kappa', group_id, error)
    call write_to_hdf5(radprop(5,1:npart), 'radiation_thick', group_id, error)
    call write_to_hdf5(radprop(6,1:npart), 'radiation_numph', group_id, error)
    call write_to_hdf5(radprop(7,1:npart), 'radiation_vorcl', group_id, error)
 endif

 ! Divergence and curl of velocity
 if (array_options%ndivcurlv >= 1) then
    call write_to_hdf5(divcurlv(1,1:npart), 'divv', group_id, error)
    if (array_options%ndivcurlv>=4) then
       call write_to_hdf5(divcurlv(2:4,1:npart), 'curlvxyz', group_id, error)
    endif
 endif

 ! Particle IDs
 call write_to_hdf5(iorig(1:npart), 'iorig', group_id, error)

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
    call write_to_hdf5(xyzmh_ptmass(15,1:nptmass),'mdotloss',group_id,error)
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
   abundance,                       &
   luminosity,                      &
   rad,                             &
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
                                abundance(:,:),    &
                                rad(:,:)
 real(kind=4),    intent(in) :: luminosity(:)
 integer(kind=1), intent(in) :: iphase(:)
 type (arrays_options_hdf5), intent(in)  :: array_options

 integer(HID_T) :: group_id
 integer :: ndusttypes

 error = 0

 ! Create particles group
 call create_hdf5group(file_id, 'particles', group_id, error)

 ! Type, position, smoothing length
 call write_to_hdf5(iphase(1:npart), 'itype', group_id, error)
 call write_to_hdf5(real(xyzh(1:3,1:npart), kind=4), 'xyz', group_id, error)
 call write_to_hdf5(real(xyzh(4,1:npart), kind=4), 'h', group_id, error)

 ! Equation of state
 if (array_options%lightcurve) then
    call write_to_hdf5(luminosity(1:npart), 'luminosity', group_id, error)
 endif

 ! MHD
 if (array_options%mhd) then
    call write_to_hdf5(real(Bxyz(:,1:npart), kind=4), 'Bxyz', group_id, error)
 endif

 ! Dust
 ndusttypes = array_options%ndustsmall + array_options%ndustlarge
 if (array_options%use_dust .and. ndusttypes > 0) then
    call write_to_hdf5(real(dustfrac(1:ndusttypes,1:npart), kind=4), 'dustfrac', group_id, error)
 endif

 ! Dustgrowth
 if (array_options%use_dustgrowth) then
    call write_to_hdf5(real(dustprop(1,1:npart), kind=4), 'grainsize', group_id, error)
    call write_to_hdf5(real(dustprop(2,1:npart), kind=4), 'graindens', group_id, error)
 endif

 ! Chemistry
 if (array_options%h2chemistry) then
    call write_to_hdf5(real(abundance(:,1:npart), kind=4), 'abundance', group_id, error)
 endif

 ! Radiation
 if (array_options%radiation) then
    call write_to_hdf5(real(rad(1,1:npart), kind=4), 'radiation_xi', group_id, error)
 endif

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
    call write_to_hdf5(real(xyzmh_ptmass(12,1:nptmass), kind=4), 'lum', group_id,error)
    call write_to_hdf5(real(xyzmh_ptmass(13,1:nptmass), kind=4), 'Teff', group_id,error)
    call write_to_hdf5(real(xyzmh_ptmass(14,1:nptmass), kind=4), 'Reff', group_id,error)
    call write_to_hdf5(real(xyzmh_ptmass(15,1:nptmass), kind=4), 'mdotloss', group_id,error)
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
   ndusttypes,               &
   iphase,                   &
   xyzh,                     &
   vxyzu,                    &
   xyzmh_ptmass,             &
   vxyz_ptmass,              &
   eos_vars,                 &
   dt_in,                    &
   alphaind,                 &
   poten,                    &
   Bxyz,                     &
   Bevol,                    &
   iorig,                    &
   dustfrac,                 &
   deltav,                   &
   dustprop,                 &
   VrelVf,                   &
   dustgasprop,              &
   abundance,                &
   pxyzu,                    &
   T_gas_cool,               &
   nucleation,               &
   dust_temp,                &
   rad,                      &
   radprop,                  &
   array_options,            &
   got_arrays)

 integer(HID_T),  intent(in)  :: file_id
 integer,         intent(in)  :: npart, nptmass, ndusttypes
 type (arrays_options_hdf5), intent(in)  :: array_options
 type (got_arrays_hdf5),     intent(out) :: got_arrays
 integer(kind=1), intent(out) :: iphase(:)
 integer(kind=8), intent(out) :: iorig(:)
 real,            intent(out) :: xyzh(:,:),         &
                                 vxyzu(:,:),        &
                                 xyzmh_ptmass(:,:), &
                                 vxyz_ptmass(:,:),  &
                                 eos_vars(:,:),     &
                                 Bxyz(:,:),         &
                                 Bevol(:,:),        &
                                 dustfrac(:,:),     &
                                 deltav(:,:,:),     &
                                 dustprop(:,:),     &
                                 dustgasprop(:,:),  &
                                 VrelVf(:),         &
                                 abundance(:,:),    &
                                 pxyzu(:,:),        &
                                 T_gas_cool(:),     &
                                 nucleation(:,:),   &
                                 dust_temp(:),      &
                                 rad(:,:),          &
                                 radprop(:,:)
 real(kind=4),    intent(out) :: dt_in(:),          &
                                 alphaind(:,:),     &
                                 poten(:)
 integer,         intent(out) :: error

 integer(HID_T) :: group_id
 logical :: got

 real(kind=4) :: rtmp(npart)
 real(kind=4) :: r2tmp(ndusttypes,npart)

 error = 0

 got_arrays%got_iphase      = .false.
 got_arrays%got_xyzh        = .false.
 got_arrays%got_vxyzu       = .false.
 got_arrays%got_dustfrac    = .false.
 got_arrays%got_deltav      = .false.
 got_arrays%got_abund       = .false.
 got_arrays%got_dt_in       = .false.
 got_arrays%got_alpha       = .false.
 got_arrays%got_poten       = .false.
 got_arrays%got_sink_data   = .false.
 got_arrays%got_sink_vels   = .false.
 got_arrays%got_Bxyz        = .false.
 got_arrays%got_psi         = .false.
 got_arrays%got_temp        = .false.
 got_arrays%got_dustprop    = .false.
 got_arrays%got_St          = .false.
 got_arrays%got_VrelVf      = .false.
 got_arrays%got_pxyzu       = .false.
 got_arrays%got_raden       = .false.
 got_arrays%got_kappa       = .false.
 got_arrays%got_Tdust       = .false.
 got_arrays%got_iorig       = .false.
 got_arrays%got_krome_mols  = .false.
 got_arrays%got_krome_gamma = .false.
 got_arrays%got_krome_mu    = .false.
 got_arrays%got_krome_T     = .false.
 got_arrays%got_x           = .false.
 got_arrays%got_z           = .false.
 got_arrays%got_mu          = .false.
 got_arrays%got_nucleation  = .false.

 ! Open particles group
 call open_hdf5group(file_id, 'particles', group_id, error)

 ! Type, position, smoothing length, velocity
 call read_from_hdf5(iphase, 'itype', group_id, got_arrays%got_iphase, error)
 call read_from_hdf5(xyzh(1:3,:), 'xyz', group_id, got, error)
 if (got) got_arrays%got_xyzh = .true.
 call read_from_hdf5(rtmp(1:npart), 'h', group_id, got, error)
 if (got) then
    xyzh(4,1:npart) = real(rtmp(1:npart))
 else
    got_arrays%got_xyzh = .false.
 endif
 call read_from_hdf5(vxyzu(1:3,:), 'vxyz', group_id, got, error)
 if (got) got_arrays%got_vxyzu = .true.

 ! Equation of state
 if (.not.array_options%isothermal) then
    call read_from_hdf5(vxyzu(4,:), 'u', group_id, got, error)
    if (.not.got) got_arrays%got_vxyzu = .false.
 endif
 if (eos_is_non_ideal(ieos)) then
    call read_from_hdf5(eos_vars(itemp,:), eos_vars_label(itemp), group_id, got_arrays%got_temp, error)
 endif

 call read_from_hdf5(eos_vars(iX,:), eos_vars_label(iX), group_id, got_arrays%got_x, error)
 call read_from_hdf5(eos_vars(iZ,:), eos_vars_label(iZ), group_id, got_arrays%got_z, error)
 call read_from_hdf5(eos_vars(imu,:), eos_vars_label(imu), group_id, got_arrays%got_mu, error)

 ! General relativity
 if (array_options%gr) then
    call read_from_hdf5(pxyzu(1:3,:), 'gr_momentum', group_id, got, error)
    call read_from_hdf5(pxyzu(4,:), 'gr_entropy', group_id, got, error)
    if (got) got_arrays%got_pxyzu = .true.
 endif

 ! Viscosity (only ever write 'first' alpha)
 if (.not. array_options%const_av) then
    call read_from_hdf5(alphaind(1,:), 'alpha', group_id, got_arrays%got_alpha, error)
 endif

 ! Individual timesteps
 if (array_options%ind_timesteps) then
    call read_from_hdf5(dt_in, 'dt', group_id, got_arrays%got_dt_in, error)
 endif

 ! Self-gravity
 if (array_options%gravity) then
    call read_from_hdf5(poten, 'poten', group_id, got_arrays%got_poten, error)
 endif

 ! MHD
 if (array_options%mhd) then
    call read_from_hdf5(Bxyz, 'Bxyz', group_id, got, error)
    if (got) got_arrays%got_Bxyz = .true.
    call read_from_hdf5(Bevol(4,:), 'psi', group_id, got_arrays%got_psi, error)
 endif

 ! Dust
 if (array_options%use_dust) then
    call read_from_hdf5(r2tmp, 'dustfrac', group_id, got, error)
    if (got) got_arrays%got_dustfrac = .true.
    dustfrac(1:ndusttypes,1:npart) = real(r2tmp(1:ndusttypes,1:npart))
 endif
 if (array_options%use_dustfrac) call read_from_hdf5(deltav, 'deltavxyz', group_id, got_arrays%got_deltav, error)

 ! Dust growth
 if (array_options%use_dustgrowth) then
    call read_from_hdf5(dustprop(1,:), 'grainsize', group_id, got_arrays%got_dustprop(1), error)
    call read_from_hdf5(dustprop(2,:), 'graindens', group_id, got_arrays%got_dustprop(2), error)
    call read_from_hdf5(VrelVf(:), 'vrel_on_vfrag', group_id, got_arrays%got_VrelVf, error)
    call read_from_hdf5(dustgasprop(3,:), 'St', group_id, got_arrays%got_St, error)
 endif

 ! Chemistry
 if (array_options%h2chemistry) then
    call read_from_hdf5(abundance, 'abundance', group_id, got, error)
    if (got) got_arrays%got_abund = .true.
 endif

 ! Chemistry (Krome)
 if (array_options%krome) then
    call read_from_hdf5(abundance, 'abundance', group_id, got, error)
    if (got) got_arrays%got_krome_mols = .true.
    call read_from_hdf5(T_gas_cool, 'T_gas_cool', group_id, got_arrays%got_krome_gamma, error)
 endif

 ! Nucleation
 if (array_options%nucleation) then
    call read_from_hdf5(nucleation(1,1:npart), 'nucleation_Jstar', group_id, got, error)
    call read_from_hdf5(nucleation(2,1:npart), 'nucleation_K0', group_id, got, error)
    call read_from_hdf5(nucleation(3,1:npart), 'nucleation_K1', group_id, got, error)
    call read_from_hdf5(nucleation(4,1:npart), 'nucleation_K2', group_id, got, error)
    call read_from_hdf5(nucleation(5,1:npart), 'nucleation_K3', group_id, got, error)
    call read_from_hdf5(nucleation(6,1:npart), 'nucleation_mu', group_id, got, error)
    call read_from_hdf5(nucleation(7,1:npart), 'nucleation_gamma', group_id, got, error)
    call read_from_hdf5(nucleation(8,1:npart), 'nucleation_S', group_id, got, error)
    call read_from_hdf5(nucleation(9,1:npart), 'nucleation_kappa', group_id, got, error)
    if (got) got_arrays%got_nucleation = .true.
 endif

 ! Radiation
 if (array_options%store_dust_temperature) then
    call read_from_hdf5(dust_temp(1:npart), 'temperature_dust', group_id, got, error)
 endif
 if (array_options%radiation) then
    call read_from_hdf5(rad(1,1:npart), 'radiation_xi', group_id, got, error)
    call read_from_hdf5(radprop(1:3,1:npart), 'radition_F', group_id, got, error)
    call read_from_hdf5(radprop(4,1:npart), 'radiation_kappa', group_id, got, error)
    call read_from_hdf5(radprop(5,1:npart), 'radiation_thick', group_id, got, error)
    call read_from_hdf5(radprop(6,1:npart), 'radiation_numph', group_id, got, error)
    call read_from_hdf5(radprop(7,1:npart), 'radiation_vorcl', group_id, got, error)
 endif

 call read_from_hdf5(iorig, 'iorig', group_id, got, error)
 if (got) got_arrays%got_orig = .true.

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
    call read_from_hdf5(xyzmh_ptmass(12,1:nptmass), 'lum', group_id, got_arrays%got_sink_data(12), error)
    call read_from_hdf5(xyzmh_ptmass(13,1:nptmass), 'Teff', group_id, got_arrays%got_sink_data(13), error)
    call read_from_hdf5(xyzmh_ptmass(14,1:nptmass), 'Reff', group_id, got_arrays%got_sink_data(14), error)
    call read_from_hdf5(xyzmh_ptmass(15,1:nptmass), 'mdotloss', group_id, got_arrays%got_sink_data(15), error)
    call read_from_hdf5(vxyz_ptmass(:,1:nptmass), 'vxyz', group_id, got, error)
    if (got) got_arrays%got_sink_vels = .true.
 endif

 ! Close the sinks group
 call close_hdf5group(group_id, error)

end subroutine read_hdf5_arrays

end module utils_dumpfiles_hdf5
