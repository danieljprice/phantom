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

 integer(HID_T), public :: hdf5_file_id

 private

contains


!--------------------------------------------------------------------
!+
!  write header
!+
!-------------------------------------------------------------------
subroutine write_hdf5_header(file_id,error,fileident,maxtypes,nblocks,isink,   &
                             nptmass,ndustlarge,ndustsmall,idust,              &
                             phantom_version_major,phantom_version_minor,      &
                             phantom_version_micro,nparttot,npartoftypetot,    &
                             iexternalforce,ieos,t,dtmax,gamma,rhozero,polyk,  &
                             hfact,tolh,C_cour,C_force,alpha,alphau,alphaB,    &
                             polyk2,qfacdisc,massoftype,Bextx,Bexty,Bextz,     &
                             xmin,xmax,ymin,ymax,zmin,zmax,get_conserv,        &
                             etot_in,angtot_in,totmom_in,mdust_in,grainsize,   &
                             graindens,udist,umass,utime,unit_Bfield)
 integer(HID_T),   intent(in) :: file_id
 integer,          intent(out):: error
 character(len=*), intent(in) :: fileident
 integer, intent(in) :: maxtypes,nblocks,isink,nptmass,ndustlarge,ndustsmall,  &
                        idust,phantom_version_major,phantom_version_minor,     &
                        phantom_version_micro,nparttot,npartoftypetot(:),      &
                        iexternalforce,ieos
 real,    intent(in) :: t,dtmax,gamma,rhozero,polyk,hfact,tolh,C_cour,C_force, &
                        alpha,alphau,alphaB,polyk2,qfacdisc,massoftype(:),     &
                        Bextx,Bexty,Bextz,xmin,xmax,ymin,ymax,zmin,zmax,       &
                        get_conserv,etot_in,angtot_in,totmom_in,mdust_in(:),   &
                        grainsize(:),graindens(:),udist,umass,utime,unit_Bfield
 integer(HID_T) :: group_id
 integer :: errors(53)

 errors(:) = 0

 ! Create header group
 call create_hdf5group(file_id,'header',group_id,errors(1))

 ! Write things to header group
 call write_to_hdf5(fileident,'fileident',group_id,errors(2))
 call write_to_hdf5(maxtypes,'ntypes',group_id,errors(3))
 call write_to_hdf5(nblocks,'nblocks',group_id,errors(4))
 call write_to_hdf5(isink,'isink',group_id,errors(5))
 call write_to_hdf5(nptmass,'nptmass',group_id,errors(6))
 call write_to_hdf5(ndustlarge,'ndustlarge',group_id,errors(7))
 call write_to_hdf5(ndustsmall,'ndustsmall',group_id,errors(8))
 call write_to_hdf5(idust,'idust',group_id,errors(9))
 call write_to_hdf5(phantom_version_major,'majorv',group_id,errors(10))
 call write_to_hdf5(phantom_version_minor,'minorv',group_id,errors(11))
 call write_to_hdf5(phantom_version_micro,'microv',group_id,errors(12))
 call write_to_hdf5(nparttot,'nparttot',group_id,errors(13))
 call write_to_hdf5(npartoftypetot(1:maxtypes),'npartoftype',group_id,errors(14))
 call write_to_hdf5(iexternalforce,'iexternalforce',group_id,errors(15))
 call write_to_hdf5(ieos,'ieos',group_id,errors(16))
 call write_to_hdf5(t,'time',group_id,errors(17))
 call write_to_hdf5(dtmax,'dtmax',group_id,errors(18))
 call write_to_hdf5(gamma,'gamma',group_id,errors(19))
 call write_to_hdf5(rhozero,'rhozero',group_id,errors(20))
 call write_to_hdf5(1.5*polyk,'RK2',group_id,errors(21))
 call write_to_hdf5(hfact,'hfact',group_id,errors(22))
 call write_to_hdf5(tolh,'tolh',group_id,errors(23))
 call write_to_hdf5(C_cour,'C_cour',group_id,errors(24))
 call write_to_hdf5(C_force,'C_force',group_id,errors(25))
 call write_to_hdf5(alpha,'alpha',group_id,errors(26))
 call write_to_hdf5(alphau,'alphau',group_id,errors(27))
 call write_to_hdf5(alphaB,'alphaB',group_id,errors(28))
 call write_to_hdf5(polyk2,'polyk2',group_id,errors(29))
 call write_to_hdf5(qfacdisc,'qfacdisc',group_id,errors(30))
 call write_to_hdf5(massoftype,'massoftype',group_id,errors(31))
 call write_to_hdf5(Bextx,'Bextx',group_id,errors(32))
 call write_to_hdf5(Bexty,'Bexty',group_id,errors(33))
 call write_to_hdf5(Bextz,'Bextz',group_id,errors(34))
 call write_to_hdf5(0.,'dum',group_id,errors(35))
 ! TODO: NEED TO FIND A WAY TO DO THIS
 ! if (iexternalforce /= 0) call write_headeropts_extern(iexternalforce,hdr,t,ierr)
 call write_to_hdf5(xmin,'xmin',group_id,errors(36))
 call write_to_hdf5(xmax,'xmax',group_id,errors(37))
 call write_to_hdf5(ymin,'ymin',group_id,errors(38))
 call write_to_hdf5(ymax,'ymax',group_id,errors(39))
 call write_to_hdf5(zmin,'zmin',group_id,errors(40))
 call write_to_hdf5(zmax,'zmax',group_id,errors(41))
 call write_to_hdf5(get_conserv,'get_conserv',group_id,errors(42))
 call write_to_hdf5(etot_in,'etot_in',group_id,errors(43))
 call write_to_hdf5(angtot_in,'angtot_in',group_id,errors(44))
 call write_to_hdf5(totmom_in,'totmom_in',group_id,errors(45))
 call write_to_hdf5(mdust_in,'mdust_in',group_id,errors(46))
 call write_to_hdf5(grainsize,'grainsize',group_id,errors(47))
 call write_to_hdf5(graindens,'graindens',group_id,errors(48))
 call write_to_hdf5(udist,'udist',group_id,errors(49))
 call write_to_hdf5(umass,'umass',group_id,errors(50))
 call write_to_hdf5(utime,'utime',group_id,errors(51))
 call write_to_hdf5(unit_Bfield,'umagfd',group_id,errors(52))

 ! Close the header group
 call close_hdf5group(group_id, errors(53))

 error = maxval(abs(errors))

end subroutine write_hdf5_header

!--------------------------------------------------------------------
!+
!  write arrays for full dump
!+
!-------------------------------------------------------------------
subroutine write_hdf5_arrays(file_id,error,npart,xyzh,vxyzu,iphase,pressure,  &
                             alphaind,dtind,poten,xyzmh_ptmass,vxyz_ptmass,   &
                             Bxyz,Bevol,divcurlB,divBsymm,eta_nimhd,dustfrac, &
                             tstop,deltav,dustprop,St,abundance,temperature,  &
                             divcurlv,luminosity,beta_pr,const_av,            &
                             ind_timesteps,gravity,nptmass,mhd,maxBevol,      &
                             ndivcurlB,mhd_nonideal,use_dust,use_dustfrac,    &
                             use_dustgrowth,h2chemistry,store_temperature,    &
                             ndivcurlv,lightcurve,prdrag,isothermal)

 integer(HID_T),  intent(in) :: file_id
 integer,         intent(out):: error
 integer,         intent(in) :: npart, nptmass,maxBevol,ndivcurlB,ndivcurlv
 real,            intent(in) :: pressure(:),dtind(:),beta_pr(:),St(:),         &
                                temperature(:),xyzh(:,:),vxyzu(:,:),Bxyz(:,:), &
                                Bevol(:,:), eta_nimhd(:,:),xyzmh_ptmass(:,:),  &
                                vxyz_ptmass(:,:),dustfrac(:,:),tstop(:,:),     &
                                dustprop(:,:),abundance(:,:),deltav(:,:,:)
 real(kind=4),    intent(in) :: poten(:),divBsymm(:),luminosity(:),            &
                                alphaind(:,:),divcurlv(:,:),divcurlB(:,:)
 integer(kind=1), intent(in) :: iphase(:)
 logical,         intent(in) :: const_av,ind_timesteps,gravity,mhd,            &
                                mhd_nonideal,use_dust,use_dustfrac,            &
                                use_dustgrowth,h2chemistry,store_temperature,  &
                                lightcurve,prdrag,isothermal

 integer(HID_T) :: group_id
 integer :: errors(44)

 errors(:) = 0

 ! Create particles group
 call create_hdf5group(file_id,'particles',group_id,errors(1))

 ! Main arrays
 call write_to_hdf5(xyzh(1:3,1:npart),'xyz',group_id,errors(2))
 ! Write smoothing length in single precision to save disc space
 call write_to_hdf5(real(xyzh(4,1:npart),kind=4),'h',group_id,errors(3))
 call write_to_hdf5(vxyzu(1:3,1:npart),'vxyz',group_id,errors(4))
 if (.not.isothermal) call write_to_hdf5(vxyzu(4,1:npart),'u',group_id,errors(5))
 call write_to_hdf5(iphase(1:npart),'itype',group_id,errors(6))
 call write_to_hdf5(pressure(1:npart),'pressure',group_id,errors(7))

 if (.not.const_av)  call write_to_hdf5(alphaind(1,1:npart),'alpha',group_id,errors(8))   ! Viscosity (only ever write 'first' alpha)
 if (ind_timesteps)  call write_to_hdf5(real(dtind(1:npart),kind=4),'dt',group_id,errors(9)) ! Individual timesteps
 if (gravity)        call write_to_hdf5(poten(1:npart),'poten',group_id,errors(10))

 ! MHD arrays
 if (mhd) then
    call write_to_hdf5(Bxyz(:,1:npart),'Bxyz',group_id,errors(11))
    if (maxBevol >= 4) then
       call write_to_hdf5(Bevol(4,1:npart),'psi',group_id,errors(12))
    endif
    if (ndivcurlB >= 1) then
       call write_to_hdf5(divcurlB(1,1:npart),'divB',group_id,errors(13))
       call write_to_hdf5(divcurlB(2:4,1:npart),'curlBxyz',group_id,errors(14))
    else
       call write_to_hdf5(divBsymm(1:npart),'divBsymm',group_id,errors(16))
    endif
    if (mhd_nonideal) then
       call write_to_hdf5(eta_nimhd(1,1:npart),'eta_{OR}',group_id,errors(17))
       call write_to_hdf5(eta_nimhd(2,1:npart),'eta_{HE}',group_id,errors(18))
       call write_to_hdf5(eta_nimhd(3,1:npart),'eta_{AD}',group_id,errors(19))
       call write_to_hdf5(eta_nimhd(4,1:npart),'ne/n'    ,group_id,errors(20))
    endif
 endif

 ! Dust arrays
 if (use_dust) then
    call write_to_hdf5(dustfrac(:,1:npart),'dustfrac',group_id,errors(21))
    call write_to_hdf5(tstop(:,1:npart),'tstop',group_id,errors(22))
 endif
 if (use_dustfrac) call write_to_hdf5(deltav(:,:,1:npart),'deltavxyz',group_id,errors(23))
 if (use_dustgrowth) then
    call write_to_hdf5(dustprop(1,1:npart),'grainsize',group_id,errors(24))
    call write_to_hdf5(dustprop(2,1:npart),'graindens',group_id,errors(25))
    call write_to_hdf5(dustprop(3,1:npart),'vrel/vfrag',group_id,errors(26))
    ! call write_to_hdf5(dustprop(4,:),'dv_dust',group_id,errors())
    call write_to_hdf5(St(1:npart),'St',group_id,errors(27))
 endif

 ! Other Arrays
 if (h2chemistry)       call write_to_hdf5(abundance(:,1:npart),'abundance',group_id,errors(28))
 if (store_temperature) call write_to_hdf5(temperature(1:npart),'T',group_id,errors(29))
 if (ndivcurlv >= 1) then
    call write_to_hdf5(divcurlv(1,1:npart),'divv',group_id,errors(30))
    if (ndivcurlv>=4) call write_to_hdf5(divcurlv(2:4,1:npart),'curlvxyz',group_id,errors(31))
 endif
 if (lightcurve) call write_to_hdf5(luminosity(1:npart),'luminosity',group_id,errors(32))
 if (prdrag)     call write_to_hdf5(real(beta_pr(1:npart),kind=4),'beta_pr',group_id,errors(33))

 ! Close the particles group
 call close_hdf5group(group_id, errors(34))

 ! Create sink group
 call create_hdf5group(file_id,'sinks',group_id,errors(35))
 if (nptmass > 0) then
    call write_to_hdf5(xyzmh_ptmass(1:3,1:nptmass),'xyz',group_id,errors, curent_error_id)
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
subroutine write_hdf5_arrays_small(file_id,error,npart,xyzh,iphase,           &
                                   xyzmh_ptmass,Bxyz,dustfrac,dustprop,St,    &
                                   abundance,luminosity,nptmass,mhd,use_dust, &
                                   use_dustgrowth,h2chemistry,lightcurve)

 integer(HID_T),  intent(in) :: file_id
 integer,         intent(out):: error
 integer,         intent(in) :: npart,nptmass
 real,            intent(in) :: St(:),xyzh(:,:),Bxyz(:,:),xyzmh_ptmass(:,:), &
                                dustfrac(:,:),dustprop(:,:),abundance(:,:)
 real(kind=4),    intent(in) :: luminosity(:)
 integer(kind=1), intent(in) :: iphase(:)
 logical,         intent(in) :: mhd,use_dust,use_dustgrowth,h2chemistry, &
                                lightcurve

 integer(HID_T) :: group_id
 integer :: errors(22)

 errors(:) = 0

 ! Create particles group
 call create_hdf5group(file_id,'particles',group_id,errors(1))

 ! Main arrays
 call write_to_hdf5(real(xyzh(1:3,1:npart),kind=4),'xyz',group_id,errors(2))
 call write_to_hdf5(real(xyzh(4,1:npart),kind=4),'h',group_id,errors(3))
 call write_to_hdf5(iphase(1:npart),'itype',group_id,errors(4))

 ! MHD arrays
 if (mhd) then
    call write_to_hdf5(real(Bxyz(:,1:npart),kind=4),'Bxyz',group_id,errors(5))
 endif

 ! Dust arrays
 if (use_dust) then
    call write_to_hdf5(real(dustfrac(:,1:npart),kind=4),'dustfrac',group_id,errors(6))
 endif
 if (use_dustgrowth) then
    call write_to_hdf5(real(dustprop(1,1:npart),kind=4),'grainsize',group_id,errors(7))
    call write_to_hdf5(real(dustprop(2,1:npart),kind=4),'graindens',group_id,errors(8))
    call write_to_hdf5(real(dustprop(3,1:npart),kind=4),'vrel/vfrag',group_id,errors(9))
    ! call write_to_hdf5(real(dustprop(4,1:npart),kind=4),'dv_dust',group_id,errors())
    call write_to_hdf5(real(St(1:npart),kind=4),'St',group_id,errors(10))
 endif

 ! Other Arrays
 if (h2chemistry) call write_to_hdf5(real(abundance(:,1:npart),kind=4),'abundance',group_id,errors(11))
 if (lightcurve)  call write_to_hdf5(luminosity(1:npart),'luminosity',group_id,errors(12))

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
subroutine read_hdf5_header(file_id,error,fileident,isink,nptmass,ndustlarge,  &
                            ndustsmall,npart,npartoftype,iexternalforce,ieos,  &
                            time,dtmax,gamma,rhozero,polyk,hfact,tolh,C_cour,  &
                            C_force,alpha,alphau,alphaB,polyk2,qfacdisc,       &
                            massoftype,Bextx,Bexty,Bextz,xmin,xmax,ymin,ymax,  &
                            zmin,zmax,get_conserv,etot_in,angtot_in,totmom_in, &
                            mdust_in,grainsize,graindens,udist,umass,utime,    &
                            unit_Bfield)

 integer(HID_T),   intent(in)  :: file_id
 character(len=*), intent(out) :: fileident

 integer, intent(out) :: error,          &
                         isink,          &
                         nptmass,        &
                         ndustlarge,     &
                         ndustsmall,     &
                         npart,          &
                         npartoftype(:), &
                         iexternalforce, &
                         ieos

 real, intent(out) :: time,          &
                      dtmax,         &
                      gamma,         &
                      rhozero,       &
                      polyk,         &
                      hfact,         &
                      tolh,          &
                      C_cour,        &
                      C_force,       &
                      alpha,         &
                      alphau,        &
                      alphaB,        &
                      polyk2,        &
                      qfacdisc,      &
                      massoftype(:), &
                      Bextx,         &
                      Bexty,         &
                      Bextz,         &
                      xmin,          &
                      xmax,          &
                      ymin,          &
                      ymax,          &
                      zmin,          &
                      zmax,          &
                      get_conserv,   &
                      etot_in,       &
                      angtot_in,     &
                      totmom_in,     &
                      mdust_in(:),   &
                      grainsize(:),  &
                      graindens(:),  &
                      udist,         &
                      umass,         &
                      utime,         &
                      unit_Bfield

 integer(HID_T) :: group_id
 integer        :: errors(53),ntypes
 real           :: rval
 logical        :: got_val

 errors(:) = 0

 ! Create header group
 call open_hdf5group(file_id,'header',group_id,errors(1))

 ! Write things to header group
 call read_from_hdf5(fileident,'fileident',group_id,got_val,errors(2))
 call read_from_hdf5(ntypes,'ntypes',group_id,got_val,errors(3))
 call read_from_hdf5(isink,'isink',group_id,got_val,errors(5))
 call read_from_hdf5(nptmass,'nptmass',group_id,got_val,errors(6))
 call read_from_hdf5(ndustlarge,'ndustlarge',group_id,got_val,errors(7))
 call read_from_hdf5(ndustsmall,'ndustsmall',group_id,got_val,errors(8))
 call read_from_hdf5(npart,'nparttot',group_id,got_val,errors(13))
 call read_from_hdf5(npartoftype(1:ntypes),'npartoftype',group_id,got_val,errors(14))
 call read_from_hdf5(iexternalforce,'iexternalforce',group_id,got_val,errors(15))
 call read_from_hdf5(ieos,'ieos',group_id,got_val,errors(16))
 call read_from_hdf5(time,'time',group_id,got_val,errors(17))
 call read_from_hdf5(dtmax,'dtmax',group_id,got_val,errors(18))
 call read_from_hdf5(gamma,'gamma',group_id,got_val,errors(19))
 call read_from_hdf5(rhozero,'rhozero',group_id,got_val,errors(20))
 call read_from_hdf5(rval,'RK2',group_id,got_val,errors(21))
 polyk = rval/1.5
 call read_from_hdf5(hfact,'hfact',group_id,got_val,errors(22))
 call read_from_hdf5(tolh,'tolh',group_id,got_val,errors(23))
 call read_from_hdf5(C_cour,'C_cour',group_id,got_val,errors(24))
 call read_from_hdf5(C_force,'C_force',group_id,got_val,errors(25))
 call read_from_hdf5(alpha,'alpha',group_id,got_val,errors(26))
 call read_from_hdf5(alphau,'alphau',group_id,got_val,errors(27))
 call read_from_hdf5(alphaB,'alphaB',group_id,got_val,errors(28))
 call read_from_hdf5(polyk2,'polyk2',group_id,got_val,errors(29))
 call read_from_hdf5(qfacdisc,'qfacdisc',group_id,got_val,errors(30))
 call read_from_hdf5(massoftype,'massoftype',group_id,got_val,errors(31))
 call read_from_hdf5(Bextx,'Bextx',group_id,got_val,errors(32))
 call read_from_hdf5(Bexty,'Bexty',group_id,got_val,errors(33))
 call read_from_hdf5(Bextz,'Bextz',group_id,got_val,errors(34))
 !call read_from_hdf5(0.,'dum',group_id,got_val,errors(35))
 ! TODO: NEED TO FIND A WAY TO DO THIS
 ! if (iexternalforce /= 0) call write_headeropts_extern(iexternalforce,hdr,t,ierr)
 call read_from_hdf5(xmin,'xmin',group_id,got_val,errors(36))
 call read_from_hdf5(xmax,'xmax',group_id,got_val,errors(37))
 call read_from_hdf5(ymin,'ymin',group_id,got_val,errors(38))
 call read_from_hdf5(ymax,'ymax',group_id,got_val,errors(39))
 call read_from_hdf5(zmin,'zmin',group_id,got_val,errors(40))
 call read_from_hdf5(zmax,'zmax',group_id,got_val,errors(41))
 call read_from_hdf5(get_conserv,'get_conserv',group_id,got_val,errors(42))
 call read_from_hdf5(etot_in,'etot_in',group_id,got_val,errors(43))
 call read_from_hdf5(angtot_in,'angtot_in',group_id,got_val,errors(44))
 call read_from_hdf5(totmom_in,'totmom_in',group_id,got_val,errors(45))
 call read_from_hdf5(mdust_in,'mdust_in',group_id,got_val,errors(46))
 call read_from_hdf5(grainsize,'grainsize',group_id,got_val,errors(47))
 call read_from_hdf5(graindens,'graindens',group_id,got_val,errors(48))
 call read_from_hdf5(udist,'udist',group_id,got_val,errors(49))
 call read_from_hdf5(umass,'umass',group_id,got_val,errors(50))
 call read_from_hdf5(utime,'utime',group_id,got_val,errors(51))
 call read_from_hdf5(unit_Bfield,'umagfd',group_id,got_val,errors(52))

 ! Close the header group
 call close_hdf5group(group_id,errors(53))

 error = maxval(abs(errors))

end subroutine read_hdf5_header

!--------------------------------------------------------------------
!+
!  read arrays for full dump
!+
!-------------------------------------------------------------------
subroutine read_hdf5_arrays(file_id,error,npart,nptmass,iphase,xyzh,vxyzu,   &
                            xyzmh_ptmass,vxyz_ptmass,dt_in,alphaind,poten,   &
                            Bxyz,Bevol,dustfrac,deltav,dustprop,tstop,St,    &
                            temperature,abundance,isothermal,const_av,       &
                            ind_timesteps,gravity,mhd,use_dust,use_dustfrac, &
                            use_dustgrowth,h2chemistry,store_temperature,    &
                            nsinkproperties,got_iphase,got_xyzh,got_vxyzu,   &
                            got_dustfrac,got_tstop,got_deltav,got_abund,     &
                            got_dt_in,got_alpha,got_poten,got_sink_data,     &
                            got_sink_vels,got_Bxyz,got_psi,got_temp,         &
                            got_dustprop,got_St)

 integer(HID_T),  intent(in)  :: file_id
 integer,         intent(in)  :: npart,nptmass,nsinkproperties
 logical,         intent(in)  :: isothermal,                     &
                                 const_av,                       &
                                 ind_timesteps,                  &
                                 gravity,                        &
                                 mhd,                            &
                                 use_dust,                       &
                                 use_dustfrac,                   &
                                 use_dustgrowth,                 &
                                 h2chemistry,                    &
                                 store_temperature
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
                                 tstop(:,:),                     &
                                 St(:),                          &
                                 temperature(:),                 &
                                 abundance(:,:)
 real(kind=4),    intent(out) :: dt_in(:),                       &
                                 alphaind(:,:),                  &
                                 poten(:)
 logical,         intent(out) :: got_iphase,                     &
                                 got_xyzh,                       &
                                 got_vxyzu,                      &
                                 got_dustfrac,                   &
                                 got_tstop,                      &
                                 got_deltav,                     &
                                 got_abund,                      &
                                 got_dt_in,                      &
                                 got_alpha,                      &
                                 got_poten,                      &
                                 got_sink_data(nsinkproperties), &
                                 got_sink_vels,                  &
                                 got_Bxyz,                       &
                                 got_psi,                        &
                                 got_temp,                       &
                                 got_dustprop(3),                &
                                 got_St
 integer,         intent(out) :: error

 integer(HID_T) :: group_id
 integer :: errors(44)
 logical :: got

 real(kind=4) :: rtmp(npart)

 errors(:) = 0

 got_iphase    = .false.
 got_xyzh      = .false.
 got_vxyzu     = .false.
 got_dustfrac  = .false.
 got_tstop     = .false.
 got_deltav    = .false.
 got_abund     = .false.
 got_dt_in     = .false.
 got_alpha     = .false.
 got_poten     = .false.
 got_sink_data = .false.
 got_sink_vels = .false.
 got_Bxyz      = .false.
 got_psi       = .false.
 got_temp      = .false.
 got_dustprop  = .false.
 got_St        = .false.

 ! Open particles group
 call open_hdf5group(file_id,'particles',group_id,errors(1))

 ! Main arrays
 call read_from_hdf5(iphase,'itype',group_id,got_iphase,errors(2))
 call read_from_hdf5(xyzh(1:3,:),'xyz',group_id,got,errors(3))
 if (got) got_xyzh = .true.
 call read_from_hdf5(rtmp,'h',group_id,got,errors(4))
 if (got) then
    xyzh(4,:) = real(rtmp)
 else
    got_xyzh = .false.
 endif
 call read_from_hdf5(vxyzu(1:3,:),'vxyz',group_id,got,errors(5))
 if (got) got_vxyzu = .true.
 if (.not.isothermal) then
    call read_from_hdf5(vxyzu(4,:),'u',group_id,got,errors(6))
    if (.not.got) got_vxyzu = .false.
 endif
 if (ind_timesteps) call read_from_hdf5(dt_in,'dt',group_id,got_dt_in,errors(7))
 if (.not.const_av) call read_from_hdf5(alphaind(1,:),'alpha',group_id,got_alpha,errors(8))
 if (gravity)       call read_from_hdf5(poten,'poten',group_id,got_poten,errors(9))

 ! MHD arrays
 if (mhd) then
    call read_from_hdf5(Bxyz,'Bxyz',group_id,got_Bxyz,errors(10))
    call read_from_hdf5(Bevol(4,:),'psi',group_id,got_psi,errors(11))
 endif

 ! Dust arrays
 if (use_dust) then
    call read_from_hdf5(dustfrac,'dustfrac',group_id,got_dustfrac,errors(12))
    call read_from_hdf5(tstop,'tstop',group_id,got_tstop,errors(13))
 endif
 if (use_dustfrac) call read_from_hdf5(deltav,'deltavxyz',group_id,got_deltav,errors(14))
 if (use_dustgrowth) then
    call read_from_hdf5(dustprop(1,:),'grainsize',group_id,got_dustprop(1),errors(15))
    call read_from_hdf5(dustprop(2,:),'graindens',group_id,got_dustprop(2),errors(16))
    call read_from_hdf5(dustprop(3,:),'vrel/vfrag',group_id,got_dustprop(3),errors(17))
    ! call read_from_hdf5(dustprop(4,:),'dv_dust',group_id,got_dv_dust,errors())
    call read_from_hdf5(St,'St',group_id,got_St,errors(18))
 endif

 ! Other Arrays
 if (h2chemistry) call read_from_hdf5(abundance,'abundance',group_id,got_abund,errors(19))
 if (store_temperature) call read_from_hdf5(temperature,'T',group_id,got_temp,errors(20))

 ! Close the particles group
 call close_hdf5group(group_id, errors(21))

 ! Open sinks group
 call open_hdf5group(file_id,'sinks',group_id,errors(22))

 ! Sink arrays
 if (nptmass > 0) then
    call read_from_hdf5(xyzmh_ptmass(1:3,1:nptmass),'xyz',group_id,got_sink_data(1),errors(23))
    got_sink_data(1:3) = got_sink_data(1)
    call read_from_hdf5(xyzmh_ptmass(4,1:nptmass),'m',group_id,got_sink_data(4),errors(24))
    call read_from_hdf5(xyzmh_ptmass(5,1:nptmass),'h',group_id,got_sink_data(5),errors(25))
    call read_from_hdf5(xyzmh_ptmass(6,1:nptmass),'hsoft',group_id,got_sink_data(6),errors(26))
    call read_from_hdf5(xyzmh_ptmass(7,1:nptmass),'maccreted',group_id,got_sink_data(7),errors(27))
    call read_from_hdf5(xyzmh_ptmass(8:10,1:nptmass),'spinxyz',group_id,got_sink_data(8),errors(28))
    got_sink_data(8:10) = got_sink_data(8)
    call read_from_hdf5(xyzmh_ptmass(11,1:nptmass),'tlast',group_id,got_sink_data(11),errors(29))
    call read_from_hdf5(vxyz_ptmass(:,1:nptmass),'vxyz',group_id,got_sink_vels,errors(30))
 endif

 ! Close the sinks group
 call close_hdf5group(group_id, errors(31))

 error = maxval(abs(errors))

end subroutine read_hdf5_arrays

end module utils_dumpfiles_hdf5
