module output_hdf5
 use hdf5,             only:HID_T,h5open_f,h5fcreate_f,H5F_ACC_TRUNC_F,h5fclose_f,h5close_f,h5gcreate_f,h5gclose_f
 use utils_outputhdf5, only:write_to_hdf5
implicit none
 public :: open_hdf5file, close_hdf5file, write_hdf5_header, write_hdf5_arrays

 integer(HID_T), public :: outputfile_id

 private

contains

subroutine open_hdf5file(filename,file_id,error)
 character(len=*), intent(in)  :: filename
 integer(HID_T),   intent(out) :: file_id
 integer,          intent(out) :: error
 integer :: errors(2)
 call h5open_f(errors(1))                                     ! Initialise Fortran h5 interfaces
 call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,errors(2)) ! Create file
 error = maxval(abs(errors))
end subroutine open_hdf5file

subroutine close_hdf5file(file_id,error)
 integer(HID_T), intent(in) :: file_id
 integer, intent(out) :: error
 integer :: errors(2)
 call h5fclose_f(file_id,errors(2)) ! Close the file
 call h5close_f(errors(1))          ! Close Fortran h5 interfaces
 error = maxval(abs(errors))
end subroutine close_hdf5file

subroutine write_hdf5_header(file_id,fileident,maxtypes,nblocks,isink,nptmass,ndustlarge,ndustsmall,idust,                  &
                            phantom_version_major,phantom_version_minor,phantom_version_micro,                    &
                            nparttot,npartoftypetot,iexternalforce,ieos,t,dtmax,gamma,rhozero,                    &
                            polyk,hfact,tolh,C_cour,C_force,alpha,alphau,alphaB,polyk2,qfacdisc,                  &
                            massoftype,Bextx,Bexty,Bextz,xmin,xmax,ymin,ymax,zmin,zmax,get_conserv,               &
                            etot_in,angtot_in,totmom_in,mdust_in,grainsize,graindens,udist,umass,utime,unit_Bfield)
 integer(HID_T),   intent(in) :: file_id
 character(len=*), intent(in) :: fileident
 integer, intent(in) :: maxtypes,              & !integer (default)
                        nblocks,               & !integer (default)
                        isink,                 & !integer (default)
                        nptmass,               & !integer (default)
                        ndustlarge,            & !integer (default)
                        ndustsmall,            & !integer (default)
                        idust,                 & !integer (default)
                        phantom_version_major, & !integer (default)
                        phantom_version_minor, & !integer (default)
                        phantom_version_micro, & !integer (default)
                        nparttot,              & !int*8
                        npartoftypetot(:),     & !int(:)*8
                        iexternalforce,        & !int*4
                        ieos                     !int*4
 real, intent(in) :: t,                     & !real (default)
                     dtmax,                 & !real (default)
                     gamma,                 & !real (default)
                     rhozero,               & !real (default)
                     polyk,                 & !real (default)
                     hfact,                 & !real (default)
                     tolh,                  & !real (default)
                     C_cour,                & !real (default)
                     C_force,               & !real (default)
                     alpha,                 & !real (default)
                     alphau,                & !real (default)
                     alphaB,                & !real (default)
                     polyk2,                & !real (default)
                     qfacdisc,              & !real (default)
                     massoftype(:),         & !real(:) (default)
                     Bextx,                 & !real (default)
                     Bexty,                 & !real (default)
                     Bextz,                 & !real (default)
                     xmin,                  & !real (default)
                     xmax,                  & !real (default)
                     ymin,                  & !real (default)
                     ymax,                  & !real (default)
                     zmin,                  & !real (default)
                     zmax,                  & !real (default)
                     get_conserv,           & !real (default)
                     etot_in,               & !real (default)
                     angtot_in,             & !real (default)
                     totmom_in,             & !real (default)
                     mdust_in(:),           & !real(:) (default)
                     grainsize(:),          & !real(:) (default)
                     graindens(:),          & !real(:) (default)
                     udist,                 & !real*8
                     umass,                 & !real*8
                     utime,                 & !real*8
                     unit_Bfield              !real*8
 integer(HID_T) :: group_id
 integer :: error

    ! Create header group
    call h5gcreate_f(file_id,'header',group_id,error)

    ! Write things to header group
    call write_to_hdf5(fileident,'fileident',group_id)
    ! call write_to_hdf5(int(nparttot),'nparttot',group_id)
    call write_to_hdf5(maxtypes,'ntypes',group_id)
    ! call write_to_hdf5(int(npartoftypetot(1:maxtypes)),'npartoftype',group_id)
    call write_to_hdf5(nblocks,'nblocks',group_id)
    call write_to_hdf5(isink,'isink',group_id)
    call write_to_hdf5(nptmass,'nptmass',group_id)
    call write_to_hdf5(ndustlarge,'ndustlarge',group_id)
    call write_to_hdf5(ndustsmall,'ndustsmall',group_id)
    call write_to_hdf5(idust,'idust',group_id)
    call write_to_hdf5(phantom_version_major,'majorv',group_id)
    call write_to_hdf5(phantom_version_minor,'minorv',group_id)
    call write_to_hdf5(phantom_version_micro,'microv',group_id)
    call write_to_hdf5(nparttot,'nparttot',group_id)
    ! call write_to_hdf5(int(maxtypes,kind=8),'ntypes',group_id)
    call write_to_hdf5(npartoftypetot(1:maxtypes),'npartoftype',group_id)
    call write_to_hdf5(iexternalforce,'iexternalforce',group_id)
    call write_to_hdf5(ieos,'ieos',group_id)
    call write_to_hdf5(t,'time',group_id)
    call write_to_hdf5(dtmax,'dtmax',group_id)
    call write_to_hdf5(gamma,'gamma',group_id)
    call write_to_hdf5(rhozero,'rhozero',group_id)
    call write_to_hdf5(1.5*polyk,'RK2',group_id)
    call write_to_hdf5(hfact,'hfact',group_id)
    call write_to_hdf5(tolh,'tolh',group_id)
    call write_to_hdf5(C_cour,'C_cour',group_id)
    call write_to_hdf5(C_force,'C_force',group_id)
    call write_to_hdf5(alpha,'alpha',group_id)
    call write_to_hdf5(alphau,'alphau',group_id)
    call write_to_hdf5(alphaB,'alphaB',group_id)
    call write_to_hdf5(polyk2,'polyk2',group_id)
    call write_to_hdf5(qfacdisc,'qfacdisc',group_id)
    call write_to_hdf5(massoftype,'massoftype',group_id) ! array
    call write_to_hdf5(Bextx,'Bextx',group_id)
    call write_to_hdf5(Bexty,'Bexty',group_id)
    call write_to_hdf5(Bextz,'Bextz',group_id)
    call write_to_hdf5(0.,'dum',group_id)
    ! if (iexternalforce /= 0) call write_headeropts_extern(iexternalforce,hdr,t,ierr)
    call write_to_hdf5(xmin,'xmin',group_id)
    call write_to_hdf5(xmax,'xmax',group_id)
    call write_to_hdf5(ymin,'ymin',group_id)
    call write_to_hdf5(ymax,'ymax',group_id)
    call write_to_hdf5(zmin,'zmin',group_id)
    call write_to_hdf5(zmax,'zmax',group_id)
    call write_to_hdf5(get_conserv,'get_conserv',group_id)
    call write_to_hdf5(etot_in,'etot_in',group_id)
    call write_to_hdf5(angtot_in,'angtot_in',group_id)
    call write_to_hdf5(totmom_in,'totmom_in',group_id)
    call write_to_hdf5(mdust_in,'mdust_in',group_id)
    call write_to_hdf5(grainsize,'grainsize',group_id)
    call write_to_hdf5(graindens,'graindens',group_id)
    call write_to_hdf5(udist,'udist',group_id)
    call write_to_hdf5(umass,'umass',group_id)
    call write_to_hdf5(utime,'utime',group_id)
    call write_to_hdf5(unit_Bfield,'umagfd',group_id)

    ! Close the header group
    call h5gclose_f(group_id, error)
end subroutine write_hdf5_header

subroutine write_hdf5_arrays(file_id,xyzh,vxyzu,iphase,pressure,alphaind,dtind,poten,xyzmh_ptmass,vxyz_ptmass, &
                             Bxyz,Bevol,divcurlB,divBsymm,eta_nimhd,dustfrac,tstop,deltav,dustprop,st,         &
                             abundance,temperature,divcurlv,luminosity,beta_pr,                                &
                             const_av,ind_timesteps,gravity,nptmass,mhd,maxBevol,ndivcurlB,mhd_nonideal,       &
                             use_dust,use_dustfrac,use_dustgrowth,h2chemistry,store_temperature,ndivcurlv,     &
                             lightcurve,prdrag,isothermal                                                      )
 integer(HID_T),         intent(in) :: file_id
 real, dimension(:),     intent(in) :: pressure,dtind,poten,beta_pr,divBsymm,st,temperature,luminosity
 real, dimension(:,:),   intent(in) :: xyzh,vxyzu,alphaind,Bxyz,Bevol,divcurlB,eta_nimhd,xyzmh_ptmass,vxyz_ptmass
 real, dimension(:,:),   intent(in) :: dustfrac,tstop,dustprop,abundance,divcurlv
 real, dimension(:,:,:), intent(in) :: deltav
 integer, intent(in) :: iphase(:)
 integer, intent(in) :: nptmass,maxBevol,ndivcurlB,ndivcurlv
 logical, intent(in) :: const_av,ind_timesteps,gravity,mhd,mhd_nonideal,use_dust,use_dustfrac,use_dustgrowth
 logical, intent(in) :: h2chemistry,store_temperature,lightcurve,prdrag,isothermal
 integer(HID_T) :: group_id
 integer :: error

 ! Create particles group
 call h5gcreate_f(file_id,'particles',group_id,error)

 ! Main arrays
 call write_to_hdf5(xyzh(1:3,:),'xyz',group_id)
 call write_to_hdf5(xyzh(4,:),'h',group_id)
 call write_to_hdf5(vxyzu(1:3,:),'vxyz',group_id)
 if (.not.isothermal) call write_to_hdf5(vxyzu(4,:),'u',group_id)
 call write_to_hdf5(iphase,'itype',group_id)
 call write_to_hdf5(pressure,'pressure',group_id)

 if (.not.const_av)  call write_to_hdf5(alphaind,'alpha',group_id) ! Viscosity
 if (ind_timesteps)  call write_to_hdf5(dtind,'dt',group_id)       ! Individual timesteps
 if (gravity)        call write_to_hdf5(poten,'poten',group_id)

 ! MHD arrays
 if (mhd) then
    call write_to_hdf5(Bxyz,'Bxyz',group_id)
    if (maxBevol >= 4) then
       call write_to_hdf5(Bevol(4,:),'psi',group_id)
    endif
    if (ndivcurlB >= 1) then
       call write_to_hdf5(divcurlB(1,:),'divB',group_id)
       call write_to_hdf5(divcurlB(2:4,:),'curlBxyz',group_id)
    else
       call write_to_hdf5(divBsymm,'divBsymm',group_id)
    endif
    if (mhd_nonideal) then
       call write_to_hdf5(eta_nimhd(1,:),'eta_{OR}',group_id)
       call write_to_hdf5(eta_nimhd(2,:),'eta_{HE}',group_id)
       call write_to_hdf5(eta_nimhd(3,:),'eta_{AD}',group_id)
       call write_to_hdf5(eta_nimhd(4,:),'ne/n'    ,group_id)
    endif
 endif

 ! Dust arrays
 if (use_dust) then
    call write_to_hdf5(dustfrac,'dustfrac',group_id)
    call write_to_hdf5(tstop,'tstop',group_id)
 endif
 if (use_dustfrac) call write_to_hdf5(deltav,'deltavxyz',group_id)
 if (use_dustgrowth) then
    call write_to_hdf5(dustprop(1,:),'grainsize' ,group_id)
    call write_to_hdf5(dustprop(2,:),'graindens' ,group_id)
    call write_to_hdf5(dustprop(3,:),'vrel/vfrag',group_id)
    call write_to_hdf5(dustprop(4,:),'dv_dust'   ,group_id)
    call write_to_hdf5(St,'St',group_id)
 endif

 ! Other Arrays
 if (h2chemistry)       call write_to_hdf5(abundance,'abundance',group_id)
 if (store_temperature) call write_to_hdf5(temperature,'T',group_id)
 if (ndivcurlv >= 1) then
    call write_to_hdf5(divcurlv(1,:),'divv',group_id)
    if (ndivcurlv>=4) call write_to_hdf5(divcurlv(2:4,:),'curlvxyz',group_id)
 endif
 if (lightcurve) call write_to_hdf5(luminosity,'luminosity',group_id)
 if (prdrag)     call write_to_hdf5(beta_pr,'beta_pr',group_id)

 ! Close the particles group
 call h5gclose_f(group_id, error)

 ! Create sink group
 call h5gcreate_f(file_id,'sinks',group_id,error)
 if (nptmass > 0) then
    call write_to_hdf5(xyzmh_ptmass(1:3,:),'xyz',group_id)
    call write_to_hdf5(xyzmh_ptmass(4,:),'m',group_id)
    call write_to_hdf5(xyzmh_ptmass(5,:),'h',group_id)
    call write_to_hdf5(xyzmh_ptmass(6,:),'hsoft',group_id)
    call write_to_hdf5(xyzmh_ptmass(7,:),'maccreted',group_id)
    call write_to_hdf5(xyzmh_ptmass(8:10,:),'spinxyz',group_id)
    call write_to_hdf5(xyzmh_ptmass(11,:),'tlast',group_id)
    call write_to_hdf5(vxyz_ptmass,'vxyz',group_id)
 endif
 ! Close the sink group
 call h5gclose_f(group_id, error)

end subroutine write_hdf5_arrays

end module output_hdf5
