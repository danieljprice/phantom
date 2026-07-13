!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine comparing time between dumps
!
! :References: None
!
! :Owner: Camille Landri
!
! :Runtime parameters: None
!
! :Dependencies: eos, hdf5, io, krome_main, krome_user, neighkdtree,
!   omp_lib, part, physcon, raytracer, units
!
 use krome_user, only:krome_nmols
 use part,       only: maxp
 use raytracer,  only: get_all_tau
 use hdf5
 use omp_lib

 implicit none
 character(len=20), parameter, public :: analysistype = 'krome'
 public :: do_analysis
 logical, allocatable :: mask(:)

 real(8), allocatable    :: abundance(:,:), abundance_prev(:,:), one(:)
 character(len=16)    :: abundance_label(krome_nmols)
 integer(8), allocatable :: iorig_old(:)
 integer, allocatable :: iprev(:)
 logical :: done_init = .false.
 real :: AuvAv = 4.65, albedo = 0.5

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,        only: isdead_or_accreted, iorig, rhoh, nptmass, xyzmh_ptmass, iReff, iboundary, igas, iphase, iamtype
 use neighkdtree, only:build_tree
 use units,       only: utime,unit_density,udist
 use physcon,     only: atomic_mass_unit
 use eos,         only: get_temperature, ieos, gamma,gmw, init_eos
 use io,          only: fatal, iverbose
 use krome_main,  only: krome_init, krome
 use krome_user,  only: krome_get_names,krome_set_user_Auv,krome_set_user_xi,&
                        krome_set_user_alb,krome_set_user_AuvAv
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real, save    :: tprev = 0.
 integer, save :: nprev = 0
 real          :: dt_cgs, rho_cgs
 real          :: numberdensity, T_gas, gammai, mui, AUV, xi
 real          :: abundance_part(krome_nmols), Y(krome_nmols), column_density(npart), xyzh_copy(4,npart)
 real          :: max_radius, radius, tstart
 integer       :: i, j, isize=0, ierr, completed_iterations, npart_copy = 0, hdferr, i_radius = 1

#ifdef __GFORTRAN__
   print*, "Setting number of threads to 1 (KROME is not thread-safe when compiled with gfortran)"
   call omp_set_num_threads(1)
#else
   print*, "running with ", omp_get_max_threads(), " threads"
#endif

 if (.not.done_init) then
    done_init = .true.

   ! initialize HDF5
   call h5open_f(hdferr)
   if (hdferr /= 0) then
      print*,'ERROR: Failed to initialize HDF5'
      return
   endif

    print*, "initialising KROME"
    call krome_init()
    print*, "Initialised KROME"

    abundance_label(:) = krome_get_names()
    allocate(abundance(krome_nmols,maxp))
    abundance = 0.
    allocate(abundance_prev(krome_nmols,maxp))
    abundance_prev = 0.
    allocate(one(maxp))
    one = 1.
    allocate(iorig_old(maxp))
    iorig_old = 0
    allocate(iprev(maxp))
    iprev = 0

    !check if a file with abundances already exists, if so read it in and use it as initial abundances
    inquire(file=trim(dumpfile)//'.h5', size=isize)
    if (isize > 0) then
         print*, "Found existing abundance file, reading in abundances"
         call read_chem(npart, dumpfile)
    else
       print*, "setting abundances"
       !$omp parallel do default(none) &
       !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev) &
       !$omp shared(abundance,abundance_prev,particlemass,unit_density) &
       !$omp shared(ieos,rho_cgs,T_gas,j) &
       !$omp private(i,abundance_part)
       do i=1, npart
          if (.not.isdead_or_accreted(xyzh(4,i))) then
             call chem_init(abundance_part)
             abundance(:,i) = abundance_part
          endif
       enddo
       call write_chem(npart, dumpfile)
    endif
    call init_eos(ieos, ierr)
    if (ierr /= 0) call fatal(analysistype, "Failed to initialise EOS")

   else
    dt_cgs = (time - tprev)*utime
    completed_iterations = 0
    print*, "not first step data, timestep = ",dt_cgs, "npart = ",npart, "nprev = ",nprev
    print*, "Building neighbour tree..."
    tstart = omp_get_wtime()
    xyzmh_ptmass(iReff,1) = 2.
    npart_copy = npart
    xyzh_copy = xyzh(:,:npart)
    call build_tree(npart_copy,npart_copy,xyzh_copy,vxyzu)
    print*, "        - Took ", omp_get_wtime() - tstart, " seconds"

    print*, "Calculating column density..."
    tstart = omp_get_wtime()
    call get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, one, 5, .false., column_density)
    max_radius = 0.0
    do i = 1, npart
       if (.not.isdead_or_accreted(xyzh(4, i))) then
          radius = sqrt(xyzh(1, i)**2 + xyzh(2, i)**2 + xyzh(3, i)**2)
          if (radius > max_radius) then
             max_radius = radius
             i_radius = i
          endif
       endif
    enddo
    column_density = column_density + rhoh(xyzh(4,i_radius),particlemass)*unit_density * max_radius * udist
    print*, "        - Took ", omp_get_wtime() - tstart, " seconds"

    print*, "Running KROME"
    tstart = omp_get_wtime()
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev,iverbose) &
    !$omp shared(abundance,abundance_prev,particlemass,unit_density,udist,iphase) &
    !$omp shared(ieos,gamma,gmw,time,completed_iterations,column_density,AuvAv,albedo) &
    !$omp private(i,j,abundance_part,Y,rho_cgs,numberdensity,T_gas,gammai,mui,AUV,xi)

    outer: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          inner: do j=1,nprev
             if (iorig(i) == iorig_old(j)) then
                iprev(i) = j
                exit inner
             endif
          enddo inner

          if (j == iprev(i)) then
             abundance_part(:) = abundance_prev(:,iprev(i))
          else
             call chem_init(abundance_part)
          endif
          if (iamtype(iphase(i)) /= iboundary .and. i > 2460) then ! 2460 is the amount of boundary particles
            !Thermodynamic quantities
             rho_cgs = rhoh(xyzh(4,i),particlemass)*unit_density
             gammai = gamma
             mui    = gmw
             numberdensity = rho_cgs / (mui * atomic_mass_unit)
             T_gas = get_temperature(ieos,xyzh(1:3, i),rhoh(xyzh(4,i),particlemass),vxyzu(:,i),gammai,mui)
             T_gas = max(T_gas,20.0d0)

             !Radiation quantities
             AUV = AuvAv * column_density(i) / (mui * atomic_mass_unit) / 1.87e21
             xi = get_xi(AUV)

             call krome_set_user_Auv(AUV)
             call krome_set_user_xi(xi)
             call krome_set_user_alb(ALBEDO)
             call krome_set_user_AuvAv(AuvAv)

             Y = abundance_part*numberdensity
             call krome(Y,T_gas,dt_cgs)
             abundance_part = Y/numberdensity
             abundance(:,i) = abundance_part
          endif
       endif
       if (iverbose > 1) then
          !$omp atomic
          completed_iterations = completed_iterations + 1
          print*, 'Completed ', completed_iterations, ' of ', npart
       endif
    enddo outer
   print*, "        - Took ", omp_get_wtime() - tstart, " seconds"

   tstart = omp_get_wtime()
   call write_chem(npart, dumpfile)
   print*, "        - Took ", omp_get_wtime() - tstart, " seconds"
 endif

 ! keep track of previous abundances for next dump
 nprev = npart
 tprev = time
 iorig_old(1:npart) = iorig(1:npart)
 abundance_prev(:,1:npart) = abundance(:,1:npart)
end subroutine do_analysis

real function get_xi(AUV)
 use physcon, only:pi
 real, intent(in) :: AUV
 real :: xi
 real :: W(6), GA(6), ceta
 integer :: i

 W(1) = 0.17132449
 W(2) = 0.36076157
 W(3) = 0.46791393
 W(4) = W(1)
 W(5) = W(2)
 W(6) = W(3)
 GA(1) = 0.93246951
 GA(2) = 0.66120939
 GA(3) = 0.23861919
 GA(4) = -GA(1)
 GA(5) = -GA(2)
 GA(6) = -GA(3)

 xi = 0.0
 do i=1,6
    ceta = (pi*GA(i)+pi)/2.0
    xi=xi+(W(i)*(sin(ceta)*exp((-AUV*ceta)/sin(ceta))))
 enddo
 xi = (pi/4.0)*xi

 get_xi = xi

end function get_xi

subroutine write_chem(npart, dumpfile)
 integer,          intent(in) :: npart
 character(len=*), intent(in) :: dumpfile
 integer                      :: hdferr
 integer(hid_t) :: file_id, dspace_id, dset_id, group_id, type_id
 integer(hsize_t) :: dims(2)
 integer(hsize_t) :: dims_labels(1)
 integer(size_t)  :: str_len

 ! create HDF5 file
 print "(1x,a)",'Writing to '//trim(dumpfile)//'.h5'
 call h5fcreate_f(trim(dumpfile)//'.h5', H5F_ACC_TRUNC_F, file_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create HDF5 file'
    return
 endif

  ! create group for particle data
 call h5gcreate_f(file_id, 'chemistry', group_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create HDF5 group'
    return
 endif

 ! create dataspace for particle datasets
 dims(1) = krome_nmols
 dims(2) = npart
 call h5screate_simple_f(2, dims, dspace_id, hdferr)
  if (hdferr /= 0) then
    print*,'ERROR: Failed to create HDF5 dataspace'
    return
 endif

 ! create one 2D dataset for all species abundances:
 call h5dcreate_f(group_id, 'abundances', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create dataset abundances'
    return
 endif

 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, abundance(:,1:npart), dims, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to write dataset abundances'
    return
 endif

 call h5dclose_f(dset_id, hdferr)

 ! store species labels so abundances row i is always identifiable
 dims_labels(1) = krome_nmols
 call h5screate_simple_f(1, dims_labels, dspace_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create HDF5 dataspace for species labels'
    return
 endif

 call h5tcopy_f(H5T_FORTRAN_S1, type_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create HDF5 string datatype for species labels'
    return
 endif

 str_len = len(abundance_label(1))
 call h5tset_size_f(type_id, str_len, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to set HDF5 string size for species labels'
    return
 endif

 call h5dcreate_f(group_id, 'species_labels', type_id, dspace_id, dset_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create dataset species_labels'
    return
 endif

 call h5dwrite_f(dset_id, type_id, abundance_label, dims_labels, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to write dataset species_labels'
    return
 endif

 call h5dclose_f(dset_id, hdferr)
 call h5tclose_f(type_id, hdferr)
 call h5sclose_f(dspace_id, hdferr)
 call h5gclose_f(group_id, hdferr)
 call h5fclose_f(file_id, hdferr)
end subroutine write_chem

subroutine read_chem(npart, dumpfile)
 integer,          intent(in) :: npart
 character(len=*), intent(in) :: dumpfile
 integer                      :: hdferr
 integer :: i
 integer(hid_t) :: file_id, dset_id, group_id, filespace_id, type_id
 integer(hsize_t) :: max_dims(2), file_dims(2)
 integer(hsize_t) :: max_dims_labels(1), file_dims_labels(1)
 integer(size_t)  :: str_len
 character(len=16) :: labels_file(krome_nmols)

 ! open HDF5 file
 print "(1x,a)",'Reading from '//trim(dumpfile)//'.h5'
 call h5fopen_f(trim(dumpfile)//'.h5', H5F_ACC_RDONLY_F, file_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to open HDF5 file'
    return
 endif

  ! open group for particle data
 call h5gopen_f(file_id, 'chemistry', group_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to open HDF5 group'
    return
 endif
 ! open matrix dataset storing all species abundances
 call h5dopen_f(group_id, 'abundance_matrix', dset_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to open dataset abundance_matrix'
    return
 endif

 call h5dget_space_f(dset_id, filespace_id, hdferr)
 call h5sget_simple_extent_dims_f(filespace_id, file_dims, max_dims, hdferr)
 if (file_dims(1) /= krome_nmols .or. file_dims(2) /= npart) then
    print*,'ERROR: abundance_matrix shape mismatch in HDF5 file'
    return
 endif

 call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, abundance(:,1:npart), file_dims, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to read dataset abundance_matrix'
    return
 endif

 call h5dclose_f(dset_id, hdferr)
 call h5sclose_f(filespace_id, hdferr)

 ! read and validate species labels to guarantee row-label mapping
 call h5dopen_f(group_id, 'species_labels', dset_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to open dataset species_labels'
    return
 endif

 call h5dget_space_f(dset_id, filespace_id, hdferr)
 call h5sget_simple_extent_dims_f(filespace_id, file_dims_labels, max_dims_labels, hdferr)
 if (file_dims_labels(1) /= krome_nmols) then
    print*,'ERROR: species_labels length mismatch in HDF5 file'
    return
 endif

 call h5tcopy_f(H5T_FORTRAN_S1, type_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create HDF5 string datatype for species_labels read'
    return
 endif

 str_len = len(labels_file(1))
 call h5tset_size_f(type_id, str_len, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to set HDF5 string size for species_labels read'
    return
 endif

 call h5dread_f(dset_id, type_id, labels_file, file_dims_labels, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to read dataset species_labels'
    return
 endif

 do i = 1, krome_nmols
    if (trim(adjustl(labels_file(i))) /= trim(adjustl(abundance_label(i)))) then
       print*,'ERROR: species_labels mismatch at index ', i
       print*,'       file: ', trim(adjustl(labels_file(i)))
       print*,'       run : ', trim(adjustl(abundance_label(i)))
       return
    endif
 enddo

 call h5dclose_f(dset_id, hdferr)
 call h5sclose_f(filespace_id, hdferr)
 call h5tclose_f(type_id, hdferr)
 call h5gclose_f(group_id, hdferr)
 call h5fclose_f(file_id, hdferr)
end subroutine read_chem

subroutine chem_init(abundance_part)
 use krome_user, only:krome_idx_H2,krome_idx_He,krome_idx_CO,krome_idx_C2H2,&
       krome_idx_HCN,krome_idx_N2,krome_idx_SiC2,krome_idx_CS,&
       krome_idx_SiS,krome_idx_SiO,krome_idx_CH4,krome_idx_H2O,&
       krome_idx_HCl,krome_idx_C2H4,krome_idx_NH3,krome_idx_HCP,&
       krome_idx_HF,krome_idx_H2S,krome_idx_e,krome_get_electrons
 real, intent(out) :: abundance_part(krome_nmols)

 ! Initial abundances for the krome model taken from Agúndez et al. (2020)
 ! H2, He, CO, C2H2, HCN, N2, SiC2, CS, SiS, SiO, CH4, H2O, HCl, C2H4, NH3, HCP, HF, H2S
 abundance_part(:)              = 0.
 abundance_part(krome_idx_H2)   = 0.5d0
 abundance_part(krome_idx_He)   = 8.5d-2
 abundance_part(krome_idx_CO)   = 4d-4
 abundance_part(krome_idx_C2H2) = 2.19d-5
 abundance_part(krome_idx_HCN)  = 2.045d-5
 abundance_part(krome_idx_N2)   = 2d-5
 abundance_part(krome_idx_SiC2) = 9.35d-6
 abundance_part(krome_idx_CS)   = 5.3d-6
 abundance_part(krome_idx_SiS)  = 2.99d-6
 abundance_part(krome_idx_SiO)  = 2.51d-6
 abundance_part(krome_idx_CH4)  = 1.75d-6
 abundance_part(krome_idx_H2O)  = 1.275d-6
 abundance_part(krome_idx_HCl)  = 1.625d-7
 abundance_part(krome_idx_C2H4) = 3.425d-8
 abundance_part(krome_idx_NH3)  = 3d-8
 abundance_part(krome_idx_HCP)  = 1.25d-8
 abundance_part(krome_idx_HF)   = 8.5d-9
 abundance_part(krome_idx_H2S)  = 2d-9
 abundance_part(krome_idx_e)    = krome_get_electrons(abundance_part(:))

end subroutine chem_init

end module analysis
