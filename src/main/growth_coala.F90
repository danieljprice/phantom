!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module growth_coala
!
! Interface to library for dust growth and fragmentation
!   using COALA solver. This module can be compiled
!   with a special rule in the Makefile to link against the
!   COALA library
!
! :References:
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon, units, part, dust
!
 use part,      only:ndusttypes,grainsize,graindens
 use physcon,   only:pi
#ifdef COALA
 use precision, only:wp
#endif
 implicit none
#ifndef COALA
 integer, parameter :: wp = kind(0.d0)
#endif
 integer :: order_growth = 0 ! order of the DG polynomials (order of scheme - 1)
 integer :: Q_coag = 10 ! number of points for Gauss quadrature
 real :: alpha_turb = 1.5
 real :: alpha_turb_disk = 1e-4
 real :: n_trans_coag = 1e9
 integer :: turb_grow = 0
 integer :: brow_grow = 0
 integer :: drift_grow = 0
 real :: K0 = 1.0
 integer :: kernel = 3

 ! COALA arrays
 real(wp), allocatable :: tabflux_coag_k0(:,:,:)
 real(wp), allocatable :: tabflux_coag(:,:,:,:,:)
 real(wp), allocatable :: tabintflux_coag(:,:,:,:,:,:)
 real(wp), allocatable :: massgrid(:)
 real(wp), allocatable :: massbins(:)
 real(wp), allocatable :: massmeanlog(:)
#ifdef COALA
 real(wp), allocatable :: vecnodes(:)
 real(wp), allocatable :: vecweights(:)
 real(wp), allocatable :: mat_coeffs_leg(:,:)
#endif

 public :: init_growth_coala, read_options_growth_coala, write_options_growth_coala, get_growth_rate_coala
 public :: check_coagflux_array, K0 ! to avoid compiler warning
 private

contains

!-----------------------------------------------------------------------
!+
!  Sanity check for coagulation flux arrays
!+
!-----------------------------------------------------------------------
subroutine check_coagflux_array(array,array_name,ierr)
 use io, only:error,fatal
 real(wp), intent(in) :: array(:)
 character(len=*), intent(in) :: array_name
 integer, intent(inout) :: ierr
 real(wp) :: val_max,val_min
 logical :: has_nan,has_inf,has_negative
 integer :: i

 ! Check for NaN using any(isnan())
 ! Note: isnan() is available in gfortran; for other compilers may need ieee_is_nan
 has_nan = any(isnan(array))
 
 ! Check for Inf (very large values)
 has_inf = any(abs(array) > huge(1.0_wp)*0.1_wp)
 
 ! Check for negative values
 has_negative = any(array < 0.0_wp)
 
 ! Find min and max (excluding NaN/Inf)
 val_max = -huge(1.0_wp)
 val_min = huge(1.0_wp)
 do i=1,size(array)
    if (array(i) == array(i) .and. abs(array(i)) < huge(1.0_wp)*0.1_wp) then
       if (array(i) > val_max) val_max = array(i)
       if (array(i) < val_min) val_min = array(i)
    endif
 enddo
 
 if (has_nan .or. has_inf) then
    print*,'ERROR: Invalid values (NaN/Inf) found in ',trim(array_name)
    call fatal('check_coagflux_array','Invalid values in '//trim(array_name))
    ierr = 1
    return
 endif
 
 print*,trim(array_name),': min = ',val_min,', max = ',val_max
 if (has_negative) then
    call error('check_coagflux_array','Negative values found in '//trim(array_name))
 endif

end subroutine check_coagflux_array

!-----------------------------------------------------------------------
!+
!  Initialize growth with COALA
!+
!-----------------------------------------------------------------------
subroutine init_growth_coala(ierr)
 use io,    only:fatal,error
 use units, only:udist
#ifdef COALA
 use io,                                only:id,master
 use coala_polynomials_legendre,        only:compute_mat_coeffs
 use coala_GQLeg_nodes_weights,         only:GQLeg_nodes,GQLeg_weights
 use coala_generate_tabflux_tabintflux, only:compute_coagtabflux_GQ_k0, &
                                             compute_coagtabflux_GQ, &
                                             compute_coagtabintflux_GQ
#endif
 integer, intent(out) :: ierr
 integer :: idust
 real(wp), dimension(1:ndusttypes+1) :: sdust
 real(wp), dimension(1:ndusttypes) :: d_grain,l_grain
 real(wp) :: ratio

 ierr = 0

 ! Allocate COALA arrays
 allocate(tabflux_coag_k0(ndusttypes,ndusttypes,ndusttypes))
 if (order_growth > 0) then
    allocate(tabflux_coag(ndusttypes,ndusttypes,ndusttypes,order_growth+1,order_growth+1))
    allocate(tabintflux_coag(ndusttypes,order_growth+1,ndusttypes,ndusttypes,order_growth+1,order_growth+1))
 endif
 allocate(massgrid(1:ndusttypes+1),stat=ierr)
 allocate(massbins(1:ndusttypes),stat=ierr)
 allocate(massmeanlog(1:ndusttypes),stat=ierr)

 massgrid = 0.0_wp
 massbins = 0.0_wp
 massmeanlog = 0.0_wp

 ! Convert grain sizes and densities from code units to cgs, then to code units
 ! In Phantom, grainsize and graindens are already in code units
 do idust=1,ndusttypes
    l_grain(idust) = grainsize(idust)
    d_grain(idust) = graindens(idust)
 enddo

 ! Build size grid (sdust) - boundaries of size bins
 ! Reconstruct from grainsize array, which contains geometric means of bin boundaries
 ! For a logarithmic grid: grainsize(i) = sqrt(grid(i)*grid(i+1))
 ! This means: grid(i+1) = grainsize(i)^2 / grid(i)
 ! For a logspace grid, the ratio between boundaries is constant: ratio = grainsize(2)/grainsize(1)
 ! Therefore: grid(1) = grainsize(1) / sqrt(ratio) = grainsize(1)^(3/2) / sqrt(grainsize(2))
 if (ndusttypes > 1) then
    ! Compute the ratio between adjacent grainsizes (should be constant for logspace)
    ratio = grainsize(2) / grainsize(1)
    
    ! First bin boundary: reconstructed from logspace relationship
    sdust(1) = grainsize(1) / sqrt(ratio)
    
    ! Middle boundaries: geometric mean of adjacent grainsizes
    do idust=2,ndusttypes
       sdust(idust) = sqrt(grainsize(idust-1)*grainsize(idust))
    end do
    
    ! Last bin boundary: use the same ratio
    sdust(ndusttypes+1) = grainsize(ndusttypes) * sqrt(ratio)
    if (id==master) then
      print "(/,1x,a)",'---------------------------------------------------------------------------'
      print "(1x,a,/)",'COALA dust growth is ON: please cite Lombart et al. (2021) MNRAS 501, 4298'
       print "(2x,a,1pg0.3,a)",'minimum grain size = ',sdust(1)*udist,' cm'
       print "(2x,a,1pg0.3,a)",'maximum grain size = ',sdust(ndusttypes+1)*udist,' cm'
       print "(2x,a,i0,a)",'number of dust types = ',ndusttypes
       print "(2x,a,i0,a)",'number of points for Gauss quadrature = ',Q_coag
       print "(2x,a,i0,a)",'order of the DG polynomials = ',order_growth
       print "(2x,a,i0,a)",'kernel = ',kernel
       print "(2x,a,1pg0.3,a)",'alpha_turb = ',alpha_turb
       print "(2x,a,1pg0.3,a)",'alpha_turb_disk = ',alpha_turb_disk
       print "(1x,a)",'---------------------------------------------------------------------------'
    endif
 else
    ! Single dust type - use a reasonable range
    sdust(1) = 0.5_wp * grainsize(1)
    sdust(2) = 1.5_wp * grainsize(1)
 endif

 ! Build massgrid and massbins
 do idust=1,ndusttypes
    massgrid(idust) = 4.0_wp/3.0_wp*pi*d_grain(idust)*(sdust(idust))**3
    massgrid(idust+1) = 4.0_wp/3.0_wp*pi*d_grain(idust)*(sdust(idust+1))**3
    massbins(idust) = 0.5_wp*(massgrid(idust+1)+massgrid(idust))
    massmeanlog(idust) = sqrt(massgrid(idust+1)*massgrid(idust))
 enddo

#ifdef COALA
 ! Allocate and compute Gauss-Legendre quadrature nodes and weights
 allocate(vecnodes(Q_coag),stat=ierr)
 allocate(vecweights(Q_coag),stat=ierr)
 vecnodes = 0.0_wp
 vecweights = 0.0_wp
 call GQLeg_nodes(Q_coag,vecnodes)
 call GQLeg_weights(Q_coag,vecweights)

 ! Ballistic kernel
 ! 0 = constant, 1 = additive, 2 = ballistic (cross section with delta v from hydro),
 ! 3 = ballistic with delta v from Brownian motion, 4 = ballistic with delta v from turbulence
 kernel = 2

 ! Normalisation of the geometrical cross-section
 ! K0 = pi*(4.0_wp/3.0_wp*pi*d_grain(1))**(-2.0_wp/3.0_wp)
 ! Using first grain density for normalization
 K0 = pi*(4.0_wp/3.0_wp*pi*d_grain(1))**(-2.0_wp/3.0_wp)

 ! Compute Legendre polynomial coefficients if needed
 allocate(mat_coeffs_leg(order_growth+1,order_growth+1))
 call compute_mat_coeffs(order_growth,mat_coeffs_leg)

 ! Precompute coagulation flux tables
 if (order_growth == 0) then
    call compute_coagtabflux_GQ_k0(kernel,K0,Q_coag,vecnodes,vecweights, &
                                    ndusttypes,order_growth,massgrid, &
                                    mat_coeffs_leg,tabflux_coag_k0)
    
    ! Sanity checks for tabflux_coag_k0
    call check_coagflux_array(reshape(tabflux_coag_k0,[size(tabflux_coag_k0)]), &
                              'tabflux_coag_k0',ierr)
    if (ierr /= 0) return
    
 else
    call compute_coagtabflux_GQ(kernel,K0,Q_coag,vecnodes,vecweights, &
                                 ndusttypes,order_growth,massgrid, &
                                 mat_coeffs_leg,tabflux_coag)
    call compute_coagtabintflux_GQ(kernel,K0,Q_coag,vecnodes,vecweights, &
                                    ndusttypes,order_growth,massgrid, &
                                    mat_coeffs_leg,tabintflux_coag)
    
    ! Sanity checks for tabflux_coag
    call check_coagflux_array(reshape(tabflux_coag,[size(tabflux_coag)]), &
                              'tabflux_coag',ierr)
    if (ierr /= 0) return
    
    ! Sanity checks for tabintflux_coag
    call check_coagflux_array(reshape(tabintflux_coag,[size(tabintflux_coag)]), &
                              'tabintflux_coag',ierr)
    if (ierr /= 0) return
    
 endif
#endif

end subroutine init_growth_coala


!-----------------------------------------------------------------------
!+
!  Main routine that updates rhodust with COALA solver
!+
!-----------------------------------------------------------------------
subroutine get_growth_rate_coala(npart,xyzh,vxyzu,fxyzu,fext,&
                                 grainsize,dustfrac,deltav,ddustevol,dt,eos_vars)
 use physcon,              only:mH=>mass_proton_cgs
 use part,                 only:rhoh,massoftype,igas,isdead_or_accreted,ics,itemp,imu,tstop
#ifdef COALA
 use coala_interface_coag, only:coala_coag_k0,coala_coag
#endif
 integer, intent(in)  :: npart
 real, intent(in)     :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real, intent(in)     :: grainsize(:),dustfrac(:,:)
 real, intent(in)     :: deltav(:,:,:),dt
 real, intent(in)     :: eos_vars(:,:)
 real, intent(inout)  :: ddustevol(:,:)

 integer :: i,idust
 real(wp) :: rhodust_old(ndusttypes),rhodust_new(ndusttypes),t_stop(ndusttypes)
 real(wp) :: m_grain(ndusttypes),d_grain(ndusttypes),l_grain(ndusttypes)
 real(wp) :: sym_dvij(ndusttypes,ndusttypes)
 real(wp) :: eps_rhodust,cs,mu_gas
 real(wp) :: rhoi,drhodust_dti,fxi,fyi,fzi,a_gas

 ! Minimum value for rhodust in code units
 eps_rhodust = 0. !1.e-30_wp

 ! compute grain properties in code units
 do idust=1,ndusttypes
    l_grain(idust) = grainsize(idust)
    d_grain(idust) = graindens(idust)  ! grain density
    m_grain(idust) = 4.0_wp/3.0_wp*pi*d_grain(idust)*l_grain(idust)**3  ! grain mass
 enddo

 ! Loop over particles
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle

    ! Get total density (gas+dust) in code units
    rhoi = rhoh(xyzh(4,i),massoftype(igas))
    if (rhoi <= 0.0) cycle

    ! Get sound speed, temperature and molecular weight from eos_vars (already in code units)
    cs = eos_vars(ics,i)
    mu_gas = real(eos_vars(imu,i),wp)

    ! Compute rhodust in code units
    do idust=1,ndusttypes
       rhodust_old(idust) = max(rhoi * dustfrac(idust,i), eps_rhodust)
       t_stop(idust) = tstop(idust,i)
    enddo

    fxi = (fxyzu(1,i) + fext(1,i))
    fyi = (fxyzu(2,i) + fext(2,i))
    fzi = (fxyzu(3,i) + fext(3,i))
    a_gas = sqrt(fxi**2 + fyi**2 + fzi**2)

    ! Compute differential velocities dvij (symmetric matrix)
    call compute_differential_velocities(ndusttypes,t_stop,cs,rhoi,m_grain,mu_gas, &
                                         deltav(:,:,i),a_gas,sym_dvij)

#ifdef COALA
    ! Call COALA coagulation routine
    if (order_growth == 0) then
       call coala_coag_k0(ndusttypes,massgrid,tabflux_coag_k0, &
                          rhodust_old,eps_rhodust, &
                          sym_dvij,dt,rhodust_new)
    else
       call coala_coag(ndusttypes,order_growth,massgrid,massbins, &
                       mat_coeffs_leg,Q_coag,vecnodes,vecweights, &
                       tabflux_coag,tabintflux_coag, &
                       rhodust_old,eps_rhodust, &
                       sym_dvij,dt,rhodust_new)
    endif

    ! Compute drhodust/dt and convert to d(sqrt(rhodust))/dt using chain rule
    ! d(sqrt(rhodust))/dt = (1/(2*sqrt(rhodust))) * drhodust/dt
    ! Note: ddustevol stores d(sqrt(rhodust))/dt in code units
    do idust=1,ndusttypes
       drhodust_dti = (rhodust_new(idust) - rhodust_old(idust)) / dt
       
       if (rhodust_old(idust) > eps_rhodust) then
          ddustevol(idust,i) = ddustevol(idust,i) + (1.0/(2.0*sqrt(rhodust_old(idust)))) * drhodust_dti
       endif
    enddo

    !print*,i,' rhodust_old = ',rhodust_old
    !print*,i,' rhodust_new = ',rhodust_new,' dt = ',dt
    !print*,i,' ddustevol = ',ddustevol(:,i)
    !read*

 enddo
#endif

end subroutine get_growth_rate_coala

!-----------------------------------------------------------------------
!+
!  Compute differential velocities between dust grain sizes
!+
!-----------------------------------------------------------------------
subroutine compute_differential_velocities(ndusttypes,t_stop,cs,rho,m_grain,mu_gas, &
                                           deltav,a_gas,sym_dvij)
 use physcon,   only:mH=>mass_proton_cgs
 use units,     only:umass
 integer,  intent(in) :: ndusttypes
 real(wp), intent(in) :: t_stop(ndusttypes),cs,rho
 real(wp), intent(in) :: m_grain(ndusttypes),mu_gas
 real,     intent(in) :: deltav(:,:)
 real(wp), intent(in) :: a_gas
 real(wp), intent(out) :: sym_dvij(ndusttypes,ndusttypes)

 integer :: idust,jdust
 real(wp) :: vdrift_turb,vdrift_brow,vdrift_hydro
 real(wp) :: f_Stokes,St1,St2,x_stokes,t_L

 ! Compute differential velocities dvij (symmetric matrix)
 sym_dvij = 0.0_wp
 do idust=1,ndusttypes
    do jdust=1,ndusttypes
       x_stokes = t_stop(jdust)/t_stop(idust)
       f_Stokes = 3.2_wp - 1.0_wp - x_stokes + 2.0_wp/(1.0_wp+x_stokes) * &
                  (1.0_wp/2.6_wp + x_stokes**3/(1.6_wp+x_stokes))

       ! Compute t_L (large eddy turnover time) - approximate from free fall time
       !t_L = sqrt(3.0_wp*pi/(32.0_wp*6.67e-8_wp*rho*unit_density)) / utime ! convert to code units
       t_L = cs / a_gas
       St1 = t_stop(idust) / t_L
       St2 = t_stop(jdust) / t_L

       ! Symmetrize
       if (jdust > idust) then
          St1 = t_stop(jdust) / t_L
          St2 = t_stop(idust) / t_L
          x_stokes = t_stop(idust)/t_stop(jdust)
          f_Stokes = 3.2_wp - 1.0_wp - x_stokes + 2.0_wp/(1.0_wp+x_stokes) * &
                     (1.0_wp/2.6_wp + x_stokes**3/(1.6_wp+x_stokes))
       endif

       ! Turbulent velocity
       vdrift_turb = 0.0_wp
       if (turb_grow > 0) then
          vdrift_turb = sqrt(alpha_turb) * cs * sqrt(f_Stokes*St1)
          if (rho > n_trans_coag) then
             vdrift_turb = sqrt(alpha_turb_disk) * cs * sqrt(f_Stokes*St1)
          endif
       endif

       ! Hydrodynamic drift velocity (from deltav)
       vdrift_hydro = 0.0_wp
       if (drift_grow > 0) then
          vdrift_hydro = sqrt(sum((deltav(:,idust) - deltav(:,jdust))**2))
       endif

       ! Brownian motion velocity
       vdrift_brow = 0.0_wp
       if (brow_grow > 0) then
          vdrift_brow = cs * sqrt(mu_gas*mH/umass) / sqrt(pi*1.4_wp/8.0_wp) * &
                        sqrt((m_grain(idust)+m_grain(jdust))/(m_grain(idust)*m_grain(jdust)))
       endif

       ! Total differential velocity
       sym_dvij(idust,jdust) = sqrt(vdrift_turb**2 + vdrift_hydro**2 + vdrift_brow**2)
    enddo
 enddo

end subroutine compute_differential_velocities

!-----------------------------------------------------------------------
!+
!  Read growth options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_growth_coala(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(order_growth,'order_growth',db,errcount=nerr,min=0,max=10)
 call read_inopt(Q_coag,'Q_coag',db,errcount=nerr,min=5,max=20)
 call read_inopt(alpha_turb,'alpha_turb',db,errcount=nerr,min=0.0,max=10.0)
 call read_inopt(alpha_turb_disk,'alpha_turb_disk',db,errcount=nerr,min=0.0,max=1.0)
 call read_inopt(n_trans_coag,'n_trans_coag',db,errcount=nerr)
 call read_inopt(turb_grow,'turb_grow',db,errcount=nerr)
 call read_inopt(brow_grow,'brow_grow',db,errcount=nerr)
 call read_inopt(drift_grow,'drift_grow',db,errcount=nerr)
 call read_inopt(kernel,'kernel',db,errcount=nerr)

end subroutine read_options_growth_coala

!-----------------------------------------------------------------------
!+
!  Write growth options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_growth_coala(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling dust growth with COALA'
 call write_inopt(order_growth,'order_growth','order of the DG polynomials (order of scheme - 1)',iunit)
 call write_inopt(Q_coag,'Q_coag','number of points for Gauss quadrature',iunit)
 call write_inopt(alpha_turb,'alpha_turb','turbulent growth parameter',iunit)
 call write_inopt(alpha_turb_disk,'alpha_turb_disk','turbulent growth parameter for disk',iunit)
 call write_inopt(n_trans_coag,'n_trans_coag','number density threshold for coagulation (cm^-3)',iunit)
 call write_inopt(turb_grow,'turb_grow','turbulent growth (0=off, 1=on)',iunit)
 call write_inopt(brow_grow,'brow_grow','Brownian growth (0=off, 1=on)',iunit)
 call write_inopt(drift_grow,'drift_grow','drift growth (0=off, 1=on)',iunit)
 call write_inopt(kernel,'kernel','kernel type (0=constant,1=additive,2=multiplicative,3=cross section)',iunit)

end subroutine write_options_growth_coala

end module growth_coala