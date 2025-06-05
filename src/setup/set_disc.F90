!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setdisc
!
! This module contains utility routines for accretion disc setups
!
! :References:
!   Lodato & Pringle (2007), MNRAS 381, 1287-1300
!   Lodato & Price (2010), MNRAS 405, 1212-1226
!   Meru & Bate (2012), MNRAS 427, 2022-2046
!
! :Owner: Daniel Mentiplay
!
! :Runtime parameters:
!   - G           : *in code units*
!   - M_disc      : *disc mass*
!   - M_star      : *mass of central star*
!   - Qmin        : *minimum Toomre Q parameter*
!   - R_c         : *characteristic radius of the exponential taper*
!   - R_in        : *inner disc boundary*
!   - R_out       : *outer disc boundary*
!   - R_ref       : *reference radius*
!   - R_warp      : *position of warp*
!   - T_in        : *temperature (K) at R=R_in*
!   - T_out       : *temperature (K) at R=R_out*
!   - T_ref       : *temperature (K) at R=R_ref*
!   - alphaSS_max : *maximum Shakura-Sunyaev alpha viscosity in disc*
!   - alphaSS_min : *minimum Shakura-Sunyaev alpha viscosity in disc*
!   - c           : *in code units*
!   - cs0         : *sound speed at R=1*
!   - n           : *number of particles in the disc*
!   - p_index     : *power law index of surface density profile*
!   - psi_max     : *maximum warp amplitude*
!   - q_index     : *power law index of sound speed profile*
!   - sig_in      : *surface density (g/cm^2) at R=R_in*
!   - sig_max     : *maximum surface density (g/cm^2)*
!   - sig_out     : *surface density (g/cm^2) at R=R_out*
!   - sig_ref     : *surface density (g/cm^2) at R=R_ref*
!   - udist       : *distance units (cgs)*
!   - umass       : *mass units (cgs)*
!   - utime       : *time units (cgs)*
!
! :Dependencies: allocutils, centreofmass, dim, eos, externalforces,
!   fileutils, grids_for_setup, infile_utils, io, mpidomain, mpiutils,
!   options, part, physcon, random, table_utils, units, vectorutils
!
 use dim,      only:maxvxyzu
 use mpidomain,only:i_belong_i4
 use io,       only:warning,error,fatal
 use mpiutils, only:reduceall_mpi
 use part,     only:igas,labeltype,maxp
 use physcon,  only:c,gg,pi
 use units,    only:umass,udist,utime,unit_angmom
 implicit none
 public :: set_disc,set_incline_or_warp,get_disc_mass,scaled_sigma

 private
 integer, parameter, public :: maxbins = 4096

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine set_disc(id,master,mixture,nparttot,npart,npart_start,rmin,rmax, &
                    rmindust,rmaxdust,phimin,phimax,indexprofile,indexprofiledust, &
                    rc,rcdust,p_index,p_indexdust,q_index,q_indexdust,HoverR,HoverRdust,gamma, &
                    disc_mass,disc_massdust,sig_norm,star_mass,xyz_origin,vxyz_origin, &
                    particle_type,particle_mass,hfact,xyzh,vxyzu,polyk, &
                    position_angle,inclination,ismooth,alpha,rwarp,warp_smoothl, &
                    e0,eindex,phiperi,eccprofile, &
                    bh_spin,bh_spin_angle,rref,enc_mass,r_grid,writefile,ierr,prefix,verbose)
 use io,   only:stdout
 use part, only:maxp,idust,maxtypes
 use centreofmass, only:get_total_angular_momentum
 use allocutils, only:allocate_array
 use grids_for_setup, only: init_grid_sigma,init_grid_ecc,deallocate_sigma,deallocate_ecc
 integer,           intent(in)    :: id,master
 integer, optional, intent(in)    :: nparttot
 integer,           intent(inout) :: npart
 integer, optional, intent(in)    :: npart_start,indexprofile,indexprofiledust
 real,              intent(in)    :: rmin,rmax
 real, optional,    intent(in)    :: rmindust,rmaxdust,p_indexdust,q_indexdust,HoverRdust,disc_massdust
 real, optional,    intent(in)    :: rc,rcdust,rref
 real, optional,    intent(in)    :: phimin,phimax
 real, optional,    intent(inout) :: alpha
 real,              intent(in)    :: p_index,q_index,HoverR,gamma,hfact
 real, optional,    intent(in)    :: disc_mass,star_mass,sig_norm
 real, optional,    intent(in)    :: xyz_origin(3),vxyz_origin(3)
 integer, optional, intent(in)    :: particle_type
 real, optional,    intent(in)    :: position_angle,inclination
 real, optional,    intent(in)    :: rwarp,warp_smoothl,bh_spin,bh_spin_angle
 real, optional,    intent(in)    :: enc_mass(maxbins),r_grid(maxbins)
 real, optional,    intent(in)    :: e0,eindex,phiperi
 integer, optional, intent(in)    :: eccprofile
 logical, optional, intent(in)    :: ismooth,mixture
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,particle_mass
 logical, optional, intent(in)    :: writefile,verbose
 integer, optional, intent(out)   :: ierr
 character(len=20), optional, intent(in) :: prefix
 integer :: itype,npart_tot,npart_start_count,i,npart_set
 integer :: sigmaprofile,sigmaprofiledust
 real    :: Q,G,cs0,cs0dust,clight
 real    :: R_in,R_out,phi_min,phi_max,H_R,H_Rdust,R_indust,R_outdust,p_inddust,q_inddust,R_c,R_c_dust
 real    :: star_m,disc_m,disc_mdust,sigma_norm,sigma_normdust,Q_tmp
 real    :: honH,alphaSS_min,alphaSS_max,rminav,rmaxav,honHmin,honHmax
 real    :: aspin,aspin_angle,posangl,incl,R_warp,H_warp,psimax
 real    :: e_0,e_index,phi_peri
 integer :: ecc_profile
 real    :: xorigini(3),vorigini(3),R_ref,L_tot(3),L_tot_mag
 real    :: enc_m(maxbins),rad(maxbins),enc_m_tmp(maxbins),rad_tmp(maxbins)
 logical :: smooth_surface_density,do_write,do_mixture
 logical :: do_verbose,exponential_taper,exponential_taper_dust
 logical :: exponential_taper_alternative,exponential_taper_dust_alternative
 logical :: use_sigma_file,use_sigmadust_file
 real, allocatable :: ecc_arr(:)
 real, allocatable :: a_arr(:)

 call allocate_array('ecc_arr', ecc_arr, maxp)
 call allocate_array('a_arr', a_arr, maxp)
 ecc_arr(:)=0.
 a_arr(:)=0.
 !
 !--set problem parameters
 !
 H_R            = HoverR       !--HoverR is H/R at R=R_ref (R_ref=R_in unless specified)
 R_in           = rmin
 R_out          = rmax
 if (present(star_mass)) then
    star_m = star_mass
 else
    star_m = 1.d0
 endif
 if (present(npart_start)) then
    npart_start_count = npart_start
 else
    npart_start_count = 1
 endif
 if (rmax < rmin) then
    if (id==master) call error('set_disc','outer radius < inner radius')
    if (present(ierr)) ierr = 2
    return
 endif
 if (present(nparttot)) then
    npart_set = nparttot
 else
    npart_set = npart
 endif
 if (npart_set <= 0) then
    if (id==master) call error('set_disc','nparttot <= 0')
    if (present(ierr)) ierr = 3
    return
 endif
 if (present(xyz_origin)) then
    xorigini = xyz_origin
 else
    xorigini = 0.
 endif
 if (present(vxyz_origin)) then
    vorigini = vxyz_origin
 else
    vorigini = 0.
 endif
 smooth_surface_density = .false.
 if (present(ismooth)) smooth_surface_density = ismooth

 exponential_taper = .false.
 exponential_taper_alternative = .false.
 use_sigma_file = .false.
 if (present(indexprofile)) then
    if (indexprofile==1) exponential_taper = .true.
    if (indexprofile==2) exponential_taper_alternative = .true.
    if (indexprofile==3) use_sigma_file = .true.
    !--loading grids for initialising sigma and eccentricity from file
    if (indexprofile==3) call init_grid_sigma(R_in,R_out)
 endif

 if (present(eccprofile)) then
    if (eccprofile==4) call init_grid_ecc(R_in,R_out)
 endif

 exponential_taper_dust = .false.
 exponential_taper_dust_alternative = .false.
 use_sigmadust_file = .false.
 if (present(indexprofiledust)) then
    if (indexprofiledust==1) exponential_taper_dust = .true.
    if (indexprofiledust==2) exponential_taper_dust_alternative = .true.
    !if (indexprofiledust==3) use_sigmadust_file = .true.
 endif

 aspin = 0.
 if (present(bh_spin)) then
    if (.not. isnan(bh_spin)) aspin = bh_spin
 endif

 aspin_angle = 0.
 if (present(bh_spin_angle)) aspin_angle = bh_spin_angle

 !--reference radius for normalisation of sigma, temperature profiles
 if (present(rref)) then
    R_ref = rref
 else
    R_ref = R_in
 endif
 do_verbose = .true.
 if (present(verbose)) do_verbose = verbose
 phi_min = 0.
 phi_max = 2.*pi
 if (present(phimin)) phi_min = phimin
 if (present(phimax)) phi_max = phimax
 G      = gg*umass*utime**2/(udist**3)
 clight = c*utime/udist
 if (id==master .and. do_verbose) then
    print "(/,a)",' Phantom: general disc setup (see .discparams file for details)'
 endif
 if (present(particle_type)) then
    itype = particle_type
 else
    itype = igas
 endif
 do_mixture = .false.
 q_inddust = 0.  ! avoid compiler warning if not set
 if (present(mixture)) then
    do_mixture = mixture
    if (do_mixture) then
       if (.not.(present(rmindust) .and. present(rmaxdust) .and. present(p_indexdust) &
           .and. present(q_indexdust) .and. present(HoverRdust) .and. present(disc_massdust))) then
          call fatal('set_disc','setup for dusty disc in the mixture is not specified')
          H_Rdust = 0. ! to prevent compiler warning
       else
          if (rmaxdust < rmindust) then
             if (id==master) call error('set_disc','dust outer radius < dust inner radius')
             if (present(ierr)) ierr = 2
             return
          endif
          if (HoverR < HoverRdust .or. q_indexdust < q_index) then
             if (id==master) call error('set_disc','dust scale height > gas scale height')
             if (present(ierr)) ierr = 2
             return
          endif
          R_indust  = rmindust
          R_outdust = rmaxdust
          p_inddust = p_indexdust
          q_inddust = q_indexdust
          H_Rdust   = HoverRdust
       endif
    endif
 endif
 if (id==master .and. do_verbose) then
    if (do_mixture) then
       print "(a,i8,a)",' Setting up disc mixture containing ',npart_set,' '// &
              trim(labeltype(itype))//'/'//trim(labeltype(itype+1))//' particles'
    else
       print "(a,i8,a)",' Setting up disc containing ',npart_set,' '//trim(labeltype(itype))//' particles'
    endif
 endif
 !
 !--set sound speed (cs0 is sound speed at R=1; H_R is at R=R_ref)
 !  and polyk
 !
 cs0 = H_R*sqrt(G*star_m/R_ref)*R_ref**q_index
 polyk = cs0**2
 if (do_mixture) cs0dust = H_Rdust*sqrt(G*star_m/R_ref)*R_ref**q_inddust   !variable to compute the dust scale height
 !
 !--set surface density profile
 !
 !    0 = power law
 !    1 = exponentially tapered power law
 !    2 = smoothed power law
 !    3 = both tapered and smoothed
 !    4 = alternative taper
 !    5 = alternative taper with smoothing
 !    6 = density profile from file
 sigmaprofile = 0
 if (exponential_taper) then
    sigmaprofile = 1
    if (present(rc)) then
       R_c = rc
    else
       R_c = R_out
    endif
 endif
 if (smooth_surface_density) sigmaprofile = 2
 if (smooth_surface_density .and. exponential_taper) sigmaprofile = 3
 if (exponential_taper_alternative) sigmaprofile = 4
 if (exponential_taper_alternative .and. smooth_surface_density) sigmaprofile = 5
 if (use_sigma_file) sigmaprofile = 6
 !--mixture
 if (do_mixture) then
    sigmaprofiledust = 0
    if (exponential_taper_dust) then
       sigmaprofiledust = 1
       if (present(rcdust)) then
          R_c_dust = rcdust
       else
          R_c_dust = R_outdust
       endif
    endif
    if (smooth_surface_density) sigmaprofiledust = 2
    if (smooth_surface_density .and. exponential_taper_dust) sigmaprofiledust = 3
    if (exponential_taper_dust_alternative) sigmaprofiledust = 4
    if (exponential_taper_dust_alternative .and. smooth_surface_density) sigmaprofiledust = 5
    if (use_sigmadust_file) sigmaprofiledust = 6
 endif
 !
 !--disc mass and sigma normalisation
 !
 if (present(r_grid)) then
    rad = r_grid
 else
    do i=1,maxbins
       rad(i) = R_in + (i-1) * (R_out-R_in)/real(maxbins-1)
       !R_in + 0.5*(R_out-R_in)/maxbins
    enddo
 endif

 if (present(sig_norm)) then
    if (present(disc_mass)) then
       call fatal('set_disc','cannot set disc_mass and sig_norm at same time')
    endif
    !--set disc mass from sigma_norm
    sigma_norm = sig_norm
    call get_disc_mass(disc_m,enc_m,rad,Q,sigmaprofile,sigma_norm, &
                       star_m,p_index,q_index,R_in,R_out,R_ref,R_c,H_R)
 elseif (present(disc_mass)) then
    !--compute trial disc mass and Toomre Q
    sigma_norm = 1.d0
    call get_disc_mass(disc_m,enc_m,rad,Q_tmp,sigmaprofile,sigma_norm, &
                       star_m,p_index,q_index,R_in,R_out,R_ref,R_c,H_R)
    sigma_norm = sigma_norm * disc_mass / disc_m
    !--recompute actual disc mass, enc_m, rad, and Toomre Q
    call get_disc_mass(disc_m,enc_m,rad,Q,sigmaprofile,sigma_norm, &
                       star_m,p_index,q_index,R_in,R_out,R_ref,R_c,H_R)
 else
    disc_m = 0.
    sigma_norm = 0.
    call fatal('set_disc','need to set disc mass directly or via sigma normalisation')
 endif
 enc_m = enc_m + star_m
 !print*, 'encm in setdisc', enc_m(1:20)
 !
 !--dust mass
 !
 sigma_normdust = 1.d0
 if (do_mixture) then
    if (present(r_grid)) then
       rad_tmp = r_grid
    else
       do i=1,maxbins
          rad_tmp(i) = R_indust + (i-1) * (R_outdust-R_indust)/real(maxbins-1)
       enddo
    endif
    !--sigma_normdust set from dust disc mass
    call get_disc_mass(disc_mdust,enc_m_tmp,rad_tmp,Q_tmp,sigmaprofiledust, &
                       sigma_normdust,star_m,p_indexdust,q_inddust, &
                       R_indust,R_outdust,R_ref,R_c_dust,H_Rdust)
    sigma_normdust = sigma_normdust*disc_massdust/disc_mdust
    disc_mdust = disc_massdust
 endif
 !
 !--set the particle mass
 !
 if (.not.do_mixture) then
    particle_mass = disc_m / dble(npart_set)
 else
    particle_mass = (disc_m+disc_mdust) / dble(npart_set)
 endif
 !
 !--count particles on this MPI thread
 !
 if (present(nparttot)) then
    npart = 0
    do i=1,npart_set
       if (i_belong_i4(i)) npart = npart + 1
    enddo
 endif
 !
 !--set particle positions and smoothing lengths
 !
 ecc_profile=0
 e_0=0.
 e_index=0.
 phi_peri=0.
 if (present(e0)) e_0=e0
 if (present(eindex)) e_index=eindex
 if (present(phiperi)) phi_peri=phiperi
 if (present(eccprofile)) ecc_profile=eccprofile

 npart_tot = npart_start_count + npart_set - 1
 if (npart_tot > maxp) call fatal('set_disc', &
    'number of particles > array size, use e.g. "./phantomsetup --maxp=10000000"')
 call set_disc_positions(npart_tot,npart_start_count,do_mixture,R_ref,R_in,R_out,&
                         R_indust,R_outdust,phi_min,phi_max,sigma_norm,sigma_normdust,&
                         sigmaprofile,sigmaprofiledust,R_c,R_c_dust,p_index,p_inddust,cs0,cs0dust,&
                         q_index,q_inddust,e_0,e_index,phi_peri,ecc_arr,a_arr,ecc_profile,&
                         star_m,G,particle_mass,hfact,itype,xyzh,honH,do_verbose)

 if (present(inclination)) then
    incl = inclination
 else
    incl = 0.
 endif
 !
 !--override enclosed mass with the correct version
 !  for the case where multiple calls to set_disc are used
 !  e.g. for discs with gas and dust
 !
 if (present(enc_mass)) then
    enc_m = enc_mass
 endif
 !
 !--set particle velocities
 !
 call set_disc_velocities(npart_tot,npart_start_count,itype,G,star_m,aspin,aspin_angle, &
                          clight,cs0,exponential_taper,p_index,q_index,gamma,R_in, &
                          rad,enc_m,smooth_surface_density,xyzh,vxyzu,incl,ecc_arr,a_arr)
 !
 !--inclines and warps
 !
 posangl = 0.
 incl = 0.
 R_warp = 0.
 H_warp = 0.
 psimax = 0.
 if (present(position_angle)) posangl = position_angle
 if (present(rwarp))           R_warp = rwarp
 if (present(warp_smoothl))    H_warp = warp_smoothl
 !
 !--work out h/H in order to set the artificial viscosity parameter to match a chosen alpha_SS
 !
 if (do_verbose) write(*,'(/,1x,"(<h>/H) per particle...: ",f9.4)') honH
 if (smooth_surface_density .and. p_index > 0.) then
    rminav = R_in*((1.+2.*p_index)**2)/(4.*p_index**2)
 else
    rminav = R_in
 endif
 rmaxav = R_out
 call get_honH(xyzh,rminav,rmaxav,honHmin,honHmax,honH,cs0,q_index,star_m,&
               npart_start_count,npart_tot,do_verbose,R_warp)

 if (present(inclination)) then
    !--incline disc at position angle and poss. warp disc
    incl = inclination
    call set_incline_or_warp(xyzh,vxyzu,npart_tot,npart_start_count,posangl,incl,&
                             R_warp,H_warp,psimax)
 endif

#ifdef DISC_VISCOSITY
 !
 !--if disc viscosity is used, set the artificial viscosity parameter
 !  in the input file so as to give the desired alpha_SS
 !
 if (present(alpha)) then
    if (do_verbose) print "(a,g11.4)", ' alphaSS requested = ', alpha
    alpha = alpha/(honH/10.0)
    !--and the min and max alphaSS present
    alphaSS_min = alpha*honHmin/10.
    alphaSS_max = alpha*honHmax/10.
    if (do_verbose) print "(a,g11.4,a)", ' Setting alpha_AV  = ',alpha,' to give alphaSS as requested'
 else
    alphaSS_min = honHmin/10.
    alphaSS_max = honHmax/10.
 endif
#else
 !
 !--if disc viscosity is not used, simply return the range of alphaSS
 !  implied in the disc by the chosen artificial viscosity parameter
 !  see Meru & Bate (2010)
 !
 alphaSS_min = honHmin*(31./525.)
 alphaSS_max = honHmax*(31./525.)
#endif
 !
 !--adjust positions and velocities so the centre of mass is at the origin
 !  also shift particles to new origin if this is not at (0,0,0)
 !
 if (present(phimax)) then
    if (do_verbose) print "(a)",'Setting up disc sector - not adjusting centre of mass'
 else
    call adjust_centre_of_mass(xyzh,vxyzu,particle_mass,npart_start_count,npart_tot,xorigini,vorigini,e_0,ecc_profile)
 endif
 ! Calculate the total angular momentum of the disc only
 call get_total_angular_momentum(xyzh,vxyzu,npart,L_tot)
 L_tot_mag = sqrt(dot_product(L_tot,L_tot))*unit_angmom
 !
 !--print out disc parameters, to file and to the screen
 !
 if (id==master) then
    if (present(writefile)) then
       do_write = writefile
    else
       do_write = .true.
    endif
    if (do_write .and. present(prefix)) then
       if (itype /= igas .and. itype > 0 .and. itype <= maxtypes) then
          open(1,file=trim(prefix)//'-'//trim(labeltype(itype))//'.discparams', &
             status='replace',form='formatted')
       else
          open(1,file=trim(prefix)//'.discparams',status='replace',form='formatted')
       endif
       call write_discinfo(1,R_in,R_out,R_ref,Q,npart,sigmaprofile,R_c,p_index,q_index, &
                           star_m,disc_m,sigma_norm,real(incl*180.0/pi),honH,cs0, &
                           alphaSS_min,alphaSS_max,R_warp,psimax,L_tot_mag,itype)
       close(1)
       if (do_mixture) then
          open(1,file=trim(prefix)//'-'//trim(labeltype(idust))//'.discparams', &
             status='replace',form='formatted')
          call write_discinfo(1,R_indust,R_outdust,R_ref,Q,npart,sigmaprofiledust, &
                              R_c_dust,p_inddust,q_inddust,star_m,disc_massdust, &
                              sigma_normdust,real(incl*180.0/pi),honH,cs0dust, &
                              alphaSS_min,alphaSS_max,R_warp,psimax,L_tot_mag,idust)
          close(1)
       endif
    endif
    !--write disc parameters to screen
    if (do_verbose) then
       call write_discinfo(stdout,R_in,R_out,R_ref,Q,npart,sigmaprofile, &
                           R_c,p_index,q_index,star_m,disc_m,sigma_norm, &
                           real(incl*180.0/pi),honH,cs0,alphaSS_min,alphaSS_max, &
                           R_warp,psimax,L_tot_mag,itype)
    endif
    if (do_mixture) then
       call write_discinfo(stdout,R_indust,R_outdust,R_ref,Q,npart,sigmaprofiledust, &
                           R_c_dust,p_inddust,q_inddust,star_m,disc_massdust, &
                           sigma_normdust,real(incl*180.0/pi),honH,cs0dust, &
                           alphaSS_min,alphaSS_max,R_warp,psimax,L_tot_mag,idust)
    endif
 endif

 if (use_sigma_file) call deallocate_sigma()
 if (ecc_profile==4) call deallocate_ecc()

 if (allocated(ecc_arr)) deallocate(ecc_arr)
 if (allocated(a_arr)) deallocate(a_arr)

 return
end subroutine set_disc

!----------------------------------------------------------------
!
! function to return the sound speed given the radius
!
!----------------------------------------------------------------
pure real function cs_func(cs0,r,q_index)
 real, intent(in) :: cs0,r,q_index

 cs_func = cs0*r**(-q_index)

end function cs_func

!---------------------------------------------------------
!
! set up the particle positions and smoothing length
! using a Monte Carlo approach
!
!---------------------------------------------------------
subroutine set_disc_positions(npart_tot,npart_start_count,do_mixture,R_ref,R_in,R_out,&
                              R_indust,R_outdust,phi_min,phi_max,sigma_norm,sigma_normdust,&
                              sigmaprofile,sigmaprofiledust,R_c,R_c_dust,p_index,p_inddust,cs0,cs0dust,&
                              q_index,q_inddust,e_0,e_index,phi_peri,ecc_arr,a_arr,ecc_profile,&
                              star_m,G,particle_mass,hfact,itype,xyzh,honH,verbose)
 use io,             only:id,master
 use part,           only:set_particle_type
 use random,         only:ran2
 use fileutils, only: load_data_file
 integer, intent(in)    :: npart_start_count,npart_tot
 real,    intent(in)    :: R_ref,R_in,R_out,phi_min,phi_max
 real,    intent(in)    :: sigma_norm,p_index,cs0,q_index,star_m,G,particle_mass,hfact
 real,    intent(in)    :: sigma_normdust,R_indust,R_outdust,R_c,R_c_dust,p_inddust,q_inddust,cs0dust
 real,    intent(in)    :: e_0,e_index,phi_peri
 real,    intent(inout) :: ecc_arr(:),a_arr(:)
 logical, intent(in)    :: do_mixture,verbose
 integer, intent(in)    :: itype,sigmaprofile,sigmaprofiledust,ecc_profile
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(out)   :: honH
 integer :: i,iseed,ninz,n_to_place
 integer :: ipart
 real    :: rand_no,randtest,R,phi,zi,ea,Mmean
 real    :: f,fr_max,fz_max,sigma,cs,omega,fmixt,distr_corr_max,distr_corr_val
 real    :: HH,HHsqrt2,z_min,z_max
 real    :: rhopart,rhoz,hpart
 real    :: xcentreofmass(3)
 real    :: dR,dphi,f_val,sigmamixt,HHdust,HHsqrt2dust,rhozmixt,csdust
 real    :: R_ecc,Rin,Rout,Rindust,Routdust,phi_perirad
 !--seed for random number generator
 iseed = -34598 + (itype - igas)
 honH = 0.
 ninz = 0
 !--n_to_place determines how many particles are initialised for each i in the do loop
 !--n_to_place=2 (assumes symmetry wrt origin), n_to_place=1 sets position individually
 n_to_place=2
 if (abs(e_0)>tiny(e_0) .or. ecc_profile  /=  0) then
    n_to_place=1
 endif

 !--converting phi_peri to radians
 phi_perirad=phi_peri*3.1415/180.

 Rin=R_in
 Rout=R_out
 !--same for the dust
 Rindust=R_indust
 Routdust=R_outdust

 !--set maximum f=R*sigma value
 dR = (Rout-Rin)/real(maxbins-1)
 fr_max = 0.
 dphi= (6.283185)/real(maxbins-1)

 distr_corr_max=0.
 distr_corr_val=0.

 do i=1,maxbins
    R = Rin + (i-1)*dR
    !---------This if cycle needed to discover maximum value of distr_ecc_corr needed for ecc geom.
    phi=0.
    distr_corr_val=distr_ecc_corr(R,phi,R_ref,e_0,e_index,phi_peri,ecc_profile)!*&
    if (e_0 > 1.) then
       call fatal('set_disc','set_disc_positions: e_0>1, set smaller eccentricity')
    elseif (distr_corr_val < 0.) then
       call fatal('set_disc','set_disc_positions: distr_corr<0, choose a shallower eccentricity profile')
    endif
    distr_corr_max=max(distr_corr_max,distr_corr_val)

    f_val = R*sigma_norm*scaled_sigma(R,sigmaprofile,p_index,R_ref,&
                                      Rin,Rout,R_c)*distr_corr_max
    !--distr_corr_max is maximum correction
    !--in distr_ecc_corr(....) for eccentric topology
    if (do_mixture) then
       if (R>=Rindust .and. R<=Routdust) then
          f_val = f_val + R*sigma_normdust*&
                  scaled_sigma(R,sigmaprofiledust,p_inddust,R_ref,&
                               Rindust,Routdust,R_c_dust)*distr_corr_max
       endif
    endif
    fr_max = max(fr_max,f_val)
 enddo

 sigmamixt = 0.
 rhozmixt = 0.
 xcentreofmass = 0.
 ipart = npart_start_count - 1

 !--loop over particles
 do i=npart_start_count,npart_tot,n_to_place
    !--get a random angle between phi_min and phi_max
    rand_no = ran2(iseed)
    phi = 0.! if circular no need of phi! phi_min + (phi_max - phi_min)*ran2(iseed)
    !--now get radius
    f = 0.
    randtest = 1.
    do while (randtest > f)
       R = Rin + (Rout - Rin)*ran2(iseed)
       !--Note that here R is the semi-maj axis, if e0=0. R=a
       ea=ecc_distrib(R,e_0,R_ref,e_index,ecc_profile)

       if ((abs(e_0) > tiny(e_0)) .or. (ecc_profile > 0)) then !-- We generate mean anomalies
          Mmean = phi_min + (phi_max - phi_min)*ran2(iseed)
          !--This is because rejection must occur on the couple (a,phi)
          !--and not only on a.
          !--we convert Mean anomaly to true anomaly, this produces right
          !--azimuthal density
          phi=m_to_f(ea,Mmean)
       else
          phi=phi_min + (phi_max - phi_min)*ran2(iseed)
       endif
       !--distr_ecc_corr is there because variable is "a" but standard density integrates on R
       !--distr_ecc_corr accounts for the change of variable RdR --> R dR/da da
       randtest = fr_max*ran2(iseed)
       f = R*sigma_norm*scaled_sigma(R,sigmaprofile,&
                                     p_index,R_ref,Rin,Rout,R_c)*&
                        distr_ecc_corr(R,phi,R_ref,e_0,e_index,phi_peri,ecc_profile)
       sigma = f/(R*distr_ecc_corr(R,phi,R_ref,e_0,e_index,phi_peri,ecc_profile))
       if (do_mixture) then
          if (R>=Rindust .and. R<=Routdust) then
             fmixt = R*sigma_normdust*scaled_sigma(R,sigmaprofiledust,&
                                                   p_inddust,R_ref,Rindust,&
                                                   Routdust,R_c_dust)*&
                       distr_ecc_corr(R,phi,R_ref,e_0,e_index,phi_peri,ecc_profile)
             f     = f + fmixt
             sigmamixt = fmixt/(R*distr_ecc_corr(R,phi,R_ref,e_0,&
                                                 e_index,phi_peri,ecc_profile))
          endif
       endif
    enddo

    !--set z
    !--first get sound speed at R
    cs = cs_func(cs0,R,q_index)
    !--then pressure scale-height - multiplied by sqrt(2)
    !  for convenience in rhoz calc.
    omega   = sqrt(G*star_m/R**3)
    HH      = cs/omega
    HHsqrt2 = sqrt(2.0d0)*HH
    if (do_mixture) then
       !-- get HH of dust at R
       csdust = cs_func(cs0dust,R,q_inddust)
       HHdust = csdust/omega
       HHsqrt2dust = sqrt(2.0d0)*HHdust
    endif
    z_min = -3.0d0*HHsqrt2
    z_max =  3.0d0*HHsqrt2

    fz_max = sigma/(HHsqrt2*sqrt(pi))
    if (do_mixture) fz_max = fz_max + sigmamixt/(HHsqrt2dust*sqrt(pi))
    f = 0.
    randtest = 1.
    do while(randtest > f)
       zi = z_min + (z_max - z_min)*ran2(iseed)
       randtest = fz_max*ran2(iseed)
       f = sigma*exp(-(zi/(HHsqrt2))**2)/(HHsqrt2*sqrt(pi))
       rhoz = f
       if (do_mixture) then
          fmixt = sigmamixt*exp(-(zi/(HHsqrt2dust))**2)/(HHsqrt2dust*sqrt(pi))
          f    = f + fmixt
          rhozmixt = fmixt
       endif
    enddo

    !--starting estimate of smoothing length for h-rho iterations
    rhopart = rhoz
    if (do_mixture) rhopart = rhopart + rhozmixt
    hpart = hfact*(particle_mass/rhopart)**(1./3.)

    if (i_belong_i4(i)) then
       ipart = ipart + 1
       !--Setting ellipse properties after MC sampling
       ecc_arr(ipart)=ea
       a_arr(ipart)=R
       if (ecc_arr(ipart)>0.99) then
          call fatal('set_disc', 'set_disc_positions: part ',i,' have ecc >1.')
       endif

       R_ecc=R*(1.-ecc_arr(ipart)**2.)/&
               (1.+ecc_arr(ipart)*cos(phi-phi_perirad))

       !--set positions -- move to origin below
       xyzh(1,ipart) = R_ecc*cos(phi)
       xyzh(2,ipart) = R_ecc*sin(phi)
       xyzh(3,ipart) = zi
       xyzh(4,ipart) = hpart
       !--set particle type
       call set_particle_type(ipart,itype)
    endif


    if (i_belong_i4(i+1) .and. i+1 <= npart_tot .and. n_to_place==2 ) then
       ipart = ipart + 1
       !--set positions -- move to origin below

       !--NB this is not redundant as we need to store ecc_arr and a_arr for each particle
       !--for initialising velocities despite n_to_place==2 using same ea and R.
       !--Defining Recc in this if avoid warning of R_ecc used uninitialised
       ecc_arr(ipart)=ea
       a_arr(ipart)=R
       R_ecc=R*(1.-ecc_arr(ipart)**2.)/&
               (1.+ecc_arr(ipart)*cos(phi-phi_perirad))
       !---------

       xyzh(1,ipart) = -R_ecc*cos(phi)
       xyzh(2,ipart) = -R_ecc*sin(phi)
       xyzh(3,ipart) = -zi
       xyzh(4,ipart) = hpart
       !--set particle type
       call set_particle_type(ipart,itype)
    endif

    !--HH is scale height
    if (zi*zi < HH*HH) then
       ninz = ninz + 1
       honH = honH + hpart/HH
    endif
    if (id==master .and.  mod(ipart,max(npart_tot/10,10))==0 .and. verbose) print*,ipart
 enddo

 !--set honH
 honH = honH/real(ninz)

end subroutine set_disc_positions

!----------------------------------------------------------------
!
! set up the particle velocities
!
!----------------------------------------------------------------
! Note that phi_peri must be added if different from 0.
subroutine set_disc_velocities(npart_tot,npart_start_count,itype,G,star_m,aspin, &
                               aspin_angle,clight,cs0,do_sigmapringle,p_index, &
                               q_index,gamma,R_in,rad,enc_m,smooth_sigma,xyzh,vxyzu,inclination,&
                               ecc_arr,a_arr)
 use externalforces, only:iext_einsteinprec
 use options,        only:iexternalforce
 use part,           only:gravity
 use dim,            only:gr
 integer, intent(in)    :: npart_tot,npart_start_count,itype
 real,    intent(in)    :: G,star_m,aspin,aspin_angle,clight,cs0,p_index,q_index
 real,    intent(in)    :: rad(:),enc_m(:),gamma,R_in
 logical, intent(in)    :: do_sigmapringle,smooth_sigma
 real,    intent(in)    :: xyzh(:,:),inclination
 real,    intent(in)    :: ecc_arr(:),a_arr(:)
 real,    intent(inout) :: vxyzu(:,:)
 real :: term,term_pr,term_bh,det,vr,vphi,cs,R,phi,a_smj,ecc
 integer :: i,itable,ipart,ierr
 real :: rg,vkep
 logical :: isecc

 isecc=any((abs(ecc_arr(:)) > tiny(ecc_arr(1))))
 print *
 print "(a)",'Setting up disc velocities'
 if (isecc) then
    print "(a)",'!!!!!!!!! Disc velocities set to be eccentric, neglecting pressure corrections !!!!!!!!!'
 endif

 ierr = 0
 ipart = npart_start_count - 1

 do i=npart_start_count,npart_tot
    if (i_belong_i4(i)) then
       ipart = ipart + 1
       !
       !--set velocities to give centrifugal balance:
       !  v_phi^2 = GM/R - term_pr - 2*vphi*term_bh
       !
       R    = sqrt(xyzh(1,ipart)**2 + xyzh(2,ipart)**2)
       phi  = atan2(xyzh(2,ipart),xyzh(1,ipart))
       ecc  = ecc_arr(i)
       a_smj= a_arr(i)
       !--term is v_phi^2, note that a_smj=R by definition in set_positions if ecc=0.
       term = G*star_m/R
       !
       !--correction for Einstein precession (assumes Rg=1)
       !
       if (iexternalforce==iext_einsteinprec) term = term*(1.0 + 6.0/R)
       !
       !--correction due to self-gravity
       !
       if (gravity) then
          itable = nint((R-rad(1))/(rad(2)-rad(1))) + 1
          term = G*enc_m(itable)/R
       endif
       !
       !--add contribution from pressure gradients
       !  pressure contribution from a powerlaw disc
       !
       select case(itype)
       case(igas)
          cs = cs_func(cs0,R,q_index)
          if (do_sigmapringle) then
             term_pr = 0.
          else
             ! NB: We do NOT correct for the smoothing of the inner disc profile in
             ! the orbital speed (as we did previously), this produces a strong response
             ! which is not desired. Instead we allow a non-zero vr in the inner disc
             term_pr = -cs**2*(1.5+p_index+q_index)
          endif
          if (term + term_pr < 0.) then
             call fatal('set_disc', 'set_disc_velocities: '// &
                'pressure correction causing -ve sqrt while setting velocities')
          endif
       case default
          cs = 0.
          term_pr = 0.
       end select
       !
       !--correction due to Lense-Thirring precession:
       !  Nealon, Nixon & Price correction for Nelson & Papaloizou v x h term
       !  this is Eq. 5.21 in Nealon (2013) multiplied by -R
       !
       term_bh = -2.*aspin*(G*star_m/R)**2/clight**3
       if (aspin_angle > tiny(aspin_angle)) then
          ierr = 1
       endif
       !
       !--now solve quadratic equation for vphi
       !
       det = term_bh**2 + 4.*(term + term_pr)
       Rg   = G*star_m/clight**2
       if (gr) then
          vkep = sqrt(G*star_m/R)
          ! Pure post-Newtonian velocity i.e. no pressure corrections
          vphi = vkep**4/clight**3 * (sqrt(aspin**2 + (R/Rg)**3) - aspin) * cos(inclination)
          vr=0.
       elseif (isecc) then
          !--if eccentric we ignore pressure correction for setup, and use a_smj in vkep
          !--eccentric velocities (Eq. 2.31-2.32 in Murray & Dermott, 1999)
          vkep = sqrt(G*star_m/a_smj)
          vphi=vkep*(1.+ecc*cos(phi))/sqrt(1.-ecc**2)
          vr = vkep*ecc*sin(phi)/sqrt(1.-ecc**2)
       else
          vphi = 0.5*(term_bh + sqrt(det))
          vr=0.
       endif
       !
       !--set velocities -- move to origin below
       !
       vxyzu(1,ipart) = -vphi*sin(phi)+ vr*cos(phi)
       vxyzu(2,ipart) = vphi*cos(phi) + vr*sin(phi)
       vxyzu(3,ipart) = 0.0d0
       !
       !--set thermal energy
       !  utherm generally should not be stored
       !  for an isothermal equation of state
       if (maxvxyzu >= 4) then
          if (itype==igas) then
             if (gamma > 1.) then
                vxyzu(4,ipart) = cs**2/(gamma - 1.)/gamma
             else
                vxyzu(4,ipart) = 1.5*cs**2
             endif
          else
             vxyzu(4,ipart) = 0.
          endif
       endif
    endif
 enddo
 if (ierr /= 0) call warning('set_disc','set_disc_velocities: '// &
    'assuming that the disc and black hole are aligned')

end subroutine set_disc_velocities

!-------------------------------------------------------------
!
! shift the particles so the centre of mass is at the origin
!
!-------------------------------------------------------------
subroutine adjust_centre_of_mass(xyzh,vxyzu,particle_mass,i1,i2,x0,v0,&
                                 e_0,ecc_profile)
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real,    intent(in)    :: particle_mass
 integer, intent(in)    :: i1,i2,ecc_profile
 real,    intent(in)    :: x0(3),v0(3),e_0
 real :: xcentreofmass(3), vcentreofmass(3)
 integer :: i,ipart
 real    :: totmass

 xcentreofmass = 0.
 vcentreofmass = 0.
 totmass       = 0.
 ipart = i1 - 1
 do i=i1,i2
    if (i_belong_i4(i)) then
       ipart = ipart + 1
       if ((abs(e_0) < tiny(e_0)) .and. (ecc_profile  /=  4)) then
          xcentreofmass = xcentreofmass + particle_mass*xyzh(1:3,ipart)
       endif
       vcentreofmass = vcentreofmass + particle_mass*vxyzu(1:3,ipart)
       totmass = totmass + particle_mass
    endif
 enddo

 totmass = reduceall_mpi('+',totmass)

 xcentreofmass = xcentreofmass/totmass
 vcentreofmass = vcentreofmass/totmass

 xcentreofmass = reduceall_mpi('+',xcentreofmass)
 vcentreofmass = reduceall_mpi('+',vcentreofmass)


 ipart = i1 - 1
 do i=i1,i2
    if (i_belong_i4(i)) then
       ipart = ipart + 1
       xyzh(1:3,ipart)  = xyzh(1:3,ipart)  - xcentreofmass + x0
       vxyzu(1:3,ipart) = vxyzu(1:3,ipart) - vcentreofmass + v0
    endif
 enddo

end subroutine adjust_centre_of_mass

!----------------------------------------------------------------
!
! This subroutine is a utility for inclining discs and
! setting up warps
!
!----------------------------------------------------------------
pure subroutine set_incline_or_warp(xyzh,vxyzu,npart_tot,npart_start,posangl,incl,&
                                    Rwarp,Hwarp,psimax)
 use vectorutils, only:rotatevec
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(out)   :: psimax
 integer, intent(in)    :: npart_start,npart_tot
 real,    intent(in)    :: posangl,incl,Rwarp,Hwarp
 real    :: R,inc,k(3),psi
 integer :: i

 !--rotation axis
 k = (/-sin(posangl), cos(posangl), 0./)

 psi = 0.
 do i=npart_start,npart_tot
    R = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    if (R < Rwarp-Hwarp) then
       inc = 0.
    elseif (R < Rwarp+Hwarp) then
       inc = asin(0.5*(1.+sin(pi/(2.*Hwarp)*(R-Rwarp)))*sin(incl))
!       inc = incl*(0.5*tanh((R-Rwarp)/1.) +0.5)
       psi = pi*Rwarp/(4.*Hwarp)*sin(incl)/sqrt(1. - (0.5*sin(incl))**2)
       psimax = max(psimax,psi)
    else
       inc = incl
    endif
    !--rotate position and velocity
    call rotatevec(xyzh(1:3,i),k,inc)
    call rotatevec(vxyzu(1:3,i),k,inc)
 enddo

end subroutine set_incline_or_warp

!------------------------------------
!
!  Compute H/R at a given radius
!
!------------------------------------
pure real function get_HonR(r,cs0,q_index,star_m,G)
 real, intent(in) :: r,cs0,q_index,star_m,G
 real :: omega,cs,HH

 omega = sqrt(G*star_m/r**3)
 cs    = cs_func(cs0,r,q_index)
 HH    = cs/omega
 get_HonR = HH/r

end function get_HonR

!-----------------------------------------------------------------------------
!
!  Print useful information about the disc to the discparams file
!
!-----------------------------------------------------------------------------
subroutine write_discinfo(iunit,R_in,R_out,R_ref,Q,npart,sigmaprofile, &
                          R_c,p_index,q_index,star_m,disc_m,sigma_norm, &
                          inclination,honH,cs0,alphaSS_min,alphaSS_max, &
                          R_warp,psimax,L_tot_mag,itype)
 use eos,          only:gmw
 use infile_utils, only:write_inopt
 use part,         only:igas
 use physcon,      only:kb_on_mh
 use units,        only:unit_velocity
 use grids_for_setup, only: init_grid_sigma,deallocate_sigma

 integer, intent(in) :: iunit,npart,itype,sigmaprofile
 real,    intent(in) :: R_in,R_out,R_ref,Q,p_index,q_index,star_m,disc_m,sigma_norm,L_tot_mag
 real,    intent(in) :: alphaSS_min,alphaSS_max,R_warp,psimax,R_c,inclination,honH,cs0
 integer :: i
 real    :: T0,T_ref,sig,dR,R

 write(iunit,"(/,a)") '# '//trim(labeltype(itype))//' disc parameters - this file is NOT read by setup'
 call write_inopt(R_in,'R_in','inner disc boundary',iunit)
 call write_inopt(R_ref,'R_ref','reference radius',iunit)
 call write_inopt(R_out,'R_out','outer disc boundary',iunit)
 if (R_warp > 0.) call write_inopt(R_warp,'R_warp','position of warp',iunit)
 if (psimax > 0.) call write_inopt(psimax,'psi_max','maximum warp amplitude',iunit)
 call write_inopt(get_HonR(R_in,cs0,q_index,star_m,1.),'H/R_in','disc aspect ratio H/R at R=R_in',iunit)
 call write_inopt(get_HonR(R_ref,cs0,q_index,star_m,1.),'H/R_ref','disc aspect ratio H/R at R=R_ref',iunit)
 call write_inopt(get_HonR(R_out,cs0,q_index,star_m,1.),'H/R_out','disc aspect ratio H/R at R=R_out',iunit)
 if (R_warp > 0.) then
    call write_inopt(get_HonR(R_warp,cs0,q_index,star_m,1.),'H/R_warp','disc aspect ratio H/R at R=R_warp',iunit)
 endif
 sig = sigma_norm*scaled_sigma(R_in,sigmaprofile,p_index,R_ref,R_in,R_out,R_c)
 sig = sig*umass/udist**2
 call write_inopt(sig,'sig_in','surface density (g/cm^2) at R=R_in',iunit)
 sig = sigma_norm*scaled_sigma(R_ref,sigmaprofile,p_index,R_ref,R_in,R_out,R_c)
 sig = sig*umass/udist**2
 call write_inopt(sig,'sig_ref','surface density (g/cm^2) at R=R_ref',iunit)
 sig = sigma_norm*scaled_sigma(R_out,sigmaprofile,p_index,R_ref,R_in,R_out,R_c)
 sig = sig*umass/udist**2
 call write_inopt(sig,'sig_out','surface density (g/cm^2) at R=R_out',iunit)
 dR = (R_out-R_in)/real(maxbins-1)
 sig = 0.
 do i=1,maxbins
    R = R_in + (i-1)*dR
    sig = max(sig,sigma_norm*scaled_sigma(R,sigmaprofile,p_index,R_ref,R_in,R_out,R_c))
 enddo
 sig = sig*umass/udist**2
 call write_inopt(sig,'sig_max','maximum surface density (g/cm^2)',iunit)
 call write_inopt(Q,'Qmin','minimum Toomre Q parameter',iunit)
 call write_inopt(npart,'n','number of particles in the disc',iunit)
 call write_inopt(p_index,'p_index','power law index of surface density profile',iunit)
 if (sigmaprofile==1 .or. sigmaprofile==3) call write_inopt(R_c,'R_c','characteristic radius of the exponential taper',iunit)
 call write_inopt(q_index,'q_index','power law index of sound speed profile',iunit)
 call write_inopt(star_m,'M_star','mass of central star',iunit)
 call write_inopt(disc_m,'M_disc','disc mass',iunit)
 call write_inopt(disc_m/star_m,'M_disc/M_star','relative disc mass',iunit)
 if (itype == igas) call write_inopt(cs0,'cs0','sound speed at R=1',iunit)

 if (itype == igas) then
    T0 = (cs_func(cs0,R_in,q_index)*unit_velocity)**2*gmw/kb_on_mh
    call write_inopt(T0,'T_in','temperature (K) at R=R_in',iunit)
    T_ref = (cs_func(cs0,R_ref,q_index)*unit_velocity)**2*gmw/kb_on_mh
    call write_inopt(T_ref,'T_ref','temperature (K) at R=R_ref',iunit)
    T0 = (cs_func(cs0,R_out,q_index)*unit_velocity)**2*gmw/kb_on_mh
    call write_inopt(T0,'T_out','temperature (K) at R=R_out',iunit)
 else
    T_ref = 0
 endif

 call write_inopt(inclination,'inc.deg','disc inclination in degrees',iunit)
 call write_inopt(honH,'<h/H>','approx. mean smoothing length over disc scale height',iunit)
 call write_inopt(alphaSS_min,'alphaSS_min','minimum Shakura-Sunyaev alpha viscosity in disc',iunit)
 call write_inopt(alphaSS_max,'alphaSS_max','maximum Shakura-Sunyaev alpha viscosity in disc',iunit)
 call write_inopt(udist,'udist','distance units (cgs)',iunit)
 call write_inopt(umass,'umass','mass units (cgs)',iunit)
 call write_inopt(utime,'utime','time units (cgs)',iunit)
 call write_inopt(gg/(udist**3/(utime**2*umass)),'G','in code units',iunit)
 call write_inopt(c/(udist/utime),'c','in code units',iunit)
 write(iunit,"(a)")

 !--print some of these diagnostics in more useful form
 if (itype == igas) then
    if (T_ref < 1.0d3) then
       write(iunit,"(a,f5.1,a,f5.1,a,f5.2,a,/)")  '# Temperature profile  = ',T_ref,'K (R/',R_ref,')^(',-2.*q_index,')'
    else
       write(iunit,"(a,es9.2,a,f5.1,a,f5.2,a,/)") '# Temperature profile  = ',T_ref,'K (R/',R_ref,')^(',-2.*q_index,')'
    endif
 endif
 if (sigmaprofile==0) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,')'
 elseif (sigmaprofile==1) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,f5.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,') exp[-(R/',R_c,')^(2-',p_index,')]'
 elseif (sigmaprofile==2) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,') (1 - sqrt(',R_in,'/R))'
 elseif (sigmaprofile==3) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,f5.1,a,f4.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,') exp[-(R/',R_c,')^(2-',p_index,')] (1 - sqrt(',R_in,'/R))'
 elseif (sigmaprofile==4) then
    write(iunit, "(a,es9.2,a,f5.1,a,f4.1,a,f5.1,a,f4.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,') (1-exp(R-',R_out,'))'
 elseif (sigmaprofile==5) then
    write(iunit, "(a,es9.2,a,f5.1,a,f4.1,a,f5.1,a,f4.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,') (1-exp(R-',R_out,')) (1 - sqrt(',R_in,'/R))'
 endif
 write(iunit,"(a,es9.2,a,/)") '# Disc total angular momentum = ',L_tot_mag,' g*cm^2/sec'

 return
end subroutine write_discinfo

!-----------------------------------------------------------------------------
!
!  Routine to compute ratio of smoothing length to disc scale height
!  Used to determine alpha viscosity in the disc when using the
!  artificial viscosity to represent Shakura-Sunyaev viscosity
!
!  Input:
!   xyzh - positions and smoothing lengths
!   rminav, rmaxav - rmin and rmax to perform calculation
!   cs0 - sound speed at R=1
!   q_index - power-law index of temperature profile
!   M_star - mass of star
!   i1,i2 - range of particle indices to use
!   rwarp - location of warp (optional)
!
!  Output:
!   honHmin, honH, honHmax - max mean and min ratio of h/H within range of R
!-----------------------------------------------------------------------------
subroutine get_honH(xyzh,rminav,rmaxav,honHmin,honHmax,honH,cs0,q_index,M_star,i1,i2,verbose,rwarp)
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: honHmax,honHmin,honH
 integer, intent(in)  :: i1,i2
 real,    intent(in)  :: cs0,q_index,M_star,rminav,rmaxav,rwarp
 logical, intent(in)  :: verbose

 integer :: i,ii,iwarp
 real :: G,rmin,rmax,dr,ri
 real :: rad(maxbins),ninbin(maxbins),h_smooth(maxbins),cs(maxbins),H(maxbins),omega(maxbins)
 integer :: ipart

 G = 1.0

 !--setup rmin and rmax for the analysis
 rmin = rminav
 rmax = rmaxav

 !--set up the radius array
 dr = (rmax-rmin)/real(maxbins-1)
 iwarp = 0
 do i=1,maxbins
    rad(i)=rmin + real(i-1)*dr
    if (rwarp > tiny(rwarp)) then
       if (rad(i) > rwarp) iwarp = i - 1
    endif
 enddo

 !--initialise arrays to zero
 ninbin = 0
 h_smooth = 0.0

 !--and thus the sound speed array
 do i=1,maxbins
    cs(i) = cs_func(cs0,rad(i),q_index)
    omega(i) = sqrt(G*M_star/rad(i)**3)
 enddo

 !--and thus the disc scale height
 do i=1,maxbins
    H(i) = cs(i)/omega(i)
 enddo

 !--loop over particles putting properties into the correct bin
 ipart = i1 - 1
 do i=i1,i2
    if (i_belong_i4(i)) then
       ipart = ipart + 1
       if (xyzh(4,ipart) > tiny(xyzh)) then ! IF ACTIVE
          ri = sqrt(dot_product(xyzh(1:3,ipart),xyzh(1:3,ipart)))
          ii = int((ri-rad(1))/dr + 1)

          if (ii > maxbins) cycle
          if (ii < 1) cycle

          !--ignoring the large smoothing length particles far from the mid-plane
          if (xyzh(3,ipart)**2 < 4.*H(ii)*H(ii)) then
             h_smooth(ii) = h_smooth(ii) + xyzh(4,ipart)
             ninbin(ii) = ninbin(ii) + 1
          endif
       endif
    endif
 enddo

 h_smooth = reduceall_mpi('+', h_smooth)
 ninbin = reduceall_mpi('+', ninbin)

 !--average h_smooth
 do i=1,maxbins
    if (ninbin(i) > 0) h_smooth(i) = h_smooth(i)/ninbin(i)
 enddo

 !--now loop over rings to calculate required quantities
 do i=1,maxbins
    if (H(i) > 0.) then
       h_smooth(i) = h_smooth(i)/H(i)
    endif
 enddo

 !--print out an average value for <h>/H and thus \alpha_SS/\alpha_AV
 honH = 0.0
 honHmin = minval(h_smooth)
 honHmax = maxval(h_smooth)
 do i=1,maxbins
    honH = honH + h_smooth(i)
 enddo
 honH = honH/real(maxbins)

 if (verbose) then
    write(*,'(1x,"(<h>/H) mean...........: ",f9.4)') honH
    write(*,'(1x,"(<h>/H) minimum........: ",f9.4)') honHmin
    write(*,'(1x,"(<h>/H) maximum........: ",f9.4)') honHmax
    if (rwarp > tiny(rwarp) .and. iwarp > 0) then
       write(*,'(1x,"(<h>/H) at R_warp......: ",f9.4)') h_smooth(iwarp)
    endif
    write(*,'(1x,"e.g. for alpha_SS = 0.1,  use alpha_AV = ",f10.5)') 0.1/(honH/10.0)
    write(*,'(1x,"   using alpha_AV = 1.0 gives alpha_SS ~ ",f8.4,"->",f8.4)') honHmin/10.0,honHmax/10.0
 endif

end subroutine get_honH

!------------------------------------------------------------------------
!
! returns surface density (sigma) at any value of R for a given profile
! with an arbitrary scaling
!
! sigmaprofile options:
!
!    0) power law
!         sigma ~ (R/Rref)^-p
!
!    1) exponentially tapered power law
!         sigma ~ (R/Rref)^-p * exp[-(R/Rc)^(2-p)]
!
!    2) smoothed power law
!         sigma ~ (R/Rref)^-p * (1 - sqrt(Rin/R))
!
!    3) both tapered and smoothed
!         sigma ~ (R/Rref)^-p * exp[-(R/Rc)^(2-p)] * (1 - sqrt(Rin/R))
!
!    4) alternative taper
!         sigma ~ (R/Rref)^-p * [1 - exp(R-Rout)]
!
!    5) alternative taper with smoothing
!         sigma ~ (R/Rref)^-p * [1 - exp(R-Rout)] * (1 - sqrt(Rin/R))
!
!    TODO: accept any surface density profile from function or file
!
!------------------------------------------------------------------------
function scaled_sigma(R,sigmaprofile,pindex,R_ref,R_in,R_out,R_c) result(sigma)
 use table_utils, only: interpolate_1d
 use grids_for_setup, only: datasigma,sigma_initialised,dsigmadx
 real,    intent(in)  :: R,R_ref,pindex
 real,    intent(in)  :: R_in,R_out,R_c
 integer, intent(in)  :: sigmaprofile
 real :: sigma

 sigma=0.

 select case (sigmaprofile)
 case (0)
    sigma = (R/R_ref)**(-pindex)
 case (1)
    sigma = (R/R_ref)**(-pindex)*exp(-(R/R_c)**(2.-pindex))
 case (2)
    sigma = (R/R_ref)**(-pindex)*(1.-sqrt(R_in/R))
 case (3)
    sigma = (R/R_ref)**(-pindex)*exp(-(R/R_c)**(2.-pindex))*(1.-sqrt(R_in/R))
 case (4)
    sigma = (R/R_ref)**(-pindex)*(1-exp(R-R_out))
 case (5)
    sigma = (R/R_ref)**(-pindex)*(1-exp(R-R_out))*(1-sqrt(R_in/R))
 case (6)
    if (sigma_initialised) then
       sigma = interpolate_1d(R,datasigma(:,1),datasigma(:,2),dsigmadx)
    else
       call fatal('set_disc', 'sigma grid not initialised, something went wrong')
    endif
 case default
    call error('set_disc','unavailable sigmaprofile; surface density is set to zero')
    sigma = 0.
 end select

end function scaled_sigma
!-------------------------------

function ecc_distrib(a,e_0,R_ref,e_index,ecc_profile) result(eccval)
 use table_utils, only: interpolate_1d
 use grids_for_setup, only: dataecc,ecc_initialised,deda
 real, intent(in) :: a,e_0,R_ref,e_index
 integer, intent(in) :: ecc_profile
 real :: eccval

 eccval=0.

 select case (ecc_profile)
 case(0)
    eccval=e_0
 case(1)
    eccval=e_0*(a/R_ref)**(-e_index)
 case(4)
    if (ecc_initialised) then
       eccval=interpolate_1d(a,dataecc(:,1),dataecc(:,2),deda)
    else
       call fatal('set_disc', 'ecc grid not initialised, something went wrong')
    endif
 case default
    call error('set_disc','unavailable eccentricity profile, eccentricity is set to zero')
    eccval = 0.
 end select

end function ecc_distrib

function deda_distrib(a,e_0,R_ref,e_index,ecc_profile) result(dedaval)
 use table_utils, only: interpolate_1d
 use grids_for_setup, only: ecc_initialised,dataecc,deda,ddeda
 real, intent(in) :: a,e_0,R_ref,e_index
 integer, intent(in) :: ecc_profile
 real :: dedaval,ea

 dedaval=0.

 select case (ecc_profile)
 case(0)
    dedaval=0.
 case(1)
    ea=e_0*(a/R_ref)**(-e_index)
    dedaval=-e_index*ea/a
 case(4)
    if (ecc_initialised) then
       dedaval=interpolate_1d(a,dataecc(:,1),deda,ddeda)
    else
       call fatal('set_disc', 'ecc grid not initialised, something went wrong')
    endif
 case default
    call error('set_disc','unavailable eccentricity profile, eccentricity is set to zero')
    dedaval = 0.
 end select

end function deda_distrib

!--This function corrects the distribution to account for eccentricity when
!--sampling a and uniform mean anomaly
function distr_ecc_corr(a,phi,R_ref,e_0,e_index,phi_peri,ecc_profile) result(distr)
 real,     intent(in) :: a,phi,R_ref,e_0,e_index,phi_peri
 integer,  intent(in) :: ecc_profile
 real :: distr,ea,deda

 ea = ecc_distrib(a,e_0,R_ref,e_index,ecc_profile) !e_0*(a/R_ref)**(-e_index)
 deda = deda_distrib(a,e_0,R_ref,e_index,ecc_profile)

 distr = 2*pi*(sqrt(1-ea**2)-(a*ea*deda)/2/sqrt(1-ea**2))
 !--distr=1 for e_0=0.

end function distr_ecc_corr

!------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------
subroutine get_disc_mass(disc_m,enc_m,rad,toomre_min,sigmaprofile,sigma_norm, &
                         star_m,pindex,qindex,R_in,R_out,R_ref,R_c,H_R)
 real,           intent(in)  :: sigma_norm,star_m,pindex,qindex,R_in,R_out,R_ref,H_R
 real, optional, intent(in)  :: R_c
 integer,        intent(in)  :: sigmaprofile
 real,           intent(in)  :: rad(:)
 real,           intent(out) :: disc_m,enc_m(:),toomre_min

 real    :: dr,dM,R,sigma,cs0,cs,kappa,G
 integer :: i

 G = gg*umass*utime**2/(udist**3)
 cs0 = H_R*sqrt(G*star_m/R_ref)*R_ref**qindex
 enc_m = 0.
 toomre_min = huge(toomre_min)
 disc_m = 0.
 dR = rad(2)-rad(1)
 do i=1,maxbins
    R = rad(i)
    sigma = sigma_norm * scaled_sigma(R,sigmaprofile,pindex,R_ref,R_in,R_out,R_c)
    !--disc mass
    dM = 2.*pi*R*sigma*dR
    disc_m = disc_m + dM
    enc_m(i) = disc_m
    !--Toomre Q
    cs = cs_func(cs0,R,qindex)
    kappa = sqrt(G*star_m/R**3)
    if (sigma > epsilon(sigma)) then
       toomre_min = min(toomre_min,real(cs*kappa/(pi*G*sigma)))
    endif
 enddo

end subroutine get_disc_mass

function m_to_f(ecc,M) result(F)
 integer          :: i
 real             :: E,A
 real             :: F
 real, intent(in)  :: ecc,M
 !--First find eccentric anomaly
 F=0.

 if (ecc < 1.) then
    if (ecc < 0.8) then
       E=M
    else
       E=pi
    endif
    A = E - ecc*sin(E) - M;
    do i=0,200
       E = E - A/(1.-ecc*cos(E));
       A = E - ecc*sin(E) - M;
       if (abs(A) < 1.E-16) then
          exit
       endif
    enddo
    !--then convert to true anomaly
    F = 2.*atan(sqrt((1.+ecc)/(1.-ecc))*tan(0.5*E))
    F=mod(2*pi + mod(f, 2*pi), 2*pi)
 else
    call fatal('set_disc', 'm_to_f: some particles have ecc >1.')
 endif

end function m_to_f

end module setdisc
