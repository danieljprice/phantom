!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setdisc
!
!  DESCRIPTION:
!  This module contains utility routines for accretion disc setups
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    G           -- in code units
!    M_disc      -- disc mass
!    M_star      -- mass of central star
!    Qmin        -- minimum Toomre Q parameter
!    R_c         -- characteristic radius of the exponential taper
!    R_in        -- inner disc boundary
!    R_out       -- outer disc boundary
!    R_ref       -- reference radius
!    R_warp      -- position of warp
!    T_in        -- temperature (K) at R=R_in
!    T_out       -- temperature (K) at R=R_out
!    T_ref       -- temperature (K) at R=R_ref
!    alphaSS_max -- maximum Shakura-Sunyaev alpha viscosity in disc
!    alphaSS_min -- minimum Shakura-Sunyaev alpha viscosity in disc
!    c           -- in code units
!    cs0         -- sound speed at R=1
!    n           -- number of particles in the disc
!    p_index     -- power law index of surface density profile
!    psi_max     -- maximum warp amplitude
!    q_index     -- power law index of sound speed profile
!    sig_in      -- surface density (g/cm^2) at R=R_in
!    sig_max     -- maximum surface density (g/cm^2)
!    sig_out     -- surface density (g/cm^2) at R=R_out
!    sig_ref     -- surface density (g/cm^2) at R=R_ref
!    udist       -- distance units (cgs)
!    umass       -- mass units (cgs)
!    utime       -- time units (cgs)
!
!  DEPENDENCIES: dim, domain, eos, externalforces, infile_utils, io,
!    mpiutils, options, part, physcon, random, units, vectorutils
!+
!--------------------------------------------------------------------------
module setdisc
 implicit none
 public :: set_disc,set_warp,scaled_sigma,scaled_discmass

 private

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine set_disc(id,master,mixture,nparttot,npart,npart_start,rmin,rmax,rmindust,rmaxdust,phimin,phimax,&
                    indexprofile,indexprofiledust,rc,rcdust,p_index,p_indexdust,q_index,HoverR,gamma,&
                    disc_Q,disc_mass,disc_massdust,sig_norm,star_mass,xyz_origin,vxyz_origin,&
                    particle_type,particle_mass,hfact,xyzh,vxyzu,polyk,position_angle,inclination,sininclination,&
                    twist,ismooth,alpha,isink,rwarp,warp_smoothl,bh_spin,rref,writefile,ierr,prefix,verbose)
 use domain,  only:i_belong
 use eos,     only:qfacdisc
 use io,      only:fatal,comment,warning,stdout,error
 use options, only:ieos
 use part,    only:maxp,igas,idust,set_particle_type,labeltype,gravity,maxtypes
 use physcon, only:pi,gg,c
 use units,   only:umass,udist,utime
 integer,                     intent(in)    :: id,master
 integer, optional,           intent(in)    :: nparttot
 integer,                     intent(inout) :: npart
 integer, optional,           intent(in)    :: npart_start,isink,indexprofile,indexprofiledust
 real,                        intent(in)    :: rmin,rmax
 real, optional,              intent(in)    :: rmindust,rmaxdust,p_indexdust,disc_massdust
 real, optional,              intent(in)    :: rc,rcdust
 real, optional,              intent(in)    :: phimin,phimax
 real, optional,              intent(inout) :: alpha
 real,                        intent(in)    :: p_index,q_index,HoverR,gamma,hfact
 real, optional,              intent(in)    :: disc_Q,disc_mass,star_mass,sig_norm
 real, optional,              intent(in)    :: xyz_origin(3),vxyz_origin(3)
 integer, optional,           intent(in)    :: particle_type
 real, optional,              intent(in)    :: position_angle,inclination,sininclination
 real, optional,              intent(in)    :: rwarp,warp_smoothl,bh_spin,rref
 logical, optional,           intent(in)    :: twist,ismooth,mixture
 real,                        intent(out)   :: xyzh(:,:)
 real,                        intent(out)   :: vxyzu(:,:)
 real,                        intent(out)   :: polyk,particle_mass
 logical, optional,           intent(in)    :: writefile,verbose
 integer, optional,           intent(out)   :: ierr
 character(len=20), optional, intent(in)    :: prefix
 integer, parameter :: maxbins = 256
 integer :: maxvxyzu,itype,npart_tot,npart_start_count,ierror,i,npart_set
 integer :: sigmaprofile,sigmaprofiledust
 real    :: Q,G,cs0,clight
 real    :: R_in,R_out,phi_min,phi_max,H_R,R_indust,R_outdust,p_inddust,R_c,R_c_dust
 real    :: star_m,disc_m,disc_mdust,rminav,rmaxav,honHmin,honHmax
 real    :: honH,alphaSS_min,alphaSS_max
 real    :: xinclination,rwarpi,hwarp,rsi,rso,psimax
 real    :: aspin,posangl,incl
 real    :: xorigini(3),vorigini(3),R_ref
 real    :: enc_m(maxbins),rad(maxbins),R,dR,dM,sigma
 logical :: do_twist,smooth_surface_density,do_write,do_mixture
 logical :: do_verbose,exponential_taper,exponential_taper_dust

 !
 !--set problem parameters
 !
 H_R            = HoverR       !--HoverR is H/R at R=R_ref (R_ref=R_in unless specified)
 R_in           = rmin
 R_out          = rmax
 R_indust       = rmin
 R_outdust      = rmax
 p_inddust      = p_index
 if (present(rmindust) .and. present(rmaxdust) .and. present(p_indexdust)) then
    R_indust    = rmindust
    R_outdust   = rmaxdust
    p_inddust   = p_indexdust
 endif
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
 if (present(disc_Q)) then
    if (present(disc_mass)) then
       if (id==master) call error('set_disc','cannot specify both disc mass and Toomre Q parameter')
       if (present(ierr)) ierr = 1
       return
    endif
    Q = disc_Q
 else
    Q = 168.d0
 endif
 if (rmax < rmin) then
    if (id==master) call error('set_disc','outer radius < inner radius')
    if (present(ierr)) ierr = 2
    return
 endif
 if (present(rmindust) .and. present(rmaxdust)) then
    if (rmaxdust < rmindust) then
       if (id==master) call error('set_disc','dust outer radius < dust inner radius')
       if (present(ierr)) ierr = 2
       return
    endif
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
    xorigini(:) = xyz_origin(:)
 else
    xorigini(:) = 0.
 endif
 if (present(vxyz_origin)) then
    vorigini(:) = vxyz_origin(:)
 else
    vorigini(:) = 0.
 endif
 if (present(inclination)) then
    xinclination = inclination
 else
    xinclination = 0.
 endif
 if (present(ismooth)) then
    smooth_surface_density = ismooth
 else
    smooth_surface_density = .false.
 endif
 if (present(indexprofile) .and. indexprofile==1) then
    exponential_taper = .true.
 else
    exponential_taper = .false.
 endif
 if (present(indexprofiledust) .and. indexprofiledust==1) then
    exponential_taper_dust = .true.
 else
    exponential_taper_dust = .false.
 endif
 if (present(bh_spin)) then
    aspin = bh_spin
 else
    aspin = 0.
 endif
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
    print*,' Phantom: general disc setup (see .discparams file for details)'
 endif
 if (present(particle_type)) then
    itype = particle_type
 else
    itype = igas
 endif
 do_mixture = .false.
 if (present(mixture)) then
    do_mixture = mixture
    if (do_mixture) then
       if (.not.(present(rmindust) .and. present(rmaxdust) .and. present(p_indexdust) .and. present(disc_massdust))) then
          call fatal('set_disc','setup for dusty disc in the mixture is not specified')
       endif
       disc_mdust = disc_massdust
    endif
 endif
 if (id==master .and. do_verbose) then
    if (do_mixture) then
       print*,' Setting up disc mixture containing ',&
                 npart_set,' '//trim(labeltype(itype))//'/'//trim(labeltype(itype+1))//' particles'
    else
       print*,' Setting up disc containing ',npart_set,' '//trim(labeltype(itype))//' particles'
    endif
 endif
 !
 !--store q_index in eos module, then gets saved to dump header
 !
 maxvxyzu = size(vxyzu(:,1))
 qfacdisc = q_index
 !
 !--use isothermal eos (ieos=3) if qfacdisc is set (and now deals with ieos=6 as well)
 !
 if (maxvxyzu < 4 .and. q_index > 0.) then
    if (present(isink)) then
       if (isink==0) then
          call comment('set_disc','setting ieos=3 in input options based on q_index setting and isink=0')
          ieos = 3
       else
          call comment('set_disc','setting ieos=6 in input options based on q_index setting and isink')
          write(*,'("isink = ",i1)') isink
          ieos = 6
       endif
    else
       call comment('set_disc','setting ieos=3 in input options based on q_index setting')
       ieos = 3
    endif
 endif
 !
 !--set sound speed (cs0 is sound speed at R=1; H_R is at R=R_ref)
 !
 cs0 = H_R*sqrt(G*star_m/R_ref)*R_ref**q_index
 !
 !--set surface density profile
 !
 !    0 = power law
 !    1 = exponentially tapered power law
 !    2 = smoothed power law
 !    3 = both tapered and smoothed
 !
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
 endif
 !
 !--get disc mass from sig_norm
 !
 if (present(sig_norm)) then
    if (present(disc_mass)) then
       call fatal('set_disc','cannot set disc_mass and sig_norm at same time')
    endif
    disc_m = sig_norm*scaled_discmass(sigmaprofile,p_index,R_in,R_out,R_ref,R_c)
 else
    disc_m = disc_mass
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
    do i = 1,npart_set
       if (i_belong(i)) npart = npart + 1
    enddo
 endif
 !
 !--set particle positions and smoothing lengths
 !
 npart_tot = npart_start_count + npart_set - 1
 if (npart_tot > maxp) call fatal('set_disc','number of particles exceeds array dimensions',var='n',ival=npart_tot)
 call set_disc_positions(npart_tot,npart_start_count,do_mixture,R_ref,R_in,R_out,&
                         R_indust,R_outdust,phi_min,phi_max,disc_m,disc_mdust,&
                         sigmaprofile,sigmaprofiledust,R_c,R_c_dust,p_index,p_inddust,cs0,&
                         q_index,star_m,G,particle_mass,hfact,itype,xyzh,honH,do_verbose)
 !
 !--set particle velocities
 !
 enc_m(:) = 0.
 dR = (R_out-R_in)/real(maxbins-1)
 do i=1,maxbins
    R = R_in + (i-1)*dR
    sigma = scaled_sigma(R,sigmaprofile,p_index,R_ref,R_in,R_c)
    dM    = 2.*pi*R*sigma*dR
    if (i>1) then
       enc_m(i) = enc_m(i-1) + dM
    else
       enc_m(i) = dM
    endif
    rad(i) = R + 0.5*dR
 enddo
 enc_m(:) = enc_m(:) * disc_m / scaled_discmass(sigmaprofile,p_index,R_in,R_out,R_ref,R_c)
 enc_m(:) = enc_m(:) + star_m
 call set_disc_velocities(npart_tot,npart_start_count,itype,G,star_m,aspin,&
                          clight,cs0,exponential_taper,p_index,q_index,gamma,R_in, &
                          maxbins,rad,enc_m,smooth_surface_density,xyzh,vxyzu,ierror)
 if (ierror==1) then
    call fatal('set_disc','error: pressure correction causing -ve sqrt while setting velocities')
 elseif (ierror /= 0) then
    call fatal('set_disc','error setting velocities')
 endif
 if (present(bh_spin) .and. aspin > 1.d-10 .and. id==master) then
    print*,'Orbital velocities corrected for black hole spin, a = ',aspin
    print*,'ASSUMING THAT THE DISC AND BLACK HOLE ARE ALIGNED'
 endif
 !
 !--set polyk
 !
 polyk = cs0**2
 !
 !--work out h/H in order to set the artificial viscosity parameter to match a chosen alpha_SS
 !
 if (do_verbose) write(*,'(/,1x,"Actual <h>/H...per particle... ",f9.4)') honH
 if (smooth_surface_density .and. p_index > 0.) then
    rminav = R_in*((1.+2.*p_index)**2)/(4.*p_index**2)
 else
    rminav = R_in
 endif
 rmaxav = R_out
 if (present(rwarp)) then
    call get_honH(xyzh,rminav,rmaxav,honHmin,honHmax,honH,npart_set,cs0,q_index,star_m,R_in,&
                  npart_start_count,npart_tot,do_verbose,rwarp)
 else
    call get_honH(xyzh,rminav,rmaxav,honHmin,honHmax,honH,npart_set,cs0,q_index,star_m,R_in,&
                  npart_start_count,npart_tot,do_verbose)
 endif

#ifdef DISC_VISCOSITY
 !
 !--if disc viscosity is used, set the artificial viscosity parameter
 !  in the input file so as to give the desired alpha_SS
 !
 if (present(alpha)) then
    if (do_verbose) print*, 'alphaSS requested = ', alpha
    alpha = alpha/(honH/10.0)
    !--and the min and max alphaSS present
    alphaSS_min = alpha*honHmin/10.
    alphaSS_max = alpha*honHmax/10.
    if (do_verbose) print*, 'Setting alpha_AV = ',alpha,' to give alphaSS as requested'
 else
    alphaSS_min = honHmin/10.
    alphaSS_max = honHmax/10.
 endif
#else
 !
 !--if disc viscosity is not used, simply return the range of alphaSS
 !  implied in the disc by the chosen artificial viscosity parameter
 !
 alphaSS_min = honHmin*(31./525.)
 alphaSS_max = honHmax*(31./525.)
#endif

 !
 !--add a warp/twist to the disc
 !
 rwarpi = 0.
 do_twist = .false.
 if (present(twist)) then
    if (twist) do_twist = .true.
 endif
 !--also twist disc if rwarp is given as an argument
 !  regardless of whether or not twist=.true.
 if (present(rwarp)) do_twist = .true.
 if (do_twist) then
    if (present(rwarp)) then
       rwarpi = rwarp
    else
       rwarpi = 0.5*(R_in+R_out)
    endif
    if (present(warp_smoothl)) then
       hwarp = warp_smoothl
    else
       hwarp = 1.5
    endif
    rsi = rwarpi - hwarp
    rso = rwarpi + hwarp
 else
    rwarpi = 0.
    hwarp  = 0.
    rsi    = 0.
    rso    = 0.
 endif

 !
 !--incline and warps
 !
 if (present(position_angle) .and. present(inclination)) then
    !--incline disc at position angle
    posangl = position_angle
    incl    = inclination
    call rotate_disc(xyzh,vxyzu,npart_tot,npart_start_count,posangl,incl)
 else
    !--call setwarp to calculate the warp
    call set_warp(npart_tot,npart_start_count,&
                  xyzh,vxyzu,inclination,sininclination,&
                  rwarpi,psimax,rsi,rso,do_twist)
 endif

 !
 !--adjust positions and velocities so the centre of mass is at the origin
 !  also shift particles to new origin if this is not at (0,0,0)
 !
 if (present(phimax)) then
    print "(a)",'Setting up disc sector - not adjusting centre of mass'
 else
    call adjust_centre_of_mass(xyzh,vxyzu,particle_mass,npart_start_count,npart_tot,xorigini,vorigini)
 endif

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
          open(1,file=trim(prefix)//'-'//trim(labeltype(itype))//'.discparams',status='replace',form='formatted')
       else
          open(1,file=trim(prefix)//'.discparams',status='replace',form='formatted')
       endif
       call write_discinfo(1,R_in,R_out,R_ref,Q,npart,sigmaprofile,&
                           R_c,p_index,q_index,star_m,disc_m,real(xinclination*180.0/pi),honH,cs0,&
                           alphaSS_min,alphaSS_max,rwarpi,psimax,itype)
       close(1)
       if (do_mixture) then
          open(1,file=trim(prefix)//'-'//trim(labeltype(idust))//'.discparams',status='replace',form='formatted')
          call write_discinfo(1,R_indust,R_outdust,R_ref,Q,npart,sigmaprofiledust,&
                              R_c_dust,p_inddust,q_index,star_m,disc_massdust,real(xinclination*180.0/pi),honH,&
                              cs0,alphaSS_min,alphaSS_max,rwarpi,psimax,idust)
          close(1)
       endif
    endif
    !--write disc parameters to screen
    if (do_verbose) then
       call write_discinfo(stdout,R_in,R_out,R_ref,Q,npart,sigmaprofile,&
                           R_c,p_index,q_index,star_m,disc_m,real(xinclination*180.0/pi),honH,cs0,&
                           alphaSS_min,alphaSS_max,rwarpi,psimax,itype)
    endif
    if (do_mixture) then
       call write_discinfo(stdout,R_indust,R_outdust,R_ref,Q,npart,sigmaprofiledust,&
                           R_c_dust,p_inddust,q_index,star_m,disc_massdust,real(xinclination*180.0/pi),honH,&
                           cs0,alphaSS_min,alphaSS_max,rwarpi,psimax,idust)
    endif
 endif

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
                              R_indust,R_outdust,phi_min,phi_max,disc_m,disc_mdust,&
                              sigmaprofile,sigmaprofiledust,R_c,R_c_dust,p_index,p_inddust,cs0,&
                              q_index,star_m,G,particle_mass,hfact,itype,xyzh,honH,verbose)
 use domain,  only:i_belong
 use io,      only:id,master
 use part,    only:igas,set_particle_type
 use physcon, only:pi
 use random,  only:ran2
 integer, intent(in)    :: npart_start_count,npart_tot
 real,    intent(in)    :: R_ref,R_in,R_out,phi_min,phi_max
 real,    intent(in)    :: disc_m,p_index,cs0,q_index,star_m,G,particle_mass,hfact
 real,    intent(in)    :: disc_mdust,R_indust,R_outdust,R_c,R_c_dust,p_inddust
 logical, intent(in)    :: do_mixture,verbose
 integer, intent(in)    :: itype,sigmaprofile,sigmaprofiledust
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(out)   :: honH
 integer :: i,j,iseed,ninz
 integer :: ipart
 real    :: rand_no,randtest,R,phi,zi
 real    :: f,fr_max,fz_max,sigma,cs,omega,fmixt
 real    :: HH,HHsqrt2,z_min,z_max
 real    :: rhopart,rhoz,hpart
 real    :: xcentreofmass(3)
 real    :: sigma_norm,sigma_normdust,dR,dRmixt
 integer, parameter :: nbins=10000
 real, dimension(nbins) :: f_vals,fmixt_vals

 !--seed for random number generator
 iseed = -34598 + (itype - igas)
 honH = 0.
 ninz = 0

 !--set maximum f=R*sigma (scaled) value
 sigma_norm = disc_m / scaled_discmass(sigmaprofile,p_index,R_in,R_out,R_ref,R_c)
 dR = (R_out-R_in)/real(nbins-1)
 do i=1,nbins
    R = R_in + (i-1)*dR
    f_vals(i) = R*scaled_sigma(R,sigmaprofile,p_index,R_ref,R_in,R_c)
 enddo
 fr_max = maxval(f_vals)
 if (do_mixture) then
    sigma_normdust = disc_mdust / scaled_discmass(sigmaprofiledust,p_inddust,R_indust,R_outdust,R_ref,R_c_dust)
    dRmixt = (R_outdust-R_indust)/real(nbins-1)
    do i=1,nbins
       R = R_indust + (i-1)*dRmixt
       fmixt_vals(i) = R*scaled_sigma(R,sigmaprofiledust,p_inddust,R_ref,R_indust,R_c_dust)
    enddo
    fr_max = fr_max + maxval(fmixt_vals)
 endif

 xcentreofmass(:) = 0.
 ipart = npart_start_count - 1

 !--loop over particles
 do i = npart_start_count,npart_tot
    if (id==master .and. mod(i,npart_tot/10)==0 .and. verbose) print*,i
    !--get a random angle between phi_min and phi_max
    rand_no = ran2(iseed)
    phi  = phi_min + (phi_max - phi_min)*ran2(iseed)
    !--now get radius
    f = 0.
    randtest = 1.
    fmixt = 0.
    do while (randtest > f)
       R = R_in + (R_out - R_in)*ran2(iseed)
       randtest = fr_max*ran2(iseed)
       f = R*scaled_sigma(R,sigmaprofile,p_index,R_ref,R_in,R_c)
       sigma = sigma_norm*f/R
       if (do_mixture) then
          if (R>=R_indust .and. R<=R_outdust) then
             fmixt = R*scaled_sigma(R,sigmaprofiledust,p_inddust,R_ref,R_indust,R_c_dust)
             f     = f + fmixt
             sigma = sigma + sigma_normdust*fmixt/R
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
    z_min = -3.0d0*HHsqrt2
    z_max =  3.0d0*HHsqrt2
    fz_max = 1.0d0
    f = 0.
    randtest = 1.
    do while(randtest > f)
       zi = z_min + (z_max - z_min)*ran2(iseed)
       randtest = fz_max*ran2(iseed)
       f = exp(-(zi/(HHsqrt2))**2)
       rhoz = f/(HHsqrt2*sqrt(pi))
    enddo

    !--starting estimate of smoothing length for h-rho iterations
    rhopart = sigma*rhoz
    hpart = hfact*(particle_mass/rhopart)**(1./3.)

    if (i_belong(i)) then
       ipart = ipart + 1
       !--set positions -- move to origin below
       xyzh(1,ipart) = R*cos(phi)
       xyzh(2,ipart) = R*sin(phi)
       xyzh(3,ipart) = zi
       xyzh(4,ipart) = hpart
       !--set particle type
       call set_particle_type(ipart,itype)
    endif

    !--HH is scale height
    if (zi*zi < HH*HH) then
       ninz = ninz + 1
       honH = honH + hpart/HH
    endif
 enddo

 !--set honH
 honH = honH/real(ninz)

end subroutine set_disc_positions

!----------------------------------------------------------------
!
! set up the particle velocities
!
!----------------------------------------------------------------
subroutine set_disc_velocities(npart_tot,npart_start_count,itype,G,star_m,aspin,&
                               clight,cs0,do_sigmapringle,p_index,q_index,gamma,R_in,&
                               maxbins,rad,enc_m,smooth_sigma,xyzh,vxyzu,ierr)
 use domain,         only:i_belong
 use externalforces, only:iext_lensethirring,iext_einsteinprec
 use options,        only:iexternalforce
 use part,           only:gravity,igas,maxvxyzu
 integer, intent(in)    :: npart_tot,npart_start_count,itype,maxbins
 real,    intent(in)    :: G,star_m,aspin,clight,cs0,p_index,q_index,gamma,R_in
 real,    intent(in)    :: rad(maxbins),enc_m(maxbins)
 logical, intent(in)    :: do_sigmapringle,smooth_sigma
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 integer, intent(out)   :: ierr
 real :: term,term_pr,term_bh,det,vr,vphi
 real :: cs,R,phi
 integer :: i,itable
 integer :: ipart

 ierr = 0
 ipart = npart_start_count -1
 do i = npart_start_count, npart_tot
    if (i_belong(i)) then
       ipart = ipart + 1
       !
       !--set velocities to give centrifugal balance:
       !  v_phi^2 = GM/R - term_pr - 2*vphi*term_bh
       !
       R    = sqrt(xyzh(1,ipart)**2 + xyzh(2,ipart)**2)
       phi  = atan2(xyzh(2,ipart),xyzh(1,ipart))
       term = G*star_m/R
       !--need to declare iexternalforce before calling setdisc
       if (iexternalforce==11) term=term*(1.0 + (6.0/R)) ! assumes Rg=1.
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
             if (smooth_sigma .and. R > R_in) then  ! R > R_in can happen because of disc shifting
                term_pr = -cs**2*(1.5+p_index+q_index - 0.5/(sqrt(R/R_in) - 1.))
             else
                term_pr = -cs**2*(1.5+p_index+q_index)
             endif
          endif
          if (term + term_pr < 0.) then
             ierr = 1
          endif
       case default
          cs = 0.
          term_pr = 0.
       end select
       !
       !--correction due to Lense-Thirring precession:
       !  Nealon, Nixon & Price correction for Nelson & Papaloizou v x h term
       !
       term_bh = 0.
       if (aspin > tiny(aspin) .and. (iexternalforce==iext_lensethirring &
                                .or. iexternalforce==iext_einsteinprec)) then
          !--this is Eq. 5.21 in Nealon (2013) multiplied by -R
          term_bh = -2.*aspin*(G*star_m/R)**2/clight**3
       endif
       !
       !--now solve quadratic equation for vphi
       !
       det = term_bh**2 + 4.*(term + term_pr)
       vphi = 0.5*(term_bh + sqrt(det))
       !
       !--radial velocities (zero in general)
       !
       vr = 0.d0
       !
       !--set velocities -- move to origin below
       !
       vxyzu(1,ipart)  = -vphi*sin(phi)+ vr*cos(phi)
       vxyzu(2,ipart)  = vphi*cos(phi) + vr*sin(phi)
       vxyzu(3,ipart)  = 0.0d0

       !
       !--set thermal energy
       !

       !--utherm generally should not be stored
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

end subroutine set_disc_velocities

!-------------------------------------------------------------
!
! shift the particles so the centre of mass is at the origin
!
!-------------------------------------------------------------
subroutine adjust_centre_of_mass(xyzh,vxyzu,particle_mass,i1,i2,x0,v0)
 use domain,     only:i_belong
 use mpiutils,   only:reduceall_mpi
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real,    intent(in)    :: particle_mass
 integer, intent(in)    :: i1,i2
 real,    intent(in)    :: x0(3),v0(3)
 real :: xcentreofmass(3), vcentreofmass(3)
 integer :: i,ipart
 real    :: totmass

 xcentreofmass = 0.
 vcentreofmass = 0.
 totmass       = 0.
 ipart = 0
 do i=i1,i2
    if (i_belong(i)) then
       ipart = ipart + 1
       xcentreofmass(:) = xcentreofmass(:) + particle_mass*xyzh(1:3,ipart)
       vcentreofmass(:) = vcentreofmass(:) + particle_mass*vxyzu(1:3,ipart)
       totmass = totmass + particle_mass
    endif
 enddo

 totmass = reduceall_mpi('+',totmass)

 xcentreofmass(:) = xcentreofmass(:)/totmass
 vcentreofmass(:) = vcentreofmass(:)/totmass

 xcentreofmass = reduceall_mpi('+',xcentreofmass)
 vcentreofmass = reduceall_mpi('+',vcentreofmass)

 ipart = 0
 do i=i1,i2
    if (i_belong(i)) then
       ipart = ipart + 1
       xyzh(1:3,ipart)  = xyzh(1:3,ipart)  - xcentreofmass(:) + x0(:)
       vxyzu(1:3,ipart) = vxyzu(1:3,ipart) - vcentreofmass(:) + v0(:)
    endif
 enddo

end subroutine adjust_centre_of_mass

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up disc warps
!
!----------------------------------------------------------------
pure subroutine set_warp(npart_tot,npart_start_count,&
                         xyzh,vxyzu,inclination,sininclination,&
                         rwarpi,psimax,rsi,rso,do_twist)
 use physcon, only:pi
 real, optional, intent(in)    :: inclination,sininclination
 real,           intent(inout) :: xyzh(:,:)
 real,           intent(inout) :: vxyzu(:,:)
 integer,        intent(in)    :: npart_start_count,npart_tot
 real,           intent(in)    :: rwarpi,rsi,rso
 logical,        intent(in)    :: do_twist
 real,           intent(out)   :: psimax
 integer :: i
 real    :: r,cosi,sini,xi,yi,zi,vxi,vyi,vzi,xinclination,r2
 real    :: psi

 psimax = 0.
 do i=npart_start_count, npart_tot
    !--inclined or warped discs
    if (present(inclination) .or. present(sininclination)) then
       !--rotate positions and velocities by required amount
       xi  = xyzh(1,i)
       yi  = xyzh(2,i)
       zi  = xyzh(3,i)
       r2  = xi*xi + yi*yi + zi*zi
       !
       !--these are two different ways of doing the same thing:
       !  1) specify inclination in radians, with an
       !  optional unsmoothed twist/warp at a specified radius
       !  2) specify sine of inclination angle, with an
       !  optional smoothed twist/warp at rwarp, smoothness
       !  specified with warp_smoothl
       !
       if (present(inclination)) then
          if (do_twist) then
             if (r2 < rwarpi**2) then
                xinclination = 0.0
             else
                xinclination = inclination
             endif
          else
             xinclination = inclination
          endif
          cosi = cos(xinclination)
          sini = sin(xinclination)
       else
          r = sqrt(r2)
          if (do_twist) then
             !--this is Equation 43 in Lodato & Price (2010)
             psi = 0.
             if (r < rsi) then
                sini = 0.
             elseif (r < rso) then
                !--sini = sininclination*exp(-abs(r-rwarpi)**2)
                sini = 0.5*(1.+sin(pi/(rso - rsi)*(r-rwarpi)))*sininclination
                psi  = pi*rwarpi/(2.*(rso - rsi))*sininclination/sqrt(1. - (0.5*sininclination)**2)
                psimax = max(psimax,psi)
             else
                sini = sininclination
             endif
          else
             sini = sininclination
          endif
          !--uncomment the following for Nelson-Papaloizou bending waves
          !if (r < rsi.or.r > rso) then
          !  sini=0.
          !else
          !  sini=ampl*sin(pi*(r-rsi)/(rso-rsi))
          !endif
          cosi = sqrt(1. - sini**2)
       endif
       vxi = vxyzu(1,i)
       vyi = vxyzu(2,i)
       vzi = vxyzu(3,i)
       !--rotate positions
       xyzh(1,i)  =  xi                  !xi*cosi + zi*sini
       xyzh(2,i)  =  yi*cosi - zi*sini   !yi
       xyzh(3,i)  =  yi*sini + zi*cosi   !-xi*sini + zi*cosi
       !--rotate velocities
       vxyzu(1,i) =  vxi                 !vxi*cosi + vzi*sini
       vxyzu(2,i) =  vyi*cosi - vzi*sini !vyi
       vxyzu(3,i) =  vyi*sini + vzi*cosi !-vxi*sini + vzi*cosi
    endif
 enddo

end subroutine set_warp

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
subroutine write_discinfo(iunit,R_in,R_out,R_ref,Q,npart,sigmaprofile,&
                          R_c,p_index,q_index,star_m,disc_m,inclination,honH,cs0,&
                          alphaSS_min,alphaSS_max,R_warp,psimax,itype)
 use dim,          only:maxvxyzu
 use eos,          only:get_temperature,init_eos,ieos
 use infile_utils, only:write_inopt
 use part,         only:labeltype,idust
 use physcon,      only:gg,c
 use units,        only:umass,utime,udist
 integer, intent(in) :: iunit,npart,itype,sigmaprofile
 real,    intent(in) :: R_in,R_out,R_ref,Q,p_index,q_index,star_m,disc_m,inclination,honH,cs0
 real,    intent(in) :: alphaSS_min,alphaSS_max,R_warp,psimax,R_c
 integer :: ierr,i
 real    :: T0,T_ref,sig,sigma_norm,dR,R
 real,    parameter :: vxyzutmp(maxvxyzu) = 0.
 integer, parameter :: nbins=10000
 real, dimension(nbins) :: sig_vals

 write(iunit,"(/,a)") '# '//trim(labeltype(itype))//' disc parameters'
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
 sigma_norm = disc_m / scaled_discmass(sigmaprofile,p_index,R_in,R_out,R_ref,R_c)
 sig = sigma_norm*scaled_sigma(R_in,sigmaprofile,p_index,R_ref,R_in,R_c)
 sig = sig*umass/udist**2
 call write_inopt(sig,'sig_in','surface density (g/cm^2) at R=R_in',iunit)
 sig = sigma_norm*scaled_sigma(R_ref,sigmaprofile,p_index,R_ref,R_in,R_c)
 sig = sig*umass/udist**2
 call write_inopt(sig,'sig_ref','surface density (g/cm^2) at R=R_ref',iunit)
 sig = sigma_norm*scaled_sigma(R_out,sigmaprofile,p_index,R_ref,R_in,R_c)
 sig = sig*umass/udist**2
 call write_inopt(sig,'sig_out','surface density (g/cm^2) at R=R_out',iunit)
 dR = (R_out-R_in)/real(nbins-1)
 do i=1,nbins
    R = R_in + (i-1)*dR
    sig_vals(i) = sigma_norm*scaled_sigma(R,sigmaprofile,p_index,R_ref,R_in,R_c)
 enddo
 sig = maxval(sig_vals)
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
 call write_inopt(cs0,'cs0','sound speed at R=1',iunit)

 call init_eos(ieos,ierr)
 T0 = get_temperature(ieos,(/R_in,0.,0./),1.,vxyzutmp)
 call write_inopt(T0,'T_in','temperature (K) at R=R_in',iunit)
 T_ref = get_temperature(ieos,(/R_ref,0.,0./),1.,vxyzutmp)
 call write_inopt(T_ref,'T_ref','temperature (K) at R=R_ref',iunit)
 T0 = get_temperature(ieos,(/R_out,0.,0./),1.,vxyzutmp)
 call write_inopt(T0,'T_out','temperature (K) at R=R_out',iunit)

 call write_inopt(inclination,'inc.deg','disc inclination in degrees',iunit)
 call write_inopt(honH,'<h/H>','approx. mean smoothing length over disc scale height',iunit)
 call write_inopt(alphaSS_min,'alphaSS_min','minimum Shakura-Sunyaev alpha viscosity in disc',iunit)
 call write_inopt(alphaSS_max,'alphaSS_max','maximum Shakura-Sunyaev alpha viscosity in disc',iunit)
 call write_inopt(udist,'udist','distance units (cgs)',iunit)
 call write_inopt(umass,'umass','mass units (cgs)',iunit)
 call write_inopt(utime,'utime','time units (cgs)',iunit)
 call write_inopt(gg/(udist**3/(utime**2*umass)),'G','in code units',iunit)
 call write_inopt(c/(udist/utime),'c','in code units',iunit)
 write(iunit,*)

 !--print some of these diagnostics in more useful form
 write(iunit,"(a,f5.1,a,f5.1,a,f4.1,a)") '# Temperature profile  = ',T_ref,'K (R/',R_ref,')^(',-2.*q_index,')'
 if (sigmaprofile==0) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,')'
 elseif (sigmaprofile==1) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,f5.1,a,f4.1,a/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,') exp[-(R/',R_c,')^(2-',p_index,')]'
 elseif (sigmaprofile==2) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,') (1 - sqrt(',R_in,'/R))'
 elseif (sigmaprofile==3) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,f5.1,a,f4.1,a,f4.1,a,/)") '# Surface density      = ',&
         sigma_norm*umass/udist**2,' g/cm^2 (R/',R_ref,')^(',-p_index,') exp[-(R/',R_c,')^(2-',p_index,')] (1 - sqrt(',R_in,'/R))'
 endif

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
!   R_in - inner edge of disc
!   i1,i2 - range of particle indices to use
!   rwarp - location of warp (optional)
!
!  Output:
!   honHmin, honH, honHmax - max mean and min ratio of h/H within range of R
!-----------------------------------------------------------------------------
subroutine get_honH(xyzh,rminav,rmaxav,honHmin,honHmax,honH,npart,cs0,q_index,M_star,R_in,i1,i2,verbose,rwarp)
 use domain,   only:i_belong
 use mpiutils, only:reduceall_mpi
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: honHmax,honHmin,honH
 integer, intent(in)  :: npart,i1,i2
 real,    intent(in)  :: cs0,q_index,M_star,rminav,rmaxav,R_in
 logical, intent(in)  :: verbose
 real,    intent(in), optional :: rwarp

 integer, parameter :: nr = 350
 integer :: i,ii,iwarp
 real :: G,rmin,rmax,dr,ri
 real :: rad(nr),ninbin(nr),h_smooth(nr),cs(nr),H(nr),omega(nr)
 integer :: ipart

 G = 1.0

 !--setup rmin and rmax for the analysis
 rmin = rminav
 rmax = rmaxav

 !--set up the radius array
 dr = (rmax-rmin)/real(nr-1)
 iwarp = 0
 do i=1,nr
    rad(i)=rmin + real(i-1)*dr
    if (present(rwarp)) then
       if (rad(i) > rwarp) iwarp = i - 1
    endif
 enddo

 !--initialise arrays to zero
 ninbin(:)=0
 h_smooth(:)=0.0

 !--and thus the sound speed array
 do i=1,nr
    cs(i) = cs_func(cs0,rad(i),q_index)
    omega(i) = sqrt(G*M_star/rad(i)**3)
 enddo

 !--and thus the disc scale height
 do i=1,nr
    H(i) = cs(i)/omega(i)
 enddo

 !--loop over particles putting properties into the correct bin
 ipart = 0
 do i = i1,i2
    if (i_belong(i)) then
       ipart = ipart + 1
       if (xyzh(4,ipart)  >  tiny(xyzh)) then ! IF ACTIVE
          ri = sqrt(dot_product(xyzh(1:3,ipart),xyzh(1:3,ipart)))
          ii = int((ri-rad(1))/dr + 1)

          if (ii > nr) cycle
          if (ii < 1)  cycle

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
 do i = 1,nr
    if (ninbin(i) > 0) h_smooth(i) = h_smooth(i)/ninbin(i)
 enddo

 !--now loop over rings to calculate required quantities
 do i = 1, nr
    if (H(i) > 0.) then
       h_smooth(i) = h_smooth(i)/H(i)
    endif
 enddo

 !--print out an average value for <h>/H and thus \alpha_SS/\alpha_AV
 honH = 0.0
 honHmin = minval(h_smooth)
 honHmax = maxval(h_smooth)
 do i=1,nr
    honH = honH + h_smooth(i)
 enddo
 honH = honH/real(nr)

 if (verbose) then
    write(*,'(1x,"Actual <h>/H.................. ",f9.4)') honH
    write(*,'(1x,"Actual (<h>/H)_min is approx.: ",f9.4)') honHmin
    write(*,'(1x,"Actual (<h>/H)_max is approx.: ",f9.4)') honHmax
    if (present(rwarp) .and. iwarp > 0) then
       write(*,'(1x,"Actual <h>/H at R_warp is    : ",f9.4)') h_smooth(iwarp)
    endif
    ! write(*,'(1x,"alpha_SS/alpha_AV is approx: ",f9.4)') honH/10.0
    write(*,'(1x,"e.g. for alpha_SS = 0.1,  use alpha_AV = ",f10.5)') 0.1/(honH/10.0)
    write(*,'(1x,"   using alpha_AV = 1.0 gives alpha_SS ~ ",f8.4,"->",f8.4)') honHmin/10.0,honHmax/10.0
 endif

end subroutine get_honH

!------------------------------------------------------------------------
!
! returns surface density (sigma) at any value of R for a given profile
! scaled by the normalisation value sigma_norm defined below
!
! sigmaprofile options:
!
!    0) power law
!         sigma = sigma_norm * (R/Rref)^-p
!
!    1) exponentially tapered power law
!         sigma = sigma_norm * (R/Rref)^-p * exp[-(R/Rc)^(2-p)]
!
!    2) smoothed power law
!         sigma = sigma_norm * (R/Rref)^-p * (1 - sqrt(Rin/R))
!
!    3) both tapered and smoothed
!         sigma = sigma_norm * (R/Rref)^-p * exp[-(R/Rc)^(2-p)] * (1 - sqrt(Rin/R))
!
!    todo: accept any surface density profile from file
!
!------------------------------------------------------------------------
function scaled_sigma(R,sigmaprofile,pindex,R_ref,R_in,R_c) result(sigma)
 use io, only:error
 real,           intent(in)  :: R,R_ref,pindex
 real, optional, intent(in)  :: R_in,R_c
 integer,        intent(in)  :: sigmaprofile
 real :: sigma

 select case (sigmaprofile)
 case (0)
    sigma = (R/R_ref)**(-pindex)
 case (1)
    sigma = (R/R_ref)**(-pindex)*exp(-(R/R_c)**(2-pindex))
 case (2)
    sigma = (R/R_ref)**(-pindex)*(1-sqrt(R_in/R))
 case (3)
    sigma = (R/R_ref)**(-pindex)*exp(-(R/R_c)**(2-pindex))*(1-sqrt(R_in/R))
 case default
    call error('set_disc','unavailable sigmaprofile; surface density is set to zero')
    sigma = 0.
 end select

end function scaled_sigma

!------------------------------------------------------------------------
!
! returns the disc mass (surface density integrated over R=R_in to R_out)
! scaled by the normalisation value sigma_norm as defined in the function
! "scaled_sigma"
!
!------------------------------------------------------------------------
function scaled_discmass(sigmaprofile,pindex,R_in,R_out,R_ref,R_c) result(mass)
 use physcon, only:pi
 real,           intent(in)  :: pindex,R_in,R_out,R_ref
 real, optional, intent(in)  :: R_c
 integer,        intent(in)  :: sigmaprofile

 integer, parameter :: nbins=10000
 real    :: dr,dM,R,mass,sigma
 integer :: i

 mass = 0.
 dR = (R_out-R_in)/real(nbins-1)
 do i=1,nbins
    R = R_in + (i-1)*dR
    sigma = scaled_sigma(R,sigmaprofile,pindex,R_ref,R_in,R_c)
    dM    = 2.*pi*R*sigma*dR
    mass  = mass + dM
 enddo

end function scaled_discmass

!------------------------------------------------------------------------
!
! incline disc around rotation axis defined by a position angle
!
!------------------------------------------------------------------------
pure subroutine rotate_disc(xyzh,vxyzu,npart_tot,npart_start,posangl,incl)
 use vectorutils, only:rotatevec
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(in)    :: posangl,incl
 integer, intent(in)    :: npart_start,npart_tot
 integer :: i
 real :: k(3)

 !--vector in direction of PA
 k = (/-sin(posangl), cos(posangl), 0./)

 do i=npart_start,npart_tot
    !--rotate positions and velocities
    call rotatevec(xyzh(1:3,i),k,incl)
    call rotatevec(vxyzu(1:3,i),k,incl)
 enddo

end subroutine rotate_disc

end module setdisc
