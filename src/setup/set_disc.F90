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
!    H_R         -- disc aspect ratio H/R at R=R_in
!    M_disc      -- disc mass
!    M_star      -- mass of central star
!    Qmin        -- minimum Toomre Q parameter
!    R_c         -- characteristic radius of the exponential taper
!    R_in        -- inner disc boundary
!    R_out       -- outer disc boundary
!    R_warp      -- position of warp
!    Rref        -- reference radius
!    Sig0        -- surface density at R=1
!    T_in        -- temperature (K) at R=R_in
!    T_out       -- temperature (K) at R=R_out
!    alphaSS_max -- maximum Shakura-Sunyaev alpha viscosity in disc
!    alphaSS_min -- minimum Shakura-Sunyaev alpha viscosity in disc
!    c           -- in code units
!    cs0         -- sound speed at R=1
!    n           -- number of particles in the disc
!    p_index     -- power law index of surface density profile
!    psi_max     -- maximum warp amplitude
!    q_index     -- power law index of sound speed profile
!    udist       -- distance units (cgs)
!    umass       -- mass units (cgs)
!    utime       -- time units (cgs)
!
!  DEPENDENCIES: dim, domain, eos, externalforces, infile_utils, io,
!    mpiutils, options, part, physcon, random, units
!+
!--------------------------------------------------------------------------
module setdisc
 implicit none
 public :: set_disc, set_warp

 private

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine set_disc(id,master,mixture,nparttot,npart,npart_start,rmin,rmax,rmindust,rmaxdust,phimin,phimax,&
                    indexprofile,indexprofiledust,rc,rcdust,p_index,p_indexdust,q_index,HoverR,gamma,&
                    disc_Q,disc_mass,disc_massdust,sig_naught,star_mass,xyz_origin,vxyz_origin,&
                    particle_type,particle_mass,hfact,xyzh,vxyzu,polyk,inclination,sininclination,&
                    twist,ismooth,alpha,isink,rwarp,warp_smoothl,bh_spin,R_ref,writefile,ierr,prefix,verbose)
 use domain,  only:i_belong
 use physcon, only:pi,gg,c
 use units,   only:umass,udist,utime
 use eos,     only:qfacdisc
 use options, only:ieos
 use part,    only:maxp,igas,idust,set_particle_type,labeltype,gravity,maxtypes
 use io,      only:fatal,warning,stdout
 integer,                     intent(in)    :: id,master
 integer, optional,           intent(in)    :: nparttot
 integer,                     intent(out)   :: npart
 integer, optional,           intent(in)    :: npart_start,isink,indexprofile,indexprofiledust
 real,                        intent(in)    :: rmin,rmax
 real, optional,              intent(in)    :: rmindust,rmaxdust,p_indexdust,disc_massdust
 real, optional,              intent(in)    :: rc,rcdust
 real, optional,              intent(in)    :: phimin,phimax
 real, optional,              intent(inout) :: alpha
 real,                        intent(in)    :: p_index,q_index,HoverR,gamma,hfact
 real, optional,              intent(in)    :: disc_Q,disc_mass,star_mass,sig_naught
 real, optional,              intent(in)    :: xyz_origin(3), vxyz_origin(3)
 integer, optional,           intent(in)    :: particle_type
 real, optional,              intent(in)    :: inclination,sininclination,rwarp,warp_smoothl,bh_spin,R_ref
 logical, optional,           intent(in)    :: twist,ismooth,mixture
 real,                        intent(out)   :: xyzh(:,:)
 real,                        intent(out)   :: vxyzu(:,:)
 real,                        intent(out)   :: polyk,particle_mass
 logical, optional,           intent(in)    :: writefile,verbose
 integer, optional,           intent(out)   :: ierr
 character(len=20), optional, intent(in)    :: prefix
 integer, parameter :: maxbins = 256
 integer :: maxvxyzu,itype,npart_tot,npart_start_count,ierror
 integer :: i
 real    :: Q,Q_tmp,G,cs0,sig0,clight,sig0dust
 real    :: R_in,R_out,phi_min,phi_max,H_R,R_indust,R_outdust,p_inddust,rc0,rc0dust
 real    :: star_M,disc_M,rminav,rmaxav,honHmin,honHmax
 real    :: honH,alphaSS_min,alphaSS_max
 real    :: xinclination,rwarpi,hwarp,rsi,rso,psimax,tempmassdust
 real    :: aspin
 real    :: enc_m(maxbins), r_c(maxbins)
 real    :: xorigini(3),vorigini(3),rref
 logical :: do_twist,smooth_surface_density,do_write,do_sigmapringle,do_sigmapringledust,do_mixture
 logical :: do_verbose
!
!  Set problem parameters
!
! HoverR is at R=R_in - not at R=1
! (if you would like it set at R=1 rather than R=R_in you must alter cs0 as described below)
! HoverR is a function of radius depending on your choice of q_index
!
 H_R            = HoverR
 R_in           = rmin
 R_out          = rmax
 R_indust       = rmin
 R_outdust      = rmax
 p_inddust      = p_index
 if (present(rmindust) .and. present(rmaxdust) .and. present(p_indexdust))then
    R_indust    = rmindust
    R_outdust   = rmaxdust
    p_inddust   = p_indexdust
 endif
 if (present(star_mass)) then
    Star_M = star_mass
 else
    Star_M = 1.d0
 endif
 if (present(npart_start)) then
    npart_start_count = npart_start
 else
    npart_start_count = 1
 endif
 if (present(disc_Q)) then
    if (present(disc_mass)) then
       if (id==master) &
         print*,' ERROR: set_disc: cannot specify both disc mass and Toomre Q parameter (use only one of these)'
       if (present(ierr)) ierr = 1
       return
    endif
    Q = disc_Q
 else
    Q = 168.d0
 endif
 if (rmax < rmin) then
    if (id==master) print*,' ERROR: outer radius < inner radius in set_disc, doing nothing'
    if (present(ierr)) ierr = 2
    return
 endif
 if (present(rmindust) .and. present(rmaxdust))then
    if (rmaxdust < rmindust) then
       if (id==master) print*,' ERROR: dust outer radius < dust inner radius in set_disc, doing nothing'
       if (present(ierr)) ierr = 2
       return
    endif
 endif
 if (nparttot <= 0) then
    if (id==master) print*,' ERROR: set_disc: nparttot <= 0 in call to set_disc, doing nothing'
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
 if (present(bh_spin)) then
    aspin = bh_spin
 else
    aspin = 0.
 endif
 ! reference radius for normalisation of sigma, temperature profiles
 if (present(R_ref)) then
    Rref = R_ref
 else
    Rref = R_in
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
       if(.not.(present(rmindust) .and. present(rmaxdust) .and. present(p_indexdust) .and. present(disc_massdust)))then
          call fatal('set_disc','setup for dusty disc in the mixture is not specified')
       endif
    endif
 endif

 if (id==master .and. do_verbose) then
    if (do_mixture) then
       print*,' Setting up disc mixture containing ',&
                 nparttot,' '//trim(labeltype(itype))//'/'//trim(labeltype(itype+1))//' particles'
    else
       print*,' Setting up disc containing ',nparttot,' '//trim(labeltype(itype))//' particles'
    endif
 endif
!
! store q_index in eos module, then gets saved to dump header
!
 maxvxyzu = size(vxyzu(:,1))
 qfacdisc = q_index
!
! this is a fudge to ensure that ieos=3 is used if qfacdisc is set and isothermal
! (and now deals with ieos=6 as well)
 if (maxvxyzu < 4 .and. q_index > 0.) then
    if (present(isink)) then
       if (isink==0) then
          call warning('set_disc','setting ieos=3 in input options based on q_index setting and isink=0')
          ieos = 3
       else
          call warning('set_disc','setting ieos=6 in input options based on q_index setting and isink')
          write(*,'("isink = ",i1)') isink
          ieos = 6
       endif
    else
       call warning('set_disc','setting ieos=3 in input options based on q_index setting')
       ieos = 3
    endif
 endif
!
! Set cs0 to give desired H/R
!
! Here is where cs0 is set such that H_R is H/R at R=R_ref
 cs0 = sqrt(1.0/Rref)*(H_R)*sqrt(G*Star_M)*(1.0/Rref)**(-q_index)
! if you would like this defined at R=1 instead, uncomment this line (and comment out the above)
! cs0 = (H_R)*sqrt(G*Star_M)
! if you would like this defined at R=15 instead, uncomment this line (and comment out the above)
! cs0 = sqrt(1.0/15.0)*(H_R)*sqrt(G*Star_M)*(1.0/15.0)**(-q_index)

! set the surface density profile: default option is the power-law
 do_sigmapringle = .false.
 rc0 = 0.
 if (present(indexprofile)) then
    if (indexprofile==1) then
       !surface density profile: power-law tapered by an exponential function (Lynden-Bell & Pringle 1974)
       do_sigmapringle = .true.
       if (present(rc))then
          rc0 = rc
       else
          rc0 = R_out
       endif
    endif
 endif
 do_sigmapringledust = .false.
 rc0dust = 0.
 if (do_mixture) then
    if (present(indexprofiledust)) then
       if (indexprofiledust==1) then
          do_sigmapringledust = .true.
          if (present(rcdust))then
             rc0dust = rcdust
          else
             rc0dust = R_outdust
          endif
       endif
    endif
 endif

 if (present(sig_naught)) then
    sig0 = sig_naught*(1./Rref)**p_index
    call get_disc_mass(sig0,do_sigmapringle,rc0,p_index,smooth_surface_density,cs0,q_index,star_M,G, &
                       r_in,r_out,Q,disc_m,maxbins,enc_m,r_c)
 else
    sig0 = 1.d0
!
! compute trial disc mass and Toomre Q
!
    call get_disc_mass(sig0,do_sigmapringle,rc0,p_index,smooth_surface_density,cs0,q_index,star_M,G, &
                    r_in,r_out,Q_tmp,disc_m,maxbins,enc_m,r_c)
!
! renormalise to give desired value
!
    if (present(disc_mass)) then
       sig0 = sig0*disc_mass/disc_m
    else
       sig0 = Q_tmp / Q
    endif
!
! recompute actual disc mass and Toomre Q
!
    call get_disc_mass(sig0,do_sigmapringle,rc0,p_index,smooth_surface_density,cs0,q_index,star_M,G, &
                       r_in,r_out,Q,disc_m,maxbins,enc_m,r_c)
 endif

 sig0dust = 0.
 if (.not.do_mixture) then
!
! set the particle mass
!
    particle_mass = disc_m / dble(nparttot)
 else
    sig0dust = 1.d0
    call get_disc_mass(sig0dust,do_sigmapringledust,rc0dust,p_inddust,smooth_surface_density,cs0,q_index,star_M,G, &
                    r_indust,r_outdust,Q_tmp,tempmassdust,maxbins,enc_m,r_c)
    sig0dust = sig0dust*disc_massdust/tempmassdust
!
! set the particle mass of the mixture
!
    particle_mass = (disc_m+disc_massdust) / dble(nparttot)
 endif
!
! count particles on this MPI thread
!
 npart = 0
 do i = 1,nparttot
    if (i_belong(i)) npart = npart + 1
 enddo
!
! set particle positions and smoothing lengths
!
 npart_tot = npart_start_count + nparttot - 1
 if (npart_tot > maxp) call fatal('set_disc','number of particles exceeds array dimensions',var='n',ival=npart_tot)

 call set_disc_positions(npart_tot,npart_start_count,do_mixture,R_in,R_out,R_indust,R_outdust,phi_min,phi_max,&
                         sig0,sig0dust,do_sigmapringle,do_sigmapringledust,rc0,rc0dust,p_index,p_inddust,cs0,&
                         q_index,star_M,G,particle_mass,hfact,smooth_surface_density,itype,xyzh,honH,do_verbose)
!
! set velocities of particles in the disc
!
 call set_disc_velocities_u(npart_tot,npart_start_count,itype,G,Star_M,aspin,&
                            clight,cs0,do_sigmapringle,p_index,q_index,gamma,r_in, &
                            maxbins,r_c,enc_m,smooth_surface_density,xyzh,vxyzu,ierror)
 if (ierror==1) then
    call fatal('set_disc','error: pressure correction causing -ve sqrt while setting velocities')
 elseif (ierror /= 0) then
    call fatal('set_disc','error setting velocities')
 endif

 if (present(bh_spin) .and. aspin > 1.d-10 .and. id==master) then
    print*,'Orbital velocities corrected for black hole spin, a = ',aspin
    print*,'ASSUMING THAT THE DISC AND BLACK HOLE ARE ALIGNED'
 endif

 polyk = cs0**2
!
! work out h/H in order to set the artificial viscosity parameter to match a chosen alpha_SS
!
 if (do_verbose) write(*,'(/,1x,"Actual <h>/H...per particle... ",f9.4)') honH
 if (smooth_surface_density .and. p_index > 0.) then
    rminav = R_in*((1.+2.*p_index)**2)/(4.*p_index**2)
 else
    rminav = R_in
 endif
 rmaxav = R_out
 if (present(rwarp)) then
    call get_honH(xyzh,rminav,rmaxav,honHmin,honHmax,honH,nparttot,H_R,q_index,Star_M,R_in,&
                  npart_start_count,npart_tot,do_verbose,rwarp)
 else
    call get_honH(xyzh,rminav,rmaxav,honHmin,honHmax,honH,nparttot,H_R,q_index,Star_M,R_in,&
                  npart_start_count,npart_tot,do_verbose)
 endif

#ifdef DISC_VISCOSITY
!
! if disc viscosity is used, set the artificial viscosity parameter
! in the input file so as to give the desired alpha_SS
!
 if(present(alpha)) then
    if (do_verbose) print*, 'alphaSS requested = ', alpha
    alpha = alpha/(honH/10.0)
! and the min and max alphaSS present
    alphaSS_min = alpha*honHmin/10.
    alphaSS_max = alpha*honHmax/10.
    if (do_verbose) print*, 'Setting alpha_AV = ',alpha,' to give alphaSS as requested'
 else
    alphaSS_min = honHmin/10.
    alphaSS_max = honHmax/10.
 endif
#else
!
! if disc viscosity is not used, simply return the range of alphaSS
! implied in the disc by the chosen artificial viscosity parameter
!
 alphaSS_min = honHmin*(31./525.)
 alphaSS_max = honHmax*(31./525.)
#endif

!
! add a warp/twist to the disc
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
!---------------------------------------------
! Call setwarp to actually calculate the warp
 call set_warp(npart_tot,npart_start_count,&
               xyzh,vxyzu,inclination,sininclination,&
               rwarpi,psimax,rsi,rso,do_twist)
!---------------------------------------------
 !
 ! adjust positions and velocities so the centre of mass is at the origin
 ! also shift particles to new origin if this is not at (0,0,0)
 !
 if (present(phimax)) then
    print "(a)",'Setting up disc sector - not adjusting centre of mass'
 else
    call adjust_centre_of_mass(xyzh,vxyzu,particle_mass,npart_start_count,npart_tot,xorigini,vorigini)
 endif

 !  Print out disc parameters, to file and to the screen
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
       call write_discinfo(1,R_in,R_out,Rref,H_R,Q,npart,do_sigmapringle,rc0,p_index,q_index,star_m,disc_m,&
             real(xinclination*180.0/pi),honH,sig0,cs0,alphaSS_min,alphaSS_max,rwarpi,psimax,itype)
       close(1)
       if (do_mixture) then
          open(1,file=trim(prefix)//'-'//trim(labeltype(idust))//'.discparams',status='replace',form='formatted')
          call write_discinfo(1,R_indust,R_outdust,Rref,H_R,Q,npart,do_sigmapringledust,rc0dust,p_inddust,q_index,&
               star_m,disc_massdust,real(xinclination*180.0/pi),honH,sig0dust,cs0,alphaSS_min,alphaSS_max,rwarpi,psimax,idust)
          close(1)
       endif
    endif
    ! write disc parameters to screen
    if (do_verbose) then
       call write_discinfo(stdout,R_in,R_out,Rref,H_R,Q,npart,do_sigmapringle,rc0,p_index,q_index,star_m,disc_m,&
          real(xinclination*180.0/pi),honH,sig0,cs0,alphaSS_min,alphaSS_max,rwarpi,psimax,itype)
    endif
    if (do_mixture) then
       call write_discinfo(stdout,R_indust,R_outdust,Rref,H_R,Q,npart,do_sigmapringledust,&
                           rc0dust,p_inddust,q_index,star_m,disc_massdust,real(xinclination*180.0/pi),&
                           honH,sig0dust,cs0,alphaSS_min,alphaSS_max,rwarpi,psimax,idust)
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

!----------------------------------------------------------------
!
! function to return the disc mass and Toomre Q parameter
! also returns table of enclosed mass enc_m(i) within a radius r_c(i)
!
!----------------------------------------------------------------
subroutine get_disc_mass(sig0,do_sigmapringle,rc0,p_index,smooth_sigma,cs0,q_index,star_M,G,&
                              R_in,R_out,Q,disc_m,maxbins,enc_m,r_c)
 use physcon, only:pi
 real,    intent(in)  :: sig0,p_index,cs0,q_index,star_M,G,R_in,R_out,rc0
 real,    intent(out) :: Q,disc_m
 integer, intent(in)  :: maxbins
 logical, intent(in)  :: do_sigmapringle,smooth_sigma
 real,    intent(out) :: enc_m(maxbins), r_c(maxbins)
 real    :: r,dr,cs,sigma,kappa,dm
 integer :: i

 Q  = 1.0d10
 disc_m = 0.
 dr = (R_out - R_in)/real(maxbins-1)
 enc_m(:) = 0.
 do i = 1, maxbins
    r = R_in + (i-1)*dr
    cs = cs_func(cs0,r,q_index)
    if (do_sigmapringle) then
       sigma=sig0*((r/rc0)**(-p_index))*exp(-(r/rc0)**(2.-p_index))
    else
       sigma = sig0*r**(-p_index)
    endif

! Use next lines if sigma tends to 0 at the inner edge
    if (smooth_sigma) then
       sigma = sigma*(1.0d0-sqrt(R_in/r))
    endif
! -----------------------------------------------------
    kappa  = sqrt(G*Star_M/r**3) ! epicyclic frequency = orbital freq if Keplerian
    if (sigma > epsilon(sigma)) then
       Q      = min(Q,real(cs*kappa/(pi*G*sigma)))
    endif
    dm     = 2.0d0*pi*r*sigma*dr
    disc_m   = disc_m + dm
    enc_m(i) = disc_m
    r_c(i)   = r+0.5d0*dr
 enddo
 enc_m(:) = enc_m(:) + star_M

end subroutine get_disc_mass

!---------------------------------------------------------
!
! Set up the particle positions and smoothing length
! using a Monte Carlo approach
!
!---------------------------------------------------------
subroutine set_disc_positions(npart_tot,npart_start_count,do_mixture,R_in,R_out,R_indust,R_outdust,phi_min,phi_max,&
                              sig0,sig0dust,do_sigmapringle,do_sigmapringledust,rc0,rc0dust,p_index,p_inddust,cs0,&
                              q_index,star_M,G,particle_mass,hfact,smooth_sigma,itype,xyzh,honH,verbose)
 use domain,  only:i_belong
 use physcon, only:pi
 use part,    only:igas,set_particle_type
 use random,  only:ran2
 use io,      only:fatal,id,master
 integer, intent(in)    :: npart_start_count,npart_tot
 real,    intent(in)    :: R_in,R_out,phi_min,phi_max
 real,    intent(in)    :: sig0,p_index,cs0,q_index,star_M,G,particle_mass,hfact
 real,    intent(in)    :: R_indust,R_outdust,sig0dust,rc0,rc0dust,p_inddust
 logical, intent(in)    :: do_sigmapringle,smooth_sigma,do_mixture,do_sigmapringledust,verbose
 integer, intent(in)    :: itype
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(out)   :: honH
 integer :: i,iseed,ninz
 integer :: ipart
 real    :: rand_no,randtest,r,phi,zi
 real    :: f,f_max,sigma,cs,omega,fmixt
 real    :: HH,HHsqrt2,z_min,z_max
 real    :: rhopart,rhoz,hpart,r_max
 real :: xcentreofmass(3)

! seed for random number generator
 iseed = -34598 + (itype - igas)
 honH = 0.
 ninz = 0

 xcentreofmass(:) = 0.
 ipart = 0
 do i = npart_start_count, npart_tot
    if (id==master .and. mod(i,npart_tot/10)==0 .and. verbose) print*,i
    !
    ! get a random angle between phi_min and phi_max
    !
    rand_no = ran2(iseed)
    phi  = phi_min + (phi_max - phi_min)*ran2(iseed)
    ! Now get radius
    if (do_sigmapringle) then
       ! power law profile tapered by an exponential function (Lynden-Bell & Pringle 1974)
       if(p_index < 1.) then
          r_max = rc0*exp(log((1.-p_index)/(2.-p_index))/(2.-p_index))
          f_max = rc0*((r_max/rc0)**(1.-p_index))*exp(-(r_max/rc0)**(2.-p_index))
       else
          f_max = rc0*((R_in/rc0)**(1.-p_index))*exp(-(R_in/rc0)**(2.-p_index))
          if (smooth_sigma) then
             f_max = rc0*((R_in/rc0)**(1.-p_index))*exp(-(R_in/rc0)**(2.-p_index))/4.
          endif
       endif
    else
       ! power law profile
       if(p_index <= 1.) then
          f_max = R_out**(-(p_index-1.))
       else
          f_max = R_in**(-(p_index-1.))
!      Use following lines if sigma = 0 at inner edge
          if (smooth_sigma) then
             f_max = R_in**(-(p_index-1.))/4.
          endif
       endif
    endif
    if (do_mixture) then
       if (do_sigmapringledust) then
          ! power law profile tapered by an exponential function (Lynden-Bell & Pringle 1974)
          if(p_inddust < 1.) then
             r_max = rc0dust*exp(log((1.-p_inddust)/(2.-p_inddust))/(2.-p_inddust))
             f_max = f_max + rc0dust*((r_max/rc0dust)**(1.-p_inddust))*exp(-(r_max/rc0dust)**(2.-p_inddust))
          else
             f_max = f_max + rc0dust*((R_indust/rc0dust)**(1.-p_inddust))*exp(-(R_indust/rc0dust)**(2.-p_inddust))
             if (smooth_sigma) then
                f_max = f_max + rc0dust*((R_indust/rc0dust)**(1.-p_inddust))*exp(-(R_indust/rc0dust)**(2.-p_inddust))/4.
             endif
          endif
       else
          ! power law profile
          if(p_inddust <= 1.) then
             f_max = f_max + R_outdust**(-(p_inddust-1.))
          else
             f_max = f_max + R_indust**(-(p_inddust-1.))
!         Use following lines if sigma = 0 at inner edge
             if (smooth_sigma) then
                f_max = f_max + R_indust**(-(p_inddust-1.))/4.
             endif
          endif
       endif
    endif
! ----------------------------------------------------
    f = 0.
    randtest = 1.
    fmixt = 0.
    do while (randtest > f)
       r = R_in + (R_out - R_in)*ran2(iseed)
       randtest = f_max*ran2(iseed)
       !--this function is R*sigma
       if (do_sigmapringle) then
          f = rc0*((r/rc0)**(1.-p_index))*exp(-(r/rc0)**(2.-p_index))
       else
          f = r**(-(p_index-1.))
       endif
       if (smooth_sigma) then
          f = f*(1.0d0-sqrt(R_in/r))
       endif
       sigma = sig0*f/r
       if (do_mixture) then
          if (r>=R_indust .and. r<=R_outdust) then
             !--this function is R*sigma
             if (do_sigmapringledust) then
                fmixt = rc0dust*((r/rc0dust)**(1.-p_inddust))*exp(-(r/rc0dust)**(2.-p_inddust))*sig0dust/sig0
             else
                fmixt = r**(-(p_inddust-1.))*sig0dust/sig0
             endif
             if (smooth_sigma) then
                fmixt = r**(-(p_inddust-1.))*(1.0d0-sqrt(R_indust/r))*sig0dust/sig0
             endif
             f = f + fmixt
             sigma = sigma+sig0dust*fmixt/r
          endif
       endif
    enddo
! ----------------------------------------------------
!
! Finally, get z
!
    ! First get sound speed at r
    cs = cs_func(cs0,r,q_index)
    ! Then pressure scale-height - multiplied by sqrt(2)
    ! for convenience in rhoz calc.
    omega   = sqrt(G*Star_M/r**3)
    HH      = cs/omega
    HHsqrt2 = sqrt(2.0d0)*HH
! -----------------------------------------------
    z_min = -3.0d0*HHsqrt2
    z_max =  3.0d0*HHsqrt2
    f_max = 1.0d0
    f = 0.
    randtest = 1.
    do while(randtest > f)
       zi = z_min + (z_max - z_min)*ran2(iseed)
       randtest = f_max*ran2(iseed)
       f = exp(-(zi/(HHsqrt2))**2)
       rhoz = f/(HHsqrt2*sqrt(pi))
    enddo

    !------------------------------------------
    !  Starting estimate of smoothing length for h-rho iterations
    rhopart = sigma*rhoz
    hpart = hfact*(particle_mass/rhopart)**(1./3.)

    if (i_belong(i)) then
       ipart = ipart + 1
       !----------------------------------------------
       ! Set positions -- move to origin below
       xyzh(1,ipart) = r*cos(phi)
       xyzh(2,ipart) = r*sin(phi)
       xyzh(3,ipart) = zi
       xyzh(4,ipart) = hpart

       !------------------------------------------
       !  Set particle type
       call set_particle_type(ipart,itype)
    endif

    ! HH is scale height
    if (zi*zi < HH*HH) then
       ninz = ninz + 1
       honH = honH + hpart/HH
    endif
 enddo
!
! Set honH
!
 honH = honH/real(ninz)

end subroutine set_disc_positions

!----------------------------------------------------------------
!
! Set up the particle positions using a Monte Carlo approach
!
!----------------------------------------------------------------
subroutine set_disc_velocities_u(npart_tot,npart_start_count,itype,G,Star_M,aspin,&
                                      clight,cs0,do_sigmapringle,p_index,q_index,gamma,r_in,&
                                      maxbins,r_c,enc_m,smooth_sigma,xyzh,vxyzu,ierr)
 use domain,         only:i_belong
 use part,           only:gravity,igas,maxvxyzu
 use options,        only:iexternalforce
 use externalforces, only:iext_lensethirring,iext_einsteinprec
 use io,             only:fatal
 integer, intent(in)    :: npart_tot,npart_start_count,itype,maxbins
 real,    intent(in)    :: G,Star_M,aspin,clight,cs0,p_index,q_index,gamma,r_in
 real,    intent(in)    :: r_c(maxbins),enc_m(maxbins)
 logical, intent(in)    :: do_sigmapringle,smooth_sigma
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 integer, intent(out)   :: ierr
 real :: term,term_pr,term_bh,det,vr,vphi
 real :: cs,r,phi
 integer :: i,itable
 integer :: ipart

 ierr = 0
 ipart = 0
 do i = npart_start_count, npart_tot
    if (i_belong(i)) then
       ipart = ipart + 1
       !----------------------------------------------
       !   Set velocities to give centrifugal balance:
       !   v_phi^2 = GM/r - term_pr - 2*vphi*term_bh
       !
       r    = sqrt(xyzh(1,ipart)**2 + xyzh(2,ipart)**2)
       phi  = atan2(xyzh(2,ipart),xyzh(1,ipart))
       term = G*Star_M/r
       ! Need to declare iexternalforce before calling setdisc
       if (iexternalforce==11) term=term*(1.0 + (6.0/r)) ! assumes Rg=1.
       !
       !  Correction due to self-gravity
       !
       if (gravity) then
          itable = nint((r-r_c(1))/(r_c(2)-r_c(1))) + 1
          term = G*enc_m(itable)/r
       endif
       !
       !  Add contribution from pressure gradients
       !  Pressure contribution from a powerlaw disc
       !
       select case(itype)
       case(igas)
          cs = cs_func(cs0,r,q_index)
          if (do_sigmapringle) then
             term_pr = 0.
          else
             if (smooth_sigma .and. r > R_in) then  ! r > R_in can happen because of disc shifting
                term_pr = -cs**2*(1.5+p_index+q_index - 0.5/(sqrt(r/R_in) - 1.))
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
       !  Correction due to Lense-Thirring precession:
       !  Nealon, Nixon & Price correction for Nelson & Papaloizou v x h term
       !
       term_bh = 0.
       if (aspin > tiny(aspin) .and. (iexternalforce==iext_lensethirring &
                                .or. iexternalforce==iext_einsteinprec)) then
          ! this is Eq. 5.21 in Nealon (2013) multiplied by -r
          term_bh = -2.*aspin*(G*Star_M/r)**2/clight**3
       endif
       !
       !  now solve quadratic equation for vphi
       !
       det = term_bh**2 + 4.*(term + term_pr)
       vphi = 0.5*(term_bh + sqrt(det))
       !
       ! radial velocities (zero in general)
       !
       vr = 0.d0

       !  --------------------------------------
       !  Set velocities -- move to origin below
       vxyzu(1,ipart)  = -vphi*sin(phi)+ vr*cos(phi)
       vxyzu(2,ipart)  = vphi*cos(phi) + vr*sin(phi)
       vxyzu(3,ipart)  = 0.0d0

       !-----------------------------------------
       !  Set thermal energy

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

end subroutine set_disc_velocities_u

!-------------------------------------------------------------
! shift the particles so the centre of mass is at the origin
!-------------------------------------------------------------
subroutine adjust_centre_of_mass(xyzh,vxyzu,particle_mass,i1,i2,x0,v0)
 use domain,      only:i_belong
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
       vxyzu(1:3,ipart)  = vxyzu(1:3,ipart)  - vcentreofmass(:) + v0(:)
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
! Rotate positions and velocities by required amount
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
             ! This is Equation 43 in Lodato & Price (2010)
             psi = 0.
             if (r < rsi) then
                sini = 0.
             elseif (r < rso) then
                !  sini = sininclination*exp(-abs(r-rwarpi)**2)
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
          !if(r < rsi.or.r > rso) then
          !  sini=0.
          !else
          !  sini=ampl*sin(pi*(r-rsi)/(rso-rsi))
          !endif
          cosi = sqrt(1. - sini**2)
       endif

       vxi = vxyzu(1,i)
       vyi = vxyzu(2,i)
       vzi = vxyzu(3,i)

! Rotate positions
       xyzh(1,i)  =  xi*cosi + zi*sini
       xyzh(2,i)  =  yi
       xyzh(3,i)  = -xi*sini + zi*cosi
! Rotate velocities
       vxyzu(1,i) =  vxi*cosi + vzi*sini
       vxyzu(2,i) =  vyi
       vxyzu(3,i) = -vxi*sini + vzi*cosi
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
subroutine write_discinfo(iunit,R_in,R_out,Rref,H_R,Q,npart,do_sigmapringle,rc0,p_index,q_index,star_m,&
                          disc_m,inclination,honH,sig0,cs0,alphaSS_min,alphaSS_max,R_warp,psimax,itype)
 use units,        only:umass,utime,udist
 use infile_utils, only:write_inopt
 use physcon,      only:gg,c
 use part,         only:labeltype,idust
 use eos,          only:get_temperature,init_eos,ieos
 use dim, only:maxvxyzu
 integer, intent(in) :: iunit,npart,itype
 real,    intent(in) :: R_in,R_out,Rref,H_R,Q,p_index,q_index,star_m,disc_m,inclination,honH,sig0,cs0
 real,    intent(in) :: alphaSS_min,alphaSS_max,R_warp,psimax,rc0
 logical, intent(in) :: do_sigmapringle
 integer :: ierr
 real :: T0
 real, parameter :: vxyzutmp(maxvxyzu) = 0.

 write(iunit,"(/,a)") '# '//trim(labeltype(itype))//' disc parameters'
 call write_inopt(R_in,'R_in','inner disc boundary',iunit)
 call write_inopt(R_out,'R_out','outer disc boundary',iunit)
 call write_inopt(Rref,'Rref','reference radius',iunit)
 if (R_warp > 0.) call write_inopt(R_warp,'R_warp','position of warp',iunit)
 if (psimax > 0.) call write_inopt(psimax,'psi_max','maximum warp amplitude',iunit)
 call write_inopt(H_R,'H_R','disc aspect ratio H/R at R=R_in',iunit)
! call write_inopt(get_HonR(R_in,cs0,q_index,star_m,1.),'H_R','disc aspect ratio H/R at R=R_in',iunit)
 call write_inopt(get_HonR(Rref,cs0,q_index,star_m,1.),'H/R (R=Rref)','disc aspect ratio H/R at R=Rref',iunit)
 call write_inopt(get_HonR(R_out,cs0,q_index,star_m,1.),'H/R (R=R_out)','disc aspect ratio H/R at R=R_out',iunit)
 if (R_warp > 0.) then
    call write_inopt(get_HonR(R_warp,cs0,q_index,star_m,1.),'H/R (R=R_warp)','disc aspect ratio H/R at R=R_warp',iunit)
 endif
 call write_inopt(Q,'Qmin','minimum Toomre Q parameter',iunit)
 call write_inopt(npart,'n','number of particles in the disc',iunit)
 call write_inopt(p_index,'p_index','power law index of surface density profile',iunit)
 if (do_sigmapringle) call write_inopt(rc0,'R_c','characteristic radius of the exponential taper',iunit)
 call write_inopt(q_index,'q_index','power law index of sound speed profile',iunit)
 call write_inopt(star_m,'M_star','mass of central star',iunit)
 call write_inopt(disc_m,'M_disc','disc mass',iunit)
 call write_inopt(disc_m/star_m,'disc m/star m','relative disc mass',iunit)
 if (do_sigmapringle) then
    call write_inopt(sig0,'Sig0','Sigma0 of the density profile Sigma = Sigma0*(R/Rc)^-p*Exp(-(R/Rc)^(2-p))',iunit)
 else
    call write_inopt(sig0,'Sig0','surface density at R=1',iunit)
 endif
 call write_inopt(cs0,'cs0','sound speed at R=1',iunit)

 call init_eos(ieos,ierr)
 T0 = get_temperature(ieos,(/R_in,0.,0./),1.,vxyzutmp)
 call write_inopt(T0,'T_in','temperature (K) at R=R_in',iunit)
 T0 = get_temperature(ieos,(/R_out,0.,0./),1.,vxyzutmp)
 call write_inopt(T0,'T_out','temperature (K) at R=R_out',iunit)
 T0 = get_temperature(ieos,(/Rref,0.,0./),1.,vxyzutmp)
! call write_inopt(T0,'T_ref','temperature (K) at R=R_ref',iunit)
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

! print some of these diagnostics in more useful form
 write(iunit,"(a,f5.1,a,f5.1,a,f4.1,a)") '# Temperature profile  = ',T0,'K (R/',Rref,')^(',-2.*q_index,')'
 if (.not.do_sigmapringle) then
    write(iunit,"(a,es9.2,a,f5.1,a,f4.1,a,/)") '# Surface density      = ',&
         sig0*(Rref)**(-p_index)*umass/udist**2,' g/cm^2 (R/',Rref,')^(',-p_index,')'
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
!   H_R - aspect ratio of disc at R=R_in
!   q_index - power-law index of temperature profile
!   M_star - mass of star
!   R_in - inner edge of disc
!   i1,i2 - range of particle indices to use
!   rwarp - location of warp (optional)
!
!  Output:
!   honHmin, honH, honHmax - max mean and min ratio of h/H within range of R
!-----------------------------------------------------------------------------
subroutine get_honH(xyzh,rminav,rmaxav,honHmin,honHmax,honH,npart,H_R,q_index,M_star,R_in,i1,i2,verbose,rwarp)
 use domain,   only:i_belong
 use mpiutils, only:reduceall_mpi
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: honHmax,honHmin,honH
 integer, intent(in)  :: npart,i1,i2
 real,    intent(in)  :: H_R,q_index,M_star,rminav,rmaxav,R_in
 logical, intent(in)  :: verbose
 real,    intent(in), optional :: rwarp

 integer, parameter :: nr = 350
 integer :: i,ii,iwarp
 real :: G,rmin,rmax,dr,cs0,ri
 real :: rad(nr),ninbin(nr),h_smooth(nr),cs(nr),H(nr),omega(nr)
 integer :: ipart

 G = 1.0

! Setup rmin and rmax for the analysis
 rmin = rminav
 rmax = rmaxav

! Set up the radius array
 dr = (rmax-rmin)/real(nr-1)
 iwarp = 0
 do i=1,nr
    rad(i)=rmin + real(i-1)*dr
    if (present(rwarp)) then
       if (rad(i) > rwarp) iwarp = i - 1
    endif
 enddo

! Initialise arrays to zero
 ninbin(:)=0
 h_smooth(:)=0.0

! Set up cs0: cs = cs0 * R^-q
 cs0 = H_R * sqrt(G*M_star) * R_in**(q_index-0.5)

! And thus the sound speed array
 do i=1,nr
    cs(i) = cs_func(cs0,rad(i),q_index)
    omega(i) = sqrt(G*M_star/rad(i)**3)
 enddo

! and thus the disc scale height
 do i=1,nr
    H(i) = cs(i)/omega(i)
 enddo

! Loop over particles putting properties into the correct bin
 ipart = 0
 do i = i1,i2
    if (i_belong(i)) then
       ipart = ipart + 1
       if (xyzh(4,ipart)  >  tiny(xyzh)) then ! IF ACTIVE
          ri = sqrt(dot_product(xyzh(1:3,ipart),xyzh(1:3,ipart)))
          ii = int((ri-rad(1))/dr + 1)

          if (ii > nr) cycle
          if (ii < 1)  cycle

          ! Ignoring the large smoothing length particles far from the mid-plane
          if (xyzh(3,ipart)**2 < 4.*H(ii)*H(ii)) then
             h_smooth(ii) = h_smooth(ii) + xyzh(4,ipart)

             ninbin(ii) = ninbin(ii) + 1
          endif
       endif
    endif
 enddo

 h_smooth = reduceall_mpi('+', h_smooth)
 ninbin = reduceall_mpi('+', ninbin)

! Average h_smooth
 do i = 1,nr
    if (ninbin(i) > 0) h_smooth(i) = h_smooth(i)/ninbin(i)
 enddo

! Now loop over rings to calculate required quantities
 do i = 1, nr
    if (H(i) > 0.) then
       h_smooth(i) = h_smooth(i)/H(i)
    endif
 enddo

! Print out an average value for <h>/H and thus \alpha_SS/\alpha_AV
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

end module setdisc
