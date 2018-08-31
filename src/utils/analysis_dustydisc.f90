!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for dustydisc
!
!  REFERENCES: None
!
!  OWNER: Mark Hutchison
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, dust, infile_utils, io, leastsquares, options, part,
!    physcon, set_dust, solvelinearsystem, units
!+
!--------------------------------------------------------------------------
module analysis
 use dim,  only:maxdusttypes
 use dust, only:grainsizecgs
 use part, only:ndusttypes
 implicit none
 character(len=20), parameter, public :: analysistype = 'dustydisc'
 public :: do_analysis

 integer, parameter :: nr = 50
 integer, parameter :: numlabels = 23
 integer, parameter :: numarrays = 7
 integer, parameter :: maxlabels = numlabels + numarrays*(maxdusttypes-1)
 integer, parameter, public :: &
          iradius    = 1,  &
          irhog      = 2,  &
          isigmagas  = 3,  &
          ih_H       = 4,  &
          ilx        = 5,  &
          ily        = 6,  &
          ilz        = 7,  &
          itilt      = 8,  &
          itwist     = 9,  &
          ipsi       = 10, &
          iH_R_init  = 11, &
          iH_R       = 12, &
          iomega     = 13, &
          ivK        = 14, &
          ics        = 15, &
          ivrgas     = 16, &
          ! initial index for arrays
          ivrdust    = 17, &
          iSt        = 18 +   (maxdusttypes-1), &
          itstop     = 19 + 2*(maxdusttypes-1), &
          ivrsigma   = 20 + 3*(maxdusttypes-1), &
          irhod      = 21 + 4*(maxdusttypes-1), &
          isigmadust = 22 + 5*(maxdusttypes-1), &
          iHdust_R   = 23 + 6*(maxdusttypes-1), &
          ! ending index for arrays
          ivrdustend    = iSt-1,        &
          iStend        = itstop-1,     &
          itstopend     = ivrsigma-1,   &
          ivrsigmaend   = irhod-1,      &
          irhodend      = isigmadust-1, &
          isigmadustend = iHdust_R-1,   &
          iHdust_R_end  = maxlabels

 real :: twist(nr),twistprev(nr)
 character(len=80) :: label(maxlabels),labelfmt,datafmt

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use dim,          only:maxp
 use io,           only:fatal
 use physcon,      only:pi,jupiterm,years,au
 use part,         only:iphase,npartoftype,igas,idust,massoftype,labeltype,dustfrac,tstop, &
                        rhoh,maxphase,iamtype,xyzmh_ptmass,vxyz_ptmass,nptmass,deltav, &
                        isdead_or_accreted,graindens
 use options,      only:use_dustfrac,iexternalforce
 use units,        only:umass,udist,utime
 use leastsquares, only:fit_slope
 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: xyzh(:,:),vxyz(:,:)
 real,             intent(in) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=80)  :: basename
 character(len=80)  :: output
 character(len=80)  :: filename
 character(len=20)  :: discprefix
 integer, parameter :: nxn = 2*maxdusttypes
 integer :: i,j,k,ir,ii,ierr,iline,ninbin(nr),ninbindust(maxdusttypes,nr),iwarp,nptmassinit
 integer :: icutgas(nr),icutdust(maxdusttypes,nr),find_ii(1),irealvisc
 integer :: itype,lu,nondustcols,dustcols,numcols
 real, allocatable :: deltavsum(:,:),dustfracisuma(:)
 real :: err,errslope,erryint,slope
 real :: Sig0,rho0,hi,rhoi
 real :: vgassol(2,nr),vdustsol(2,maxdusttypes,nr)
 real :: flat_cut,scale_cut,eps_dust_cut
 real :: stan_dev(maxdusttypes,nr)
 real :: zeta(nr+2),dzetadr(nr),Pr(nr+1),dPrdr(nr)
 real :: St_mid(maxdusttypes,nr),St_from_tstop(maxdusttypes,nr)
 real :: rhogmid(nr),rhodmid(maxdusttypes,nr)
 real :: rhog(npart),rhod(maxdusttypes,npart),rhogbin(npartoftype(igas),nr),rhodbin(maxdusttypes,npartoftype(igas),nr)
 real :: vK(nr),etabin(npartoftype(igas),nr),meaneta(nr),nuvisc(nr),shearvisc,alphaAV
 real :: vrgasbin(npartoftype(igas),nr),vrdustbin(maxdusttypes,npartoftype(igas),nr)
 real :: meanvrgas(nr),meanvrdust(maxdusttypes,nr)
 real :: meanrhog(nr),meanrhod(maxdusttypes,nr)
 real :: vgas(3),vdust(3,maxdusttypes),vrgas(npart),vrdust(maxdusttypes,npart)
 real :: R_in,R_out,R_ref,R_warp
 real :: H_R_in,H_R_out,H_R_ref,p_index,q_index,M_star,M_disc,R_c,R_cdust
 real :: R_in_dust,R_out_dust,R_ref_dust,R_warp_dust
 real :: H_R_in_dust,H_R_out_dust,H_R_ref_dust,p_index_dust,M_star_dust,M_disc_dust
 real :: G,rmin,rmax,cs0,angx,angy,angz,ri,area !Hi_part
 real :: dreven,log_dr,drlog(nr),log_grid(nr+1),grid(nr+1)
 real :: angtot,Ltot,tilt,dtwist,Li(3)
 real :: rad(nr),Lx(nr),Ly(nr),Lz(nr),h_smooth(nr),sigmagas(nr),sigmadust(maxdusttypes,nr),cs(nr),H(nr),omega(nr)
 real :: zsetgas(npartoftype(igas),nr),hgas(nr),meanzgas(nr)
 real :: dustfraci_bin(maxdusttypes,npartoftype(igas),nr),meandustfraci(maxdusttypes,nr)
 real :: tstopbin(maxdusttypes,npartoftype(igas),nr),meantstop(maxdusttypes,nr)
 real :: dustfracisum_bin(npartoftype(igas),nr),meandustfracisum(nr)
 real :: d2g_ratio_bin(maxdusttypes,npartoftype(igas),nr),meand2g_ratio(maxdusttypes,nr)
 real :: unitlx(nr),unitly(nr),unitlz(nr),tp(nr)
 real :: zsetdust(maxdusttypes,max(npartoftype(idust),npartoftype(igas)),nr)
 real :: hdust(maxdusttypes,nr),meanzdust(maxdusttypes,nr)
 real :: psi_x,psi_y,psi_z,psi,Mdust,Mgas,Mtot,Macc,pmassi,pgasmass,pdustmass(maxdusttypes)
 real :: dustfraci(maxdusttypes),dustfracisum,rhoeff(maxdusttypes)
 real :: ri_mid
 real :: l_planet(3),bigl_planet,rad_planet,inc,planet_mass
 real, save :: Mtot_in,Mgas_in,Mdust_in
 logical, save :: init  = .false.

 integer :: scale_vel   = 0
 logical :: fit_sigma   = .false.
 logical :: fixslope    = .false.
 logical :: use_log_r   = .false.
 logical :: logdata     = .false.
 logical :: comparedata = .false.
 logical :: solve_baistone = .false.
 logical :: print_part_in_bin = .false.
 logical :: check_radial_migration = .false.

 integer, parameter :: iparams = 10
 integer, parameter :: iplanet = 23
 integer, parameter :: iprec   = 24
 integer, parameter :: isplash = 33
 integer, parameter :: isol    = 34
 integer, parameter :: iunit1  = 22
 logical :: do_precession,ifile

 if (use_dustfrac .and. ndusttypes > 1) then
    fit_sigma   = .false.
    fixslope    = .true.
    comparedata = .true.
    logdata     = .true.
    !use_log_r   = .true.
    scale_vel   = 2 ! 0=no scaling, 1=scale with meaneta, 2=scale with v_P
    check_radial_migration = .true.
    print_part_in_bin = .true.

    flat_cut     = 100.
    scale_cut    = 1.
    eps_dust_cut = 1.e-5
 else
    flat_cut     = huge(flat_cut)
    scale_cut    = huge(scale_cut)
    eps_dust_cut = tiny(eps_dust_cut)
 endif

 !--allocate memory for deltav calculations (only needed for one-fluid multigrain)
 if (.not. allocated(deltavsum)) allocate(deltavsum(3,npart))

 if (use_dustfrac) then
    allocate(dustfracisuma(npart))
    dustfracisuma(:) = sum(dustfrac(:,1:npart),1)
    where (dustfracisuma > 0.)
       deltavsum(1,:)   = sum(dustfrac(:,1:npart)*deltav(1,:,1:npart),1)/dustfracisuma
       deltavsum(2,:)   = sum(dustfrac(:,1:npart)*deltav(2,:,1:npart),1)/dustfracisuma
       deltavsum(3,:)   = sum(dustfrac(:,1:npart)*deltav(3,:,1:npart),1)/dustfracisuma
    elsewhere (dustfracisuma == 0.)
       deltavsum(1,:)   = 0.
       deltavsum(2,:)   = 0.
       deltavsum(3,:)   = 0.
    endwhere
    deallocate(dustfracisuma)
 else
    deltavsum = 0.
    deltav    = 0.
 endif

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 if (comparedata) then
    iline = index(dumpfile,'_')
    discprefix = dumpfile(1:iline-1)
    write(basename,"(a,i5.5)") trim(discprefix)//'_',numfile
    write(output,"(a)") trim(basename)//'.analysis'
 else
    write(output,"(a4,i5.5)") 'angm',numfile
 endif
 write(*,'("Output file name is ",A)') output

 if (use_dustfrac) then
    write(*,'("one-fluid model")')
 else
    write(*,'("two-fluid model")')
 endif

 iline = index(dumpfile,'_')
 discprefix = dumpfile(1:iline-1)
 write(filename,"(a)") trim(discprefix)//'.discparams'
 inquire(file=filename, exist=ifile)
 if (.not.ifile) write(filename,"(a)") 'discparams.list'
 call read_discparams(filename,R_in,R_out,R_ref,R_warp,H_R_in,H_R_out,H_R_ref, &
                           p_index,R_c,q_index,G,M_star,M_disc,iparams,ierr,cs0,Sig0)
 if (ierr /= 0) call fatal('analysis','could not open/read '//trim(filename))

 iline = index(dumpfile,'_')
 discprefix = dumpfile(1:iline-1)
 write(filename,"(a)") trim(discprefix)//'.in'
 inquire(file=filename, exist=ifile)
 if (ifile) then
    call read_in(filename,irealvisc,alphaAV,shearvisc,iunit,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read .in file')

    if (comparedata) then
       if (alphaAV == 0. .or. solve_baistone) then
          write(output,"(a)") 'baistone'//trim(adjustl(basename))//'.analysis'
       else
          if (irealvisc == 1) print*,'WARNING!!! Dipierro solution for irealvisc == 1 needs to be corrected'
          write(output,"(a)") 'dipierro'//trim(adjustl(basename))//'.analysis'
       endif
    endif
 endif

! Print out the parameters of gas disc
 write(*,*)
 write(*,'("Gas disc parameters are:")')
 write(*,*) 'R_in    = ',R_in
 write(*,*) 'R_out   = ',R_out
 write(*,*) 'R_ref   = ',R_ref
 write(*,*) 'H_R_in  = ',H_R_in
 write(*,*) 'H_R_out = ',H_R_out
 write(*,*) 'H_R_ref = ',H_R_ref
 write(*,*) 'Sig0    = ',Sig0
 write(*,*) 'cs0     = ',cs0*udist/utime
 if(R_warp /= 0.) &
 write(*,*) 'Rwarp   = ',R_warp
 write(*,*) 'p_index = ',p_index
 if(R_c /= 0.) &
 write(*,*) 'R_c     = ',R_c
 write(*,*) 'q_index = ',q_index
 write(*,*) 'G       = ',G
 write(*,*) 'M_star  = ',M_star
 write(*,*) 'M_disc  = ',M_disc
 write(*,*)
 write(*,*)

 write(filename,"(a)") trim(discprefix)//'-'//trim(labeltype(idust))//'.discparams'
 inquire(file=filename, exist=ifile)
 if (ifile) then
    call read_discparams(filename,R_in_dust,R_out_dust,R_ref_dust,R_warp_dust, &
                         H_R_in_dust,H_R_out_dust,H_R_ref_dust,p_index_dust,R_cdust, &
                         q_index,G,M_star_dust,M_disc_dust,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read -'//trim(filename))

    ! Print out the parameters of dust disc
    write(*,*)
    write(*,'("Dust disc parameters are:")')
    write(*,*) 'R_in    = ',R_in_dust
    write(*,*) 'R_out   = ',R_out_dust
    write(*,*) 'H_R_in  = ',H_R_in_dust
    write(*,*) 'H_R_out = ',H_R_out_dust
    write(*,*) 'H_R_ref = ',H_R_ref_dust
    if(R_warp_dust/=0.) &
    write(*,*) 'Rwarp   = ',R_warp_dust
    write(*,*) 'p_index = ',p_index_dust
    if(R_cdust/=0.) &
    write(*,*) 'R_c     = ',R_cdust
    write(*,*) 'G       = ',G
    write(*,*) 'M_disc  = ',M_disc_dust
    write(*,*)
    write(*,*)
 endif

! Setup rmin and rmax for the analysis
 if (npartoftype(idust) > 0) then
    rmin = min(R_in,R_in_dust)
    rmax = max(R_out,R_out_dust)
 else
    !rmin = R_in
    !rmax = R_out
    rmin = max(R_in,R_in_dust)
    rmax = min(R_out,R_out_dust)
 endif

! Set up the radius array
 if (use_log_r) then
    !--Create a uniform grid with N+1 points between rmin and rmax (inclusive)
    log_dr = log10(rmax/rmin)/real(nr)
    do i = 1,nr+1
       log_grid(i) = log10(rmin) + (i-1)*log_dr
    enddo

    !--Convert grid coordinates back to real space
    grid = 10.**log_grid

    !--Find representative point in each cell
    !  (skewed towards left edge...no real reason why except that is how it is done for grainsizes)
    do i = 1,nr
       rad(i) = sqrt(grid(i)*grid(i+1))
       drlog(i) = grid(i+1)-grid(i)
    enddo
 else
    ! Set up an evenly spaced radius array
    dreven = (rmax-rmin)/real(nr-1)
    do i=1,nr
       rad(i) = rmin + real(i-1)*dreven
    enddo
 endif
 rad = rad*au/udist

! Initialise arrays to zero
 icutgas(:)      = 0
 icutdust(:,:)   = 0

 ninbin(:)       = 0
 ninbindust(:,:) = 0

 lx(:) = 0.
 ly(:) = 0.
 lz(:) = 0.

 angx  = 0.
 angy  = 0.
 angz  = 0.

 H(:)      = 0.
 cs(:)     = 0.
 vK(:)     = 0.
 omega(:)  = 0.
 nuvisc(:) = 0.

 hgas(:)     = 0.
 hdust(:,:)  = 0.
 h_smooth(:) = 0.

 rhog(:)        = 0.
 rhod(:,:)      = 0.
 rhogmid(:)     = 0.
 rhodmid(:,:)   = 0.
 rhogbin(:,:)   = 0.
 rhodbin(:,:,:) = 0.

 sigmagas(:)     = 0.
 sigmadust(:,:)  = 0.

 tstopbin(:,:,:) = 0.

 zsetgas(:,:)    = 0.
 zsetdust(:,:,:) = 0.

 vrgasbin(:,:)    = 0.
 vrdustbin(:,:,:) = 0.

 meanrhog(:)      = 0.
 meanzgas(:)      = 0.
 meanvrgas(:)     = 0.
 meanrhod(:,:)    = 0.
 meantstop(:,:)   = 0.
 meanzdust(:,:)   = 0.
 meanvrdust(:,:)  = 0.
 stan_dev(:,:)    = 0.

 meandustfraci(:,:)  = 0.
 meandustfracisum(:) = 0.
 meand2g_ratio(:,:)  = 0.

 d2g_ratio_bin(:,:,:)  = 0.
 dustfraci_bin(:,:,:)  = 0.
 dustfracisum_bin(:,:) = 0.

! sound speed, orbital frequency/velocity, and gas scale-height
 do i = 1,nr
    cs(i)    = cs0*rad(i)**(-q_index)
    omega(i) = sqrt(G*M_star/rad(i)**3)
    vK(i)    = sqrt(G*M_star/rad(i))
    H(i)     = cs(i)/omega(i)
 enddo

 if (R_warp/=0.)then
    iwarp=3
 else
    iwarp=2
 endif

! if the central star is represented by a sink (change to 2 if you set a binary)
 nptmassinit = 1

! if the central star is represented by by external force
 if(iexternalforce /= 0) nptmassinit = 0

 if(nptmass > nptmassinit)then
    do i = nptmassinit+1,nptmass
       write(*,*)"Planet",i-nptmassinit,"mass",xyzmh_ptmass(4,i)*umass/jupiterm,"Jupiter mass, radius",&
     sqrt(dot_product(xyzmh_ptmass(1:iwarp,i),xyzmh_ptmass(1:iwarp,i)))*udist/au,"au"," coords xy ",&
     xyzmh_ptmass(1,i),xyzmh_ptmass(2,i)
    enddo
 endif

! Loop over gas particles putting properties into the correct bin
 Mtot   = 0.
 Mgas   = 0.
 Mdust  = 0.
 Macc   = 0.
 pmassi = massoftype(igas)
 rhoeff = graindens*sqrt(pi/8.)
 dustfracisum = 0.
 dustfraci(:) = 0.

 do i = 1,npart
    hi = xyzh(4,i)
    if (maxphase==maxp) then
       itype = iamtype(iphase(i))
       pmassi = massoftype(itype)
    endif
    Mtot = Mtot + pmassi
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (itype==igas) then
          if (use_dustfrac) then
             dustfraci = dustfrac(:,i)
             pdustmass = pmassi*dustfraci
          else
             dustfraci = 0.
             pdustmass = 0.
          endif
          dustfracisum = sum(dustfraci)
          pgasmass     = pmassi*(1. - dustfracisum)
          Mgas  = Mgas  + pgasmass
          Mdust = Mdust + pmassi*dustfracisum
       elseif (itype==idust) then
          Mdust = Mdust + pmassi
       endif
    else
       Macc  = Macc + pmassi
    endif

    !--Calculate vr for the gas and the N dust phases
    ri_mid = sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
    do j=1,ndusttypes
       if (use_dustfrac) then
          rhoi = rhoh(hi,pmassi)
          rhog(i)    = (1.-dustfracisum)*rhoi
          rhod(j,i)  = dustfraci(j)*rhoi
          if (ndusttypes > 1) then
             do k=1,ndusttypes
                if (isnan(deltavsum(k,i))) then
                   deltavsum(k,i) = 0.
                   deltav(:,k,i) = 0.
                endif
             enddo
             vgas(:)    = vxyz(:,i) - dustfracisum*deltavsum(:,i)
             vdust(:,j) = vxyz(:,i) + deltav(:,j,i) - dustfracisum*deltavsum(:,i)
          else
             vgas(:)    = vxyz(:,i) - dustfraci(j)*deltav(:,j,i)
             vdust(:,j) = vxyz(:,i) + (1. - dustfraci(j))*deltav(:,j,i)
          endif
       else
          select case(itype)
          case(igas)
             rhog(i)    = rhog(i)   + rhoh(hi,pmassi)
             vgas(:)    = vxyz(:,i)
             vdust(:,j) = 0.
          case(idust)
             rhod(j,i)  = rhod(j,i) + rhoh(hi,pmassi)
             vgas(:)    = 0.
             vdust(:,j) = vxyz(:,i)
          case default
             vgas(:)    = 0.
             vdust(:,j) = 0.
          end select
       endif
       vrgas(i)    = 0.5*(2*xyzh(1,i)*vgas(1)    + 2*xyzh(2,i)*vgas(2))   /ri_mid
       vrdust(j,i) = 0.5*(2*xyzh(1,i)*vdust(1,j) + 2*xyzh(2,i)*vdust(2,j))/ri_mid
    enddo

    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       ri = sqrt(dot_product(xyzh(1:iwarp,i),xyzh(1:iwarp,i)))
       !Hi_part = cs0*ri**(-q_index)/(sqrt(G*M_star/ri**3))
       if (use_log_r) then
          find_ii = minloc(ri-rad(:),ri-rad(:) > 0.)
          ii = find_ii(1)
       else
          ii = int((ri-rad(1))/dreven + 1)
       endif

       if (ii > nr) cycle
       if (ii < 1)  cycle

       if (use_log_r) then
          area = (pi*(grid(ii+1)**2-grid(ii)**2))
       else
          area = (pi*((rad(ii) + dreven/2.)**2 - (rad(ii) - dreven/2.)**2))
       endif

       if(iphase(i) == igas) then

          sigmagas(ii) = sigmagas(ii) + pgasmass/area
          if (use_dustfrac) sigmadust(:,ii) = sigmadust(:,ii) + pdustmass(:)/area

          Li(1)  = pgasmass*(xyzh(2,i)*vxyz(3,i)-xyzh(3,i)*vxyz(2,i))
          Li(2)  = pgasmass*(xyzh(3,i)*vxyz(1,i)-xyzh(1,i)*vxyz(3,i))
          Li(3)  = pgasmass*(xyzh(1,i)*vxyz(2,i)-xyzh(2,i)*vxyz(1,i))

          Lx(ii) = Lx(ii) + Li(1)
          Ly(ii) = Ly(ii) + Li(2)
          Lz(ii) = Lz(ii) + Li(3)

          h_smooth(ii) = h_smooth(ii) + xyzh(4,i)

          ninbin(ii) = ninbin(ii) + 1
          zsetgas(ninbin(ii),ii) = xyzh(3,i)

          etabin(ninbin(ii),ii)  = 0.5*(cs(ii)/vK(ii))**2*(p_index + q_index + 1.5 - &
                                   (1.5 - q_index)*(xyzh(3,i)/H(ii))**2)
          !etabin(ninbin(ii),ii)  = 0.5*(Hi_part/ri)**2*(p_index + q_index + 1.5 - &
          !                         (1.5 - q_index)*(xyzh(3,i)/Hi_part)**2)

          if (abs(xyzh(3,i)) < min(flat_cut,scale_cut*H(ii))) then
             rhogbin(ninbin(ii),ii)  = rhog(i)
             vrgasbin(ninbin(ii),ii) = vrgas(i)
             if (use_dustfrac) then
                vrdustbin(:,ninbin(ii),ii) = vrdust(:,i)
                tstopbin(:,ninbin(ii),ii)  = tstop(:,i)

                dustfracisum_bin(ninbin(ii),ii) = dustfracisum
                dustfraci_bin(:,ninbin(ii),ii)  = dustfraci(:)
                d2g_ratio_bin(:,ninbin(ii),ii)  = dustfraci(:)/(1. - dustfracisum)
                do j=1,ndusttypes
                   if (use_dustfrac .and. dustfraci(j) > eps_dust_cut) then
                      ninbindust(j,ii) = ninbindust(j,ii) + 1
                      rhodbin(j,ninbindust(j,ii),ii)  = rhod(j,i)
                      zsetdust(j,ninbindust(j,ii),ii) = xyzh(3,i)
                   else
                      icutdust(j,ii) = icutdust(j,ii) + 1
                   endif
                enddo
             endif
          else
             icutgas(ii) = icutgas(ii) + 1
          endif
       elseif(iphase(i) == idust) then
          do j=1,ndusttypes
             sigmadust(j,ii) = sigmadust(j,ii) + pdustmass(j)/area

             ninbindust(j,ii) = ninbindust(j,ii) + 1
             if (abs(xyzh(3,i)) < min(flat_cut,scale_cut*H(ii))) then
                ! Important for two-fluid method:
                !  <>  dustfrac = dust-to-gas ratio for iphase = igas
                !  <>  dustfrac = gas-to-dust ratio for iphase = idust
                d2g_ratio_bin(j,ninbindust(j,ii),ii) = 1./dustfrac(j,i)
                tstopbin(j,ninbindust(j,ii),ii)  = tstop(j,i)
                rhodbin(j,ninbindust(j,ii),ii)   = rhod(j,i)
                vrdustbin(j,ninbindust(j,ii),ii) = vrdust(j,i)
                zsetdust(j,ninbindust(j,ii),ii)  = xyzh(3,i)
             else
                icutdust(j,ii) = icutdust(j,ii) + 1
             endif
          enddo
       endif
    elseif (xyzh(4,i) < -tiny(xyzh) .and. iphase(i)==igas) then !ACCRETED
       angx = angx + pgasmass*(xyzh(2,i)*vxyz(3,i) - xyzh(3,i)*vxyz(2,i))
       angy = angy + pgasmass*(xyzh(3,i)*vxyz(1,i) - xyzh(1,i)*vxyz(3,i))
       angz = angz + pgasmass*(xyzh(1,i)*vxyz(2,i) - xyzh(2,i)*vxyz(1,i))
    endif
 enddo
 write(*,*)"Dust mass: ",Mdust
 write(*,*)"Gas mass: ",Mgas

 numcols = 5 ! # of total columns
 if (.not.init) then
    open(newunit=lu,file='dustmass.ev',status='replace')
    write(labelfmt,'(I10)') numcols
    write(labelfmt,'(A)') '(''#'','//trim(adjustl(labelfmt))//'(1x,''['',i2.2,1x,a12,'']'',1x))'
    write(lu,labelfmt) 1,'time',2,'Mtot',3,'Mgas',4,'Mdust',5,'Macc'
    init = .true.
    Mgas_in  = Mgas
    Mdust_in = Mdust
    Mtot_in  = Mtot
 else
    open(newunit=lu,file='dustmass.ev',status='old',position='append')
 endif

 print*,' Mtot = ',Mtot,' Mgas = ',Mgas,' Mdust = ',Mdust,' Macc = ',Macc
 print "(4(/,a,2pf6.2,'%'))",' dMtot  = ',(Mtot-Mtot_in)/Mtot_in,&
       ' dMgas  = ',(Mgas-Mgas_in)/Mtot_in,&
       ' dMdust = ',(Mdust-Mdust_in)/Mtot_in,&
       ' dMacc  = ',Macc/Mtot_in
 write(datafmt,'(I10)') numcols
 write(datafmt,'(A)') '('//trim(adjustl(datafmt))//'(es18.10,1X))'
 write(lu,datafmt) time,Mtot,Mgas,Mdust,Macc
 close(lu)

 ! Computing Hgas, Hdust, vrgas, vrdust
 do i = 1,nr
    meaneta(i)  = sum(etabin(1:ninbin(i),i))/real(ninbin(i))
    if(ninbin(i) > 1)then
       if (print_part_in_bin) print*,'The # of particles in bin',i,'is',ninbin(i)-icutgas(i)
       meanzgas(i)  = sum(zsetgas(1:ninbin(i),i))/real(ninbin(i))
       meanrhog(i)  = sum(rhogbin(1:ninbin(i),i))/(real(ninbin(i))-icutgas(i))
       meanvrgas(i) = sum(vrgasbin(1:ninbin(i),i))/(real(ninbin(i))-icutgas(i))
       if (use_dustfrac) then
          do j=1,ndusttypes
             meanvrdust(j,i) = sum(vrdustbin(j,1:ninbin(i),i))/(real(ninbin(i))-icutgas(i))
             stan_dev(j,i)   = sqrt(sum((vrdustbin(j,1:ninbin(i),i)-meanvrdust(j,i))**2) &
                               /(real(ninbin(i)) - icutgas(i) - 1))
             meandustfraci(j,i) = sum(dustfraci_bin(j,1:ninbin(i),i))/(real(ninbin(i))-icutgas(i))
             meand2g_ratio(j,i) = sum(d2g_ratio_bin(j,1:ninbin(i),i))/(real(ninbin(i))-icutgas(i))
             meantstop(j,i) = sum(tstopbin(j,1:ninbin(i),i))/(real(ninbin(i))-icutgas(i))
          enddo
          meandustfracisum(i) = sum(dustfracisum_bin(1:ninbin(i),i))/(real(ninbin(i))-icutgas(i))
       endif
       hgas(i) = sqrt(sum(((zsetgas(1:ninbin(i),i)-meanzgas(i))**2)/(real(ninbin(i)-1))))
    endif
    do j=1,ndusttypes
       if (ninbindust(j,i) > 1) then
          meanzdust(j,i)  = sum(zsetdust(j,1:ninbindust(j,i),i))/(real(ninbindust(j,i)))
          meanrhod(j,i)   = sum(rhodbin(j,1:ninbindust(j,i),i))/(real(ninbindust(j,i))-icutdust(j,i))
          hdust(j,i)      = sqrt(sum(((zsetdust(j,1:ninbindust(j,i),i)-meanzdust(j,i))**2)/(real(ninbindust(j,i)-1))))
          if (.not. use_dustfrac) then
             meanvrdust(j,i)    = sum(vrdustbin(j,1:ninbindust(j,i),i))/(real(ninbindust(j,i))-icutdust(j,i))
             stan_dev(j,i)   = sqrt(sum((vrdustbin(j,1:ninbindust(j,i),i)-meanvrdust(j,i))**2) &
                               /(real(ninbindust(j,i)) - icutdust(j,i) - 1))
             meand2g_ratio(j,i) = sum(d2g_ratio_bin(j,1:ninbindust(j,i),i))/(real(ninbindust(j,i))-icutdust(j,i))
             meantstop(j,i)     = sum(tstopbin(j,1:ninbindust(j,i),i))/(real(ninbindust(j,i))-icutdust(j,i))
          endif
       endif
    enddo
 enddo

 ! Print angular momentum of accreted particles
 angtot = sqrt(angx*angx + angy*angy + angz*angz)
 print*,' angular momentum of accreted particles = ',angtot

! Convert total angular momentum into a unit vector, and average h_smooth
 do i = 1,nr
    Ltot = sqrt(Lx(i)*Lx(i) + Ly(i)*Ly(i) + Lz(i)*Lz(i))

    if(Ltot/=0.) then
       unitlx(i) = Lx(i)/Ltot
       unitly(i) = Ly(i)/Ltot
       unitlz(i) = Lz(i)/Ltot
    else
       unitlx(i) = 0.
       unitly(i) = 0.
       unitlz(i) = 0.
    endif

    if (ninbin(i) > 0) h_smooth(i) = h_smooth(i)/ninbin(i)
 enddo

! Now loop over rings to calculate required quantities
 do i = 1, nr
    if(ninbin(i) == 0) then
       lx(i) = 0.
       ly(i) = 0.
       lz(i) = 0.
       sigmagas(i) = 0.
       h_smooth(i) = 0.
    else
       h_smooth(i) = h_smooth(i)/H(i)
    endif
 enddo

 print*,' '
 print*,'Warning: removed a total of',sum(icutgas(:)),'gas particles from the velocity'
 do j=1,ndusttypes
    print*,'             and a total of',sum(icutdust(j,:)),'dust particles from the velocity'
 enddo
 print*,' '

 !--if Sig0 does not match the binned sigma calculated above so we can do a least squares
 if (fit_sigma) then
    print*,' '
    print*,'Fitting the power-law surface density profile for the analytic solution...'
    slope = -p_index
    call fit_slope(nr,rad,sigmagas,slope,Sig0,err,errslope,erryint,R_in_dust,R_out_dust,logdata,fixslope)
    if (logdata) Sig0 = 10.**Sig0
    rho0 = Sig0/(sqrt(2.*pi)*H_R_ref)
    print*,'   Corrected value for Sig0 = ',Sig0*umass/udist**2
    print*,'   Corrected value for rho0 = ',rho0
    print*,'...Done'
    print*,' '

    rhogmid(:) = rho0*rad(:)**(-(p_index-q_index+1.5))
    do j=1,ndusttypes
       rhodmid(j,:) = rhogmid(:)*meandustfraci(j,:)
    enddo
    !sigmagas(i) = Sig0*rad(i)**(-p_index)
    sigmagas(:) = Sig0*(rad(:)/R_ref)**(-p_index)*(1.-sqrt(R_in/rad(:)))
 else
    rhogmid  = meanrhog
    rhodmid  = meanrhod
    sigmagas = sigmagas(:)
 endif

! and thus the Stokes Number
! St = tstop*omegaK = rhoeff*s/(cs*rhogmid)*omegaK
!                   = rhoeff*s*sqrt(2*pi)/(sigma*exp(-0.5*(z/Hg)**2)))
!                   = rhoeff*s*sqrt(2*pi)/sigma    ! if z = 0
 do i = 1,nr
    ! Note that for the Dipierro et al. (2018) solution we ended up scaling St to represent
    ! the values that were closer to the mid-plane. However, we did not need this for the
    ! Bai & Stone (2010) test. This could indicate that there is a bug somewhere in this analysis
    ! tool or in the code itself.
    St_from_tstop(:,i) = meantstop(:,i)*omega(i)
    St_mid(:,i) = St_from_tstop(:,i)

!    St_mid(:,i) = rhoeff(:)*grainsize(:)*(sqrt(2.*pi)/sigmagas(i)) ! St_mid = Lambda_k
!    meantstop(:,i) = St_mid(:,i)/omega(i)
    if (irealvisc == 0) then
       nuvisc(i) = (0.1*alphaAV*h_smooth(i))*cs(i)*H(i)
    elseif (irealvisc == 1) then
       nuvisc(i) = shearvisc
    elseif (irealvisc == 2) then
       nuvisc(i) = shearvisc*cs(i)*H(i)
    endif

    ! Variables needed to calculate v_nu and v_P in an evolving disc
    zeta(i+1) = nuvisc(i)*rad(i)*rhogmid(i)*vK(i)
    Pr(i+1) = cs(i)**2*rhogmid(i)
 enddo

 ! Extrapolate out to populate the ghost cells
 zeta(1) = 3.*zeta(2) -3.*zeta(3) + zeta(4)
 zeta(nr+2) = zeta(nr-1) - 3.*zeta(nr)  + 3.*zeta(nr+1)
 Pr(1) = 3.*Pr(2) -3.*Pr(3) + Pr(4)
 Pr(nr+2) = Pr(nr-1) - 3.*Pr(nr)  + 3.*Pr(nr+1)

 ! Calculate the derivatives
 do i = 2,nr+1
    dzetadr(i-1) = (zeta(i+1) - zeta(i-1))/(2.*dreven)
    dPrdr(i-1) = (Pr(i+1) - Pr(i-1))/(2.*dreven)
    dPrdr(i-1) = dPrdr(i-1)/(rhogmid(i-1)*omega(i-1))
 enddo

 ! Scale the velocities
 if (scale_vel == 1) then
    do i = 1,nr
       meanvrgas(i)    = meanvrgas(i)/(meaneta(i)*vK(i))
       meanvrdust(:,i) = meanvrdust(:,i)/(meaneta(i)*vK(i))
       stan_dev(:,i)   = stan_dev(:,i)/(meaneta(i)*vK(i))
    enddo
 elseif (scale_vel == 2) then
    do i = 1,nr
       meanvrgas(i)    = meanvrgas(i)/abs(dPrdr(i))
       meanvrdust(:,i) = meanvrdust(:,i)/abs(dPrdr(i))
       stan_dev(:,i)   = stan_dev(:,i)/abs(dPrdr(i))
    enddo
 endif

 ! Make labels for header of output file
 label(iradius)    = 'radius'
 label(irhog)      = 'rhogmid'
 label(irhod)      = 'rhodmid'
 label(isigmagas)  = 'sigmagas'
 label(isigmadust) = 'sigmadust'
 label(ih_H)       = '<h>/H'
 label(ilx)        = 'lx'
 label(ily)        = 'ly'
 label(ilz)        = 'lz'
 label(itilt)      = 'tilt'
 label(itwist)     = 'twist'
 label(ipsi)       = 'psi'
 label(iH_R_init)  = 'H/R init'
 label(iH_R)       = 'H/R'
 label(iHdust_R)   = 'Hdust/R'
 label(iomega)     = 'OmegaK'
 label(ivK)        = 'vK'
 label(ics)        = 'cs'
 label(ivrgas)     = '<vrgas>'

 call make_output_labels(ivrdust,   ivrdustend,   '<vrdust','>')
 call make_output_labels(iSt,       iStend,       'St',     '' )
 call make_output_labels(itstop,    itstopend,    'tstop',  '' )
 call make_output_labels(ivrsigma,  ivrsigmaend,  'vrsigma','' )
 call make_output_labels(irhod,     irhodend,     'rhodmid','' )
 call make_output_labels(isigmadust,isigmadustend,'sigmad', '' )
 call make_output_labels(iHdust_R,  iHdust_R_end, 'Hdust/R','' )

 ! Make splash.columns file using the labels above
 write(filename,"(a)") 'dustanal.columns'
 open(unit=isplash,file=filename,status="replace")
 do i = 1,maxlabels
    write(isplash,'(A)') trim(adjustl(label(i)))
 enddo
 close(unit=isplash)

 ! Make format for write statements
 write(labelfmt,'(I10)') maxlabels
 write(labelfmt,'(A)') '(''#'','//trim(adjustl(labelfmt))//'(1x,''['',i2.2,1x,a11,'']'',2x))'
 write(datafmt,'(I10)') maxlabels
 write(datafmt,'(A)') '('//trim(adjustl(datafmt))//'(es18.10,1X))'

 ! Write the header to the file
 open(iunit,file=output)
 write(iunit,'("# time:             time unit ( yrs)")')
 write(iunit,'("# ",es20.12,es20.12)') time,utime
 write(iunit,labelfmt) (j,trim(adjustl(label(j))),j=1,maxlabels)

 do_precession = .false.

! Write the data to the file
 do i=1,nr
    if(i /= 1.and.i /= nr) then
       psi_x = (unitlx(i+1)-unitlx(i-1))/(rad(i+1)-rad(i-1))
       psi_y = (unitly(i+1)-unitly(i-1))/(rad(i+1)-rad(i-1))
       psi_z = (unitlz(i+1)-unitlz(i-1))/(rad(i+1)-rad(i-1))
       psi   = sqrt(psi_x**2 + psi_y**2 + psi_z**2)*rad(i)
    else
       psi = 0.
    endif

    if (ninbin(i) > 0) then
       tilt     = acos(unitlz(i))
       twist(i) = atan2(unitly(i),unitlx(i))
       if (i == 1 .or. time == 0.) then
          twistprev(i) = 0.
       endif
       ! Taking into account negative twist
       if (twist(i) < 0) then
          twistprev(i) = 2.*pi + twist(i)
       else
          twistprev(i) = twist(i) !cumulative twist
       endif
    else
       tilt   = 0.
       twist  = 0.
       dtwist = 0.
    endif

! Calculate the precession time
    if (twist(i) > tiny(twist(i))) then
       tp(i) = time*2.*pi/twist(i)
    else
       tp(i) = 0.
    endif

    write(iunit,datafmt)  &
          rad(i),         &
          rhogmid(i),     &
          sigmagas(i),    &
          h_smooth(i),    &
          unitlx(i),      &
          unitly(i),      &
          unitlz(i),      &
          tilt,           &
          twist(i),       &
          psi,            &
          H(i)/rad(i),    &
          hgas(i)/rad(i), &
          omega(i),       &
          vK(i),          &
          cs(i),          &
          meanvrgas(i),   &
          (meanvrdust(j,i),  j=1,ndusttypes), &
          (St_mid(j,i),      j=1,ndusttypes), &
          (meantstop(j,i),   j=1,ndusttypes), &
          (stan_dev(j,i),    j=1,ndusttypes), &
          (rhodmid(j,i),     j=1,ndusttypes), &
          (sigmadust(j,i),   j=1,ndusttypes), &
          (Hdust(j,i)/rad(i),j=1,ndusttypes)

! Printing time and twist for each radius bin
    if (do_precession) then
       write(filename,"(a,i3.3)")"precess",i
       if (numfile==0) then
          open(unit=iprec,file=filename,status="replace")
          write(iprec,'("# tilt and twist with time for r = ",es18.10)') rad(i)
          write(iprec,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
               1,'rad', &
               2,'time', &
               3,'tilt', &
               4,'twist', &
               5,'tot twist', &
               6,'tp'
       else
          open(unit=iprec,file=filename,status="old",position="append")
       endif
       write(iprec,'(6(es18.10,1X))') rad(i),time,tilt,twist(i),twistprev(i),tp(i)
       close(unit=iprec)
    endif

 enddo

 close(iunit)

 if (check_radial_migration) then
    print *,' '
    if (alphaAV == 0. .or. solve_baistone) then
       print *,'Solving the Bai and Stone (2010) analytic solution for radial drift...'
       call solve_bai_stone_2010(meand2g_ratio,nxn,meaneta,vK,vgassol,vdustsol,St_mid)
    else
       print *,'Solving the Dipierro et al. (2018) analytic solution for radial drift...'
       call solve_dipierro_2018(irealvisc,vgassol,vdustsol,meand2g_ratio,rad,cs,vK, &
                                 nuvisc,p_index,q_index,St_mid,zeta,dzetadr,dPrdr)
    endif

    if (scale_vel == 1) then
       do i = 1,nr
          vgassol(:,i)    = vgassol(:,i)/(meaneta(i)*vK(i))
          vdustsol(:,:,i) = vdustsol(:,:,i)/(meaneta(i)*vK(i))
       enddo
    elseif (scale_vel == 2) then
       do i = 1,nr
          vgassol(:,i)    = vgassol(:,i)/abs(dPrdr(i))
          vdustsol(:,:,i) = vdustsol(:,:,i)/abs(dPrdr(i))
       enddo
    endif

    print*,'...Done'
    print *,' '

    ! Print solution to file
    write(filename,"(a)") trim(output)//'.sol'
    open(unit=isol,file=filename,status="replace")
    nondustcols = 3                  ! # of non-dust columns
    dustcols = iStend - ivrdust + 1  ! # of dust columns
    numcols = nondustcols + dustcols ! # of total columns
    write(labelfmt,'(I10)') numcols
    write(labelfmt,'(A)') '(''#'','//trim(adjustl(labelfmt))//'(1x,''['',i2.2,1x,a11,'']'',2x))'
    write(datafmt,'(I10)') numcols
    write(datafmt,'(A)') '('//trim(adjustl(datafmt))//'(es18.10,1X))'
    write(isol,labelfmt) 1,trim(adjustl(label(iradius))), &
                         2,trim(adjustl(label(isigmagas))), &
                         3,trim(adjustl(label(ivrgas))), &
                        (4 + j-ivrdust,trim(adjustl(label(j))),j=ivrdust,iStend)
    do ir = 1,nr
       if (sigmagas(ir) == 0.) then
          vgassol(1,ir) = 0.
          vdustsol(1,:,ir) = 0.
          St_mid(:,ir) = 0.
       endif
       write(isol,datafmt)  &
             rad(ir),       &
             sigmagas(ir),  &
             vgassol(1,ir), &
             (vdustsol(1,j,ir), j=1,ndusttypes), &
             (St_mid(j,ir),j=1,ndusttypes)
    enddo
    close(unit=isol)

    ! Make splash.columns file using the labels above
    write(filename,"(a)") 'sol.columns'
    open(unit=isplash,file=filename,status="replace")
    write(isplash,'(A)') trim(adjustl(label(iradius)))
    write(isplash,'(A)') trim(adjustl(label(isigmagas)))
    do i = ivrdust,ivrdustend
       write(isplash,'(A)') trim(adjustl(label(i)))
    enddo
    do i = iSt,iStend
       write(isplash,'(A)') trim(adjustl(label(i)))
    enddo
    close(unit=isplash)
 endif

 write(*,*)
 write(*,*) '--------------------------------------------------------------------------------'

! Printing the information for the planet as a function of time
 if(nptmass>nptmassinit)then
    do i=nptmassinit+1,nptmass
       write(filename,"(a,i3.3)")"planet_",i-1
       if (numfile==0) then
          open(iplanet,file=filename,status="replace")
          write(iplanet,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
               1,'time', &
               2,'x', &
               3,'y', &
               4,'z', &
               5,'lx', &
               6,'ly', &
               7,'lz', &
               8,'tilt', &
               9,'rad_sph'
       else
          open(iplanet,file=filename,status="old",position="append")
       endif
       planet_mass = xyzmh_ptmass(4,i)
       rad_planet = sqrt((xyzmh_ptmass(1,i)-xyzmh_ptmass(1,1))**2 + &
                    (xyzmh_ptmass(2,i)-xyzmh_ptmass(2,1))**2 + (xyzmh_ptmass(3,i)-xyzmh_ptmass(3,1))**2)
       l_planet(1) = planet_mass*((xyzmh_ptmass(2,i)*vxyz_ptmass(3,i)) - (xyzmh_ptmass(3,i)*vxyz_ptmass(2,i)))
       l_planet(2) = planet_mass*((xyzmh_ptmass(3,i)*vxyz_ptmass(1,i)) - (xyzmh_ptmass(1,i)*vxyz_ptmass(3,i)))
       l_planet(3) = planet_mass*((xyzmh_ptmass(1,i)*vxyz_ptmass(2,i)) - (xyzmh_ptmass(2,i)*vxyz_ptmass(1,i)))
       bigl_planet = sqrt(dot_product(l_planet,l_planet))
       l_planet = l_planet/bigl_planet
       inc = acos(abs(l_planet(3)))*180./pi
       write(iplanet,'(9(es18.10,1X))') time, xyzmh_ptmass(1,i), xyzmh_ptmass(2,i), xyzmh_ptmass(3,i),&
            l_planet(1),l_planet(2),l_planet(3),inc,rad_planet
       close(iplanet)
    enddo
 endif

 if (allocated(deltavsum)) deallocate(deltavsum)

 return
end subroutine do_analysis


!----------------------------------------------------------------
!+
!  Solve for the radial velocities in Bai & Stone 2010b test
!+
!----------------------------------------------------------------
subroutine solve_bai_stone_2010(d2g_ratio,nxn,eta,vK,vgassol,vdustsol,St_mid)
 use solvelinearsystem, only:rowk,dple
 integer, intent(in)  :: nxn
 real,    intent(in)  :: d2g_ratio(:,:),eta(:),vK(:),St_mid(:,:)
 real,    intent(out) :: vgassol(:,:),vdustsol(:,:,:)
 integer, parameter   :: dp = selected_real_kind(14, 60)
 integer  :: i,ir,ierr
 real(dp) :: soln(nxn)
 real :: Imat(maxdusttypes,maxdusttypes)
 real :: Lambda(maxdusttypes,maxdusttypes)
 real :: Gamma(maxdusttypes,maxdusttypes)
 real :: Amat(nxn,nxn),Bmat(nxn)

 vdustsol(:,:,:) = 0.
 do ir = 1,nr
    Bmat(1:ndusttypes) = 0.
    Bmat(ndusttypes+1:nxn) = 1.

    Imat(:,:) = 0.
    do i=1,ndusttypes
       Imat(i,i) = 1.
    enddo

    Lambda(:,:) = 0.
    do i=1,ndusttypes
       Lambda(i,i) = St_mid(i,ir)
       Gamma(i,:)  = d2g_ratio(:,ir)
    enddo

    Amat(:,:) = 0.
    do i=1,ndusttypes
       Amat(i,           :) = [ (Imat(i,:)+Gamma(i,:)), -2.*Lambda(i,:)        ]
       Amat(i+ndusttypes,:) = [ 0.5*Lambda(i,:)       , (Imat(i,:)+Gamma(i,:)) ]
    enddo

    call dple(rowk,nxn,Amat,Bmat,soln,ierr)
    if (ierr /= 0) then
       write(*, *) ' Error = ', ierr
    endif

    ! Save the solution to an array
    vdustsol(1,:,ir) = -eta(ir)*vK(ir)*soln(1:ndusttypes)
    vdustsol(2,:,ir) = -eta(ir)*vK(ir)*soln(ndusttypes+1:nxn)

    vgassol(1,ir) = -sum(d2g_ratio(:,ir)*vdustsol(1,:,ir))
    vgassol(2,ir) = -sum(d2g_ratio(:,ir)*vdustsol(2,:,ir)) - eta(ir)*vK(ir)
 enddo

 return
end subroutine solve_bai_stone_2010


!----------------------------------------------------------------
!+
!  Solve for the radial velocities in Dipierro et al. 2018 test
!+
!----------------------------------------------------------------
subroutine solve_dipierro_2018(irealvisc,vgassol,vdustsol,d2g_ratio,r,cs,vK,nu,p,q, &
                               St_mid,zeta,dzetadr,dPrdr)

 integer, intent(in)  :: irealvisc
 real,    intent(out) :: vgassol(:,:)
 real,    intent(out) :: vdustsol(:,:,:)
 real,    intent(in)  :: d2g_ratio(:,:)
 real,    intent(in)  :: r(:),cs(:),vK(:),nu(:)
 real,    intent(in)  :: p,q
 real,    intent(in)  :: St_mid(:,:),zeta(:),dzetadr(:),dPrdr(:)

 integer :: i
 real    :: denom2(maxdusttypes)
 real    :: lambda0,lambda1,v_P,v_nu
 real    :: denom1

 do i = 1,nr
    lambda0 = sum(d2g_ratio(:,i)/(1. + St_mid(:,i)**2))
    lambda1 = sum(d2g_ratio(:,i)*St_mid(:,i)/(1. + St_mid(:,i)**2))

    if (irealvisc == 0) then
       v_P = dPrdr(i)
       v_nu = -3.*nu(i)*dzetadr(i)/zeta(i+1)
       !v_nu = nu(i)*(2.*p + q + 1.5)/r(i)
    elseif (irealvisc == 1) then
       v_P = -(p + q + 1.5)*cs(i)**2/vK(i)
       v_nu = 3.*nu(i)*(p - q + 1.)/r(i)
    elseif (irealvisc == 2) then
       v_P = -(p + q + 1.5)*cs(i)**2/vK(i)
       v_nu = 3.*nu(i)*(p + q - 0.5)/r(i)
    endif

    denom1          = (1. + lambda0)**2 + lambda1**2
    denom2(:)       = denom1*(1. + St_mid(:,i)**2)
    vgassol(1,i)    = (-lambda1*v_P + (1. + lambda0)*v_nu)/denom1
    vdustsol(1,:,i) = (v_P*((1. + lambda0)*St_mid(:,i) - lambda1) +  &
                       v_nu*(1. + lambda0 + St_mid(:,i)*lambda1))/denom2(:)
    vgassol(2,i)    = 0.5*(v_P*(1. + lambda0) + v_nu*lambda1)/denom1
    vdustsol(2,:,i) = 0.5*(v_P*(1. + lambda0 + St_mid(:,i)*lambda1) -  &
                       v_nu*((1. + lambda0)*St_mid(:,i) - lambda1))/denom2(:)
 enddo

 return
end subroutine solve_dipierro_2018


!----------------------------------------------------------------
!+
!  Read disc information from discparams.list file
!+
!----------------------------------------------------------------
subroutine read_discparams(filename,R_in,R_out,R_ref,R_warp,H_R_in,H_R_out,H_R_ref, &
                           p_index,R_c,q_index,G,M_star,M_disc,iunit,ierr,cs0,Sig0)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: R_in,R_out,R_ref,R_warp,R_c
 real,             intent(out) :: H_R_in,H_R_out,H_R_ref
 real,             intent(out) :: p_index,q_index
 real,             intent(out) :: G,M_star,M_disc
 real, optional,   intent(out) :: cs0,Sig0
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 real :: sig_in,sig_ref,sig_out,sig_max
 type(inopts), allocatable :: db(:)

! Read in parameters from the file discparams.list
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(R_in,'R_in',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_out,'R_out',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_ref,'R_ref',db,ierr)
 call read_inopt(R_warp,'R_warp',db,ierr)
 call read_inopt(H_R_in,'H/R_in',db,ierr)
 call read_inopt(H_R_out,'H/R_out',db,ierr)
 call read_inopt(H_R_ref,'H/R_ref',db,ierr)
 call read_inopt(sig_max,'sig_max',db,ierr)
 if (ierr /= 0) return
 call read_inopt(p_index,'p_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_c,'R_c',db,ierr)
 call read_inopt(q_index,'q_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_star,'M_star',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_disc,'M_disc',db,ierr)
 if (ierr /= 0) return
 call read_inopt(G,'G',db,ierr)
 if (ierr /= 0) then
    ! Assuming G=1
    write(*,'("cannot find G......assuming G==1")')
    G = 1.0
 endif
 if (present(cs0)) then
    call read_inopt(cs0,'cs0',db,ierr)
    if (ierr /= 0) then
       cs0 = H_R_ref*sqrt(G*M_star/R_ref)*R_ref**q_index
    endif
 endif
 if (present(Sig0)) then
    call read_inopt(sig_in,'sig_in',db,ierr)
    call read_inopt(sig_ref,'sig_ref',db,ierr)
    call read_inopt(sig_out,'sig_out',db,ierr)
    Sig0 = sig_ref
 endif

 call close_db(db)

 return
end subroutine read_discparams


!----------------------------------------------------------------
!+
!  read disc information from .in file
!+
!----------------------------------------------------------------
subroutine read_in(filename,irealvisc,alphaAV,shearvisc,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: alphaAV,shearvisc
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr,irealvisc
 type(inopts), allocatable :: db(:)

! read in parameters from the .setup file
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(alphaAV,'alpha',db,ierr)
 if (ierr /= 0) return
 call read_inopt(irealvisc,'irealvisc',db,ierr)
 if (ierr /= 0) return
 call read_inopt(shearvisc,'shearparam',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_in


!----------------------------------------------------------------
!+
!  make tags for the dump file
!+
!----------------------------------------------------------------
subroutine make_output_labels(istart,iend,prestring,poststring)
 integer,          intent(in) :: istart,iend
 character(len=*), intent(in) :: prestring,poststring
 integer :: i
 character(len=20) :: istring

 do i = istart,iend
    write(istring,'(i10)') i-(istart-1)
    write(istring,'(a)') prestring//trim(adjustl(istring))//poststring
    label(i) = istring
 enddo

 return
end subroutine make_output_labels

end module analysis
