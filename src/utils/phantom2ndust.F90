!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
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
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, infile_utils, io, options, part, physcon, units
!+
!--------------------------------------------------------------------------
module analysis
 use dim,  only:ndusttypes
 use dust, only:grainsizecgs
 implicit none
 character(len=20), parameter, public :: analysistype = 'dustydisc'
 public :: do_analysis

 integer, parameter :: nr = 10
 integer, parameter :: nxn = 2*ndusttypes
 integer, parameter :: numlabels = 21
 integer, parameter :: numarrays = 4
 integer, parameter :: maxlabels = numlabels + numarrays*(ndusttypes-1)
 integer, parameter, public :: &
          iradius    = 1, &
          irhog      = 2, &
          isigma     = 3, &
          isigmadust = 4, &
          ih_H       = 5, &
          ilx        = 6, &
          ily        = 7, &
          ilz        = 8, &
          itilt      = 9, &
          itwist     = 10, &
          ipsi       = 11, &
          iH_R_init  = 12, &
          iH_R       = 13, &
          iomega     = 14, &
          ivK        = 15, &
          ics        = 16, &
          ivrgas     = 17, &
          ! initial index for arrays
          ivrdust    = 18, &
          iSt        = 19 +   (ndusttypes-1), &
          itstop     = 20 + 2*(ndusttypes-1), &
          ivrerr     = 21 + 3*(ndusttypes-1), &
          ! ending index for arrays
          ivrdustend = iSt-1, &
          iStend     = itstop-1, &
          itstopend  = ivrerr-1, &
          ivrerrend  = maxlabels

 real :: St(ndusttypes,nr),St_smooth(ndusttypes,nr)
 real :: twist(nr),twistprev(nr)
 real :: Amat(nxn,nxn)
 character(len=80) :: label(maxlabels),labelfmt,datafmt

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,deltavsum,deltav, &
                       particlemass,npart,time,iunit)
 use dim,          only:use_dustfrac,maxp
 use io,           only:fatal
 use physcon,      only:pi,jupiterm,years,au
 use part,         only:iphase,npartoftype,igas,idust,massoftype,labeltype,dustfrac, &
                        maxphase,iamtype,xyzmh_ptmass,nptmass
 use options,      only:iexternalforce
 use units,        only:umass,udist!,utime
 use dust,         only:graindens,grainsize
 use leastsquares, only:fit_slope
 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: deltavsum(:,:),deltav(:,:,:)
 real,             intent(in) :: particlemass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=14) :: output
 character(len=20) :: filename
 character(len=8)  :: basename = 'dustanal'
 character(len=20) :: discprefix
 character(len=20) :: vrdust_string,St_string,tstop_string,vrerr_string
 integer :: i,j,k,ir,ii,find_ii(1),ierr,iline,ninbin(nr),ninbindust(nr),iwarp,nptmassinit
 real :: St_from_tstop(ndusttypes,nr),St_from_sigma(ndusttypes,nr)
 real :: err,errslope,erryint,slope
 real :: R_in,R_out,R_warp,H_R,p_index,q_index,M_star,M_disc,R_c,R_cdust,Sig0,Sig0dust
 real :: R_in_dust,R_out_dust,R_warp_dust,H_R_dust,p_index_dust,M_star_dust,M_disc_dust
 real :: G,rmin,rmax,dreven,dr(nr),cs0,rho0,angx,angy,angz,ri,ri_mid,area
 real :: angtot,Ltot,tilt,dtwist
 real :: Li(3),pgasmass
 real :: rad(nr),Lx(nr),Ly(nr),Lz(nr),h_smooth(nr),sigma(nr),Sigma_smooth(nr)
 real :: cs(nr),rhog(nr),H(nr),omega(nr),vK(nr),eta(nr),tstop(ndusttypes,nr)
 real :: zsettlgas(npartoftype(igas),nr),hgas(nr),meanzgas(nr)
 real :: vrgasbin(npartoftype(igas),nr),vrdustbin(ndusttypes,npartoftype(igas),nr)
 real :: dustfracsum(ndusttypes,nr),dust_fraction(ndusttypes,nr)
 real :: unitlx(nr),unitly(nr),unitlz(nr),tp(nr)
 real :: sigmadust(nr),zsettldust(npartoftype(idust),nr),hdust(nr),meanzdust(nr)
 real :: meanvrgas(nr),meanvrdust(ndusttypes,nr)
 real :: vgas(3),vdust(3,ndusttypes),vrgas(npart),vrdust(ndusttypes,npart)
 real :: stan_dev(ndusttypes,nr)
 real :: vrsol(nxn,nr),tol = 1.e-5
 real :: psi_x,psi_y,psi_z,psi,Mdust,Mgas,Mtot,Macc,pmassi
 real :: dustfraci(ndusttypes),dustfracisum,rhoeff
 !real :: beta,phi,Q,tstop1,tstop2
 real :: d2g_ratio
 real :: log_dr,log_grid(nr+1),grid(nr+1)
 real :: cut_fact = 1000.
 real, save :: Mtot_in,Mgas_in,Mdust_in
 integer :: itype,lu,icut(nr)
 logical, save :: init = .false.
 logical :: use_log_r  = .true.
 logical :: scale_vel  = .true.
 logical :: fit_sigma  = .true.
 logical :: fixslope   = .true.
 logical :: logdata    = .true.
 integer, parameter :: isetupparams = 2
 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24
 integer, parameter :: isplash = 33
 integer, parameter :: isol    = 34
 logical :: do_precession,ifile,isetupfile

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a9,i5.5)") basename//'_',numfile
 write(*,'("Output file name is ",A)') output

 if (use_dustfrac) write(*,'("one-fluid model")')

! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G==1")')
 G = 1.0

 iline = index(dumpfile,'_')
 discprefix = dumpfile(1:iline-1)
 inquire(file=trim(discprefix)//'.discparams', exist=ifile)
 if (ifile) then
    call read_discparams(trim(discprefix)//'.discparams', &
         R_in,R_out,R_warp,H_R,p_index,R_c,q_index,M_star,M_disc,Sig0,cs0,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read .discparams file')
 else
    call read_discparams('discparams.list', &
         R_in,R_out,R_warp,H_R,p_index,R_c,q_index,M_star,M_disc,Sig0,cs0,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read discparams.list')
 endif
 
 iline = index(dumpfile,'_')
 discprefix = dumpfile(1:iline-1)
 inquire(file=trim(discprefix)//'.setup', exist=isetupfile)
 if (isetupfile) then
    call read_setup(trim(discprefix)//'.setup',d2g_ratio,grainsizecgs(:),isetupparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read .setup file')
 endif
 
 ! Print out the parameters of gas disc
 write(*,*)
 write(*,'("Gas disc parameters are:")')
 write(*,*) 'R_in    = ',R_in
 write(*,*) 'R_out   = ',R_out
 write(*,*) 'H_R     = ',H_R
 write(*,*) 'Sig0    = ',Sig0
 write(*,*) 'cs0     = ',cs0
 if(R_warp/=0.) write(*,*) 'Rwarp     = ',R_warp
 write(*,*) 'p_index = ',p_index
 if(R_c/=0.) write(*,*) 'R_c     = ',R_c
 write(*,*) 'q_index = ',q_index
 write(*,*) 'M_star  = ',M_star
 write(*,*) 'M_disc  = ',M_disc
 write(*,*)
 write(*,*)

 inquire(file=trim(discprefix)//'-'//trim(labeltype(idust))//'.discparams', exist=ifile)
 if (ifile) then
    call read_discparams(trim(discprefix)//'-'//trim(labeltype(idust))//'.discparams',&
    R_in_dust,R_out_dust,R_warp_dust,H_R_dust,p_index_dust,R_cdust,q_index,M_star_dust,M_disc_dust,Sig0dust,cs0,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read -'//trim(labeltype(idust))//' .discparams file')

    ! Print out the parameters of dust disc
    write(*,*)
    write(*,'("Dust disc parameters are:")')
    write(*,*) 'R_in    = ',R_in_dust
    write(*,*) 'R_out   = ',R_out_dust
    write(*,*) 'H_R     = ',H_R_dust
    if(R_warp_dust/=0.) write(*,*) 'Rwarp     = ',R_warp_dust
    write(*,*) 'p_index = ',p_index_dust
    if(R_cdust/=0.) write(*,*) 'R_c     = ',R_cdust
    write(*,*) 'M_disc  = ',M_disc_dust
    write(*,*)
    write(*,*)
 endif

! Setup rmin and rmax for the analysis
 if (npartoftype(idust) > 0) then
    rmin = min(R_in,R_in_dust)
    rmax = max(R_out,R_out_dust)
 else
    rmin = R_in
    rmax = R_out
 endif

 if (use_log_r) then
    !--Create a uniform grid with N+1 points between smax and smin (inclusive)
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
       dr(i) = grid(i+1)-grid(i)
    enddo
 else
    ! Set up the radius array
    dreven = (rmax-rmin)/real(nr-1)
    do i=1,nr
       rad(i) = rmin + real(i-1)*dreven
    enddo
 endif

! Initialise arrays to zero
 ninbin(:)=0
 lx(:)=0.0
 ly(:)=0.0
 lz(:)=0.0
 h_smooth(:)=0.0
 sigma(:)=0.0
 Sigma_smooth(:)=0.0
 dustfracsum(:,:)=0.0
 dust_fraction(:,:)=0.0
 sigmadust(:)=0.0
 ninbindust(:)=0
 hgas(:)=0.0
 hdust(:)=0.0
 rhog(:)=0.0
 cs(:)=0.0
 omega(:)=0.0
 vK(:)=0.0
 H(:)=0.0
 tstop(:,:)=0.0
 meanzgas(:)=0.0
 meanzdust(:)=0.0
 zsettlgas(:,:)=0.0
 zsettldust(:,:)=0.0
 vrgasbin(:,:)=0.0
 vrdustbin(:,:,:)=0.0
 meanvrgas(:)=0.0
 meanvrdust(:,:)=0.0
 icut(:) = 0

 angx = 0.0
 angy = 0.0
 angz = 0.0

 if (R_warp/=0.)then
    iwarp=3
 else
    iwarp=2
 endif

 !if the central star is represented by a sink (change to 2 if you set a binary)
 nptmassinit = 1

 !if the central star is represented by by external force
 if(iexternalforce/=0) nptmassinit = 0

 if(nptmass>nptmassinit)then
    do i=nptmassinit+1,nptmass
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
    if (maxphase==maxp) then
       itype = iamtype(iphase(i))
       pmassi = massoftype(itype)
    endif
    Mtot = Mtot + pmassi
    if (xyzh(4,i)  >  tiny(xyzh)) then
       if (use_dustfrac) dustfraci(:) = dustfrac(:,i)
       dustfracisum = sum(dustfraci)
       pgasmass = pmassi*(1. - dustfracisum)
       Mgas  = Mgas  + pgasmass
       Mdust = Mdust + pmassi*dustfracisum
    else
       pgasmass = 0.
       Macc  = Macc + pmassi
    endif

    !--Calculate vr for the gas and the N dust phases
    ri_mid = sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
    do j = 1,ndusttypes
       if (ndusttypes > 1) then
          do k = 1,ndusttypes
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
       vrgas(i)    = 0.5*(2*xyzh(1,i)*vgas(1)    + 2*xyzh(2,i)*vgas(2))   /ri_mid
       vrdust(j,i) = 0.5*(2*xyzh(1,i)*vdust(1,j) + 2*xyzh(2,i)*vdust(2,j))/ri_mid
    enddo

    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       ri = sqrt(dot_product(xyzh(1:iwarp,i),xyzh(1:iwarp,i)))
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
          area = (pi*((rad(ii)+dreven/2.)**2-(rad(ii)- dreven/2.)**2))
       endif

       if(iphase(i)==igas)then

          sigma(ii) = sigma(ii) + pgasmass/area

          Li(1) = pgasmass*(xyzh(2,i)*vxyz(3,i)-xyzh(3,i)*vxyz(2,i))
          Li(2) = pgasmass*(xyzh(3,i)*vxyz(1,i)-xyzh(1,i)*vxyz(3,i))
          Li(3) = pgasmass*(xyzh(1,i)*vxyz(2,i)-xyzh(2,i)*vxyz(1,i))

          Lx(ii)=Lx(ii)+Li(1)
          Ly(ii)=Ly(ii)+Li(2)
          Lz(ii)=Lz(ii)+Li(3)

          h_smooth(ii) = h_smooth(ii) + xyzh(4,i)

          ninbin(ii) = ninbin(ii) + 1
          zsettlgas(ninbin(ii),ii)=xyzh(3,i)
          !if (abs(xyzh(3,i)) < cut_fact*H_R*ri**(1.5-q_index)) then
          if (abs(xyzh(3,i)) < min(cut_fact,H_R*ri**(1.5-q_index))) then
             vrgasbin(ninbin(ii),ii)=vrgas(i)
             if (use_dustfrac) then
                do j = 1,ndusttypes
                   vrdustbin(j,ninbin(ii),ii) = vrdust(j,i)
                enddo
             endif
          else
             icut(ii) = icut(ii) + 1
          endif
          if (use_dustfrac) dustfracsum(:,ii) = dustfracsum(:,ii) + dustfrac(:,i)
       elseif(iphase(i)==idust) then
          sigmadust(ii) = sigmadust(ii) + massoftype(iphase(i))/area

          ninbindust(ii) = ninbindust(ii) + 1
          zsettldust(ninbindust(ii),ii)=xyzh(3,i)
       endif
    elseif (xyzh(4,i) < -tiny(xyzh) .and. iphase(i)==igas) then !ACCRETED
       angx = angx + pgasmass*(xyzh(2,i)*vxyz(3,i) - xyzh(3,i)*vxyz(2,i))
       angy = angy + pgasmass*(xyzh(3,i)*vxyz(1,i) - xyzh(1,i)*vxyz(3,i))
       angz = angz + pgasmass*(xyzh(1,i)*vxyz(2,i) - xyzh(2,i)*vxyz(1,i))

    endif
 enddo
 write(*,*)"Massa della polvere: ",Mdust
 write(*,*)"Massa del gas: ",Mgas

 if (.not.init) then
    open(newunit=lu,file='dustmass.ev',status='replace')
    write(lu,"('# ',5('[',i2.2,1x,a12,']',1x))") 1,'time',2,'Mtot',3,'Mgas',4,'Mdust',5,'Macc'
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
 write(lu,*) time,Mtot,Mgas,Mdust,Macc
 close(lu)

 !--Usually Sig0 does not match the binned sigma calculated above so we do a least squares
 !  fit in order to ensure that the assumptions for the analytic solution matches the data
 if (fit_sigma) then
    print*,' '
    print*,'Fitting the power-law surface density profile for the analytic solution...'
    slope = -p_index
    call fit_slope(nr,rad,sigma,slope,Sig0,err,errslope,erryint,rmin,rmax,logdata,fixslope)
    if (logdata) Sig0 = 10.**Sig0
    rho0 = Sig0/(sqrt(2.*pi)*H_R)
    print*,'   Corrected value for Sig0 = ',Sig0
    print*,'   Corrected value for rho0 = ',rho0
    print*,'...Done'
    print*,' '
 endif

 ! Calculate other disc quantities
 do i=1,nr
    cs(i)    = cs0*rad(i)**(-q_index)
    rhog(i)  = rho0*rad(i)**(-(p_index-q_index+1.5))
    omega(i) = sqrt(G*M_star/rad(i)**3)
    vK(i)    = sqrt(G*M_star/rad(i))
    H(i)     = cs(i)/omega(i)
    Sigma_smooth(i) = Sig0*rad(i)**(-p_index)
 enddo

! and thus the Stokes Number
! St = tstop*omegaK = rhoeff*s/(cs*rhog)*omegaK
!                   = rhoeff*s/(sqrt(2*pi)*sigma*exp(-0.5*(z/Hg)**2)))
!                   = rhoeff*s/(sqrt(2*pi)*sigma)    ! if z = 0
 do i=1,nr
       if (ndusttypes == 1) then
          tstop(:,i) = rhoeff*grainsize(:)/(cs(i)*rhog(i))

          St_from_tstop(:,i) = tstop(:,i)*omega(i)
          St_from_sigma(:,i) = rhoeff*grainsize(:)*(sqrt(2.*pi)/Sigma_smooth(i))
          if (all(abs((St_from_tstop(:,i)-St_from_sigma(:,i))/St_from_sigma(:,i))>tol)) then
             print*,' '
             print*,'WARNING!!! There is a problem in the scaling somewhere!'
             print*,'   St calculated from tstop is not the same as when using Sigma'
             print*,' '
          endif
          St_smooth(:,i) = St_from_sigma(:,i)

          St(:,i) = rhoeff*grainsize(:)*(sqrt(2.*pi)/sigma(i))
       elseif (ndusttypes == 2) then
          St_smooth(:,i) = rhoeff*grainsize(:)*(sqrt(2.*pi)/Sigma_smooth(i))
          St(:,i) = rhoeff*grainsize(:)*(sqrt(2.*pi)/sigma(i))

          tstop(:,i) = St_smooth(:,i)/omega(i)

          !!--Actual *physical* stopping time for 2 coupled dust phases
          !!  as opposed to the stopping time of a single grain neglecting coupling
          !!  These equations are not relevant for the Bai and Stone solution though
          !!-----------------------------------
          !! Find tstop and St using Sigma_smooth
          !tstop(:,i) = rhoeff*grainsize(:)*(sqrt(2.*pi)/Sigma_smooth(i))/omega(i)

          !beta = tstop(1,i)/tstop(2,i)
          !phi = dustfraci(1)/dustfracisum
          !
          !Q = (4.*beta*(1.-dustfracisum))/((1.-phi)*(1.-dustfracisum*(1.-phi)) + &
          !    beta*phi*(1.-dustfracisum*phi))**2
          !
          !tstop1 = 1./(0.5*(1./tstop(1,i)+1./tstop(2,i))*(1.+sqrt(1.-Q)))
          !tstop2 = 1./(0.5*(1./tstop(1,i)+1./tstop(2,i))*(1.-sqrt(1.-Q)))

          !St_smooth(1,i) = tstop1*omega(i)
          !St_smooth(2,i) = tstop2*omega(i)

          !!-----------------------------------
          !! Find tstop and St using sigma
          !tstop(:,i) = rhoeff*grainsize(:)*(sqrt(2.*pi)/sigma(i))/omega(i)

          !beta = tstop(1,i)/tstop(2,i)
          !phi = dustfraci(1)/dustfracisum
          !
          !Q = (4.*beta*(1.-dustfracisum))/((1.-phi)*(1.-dustfracisum*(1.-phi)) + &
          !    beta*phi*(1.-dustfracisum*phi))**2
          !
          !tstop1 = 1./(0.5*(1./tstop(1,i)+1./tstop(2,i))*(1.+sqrt(1.-Q)))
          !tstop2 = 1./(0.5*(1./tstop(1,i)+1./tstop(2,i))*(1.-sqrt(1.-Q)))

          !St(1,i) = tstop1*omega(i)
          !St(2,i) = tstop2*omega(i)
          !!-----------------------------------
       else
          stop 'ndusttypes > 2 have not been coded...'
       endif

       eta(i)  = ((1.5+p_index+q_index))*(cs(i)/vK(i))**2
 enddo

 ! Computing Hgas, Hdust, vrgas, vrdust
 do i=1,nr
    if(ninbin(i)>1)then
       print*,'The # of particles in bin',i,'is',ninbin(i)-icut(i)
       meanzgas(i) =sum(zsettlgas(1:ninbin(i),i))/real(ninbin(i))
       meanvrgas(i) = sum(vrgasbin(1:ninbin(i),i))/(real(ninbin(i))-icut(i))
       if (use_dustfrac) then
          do j = 1,ndusttypes
             meanvrdust(j,i) = sum(vrdustbin(j,1:ninbin(i),i))/(real(ninbin(i))-icut(i))
             stan_dev(j,i) = sqrt(sum(((vrdustbin(j,1:ninbin(i),i)-meanvrdust(j,i)) &
                             /(eta(i)*vK(i)))**2)/(real(ninbin(i)) - icut(i) - 1))
             if (scale_vel) then
                meanvrdust(j,i) = meanvrdust(j,i)/(eta(i)*vK(i))
             endif
          enddo
       endif
       hgas(i)=sqrt(sum(((zsettlgas(1:ninbin(i),i)-meanzgas(i))**2)/(real(ninbin(i)-1))))
       if (use_dustfrac) dust_fraction(:,i)=dustfracsum(:,i)/real(ninbin(i))
    endif
    if(ninbindust(i)>1)then
       meanzdust(i)=sum(zsettldust(1:ninbindust(i),i))/real(ninbindust(i))
       hdust(i)=sqrt(sum(((zsettldust(1:ninbindust(i),i)-meanzdust(i))**2)/(real(ninbindust(i)-1))))
    endif
 enddo
 print*,' '
 print*,'Warning: removed a total of',sum(icut),'particles from the velocity'
 print*,' '

! Plotting Hdust/Hgas
! print*,' Hdust = ',hdust,'Hgas = ',hgas, 'Hdust/Hgas = ',hdust/hgas

! Print angular momentum of accreted particles
 angtot = sqrt(angx*angx + angy*angy + angz*angz)
 print*,' angular momentum of accreted particles = ',angtot


! Convert total angular momentum into a unit vector, and average h_smooth
 do i = 1,nr
    Ltot = sqrt(Lx(i)*Lx(i) + Ly(i)*Ly(i) + Lz(i)*Lz(i))

    unitlx(i) = Lx(i)/Ltot
    unitly(i) = Ly(i)/Ltot
    unitlz(i) = Lz(i)/Ltot

    if (ninbin(i) > 0) h_smooth(i) = h_smooth(i)/ninbin(i)
 enddo

! Now loop over rings to calculate required quantities
 do i = 1, nr
    if(ninbin(i)==0) then
       lx(i)=0.0
       ly(i)=0.0
       lz(i)=0.0
       sigma(i)=0.0
       h_smooth(i) = 0.0
    else
       h_smooth(i) = h_smooth(i)/H(i)
    endif
 enddo

 ! Make labels for header of output file
 label(iradius)    = 'radius'
 label(irhog)      = 'rhog'
 label(isigma)     = 'sigma'
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
 label(iomega)     = 'OmegaK'
 label(ivK)        = 'vK'
 label(ics)        = 'cs'
 label(ivrgas)     = '<vrgas>'

 ! Make N vrdust labels
 do i = ivrdust,ivrdustend
    write(vrdust_string,'(I10)') i-(ivrdust-1)
    write(vrdust_string,'(A)') '<vrdust'//trim(adjustl(vrdust_string))//'>'
    label(i) = vrdust_string
 enddo

 ! Make N Stokes labels
 do i = iSt,iStend
    write(St_string,'(I10)') i-(iSt-1)
    write(St_string,'(A)') 'St'//trim(adjustl(St_string))
    label(i) = St_string
 enddo

 ! Make N tstop labels
 do i = itstop,itstopend
    write(tstop_string,'(I10)') i-(itstop-1)
    write(tstop_string,'(A)') 'tstop'//trim(adjustl(tstop_string))
    label(i) = tstop_string
 enddo

 ! Make N vr error labels
 do i = ivrerr,ivrerrend
    write(vrerr_string,'(I10)') i-(ivrerr-1)
    write(vrerr_string,'(A)') 'vrerr'//trim(adjustl(vrerr_string))
    label(i) = vrerr_string
 enddo

 ! Make splash.columns file using the labels above
 write(filename,"(a16)") basename//'.columns'
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
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 if (npartoftype(idust)==0)then
    if (use_dustfrac) then
       write(iunit,labelfmt) (j,trim(adjustl(label(j))),j=1,maxlabels)
    else
       stop '...still needs work...'
    endif
 else
    stop '...still needs work...'
 endif

 do_precession = .false.

 ! Write the data to the file
do i=1,nr
    if(i /= 1.and.i /= nr) then
       psi_x=(unitlx(i+1)-unitlx(i-1))/(rad(i+1)-rad(i-1))
       psi_y=(unitly(i+1)-unitly(i-1))/(rad(i+1)-rad(i-1))
       psi_z=(unitlz(i+1)-unitlz(i-1))/(rad(i+1)-rad(i-1))
       psi=sqrt(psi_x**2 + psi_y**2 + psi_z**2)*rad(i)
    else
       psi=0.
    endif

    if (ninbin(i) > 0) then
       tilt  = acos(unitlz(i))
       twist(i) = atan2(unitly(i),unitlx(i))
       if (i==1 .or. time==0.0) then
          twistprev(i) = 0.0
       endif
       ! Taking into account negative twist
       if (twist(i) < 0) then
          twistprev(i) = 2.*pi + twist(i)
       else
          twistprev(i) = twist(i) !cumulative twist
       endif
    else
       tilt = 0.0
       twist = 0.0
       dtwist = 0.0
    endif

! Calculate the precession time
    if (twist(i) > tiny(twist(i))) then
       tp(i) = time*2.*pi/twist(i)
    else
       tp(i) = 0.0
    endif

    if (npartoftype(idust)==0)then
       if (ninbin(i) > 0) then
          if (use_dustfrac) then
             write(iunit,datafmt) &
                   rad(i), &
                   rhog(i), &
                   sigma(i), &
                   sigma(i)*sum(dust_fraction(:,i))/(1.-sum(dust_fraction(:,i))), &
                   !sigma(i)*(1.-sum(dust_fraction(:,i))), &
                   !sigma(i)*sum(dust_fraction(:,i)), &
                   h_smooth(i), &
                   unitlx(i), &
                   unitly(i), &
                   unitlz(i), &
                   tilt, &
                   twist(i),  &
                   psi, &
                   H(i)/rad(i), &
                   hgas(i)/rad(i), &
                   omega(i), &
                   vK(i), &
                   cs(i), &
                   meanvrgas(i),  &
                   (meanvrdust(j,i),j=1,ndusttypes), &
                   (St(j,i),j=1,ndusttypes), &
                   (tstop(j,i),j=1,ndusttypes), &
                   (stan_dev(j,i),j=1,ndusttypes)
          else
             stop '...still needs work...'
          endif
       endif
    else
       stop '...still needs work...'
    endif

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

 print *,' '
 print *,'Solving the Bai and Stone (2010) analytic solution for radial drift...'
 call solve_bai_stone_2010(d2g_ratio,eta,vK,vrsol,scale_vel)
 print*,'...Done'
 print *,' '

 ! Print solution to file
 write(filename,"(a12)") basename//'.sol'
 open(unit=isol,file=filename,status="replace")
 write(labelfmt,'(I10)') iStend-ivrdust+3
 write(labelfmt,'(A)') '(''#'','//trim(adjustl(labelfmt))//'(1x,''['',i2.2,1x,a11,'']'',2x))'
 write(datafmt,'(I10)') iStend-ivrdust+3
 write(datafmt,'(A)') '('//trim(adjustl(datafmt))//'(es18.10,1X))'
 write(isol,labelfmt) 1,trim(adjustl(label(iradius))), &
                      2,trim(adjustl(label(isigma))), &
                      (j-ivrdust+3,trim(adjustl(label(j))),j=ivrdust,iStend)
 do ir = 1,nr
    write(isol,datafmt) rad(ir),Sigma_smooth(ir),(vrsol(j,ir),j=1,ndusttypes),(St_smooth(j,ir),j=1,ndusttypes)
 enddo
 close(unit=isol)

 ! Make splash.columns file using the labels above
 write(filename,"(a11)") 'sol.columns'
 open(unit=isplash,file=filename,status="replace")
 write(isplash,'(A)') trim(adjustl(label(iradius)))
 write(isplash,'(A)') trim(adjustl(label(isigma)))
 do i = ivrdust,ivrdustend
    write(isplash,'(A)') trim(adjustl(label(i)))
 enddo
 do i = iSt,iStend
    write(isplash,'(A)') trim(adjustl(label(i)))
 enddo
 close(unit=isplash)

 return
end subroutine do_analysis

subroutine solve_bai_stone_2010(d2g_ratio,eta,vK,vrsol,scale_vel)
 use dim, only:ndusttypes
 use dust, only:set_dustfrac,smincgs,smaxcgs,sindex
 integer, parameter :: dp = selected_real_kind(14, 60)
 real (dp) :: soln(nxn)
 real, intent(in)  :: d2g_ratio,eta(:),vK(:)
 real, intent(out) :: vrsol(:,:)
 logical, intent(in) :: scale_vel
 real :: dustfraci(ndusttypes)
 real :: Imat(ndusttypes,ndusttypes)
 real :: Lambda(ndusttypes,ndusttypes)
 real :: Gamma(ndusttypes,ndusttypes)
 real :: Bmat(nxn)
 real :: dustfracsum
 integer :: i,ir,ierr

 vrsol(:,:) = 0.
 do ir = 1,nr
    Bmat(1:ndusttypes) = 0.
    Bmat(ndusttypes+1:nxn) = 1.

    call set_dustfrac(d2g_ratio,dustfraci,smincgs,smaxcgs,sindex)

    dustfracsum = sum(dustfraci)

    Imat(:,:) = 0.
    do i = 1,ndusttypes
       Imat(i,i) = 1.
    enddo

    Lambda(:,:) = 0.
    do i = 1,ndusttypes
       Lambda(i,i) = St_smooth(i,ir)
       Gamma(i,:)  = dustfraci(:)/(1.-dustfracsum)
    enddo

    Amat(:,:) = 0.
    do i = 1,ndusttypes
       Amat(i,           :) = [ (Imat(i,:)+Gamma(i,:)), -2.*Lambda(i,:)        ]
       Amat(i+ndusttypes,:) = [ 0.5*Lambda(i,:)       , (Imat(i,:)+Gamma(i,:)) ]
    enddo

    call dple(rowk,nxn,Bmat,soln,ierr)
    if (ierr /= 0) then
       write(*, *) ' Error = ', ierr
    end if

    ! Save the solution to an array
    vrsol(:,ir) = -eta(ir)*vK(ir)*soln(:)
    if (scale_vel) then
       vrsol(:,ir) = vrsol(:,ir)/(eta(ir)*vK(ir))
    endif
 enddo

 print*,'   The solution assumes the following dust properties:'
 print*,'   (Warning! Some of these need to be changed manually)'
 print*,'   Total dust-to-gas ratio: ',d2g_ratio
 print*,'   Total dust fraction: ',dustfracsum
 if (ndusttypes>1) then
    print*,'   smin [cm] = ',smincgs
    print*,'   smax [cm] = ',smaxcgs
    print*,'   sindex = ',sindex
 endif
 do i = 1,ndusttypes
    print*,i,': dustfrac = ',dustfraci(i),';  s [cm] = ',grainsizecgs(i)
 enddo

 return
end subroutine solve_bai_stone_2010


subroutine rowk(n, k, r)
 use dim,  only:ndusttypes
 integer, parameter  :: dp = selected_real_kind(14, 60)
 integer, intent(in)     :: n, k
 real (dp), intent(out)  :: r(:)

 r(:) = Amat(k,:)

 return
end subroutine rowk





!***********************************************************************
!* Solve AX = B using a partial pivoting algorithm and reduced storage *
!* ------------------------------------------------------------------- *
!* SAMPLE RUN:                                                         *
!*                                                                     *
!* System to solve:                                                    *
!*   2.0000  -1.0000   1.0000   7.0000 -12.5400    5.0000              *
!*   1.0000   5.0000  -2.0000  -8.0000 100.0000    1.0000              *
!*   3.0000  -2.0000   3.0000  45.0000  27.3333    3.0000              *
!*  11.0000   0.5500  -2.0000  -4.0000   1.0000    4.0000              *
!*  33.0000   2.0000  -3.0000   5.0000   7.3333  -10.0000              *
!*                                                                     *
!* Solution is:                                                        *
!*   2.11149597961869                                                  *
!*  -25.8290267820056                                                  *
!*   8.17423194407132                                                  *
!*  -2.52146730210577                                                  *
!*   1.24210363401706                                                  *
!*                                                                     *
!* ------------------------------------------------------------------- *
!* Ref.: "Wassyng, A. - Solving Ax = b: A method with reduced storage  *
!*        requirements, SIAM J. Numerical Analysis, vol.19 (1982),     *
!*        pp. 197-204".                                                *
!*                                                                     *
!*                                  F90 Release By J-P Moreau, Paris.  *
!*                                         (www.jpmoreau.fr)           *
!***********************************************************************
subroutine dple(rowk, n, b, c, ierr)

! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-16  Time: 12:26:32

! ******************************************************************
!        SOLUTION OF LINEAR EQUATIONS WITH REDUCED STORAGE
! ******************************************************************

! Uses the Henderson-Wassyng partial pivot algorithm.
! Wassyng, A. 'Solving Ax = b: A method with reduced storage requirements',
! SIAM J. Numerical Analysis, vol.19 (1982), pp. 197-204.

! The user must provide a routine ROWK to return the requested row of the
! matrix A.

! N.B. Arguments D and IP have been removed.

implicit none
integer, parameter  :: dp = selected_real_kind(14, 60)

integer, intent(in)     :: n
real (dp), intent(in)   :: b(n)
real (dp), intent(out)  :: c(n)
integer, intent(out)    :: ierr

! external rowk
interface
  subroutine rowk(n, k, r)
    implicit none
    integer, parameter  :: dp = selected_real_kind(14, 60)
    integer, intent(in)     :: n, k
    real (dp), intent(out)  :: r(:)
  end subroutine rowk
end interface

! Local variables
real (dp)  :: bk, cj, ck, c1, dkj
real (dp), parameter  :: zero = 0.0_dp
real (dp)  :: wk(n*n/4 + n + 3)
integer    :: i, iflag, ij, ijold, ik, iwk(n), j, k, kjold, km1, kp1,   &
              last, lastm1, lcol, lcolp1, m, maxwk, mjold, nm1, np1

! Set the necessary constants
ierr = 0
maxwk = n * n / 4 + n + 3
np1 = n + 1
k = 1
iflag = -1

! Get the first column of the transposed system
call rowk(n, 1, c)
bk = b(1)

if (n <= 1) then
  if (c(1) == zero) GO TO 130
  c(1) = bk / c(1)
  return
end if

! Find the pivot for column 1
m = 1
do  i = 2, n
  if (abs(c(m)) < abs(c(i))) m = i
end do

iwk(1) = m
c1 = c(m)
c(m) = c(1)
c(1) = c1
if (c(1) /= zero) then

! Find the first elementary matrix and store it in d
  do  i = 2, n
    wk(i-1) = -c(i) / c(1)
  end do
  wk(n) = bk / c(1)

! k loop - each k for a new column of the transposed system
  do  k = 2, n
    kp1 = k + 1
    km1 = k - 1

! Get column k
    call rowk(n, k, c)
    do  j = 1, km1
      m = iwk(j)
      cj = c(j)
      c(j) = c(m)
      c(m) = cj
    end do
    bk = b(k)

    iflag = -iflag
    lcol = np1 - k
    lcolp1 = lcol + 1
    lastm1 = 1
    last = maxwk - n + k
    if (k /= 2) then

      lastm1 = maxwk - n + km1
      if (iflag < 0) last = last - n + k - 2
      if (iflag > 0) lastm1 = lastm1 - n + k - 3
    end if

! j loop - effect of columns 1 to k-1 of l-inverse
    do  j = 1, km1
      cj = c(j)
      ij = (j-1) * lcolp1
      if (j == km1) ij = lastm1 - 1

! i loop - effect of l-inverse on rows k to n+1
      do  i = k, n
        ij = ij + 1
        c(i) = c(i) + wk(ij) * cj
      end do
      bk = bk - wk(ij+1) * cj
    end do

! k=n case
    m = k
    if (k >= n) then
      if (c(k) == zero) GO TO 130
      wk(last) = bk / c(k)
    else

! Find the pivot
      do  i = kp1, n
        if (abs(c(m)) < abs(c(i))) m = i
      end do

      iwk(k) = m
      ck = c(m)
      c(m) = c(k)
      c(k) = ck
      if (c(k) == zero) GO TO 130

! Find the k-th elementary matrix
      ik = last
      do  i = kp1, n
        wk(ik) = -c(i) / c(k)
        ik = ik + 1
      end do
      wk(ik) = bk / c(k)
    end if

! Form the product of the elementary matrices
    do  j = 1, km1
      kjold = j * lcolp1 + k - np1
      mjold = kjold + m - k
      ij = (j-1) * lcol
      ijold = ij + j
      if (j == km1) then

        kjold = lastm1
        mjold = lastm1 + m - k
        ijold = lastm1
      end if

      ik = last - 1
      dkj = wk(mjold)
      wk(mjold) = wk(kjold)
      do  i = kp1, np1
        ij = ij + 1
        ijold = ijold + 1
        ik = ik + 1
        wk(ij) = wk(ijold) + wk(ik) * dkj
      end do
    end do
  end do

  last = maxwk
  if (iflag < 0) last = maxwk - 2
  wk(n) = wk(last)

! Insert the solution in c

  c(1:n) = wk(1:n)

  nm1 = n - 1
  do  i = 1, nm1
    k = n - i
    m = iwk(k)
    ck = c(k)
    c(k) = c(m)
    c(m) = ck
  end do
  return
end if

! The system is singular
130 ierr = k

return
end subroutine dple

!----------------------------------------------------------------
!+
!  Read disc information from discparams.list file
!+
!----------------------------------------------------------------
subroutine read_discparams(filename,R_in,R_out,R_warp,H_R,p_index,R_c,q_index, &
                           M_star,M_disc,Sig0,cs0,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: R_in,R_out,R_warp,H_R,p_index,q_index
 real,             intent(out) :: M_star,M_disc,R_c,Sig0,cs0
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)

! Read in parameters from the file discparams.list
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(R_in,'R_in',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_out,'R_out',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_warp,'R_warp',db,ierr)
 call read_inopt(H_R,'H_R',db,ierr)
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
 call read_inopt(Sig0,'Sig0',db,ierr)
 if (ierr /= 0) return
 call read_inopt(cs0,'cs0',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_discparams

!----------------------------------------------------------------
!+
!  Read disc information from discparams.list file
!+
!----------------------------------------------------------------
subroutine read_setup(filename,d2g_ratio,grainsize,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: d2g_ratio,grainsize(:)
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 real :: sgrain
 type(inopts), allocatable :: db(:)
 
! Read in parameters from the .setup file
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(d2g_ratio,'dust_to_gas_ratio',db,ierr)
 if (ierr /= 0) return
 call read_inopt(sgrain,'grainsizeinp',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

 grainsize(:) = sgrain

end subroutine read_setup

end module analysis

!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom2ndust
!
!  DESCRIPTION: This program is a post-processing tool to calculate dust
!               properties for the ndusttypes
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantom2ndust dumpfile(s)
!
!  DEPENDENCIES: deriv, dim, initial, io, kernel, part, readwrite_dumps
!+
!--------------------------------------------------------------------------
program phantom2ndust
 use dim,             only:maxp,tagline,ndusttypes
 use part,            only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol, &
                           hfact,rhoh,dhdrho,igas,isetphase,iphase,maxphase,&
                           dustfrac,ddustfrac,deltav,massoftype
 use io,              only:set_io_unit_numbers,iprint,idisk1,ievfile
 use initial,         only:initialise
 use readwrite_dumps, only:read_dump,write_fulldump
 use deriv,           only:derivs
 use kernel,          only:hfact_default
 use dust,            only:init_drag
 use eos,             only:extract_eos_from_hdr
 use fileutils,       only:numfromfile
 use analysis,        only:do_analysis
 implicit none
 integer :: nargs
 integer :: ierr,iarg
 character(len=120) :: dumpfile
 real :: time,dtdum
 real, allocatable :: dustfracisum1(:),deltavsum(:,:)

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a)",trim(tagline)
    print "(a)",' Usage: phantom2ndust dumpfile(s)'
    stop
 endif

 print "(/,a,/)",' Phantom2ndust: dust thou art, and unto dust shalt thou return'

 call initialise()
! if (ndivcurlv < 1) stop 'error: need ndivcurlv=1 for this to do anything'
! if (ndivcurlv < 4) print "(a)",' WARNING: need ndivcurlv=4 in dim file to get curl v as well as div v'

 over_args: do iarg=1,nargs

    call get_command_argument(iarg,dumpfile)
!
!--read particle setup from dumpfile
!
    extract_eos_from_hdr = .true.
    call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) stop 'error reading dumpfile'
    if (hfact<epsilon(hfact)) hfact=hfact_default
!
!--initialise the dust quantities
!
    call init_drag(ierr)
!
!--allocate memory for deltav calculations
!
    if (.not. allocated(deltavsum)) allocate(deltavsum(3,npart))
    if (.not. allocated(dustfracisum1)) allocate(dustfracisum1(npart))
!
!--calculate derivatives including deltav
!
    if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustfrac,ddustfrac,0.,0.,dtdum)

    dustfracisum1(:) = 1./sum(dustfrac(:,1:npart),1)
    deltavsum(1,:)   = dustfracisum1*sum(dustfrac(:,1:npart)*deltav(1,:,1:npart),1)
    deltavsum(2,:)   = dustfracisum1*sum(dustfrac(:,1:npart)*deltav(2,:,1:npart),1)
    deltavsum(3,:)   = dustfracisum1*sum(dustfrac(:,1:npart)*deltav(3,:,1:npart),1)

    call do_analysis(trim(dumpfile),numfromfile(dumpfile),xyzh,vxyzu,deltavsum,deltav, &
                     massoftype(1),npart,time,ievfile)

 enddo over_args
 print "(/,a,/)",' Phantom2ndust: To be honest, I think I would prefer the sandbox'

end program phantom2ndust
