!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module discanalysisutils
!
! Routine to calculate azimuithally averaged properties in a disc
!  Can handle gas disc, gas disc + sinks and warped discs
!
! :References:
!
! :Owner: Bec Nealon
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, externalforces, options, physcon, prompting,
!   vectorutils
!
 implicit none

 character(len=20), parameter, public :: analysistype = 'disc'
 public :: disc_analysis

 private

contains

!----------------------------------------------------------------
!+
!  Disc analysis routine - so that this can be called externally
!  Adapted from a routine written by Chris Nixon
!+
!----------------------------------------------------------------

subroutine disc_analysis(xyzh,vxyz,npart,pmass,time,nbin,rmin,rmax,G,M_star,&
                     tilt,tilt_acc,twist,twistprev,psi,H,bin,h_smooth,sigma,unitlx,unitly,unitlz,Lx,Ly,Lz,&
                     ecc,ninbin,assume_Ltot_is_same_as_zaxis,xyzmh_ptmass,vxyz_ptmass,nptmass)
 use physcon,        only:pi
 use centreofmass,   only:get_total_angular_momentum,reset_centreofmass
 use externalforces, only:iext_einsteinprec
 use options,        only:iexternalforce
 use vectorutils,    only:rotatevec
 use prompting,      only:prompt
 real, intent(inout)              :: xyzh(:,:),vxyz(:,:),pmass,time
 integer, intent(in)              :: nbin,npart
 real, intent(in)                 :: rmin,rmax,G,M_star
 logical, intent(in)              :: assume_Ltot_is_same_as_zaxis
 real, optional, intent(inout)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, optional, intent(inout) :: nptmass
 real, intent(out)                :: Lx(nbin),Ly(nbin),Lz(nbin),bin(nbin),unitly(nbin)
 real, intent(out)                :: tilt(nbin),tilt_acc(nbin),twistprev(nbin)
 real, intent(out)                :: psi(nbin),H(nbin),ecc(nbin),unitlz(nbin)
 real, intent(out)                :: sigma(nbin),h_smooth(nbin),unitlx(nbin)
 integer, intent(out)             :: ninbin(nbin)
 real                             :: dbin,angx,angy,angz,unitangz
 real                             :: angtot,Ltot
 real                             :: rsphi,rcyli,area,Ei,mu,term,ecci
 real                             :: Li(3),xi(3),vi(3),Limag,dtwist
 real                             :: psi_x,psi_y,psi_z,tp(nbin)
 real                             :: L_tot(3),L_tot_mag,temp(3),temp_mag
 real                             :: rotate_about_z,rotate_about_y
 real                             :: zdash,twist(nbin),ai
 real, allocatable                :: zsetgas(:,:),myz(:)
 integer                          :: mybin(npart)
 integer                          :: i,ii,sorting_choice,iallocerr
 logical                          :: rotate,acceptance

! Options
! Rotation options for tilt and twist calculations
 if (assume_Ltot_is_same_as_zaxis) then
    rotate = .false.
 else
    rotate = .true.
 endif

! Sorting: Particles are sorted by one of
! 1 = cylindrical radius (this is the default option)
! 2 = semi-major axis
 sorting_choice = 1

! Set up the radius array
 dbin = (rmax-rmin)/real(nbin-1)
 do i=1,nbin
    bin(i)=rmin + real(i-1)*dbin
 enddo

! Initialise everything
 ninbin(:)=0
 lx(:)=0.0
 ly(:)=0.0
 lz(:)=0.0
 h_smooth(:)=0.0
 sigma(:)=0.0
 ecc=0.0
 angx = 0.0
 angy = 0.0
 angz = 0.0
 twist = 0.0
 twistprev = 0.0
 mu = G*M_star

 allocate(zsetgas(npart,nbin),stat=iallocerr)
 ! If you don't have enough memory to allocate zsetgas, then calculate H the slow way with less memory.
 if (iallocerr/=0) then
    write(*,'(/,a)') ' WARNING: Could not allocate memory for array zsetgas!'
    write(*,'(a)')   '          (It possibly requires too much memory)'
    write(*,'(a,/)') '          Height of the disc, H, will be calculated the slow way.'
    if (allocated(zsetgas)) deallocate(zsetgas)
    allocate(myz(npart))
 endif

! Move everything so that the centre of mass is at the origin
! and run everything in CoM frame
! NB: this is not suitable for discs/flybys
 if (nptmass > 0) then
    call reset_centreofmass(npart,xyzh,vxyz,nptmass,xyzmh_ptmass,vxyz_ptmass)
 endif

! Loop over particles putting properties into the correct bin
 do i = 1,npart

    ! i for the particle number, ii for the bin number
    xi = xyzh(1:3,i)
    vi = vxyz(1:3,i)

    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE

       rsphi = sqrt(dot_product(xi(1:3),xi(1:3)))
       rcyli = sqrt(dot_product(xi(1:2),xi(1:2)))

       Li(1) = pmass*(xi(2)*vi(3) - xi(3)*vi(2))
       Li(2) = pmass*(xi(3)*vi(1) - xi(1)*vi(3))
       Li(3) = pmass*(xi(1)*vi(2) - xi(2)*vi(1))

       Limag = sqrt(dot_product(Li,Li))/pmass

       ! NB: No internal energy as isothermal
       if (iexternalforce==iext_einsteinprec) then
          Ei = 0.5*dot_product(vi,vi) - G*M_star/rsphi -3.*G*M_star/(rsphi**2)
          term = 2.*Ei*(Limag**2 - 6.*mu*mu)/(mu**2)
       else
          Ei = 0.5*dot_product(vi,vi) - G*M_star/rsphi
          term = 2.*Ei*Limag**2/(mu**2)
       endif

       ai = -M_star/(2.*Ei)

       ! Now choose and store the bin
       if (sorting_choice==2) then
          ii = int((ai-bin(1))/dbin + 1)
       else
          ii = int((rcyli - bin(1))/dbin + 1)
       endif
       mybin(i) = ii

       ! If it's not in the range for the analysis, cycle
       if (ii > nbin .or. ii < 1) cycle

       area = (pi*((bin(ii)+dbin/2.)**2-(bin(ii)- dbin/2.)**2))
       sigma(ii) = sigma(ii) + pmass/area

       ecci = sqrt(1. + term)

       Lx(ii)=Lx(ii)+Li(1)
       Ly(ii)=Ly(ii)+Li(2)
       Lz(ii)=Lz(ii)+Li(3)
       ecc(ii) = ecc(ii) + ecci
       h_smooth(ii) = h_smooth(ii) + xyzh(4,i)

       ninbin(ii) = ninbin(ii) + 1

    elseif (xyzh(4,i) < -tiny(xyzh)) then !ACCRETED
       angx = angx + pmass*(xi(2)*vi(3) - xi(3)*vi(2))
       angy = angy + pmass*(xi(3)*vi(1) - xi(1)*vi(3))
       angz = angz + pmass*(xi(1)*vi(2) - xi(2)*vi(1))
    endif
 enddo

! Convert total angular momentum into a unit vector, and average h_smooth
 do i = 1,nbin
    Ltot = sqrt(Lx(i)*Lx(i) + Ly(i)*Ly(i) + Lz(i)*Lz(i))

    if (Ltot<tiny(Ltot)) then
       unitlx(i) = 0.
       unitly(i) = 0.
       unitlz(i) = 0.
    else
       unitlx(i) = Lx(i)/Ltot
       unitly(i) = Ly(i)/Ltot
       unitlz(i) = Lz(i)/Ltot
    endif

    if (ninbin(i) > 0) then
       h_smooth(i) = h_smooth(i)/ninbin(i)
       ecc(i) = ecc(i)/ninbin(i)
    endif
 enddo

 ! Now go through and get the scale height right
 ! This has to be done in a separate loop because
 ! the unit angular momentum vector of the bin
 ! must already be known

 ninbin = 0
 do i = 1,npart
    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       xi = xyzh(1:3,i)
       ii = mybin(i)

       if (ii > nbin .or. ii < 1) cycle
       ninbin(ii) = ninbin(ii) + 1

       ! get vertical height above disc midplane == z if disc is not warped
       zdash = unitlx(ii)*xi(1) + unitly(ii)*xi(2) + unitlz(ii)*xi(3)
       if (iallocerr==0) then
          zsetgas(ninbin(ii),ii) = zdash
       else
          myz(i) = zdash
       endif
    endif
 enddo

! Calculate H from the particle positions
 if (iallocerr==0) then
    call calculate_H_fast(nbin,H,zsetgas,ninbin)
    deallocate(zsetgas) ! clean up
 else
    call calculate_H_slow(nbin,npart,H,mybin,ninbin,myz)
    deallocate(myz) ! clean up
 endif

! Print angular momentum of accreted particles
 angtot = sqrt(angx*angx + angy*angy + angz*angz)

! For unit angular momentum accreted, z component
 if (angtot<tiny(angtot)) then
    unitangz = 0.
 else
    unitangz = angz/angtot
 endif

 print*,' angular momentum of accreted particles = ',angtot!,angx,angy,angz,unitangz

! Now loop over rings to calculate required quantities
 do i = 1, nbin
    if (ninbin(i)==0 .or. ninbin(i)==1) then
       lx(i)=0.0
       ly(i)=0.0
       lz(i)=0.0
       sigma(i)=0.0
       h_smooth(i) = 0.0
       H(i) = 0.
    else
       h_smooth(i) = h_smooth(i)/H(i)
    endif
 enddo

 ! Calculate the total angular momentum vector and rotate unitl[x,y,z] if required
 if (rotate) then
    if (nptmass /= 0) then
       call get_total_angular_momentum(xyzh,vxyz,npart,L_tot,xyzmh_ptmass,vxyz_ptmass,nptmass)
    else
       call get_total_angular_momentum(xyzh,vxyz,npart,L_tot)
    endif

    temp = (/L_tot(1),L_tot(2),0./)
    temp_mag = sqrt(dot_product(temp,temp))
    rotate_about_z = -acos(dot_product((/1.,0.,0./),temp/temp_mag))*temp(2)/abs(temp(2))

    ! Rotate second about y-axis
    L_tot_mag = sqrt(dot_product(L_tot,L_tot))
    rotate_about_y = -acos(dot_product((/0.,0.,1./),L_tot/L_tot_mag))

    call rotatevec(L_tot,(/0.,0.,1.0/),rotate_about_z)
    call rotatevec(L_tot,(/0.,1.0,0./),rotate_about_y)
    do i=1,nbin
       temp(1) = unitlx(i)
       temp(2) = unitly(i)
       temp(3) = unitlz(i)
       call rotatevec(temp,(/0.,0.,1.0/),rotate_about_z)
       call rotatevec(temp,(/0.,1.0,0./),rotate_about_y)
       unitlx(i) = temp(1)
       unitly(i) = temp(2)
       unitlz(i) = temp(3)
    enddo
 endif

 do i=1,nbin
    if (i /= 1.and.i /= nbin) then
       psi_x=(unitlx(i+1)-unitlx(i-1))/(bin(i+1)-bin(i-1))
       psi_y=(unitly(i+1)-unitly(i-1))/(bin(i+1)-bin(i-1))
       psi_z=(unitlz(i+1)-unitlz(i-1))/(bin(i+1)-bin(i-1))
       psi(i)=sqrt(psi_x**2 + psi_y**2 + psi_z**2)*bin(i)
    else
       psi(i)=0.
    endif
    if (ninbin(i) > 0) then
       tilt(i)  = acos(unitlz(i))
       tilt_acc(i) = acos(unitangz)
       twist(i) = atan2(unitly(i),unitlx(i))
       if (i==1 .or. time==0.0) then
          twistprev(i) = 0.0
       endif
       ! Taking into account negative twist
       if (twist(i) < 0.) then
          twistprev(i) = 2.*pi + twist(i)
       else
          twistprev(i) = twist(i) !cumulative twist
       endif
    else
       tilt(i) = 0.0
       twist(i) = 0.0
       dtwist = 0.0
    endif

! Calculate the precession time
    if (twistprev(i) > tiny(twistprev(i))) then
       tp(i) = time*2.*pi/twistprev(i)
    else
       tp(i) = 0.0
    endif
 enddo

 if (time == 0. .and. sorting_choice > 1) then
    print "(/,a)",' + ----------------------------------------------------- +'
    print "(a)",  ' |                                                       |'
    print "(a)",  ' |   You are sorting particles by something other        |'
    print "(a)",  ' |   than cylindrical radius. Do you really,             |'
    print "(a)",  ' |   definitely want to do this?                         |'
    print "(a)",  ' |                                                       |'
    print "(a,/)",' + ----------------------------------------------------- +'
    call prompt('Enter yes to continue:',acceptance)
    if (.not.acceptance) call exit(1)
 endif

end subroutine disc_analysis

subroutine calculate_H_fast(nbin,H,zsetgas,ninbin)
 integer, intent(in)  :: nbin
 real,    intent(out) :: H(:)
 real,    intent(in)  :: zsetgas(:,:)
 integer, intent(in)  :: ninbin(:)
 integer :: ii
 real    :: meanzii

 do ii = 1,nbin
    if (ninbin(ii)==0) then
       meanzii = 0.
    else
       meanzii = sum(zsetgas(1:ninbin(ii),ii))/real(ninbin(ii))
    endif
    H(ii) = sqrt(sum(((zsetgas(1:ninbin(ii),ii)-meanzii)**2)/(real(ninbin(ii)-1))))
 enddo

end subroutine calculate_H_fast

subroutine calculate_H_slow(nbin,npart,H,mybin,ninbin,myz)
 integer, intent(in)  :: nbin,npart
 real,    intent(out) :: H(:)
 integer, intent(in)  :: mybin(:),ninbin(:)
 real,    intent(in)  :: myz(:)
 integer :: ii,i
 real :: zsum,hsum,zmeanii

 !$omp parallel do default(none) &
 !$omp private(ii,i,zsum,zmeanii,hsum) &
 !$omp shared(nbin,npart,mybin,myz,ninbin,H)
 do ii = 1,nbin
    zsum = 0.
    ! For each bin, loop over all the particles and see if they have the correct bin
    do i = 1,npart
       if (mybin(i)==ii) then
          zsum = zsum + myz(i)
       endif
    enddo
    ! Calculat the mean z for the bin
    if (ninbin(ii)==0) then
       zmeanii = 0.
    else
       zmeanii = zsum/real(ninbin(ii))
    endif
    ! Now that we have the mean z for the bin, calculate H for each bin
    hsum = 0.
    do i = 1,npart
       if (mybin(i)==ii) then
          hsum = hsum + (myz(i) - zmeanii)**2
       endif
    enddo
    if (ninbin(ii)<1) then
       H(ii) = 0.
    else
       H(ii) = sqrt(hsum/real(ninbin(ii)-1))
    endif
 enddo
 !omp end parallel do

end subroutine calculate_H_slow

end module discanalysisutils
