!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: discanalysisutils
!
!  DESCRIPTION:
!  Routine to calculate azimuithally averaged properties in a disc
!  Can handle gas disc, gas disc + sinks, warped disc & eccentric disc
!
!  REFERENCES:
!
!  OWNER: Bec Nealon
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, externalforces, options, physcon,
!    vectorutils
!+
!--------------------------------------------------------------------------
module discanalysisutils
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

subroutine disc_analysis(xyzh,vxyz,npart,pmass,time,nr,rmin,rmax,H_R,G,M_star,q_index,&
                     tilt,tilt_acc,twistprev,psi,H,a,h_smooth,sigma,unitlx,unitly,unitlz,Lx,Ly,Lz,&
                     ecc,ninbin,assume_Ltot_is_same_as_zaxis,xyzmh_ptmass,vxyz_ptmass,nptmass)
 use physcon,      only:pi
 use centreofmass, only:get_total_angular_momentum,reset_centreofmass
 use externalforces, only:iext_einsteinprec
 use options,      only:iexternalforce
 use vectorutils, only:rotatevec
 real, intent(inout)              :: xyzh(:,:),vxyz(:,:),pmass,time
 integer, intent(in)              :: nr,npart
 real, intent(in)                 :: rmin,rmax,H_R,G,M_star,q_index
 logical, intent(in)              :: assume_Ltot_is_same_as_zaxis
 real, optional, intent(inout)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, optional, intent(inout) :: nptmass
 real, intent(out)                :: Lx(nr),Ly(nr),Lz(nr),a(nr),unitly(nr)
 real, intent(out)                :: tilt(nr),tilt_acc(nr),twistprev(nr)
 real, intent(out)                :: psi(nr),H(nr),ecc(nr),unitlz(nr)
 real, intent(out)                :: sigma(nr),h_smooth(nr),unitlx(nr)
 integer, intent(out)             :: ninbin(nr)
 real                             :: da,cs0,angx,angy,angz,unitangz
 real                             :: cs(nr),omega(nr),angtot,Ltot
 real                             :: ri,area,Ei,mu,term,ecci
 real                             :: Li(3),xi(3),vi(3),Limag,dtwist
 real                             :: psi_x,psi_y,psi_z,tp(nr)
 real                             :: L_tot(3),L_tot_mag,temp(3),temp_mag
 real                             :: rotate_about_z,rotate_about_y
 real                             :: meanzgas(nr),zdash,twist(nr),ai
 real, allocatable                :: zsetgas(:,:)
 integer                          :: i,ii
 logical                          :: rotate

! Options
 if (assume_Ltot_is_same_as_zaxis) then
    rotate = .false.
 else
    rotate = .true.
 endif

 mu = G*M_star

! Set up the radius array
 da = (rmax-rmin)/real(nr-1)
 do i=1,nr
    a(i)=rmin + real(i-1)*da
 enddo

! Initialise arrays to zero
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

! Set up cs0: cs = cs0 * R^-q
 cs0 = H_R * sqrt(G*M_star) * rmin**(q_index-0.5)

! And thus the sound speed array
 do i=1,nr
    cs(i) = cs0 * a(i)**(-q_index)
    omega(i) = sqrt(G*M_star/a(i)**3)
 enddo

 allocate(zsetgas(npart,nr))

! Move everything so that the centre of mass is at the origin
  if (nptmass > 0) then
    call reset_centreofmass(npart,xyzh,vxyz,nptmass,xyzmh_ptmass,vxyz_ptmass)
  endif

! Loop over particles putting properties into the correct bin

 do i = 1,npart

    xi = xyzh(1:3,i)
    vi = vxyz(1:3,i)

    ! i for the particle number, ii for the bin number
    ! Particles are sorted by their semi-major axis, not radius
    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE

       ri = sqrt(dot_product(xi(1:3),xi(1:3)))

       Li(1) = pmass*(xi(2)*vi(3)-xi(3)*vi(2))
       Li(2) = pmass*(xi(3)*vi(1)-xi(1)*vi(3))
       Li(3) = pmass*(xi(1)*vi(2)-xi(2)*vi(1))

       Limag = sqrt(dot_product(Li,Li))/pmass

       ! NB: No internal energy as isothermal
       if (iexternalforce==iext_einsteinprec) then
          Ei = 0.5*dot_product(vi,vi) - G*M_star/ri -3.*G*M_star/(ri**2)
          term = 2.*Ei*(Limag**2 - 6.*mu*mu)/(mu**2)
       else
          Ei = 0.5*dot_product(vi,vi) - G*M_star/ri
          term = 2.*Ei*Limag**2/(mu**2)
       endif

       ai = -M_star/(2.*Ei)
       ii = int((ai-a(1))/da + 1)

       ! If it's not in the range for the analysis, cycle
       if (ii > nr) cycle
       if (ii < 1)  cycle

       area = (pi*((a(ii)+da/2.)**2-(a(ii)- da/2.)**2))
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
 do i = 1,nr
    Ltot = sqrt(Lx(i)*Lx(i) + Ly(i)*Ly(i) + Lz(i)*Lz(i))

    unitlx(i) = Lx(i)/Ltot
    unitly(i) = Ly(i)/Ltot
    unitlz(i) = Lz(i)/Ltot

    if (ninbin(i) > 0) then
       h_smooth(i) = h_smooth(i)/ninbin(i)
       ecc(i) = ecc(i)/ninbin(i)
    endif
 enddo

 ! Now go through and get the scale height right
 ninbin = 0
 do i = 1,npart
    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE

    xi = xyzh(1:3,i)
    vi = vxyz(1:3,i)
    ri = sqrt(dot_product(xi(1:3),xi(1:3)))

    ! NB: No internal energy as isothermal
    if (iexternalforce==iext_einsteinprec) then
       Ei = 0.5*dot_product(vi,vi) - G*M_star/ri -3.*G*M_star/(ri**2)
       term = 2.*Ei*(Limag**2 - 6.*mu*mu)/(mu**2)
    else
       Ei = 0.5*dot_product(vi,vi) - G*M_star/ri
       term = 2.*Ei*Limag**2/(mu**2)
    endif

    ai = -M_star/(2.*Ei)

    ii = int((ai-a(1))/da + 1)

    if (ii > nr .or. ii < 1) cycle
    ninbin(ii) = ninbin(ii) + 1

    ! get vertical height above disc midplane == z if disc is not warped
    zdash = unitlx(ii)*xi(1) + unitly(ii)*xi(2) + unitlz(ii)*xi(3)
    zsetgas(ninbin(ii),ii) = zdash
    endif
 enddo

! Calculate H from the particle positions
 do i = 1,nr
    meanzgas(i)  = sum(zsetgas(1:ninbin(i),i))/real(ninbin(i))
    H(i) = sqrt(sum(((zsetgas(1:ninbin(i),i)-meanzgas(i))**2)/(real(ninbin(i)-1))))
 enddo

! clean up
 deallocate(zsetgas)

! Print angular momentum of accreted particles
 angtot = sqrt(angx*angx + angy*angy + angz*angz)

! For unit angular momentum accreted, z component
 unitangz = angz/angtot
 print*,' angular momentum of accreted particles = ',angtot!,angx,angy,angz,unitangz

! Now loop over rings to calculate required quantities
 do i = 1, nr
    if(ninbin(i)==0 .or. ninbin(i)==1) then
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
 if(rotate) then
    if (nptmass /= 0) then
       call get_total_angular_momentum(xyzh,vxyz,npart,L_tot,xyzmh_ptmass,vxyz_ptmass,nptmass)
    else
       call get_total_angular_momentum(xyzh,vxyz,npart,L_tot)
    endif

    temp = (/L_tot(1),L_tot(2),0./)
    temp_mag = sqrt(dot_product(temp,temp))
    rotate_about_z = acos(dot_product((/1.,0.,0./),temp/temp_mag))

    ! Rotate second about y-axis
    L_tot_mag = sqrt(dot_product(L_tot,L_tot))
    rotate_about_y = -acos(dot_product((/0.,0.,1./),L_tot/L_tot_mag))

    call rotatevec(L_tot,(/0.,0.,1.0/),-rotate_about_z)
    call rotatevec(L_tot,(/0.,1.0,0./),rotate_about_y)
    do i=1,nr
       temp(1) = unitlx(i)
       temp(2) = unitly(i)
       temp(3) = unitlz(i)
       call rotatevec(temp,(/0.,0.,1.0/),-rotate_about_z)
       call rotatevec(temp,(/0.,1.0,0./),rotate_about_y)
       unitlx(i) = temp(1)
       unitly(i) = temp(2)
       unitlz(i) = temp(3)
    enddo
 endif

 do i=1,nr
    if(i /= 1.and.i /= nr) then
       psi_x=(unitlx(i+1)-unitlx(i-1))/(a(i+1)-a(i-1))
       psi_y=(unitly(i+1)-unitly(i-1))/(a(i+1)-a(i-1))
       psi_z=(unitlz(i+1)-unitlz(i-1))/(a(i+1)-a(i-1))
       psi(i)=sqrt(psi_x**2 + psi_y**2 + psi_z**2)*a(i)
    else
       psi=0.
    endif
    if (ninbin(i) > 0) then
       tilt(i)  = acos(unitlz(i))
       tilt_acc(i) = acos(unitangz)
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

end subroutine disc_analysis

end module discanalysisutils
