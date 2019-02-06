!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for discs by RN, adapted from a routine by CJN
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, eos, infile_utils, io, part, physcon
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'MRI'
 public :: do_analysis

 integer, parameter :: nr = 300
 real,dimension(nr) :: twist,twistprev

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use dim,          only:maxp
 use physcon, only:pi
 use part, only: rhoh,Bxyz,massoftype,iphase,iamtype,igas,maxphase,mhd
 use eos,  only: equationofstate
 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: xyzh(4,npart),vxyz(3,npart)
 real,             intent(in) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 integer,parameter::nmaganalysis = 5
 character(len=9) :: output
 character(len=20) :: filename
 character(len=20) :: discprefix
 integer :: i,ii,ierr,iline
 real :: R_in,R_out,H_R,p_index,q_index,M_star
 real :: G,rmin,rmax,dr,cs0,angx,angy,angz,ri,area
 real :: angtot,Ltot,tilt,dtwist
 real :: Li(3)
 real :: rad(nr),ninbin(nr),Lx(nr),Ly(nr),Lz(nr),h_smooth(nr),sigma(nr),cs(nr),H(nr),omega(nr)
 real :: unitlx(nr),unitly(nr),unitlz(nr),tp(nr)
 real :: psi_x,psi_y,psi_z,psi,Bx,By,Bz
 real :: B_x(nr),B_y(nr),B_z(nr),B_r(nr),B_phi(nr),B_theta(nr),B_mag(nr)
 real :: rcyli,pressure,Bx_net,By_net,Bz_net
 real :: Br,Bphi,Bmag,delta_vphi,delta_vr
 real :: den_sum(nmaganalysis),num_sum(nmaganalysis)
 real :: alpha_estimates(nmaganalysis),cssqrd,ponrhoi,rho

 integer, parameter :: iparams = 10
 integer, parameter :: ialphamag   = 24
 logical :: do_alphamag,ifile

 ! Since calculating multiple versions of alpha, use
 ! PB or 1 and 2 = Parkin & Bicknell 2013 Equations 17 and 18
 ! A or 4 = Armitage et al 2001 Equation 3
 ! FN or 3 = Fromang & Nelson 2006 Equation 13
 ! G or 5 = Gaburov et al 2011 bottom of page 147
 ! NB: These only use weighted mass averages, not volume weighted (too hard basket)
 ! Calculate values across the radial bins, then sum and/or average as required, then
 ! take ratios to find alphas.


 if (mhd) then
    do_alphamag = .true.
 else
    call fatal('analysis','Cannot do MRI analysis without MHD. Use with MHD=yes.')
 endif

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'mag',numfile
 write(*,'("Output file name is ",A)') output

! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G==1")')
 G = 1.0

 iline = index(dumpfile,'_')
 discprefix = dumpfile(1:iline-1)
 inquire(file=trim(discprefix)//'.discparams', exist=ifile)
 if (ifile) then
    call read_discparams(trim(discprefix)//'.discparams',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read .discparams file')
 else
    call read_discparams('discparams.list',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read discparams.list')
 endif

! Print out the parameters
 write(*,*)
 write(*,'("Parameters are:")')
 write(*,*) 'R_in    = ',R_in
 write(*,*) 'R_out   = ',R_out
 write(*,*) 'H_R     = ',H_R
 write(*,*) 'p_index = ',p_index
 write(*,*) 'q_index = ',q_index
 write(*,*) 'M_star  = ',M_star
 write(*,*)
 write(*,*)

! Setup rmin and rmax for the analysis
 rmin = R_in
 rmax = R_out

! Extended setup for magnetic discs - prevents NaNs in l (because the disc spreads?)
 if (do_alphamag) then
    rmin = 0.0
    rmax = 2.*R_out
 endif

! Setup rmin and rmax for the analysis
 rmin = R_in
 rmax = R_out

! Set up the radius array
 dr = (rmax-rmin)/real(nr-1)
 do i=1,nr
    rad(i)=rmin + real(i-1)*dr
 enddo

! Initialise arrays to zero
 ninbin(:)=0
 lx(:)=0.0
 ly(:)=0.0
 lz(:)=0.0
 h_smooth(:)=0.0
 sigma(:)=0.0
 B_x(:) = 0.0
 B_y(:) = 0.0
 B_z(:) = 0.0
 B_r(:) = 0.0
 B_phi(:) = 0.0
 B_theta(:) = 0.0
 B_mag(:) = 0.0
 Bx_net = 0.0
 By_net = 0.0
 Bz_net = 0.0

! Set up cs0: cs = cs0 * R^-q
 cs0 = H_R * sqrt(G*M_star) * R_in**(q_index-0.5)

! And thus the sound speed array
 do i=1,nr
    cs(i) = cs0 * rad(i)**(-q_index)
    omega(i) = sqrt(G*M_star/rad(i)**3)
 enddo

! and thus the disc scale height
 do i=1,nr
    H(i) = cs(i)/omega(i)
 enddo

 angx = 0.0
 angy = 0.0
 angz = 0.0

! Loop over particles putting properties into the correct bin
 do i = 1,npart

! i refers to particle, ii refers to bin
    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       ri = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
       rcyli = sqrt(dot_product(xyzh(1:2,i),xyzh(1:2,i)))
       ii = int((ri-rad(1))/dr + 1)

       if (ii > nr) cycle
       if (ii < 1)  cycle

       area = (pi*((rad(ii)+dr/2.)**2-(rad(ii)- dr/2.)**2))
       sigma(ii) = sigma(ii) + pmass/area

       Li(1) = pmass*(xyzh(2,i)*vxyz(3,i)-xyzh(3,i)*vxyz(2,i))
       Li(2) = pmass*(xyzh(3,i)*vxyz(1,i)-xyzh(1,i)*vxyz(3,i))
       Li(3) = pmass*(xyzh(1,i)*vxyz(2,i)-xyzh(2,i)*vxyz(1,i))

       Lx(ii)=Lx(ii)+Li(1)
       Ly(ii)=Ly(ii)+Li(2)
       Lz(ii)=Lz(ii)+Li(3)

       h_smooth(ii) = h_smooth(ii) + xyzh(4,i)

       Bx = Bx + Bxyz(1,i)
       By = Bxyz(2,i)
       Bz = Bxyz(3,i)
       Bmag = sqrt(Bx**2 + By**2 + Bz**2)

       B_x(ii) = B_x(ii) + Bx
       B_y(ii) = B_y(ii) + By
       B_z(ii) = B_z(ii) + Bz
       B_mag(ii) = B_mag(ii) + Bmag

       B_r(ii) = B_r(ii) + (Bx*xyzh(1,i)/ri) + (By*xyzh(2,i)/ri) + (Bz*xyzh(3,i)/ri)
       B_phi(ii) = B_phi(ii) - (Bx*xyzh(2,i)/rcyli) + (By*xyzh(1,i)/rcyli)
       B_theta(ii) = B_theta(ii) + (Bx*xyzh(1,i)*xyzh(3,i)/(ri*rcyli))  &
                     + (By*xyzh(2,i)*xyzh(3,i)/(ri*rcyli)) - (Bz*rcyli/ri)

       ninbin(ii) = ninbin(ii) + 1

    elseif (xyzh(4,i) < -tiny(xyzh)) then !ACCRETED
       angx = angx + pmass*(xyzh(2,i)*vxyz(3,i) - xyzh(3,i)*vxyz(2,i))
       angy = angy + pmass*(xyzh(3,i)*vxyz(1,i) - xyzh(1,i)*vxyz(3,i))
       angz = angz + pmass*(xyzh(1,i)*vxyz(2,i) - xyzh(2,i)*vxyz(1,i))
    endif
 enddo

! Print angular momentum of accreted particles
 angtot = sqrt(angx*angx + angy*angy + angz*angz)
 print*,' angular momentum of accreted particles = ',angtot

! Convert total angular momentum into a unit vector, and average h_smooth
 do i = 1,nr
    Ltot = sqrt(Lx(i)*Lx(i) + Ly(i)*Ly(i) + Lz(i)*Lz(i))

    unitlx(i) = Lx(i)/Ltot
    unitly(i) = Ly(i)/Ltot
    unitlz(i) = Lz(i)/Ltot

    if (ninbin(i) > 0) then
       h_smooth(i) = h_smooth(i)/ninbin(i)
       B_x(i) = B_x(i)/ninbin(i)
       B_y(i) = B_y(i)/ninbin(i)
       B_z(i) = B_z(i)/ninbin(i)
       B_r(i) = B_r(i)/ninbin(i)
       B_phi(i) = B_phi(i)/ninbin(i)
       B_theta(i) = B_theta(i)/ninbin(i)
       B_mag(i) = B_mag(i)/ninbin(i)
    endif
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

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',17(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'radius', &
       2,'sigma', &
       3,'<h>/H', &
       4,'lx', &
       5,'ly', &
       6,'lz', &
       7,'tilt', &
       8,'twist', &
       9,'psi', &
       10,'H/R', &
       11,'Bx', &
       12,'By', &
       13,'Bz', &
       14,'Br', &
       15,'Bphi', &
       16,'Btheta', &
       17,'|B|'

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

    if (ninbin(i) > 0) then
       write(iunit,'(17(es18.10,1X))') rad(i),sigma(i),h_smooth(i),unitlx(i),unitly(i),unitlz(i),&
                                       tilt,twist(i),psi,H(i)/rad(i),B_x(i),B_y(i),B_z(i),B_r(i),&
                                       B_phi(i),B_theta(i),B_mag(i)
    endif

 enddo

! Calculate the required properties for each particle, sum to the necessary place
 den_sum(:) = 0.
 num_sum(:) = 0.

 do i=1,npart
    ri = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
    rcyli = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
    Bx = Bxyz(1,i)
    By = Bxyz(2,i)
    Bz = Bxyz(3,i)
    Bmag = sqrt(Bx**2 + By**2 + Bz**2)

    Bx_net = Bx_net + Bx
    By_net = By_net + By
    Bz_net = Bz_net + Bz

    rho = rhoh(xyzh(4,i),pmass)
    call equationofstate(3,ponrhoi,cssqrd,rho,xyzh(1,i),xyzh(2,i),xyzh(3,i))
    pressure = cssqrd*rho

    Br = (Bx*xyzh(1,i)/ri) + (By*xyzh(2,i)/ri) + (Bz*xyzh(3,i)/ri)
    Bphi = (-Bx*xyzh(2,i)/rcyli) + (By*xyzh(1,i)/rcyli)

    delta_vr = (vxyz(1,i)*xyzh(1,i)/ri) + (vxyz(2,i)*xyzh(2,i)/ri) + (vxyz(3,i)*xyzh(3,i)/ri)
    delta_vphi = (-vxyz(1,i)*xyzh(2,i)/rcyli) + (vxyz(2,i)*xyzh(1,i)/rcyli) - (1./sqrt(ri))

    num_sum(1) = num_sum(1) + pmass*(rho*delta_vr*delta_vphi - Br*Bphi)
    den_sum(1) = den_sum(1) + pmass*pressure

    num_sum(2) = num_sum(2) - 2.*pmass*Br*Bphi
    den_sum(2) = den_sum(2) + pmass*Bmag**2

    num_sum(3) = num_sum(3) + pmass*(delta_vphi*delta_vr - (Bphi*Br/rho))
    den_sum(3) = den_sum(3) + pmass*pressure

    num_sum(4) = num_sum(4) - 2./3.*Br*Bphi/pressure
    den_sum(4) = den_sum(4) + pmass

 enddo

 num_sum(5) = num_sum(2)
 den_sum(5) = den_sum(1)

 do i=1,nmaganalysis
    if (den_sum(i) > tiny(den_sum(i))) then
       alpha_estimates(i) = num_sum(i)/den_sum(i)
    else
       alpha_estimates(i) = 7.0
    endif
 enddo

 if (do_alphamag) then
    write(filename,"(a,i3.3)")"alpha_mag"
    if (numfile==0) then
       open(unit=ialphamag,file=filename,status="replace")
       write(ialphamag,'("# various alpha_mag measures averaged over disc",es18.10)')
       write(ialphamag,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
               1,'time', &
               2,'PB alpha_p', &
               3,'PB alpha_m', &
               4,'FN06', &
               5,'Armitage', &
               6,'Guburov'
    else
       open(unit=ialphamag,file=filename,status="old",access="append")
    endif
    write(ialphamag,'(6(es18.10,1X))') time,alpha_estimates(1),alpha_estimates(2),&
                                          alpha_estimates(3),alpha_estimates(4),alpha_estimates(5)
    close(unit=ialphamag)
 endif

 close(iunit)

 ! Check that the net magnetic field is roughly zero (less than 10% of field)
 if (abs(Bx_net)/npart > 0.1*maxval(abs(B_x)) .or. abs(By_net)/npart > 0.1*maxval(abs(B_y)) &
      .or. abs(Bz_net)/npart > 0.1*maxval(abs(B_z))) then
    print*,'WARNING! Net field is > 10% of max.'
    print*,'In x, net is ',abs(Bx_net)/npart,'compared to max of',maxval(abs(B_x))
    print*,'In y, net is ',abs(By_net)/npart,'compared to max of',maxval(abs(B_y))
    print*,'In z, net is ',abs(Bz_net)/npart,'compared to max of',maxval(abs(B_z))
 endif

end subroutine do_analysis

!----------------------------------------------------------------
!+
!  Read disc information from discparams.list file
!+
!----------------------------------------------------------------
subroutine read_discparams(filename,R_in,R_out,H_R,p_index,q_index,M_star,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: R_in,R_out,H_R,p_index,q_index,M_star
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
 call read_inopt(H_R,'H_R',db,ierr)
 if (ierr /= 0) return
 call read_inopt(p_index,'p_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(q_index,'q_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_star,'M_star',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_discparams

end module analysis

