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
!  Analysis routine for discs
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: infile_utils, io, physcon
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'CJN'
 public :: do_analysis, disc_analysis

 integer, parameter :: nr = 300
 real,dimension(nr) :: twist,twistprev

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=9) :: output
 character(len=20) :: filename
 character(len=20) :: discprefix
 integer :: i,ierr,iline
 real :: R_in,R_out,H_R,p_index,q_index,M_star
 real :: G,rmin,rmax
 real :: tilt(nr)
 real :: rad(nr),ninbin(nr),h_smooth(nr),sigma(nr),H(nr)
 real :: unitlx(nr),unitly(nr),unitlz(nr),tp(nr),ecc(nr)
 real :: psi(nr),tilt_acc(nr)

 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24
 logical :: do_precession,ifile

 do_precession = .false.

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'angm',numfile
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
 write(*,*) 'H/R_ref = ',H_R
 write(*,*) 'p_index = ',p_index
 write(*,*) 'q_index = ',q_index
 write(*,*) 'M_star  = ',M_star
 write(*,*)
 write(*,*)

! Setup rmin and rmax for the analysis
 rmin = R_in
 rmax = R_out

 call disc_analysis(xyzh,vxyz,npart,pmass,time,nr,rmin,rmax,H_R,G,M_star,q_index,&
                     tilt,tilt_acc,tp,psi,H,rad,h_smooth,sigma,unitlx,unitly,unitlz,ecc,ninbin)

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
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
       11,'|e|'

 do i = 1,nr
    if (ninbin(i) > 0) then
       write(iunit,'(13(es18.10,1X))') rad(i),sigma(i),h_smooth(i),unitlx(i),unitly(i),unitlz(i),&
                                         tilt(i),twist(i),psi(i),H(i)/rad(i),ecc(i)
    endif



! Printing time and twist for each radius bin
    if (do_precession) then
       write(filename,"(a,i3.3)")"precess",i
       if (numfile==0) then
          open(unit=iprec,file=filename,status="replace")
          write(iprec,'("# tilt and twist with time for r = ",es18.10)') rad(i)
          write(iprec,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
               1,'rad', &
               2,'time', &
               3,'tilt', &
               4,'twist', &
               5,'tot twist', &
               6,'tp', &
               7,'|e|'
       else
          open(unit=iprec,file=filename,status="old",position="append")
       endif
       write(iprec,'(7(es18.10,1X))') rad(i),time,tilt(i),twist(i),twistprev(i),tp(i),ecc(i)
       !      write(iprec,'(7(es18.10,1X))') rad(i),time,tilt,tilt_acc,twistprev(i),tp(i),ecc(i)
       !      print*,tilt,twist(i),twistprev(i)
       close(unit=iprec)
    endif

 enddo

 close(iunit)

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
 call read_inopt(H_R,'H/R_ref',db,ierr)
 if (ierr /= 0) return
 call read_inopt(p_index,'p_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(q_index,'q_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_star,'M_star',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_discparams

!----------------------------------------------------------------
!+
!  Disc analysis routine - so that this can be called externally
!  Adapted from a routine written by Chris Nixon
!+
!----------------------------------------------------------------

subroutine disc_analysis(xyzh,vxyz,npart,pmass,time,nr,rmin,rmax,H_R,G,M_star,q_index,&
                     tilt,tilt_acc,tp,psi,H,rad,h_smooth,sigma,unitlx,unitly,unitlz,ecc,ninbin)
 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use physcon,      only:pi
 use centreofmass, only:get_total_angular_momentum
 use options,      only:iexternalforce
 use vectorutils, only:rotatevec
 real, intent(inout) :: xyzh(:,:)
 real, intent(inout) :: vxyz(:,:)
 real, intent(inout) :: pmass,time
 integer, intent(in) :: nr,npart
 real, intent(in)    :: rmin,rmax,H_R,G,M_star,q_index
 real                :: dr,cs0,angx,angy,angz,unitangz
 real                :: Lx(nr),Ly(nr),Lz(nr)
 real                :: cs(nr),omega(nr),angtot,Ltot
 real                :: ri,area,Ei,mu,term,ecci
 real                :: Li(3),xi(3),vi(3),Limag,dtwist,psi_x,psi_y,psi_z
 real, intent(out)   :: tilt(nr),tilt_acc(nr),tp(nr),psi(nr),H(nr),ecc(nr),rad(nr)
 real, intent(out)   :: sigma(nr),h_smooth(nr),unitlx(nr),unitly(nr),unitlz(nr),ninbin(nr)
 real :: L_tot(3),L_tot_mag,temp(3),temp_mag,rotate_about_z,rotate_about_y
 integer             :: i,ii

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
 ecc=0.0

! Set up cs0: cs = cs0 * R^-q
 cs0 = H_R * sqrt(G*M_star) * rmin**(q_index-0.5)

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

    if (xyzh(4,i)  >  tiny(xyzh)) then ! IF ACTIVE
       ri = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
       ii = int((ri-rad(1))/dr + 1)
       xi = xyzh(1:3,i)
       vi = vxyz(1:3,i)

       if (ii > nr) cycle
       if (ii < 1)  cycle

       area = (pi*((rad(ii)+dr/2.)**2-(rad(ii)- dr/2.)**2))
       sigma(ii) = sigma(ii) + pmass/area

       Li(1) = pmass*(xyzh(2,i)*vxyz(3,i)-xyzh(3,i)*vxyz(2,i))
       Li(2) = pmass*(xyzh(3,i)*vxyz(1,i)-xyzh(1,i)*vxyz(3,i))
       Li(3) = pmass*(xyzh(1,i)*vxyz(2,i)-xyzh(2,i)*vxyz(1,i))

       Limag = sqrt(dot_product(Li,Li))/pmass
       ! Energy for standard potential, if using Einstein precession uncomment last bit and change 'term'
       ! NB: No internal energy as isothermal
       Ei = 0.5*dot_product(vi,vi) - G*M_star/ri! -3.*G*M_star/(ri**2)
       mu = G*M_star
       term = 2.*Ei*Limag**2/(mu**2)
       !term = 2.*Ei*(Limag**2 - 6.*mu*mu)/(mu**2)

       ecci = sqrt(1. + term)

       Lx(ii)=Lx(ii)+Li(1)
       Ly(ii)=Ly(ii)+Li(2)
       Lz(ii)=Lz(ii)+Li(3)
       ecc(ii) = ecc(ii) + ecci
       h_smooth(ii) = h_smooth(ii) + xyzh(4,i)

       ninbin(ii) = ninbin(ii) + 1

    elseif (xyzh(4,i) < -tiny(xyzh)) then !ACCRETED
       angx = angx + pmass*(xyzh(2,i)*vxyz(3,i) - xyzh(3,i)*vxyz(2,i))
       angy = angy + pmass*(xyzh(3,i)*vxyz(1,i) - xyzh(1,i)*vxyz(3,i))
       angz = angz + pmass*(xyzh(1,i)*vxyz(2,i) - xyzh(2,i)*vxyz(1,i))
    endif
 enddo

! Print angular momentum of accreted particles
 angtot = sqrt(angx*angx + angy*angy + angz*angz)
! For unit angular momentum accreted, z component
 unitangz = angz/angtot
 print*,' angular momentum of accreted particles = ',angtot!,angx,angy,angz,unitangz

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

 ! Calculate the total angular momentum vector and rotate unitl[x,y,z] if required
 if(iexternalforce == 0) then
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
       psi_x=(unitlx(i+1)-unitlx(i-1))/(rad(i+1)-rad(i-1))
       psi_y=(unitly(i+1)-unitly(i-1))/(rad(i+1)-rad(i-1))
       psi_z=(unitlz(i+1)-unitlz(i-1))/(rad(i+1)-rad(i-1))
       psi(i)=sqrt(psi_x**2 + psi_y**2 + psi_z**2)*rad(i)
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

end module analysis

