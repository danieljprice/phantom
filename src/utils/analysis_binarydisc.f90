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
!  Analysis routine for discs in binary systems by CJN
!
! Method:
! "primary(secondary)" implies tag=1(2) not anything to do with mass or disc hosting etc.
! Circumprimary   disc: only particles inside the primary   Roche-lobe & bound
! Circumsecondary disc: only particles inside the secondary Roche-lobe & bound
! Circumbinary    disc: only particles at R > a (semi-major axis, not separation) & bound
!
! Creates: binary.dat      for binary params
!          circumprimary_*.dat   for circumprimary disc params
!          circumsecondary_*.dat for circumsecondary disc params
!          circumbinary_*.dat    for circumbinary disc params
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, infile_utils, io, options, part, physcon, setbinary
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'binaryanalysis'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

 use io,       only:fatal
 use part,     only:xyzmh_ptmass,vxyz_ptmass,ihacc
 use physcon,  only:pi
 use setbinary,only:Rochelobe_estimate
 use eos,      only:gamma
 use options,  only:ieos

 integer, parameter :: ngrid = 100
 integer, parameter :: iout  = 149
 integer, parameter :: iparams = 150
 integer, parameter :: itime = 151
 character(len=*), intent(in) :: dumpfile
 real,dimension(:,:), intent(in) :: xyzh,vxyzu
 real, intent(in) :: pmass,time
 integer, intent(in) :: npart,iunit,numfile

 character(len=50) :: output, timeoutput, filename
 character(len=15),dimension(3) :: disc_type

 logical :: use_a = .false.
 logical :: exists

 integer :: idot
 integer :: nsinks,i,j,ipri,isec,ibin,nunbound,icell,ierr
 integer,dimension(3) :: n
 integer,dimension(npart) :: imysink
 integer :: check
 integer :: i1,i2,i3

 real,parameter :: G = 1.0
 real :: ui,ri,Ei,ai,a,ecc,dr,csi,ecci,xprir,xsecr,Limag,Hi,Ltot,omegai
 real :: psi_x,psi_y,psi_z,psi,rmax,rmin,tilt,twist,Lmagi,ai_cell
 real,dimension(3) :: xi,vi,Li,mptmass,xcom,xpricom,xseccom,rtest
 real,dimension(ngrid) :: rad,sigma,h_smooth,ecc_cell,ecc2_cell
 real,dimension(ngrid) :: ninbin,cs_cell,omega_cell,H_cell,area,E_cell,zgas,zgas2,meanzgas,hgas
 real,dimension(3,ngrid) :: L_cell,unitL_cell
 real,dimension(3,3) :: xptmass,vptmass
 real :: ecc1,ecc2,ecc3
 real :: inc1,inc2,inc3
 real :: r1,r2, r3

! Use two variables solely to remove compiler warnings...
 write(*,'("Performing analysis on ",a,"... which is unit ",i5,"...")') trim(dumpfile),iunit
 idot = index(dumpfile,'_') - 1
 filename = dumpfile(1:idot)  !create filename

! read the .in file
 call read_dotin(''//trim(filename)//'.in',ieos,iparams,ierr)
 if (ierr /= 0) call fatal('analysis','could not open/read disc.in')

 if (gamma < tiny(gamma)) call fatal(analysistype,'gamma not set...',var='gamma',val=gamma)

! Check nsinks is >=2
 nsinks = size(xyzmh_ptmass(1,:))
 if (nsinks < 2) call fatal(analysistype,'Not enough sinks...')

! Labels for binary - in principle allows more than 2 sinks
 ipri = 1
 isec = 2
 ibin = 3

 ecc1 = -1.0
 ecc2 = -1.0
 ecc3 = -1.0
 inc1 = -1.0
 inc2 = -1.0
 inc3 = -1.0

 write(disc_type(1),'("circumprimary")')
 write(disc_type(2),'("circumsecondary")')
 write(disc_type(3),'("circumbinary")')

! Calculate binary params (and write to binary.dat):
 call get_binary_params(ipri,isec,xyzmh_ptmass,vxyz_ptmass,time,a,ecc,G)
 mptmass(ipri) = xyzmh_ptmass(4,ipri)
 mptmass(isec) = xyzmh_ptmass(4,isec)
 mptmass(ibin) = mptmass(ipri) + mptmass(isec)

 xptmass(:,ipri) = xyzmh_ptmass(1:3,ipri)
 xptmass(:,isec) = xyzmh_ptmass(1:3,isec)
 xptmass(:,ibin) = (mptmass(ipri)*xptmass(:,ipri) + mptmass(isec)*xptmass(:,isec))/mptmass(ibin)
 a=sqrt(dot_product((xptmass(:,ipri)-xptmass(:,isec)),(xptmass(:,ipri)-xptmass(:,isec))))

 vptmass(:,ipri) = vxyz_ptmass(1:3,ipri)
 vptmass(:,isec) = vxyz_ptmass(1:3,isec)
 vptmass(:,ibin) = (mptmass(ipri)*vptmass(:,ipri) + mptmass(isec)*vptmass(:,isec))/mptmass(ibin)

! location of centre of mass of binary
 xcom(:) = xptmass(:,ibin)

! Calculate Roche lobe size of primary - based on semi-major axis, not instantaneous separation
 xpricom(:) = xptmass(1:3,ipri) - xcom(:)
 xseccom(:) = xptmass(1:3,isec) - xcom(:)
 xprir = sqrt(dot_product(xpricom,xpricom)) ! distance from primary to c. o. m.
 xsecr = sqrt(dot_product(xseccom,xseccom)) ! distance from secondary to c. o. m.
 rtest(ipri) = Rochelobe_estimate(mptmass(isec),mptmass(ipri),a) ! max radius of circumprimary disc
 rtest(isec) = Rochelobe_estimate(mptmass(ipri),mptmass(isec),a) ! max radius of circumsecondary disc
 rtest(ibin) = max(a,xprir+rtest(ipri),xsecr+rtest(isec)) ! this is an underestimate of the circumbinary inner edge...


 write(*,'("Using ieos: ",i2.1)') ieos
 write(*,'("Roche-lobe of primary is: ",es17.10)') Rochelobe_estimate(mptmass(isec),mptmass(ipri),a)
 write(*,'("Roche-lobe of secondary is: ",es17.10)') Rochelobe_estimate(mptmass(ipri),mptmass(isec),a)

 write(*,'("Rejecting particles that are a distance ",es17.10, " from the primary")') rtest(ipri)
 write(*,'("Rejecting particles that are a distance ",es17.10, " from the secondary")') rtest(isec)

 n(:) = 0 !number of particles in each 'disc'
 nunbound = 0
 imysink(:) = 0

! FIRST LOOP: tag particles for each disc
 do j=1,3 ! SINK LOOP
    write(*,'("Analysing ",a15)') disc_type(j)

    sigma(:)      = 0.0
    L_cell(:,:)   = 0.0
    h_smooth(:)   = 0.0
    cs_cell(:)    = 0.0
    omega_cell(:) = 0.0
    H_cell(:)     = 0.0
    ecc_cell(:)   = 0.0
    E_cell(:)     = 0.0
    zgas(:)       = 0.0
    zgas2(:)      = 0.0
    meanzgas(:)   = 0.0
    ninbin(:) = 0

    rmin = huge(rmin)
    rmax = tiny(rmax)

    do i=1,npart ! PARTICLE LOOP
! First skip dead particles
       if (xyzh(4,i) < tiny(xyzh)) cycle

! Second check if particles in right place
       xi(1:3) = xyzh(1:3,i) - xptmass(1:3,j)
       ri = sqrt(dot_product(xi,xi))
       if (j /= ibin .and. ri > rtest(j)) cycle ! If pri/sec then check particles are inside RL
       if (j == ibin .and. ri < rtest(j)) cycle ! If bin then check particles are outside binary

! Local storage of particle velocity and thermal energy
       vi(1:3) = vxyzu(1:3,i) - vptmass(1:3,j)
       if (size(vxyzu(:,1)) == 4) then
          ui = vxyzu(4,i)
       else
          call get_utherm(ieos,xi(1),xi(2),xi(3),gamma,ui,csi)
       endif

! By this point, particles are interesting (i.e. in the right area).
! Now check they are also bound to their sink (i.e. specific energy negative):
       Ei = get_particle_energy(G,mptmass(j),ri,vi,ui)
       if (Ei > 0.0) nunbound = nunbound + 1
       if (Ei > 0.0) cycle !unbound

! Calculate particle eccentricity and semi-major axis
       call cross(xi,vi,Li)
       Limag = sqrt(dot_product(Li,Li))
       if (abs(Ei) < tiny(Ei)) stop 'initial particle energy problem'
       call get_ae(Limag,Ei,mptmass(j),pmass,ai,ecci)

       imysink(i) = j
       n(j) = n(j) + 1

! min and max of semi-major axis for grid
       if (use_a) then
          rmin = min(rmin,ai)
          rmax = max(rmax,ai)
       else
          rmin = min(rmin,ri)
          rmax = max(rmax,ri)
       endif

    enddo ! END PARTICLE LOOP
    rmin = 0.99*rmin
    rmax = 1.01*rmax

    if (j == ipri) then
       r1=0.1*a
       r2=0.2*a
       r3=0.3*a

       i1=huge(i1)
       i2=huge(i2)
       i3=huge(i3)
    endif

! Set up the radius array
    dr = (rmax-rmin)/real(ngrid-1)
    do i=1,ngrid
       rad(i)=rmin + real(i-1)*dr
       area(i) = (pi*((rad(i)+dr/2.0)**2-(rad(i)- dr/2.0)**2))

       if (j == ipri .and. i /= 1) then
          if (rad(i-1) < r1 .and. rad(i) > r1) i1=i
          if (rad(i-1) < r2 .and. rad(i) > r2) i2=i
          if (rad(i-1) < r3 .and. rad(i) > r3) i3=i
       endif

    enddo

!    write(*,'("Setting up grid from ",es17.10," to ",es17.10," with ",i3," points...")') rmin,rmax,ngrid


! SECOND LOOP: calculate particle properties and place into bins
    do i=1,npart ! SECOND PARTICLE LOOP
       if (imysink(i) /= j) cycle

       xi(1:3) = xyzh(1:3,i) - xptmass(1:3,j)
       vi(1:3) = vxyzu(1:3,i) - vptmass(1:3,j)
       ri = sqrt(dot_product(xi,xi))
       if (ri < rmin .or. ri > rmax) cycle
       if (size(vxyzu(:,1)) == 4) then
          ui = vxyzu(4,i)
          csi = sqrt(gamma*(gamma-1)*ui)
       else
          call get_utherm(ieos,xi(1),xi(2),xi(3),gamma,ui,csi)
       endif
       Ei = get_particle_energy(G,mptmass(j),ri,vi,ui)
       if (abs(Ei) < tiny(Ei)) stop 'Ei is zero...'

! Calculate particle eccentricity and semi-major axis
       call cross(xi,vi,Li)
       Limag = sqrt(dot_product(Li,Li))
       call get_ae(Limag,Ei,mptmass(j),pmass,ai,ecci)

       if (use_a) then
          icell = int((ai-rad(1))/dr) + 1
       else
          icell = int((ri-rad(1))/dr) + 1
       endif
       if (icell > ngrid .or. icell < 1) call fatal(analysistype,'particles where they shouldn''t be...')

! Angular momentum with pmass for later
       Li(:) = pmass*Li(:)

! Omegai and hence Hi
       if (use_a) then
          omegai = sqrt(G*mptmass(j)/ai**3)
       else
          omegai = sqrt(G*mptmass(j)/ri**3)
       endif
       Hi     = csi/omegai

! Add particle to icell
       sigma(icell)      = sigma(icell)      + pmass/area(icell)
       L_cell(:,icell)   = L_cell(:,icell)   + Li(:)
       h_smooth(icell)   = h_smooth(icell)   + xyzh(4,i)
       cs_cell(icell)    = cs_cell(icell)    + csi
       omega_cell(icell) = omega_cell(icell) + omegai
       H_cell(icell)     = H_cell(icell)     + Hi
       ecc_cell(icell)   = ecc_cell(icell)   + ecci
       E_cell(icell)     = E_cell(icell)     + Ei
       zgas(icell)       = zgas(icell)       + xi(3)
       zgas2(icell)      = zgas2(icell)      + xi(3)**2

       ninbin(icell) = ninbin(icell) + 1
    enddo

! Average arrays that need it
    do i=1,ngrid
       if (ninbin(i) /= 0) then
          h_smooth(i)    = h_smooth(i)/ninbin(i)
          cs_cell(i)     = cs_cell(i)/ninbin(i)
          omega_cell(i)  = omega_cell(i)/ninbin(i)
          H_cell(i)      = H_cell(i)/ninbin(i)
          ecc_cell(i)    = ecc_cell(i)/ninbin(i)
          E_cell(i)    = E_cell(i)/ninbin(i)
       endif
    enddo

!-------------------------------------
! Compute disc thickness using std dev
!-------------------------------------


    do i=1,ngrid
       if (ninbin(i)/=0) then
          meanzgas(i)=zgas(i)/real(ninbin(i))
          hgas(i)=sqrt((zgas2(i)-2*meanzgas(i)*zgas(i)+ninbin(i)*meanzgas(i)**2)/(real(ninbin(i)-1)))
       endif
    enddo

    if (j == ipri) then
       write(output,'(a13,"_",i5.5,".dat")') disc_type(j),numfile
    elseif (j == isec) then
       write(output,'(a15,"_",i5.5,".dat")') disc_type(j),numfile
    elseif (j == ibin) then
       write(output,'(a12,"_",i5.5,".dat")') disc_type(j),numfile
    else
       call fatal(analysistype,'j is some funny number...')
    endif

    open(iout,file=trim(output))
    write(iout,'("# Analysis data at t = ",es20.12)') time
    write(iout,"('#',13(1x,'[',i2.2,1x,a11,']',2x))") &
         1,'radius', &
         2,'sigma', &
         3,'<h>/H', &
         4,'lx', &
         5,'ly', &
         6,'lz', &
         7,'tilt', &
         8,'twist', &
         9,'psi', &
         10,'ecc', &
         11,'ecc2', &
         12,'H_R', &
         13,'ninbin'

! Work out unitL array
    do i=1,ngrid
       if (ninbin(i) /= 0) then
          Ltot = sqrt(dot_product(L_cell(:,i),L_cell(:,i)))
          unitL_cell(:,i) = L_cell(:,i)/Ltot
       endif
    enddo

! Loop, calculate final quantities and then print to file
    do i=1,ngrid
       if (ninbin(i) /= 0) then
          if (i /= 1 .and. i /= ngrid) then
             if (ninbin(i-1) /= 0 .and. ninbin(i+1) /= 0) then
                psi_x=(unitL_cell(1,i+1)-unitL_cell(1,i-1))/(rad(i+1)-rad(i-1))
                psi_y=(unitL_cell(2,i+1)-unitL_cell(2,i-1))/(rad(i+1)-rad(i-1))
                psi_z=(unitL_cell(3,i+1)-unitL_cell(3,i-1))/(rad(i+1)-rad(i-1))
                psi=sqrt(psi_x**2 + psi_y**2 + psi_z**2)*rad(i)
             else
                psi=0.
             endif
          else
             psi=0.
          endif
          if (psi /= psi) then
             print*, 'Error with psi... at i = ',i
             print*, psi_x,psi_y,psi_z
             print*, rad(i-1),rad(i),rad(i+1)
             print*, unitL_cell(1,i-1),unitL_cell(2,i-1),unitL_cell(3,i-1)
             print*, unitL_cell(1,i+1),unitL_cell(2,i+1),unitL_cell(3,i+1)
             stop
          endif
       else
          unitL_cell(:,i) = 0.0
          psi = 0.0
       endif

       if (ninbin(i) > 0) then
          tilt  = acos(unitL_cell(3,i))
          twist = atan2(unitL_cell(2,i),unitL_cell(1,i))
       else
          tilt = 0.0
          twist = 0.0
       endif

       if (ninbin(i) > 0) then
          Lmagi = sqrt(dot_product(L_cell(:,i),L_cell(:,i)))/(pmass*real(ninbin(i)))
          call get_ae(Lmagi,E_cell(i),mptmass(j),pmass,ai_cell,ecc2_cell(i))
       else
          ecc2_cell(i) = 0.0
       endif

! PUT A CHECK HERE THAT ai_cell is consistent with rad(i)...

       if (ninbin(i) > 0) then
          write(iout,'(13(es18.10,1X))') rad(i),sigma(i),h_smooth(i)/H_cell(i), &
               unitL_cell(1,i),unitL_cell(2,i),unitL_cell(3,i), &
               tilt,twist,psi,ecc_cell(i),ecc2_cell(i),hgas(i)/rad(i),real(ninbin(i))
       endif

       if (j == ipri .and. i1 /= huge(i1)) then
          ecc1=ecc_cell(i1)
          inc1=acos(unitL_cell(3,i1))
       endif
       if (j == ipri .and. i2 /= huge(i2)) then
          ecc2=ecc_cell(i2)
          inc2=acos(unitL_cell(3,i2))
       endif
       if (j == ipri .and. i3 /= huge(i3)) then
          ecc3=ecc_cell(i3)
          inc3=acos(unitL_cell(3,i3))
       endif

    enddo
    close(iout)

 enddo ! END SINK LOOP

 write(timeoutput,'("time.dat")')
 if (time <= tiny(time)) then
    open(itime,file=timeoutput,status='replace',action='write',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open time.dat file at t=0.0')
    write(itime,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
         1,'time', &
         2,'e_1', &
         3,'e_2', &
         4,'e_3', &
         5,'i_1', &
         6,'i_2', &
         7,'i_3'
 else
    inquire(file=timeoutput,exist=exists)
    if (.not. exists) call fatal(analysistype,'t /= 0.0, but the time analysis output file does not exist...')
    open(itime,file=timeoutput,status='old',action='write',position='append',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open time.dat file during run')
 endif
 write(itime,'(7(ES18.10,1X))') time,ecc1,ecc2,ecc3,inc1,inc2,inc3
 close(itime)

 write(*,'("Number of particles in primary,secondary and binary discs: ",i9,", ",i9," and ",i9)') n(1),n(2),n(3)
 write(*,'("Number of unbound particles: ",i6)') nunbound
 write(*,'("Finished analysis.")')
 write(*,'(" ")')
 write(*,'(" ")')

end subroutine do_analysis
!----------------------------------------------------------------
!+
!  Read information from .in file
!+
!----------------------------------------------------------------
subroutine read_dotin(filename,ieos,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 implicit none
 character(len=*), intent(in)  :: filename
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ieos,ierr
 type(inopts), dimension(:), allocatable :: db

! Read in parameters from the .in file
 open(unit=iunit,file=filename,status='old',iostat=ierr,form='formatted')
 if (ierr /= 0) return

 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(ieos,'ieos',db,ierr)
 call close_db(db)
 close(iunit)

end subroutine read_dotin
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
! Calculate the binary paramters
!-----------------------------------------------------------------------
subroutine get_binary_params(ipri,isec,xyzmh_ptmass,vxyz_ptmass,time,a,ecc,G)
!-----------------------------------------------------------------------
 use io, only:fatal

 implicit none

 integer,intent(in) :: ipri,isec
 real,intent(in) :: time,G
 real,dimension(:,:),intent(in) :: xyzmh_ptmass,vxyz_ptmass
 real,intent(out) :: a,ecc

 logical :: exists
 character(len=10) :: output
 integer,parameter :: iunit = 150
 integer :: check
 real :: rbin,mpri,msec,E,Lmag
 real,dimension(3) :: xpri,vpri,xsec,vsec,dr,dv,L

 write(output,"(a10)") 'binary.dat'

 mpri = xyzmh_ptmass(4,ipri)
 msec = xyzmh_ptmass(4,isec)

 xpri(:) = xyzmh_ptmass(1:3,ipri)
 vpri(:) = vxyz_ptmass(1:3,ipri)
 xsec(:) = xyzmh_ptmass(1:3,isec)
 vsec(:) = vxyz_ptmass(1:3,isec)
 dr(:) = xpri(:) - xsec(:)
 dv(:) = vpri(:) - vsec(:)
 rbin  = sqrt(dot_product(dr,dr))

! Calculate the binary specific relative ang. mom and energy
 call cross(dr,dv,L)
 Lmag = sqrt(dot_product(L,L))
 E = 0.5*dot_product(dv,dv) - G*(mpri+msec)/rbin

 if (abs(E) < tiny(E)) stop 'binary energy problem'
 call get_ae(Lmag,E,mpri,msec,a,ecc)

 if (time <= tiny(time)) then
    open(iunit,file=output,status='replace',action='write',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open binary.dat file at t=0.0')
    write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
         1,'time', &
         2,'a', &
         3,'eccen'
 else
    inquire(file=output,exist=exists)
    if (.not. exists) call fatal(analysistype,'t /= 0.0, but the analysis output file does not exist...')
    open(iunit,file=output,status='old',action='write',position='append',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open binary.dat file during run')
 endif
 write(iunit,'(3(ES18.10,1X))') time,a,ecc
 close(iunit)
!-----------------------------------------------------------------------
end subroutine get_binary_params
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
subroutine get_ae(Lmag,E,m1,m2,a,ecc)
!-----------------------------------------------------------------------
! Return the semi-major axis and eccentricity between two objects
!-----------------------------------------------------------------------
 implicit none
 real,intent(out) :: a,ecc
 real,intent(in) :: Lmag,E,m1,m2

 if (Lmag < tiny(Lmag)) stop 'Lmag is zero in get_ae'
 if (abs(E) < tiny(E)) stop 'E is zero in get_ae'

! Hence obtain the binary eccentricity
 ecc = sqrt(1.0 + (2.0*E*Lmag**2)/((m1+m2)**2))

! and semi-major axis
 a = Lmag*Lmag/((m1+m2)*(1.0-ecc*ecc))
!-----------------------------------------------------------------------
end subroutine get_ae
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
subroutine cross(a,b,c)
!-----------------------------------------------------------------------
! Return the vector cross product of two 3d vectors
!-----------------------------------------------------------------------
 implicit none
 real,intent(in),dimension(3)  :: a,b
 real,intent(out),dimension(3) :: c

 c(1) = a(2)*b(3)-b(2)*a(3)
 c(2) = a(3)*b(1)-b(3)*a(1)
 c(3) = a(1)*b(2)-b(1)*a(2)

!-----------------------------------------------------------------------
end subroutine cross
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
! Function to calculate particle energy
!-----------------------------------------------------------------------
real(kind=8) function get_particle_energy(G,msink,ri,vi,ui)
 implicit none
 real,intent(in) :: G,msink,ri,ui
 real,dimension(3),intent(in) :: vi

 get_particle_energy = -G*msink/ri + 0.5*dot_product(vi,vi) + ui

!-----------------------------------------------------------------------
end function get_particle_energy
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
! Subroutine to return particle internal energy and sound speed
!-----------------------------------------------------------------------
subroutine get_utherm(ieos,xi,yi,zi,gamma,ui,csi)
 use eos, only:equationofstate
 implicit none

 integer,intent(in) :: ieos
 real,intent(in) :: xi,yi,zi,gamma
 real,intent(out) :: ui,csi

 real :: rhoi = 1.0 ! this is essentially a dummy variable, not needed here
!(only needed if adiabatic, but this routine is not called in that case...)
 real :: ponrhoi

 call equationofstate(ieos,ponrhoi,csi,rhoi,xi,yi,zi)

 if (gamma == 1.0) then
    ui = ponrhoi
 else
    ui = ponrhoi/(gamma-1.0)
 endif

!-----------------------------------------------------------------------
end subroutine get_utherm
!-----------------------------------------------------------------------
end module analysis

