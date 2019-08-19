!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION: analysis programmes used to analys a neutron star merger.
!  This is an interactive routine that includes multiple analysis options
!  Note: all outputs are in code units unless explicitly stated
!  Author: Bernard Field & Madeline Marshall (supervisors: James Wurster & Paul Lasky)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, extern_gwinspiral, io, part, physcon,
!    prompting, readwrite_dumps, units
!+
!--------------------------------------------------------------------------
module analysis
 use io,              only: fatal
 use part,            only: rhoh
 use physcon,         only: pi
 use centreofmass,    only: get_centreofmass
 use readwrite_dumps, only: opened_full_dump
 use extern_gwinspiral, only:Nstar
 implicit none
 character(len=20), parameter, public :: analysistype = 'NSmerger'
 !
 integer, parameter, private :: nana_opts          = 20
 real,               private :: density_cutoff_cgs = 5.0d-5     ! The density threshhold in code units (opts 2,3,4)
 real,               private :: thickness          = 2.0
 real,               private :: dtheta             = 5.*pi/180. ! width of slice in the directions of the minor and major axes (opt 4)
 logical,            private :: firstcall          = .true.
 !
 integer,            private :: choice
 real,               private :: density_cutoff
 real,               private :: com(3),vcom(3),com1(3),com2(3),vcom1(3),vcom2(3)
 logical,            private :: iexist
 character(len=200), private :: fileout,analysis_opt(nana_opts)
 !
 private :: trace_com,calculate_TW,calculate_I,calculate_midplane_profile
 private :: get_momentofinertia,jacobi
 public  :: do_analysis

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use prompting,      only: prompt
 use units,          only: unit_density
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i
 !
 !--Determine which analysis to run if this is the first step
 !
 if ( firstcall ) then
    analysis_opt(:) = 'none'
    analysis_opt(1) = 'Trace centre of mass of each star and the system'
    analysis_opt(2) = 'Brute force calculation of T/W'
    analysis_opt(3) = 'Calculate the moment of inertia tensor'
    analysis_opt(4) = 'Calculate the radial profile of a slice through the midplane'
    write(*,"(a)") 'Analysis options: '
    do i = 1, nana_opts
       if (trim(analysis_opt(i)) /= 'none') write(*,"(a5,i2,1x,a60)") 'Case ', i, analysis_opt(i)
    enddo
    choice = 1
    call prompt('Enter analysis choice',choice,1,nana_opts)
    !
    !--Get the index range for each star
    ! (Note from DJP: This should now be automatically read from the dump header)
    if (Nstar(1) <= 0) call fatal('analysis_NSmerger','Require Nstar(1) > 0 in header of dump file')
    if (Nstar(2) <= 0) call fatal('analysis_NSmerger','Require Nstar(2) > 0 in header of dump file')
    !
    !--Prompt for the density_cut off
    if (choice > 1) then
       call prompt('Enter cutoff density (cgs):',density_cutoff_cgs,0.)
       density_cutoff = density_cutoff_cgs / unit_density
    endif
    !--Prompt for thickness
    if (choice == 4) then
       call prompt('Enter the thickness of the mid-plane slice (code units):',thickness,0.)
    endif
 endif
 !
 !--Calculate the centre of mass and velocity of the system
 !
 call get_centreofmass(com,vcom,npart,xyzh,vxyzu)
 !
 !--Run the analysis
 !
 select case(choice)
 case(1)
    !--Trace the centre of masses
    call trace_com(dumpfile,xyzh,vxyzu,time,npart,iunit)
 case(2)
    !--Calculate T/W using brute force
    call calculate_TW(dumpfile,xyzh,vxyzu,time,npart,iunit,particlemass)
 case(3)
    !--Calculate the moment of inertia tensor
    call calculate_I(dumpfile,xyzh,time,npart,iunit,particlemass)
 case(4)
    !--Calculate the radial profile of a slice through the midplane'
    call calculate_midplane_profile(dumpfile,xyzh,vxyzu,npart,iunit,particlemass)
 end select
 !
 close(iunit)
 firstcall = .false. ! done here since this logical is required for opening files
end subroutine do_analysis
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Trace centre of mass of each star and the system
!+
!-----------------------------------------------------------------------
subroutine trace_com(dumpfile,xyzh,vxyzu,time,npart,iunit)
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: time
 real                         :: rad
 !
 !--Open file (appendif exists)
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_orbit.dat'
 inquire(file=fileout,exist=iexist)
 if ( firstcall .or. .not.iexist ) then
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',11(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'time',&
        2,'x',   &
        3,'y',   &
        4,'z',   &
        5,'x1',  &
        6,'y1',  &
        7,'z1',  &
        8,'x2',  &
        9,'y2',  &
       10,'z2',  &
       11,'r'
 else
    open(iunit,file=fileout,position='append')
 endif
 !
 !--Get centre of masses of the stars
 call get_centreofmass(com1,vcom1,nstar(1),xyzh(:,1:nstar(1)),vxyzu(:,1:nstar(1)))
 call get_centreofmass(com2,vcom2,nstar(2),xyzh(:,nstar(1)+1:npart),vxyzu(:,nstar(1)+1:npart))
 !
 rad = sqrt( (com1(1)-com2(1))**2 + (com1(2)-com2(2))**2 + (com1(3)-com2(3))**2 )
 !
 !--Write results to file
 write(iunit,'(11(Es18.10,1x))') time,com,com1,com2,rad
 !
end subroutine trace_com
!-----------------------------------------------------------------------
!+
! Calculate T/W using brute force
!+
!-----------------------------------------------------------------------
subroutine calculate_TW(dumpfile,xyzh,vxyzu,time,npart,iunit,particlemass)
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i,j,npartmeasured
 real                         :: rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2,rad2,rmax2
 real                         :: erot,erotx,eroty,erotz,grav
 real                         :: r(3),v(3)
 !
 !--Skip if not a full dump
 if (.not.opened_full_dump) return
 !
 !--Open file (appendif exists)
 fileout = trim(dumpfile(1:index(dumpfile,'_')-1))//'_TW.dat'
 inquire(file=fileout,exist=iexist)
 if ( firstcall .or. .not.iexist ) then
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',       &
          2,'T/W',     &
          3,'T',     &
          4,'W',    &
          5,'Star mass', &
          6,'Star radius'
 else
    open(iunit,file=fileout,position='append')
 endif
 !
 !--Calculate rotational kinetic energy
 npartmeasured = 0
 rmax2 = 0.
 erot  = 0.
 erotx = 0.
 eroty = 0.
 erotz = 0.
!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,particlemass,com,vcom,density_cutoff) &
!$omp private(i,r,v,rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2,rad2) &
!$omp reduction(+:erotx,eroty,erotz,npartmeasured) &
!$omp reduction(max:rmax2)
!$omp do
 do i=1,npart
    if (rhoh(xyzh(4,i),particlemass) > density_cutoff) then
       r = xyzh(1:3,i)  - com
       v = vxyzu(1:3,i) - vcom
       ! r cross v
       rcrossvx = (r(2)*v(3) - r(3)*v(2))
       rcrossvy = (r(3)*v(1) - r(1)*v(3))
       rcrossvz = (r(1)*v(2) - r(2)*v(1))
       ! rotational energy around each axis through the origin
       radxy2 = r(1)*r(1) + r(2)*r(2)
       radyz2 = r(3)*r(3) + r(2)*r(2)
       radxz2 = r(1)*r(1) + r(3)*r(3)
       if (radyz2 > 0.) erotx = erotx + particlemass*rcrossvx*rcrossvx/radyz2
       if (radxz2 > 0.) eroty = eroty + particlemass*rcrossvy*rcrossvy/radxz2
       if (radxy2 > 0.) erotz = erotz + particlemass*rcrossvz*rcrossvz/radxy2
       !
       ! size of the star
       rad2 = dot_product(r,r)
       rmax2 = max(rmax2,rad2)
       !
       ! additional bookkeeping
       npartmeasured = npartmeasured + 1
    endif
 enddo
!$omp enddo
!$omp end parallel
 erotx = 0.5*erotx
 eroty = 0.5*eroty
 erotz = 0.5*erotz
 erot  = sqrt(erotx**2 + eroty**2 + erotz**2)
 !
 !--Calculate gravitational potential energy
 grav = 0.
!$omp parallel default(none) &
!$omp shared(npart,xyzh,particlemass,density_cutoff) &
!$omp private(i,j,r,rad2) &
!$omp reduction(+:grav)
!$omp do
 do i=1,npart
    if (rhoh(xyzh(4,i),particlemass) > density_cutoff) then
       do j=i+1,npart
          if (rhoh(xyzh(4,j),particlemass) > density_cutoff) then
             r    = xyzh(1:3,i) - xyzh(1:3,j)
             rad2 = dot_product(r,r)
             if (rad2 > 0.0) grav = grav + 1.0/sqrt(rad2)
          endif
       enddo
    endif
 enddo
!$omp enddo
!$omp end parallel
 grav = -grav*particlemass*particlemass
 !
 !--Write results
 print *, "time:",time," T/|W|:", erot/abs(grav), " T:", erot, " W:", grav, "ignored particles:", npart - npartmeasured, &
  "star mass:", npartmeasured*particlemass, "rstar:",sqrt(rmax2)
 write(iunit,'(6(es18.10,1x))') time,erot/abs(grav),erot,grav,npartmeasured*particlemass,sqrt(rmax2)
 !
end subroutine calculate_TW
!-----------------------------------------------------------------------
!+
!  Determines the moment of inertia tensor, in diagonalised form.
!  Can exclude particles which are below a specified cut-off density.
!  Will output the mass and radius of the measured area.
!+
!-----------------------------------------------------------------------
subroutine calculate_I(dumpfile,xyzh,time,npart,iunit,particlemass)
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: npart,iunit
 real,             intent(in) :: xyzh(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: i,npartused
 real                         :: rmax,bigI,medI,smallI
 real                         :: principle(3),evectors(3,3),ellipticity(2)
 !
 !--Open file (appendif exists)
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_inertia.dat'
 inquire(file=fileout,exist=iexist)
 if ( firstcall .or. .not.iexist ) then
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',18(1x,'[',i2.2,1x,a11,']',2x))") &
           1,'time',&
           2,'I1',  &
           3,'I2',  &
           4,'I3',  &
           5,'e1',  &
           6,'e2',  &
           7,'v1,1',&
           8,'v1,2',&
           9,'v1,3',&
          10,'v2,1',&
          11,'v2,2',&
          12,'v2,3',&
          13,'v3,1',&
          14,'v3,2',&
          15,'v3,3',&
          16,'excluded parts',&
          17,'mstar',     &
          18,'rstar'
 else
    open(iunit,file=fileout,position='append')
 endif
 !
 !--Calculate the tensor
 call get_momentofinertia(xyzh,npart,npartused,principle,evectors,particlemass,rmax)
 !
 !--Sort the principle moments, since ellipticity depends on it.
 bigI   = maxval(principle)
 smallI = minval(principle)
 medI   = 0.5*(bigI+smallI)  ! to avoid compiler warning
 do i=1,3
    if (smallI < principle(i) .and. principle(i) < bigI) medI = principle(i)
 enddo
 ellipticity(1) = sqrt(2.0*(bigI-smallI)/smallI)
 ellipticity(2) = sqrt(2.0*(bigI-medI ) /medI  )
 !
 !--Write to file
 write(iunit,'(18(es18.10,1x))') &
     time,principle(1),principle(2),principle(3),&
     ellipticity(1),ellipticity(2),&
     evectors(1,1),evectors(2,1),evectors(3,1),&
     evectors(1,2),evectors(2,2),evectors(3,2),&
     evectors(1,3),evectors(2,3),evectors(3,3),&
     real(npart-npartused),npartused*particlemass,rmax
 !
end subroutine calculate_I
!-----------------------------------------------------------------------
!+
!  Determines the radial profile of the midplane (of a given thickness)
!+
!-----------------------------------------------------------------------
subroutine calculate_midplane_profile(dumpfile,xyzh,vxyzu,npart,iunit,particlemass)
 use part, only: alphaind
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass
 integer,          parameter  :: nbins = 800
 integer                      :: i,j,npartused,zloc,majorloc,minorloc
 real                         :: rmax,thetamajor,thetaminor
 integer                      :: bincountmaj(nbins),bincountmin(nbins),bincountavg(nbins)
 real                         :: rtocm(npart),theta(npart),vtheta(npart),angv(npart)
 real                         :: avvinbinmaj(nbins),avvinbinmin(nbins),avvinbinavg(nbins)
 real                         :: vinbinmaj(nbins),  vinbinmin(nbins),  vinbinavg(nbins)
 real                         :: alphabinmaj(nbins),alphabinmin(nbins),alphabinavg(nbins)
 real                         :: partdensmaj(nbins),partdensmin(nbins),partdensavg(nbins)
 real                         :: radbin(nbins),vol(nbins)
 real                         :: principle(3),principlenew(2),evectors(3,3),evectorsnew(3,2)
 real                         :: major(3),minor(3)
 !
 !--Skip if not a full dump
 if (.not.opened_full_dump)return
 !
 !--Open file
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_rotataxesprofile'//trim(dumpfile(INDEX(dumpfile,'_'):))//'.dat'
 open(iunit,file=fileout,status='replace')
 write(iunit,"('#',10(1x,'[',i2.2,1x,a11,']',2x))") &
    1,'outer bin rad',&
    2,'density_min',  &
    3,'ang v_minor',  &
    4,'alpha_minor',  &
    5,'density_maj',  &
    6,'ang v_major',  &
    7,'alpha_major',  &
    8,'density_tot',  &
    9,'ang v_tot',    &
   10,'alpha_tot'
 !
 !--Initialise variables
 !
 theta       = 0.0
 rtocm       = 0.0
 vtheta      = 0.0
 angv        = 0.0
 vol         = 0.0
 bincountmaj = 0
 vinbinmaj   = 0.0
 avvinbinmaj = 0.0
 alphabinmaj = 0.0
 partdensmaj = 0.0
 bincountmin = 0
 vinbinmin   = 0.0
 avvinbinmin = 0.0
 alphabinmin = 0.0
 partdensmin = 0.0
 bincountavg = 0
 vinbinavg   = 0.0
 avvinbinavg = 0.0
 alphabinavg = 0.0
 partdensavg = 0.0
 !
 !--Calculate radius and angle of particle from CoM coordinates in x-y plane
!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,com,vcom,rtocm,theta,vtheta,angv) &
!$omp private(i)
!$omp do
 do i=1,npart
    rtocm (i) = sqrt((xyzh(1,i)-com(1))**2 + (xyzh(2,i)-com(2))**2)
    theta (i) = atan2((xyzh(2,i)-com(2)),(xyzh(1,i)-com(1)))
    vtheta(i) = -(vxyzu(1,i)-vcom(1))*sin(theta(i)) + (vxyzu(2,i)-vcom(2))*cos(theta(i))
    angv  (i) = abs(vtheta(i))/rtocm(i)
 enddo
!$omp enddo
!$omp end parallel
 !
 !--Calculate moment of inertia
 call get_momentofinertia(xyzh,npart,npartused,principle,evectors,particlemass,rmax)
 !
 !--Find location of major and minor axes
 zloc = maxloc(evectors(3,:),1)
 j = 1
 do i = 1,3
    if (i/=zloc) then
       principlenew( j) =principle( i)
       evectorsnew(:,j) =evectors(:,i)
       j = j+1
    endif
 enddo
 majorloc   = minloc(principlenew,1)
 minorloc   = maxloc(principlenew,1)
 major      = evectorsnew(:,majorloc)
 minor      = evectorsnew(:,minorloc)
 thetamajor = atan2(major(2),major(1))
 thetaminor = atan2(minor(2),minor(1))
 write(*,*) 'evectors',evectors
 write(*,*) 'evectorsnew',evectorsnew
 write(*,*) 'major',major
 write(*,*) 'minor',minor
 write(*,*) 'thetamajor',thetamajor
 write(*,*) 'thetaminor',thetaminor
 !
 !Set radii and calculate volume of slice bins
 rmax = maxval(rtocm)
 do i = 1,nbins
    radbin(i) = rmax*float(i)/float(nbins)
    if (i==1) then
       vol(i) = thickness*dtheta*radbin(1)**2
    else
       vol(i) = thickness*dtheta*(radbin(i)**2-radbin(i-1)**2)
    endif
 enddo
 !
 !--Sort particles into bins, find total angular velocity and total alpha in each
 do i=1,npart
    do j=1,nbins
       if (rtocm(i)<radbin(j)) then
          if (abs(xyzh(3,i)-com(3)) < 0.5*thickness) then
             bincountavg(j) = bincountavg(j) + 1
             vinbinavg  (j) = vinbinavg  (j) + angv(i)
             alphabinavg(j) = alphabinavg(j) + alphaind(1,i)
             if (theta(i) < (thetamajor+dtheta) .and. theta(i) > (thetamajor-dtheta)) then
                bincountmaj(j) = bincountmaj(j) + 1
                vinbinmaj  (j) = vinbinmaj  (j) + angv(i)
                alphabinmaj(j) = alphabinmaj(j) + alphaind(1,i)
             elseif (theta(i) < (thetaminor+dtheta) .and. theta(i) > (thetaminor-dtheta)) then
                bincountmin(j) = bincountmin(j) + 1
                vinbinmin  (j) = vinbinmin  (j) + angv(i)
                alphabinmin(j) = alphabinmin(j) + alphaind(1,i)
             endif
          endif
          exit
       endif
    enddo
 enddo
 !
 !--Convert totals to averages for each bin
 do i = 1,nbins
    if (bincountmaj(i) > 0) then
       avvinbinmaj(i) = vinbinmaj(i)  /float(bincountmaj(i))
       alphabinmaj(i) = alphabinmaj(i)/float(bincountmaj(i))
       partdensmaj(i) = float(bincountmaj(i))*particlemass/vol(i)
    endif
    if (bincountmin(i) > 0) then
       avvinbinmin(i) = vinbinmin(i)  /float(bincountmin(i))
       alphabinmin(i) = alphabinmin(i)/float(bincountmin(i))
       partdensmin(i) = float(bincountmin(i))*particlemass/vol(i)
       print*, partdensmin(i) ,float(bincountmin(i)),particlemass,vol(i)
    endif
    if (bincountavg(i) > 0) then
       avvinbinavg(i) = vinbinavg(i)  /float(bincountavg(i))
       alphabinavg(i) = alphabinavg(i)/float(bincountavg(i))
       partdensavg(i) = float(bincountavg(i))*particlemass/(vol(i)*pi/dtheta)
    endif
 enddo
 !
 !--Write results to file
 do i=1,nbins
    write(iunit,'(10(es18.10,1x))')  radbin(i), partdensmin(i), avvinbinmin(i), alphabinmin(i), &
    partdensmaj(i), avvinbinmaj(i), alphabinmaj(i), partdensavg(i), avvinbinavg(i), alphabinavg(i)
 enddo
!
end subroutine calculate_midplane_profile
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Calculates the moment of inertia
! This is done about the coordinate axes whose origin is at the
! centre of mass
!+
!-----------------------------------------------------------------------
subroutine get_momentofinertia(xyzh,npart,npartused,principle,evectors,particlemass,rmax)
 integer,          intent(in)  :: npart
 integer,          intent(out) :: npartused
 real,             intent(in)  :: xyzh(:,:)
 real,             intent(in)  :: particlemass
 real,             intent(out) :: principle(3), evectors(3,3),rmax
 integer                       :: i
 real                          :: inertia(3,3)
 real                          :: x,y,z,r2,rmax2
 !
 inertia   = 0.
 npartused = 0
 rmax2     = 0.0
 do i = 1,npart
    if (rhoh(xyzh(4,i),particlemass) > density_cutoff) then
       x = xyzh(1,i) - com(1)
       y = xyzh(2,i) - com(2)
       z = xyzh(3,i) - com(3)
       inertia(1,1) = inertia(1,1) + y**2 + z**2
       inertia(2,2) = inertia(2,2) + x**2 + z**2
       inertia(3,3) = inertia(3,3) + x**2 + y**2
       inertia(1,2) = inertia(1,2) - x*y
       inertia(1,3) = inertia(1,3) - x*z
       inertia(2,3) = inertia(2,3) - y*z
       ! Additional useful values
       npartused    = npartused + 1
       r2           = x*x + y*y + z*z
       rmax2        = max(rmax2,r2)
    endif
 enddo
 rmax = sqrt(rmax2)
 !--The symmetric components
 inertia(2,1) = inertia(1,2)
 inertia(3,1) = inertia(1,3)
 inertia(3,2) = inertia(2,3)
 !--Multiply in constant
 inertia      = inertia*particlemass
 !
 !--Find the eigenvectors
 !  note: i is a dummy out-integer that we don't care about
 call jacobi(inertia,3,3,principle,evectors,i)
 !
end subroutine get_momentofinertia
!-----------------------------------------------------------------------
!+
! Calculates the Jacobian
! Source: http://www.fing.edu.uy/if/cursos/fiscomp/extras/numrec/book/f11.pdf
!+
!-----------------------------------------------------------------------
subroutine jacobi(a,n,np,d,v,nrot)
 integer, intent(in)    :: n,np
 integer, intent(out)   :: nrot
 real,    intent(inout) :: a(np,np)
 real,    intent(out)   :: d(np),v(np,np)
 integer, parameter :: nmax = 500
!
! Computes all eigenvalues and eigenvectors of a real symmetric matrix, a,
! whichisofsize n by n, stored in a physical np by np array.
! On output, elements of a above the diagonal are destroyed.
! d returns the eigenvalues of a in its first n elements.
! v is a matrix with the same logical  and  physical  dimensions  as a,
! whose  columns  contain,  on  output,  the  normalized eigenvectors of a.
! nrot returns the number  of Jacobi rotations that were required.
!
 integer :: i,ip,iq,j
 real ::  c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
 do 12, ip=1,n  !Initialize  to  the  identity  matrix.
    do 11, iq=1,n
       v(ip,iq)=0.
11  enddo
    v(ip,ip)=1.
12 enddo
 do 13, ip=1,n
    b(ip)=a(ip,ip)
!Initialize b and d to the diagonal of a.
    d(ip)=b(ip)
    z(ip)=0.
!This  vector  will  accumulate  terms  of  the  form tapq as  in equation  (11.1.14).
13 enddo

 nrot=0
 do 24,i=1,50
    sm=0.
    do 15,ip=1,n-1
!Sum  off-diagonal elements.
       do 14,iq=ip+1,n
          sm=sm+abs(a(ip,iq))
14     enddo
15  enddo
    if (sm==0.)return
!The normal return, which relies on quadratic convergence to machine  underflow.
    if (i < 4) then
       tresh=0.2*sm/n**2
!...on the first  three sweeps.
    else
       tresh=0.
!...thereafter.
    endif
    do 22,ip=1,n-1
       do 21,iq=ip+1,n
          g=100.*abs(a(ip,iq))
!After four sweeps, skip the rotation if the off-diagonal element is small.
          if ((i > 4).and.(abs(d(ip))+g==abs(d(ip))).and.(abs(d(iq))+g==abs(d(iq)))) then
             a(ip,iq)=0.
          elseif (abs(a(ip,iq)) > tresh) then
             h=d(iq)-d(ip)
             if (abs(h)+g==abs(h)) then
                t=a(ip,iq)/h
!t=1/(2(theta))
             else
                theta=0.5*h/a(ip,iq)
!Equation  (11.1.10).
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if (theta < 0.)t=-t
             endif
             c=1./sqrt(1+t**2)
             s=t*c
             tau=s/(1.+c)
             h=t*a(ip,iq)
             z(ip)=z(ip)-h
             z(iq)=z(iq)+h
             d(ip)=d(ip)-h
             d(iq)=d(iq)+h
             a(ip,iq)=0.
             do 16,j=1,ip-1
!Case of rotations 1<=j<p.
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16           enddo
             do 17,j=ip+1,iq-1
!Case of rotations p<j<q.
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17           enddo
             do 18,j=iq+1,n
!Case of rotations q<j<=n.
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18           enddo
             do 19,j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19           enddo
             nrot=nrot+1
          endif
21     enddo
22  enddo
    do 23,ip=1,n
       b(ip)=b(ip)+z(ip)
       d(ip)=b(ip)
!Update d with the  sum of tapq,
       z(ip)=0.
!and  reinitialize z.
23  enddo
24 enddo
 return
end subroutine jacobi

end module
