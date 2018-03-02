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
!  Analysis routine calculating to determine the radial profile of a sphere
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, eos, infile_utils, part, physcon,
!    rho_profile
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'polytropes'
 public :: do_analysis
 !--Free parameters
 logical, private :: firstcall = .true.
 logical, private :: binary    = .false. ! The model is of a binary star
 real,    private :: rthresh   = 0.5     ! Radius within which the L1 error will be calculated
 integer, private :: frequency = 10      ! Will determine the density profile every frequency-th dump;
 ! keep this since, if binary, then the CoM will be required
 ! more frequently
 !--stored values to track orbital period of binary
 integer          :: iperiod(4)
 real             :: period(8,32)
 logical          :: nextP  (4)

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use centreofmass,   only: get_centreofmass,reset_centreofmass
 use rho_profile,    only: rho_polytrope
 use physcon,        only: pi
 use part,           only: rhoh
 use eos,            only: gamma,polyk
 use infile_utils,   only: open_db_from_file,inopts,close_db,read_inopt
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time
 integer,          parameter     :: ng   = 2048
 integer,          parameter     :: bins =  128
 integer,          parameter     :: lu   = 21
 integer                         :: i,j,k,npts,icount,nstar,is,ie,imin,idr1,nerr,ierr,nL1
 real                            :: r(ng),den(ng),profile(4,npart),profileH(4,bins)
 real                            :: totmass,rad,rho_actual,rho_exact,pd,L1error
 real                            :: deltar,rmin,rmax,aver(4),proftmp(4)
 real(kind=8)                    :: xposA(3),xposB(3),vpos(3),sep
 real                            :: xi,yi,zi
 logical                         :: iexist,fill_next,calc_average
 type(inopts), allocatable       :: db(:)
 character(len=200)              :: fileout
 !
 !--from .setup, determine if binary or not
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'.setup'
 inquire(file=fileout,exist=iexist)
 if ( iexist ) then
    write(*,'(2a)') "reading setup file: ",trim(fileout)
    call open_db_from_file(db,fileout,lu,ierr)
    call read_inopt(binary,'binary',db,errcount=nerr)
    write(*,*) "From ",trim(fileout),", binary = ",binary
    call close_db(db)
 endif
 !
 !
 !--Reset centre of mass
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 if ( binary ) then
    !
    !--Open file (appendif exists)
    fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_centres.dat'
    inquire(file=fileout,exist=iexist)
    if ( .not.iexist .or. firstcall ) then
       open(iunit,file=fileout,status='replace')
       write(iunit,"('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
             1,'time', &
             2,'x_A',  &
             3,'y_A',  &
             4,'z_A',  &
             5,'x_B',  &
             6,'y_B',  &
             7,'z_B',  &
             8,'d'
       period  = 0.0
       iperiod = 1
       nextP   = .false.
    else
       open(iunit,file=fileout,position='append')
    endif
    !
    !--Determine the centre of mass of each star
    call get_centreofmass(xposA,vpos,npart/2,xyzh(:,1:npart/2),vxyzu(:,1:npart/2))
    call get_centreofmass(xposB,vpos,npart/2,xyzh(:,npart/2+1:npart),vxyzu(:,npart/2+1:npart))
    sep = sqrt(dot_product(xposA-xposB,xposA-xposB))
    write(iunit,'(8(es18.10,1x))') time,xposA,xposB,sep
    close(iunit)
    nstar     = 2
    !--to track period
    call track_period(time,xposA(1),1)
    call track_period(time,xposA(2),2)
    call track_period(time,xposB(1),3)
    call track_period(time,xposB(2),4)
 else
    nstar     = 1
    frequency = 1
    xposA     = 0.0
 endif
 !
 !
 L1error = 0.0
 nL1     = 0
 if (mod(num,frequency)==0) then
    !
    !--Calculate the exact polytropic solution
    totmass = npart*particlemass/nstar
    print*, "polyk, gamma = ", polyk,gamma
    print*, "total mass = ", totmass
    call rho_polytrope(gamma,polyk,totmass,r,den,npts)
    !
    !--For each particle determine its actual density, and determine what its density should be
    do i = 1,npart
       if (i > npart/nstar) then
          xi = xyzh(1,i) - xposB(1)
          yi = xyzh(2,i) - xposB(2)
          zi = xyzh(3,i) - xposB(3)
       else
          xi = xyzh(1,i) - xposA(1)
          yi = xyzh(2,i) - xposA(2)
          zi = xyzh(3,i) - xposA(3)
       endif
       rad        = sqrt( xi*xi + yi*yi + zi*zi )
       rho_actual = rhoh(xyzh(4,i),particlemass)
       j = 2
       do while(j < npts .and. rad > r(j))
          j = j + 1
       enddo
       rho_exact = den(j-1) +  (den(j)-den(j-1))/(r(j)-r(j-1)) * (rad-r(j-1))
       pd = abs(rho_actual-rho_exact)/rho_exact * 100.0
       !
       profile(1,i) = rad         ! radial position
       profile(2,i) = rho_actual  ! actual density
       profile(3,i) = rho_exact   ! exact density
       profile(4,i) = pd          ! percent difference
       if (rad < rthresh) then
          L1error = L1error + abs(rho_actual-rho_exact)   ! L1 error
          nL1     = nL1     + 1
       endif
    enddo
    !
    !--Brute force sort of profile by profile(1,:)
    !
    do k = 1,nstar
       if (k==1) then
          is = 1
          ie = npart/nstar
       else
          is = npart/nstar+1
          ie = npart
       endif
       do i = is,ie
          rmin = profile(1,i)
          imin = i
          do j = i+1,ie
             if (profile(1,j) < rmin) then
                rmin = profile(1,j)
                imin = j
             endif
          enddo
          proftmp         = profile(:,i   )
          profile(:,i   ) = profile(:,imin)
          profile(:,imin) = proftmp
       enddo
    enddo
    !
    !--Average values over deltar
    !
    profileH = 0.0
    j        = 0
    deltar   = 1.0/float(bins)*nstar
    idr1     = npart
    do k = 1,nstar
       aver   = 0.0
       rmin   = 0.0
       rmax   = deltar
       icount = 0
       calc_average = .false.
       if (k==1) then
          is = 1
          ie = npart/nstar
       else
          is = npart/nstar+1
          ie = npart
       endif
       do i = is,ie
          if (profile(1,i) < rmax) then
             aver   = aver + profile(:,i)
             icount = icount + 1
          else
             calc_average = .true.
          endif
          if (i==ie) calc_average = .true.
          if (calc_average) then
             ! calculated the averaged values
             j = j + 1
             calc_average    = .false.
             profileH(  1,j) = 0.5*(rmax+rmin)
             if (icount > 0)    profileH(2:3,j) = aver(2:3)/float(icount)
             if (aver(3) > 0.0) profileH(  4,j) = abs(aver(2)-aver(3))/aver(3)*100.0
             ! rest values with new properties
             fill_next = .true.
             rmin = rmax
             do while (fill_next)
                rmax = rmax + deltar
                aver = profile(:,i)
                icount = 1
                if (profile(1,i) < rmax) fill_next = .false.
             enddo
             if (k==1 .and. i==ie) idr1 = j
          endif
       enddo
    enddo
    !
    !--Write results to file
    write(fileout,'(3a)') 'analysisout_',trim(dumpfile),'.dat'
    fileout=trim(fileout)
    open(iunit,file=fileout)
    write(iunit,"('#',4(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'r',       &
          2,'rho_act', &
          3,'rho_exa', &
          4,'p.diff'
    do i = 1,npart
       write(iunit,'(4(1pe18.10,1x))') profile(:,i)
       if (i==npart/nstar .and. nstar==2) write(iunit,'(a)') " "
    enddo
    write(iunit,'(a)') " "
    do i = 1,j
       write(iunit,'(4(1pe18.10,1x))') profileH(:,i)
       if (i==idr1) write(iunit,'(a)') " "
    enddo
    close(iunit)
    !
    !--print period tracking to file (overwriting anything in existance)
    if ( binary ) then
       fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'_period.dat'
       fileout=trim(fileout)
       open(iunit,file=fileout)
       write(iunit,"('#',4(1x,'[',i2.2,1x,a11,']',2x))") &
             1,'P_xposA',       &
             2,'P_yposA', &
             3,'P_xposB', &
             4,'P_yposB'
       do i = 2,min(iperiod(1),iperiod(2),iperiod(3),iperiod(4))-1
          write(iunit,'(4(1pe18.10,1x))') period(5:8,i)-period(5:8,i-1)
       enddo
       close(iunit)
    endif
    L1error = L1error/float(nL1)
    write(*,*) "npart,N(Rthresh), L_1 error = ",npart,nL1,L1error
 endif
 !
 if (firstcall) firstcall = .false. ! performed here since multiple commands require this knowledge
 !
end subroutine do_analysis
!
subroutine track_period(time,pos,ival)
 implicit none
 real,         intent(in) :: time
 real(kind=8), intent(in) :: pos
 integer,      intent(in) :: ival
 !
 if (pos < period(ival,iperiod(ival))) then
    period(ival  ,iperiod(ival)) = pos
    period(ival+4,iperiod(ival)) = time
    nextP  (ival) = .true.
 else if (pos > 0.0 .and. nextP  (ival)) then
    iperiod(ival) = iperiod(ival) + 1
    nextP  (ival) = .false.
 endif
 !
end subroutine track_period
!-----------------------------------------------------------------------
!
end module
