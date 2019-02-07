!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: lombperiod
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: lombperiod time1 time2 *-colxxxxx.mflow file
!
!  DEPENDENCIES: powerspectrums
!+
!--------------------------------------------------------------------------
program lombperiod

 use powerspectrums, only: powerspectrum

 implicit none

 integer           :: nargs,iunitin=42,iunitout=43,ierr=0,i=0,k=1,nfreq,imf
 character(len=120):: mflowfile,outname,num1,num2
 real              :: time1,time2
 real, allocatable :: time(:),dat(:)
 !--physical variables
 real              :: tbin,binarysep,maxmultiple=8,maxfreq
 real, allocatable :: freq(:), power(:)

 binarysep=1 ! change to get appropriate time normalization
 nfreq=5000

 tbin=1!2*pi*binarysep**(1.5)
 maxfreq=maxmultiple*tbin

 allocate(freq(nfreq),power(nfreq))
 power(:)=0
 freq(:)=0

 !--frequencies
 do i=1, nfreq
    freq(i)=maxfreq*i/nfreq
! print*,freq(i),i, maxfreq, maxmultiple,nfreq
 enddo

 i=1

 nargs = command_argument_count()
 if (nargs < 3 .or. nargs>3 ) then
    print "(a,/)",' Lomb-Scargle periodogram'
    call get_command_argument(0,mflowfile)
    print "(a)",' Usage: '//trim(mflowfile)//' time1 time2 *-colxxxxx.mflow file'
    stop
 endif

 call get_command_argument(1,num1)
 call get_command_argument(2,num2)
 call get_command_argument(3,mflowfile)

 !--Creating output name

 read(num1,*) time1
 read(num2,*) time2

 imf = index(mflowfile,'-col') + 1
 if (imf <= 0) then
    outname = trim(mflowfile)//'.lsp'//trim(num1)//'_'//trim(num2)
 else
    outname = trim(mflowfile(imf:len_trim(mflowfile)))//'.lsp'//trim(num1)//'_'//trim(num2)
 endif


 print*,outname
 !--Open input and output

 open(unit=iunitin,file=trim(mflowfile),status='old',form='formatted',iostat=ierr)
 open(unit=iunitout,file=trim(outname),status='replace',form='formatted')

 !--Check if inputfile present
 if(ierr/=0) then
    print*, "Couldn't find input file!"
    stop
 endif

 !--Count how many lines

 do while(ierr==0)
    read (iunitin,*,iostat=ierr)
    i=i+1
 enddo

 ierr=0
 close(iunitin)

 !--Allocating vectors
 allocate(time(i),dat(i))

 time(:)=0
 dat(:)=0

 open(unit=443,file=trim(mflowfile),status='old',form='formatted',iostat=ierr)
 !--Cycle over times
 do while(ierr==0)

    read(443,*,iostat=ierr)time(k),dat(k)
    if(time1>time(k)) then
       cycle
    elseif(time1<time(k) .and. time2>time(k)) then
       k=k+1
    endif
 enddo

 k=k-1 !otherwise count also last term in input file

! print*,freq, power
 time=time*tbin

 call powerspectrum(k,time,dat,nfreq,freq,power,.true.)


 !--write output

 do i=1, nfreq
    write(iunitout,"(2(1x,es18.10,1x))")freq(i)/(tbin),power(i)/k

 enddo


 print*, "i k :",i,k
 print*, mflowfile
 print*, outname
 print*, time1
 print*, time2
! print*, time
! print*, dat

 close(443)
 close(iunitout)

end program lombperiod



