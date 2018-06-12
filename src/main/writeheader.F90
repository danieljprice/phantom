!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: writeheader
!
!  DESCRIPTION: writes runtime header to the logfile
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dust, eos, growth, io, kernel, options,
!    part, physcon, readwrite_infile, units, viscosity
!+
!--------------------------------------------------------------------------
module writeheader
 implicit none
 public :: write_header,write_codeinfo

 private

contains

subroutine write_codeinfo(iunit)
 integer, intent(in) :: iunit
!
!--write out code name, version and time
!
 write(iunit,10) '1.0, released 13th March 2018'

10 format(/, &
   "  _ \  |                 |                    ___|   _ \  |   |",/, &
   " |   | __ \   _` | __ \  __|  _ \  __ `__ \ \___ \  |   | |   |",/, &
   " ___/  | | | (   | |   | |   (   | |   |   |      | ___/  ___ |",/, &
   "_|    _| |_|\__,_|_|  _|\__|\___/ _|  _|  _|_____/ _|    _|  _|",/, &
   "              ___                         , _               ",/, &
   " |)          (|  \  _,        o  _ |\    /|/ \ ,_  o  _   _ ",/, &
   " |/\_|  |     |   |/ |  /|/|  | |/ |/     |__//  | | /   |/ ",/, &
   "  \/  \/|/   (\__/ \/|_/ | |_/|/|_/|_/    |      |/|/\__/|_/",/, &
   "       (|                                                   ",//,  &
   " Version: ",a)

 return
end subroutine write_codeinfo

!-----------------------------------------------------------------
!+
subroutine write_header(icall,infile,evfile,logfile,dumpfile,ntot)
!
!  This subroutine writes header to logfile / screen
!  with the main parameters used for the run
!
!  icall = 1 (before particle have been set up)
!  icall = 2 (after particle setup)
!+
!-----------------------------------------------------------------
 use dim,              only:maxp,maxvxyzu,maxalpha,ndivcurlv,mhd_nonideal,nalpha,use_dustgrowth
 use io,               only:iprint
 use boundary,         only:xmin,xmax,ymin,ymax,zmin,zmax
 use options,          only:tolh,alpha,alphau,alphaB,ieos,alphamax,use_dustfrac
 use part,             only:hfact,massoftype,mhd,maxBevol,&
                            gravity,h2chemistry,periodic,npartoftype,massoftype,&
                            igas,idust,iboundary,istar,idarkmatter,ibulge
 use eos,              only:eosinfo
 use readwrite_infile, only:write_infile
 use physcon,          only:pi
 use kernel,           only:kernelname,radkern
 use viscosity,        only:irealvisc,viscinfo
 use units,            only:print_units
#ifdef DUST
 use dust,             only:print_dustinfo
#ifdef DUSTGROWTH
 use growth,                   only:print_growthinfo
#endif
#endif
 integer                      :: Nneigh,i
 integer,          intent(in) :: icall
 character(len=*), intent(in) :: infile,evfile,logfile,dumpfile
 integer(kind=8),  intent(in), optional :: ntot
 character(len=10) :: startdate, starttime
 character(len=11) :: parttype
! real :: have,hmin,hmax,v2i,B2i,pri,ponrhoi,spsoundi,rhoi

!-----------------------------------------------------------------------
! 1st header after options have been read, but before particle setup
!-----------------------------------------------------------------------

 if (icall==1) then

    call date_and_time(startdate,starttime)
    startdate = startdate(7:8)//'/'//startdate(5:6)//'/'//startdate(1:4)
    starttime = starttime(1:2)//':'//starttime(3:4)//':'//starttime(5:)
    write(iprint,"(' Run started on ',a,' at ',a)") startdate,starttime

    write(iprint, 20) trim(infile),trim(evfile),trim(logfile)
    if (iprint /= 6) write(*, 20) trim(infile),trim(evfile),trim(logfile)
20  format(/,                                &
      ' Read input from   : ',a,/,        &
      ' Writing energy to : ',a,/,        &
      ' Writing log to    : ',a,/)
!
!--write the input file into the run log
!
    call write_infile(infile,logfile,evfile,dumpfile,iprint,iprint)

!----------------------------------------------
! 2nd header after particles have been setup
!----------------------------------------------

 elseif (icall==2) then
!
!--number of particles
!
    if (present(ntot)) then
       write(iprint,"(/,' Number of particles = ',i12)") ntot
       do i = 1,6
          if (npartoftype(i) > 0) then
             if (i==igas) then
                parttype = "gas"
             else if (i==idust) then
                parttype = "dust"
             else if (i==iboundary) then
                parttype = "boundary"
             else if (i==istar) then
                parttype = "star"
             else if (i==idarkmatter) then
                parttype = "dark matter"
             else if (i==ibulge) then
                parttype = "bulge star"
             endif
             write(iprint,"(1x,3a,i12,a,es14.6)") &
                "Number & mass of ",parttype," particles: ", npartoftype(i),", ",massoftype(i)
          endif
       enddo
       write(iprint,"(a)") " "
    endif
    if (periodic) then
       write(iprint,"(1x,a)") 'Periodic boundaries: '
       if (abs(xmin) > 1.0d4 .or. abs(xmax) > 1.0d4 .or. &
           abs(ymin) > 1.0d4 .or. abs(ymax) > 1.0d4 .or. &
           abs(zmin) > 1.0d4 .or. abs(zmax) > 1.0d4      ) then
          write(iprint,"(2x,2(a,es14.6))") 'xmin = ',xmin,' xmax = ',xmax
          write(iprint,"(2x,2(a,es14.6))") 'ymin = ',ymin,' ymax = ',ymax
          write(iprint,"(2x,2(a,es14.6))") 'zmin = ',zmin,' zmax = ',zmax
       else
          write(iprint,"(2x,2(a,f10.6))")  'xmin = ',xmin,' xmax = ',xmax
          write(iprint,"(2x,2(a,f10.6))")  'ymin = ',ymin,' ymax = ',ymax
          write(iprint,"(2x,2(a,f10.6))")  'zmin = ',zmin,' zmax = ',zmax
       endif
    else
       write(iprint,"(a)") ' No boundaries set '
    endif

    write(iprint,"(/,a)") ' Using '//trim(kernelname)//' kernel'

    Nneigh = nint(4./3.*pi*(radkern*hfact)**3)
    write(iprint,50) hfact, massoftype(1), tolh, Nneigh
50  format(/,' Variable smoothing length: ',/, &
      6x,' h = ',f5.2,'*[',es9.2,'/rho]^(1/3); htol = ',es9.2,/ &
      6x,' Number of neighbours = ',i4)

!
!--MHD compile time options
!
    if (mhd) then
       if (maxBevol==4) then
          write(iprint,60) 'B/rho with cleaning'
       else
          write(iprint,60) 'B/rho'
       endif
60     format(/,' Magnetic fields are ON, evolving ',a)
    endif
    if (gravity)     write(iprint,"(1x,a)") 'Self-gravity is ON'
    if (h2chemistry) write(iprint,"(1x,a)") 'H2 Chemistry is ON'
    if (use_dustfrac) write(iprint,"(1x,a)") 'One-fluid dust is ON'
    if (use_dustgrowth) write(iprint,"(1x,a)") 'Dust growth is ON'

    call eosinfo(ieos,iprint)

    if (maxalpha==maxp) then
       if (nalpha >= 2) then
          write(iprint,"(2(a,f10.6))") ' Art. viscosity w/Cullen & Dehnen switch    : alpha  = ',alpha,' ->',alphamax
       else
          write(iprint,"(2(a,f10.6))") ' Art. visc. w/Morris & Monaghan switch      : alpha  = ',alpha,' ->',alphamax
       endif
    else
       write(iprint,"(a,f10.6)") ' Artificial viscosity                       : alpha  = ',alpha
    endif
    if (mhd) then
       write(iprint,"(a,f10.6)") ' Artificial resistivity, vsig=|vab x rab|   : alphaB = ',alphaB
    endif
    if (maxvxyzu >= 4) then
       if (gravity) then
          write(iprint,"(a,f10.6)") ' Art. conductivity w/divv switch (gravity)  : alphau = ',alphau
       else
          write(iprint,"(a,f10.6)") ' Art. conductivity w/Price 2008 switch      : alphau = ',alphau
       endif
    endif
    write(iprint,*)

!
!--if physical viscosity is set, print info about this
!
    call viscinfo(irealvisc,iprint)

#ifdef DUST
    call print_dustinfo(iprint)
#ifdef DUSTGROWTH
    call print_growthinfo(iprint)
#endif
#endif
!
!  print units information
!
    call print_units

 endif

 return
end subroutine write_header

end module writeheader
