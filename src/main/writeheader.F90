!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module writeheader
!
! writes runtime header to the logfile
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, boundary_dyn, cooling, dim, dust, eos, gitinfo,
!   growth, io, kernel, metric_tools, mpiutils, options, part, physcon,
!   readwrite_infile, units, viscosity
!
 implicit none
 public :: write_header,write_codeinfo

 private

contains

!-----------------------------------------------------------------
!+
!  pretty print code header
!+
!-----------------------------------------------------------------
subroutine write_codeinfo(iunit)
 use dim,     only:phantom_version_string
 use gitinfo, only:get_and_print_gitinfo
 integer, intent(in) :: iunit
!
!--write out code name, version and time
!
 write(iunit,10) trim(phantom_version_string)

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
!
!--write info on latest git commit
!
#ifdef MCFOST
 write(*,*) ""
 write(*,*) "--------------------------"
 write(*,*) "| This is Phantom+mcfost |"
 write(*,*) "--------------------------"
 write(*,*) ""
#endif

 call get_and_print_gitinfo(iunit)

end subroutine write_codeinfo

!-----------------------------------------------------------------
!+
!  This subroutine writes header to logfile / screen
!  with the main parameters used for the run
!
!  icall = 1 (before particle have been set up)
!  icall = 2 (after particle setup)
!+
!-----------------------------------------------------------------
subroutine write_header(icall,infile,evfile,logfile,dumpfile,ntot)
 use dim,              only:maxp,maxvxyzu,maxalpha,ndivcurlv,mhd_nonideal,nalpha,use_dust,&
        use_dustgrowth,gr,h2chemistry
 use io,               only:iprint
 use boundary,         only:xmin,xmax,ymin,ymax,zmin,zmax
 use boundary_dyn,     only:dynamic_bdy,rho_thresh_bdy,width_bkg
 use options,          only:tolh,alpha,alphau,alphaB,ieos,alphamax,use_dustfrac
 use part,             only:hfact,massoftype,mhd,gravity,periodic,massoftype,npartoftypetot,&
                            labeltype,maxtypes
 use mpiutils,         only:reduceall_mpi
 use eos,              only:eosinfo
 use cooling,          only:cooling_in_step,Tfloor,ufloor
 use readwrite_infile, only:write_infile
 use physcon,          only:pi,pc
 use kernel,           only:kernelname,radkern
 use viscosity,        only:irealvisc,viscinfo
 use units,            only:print_units,unit_density,udist,unit_ergg
 use dust,             only:print_dustinfo,drag_implicit
 use growth,           only:print_growthinfo
#ifdef GR
 use metric_tools,     only:print_metricinfo
#endif
 integer                      :: Nneigh,i
 integer,          intent(in) :: icall
 character(len=*), intent(in) :: infile,evfile,logfile,dumpfile
 integer(kind=8),  intent(in), optional :: ntot
 character(len=10) :: startdate, starttime

!-----------------------------------------------------------------------
! 1st header after options have been read, but before particle setup
!-----------------------------------------------------------------------

 if (icall==1) then

    call date_and_time(startdate,starttime)
    startdate = startdate(7:8)//'/'//startdate(5:6)//'/'//startdate(1:4)
    starttime = starttime(1:2)//':'//starttime(3:4)//':'//starttime(5:)
    write(iprint,"(' Started on ',a,' at ',a)") startdate,starttime

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
       do i = 1,maxtypes
          if (npartoftypetot(i) > 0) then
             write(iprint,"(1x,3a,i12,a,es14.6)") &
                "Number & mass of ",labeltype(i)," particles: ", npartoftypetot(i),", ",massoftype(i)
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
          write(iprint,"(2x,2(a,g12.5))")  'xmin = ',xmin,' xmax = ',xmax
          write(iprint,"(2x,2(a,g12.5))")  'ymin = ',ymin,' ymax = ',ymax
          write(iprint,"(2x,2(a,g12.5))")  'zmin = ',zmin,' zmax = ',zmax
       endif
    else
       write(iprint,"(a)") ' No boundaries set '
    endif
    if (dynamic_bdy) then
       write(iprint,"(a)") ' Using dynamic boundaries '
       write(iprint,"(2x,a,Es18.6)") ' Min density of relevant gas (cgs): ',rho_thresh_bdy*unit_density
       write(iprint,"(2x,a,2f10.5)") ' dx/pc on both sides of relevant gas: ',width_bkg(1,:)*udist/pc
       write(iprint,"(2x,a,2f10.5)") ' dy/pc on both sides of relevant gas: ',width_bkg(2,:)*udist/pc
       write(iprint,"(2x,a,2f10.5)") ' dz/pc on both sides of relevant gas: ',width_bkg(3,:)*udist/pc
    endif

    write(iprint,"(/,a)") ' Using '//trim(kernelname)//' kernel'

    Nneigh = nint(4./3.*pi*(radkern*hfact)**3)
    write(iprint,50) hfact, massoftype(1), tolh, Nneigh
50  format(6x,' h = ',f5.2,'*[',es9.2,'/rho]^(1/3); htol = ',es9.2,/ &
           6x,' Number of neighbours = ',i4)

    if (mhd)              write(iprint,"(1x,a)") 'Magnetic fields are ON, evolving B/rho with cleaning'
    if (gravity)          write(iprint,"(1x,a)") 'Self-gravity is ON'
    if (h2chemistry)      write(iprint,"(1x,a)") 'H2 Chemistry is ON'
    if (use_dustfrac) then
       write(iprint,"(1x,a)") 'One-fluid dust is ON'
    else
       if (drag_implicit) then
          write(iprint,"(1x,a)") 'Two-fluid dust implicit scheme is ON'
       else
          write(iprint,"(1x,a)") 'Two-fluid dust explicit scheme is ON'
       endif
    endif
    if (use_dustgrowth)   write(iprint,"(1x,a)") 'Dust growth is ON'
    if (cooling_in_step)  then
       write(iprint,"(1x,a)") 'Cooling is calculated in step'
    else
       write(iprint,"(1x,a)") 'Cooling is explicitly calculated in force'
    endif
    if (ufloor > 0.) then
       write(iprint,"(3(a,Es10.3),a)") ' WARNING! Imposing temperature floor of = ',Tfloor,' K = ', &
       ufloor*unit_ergg,' erg/g = ',ufloor,' code units'
    endif
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
       if (gr) then
          write(iprint,"(a,f10.6)") ' Art. conductivity                          : alphau = ',alphau
       elseif (gravity .and. .not. gr) then
          write(iprint,"(a,f10.6)") ' Art. conductivity w/divv switch (gravity)  : alphau = ',alphau
       else
          write(iprint,"(a,f10.6)") ' Art. conductivity w/Price 2008 switch      : alphau = ',alphau
       endif
    endif
    if (gr) write(iprint,"(a)") '    GR --- See Liptai & Price (2018) for implementaion of shock dissipation terms'
    write(iprint,*)

!
!--if physical viscosity is set, print info about this
!
    call viscinfo(irealvisc,iprint)

    if (use_dust) call print_dustinfo(iprint)
    if (use_dustgrowth) call print_growthinfo(iprint)

#ifdef GR
    call print_metricinfo(iprint)
#endif
!
!  print units information
!
    call print_units

 endif

end subroutine write_header

end module writeheader
