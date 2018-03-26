!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: multirun
!
!  DESCRIPTION: This program generates a series of input files
! varying a particular parameter in each
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: multirun infile
!
!  DEPENDENCIES: dim, io, options, readwrite_infile, timestep, viscosity
!+
!--------------------------------------------------------------------------
program multirun
 use dim,              only:tagline
 use readwrite_infile, only:read_infile,write_infile
 use options,          only:alpha,set_default_options
 use io,               only:set_io_unit_numbers,iprint,iwritein,formatreal
 use timestep,         only:C_cour,C_force
 use viscosity,        only:irealvisc,shearparam
 implicit none
 integer :: i,nargs,nruns,np,iampl,idot
 character(len=20) :: infile,logfile,evfile,dumpfile
 character(len=60) :: outfile
 character(len=12) :: string
 real :: fac,alphamin,dalpha,alphaSS,C_courin,C_forcein
 real :: HonR,facHonR

 call set_io_unit_numbers
 call set_default_options
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs /= 1) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: multirun infile'
    stop
 endif

 call get_command_argument(1,infile)
 call read_infile(infile,logfile,evfile,dumpfile)

 C_courin = C_cour
 C_forcein = C_force

 print "(a)",' Enter resolution for runs (in units of 10^6 particles) '
 read*,np
 fac = (np/2.)**(1./3.)

 print "(a)",' enter which series to run (1=LP10 seq hires only, 2=LP10 seq A=0.5, 3=LP10 seq AV, 4=LP10 seq flebbe, 5=psi/alpha)'
 read*,iampl

 print "(a)",' Enter H/R (default is 0.02)'
 read*,HonR
 facHonR = HonR/0.02
 alphaSS = 0.

 nruns = 15
 do i=1,nruns
    select case(iampl)
    case(1)
       select case(i)
       case(1)
          alpha = 1.7*fac  ! alpha = 0.07
       case(2)
          alpha = 2.5*fac  ! alpha = 0.14
       case(3)
          alpha = 5.0*fac  ! alpha = 0.23
       case(4)
          alpha = 7.5*fac  ! alpha = 0.28
       end select
    case(2)
       select case(i)
       case(1)
          alpha = 2.5*fac
       case(2)
          alpha = 5.0*fac
       end select
    case(3,4)
       alphamin = 0.03
       dalpha = 0.02
       alphaSS = alphamin + (i-1)*dalpha
       if (iampl==4) then
          irealvisc = 2
          shearparam = alphaSS
          alpha = 0.
       else
          irealvisc = 0
          shearparam = 0.
          alpha = alphaSS*(5./0.23)*fac
       endif
       if (alpha > 5.) then
          C_cour = C_courin/(alpha/5.)
          C_force = C_forcein/(alpha/5.)
          print*,' C_cour = ',C_cour,' C_force = ',C_force
       else
          C_cour = C_courin
          C_force = C_forcein
       endif
    case(5)
       alphamin = 0.01
       dalpha = 0.01
       alphaSS = alphamin + (i-1)*dalpha
       alpha = alphaSS/0.064537877668*fac*facHonR
    case default
       stop 'no parameter set known for this ampl'
    end select

    idot = index(infile,'.in')
    if (idot <= 0) idot = len_trim(infile)

    if (alphaSS > 0.) then
       if (iampl==5) then
          write(string,"(f5.3)") alphaSS
       else
          call formatreal(alphaSS,string)
       endif
       outfile = infile(1:idot-1)//'_alphaSS'//trim(string)//'.in'
    else
       call formatreal(alpha,string)
       outfile = infile(1:idot-1)//'_alpha'//trim(string)//'.in'
    endif

    print *,trim(outfile),' : alpha = ',alpha
    call write_infile(trim(outfile),logfile,evfile,dumpfile,iwritein,iprint)
 enddo

 print "(a)",' wishing you a pleasant run time'

end program multirun
