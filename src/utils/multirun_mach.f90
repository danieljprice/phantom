!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
!  DEPENDENCIES: dim, forcing, io, options, readwrite_infile, timestep
!+
!--------------------------------------------------------------------------
program multirun
 use dim,              only:tagline
 use readwrite_infile, only:read_infile,write_infile
 use options,          only:set_default_options
 use timestep,         only:dtmax,tmax
 use io,               only:set_io_unit_numbers,iprint,iwritein,formatreal
 use forcing,          only:st_decay,st_energy,st_dtfreq,st_stirmin,st_stirmax,st_solweight
 implicit none
 integer :: i,nargs,nruns,idot
 character(len=20) :: infile,logfile,evfile,dumpfile
 character(len=20) :: machstring
 real :: xmach,tdyn

 call set_io_unit_numbers
 call set_default_options
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1 .or. nargs > 2) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: multirun infile'
    stop
 endif

 call get_command_argument(1,infile)
 call read_infile(infile,logfile,evfile,dumpfile)

 nruns = 10
 xmach = 0.
 do i=1,nruns
    xmach = xmach + 2.0d0

    tdyn = 1./(2.*xmach)
    st_decay  = tdyn
    st_dtfreq = 0.2*st_decay
    st_energy = 8.0*(xmach/10.)**2
    st_stirmin = 6.28
    st_stirmax = 18.86
    st_solweight = 1.0

    print*
    print*,' Mach number = ',xmach
    print*,' st_energy   = ',st_energy
    print*,' st_decay    = ',st_decay
    print*,' st_dtfreq   = ',st_dtfreq

    tmax  = 10.*tdyn
    dtmax = tmax/200.
    print*,'        tmax = ',tmax
    print*,'       dtmax = ',dtmax

    call formatreal(xmach,machstring)
    idot = index(infile,'.in')
    if (idot==0) idot = len_trim(infile)

    print*,' writing '//infile(1:idot-1)//'_mach'//trim(machstring)//'.in'
    call write_infile(infile(1:idot-1)//'_mach'//trim(machstring)//'.in',logfile,evfile,dumpfile,iwritein,iprint)
 enddo

 print "(a)",' wishing you a pleasant run time'

end program multirun
