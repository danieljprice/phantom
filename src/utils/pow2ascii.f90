!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: pow2ascii
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: pow2ascii powfile
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
program pow2ascii
 implicit none
 integer :: nargs,i,mk,n,nbins,ierr,ibin
 integer, parameter :: iunitr = 1, iunitw = 2, mlabel = 60, mfile = 60
 integer, parameter :: maxk = 1024
 character(len=120) :: powfile,outfile
 character(len=mfile) :: origin
 character(len=mlabel) :: label
 real :: rho_power,ptot
 integer :: nk(maxk)
 real :: xk(maxk),power(maxk)

 namelist /size/mk,n
 namelist /describe/label
 namelist /files/origin
 namelist /component/rho_power,ptot

 nargs = command_argument_count()
 if (nargs < 1) then
    stop 'usage: pow2ascii powfile'
 endif

 do i=1,nargs
    call get_command_argument(i,powfile)

    open(unit=iunitr,file=powfile,form='formatted',status='old')
    print*,'reading file '//trim(powfile)
    read(iunitr,NML=size)
    nbins = mk
    print*,'numk = ',nbins
    read(iunitr,NML=files)
    print*,'origin = '//trim(origin)
    read(iunitr,NML=describe)
    print*,'label = '//trim(label)
    read(iunitr,*) xk(1:nbins)
    read(iunitr,NML=describe)
    print*,'label = '//trim(label)
    read(iunitr,*) nk(1:nbins)
    ierr = 0
    do while (ierr == 0)
       read(iunitr,NML=describe,iostat=ierr)
       print*,'label = '//trim(label)
       read(iunitr,NML=component,iostat=ierr)
       print*,'rho power = ',rho_power, ' ptot = ',ptot
       read(iunitr,*,iostat=ierr) power(1:nbins)
       write(outfile,"(a,f4.2,a)") powfile(1:index(powfile,'.pow')-1)//'_'// &
                                   trim(label)//'_rho',rho_power,'.pspec'
       print*,'writing to '//trim(outfile)
       open(iunitw,file=outfile,form='formatted',status='replace')
       write(iunitw,*) '# ',nbins,ptot,rho_power
       do ibin=1,nbins
          write(iunitw,*) xk(ibin),power(ibin),power(ibin)*xk(ibin), &
                          power(ibin)*xk(ibin)**(5./3.),power(ibin)*xk(ibin)**2, &
                          power(ibin)*xk(ibin)**1.5
       enddo
       close(iunitw)
    enddo
    close(iunitr)
 enddo

end program pow2ascii
