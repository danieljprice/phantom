!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! moddump addplanets routine: add some planets in the disc
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - nplanets : *number of planets to add*
!
! :Dependencies: centreofmass, infile_utils, moddump_utils, part, physcon,
!   prompting, units
!
 use moddump_utils, only:prompt_for_params,write_moddump_header
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer, parameter :: maxplanets = 9
 character(len=*), dimension(maxplanets), parameter :: planets = &
  (/'1','2','3','4','5','6','7','8','9' /)

 integer :: nplanets = 1
 real    :: mplanet(maxplanets), rplanet(maxplanets), accrplanet(maxplanets)

 public :: init_moddump,read_moddump,write_moddump
 logical, parameter :: moddump_interactive = .true.
 public :: moddump_interactive

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:nptmass,xyzmh_ptmass,vxyz_ptmass,igas,ihacc,ihsoft
 use units,             only:umass,utime,udist,print_units
 use physcon,        only:au,solarm,jupiterm,pi,years
 use centreofmass,      only:reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,j
 real    :: phi,vphi,sinphi,cosphi,omega,r2,disc_m_within_r,star_m
 real    :: rsinkmass(maxplanets+1)

 print "(1x,45('-'))"
 print "(1x,45('-'))"

 call print_units

 print*, 'Current number of sink masses: ', nptmass

 do i=1,nptmass
!  print *,' sink mass: ',i,', xyzmh = ',xyzmh_ptmass(:,i)
    rsinkmass(i) = sqrt( xyzmh_ptmass(1,i)**2 + xyzmh_ptmass(2,i)**2 + xyzmh_ptmass(3,i)**2)
    if (i==1) then
       print *,' -Sink mass ',i,' (central star):'
!   print *,'    distance from the origin: ',rsinkmass(i)
       print *,'    mass: ', xyzmh_ptmass(4,i)
       print *,'    accretion radius: ',xyzmh_ptmass(5,i)
       print *,' '
    else
       print *,' -Planet ',i-1,' :'
       print *,'    distance from the origin: ',rsinkmass(i)
       print *,'    mass: ', xyzmh_ptmass(4,i)
       print *,'    accretion radius: ',xyzmh_ptmass(5,i)
       print *,' '
    endif

 enddo

 star_m = xyzmh_ptmass(4,1)
! print*,'Mass of central star: ', star_m
 !
 !--prompt only when the driver did not find a parameter file
 !
 if (prompt_for_params) call read_interactive_moddumpfile()

 print "(a,i2,a)",' --------- added ',nplanets,' planets ------------'
 do i=1,nplanets
    nptmass = nptmass + 1
    phi = (i-1)*180.*pi/180.   ! phi = 0.*pi/180.
    cosphi = cos(phi)
    sinphi = sin(phi)
    disc_m_within_r = 0.
    do j=1,npart
       r2 = xyzh(1,j)**2 + xyzh(2,j)**2 + xyzh(3,j)**2
       if (r2 < rplanet(i)**2) then
          disc_m_within_r = disc_m_within_r + massoftype(igas)
       endif
    enddo
    xyzmh_ptmass(1:3,nptmass)   = (/rplanet(i)*cosphi,rplanet(i)*sinphi,0./)
    xyzmh_ptmass(4,nptmass)     = mplanet(i)
    xyzmh_ptmass(ihacc,nptmass)  = accrplanet(i) ! 0.25*au/udist
    xyzmh_ptmass(ihsoft,nptmass) = accrplanet(i) ! 0.25*au/udist
    vphi = sqrt((star_m + disc_m_within_r)/rplanet(i))
    vxyz_ptmass(1:3,nptmass)    = (/-vphi*sinphi,vphi*cosphi,0./)
    print "(a,i2,a)",       ' planet ', nptmass-1 ,':'!,i,':'
    print "(a,g10.3,a)",    ' radius: ',rplanet(i)*udist/au,' AU'
    print "(a,g10.3,a,2pf7.3,a)",    ' M(<R) : ',(disc_m_within_r + star_m)*umass/solarm, &
          ' MSun, disc mass correction is ',disc_m_within_r/star_m,'%'
    print "(a,2(g10.3,a))", ' mass  : ',mplanet(i)*umass/jupiterm,' MJup, or ',mplanet(i)*umass/solarm,' MSun'
    print "(a,2(g10.3,a))", ' period: ',2.*pi*rplanet(i)/vphi*utime/years,' years or ',2*pi*rplanet(i)/vphi,' in code units'
    omega = vphi/rplanet(i)
    print "(a,g10.3,a)",   ' resonances: 3:1: ',(sqrt(star_m)/(3.*omega))**(2./3.),' AU'
    print "(a,g10.3,a)",   '             4:1: ',(sqrt(star_m)/(4.*omega))**(2./3.),' AU'
    print "(a,g10.3,a)",   '             5:1: ',(sqrt(star_m)/(5.*omega))**(2./3.),' AU'
    print "(a,g10.3,a)",   '             9:1: ',(sqrt(star_m)/(9.*omega))**(2./3.),' AU'
 enddo
 print "(1x,45('-'))"
 print "(1x,45('-'))"

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

end subroutine modify_dump

!
!---Interactively set the moddump parameters--------------------------------
!
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt
 use part,      only:nptmass
 integer :: i

 call prompt('Enter the number of planet you want to add: ', nplanets, 0, maxplanets-nptmass+1)
 do i=1,nplanets
    call prompt('Enter mass (code units) of planet '//trim(planets(i))//' :',mplanet(i),0.)
    call prompt('Enter distance from the central star (code units) of planet '//trim(planets(i))//' :',rplanet(i),0.)
    call prompt('Enter accretion radius (code units) of planet '//trim(planets(i))//' :',accrplanet(i),0.)
 enddo

end subroutine read_interactive_moddumpfile

!
!---Write the moddump parameter file----------------------------------------
!
subroutine write_moddump(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer :: i

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 call write_inopt(nplanets,'nplanets','number of planets to add',iunit)
 do i=1,nplanets
    call write_inopt(mplanet(i),'mplanet'//trim(planets(i)),'mass of planet [code units]',iunit)
    call write_inopt(rplanet(i),'rplanet'//trim(planets(i)),'distance from central star [code units]',iunit)
    call write_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),'accretion radius of planet [code units]',iunit)
 enddo
 close(iunit)

end subroutine write_moddump

!
!---Read the moddump parameter file-----------------------------------------
!
subroutine read_moddump(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr,i
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(nplanets,'nplanets',db,min=0,max=maxplanets,errcount=nerr)
 do i=1,nplanets
    call read_inopt(mplanet(i),'mplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
    call read_inopt(rplanet(i),'rplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
    call read_inopt(accrplanet(i),'accrplanet'//trim(planets(i)),db,min=0.,errcount=nerr)
 enddo
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddump

subroutine init_moddump()
 use physcon, only:au
 use units, only:udist
 integer :: i

 do i=1,maxplanets
    mplanet(i)    = 0.001
    rplanet(i)    = 10.*i*au/udist
    accrplanet(i) = 0.25*au/udist
 enddo
end subroutine init_moddump

end module moddump

