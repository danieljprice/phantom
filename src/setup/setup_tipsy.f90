!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, io, part, physcon, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup from Tipsy file
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,      only:maxp,tagline
 use boundary, only:set_boundary,xmin,ymin,zmin,dxbound,dybound,dzbound
 use io,       only:master
 use part,     only:igas
 use units,    only:set_units
 use physcon,  only:pc,solarm
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: nargs,ngas,nptmass
 character(len=20) :: tipsyfile
 real :: pmassprev,utherm,rhoi
 integer :: i,ndim
!
!--boundaries
!
 call set_boundary(0.,0.29,0.,0.29,0.,0.29)
!
!--units
!
 call set_units(dist=pc,mass=solarm,time=3.640914d12)
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 2) then
    if (id==master) then
       print "(a,/)",trim(tagline)
       print "(a)",' Usage: phantomsetup blah tipsyfile '
    endif
    stop
 endif
 call get_command_argument(2,tipsyfile)

 open(unit=15,file=tipsyfile,form = 'formatted')
 read(15,*) npart,ngas,nptmass! ntotal,ngas,nstars
 read(15,*) ndim                 ! ndimensions
 if (ndim /= 3) stop 'ndim /= 3'
 read(15,*) time              ! time

 print*,'npart = ',npart,' ngas = ',ngas,' nptmass = ',nptmass
 if (nptmass > 0) stop 'non-zero point masses'
 if (npart > maxp) stop 'npart > maxp'

 npartoftype(igas) = ngas
 pmassprev = 0.
 do i=1,npart
    read(15,*) massoftype(igas) !!pmass(i)
    if (abs(massoftype(igas)-pmassprev) > tiny(pmassprev)) then
       print*,' mass = ',massoftype(igas)
    endif
    pmassprev = massoftype(igas)
 enddo
 print*,'x'
 do i=1,npart
    read(15,*) xyzh(1,i)
 enddo
 print*,'y'
 do i=1,npart
    read(15,*) xyzh(2,i)
 enddo
 print*,'z'
 do i=1,npart
    read(15,*) xyzh(3,i)
 enddo
 print*,'vx'
 do i=1,npart
    read(15,*) vxyzu(1,i)
 enddo
 print*,'vy'
 do i=1,npart
    read(15,*) vxyzu(2,i)
 enddo
 print*,'vz'
 do i=1,npart
    read(15,*) vxyzu(3,i)
 enddo
 do i=1,npart
    read(15,*) rhoi !rho(i)\
    !xyzh(4,i) = hfact*(massoftype(1)/rhoi)**(1./3.)
    if (i < 10) print*,'rhoi = ',rhoi
 enddo
 do i=1,npart
    read(15,*) utherm !!vxyzu(4,i)
    if (i < 10) print*,'utherm = ',utherm
 enddo
 do i=1,npart
    read(15,*) xyzh(4,i)
 enddo
 close(15)

 print*,'uthermlast = ',utherm
 polyk = 2./3.*utherm
 print*,'polyk = ',polyk
! polyk = 0.
 hfact = 1.2
 gamma = 1.

end subroutine setpart

end module setup

