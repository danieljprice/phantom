!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for Kelvin-Helmholtz instability from Robertson et al. (2010)
!
! :References:
!   Robertson et al. (2010), MNRAS 401, 2463-2476
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, io, mpidomain, mpiutils, options, part, physcon,
!   prompting, setup_params, timestep, unifdis
!
 implicit none
 public :: setpart

 private
 !--Hard-coded input parameters
 real, parameter :: rho1   =  1.   ! density  of medium 1
 real, parameter :: rho2   =  2.   ! density  of medium 2
 real            :: v1     = -0.5  ! velocity of medium 1
 real            :: v2     =  0.5  ! velocity of medium 2
 real            :: przero =  2.5  ! initial constant pressure
 real            :: xsize  =  1.0  ! size of the box in the x-direction
 real            :: ysize  =  1.0  ! size of the box in the y-direction
 real            :: dy2    =  0.5  ! Width of medium 2 (i.e. central medium)
contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:npart_total
 use io,           only:master
 use options,      only:nfulldump
 use unifdis,      only:set_unifdis,rho_func
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,&
                        dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:igas,periodic
 use prompting,    only:prompt
 use physcon,      only:pi
 use timestep,     only:dtmax,tmax
 use mpidomain,    only:i_belong
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=26)                :: filename
 logical :: iexist
 integer :: i,maxp,maxvxyzu,npartx
 real    :: totmass,deltax
 procedure(rho_func), pointer :: density_func
!
!--general parameters
!
 time  = 0.
 gamma = 5./3
 polyk = 0.
 filename= trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then
    tmax      = 2.00
    dtmax     = 0.1
    nfulldump = 1
 endif
!
!--set particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 if (id==master) then
    npartx = 64
    call prompt('enter number of particles in x direction ',npartx,1,nint(sqrt(maxp/12.)))
 endif
 call bcast_mpi(npartx)
 deltax = dxbound/npartx
!
!--boundary
!
 call set_boundary(0.,xsize,0.,ysize,-2.*sqrt(6.)/npartx,2.*sqrt(6.)/npartx)

 npart = 0
 npart_total = 0
 density_func => rhofunc
 call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,periodic,nptot=npart_total,&
                  rhofunc=density_func,dir=2,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart
 print*,' npart = ',npart,npart_total

 totmass = dy2*dxbound*dzbound*rho2 + (dybound-dy2)*dxbound*dzbound*rho1
 massoftype(igas) = totmass/npart_total
 print*,' particle mass = ',massoftype(igas)

 do i=1,npart
    vxyzu(1,i) = v1 + Rfunc(xyzh(2,i))*(v2 - v1)
    vxyzu(2,i) = 0.1*sin(2.*pi*xyzh(1,i))
    vxyzu(3,i) = 0.
    if (maxvxyzu > 3) then
       vxyzu(4,i) = przero/((gamma - 1.)*rhofunc(xyzh(2,i)))
    endif
 enddo

end subroutine setpart

!------------------------------------------
!+
!  desired density profile in y direction
!+
!------------------------------------------
real function rhofunc(y)
 real, intent(in) :: y

 rhofunc = rho1 + Rfunc(y)*(rho2 - rho1)

end function rhofunc

!---------------------------------------------------
!+
!  smoothing function, as per Robertson et al paper
!+
!---------------------------------------------------
real function Rfunc(y)
 real, parameter  :: delta = 0.05
 real, intent(in) :: y
 real :: fac1,fac2,yedgel,yedger

 yedgel = 0.5*ysize - 0.5*dy2
 yedger = 0.5*ysize + 0.5*dy2
 fac1 = (1. - 1./(1. + exp(2.*(y-yedgel)/delta)))
 fac2 = (1. - 1./(1. + exp(2.*(yedger-y)/delta)))
 Rfunc = fac1*fac2

end function Rfunc

end module setup
