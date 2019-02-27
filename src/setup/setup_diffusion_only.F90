!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Setup of a linear sound wave in a box
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dust, io, kernel, mpiutils, options, part,
!    physcon, prompting, set_dust, setup_params, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:labeltype,set_particle_type,igas,idust,radenergy
 use physcon,      only:pi
 use kernel,       only:radkern
 use dim,          only:maxvxyzu,use_dust,maxp
 use options,      only:nfulldump,use_dustfrac
 use timestep,     only:dtmax,tmax
 use prompting,    only:prompt
 use dust,         only:K_code,idrag
 use set_dust,     only:set_dustfrac
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: totmass,fac,deltax,deltay,deltaz
 integer :: i
 integer :: itype,itypes,ntypes,npartx
 integer :: npart_previous,dust_method
 logical, parameter :: ishift_box =.true.
 real, parameter    :: dust_shift = 0.
 real    :: xmin_dust,xmax_dust,ymin_dust,ymax_dust,zmin_dust,zmax_dust
 real    :: kwave,denom,length,uuzero,przero !,dxi
 real    :: xmini,xmaxi,ampl,cs,dtg,massfac

  npartx  = 64
  ntypes  = 1
  rhozero = 1.
  massfac = 1.
  cs      = 1.
  ! ampl    = 1.d-4
  use_dustfrac = .false.
  if (id==master) then
    itype = 1
    print "(/,a,/)",'  >>> Setting up particles for diffusion only test <<<'
    call prompt(' enter number of '//trim(labeltype(itype))//' particles in x ',npartx,8,int(maxp/144.))
  endif
  call bcast_mpi(npartx)
  !
  ! boundaries
  !
  xmini = -0.5
  xmaxi =  0.5
  length = xmaxi - xmini
  deltax = length/npartx
  ! try to give y boundary that is a multiple of 6 particle spacings in the low density part
  fac = 6.*(int((1.-epsilon(0.))*radkern/6.) + 1)
  deltay = fac*deltax*sqrt(0.75)
  deltaz = fac*deltax*sqrt(6.)/3.
  call set_boundary(xmin,xmax,-deltay,deltay,-deltaz,deltaz)

  gamma = 5./3.

  npart = 0
  npart_total = 0
  npartoftype(:) = 0

  itype = igas
  rhozero = 1.0
  cs = 1.0
  polyk = 0.

  nfulldump = 1
  dtmax     = 0.1
  tmax      = 1.0
  call bcast_mpi(rhozero)

  call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                    hfact,npart,xyzh,nptot=npart_total,rhofunc=rhofunc)

 !--set which type of particle it is
  do i=1,npart
    call set_particle_type(i,itype)
    vxyzu(1:3,i) = 0.
    vxyzu(4,i) = cs**2/(gamma*(gamma-1.))
    if (xyzh(1,i) < 0.0) then
      radenergy(i) = 1.0
    else
      radenergy(i) = 2.0
    end if
  enddo

  npartoftype(itype) = npart - npart_previous
  if (id==master) print*,' npart = ',npart,npart_total

  totmass = massfac*rhozero*dxbound*dybound*dzbound
  if (id==master) print*,' box volume = ',dxbound*dybound*dzbound,' rhozero = ',rhozero

  massoftype(itype) = totmass/npartoftype(itype)
  if (id==master) print*,' particle mass = ',massoftype(itype)

contains

real function rhofunc(x)
  real, intent(in) :: x

  rhofunc = 1.
end function rhofunc

end subroutine setpart

end module setup
