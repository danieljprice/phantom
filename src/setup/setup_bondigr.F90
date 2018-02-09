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
!  DEPENDENCIES: centreofmass, io, options, part, physcon, setup_params,
!    spherical, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 use dim,            only:gr
 use physcon,        only:pi
 use externalforces, only:accradius1
 use metric,         only:mass1
 implicit none

 public :: setpart

 private

 real :: rcrit, C1,C2,n,Tc
 real :: wind_gamma = 5./3.

contains

!----------------------------------------------------------------
!+
!  setup for bondi accretion
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master,fatal
 use spherical,    only:set_sphere
 use options,      only:ieos,iexternalforce,nfulldump
 use timestep,     only:tmax,dtmax
 use centreofmass, only:reset_centreofmass
 use units,        only:udist,umass,utime,set_units
 use physcon,      only:pc,solarm,gg
 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,igas,set_particle_type,iboundary
 use stretchmap,   only:get_mass_r
 use kernel,       only:radkern
 use externalforces,only:accradius1_hard
 use prompting,    only:prompt
 use metric,       only:imetric
 use metric_tools, only:imet_schwarzschild
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: psep,tff
 real    :: vol,rmax,rmin
 real    :: gcode
 integer :: i,np,nx,ilattice,maxvxyzu,nbound
 real    :: r,pos(3),cs2
 real    :: totmass
 real    :: approx_m,approx_h
 real :: rho,v,u

!-- Set code units
 ! call set_units(mass=solarm,G=1.d0,c=1.d0)
 call set_units(G=1.,c=1.)

 rcrit = 8. ! Default
 call prompt(' Enter the critical point rcrit in units of central mass M: ',rcrit,2.+1e-5)
 rcrit = rcrit*mass1

 call compute_constants(rcrit)

 gcode = gg*umass*utime**2/udist**3
 print*,' gcode = ',gcode

 maxvxyzu = size(vxyzu(:,1))

!--Set general parameters
 time  = 0.
 gamma = wind_gamma
 polyk = 1.


!--- Setup particles
 rmax = 20.*mass1
 np   = 10*1000
 vol  = 4./3.*pi*rmax**3
 nx   = int(np**(1./3.))
 psep = vol**(1./3.)/real(nx)

 accradius1 = 7.*mass1
 approx_m    = get_mass_r(rhofunc,rmax,accradius1)/np
 approx_h    = hfact*(approx_m/rhofunc(accradius1))**(1./3.)
 ! accradius1_hard = (accradius1 - radkern*approx_h)
 accradius1_hard = 2.5*mass1

 rmin  = 2.6!accradius1_hard
 if(rmin<=2.*mass1) STOP 'rmin is less than Schwarzschild radius!'
 totmass = get_mass_r(rhofunc,rmax,rmin)
 rhozero = totmass/vol
 tff     = sqrt(3.*pi/(32.*rhozero))

 tmax = 10.*tff
 dtmax = tmax/500.

 print*,''
 print*,' Setup for gas: '
 print*,' min,max radius = ',rmin,rmax
 print*,' volume         = ',vol        ,' particle separation = ',psep
 print*,' vol/psep**3    = ',vol/psep**3,' totmass             = ',totmass
 print*,' free fall time = ',tff        ,' tmax                = ',tmax
 print*,''

 ilattice       = 2
 ! nfulldump      = 1
 iexternalforce = 1

 if (.not.gr) call fatal('setup_bondiwind','This setup only works with GR on')
 if (imetric/=imet_schwarzschild) call fatal('setup_bondiwind','This setup is meant for use with the Schwarzschild metric')


!--- Add stretched sphere
 npart = 0
 npart_total = 0
 call set_sphere('closepacked',id,master,rmin,rmax,psep,hfact,npart,xyzh,rhofunc=rhofunc,nptot=npart_total)
 massoftype(:) = totmass/npart
 print*,' npart = ',npart
 print*,''

 nbound = 0
 do i=1,npart
    pos = xyzh(1:3,i)
    r = sqrt(dot_product(pos,pos))
    call get_solution(rho,v,u,r)
    vxyzu(4,i)   = u
    vxyzu(1:3,i) = -v*pos/r
    ! attempt at outer boundary
    ! if (r + radkern*xyzh(4,i)>rmax) then
    !    call set_particle_type(i,iboundary)
    !    nbound = nbound + 1
    ! else
    !    call set_particle_type(i,igas)
    ! endif
 enddo

!--- Reset centre of mass to the origin
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 npartoftype(:) = 0
 npartoftype(igas) = npart_total-nbound
 npartoftype(iboundary) = nbound

end subroutine setpart
!======================================================================================================
subroutine compute_constants(rcrit)
 real, intent(in) :: rcrit
 real :: uc2,vc2

 n   = 1./(wind_gamma-1.)
 uc2 = mass1/(2.*rcrit)
 vc2 = uc2/(1.-3.*uc2)
 Tc  = vc2*n/(1.+n-vc2*n*(1.+n))

 C1  = sqrt(uc2) * Tc**n * rcrit**2
 C2  = (1. + (1.+n)*Tc)**2 * (1. - 2.*mass1/rcrit + C1**2/(rcrit**4*Tc**(2.*n)))
 print*,'Constants have been set'
 print*,'C1 = ',C1
 print*,'C2 = ',C2

end subroutine compute_constants

subroutine get_solution(rho,v,u,r)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r
 real, parameter :: adiabat = 1.
 real :: T,uvel,term,u0,dens,sqrtg

 ! Given an r, solve eq 76 for T numerically
 call Tsolve(T,r)

 uvel = C1/(r**2 * T**n)
 dens = adiabat*T**n
 u = T*n

 !get u0 at r
 term = 1./(1.-2.*mass1/r)
 u0  = sqrt(term*(1.+term*uvel**2))
 v   = uvel/u0

 sqrtg = 1. !???? FIX
 rho = sqrtg*u0*dens

end subroutine get_solution

! Newton Raphson
subroutine Tsolve(T,r)
   real, intent(in) :: r
   real, intent(out) :: T
   real :: Tnew, diff
   logical :: converged
   integer :: its
   integer, parameter :: itsmax = 100
   real, parameter :: tol = 1.e-5

   T = Tc !Guess

   converged = .false.
   its = 0
   do while (.not.converged .and. its<itsmax)
      Tnew = T - ffunc(T,r)/df(T,r)
      diff = abs(Tnew - T)/abs(T)
      converged = diff < tol
      T = Tnew
      its = its+1
   enddo

   if (.not. converged) print*,'not converged for r =',r,'. diff = ',diff

   ! print*,'Found T to be:',T,' with ',its,' iterations'

end subroutine Tsolve

real function ffunc(T,r)
  real, intent(in) :: T,r
  ffunc = (1. + (1. + n)*T)**2*(1. - (2.*mass1)/r + C1**2/(r**4*T**(2.*n))) - C2
end function ffunc

real function df(T,r)
   real, intent(in) :: T,r
   df = (2.*(1. + T + n*T)*((1. + n)*r**3*(-2.*mass1 + r) - C1**2*T**(-1. - 2.*n)*(n + (-1. + n**2)*T)))/r**4
end function df

real function rhofunc(r)
 real, intent(in) :: r
 real :: rho,v,u
 call get_solution(rho,v,u,r)
 rhofunc = rho
end function rhofunc


end module setup
