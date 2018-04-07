!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testbinary
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: physcon, setbinary
!+
!--------------------------------------------------------------------------
module testbinary
 implicit none
 real, parameter :: pi = 4.*atan(1.)

contains

!----------------------------------------------------
!+
!  Numerically integrate a binary orbit given
!  input orbital elements. Print resulting 3D
!  track (in cartesian coordinates) to file
!+
!----------------------------------------------------
subroutine test_binary(m1,m2,a,e,inc,o,w,f,jfile,itex)
 use setbinary, only:set_binary
 use physcon,   only:au,solarm,gg,years
 real, intent(in) :: m1,m2,a,e,inc,o,w,f
 integer, intent(in) :: jfile,itex
 real :: xyz(3,2), vxyz(3,2), fxyz(3,2), xyzmh(6,2)
 real :: mrat,h1,h2
 real :: t,dt,period,period_yrs,proj_max,tsep,dx_proj(2),r_proj
 real(8) :: utime,udist,umass
 integer :: i,n,nsteps,lu
 character(len=48) :: filename

 mrat = m2/m1
 h1 = 0.01
 h2 = 0.01
 n = 0
 call set_binary(m1,mrat,a,e,h1,h2,xyzmh,vxyz,n,posang_ascnode=o,arg_peri=w,incl=inc,f=f,verbose=.false.)
 xyz(1:3,1) = xyzmh(1:3,1)
 xyz(1:3,2) = xyzmh(1:3,2)

 t = 0.
 period = sqrt(4.*pi**2*a**3/(m1+m2))
 nsteps = 2000
 dt = period/nsteps
 udist = au
 umass = solarm
 utime = sqrt(udist**3/(gg*umass))
 period_yrs = period*utime/years
 if (itex /= 0) then ! write line to orbits.tex file
    write(itex,"(6(f6.2,1x,'&',1x),f8.1,1x)") a,e,inc,o,w,f,period_yrs
 endif
!  write(*,"(6(a,'=',f6.2,1x),a,'=',f8.1,1x)") 'a',a,'e',e,'i',inc,'o',o,'w',w,'f',f,'P',period_yrs
 !print*,'period = ',period
 call get_f(m1,m2,xyz,fxyz)

 ! write trajectory of orbit to file
 write(filename,"(a,i3.3,a)") 'orbit',jfile,'.out'
 print "(a,6(f6.2,1x))",'writing to '//trim(filename)//' aeiowf = ',a,e,inc,o,w,f
 open(newunit=lu,file=filename,status='replace')
 tsep = 0.
 proj_max = 20.
 do i=1,nsteps
    t = (i-1)*dt
    vxyz = vxyz + 0.5*dt*fxyz
    xyz  = xyz  + dt*vxyz
    call get_f(m1,m2,xyz,fxyz)
    vxyz = vxyz + 0.5*dt*fxyz
    dx_proj = xyz(1:2,2) - xyz(1:2,1)
    r_proj  = sqrt(dot_product(dx_proj,dx_proj))
    if (r_proj < proj_max) then
       tsep = tsep + dt
    endif
    !print*,t
    write(lu,*) xyz(1:3,2) - xyz(1:3,1)
 enddo
 close(lu)
 !print*,'% of time spent at projected separation < ',proj_max,' = ',tsep/period

end subroutine test_binary

!----------------------------------------------------
!+
!  get acceleration for integration of binary orbit
!+
!----------------------------------------------------
subroutine get_f(m1,m2,x,fx)
 real, intent(in) :: m1,m2,x(3,2)
 real, intent(out) :: fx(3,2)
 real :: dx(3),r,r2

 dx = x(:,1) - x(:,2)
 r2 = dot_product(dx,dx)
 r  = sqrt(r2)
 fx(:,1) = -m2*dx/(r2*r)
 fx(:,2) =  m1*dx/(r2*r)

end subroutine get_f

end module testbinary
