!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
program plot_kernel
!
! Plots the kernel functions
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: plot_kernel [no arguments]
!
! :Dependencies: giza, kernel, plotk
!
 implicit none

contains

real function wfunc(q)
 use kernel, only:wkern
 real, intent(in) :: q
 real :: q2

 q2 = q*q
 wfunc = wkern(q2,q)

end function wfunc

real function potk(q)
 real, intent(in) :: q

 if (q > 1.e-6) then
    potk = 1./q
 else
    potk = 1.e10
 endif

end function potk

real function fk(q)
 real, intent(in) :: q

 if (q > 1.e-6) then
    fk = 1./(q*q)
 else
    fk = 1.e10
 endif

end function fk

end program plot_kernel

program plot_kernel
 use giza
 use plotk,   only:potk,fk
 use kernel,  only:kernelname,radkern,get_kernel_grav1,get_kernel_grav2,cnormk
 integer, parameter :: npts = 1000
 real, parameter :: xmin = 0., xmax = 3.99
 real :: x(npts),w(npts),grw(npts),dphidh(npts),pot(npts),fsoft(npts)
 real :: xi
 integer :: id,i

 do i=1,npts
    xi = 2.*(i-1)*radkern/(npts-1)
    x(i) = xi
    call get_kernel_grav1(xi*xi,xi,w(i),grw(i),dphidh(i))
    call get_kernel_grav2(xi*xi,xi,grw(i),pot(i),fsoft(i))
 enddo
 w = cnormk*w
 grw = cnormk*grw

 id = giza_open_device('?','kernel')
 call giza_set_environment(xmin,xmax,-1.,1.5,0,0)
 call giza_label('r/h','W',trim(kernelname))

 call giza_line(npts,x,w)
! call giza_set_line_style(2)
! call giza_line(npts,x,grw)
 call giza_set_line_style(3)
 call giza_line(npts,x,fsoft)
 call giza_function_x(fk,npts,xmin,xmax,1)
 call giza_set_line_style(4)
 call giza_line(npts,x,-pot)
 call giza_function_x(potk,npts,xmin,xmax,1)
! call giza_set_line_style(5)
! call giza_line(npts,x,dphidh)

 call giza_close()

contains

end program plot_kernel
