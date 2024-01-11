!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine reconstructing the dust grain size distribution
! from the moments evolved in the dust nucleation network
!
! :References: 
!    Siess, Homan, Toupin & Price (2022), A&A 667, A75
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, fileutils, part, physcon, units
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'dustmoments'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use dim,       only:do_nucleation
 use part,      only:nucleation,idK0,idK3,rhoh
 use units,     only:udist,unit_density
 use physcon,   only:micron,amu,pi
 use dust_formation, only:mass_per_H
 use reconstruct_from_moments, only:reconstruct_maxent,&
                                    integrand,fsolve_error
 use integrate,                only:integrate_trap_log,integrate_trap
 use table_utils,              only:logspace
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,j,k,ierr
 integer, parameter :: ngrid = 4001
 real :: khat(0:3),ki(0:3),kgot(0:3)
 real :: mu(4),x(ngrid),lambsol(4),f(ngrid),ftmp(ngrid),lambguess(4)
 real :: xmin,xmax,scalefac,sigma_log,factor,logxmean,a0,rhoi,prefac
 real, parameter :: Acarb = 12.011, rhocarb = 2.25
 character(len=20) :: filename
 logical, parameter :: debug = .false.

 if (.not.do_nucleation) then
    stop 'ERROR: need DUST_NUCLEATION=yes for this analysis type'
 endif

 xmin = 1e4  ! critical value of N, number of monomers
 xmax = 1e20  ! want to integrate to infinity, but this is good enough
 sigma_log = 0.1

 call logspace(x,xmin,xmax)
 print*,' xmin, xmax  = ',xmin,xmax,' micron in code units = ',micron/udist

 logxmean = 10.
 scalefac = 100.
 ftmp = scalefac/(x*sqrt(2.*pi*sigma_log**2))*exp(-0.5*(log(x)-logxmean)**2/sigma_log**2)
 ftmp = exp(log(100./(x*sqrt(2.*pi*sigma_log**2))) -0.5*(log(x)-logxmean)**2/sigma_log**2)

 prefac = scalefac/sqrt(2.*pi*sigma_log**2)
 factor = -0.5/sigma_log**2
 ftmp = exp(log(prefac) + factor*logxmean**2 + (factor*(-2*logxmean) - 1.)*log(x) + factor*(log(x)**2))

 lambguess = [log(prefac) + factor*logxmean**2,factor*(-2.*logxmean) - 1.,factor,0.]
 ftmp = exp(lambguess(1) + lambguess(2)*log(x) + lambguess(3)*log(x)**2 + lambguess(4)*log(x)**3)
 print*,' integrated log normal Gaussian=',integrate_trap_log(ngrid,x,ftmp), ' should be ',scalefac

 a0 = (3.*Acarb*amu/(4.*pi*rhocarb))**(1./3.)
 print*,' mass per H = ',mass_per_H,' g'
 print*,' a0 = ',a0,' cm = ',a0/micron,' micron'

 ierr = 0
 do i=1,npart
    rhoi = rhoh(xyzh(4,i),particlemass)
    khat(0:3) = nucleation(idK0:idK3,i)
    ki(0:3) = khat(0:3)*rhoi*unit_density/mass_per_H 
    if (debug) print "(/,a,i0)",' reconstructing on particle ',i
    if (debug) call print_moments(ki,a0,rhoi)

    mu = ki

    ! try a log normal as the initial guess
    if (.true.) then
       prefac = mu(1)/sqrt(2.*pi*sigma_log**2)
       logxmean = log(ki(3)/ki(0))
       sigma_log = 1.
       if (debug) print*,' guessing normalisation = ',prefac,' logxmean = ',logxmean, ' sigma_log = ',sigma_log
       factor = -0.5/sigma_log**2
       lambguess = [log(prefac) + factor*logxmean**2,factor*(-2.*logxmean) - 1.,factor,0.]
    else
       lambguess = lambsol
    endif

    call reconstruct_maxent(mu,x,f,lambsol,ierr,use_log=.true.,lambguess=lambguess)
    if (ierr /= 1) print*,' INFO: '//fsolve_error(ierr)//' on particle ',i

    !lambsol = lambguess

    !print*,' got lambsol = ',lambsol
    f = integrand(x, lambsol, 0)

    do k=0,size(mu)-1
       ftmp = f*x**(k/3.)
       kgot(k) = integrate_trap_log(ngrid,x,ftmp)
    enddo

    if (debug .or. ierr /= 1) then
       print*,' MOMENTS IN  : ',mu
       print*,' MOMENTS OUT : ',kgot
    endif
    if (debug) call print_moments(kgot,a0,rhoi)

     ! print results
    if (debug) then
       write(filename,"(a,i5.5,a)") 'recon',i,'.out'
       open(unit=1,file=filename,status='replace',action='write')
       write(1,"(a)") '# N,f(N)'
       do j=1,ngrid
          write(1,*) x(j),f(j)
       enddo
       close(1)
    endif
 enddo

end subroutine do_analysis

subroutine print_moments(ki,a0,rhoi)
 use physcon, only:micron,pi
 real, intent(in) :: ki(0:3),a0,rhoi

 print "(2x,a,1pg0.4,a)",         'a) number density of dust = ',ki(0),' cm^-3'
 print "(2x,a,1pg0.4,a,1pg0.2,a)",'b) average grain radius   = ',ki(1)/ki(0),' times a0, or ',ki(1)/ki(0)*a0/micron, ' micron'
 print "(2x,a,1pg0.4,a,1pg0.2,a)",'c) average grain area     = ',4.*pi*a0**2*ki(2)/ki(0)/micron**2,' micron^2'
 print "(2x,a,1pg0.4,a,1pg0.2,a)",'d) average particle size  = ',ki(3)/ki(0),' monomers'
 print "(2x,a,1pg0.4,a,1pg0.2,a)",'e) opacity at Td=100      = ',ki(3)*pi*a0**3/rhoi*6.7*100.  ! Eq (42-43) of Siess et al. 2022

end subroutine print_moments

end module analysis
