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
 use units,     only:unit_density
 use physcon,   only:micron,amu,pi
 use dust_formation, only:mass_per_H
 use reconstruct_from_moments, only:integrand,fsolve_error,reconstruct_gamma_dist,gamma_func,&
                                    gamma_func_moment,gamma_func_from_moments
 use integrate,                only:integrate_trap_log,integrate_trap,&
                                    gauss_legendre_nodes_weights_log,integrate_gauss_legendre
 use table_utils,              only:logspace
 use timing,                   only:getused,printused
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,j,k,ierr,n,iu,nbad
 integer, parameter :: ngrid = 2001
 real :: khat(0:3),ki(0:3),kgot(0:3)
 real :: mu(4),x(ngrid),w(ngrid),lambsol(4),f(ngrid),ftmp(ngrid),lambguess(4),err(4)
 real :: xmin,xmax,scalefac,sigma_log,factor,logxmean,a0,rhoi,prefac,errmax(4),d_on_p,theta,p
 real(4) :: t1
 real, parameter :: Acarb = 12.011, rhocarb = 2.25
 character(len=64) :: filename,string
 logical, parameter :: debug = .false., check_nans = .false.
 logical :: flag_particle,bad_soln

 if (.not.do_nucleation) then
    stop 'ERROR: need DUST_NUCLEATION=yes for this analysis type'
 endif

 xmin = 1e1  ! critical value of N, number of monomers
 xmax = 1e20  ! want to integrate to infinity, but this is good enough
 sigma_log = 0.1

 print*,' TESTING INTEGRATION METHODS...'
 call logspace(x,xmin,xmax)
 print*,' xmin, xmax  = ',xmin,xmax

 logxmean = 10.
 scalefac = 100.
 ftmp = scalefac/(x*sqrt(2.*pi*sigma_log**2))*exp(-0.5*(log(x)-logxmean)**2/sigma_log**2)
 ftmp = exp(log(100./(x*sqrt(2.*pi*sigma_log**2))) -0.5*(log(x)-logxmean)**2/sigma_log**2)

 prefac = scalefac/sqrt(2.*pi*sigma_log**2)
 factor = -0.5/sigma_log**2
 ftmp = exp(log(prefac) + factor*logxmean**2 + (factor*(-2*logxmean) - 1.)*log(x) + factor*(log(x)**2))

 lambguess = [log(prefac) + factor*logxmean**2,factor*(-2.*logxmean) - 1.,factor,0.]
 ftmp = exp(lambguess(1) + lambguess(2)*log(x) + lambguess(3)*log(x)**2 + lambguess(4)*log(x)**3)
 print*,' integrated log normal on log-spaced grid=',integrate_trap_log(ngrid,x,ftmp), ' should be ',scalefac

 call gauss_legendre_nodes_weights_log(x,xmin,xmax,w)
 ftmp = exp(lambguess(1) + lambguess(2)*log(x) + lambguess(3)*log(x)**2 + lambguess(4)*log(x)**3) 
 print*,' integration with Gauss-Legendre = ',integrate_gauss_legendre(ngrid,w,ftmp)

 a0 = (3.*Acarb*amu/(4.*pi*rhocarb))**(1./3.)
 print*
 print*,' mass per H = ',mass_per_H,' g'
 print*,' a0 = ',a0,' cm = ',a0/micron,' micron'
 print "(/a)",' Reconstructing from moments...'

 ierr = 0
 n = 0
 nbad = 0
 errmax = 0.
 call getused(t1)

 !$omp parallel do default (none) &
 !$omp private (i,j,k,iu,ierr,khat,ki,kgot,mu,lambsol,f,ftmp,lambguess,bad_soln) &
 !$omp private (rhoi,prefac,factor,logxmean,sigma_log,filename,err,flag_particle,string,scalefac,d_on_p,theta,p) &
 !$omp shared(xyzh,nucleation,npart,particlemass,unit_density,mass_per_H,a0,t1,x,n,w,nbad) &
 !$omp reduction(max:errmax) &
 !$omp schedule (dynamic)
 do i=1,npart
    rhoi = rhoh(xyzh(4,i),particlemass)
    khat(0:3) = nucleation(idK0:idK3,i)
    ki(0:3) = khat(0:3)*rhoi*unit_density/mass_per_H 
    if (debug) print "(/,a,i0)",' reconstructing on particle ',i
    if (debug) call print_moments(ki,a0,rhoi)

    mu = ki !khat

    ! try a log normal as the initial guess
    if (.false.) then
       sigma_log = 1.
       prefac = mu(1)/sqrt(2.*pi*sigma_log**2)
       logxmean = log(ki(3)/ki(0))
       if (debug) print*,' guessing normalisation = ',prefac,' logxmean = ',logxmean, ' sigma_log = ',sigma_log
       factor = -0.5/sigma_log**2
       lambguess = [log(prefac) + factor*logxmean**2,factor*(-2.*logxmean) - 1.,factor,0.]
    elseif (.true.) then
       lambguess = [mu(1),1.,1.7,1.]   ! normalisation, theta, d, p
       d_on_p = 2.0
       theta = (mu(2)/mu(1) * Gamma(d_on_p) / Gamma(d_on_p + 1./3.))**3
    else
       lambguess = lambsol
    endif

    call reconstruct_gamma_dist(mu,lambsol,err,ierr,lambguess=lambguess,verbose=debug)

    ! now check by numerical integration
    d_on_p = lambsol(1)
    p = lambsol(2)
    theta = (mu(2)/mu(1) * Gamma(d_on_p) / Gamma(d_on_p + 1./(3.*p)))**3
    lambsol = [mu(1),theta,d_on_p,p]

    !lambguess = lambsol
    !print*,' lambsol from recon_gamma_dist = ',lambsol
    !call reconstruct_maxent(mu,x,f,lambsol,err,ierr,use_log=.true.,lambguess=lambguess) !,weights=w)
    !lambsol = lambguess
    !
    !--record max fractional error in each moment
    !
    flag_particle = .false.
    if (any(lambsol(3:4) < 0.)) flag_particle = .true.
    if ((any(abs(err) > 0.1 .and. mu(1) > tiny(0.)))) flag_particle = .true.
 
    if (ierr /= 1 .and. (debug .and. flag_particle)) print*,' INFO: '//trim(fsolve_error(ierr))//&
       ' on particle ',i,' k3=',mu(4),' got ',gamma_func_moment(2,lambsol(3:4),mu,3)

    if (ierr /= 1 .and. (debug .and. flag_particle)) then
       f = integrand(x, lambsol, 0)
       !do k=1,ngrid
       !   f(k) = gamma_func_from_moments(x(k),mu(1:2),lambsol(3:4))
       !enddo
       do k=0,size(mu)-1
          ftmp = f*x**(k/3.)
          !kgot(k) = integrate_trap_log(ngrid,x,ftmp)
          kgot(k) = integrate_gauss_legendre(ngrid,w,ftmp)
       enddo
       print*,' lambsol = ',lambsol
       print*,' MOMENTS IN             : ',mu
       print*,' MOMENTS OUT (numerical): ',kgot
       scalefac = mu(1) * lambsol(4) / ((lambsol(2)) * Gamma(lambsol(3)))
       print*,' Moments OUT (analytic ): ',(scalefac * (lambsol(2)) ** (1.+k/3.) / &
                              lambsol(4) * Gamma(lambsol(3) + k/(3.*lambsol(4))),k=0,size(mu)-1)
       if (debug) call print_moments(kgot,a0,rhoi)
    endif

    if (flag_particle) then
       !$omp atomic
       nbad = nbad + 1
    endif

    !$omp atomic
    n = n + 1
    if (mod(n,10000)==0) then
       write(string,"(i0,a,i0,a)") n,' [',nbad,' with errors > 10%] completed in'
       call printused(t1,string)
    endif

    bad_soln = .false.
    if (check_nans) then
       f = integrand(x,lambsol,0)
       if (any(isnan(f))) then
          print*,i,' grain size function contains NaNs'
          bad_soln = .true.
       endif
       if (f(ngrid) > exp(-50.)) then
          print*,i,' function is exploding at the large N boundary',f(ngrid)
          bad_soln = .true.
       endif
    endif
    ! print results
!    if (.true.) then
    if (debug .and. (flag_particle .or. bad_soln)) then
       f = integrand(x, lambsol, 0)
       write(filename,"(a,i7.7,a)") 'recon',i,'.out'
       open(newunit=iu,file=filename,status='replace',action='write')
       write(iu,"(a)") '# N,f(N)'
       do j=1,ngrid
          write(iu,*) x(j),f(j)
       enddo
       close(iu)
    endif
    if (bad_soln) stop

    if (any(lambsol(3:4) < 0)) then
       print*,i,'ERROR: d_on_p,p = ',lambsol(3:4),' should be > 0'
       stop
    endif

 enddo
 !$omp end parallel do
 call printused(t1)

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
