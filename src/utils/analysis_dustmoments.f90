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
 use integrate,                only:integrate_trap_log
 use table_utils,              only:linspace,logspace
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,j,k,ierr
 integer, parameter :: ngrid = 4001
 real :: khat(0:3),ki(0:3),kgot(0:3)
 real :: mu(4),x(ngrid),lambsol(4),f(ngrid),ftmp(ngrid),lambguess(4)
 real :: xmin,xmax,scalefac,sigma_log,factor,logxmean,a0,rhoi
 real, parameter :: Acarb = 12.011, rhocarb = 2.25
 character(len=20) :: filename

 if (.not.do_nucleation) then
    stop 'ERROR: need DUST_NUCLEATION=yes for this analysis type'
 endif

 xmin = 1e4  ! critical value of N, number of monomers
 xmax = 1e20  ! want to integrate to infinity, but this is good enough
 sigma_log = 1.

 call logspace(x,xmin,xmax)
 print*,' xmin, xmax  = ',xmin,xmax,' micron in code units = ',micron/udist


 a0 = (3.*Acarb*amu/(4.*pi*rhocarb))**(1./3.)
 print*,' mass per H = ',mass_per_H,' g'
 print*,' a0 = ',a0,' cm = ',a0/micron,' micron'
 !lambguess = [78.148545951275779 ,      -79.077042328656518 ,       16.875452444253007  ,     -1.1286214050246528 ]
 do i=1,10
    rhoi = rhoh(xyzh(4,i),particlemass)
    khat(0:3) = nucleation(idK0:idK3,i)
    ki(0:3) = khat(0:3)*rhoi*unit_density/mass_per_H 
    print "(/,a,i0)",' reconstructing on particle ',i
    call print_moments(ki,a0,rhoi)

    mu = ki
    print*,' MOMENTS IN  : ',mu

    ! try a log normal as the initial guess
    logxmean = log(ki(3)/ki(0))
    !print*,' guessing logxmean = ',logxmean
    factor = -0.5/sigma_log**2
    lambguess = factor*[logxmean**2,-logxmean,1.,0.]
    !lambguess = [100.,-3.5,0.,0.]

    !lambguess = [1./sqrt(2.*pi*sigma_log),ki(3)/ki(0),1000.,0.]
    !lambguess =[-10., 0., 0., 0.]
    !lambguess = [14.488032400232779,       -15.263198315586006,        1.2546701587241598 ,      -2.9858628177166273E-002]
    !lambguess = [-12.672479403225283,       -12.291129733347596,        1.1116072069896121,       -2.6731093488116833E-002]
    !lambguess = [-61.535685209145377,       -4.5732464265546762,       0.71087277282226946,       -1.9962892974256529E-002]
    !lambguess=[-348.08538660805169,        39.219752218712493,       -1.4923716926995252,        1.6675727027061004E-002]
    lambguess = [-20.449156136430453,       -10.860126886269576 ,       1.0326425621279325,       -2.5382727424015664E-002]
    call reconstruct_maxent(mu,x,f,lambsol,ierr,use_log=.true.,lambguess=lambguess)
    if (ierr /= 1) print*,' INFO: '//fsolve_error(ierr)

    print*,' got lambsol = ',lambsol
    f = integrand(x, lambsol, 0)

    do k=0,size(mu)-1
       ftmp = f*x**(k/3.)
       kgot(k) = integrate_trap_log(ngrid,x,ftmp)
    enddo
    print*,' MOMENTS OUT : ',kgot
    call print_moments(kgot,a0,rhoi)

     ! print results
    write(filename,"(a,i5.5,a)") 'recon',i,'.out'
    open(unit=1,file=filename,status='replace',action='write')
    write(1,"(a)") '# N,f(N)'
    do j=1,ngrid
       write(1,*) x(j),f(j)
    enddo
    close(1)

    read*
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
