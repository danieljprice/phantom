!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: testeos
!
!  DESCRIPTION:
!  Unit tests of the equation of state module
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, io, physcon, testutils, units
!+
!--------------------------------------------------------------------------
module testeos
 implicit none
 public :: test_eos

 private

contains

subroutine test_eos(ntests,npass)
 use eos,       only:maxeos,equationofstate,rhocrit1cgs,polyk,polyk2,eosinfo,init_eos,isink
 use io,        only:id,master,stdout
 use testutils, only:checkval,checkvalbuf,checkvalbuf_start,checkvalbuf_end
 use physcon,   only:pi,solarm,pc
 use units,     only:unit_density,set_units
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(maxeos),ncheck(maxeos)
 integer :: i,maxpts,ierr,ieos
 real    :: rhoi,xi,yi,zi,ponrhoi,spsoundi,ponrhoprev,spsoundprev
 real    :: errmax
 character(len=20) :: pdir
 logical :: got_phantom_dir

 if (id==master) write(*,"(/,a,/)") '--> TESTING EQUATION OF STATE MODULE'

 call set_units(mass=solarm,dist=1.d16,G=1.d0)

 if (id==master) write(*,"(/,a)") '--> testing equation of state initialisation'
 nfailed = 0
 rhoi = 1.e-6*rhocrit1cgs/unit_density
 polyk = 0.1
 polyk2 = 0.1
 isink  = 1
 ntests = ntests + 1
 call get_environment_variable('PHANTOM_DIR',pdir) ! MESA EOS cannot initialise if this not set
 got_phantom_dir = (len_trim(pdir) > 0)
 do ieos=1,maxeos
    call init_eos(ieos,ierr)
    if (ieos==10 .and. .not. got_phantom_dir) cycle
    call checkval(ierr,0,0,nfailed(ieos),'eos initialisation')
 enddo
 if (all(nfailed==0)) npass = npass + 1

 if (id==master) write(*,"(/,a)") '--> testing barotropic equation of state'
 nfailed = 0
 call eosinfo(8,stdout)
 call checkvalbuf_start('equation of state is continuous')
 ncheck(:) = 0
 maxpts = 5000
 errmax = 0.
 do i=1,maxpts
    rhoi = 1.01*rhoi
    call equationofstate(8,ponrhoi,spsoundi,rhoi,xi,yi,zi)
    write(1,*) rhoi*unit_density,ponrhoi,ponrhoi*rhoi,spsoundi
    if (i > 1) call checkvalbuf(ponrhoi,ponrhoprev,1.e-2,'p/rho is continuous',nfailed(1),ncheck(1),errmax)
    !if (i > 1) call checkvalbuf(spsoundi,spsoundprev,1.e-2,'cs is continuous',nfailed(2),ncheck(2),errmax)
    ponrhoprev = ponrhoi
    spsoundprev = spsoundi
 enddo
 call checkvalbuf_end('p/rho is continuous',ncheck(1),nfailed(1),0,0,maxpts)
 !call checkvalbuf_end('cs is continuous',ncheck(2),nfailed(2),0,0,maxpts)
 ntests = ntests + 1
 if (nfailed(1)==0) npass = npass + 1

 if (id==master) write(*,"(/,a)") '<-- EQUATION OF STATE TEST COMPLETE'

end subroutine test_eos

end module testeos
