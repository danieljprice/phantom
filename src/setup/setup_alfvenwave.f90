!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Setup for MHD wave tests in 3D
!
!  REFERENCES:
!    Toth G., 2000, J. Comp. Phys., 161, 605
!    Stone J. M., et al., 2008, ApJS, 178, 137
!    Gardiner T. A., Stone J. M., 2008, J. Comp. Phys., 227, 4123
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    gamma   -- adiabatic index
!    iselect --  which wave test to run
!    nx      -- resolution (number of particles in x) for -xleft < x < xshock
!    rotated --  rotate wave vector?
!
!  DEPENDENCIES: boundary, dim, geometry, infile_utils, io, mpiutils, part,
!    physcon, prompting, setup_params, timestep, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart
!
! runtime options and default settings
!
 integer :: iselect
 integer :: nx
 logical :: rotated
 real    :: wavelength, ampl

 integer, parameter :: maxwaves = 5
 character(len=*), parameter :: wavetype(0:maxwaves-1) = &
      (/'non-linear circularly polarised Alfven wave', &
        'linear fast wave                           ', &
        'linear Alfven wave                         ', &
        'linear slow wave                           ', &
        'linear entropy wave                        '/)

 private

contains

!----------------------------------------------------------------
!+
!  setup for MHD wave tests
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxp,maxvxyzu
 use setup_params, only:rhozero,ihavesetupB
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd
 use io,           only:master
 use prompting,    only:prompt
 use mpiutils,     only:bcast_mpi
 use physcon,      only:pi
 use geometry,     only:igeom_rotated,igeom_cartesian,set_rotation_angles,coord_transform
 use timestep,     only:tmax,dtmax
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(in)    :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: deltax,totmass
 integer :: i,ierr,igeom
 real :: przero,uuzero,Bvec(3),vvec(3),Bnew(3),Bzero(3),vzero(3)
 real :: sina,sinb,cosa,cosb,wk,sinx1,cosx1,rvec(8),uui
 real :: runit(3),gam1,x0(3),xdot0,rhoi,x1
 real :: drho,dv(3),dB(3),du
 real :: q0(8),q(8),vwave,denom,dxi
 character(len=len(fileprefix)+6) :: setupfile
!
!--general parameters
!
 time = 0.
 gamma = 5./3.
!
!--default settings
!
 rotated = .false.
 iselect = 0
 wavelength = 1.
 !
 ! read setup parameters from the .setup file.
 ! if file does not exist, then ask for user input
 !
 setupfile = trim(fileprefix)//'.setup'
 call read_setupfile(setupfile,gamma,ierr)
 if (ierr /= 0 .and. id==master) then
    call interactive_setup()
    call write_setupfile(setupfile,gamma)
 endif
 gamma = 5./3.
!
!--setup parameters
!
 if (rotated) then
    sina = 2./3.
    sinb = 2./sqrt(5.)
    cosa = sqrt(5.)/3.
    cosb = 1./sqrt(5.)
    igeom = igeom_rotated
    call set_rotation_angles(sin_a=sina,sin_b=sinb,cos_a=cosa,cos_b=cosb)
 else
    sina = 0.
    sinb = 0.
    cosa = 1.
    cosb = 1.
    igeom = igeom_cartesian
 endif
 runit = (/cosa*cosb,cosa*sinb,sina/)
 rhozero  = 1.
 rvec     = 0.
 du = 0.
 drho = 0.
 dv = 0.
 dB = 0.
 gam1 = gamma - 1.
 select case(iselect)
 case(1:4)
    ampl   = 1.e-4
    przero = 3./5. ! ! 1./gamma
    Bzero  = (/1.,1.5,0./)
    vzero  = 0.
    uuzero = przero/(gam1*rhozero)
    if (iselect==4) vzero(1) = 1.
    call get_eigenvector(iselect,rvec)
    call get_amplitudes(iselect,Bzero,sqrt(gamma*przero/rhozero),rhozero,uuzero,przero,ampl,drho,dv,dB,du,vwave)
 case default
    ampl   = 0.1
    przero = 0.1
    Bzero  = (/1.,0.,0./)
    vzero  = 0.
    uuzero = przero/(gam1*rhozero)
 end select

 call prim_to_cons(rhozero,vzero,Bzero,uuzero,q0)

 print "(/,a)",' MHD wave setup : '//trim(wavetype(iselect))
 print 10, sina,sinb,rhozero,przero,vwave,wavelength/vwave
10 format(/,' sin(alpha) = ',f6.3,',   sin(beta) = ',f6.3,/,&
            ' density = ',f6.3,',       pressure = ',f6.3,/, &
            '   vwave = ',f6.3,',         period = ',f6.3,/)

 call print_amplitudes(rhozero,drho,vzero,dv,Bzero,dB,uuzero,du)
 print*,' rhozero = ',rhozero,'dv = ',dv,' dB = ',dB
 print*,drho,rhozero*dv,dB,du,du + dot_product(vzero,dv) + dot_product(Bzero,dB)/rhozero

 if (maxvxyzu < 4) then
    polyk = przero/rhozero**gamma
 else
    polyk = 0.
 endif

 call bcast_mpi(nx)
!
!--boundaries
!
 call set_boundary(-1.5,1.5,-0.75,0.75,-0.75,0.75)
 deltax = dxbound/nx

 wk = 2.*pi/wavelength
 x0 = (/xmin,ymin,zmin/)
 xdot0 = dot_product(x0,runit)

 if (rvec(1) > 0.) then
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,&
                     rhofunc=rhofunc,geom=igeom)
 else
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh)
 endif
 npartoftype(:) = 0
 npartoftype(1) = npart

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart
 print*,'npart = ',npart,' particle mass = ',massoftype(1)


 do i=1,npart
    x1    = dot_product(xyzh(1:3,i),runit)
    sinx1 = sin(wk*(x1-xdot0))
    cosx1 = cos(wk*(x1-xdot0))

    !denom = (xmax-xmin) - ampl/wk*(cos(wk*(xmax-xmin))-1.0)
    !call set_perturbation(xyzh(1,i),xmin,(xmax-xmin),wk,ampl,denom,dxi)
    !xyzh(1,i) = xmin + dxi

    select case(iselect)
    case(1:4)
       vvec = vzero + dv*cosx1
       Bvec = Bzero + dB*cosx1
!       if (iselect==4) then
!          uui = przero/(gam1*rhofunc(x1))
!       else
       uui = uuzero + du*cosx1
!       endif

       !q = q0 + ampl*rvec*cosx1
       !call cons_to_prim(q,rhoi,vvec,Bvec,uui)
    case default
       vvec = vzero + ampl*(/0.,sinx1,cosx1/)
       bvec = Bzero + ampl*(/0.,sinx1,cosx1/)
       uui  = uuzero
    end select

    call transform_vec(vvec,vxyzu(1:3,i),sina,sinb,cosa,cosb)
    call transform_vec(Bvec,Bnew,sina,sinb,cosa,cosb)

    if (maxvxyzu >= 4) vxyzu(4,i) = uui
    if (mhd) Bxyz(1:3,i) = Bnew
 enddo

 if (mhd) ihavesetupB = .true.

 tmax = wavelength/vwave
 dtmax = 0.1*tmax

contains

real function rhofunc(xi)
 real, intent(in) :: xi

 rhofunc = rhozero + drho*cos(wk*(xi - xdot0))
 !rhofunc = rhozero + drho*sin(wk*(xi - xdot0))
 !rhofunc = q0(1) + ampl*rvec(1)*cos(wk*(xi - xdot0))

end function rhofunc

end subroutine setpart

!----------------------------------------------------------
subroutine set_perturbation(xi,xmin,length,kwave,ampl,denom,dxi)
 use io, only:fatal
 real, intent(in)  :: xi,xmin,length,kwave,ampl,denom
 real, intent(out) :: dxi
 integer, parameter :: itsmax = 20
 real, parameter    :: tol = 1.d-10
 integer :: its
 real    :: dxprev, xmassfrac, func, fderiv

 dxi = xi-xmin
 dxprev = length*2.
 xmassfrac = dxi/length

 its = 0
 do while ((abs(dxi-dxprev) > tol).and.(its < itsmax))
    dxprev = dxi
    func = xmassfrac*denom - (dxi - ampl/kwave*(cos(kwave*dxi)-1.0))
    fderiv = -1.0 - ampl*sin(kwave*dxi)
    dxi = dxi - func/fderiv  ! Newton-Raphson iteration
    its = its + 1
 enddo

 if (its >= itsmax) then
    print*,'Error: soundwave - too many iterations'
    call fatal('setup_dustywave','Error: soundwave - too many iterations')
 endif

end subroutine set_perturbation

subroutine get_amplitudes(iwave,B,cs,rho,u,p,amp,drho,dv,dB,du,vwave)
 integer, intent(in) :: iwave
 real, intent(in)  :: B(3),cs,rho,u,p,amp
 real, intent(out) :: drho,dv(3),dB(3),du,vwave
 real :: va2,terma,termb,vfast,vslow,valfven,term

 ! determine wave speeds
 va2   = dot_product(B,B)/rho
 terma = cs**2 + va2
 termb = terma**2 - 4.*(cs*B(1))**2/rho
 ! define MHD wave speeds
 vfast = sqrt(0.5*(terma + sqrt(termb)))
 vslow = sqrt(0.5*(terma - sqrt(termb)))
 valfven = sqrt(va2)
 print*,' vfast = ',vfast,' vslow = ',vslow,' valfven = ',valfven

 select case(iwave)
 case(1)
    vwave = vfast
 case(2)
    vwave = valfven
 case(3)
    vwave = vslow
 case(4)
    vwave = 1.
 case default
    vwave = valfven
 end select

 ! now give perturbation amplitudes for v and B
 select case(iwave)
 case(1,3)     ! fast and slow MHD waves
    term = vwave**2 - B(1)**2/rho
    drho  = amp
    dv(1) = vwave*amp
    dv(2) = -dv(1)*B(1)*B(2)/(rho*term)
    dv(3) = -dv(1)*B(1)*B(3)/(rho*term)
    dB(1) = 0.
    dB(2) = vwave*B(2)*dv(1)/term
    dB(3) = vwave*B(3)*dv(1)/term
    du    = p/rho**2*drho
 case(4)       ! entropy mode
    drho = amp
    dv(1) = vwave*amp
    dv(2:3) = 0.
    dB = 0.
    du = -u*drho/rho !0.5*amp !-u*drho/rho
 case default  ! Alfven waves
    drho  = 0.
    dv(1) = 0.
    dv(2) = amp
    dv(3) = 0.
    dB(1) = 0.
    dB(2) = amp
    dB(3) = 0.
    du    = 0.
 end select

end subroutine get_amplitudes

subroutine print_amplitudes(rho,drho,v,dv,B,dB,u,du)
 real, intent(in) :: rho,drho,v(3),dv(3),B(3),dB(3),u,du

 write(*,"('    |',8(a10,1x,'|'))") 'rho','v1','v2','v3','B1','B2','B3','u'
 write(*,"(' q0 |',8(1pg10.2,1x,'|'))") rho,v,B,u
 write(*,"(' dq |',8(1pg10.2,1x,'|'))") drho,dv,dB,du

end subroutine print_amplitudes

!------------------------------------------------
!+
!  Right eigenvectors for the various MHD waves
!  from the Appendix to Gardiner & Stone (2008)
!+
!------------------------------------------------
pure subroutine get_eigenvector(iwave,rvec)
 integer, intent(in)  :: iwave
 real,    intent(out) :: rvec(8)
 real :: dir

 dir = 1. ! or -1
 rvec = 0.
 select case(iwave)
 case(1)
    ! fast wave
    rvec = (/2., 4.*dir, -2.*dir, 0., 0., 4., 0., 9./)/(2.*sqrt(5.))
 case(2)
    ! Alfven wave
    rvec = (/0., 0., 0., -1.*dir, 0., 0., 1., 0./)
 case(3)
    ! slow wave
    rvec = (/4., 2.*dir, 4.*dir, 0., 0., -2., 0., 3./)/(2.*sqrt(5.))
 case(4)
    ! entropy wave
    rvec = 0.5*(/2., 2., 0., 0., 0., 0., 0., 1./)
 end select

end subroutine get_eigenvector

!------------------------------------------------
!+
!  primitive to conservative variable transform
!+
!------------------------------------------------
pure subroutine prim_to_cons(rho,v,B,u,q)
 real, intent(in)  :: rho,v(3),B(3),u
 real, intent(out) :: q(8)

 q = (/rho,rho*v,B,rho*(u + 0.5*dot_product(v,v))/)

end subroutine prim_to_cons

!------------------------------------------------
!+
!  conservative to primitive variable transform
!+
!------------------------------------------------
pure subroutine cons_to_prim(q,rho,v,B,u)
 real, intent(in)  :: q(8)
 real, intent(out) :: rho,v(3),B(3),u

 rho = q(1)
 v   = q(2:4)/rho
 B   = q(5:7)
 u   = q(8)/rho - 0.5*dot_product(v,v)

end subroutine cons_to_prim

!------------------------------------------------------
!+
!  transform vectors from rotated to unrotated coords
!+
!------------------------------------------------------
pure subroutine transform_vec(xvec,x,sina,sinb,cosa,cosb)
 real, intent(in)  :: xvec(3)
 real, intent(out) :: x(3)
 real, intent(in)  :: sina,sinb,cosa,cosb

 x(1) = xvec(1)*cosa*cosb - xvec(2)*sinb - xvec(3)*sina*cosb
 x(2) = xvec(1)*cosa*sinb + xvec(2)*cosb - xvec(3)*sina*sinb
 x(3) = xvec(1)*sina + xvec(3)*cosa

end subroutine transform_vec

!------------------------------------------------------
!+
!  transform vectors from unrotated to rotated coords
!+
!------------------------------------------------------
pure subroutine untransform_vec(x,xvec,sina,sinb,cosa,cosb)
 real, intent(in)  :: x(3)
 real, intent(out) :: xvec(3)
 real, intent(in)  :: sina,sinb,cosa,cosb

 xvec(1) = x(3)*sina + cosa*(x(1)*cosb + x(2)*sinb)
 xvec(2) = x(2)*cosb - x(1)*sinb
 xvec(3) = x(3)*cosa - sina*(x(1)*cosb + x(2)*sinb)

end subroutine untransform_vec

!------------------------------------------
!+
!  Prompt user for setup options
!+
!------------------------------------------
subroutine interactive_setup()
 use prompting, only:prompt
 integer :: i

 print "(5(/,i2,' : ',a))",(i,trim(wavetype(i)),i=0,maxwaves-1)
 call prompt('Select which problem to run ',iselect,0,maxwaves-1)

 nx = 128
 call prompt('Enter resolution (number of particles in x)',nx,8)

end subroutine interactive_setup

!------------------------------------------
!+
!  Write setup parameters to input file
!+
!------------------------------------------
subroutine write_setupfile(filename,gamma)
 use infile_utils, only:write_inopt
 use dim,          only:tagline,maxvxyzu
 character(len=*), intent(in) :: filename
 real,             intent(in) :: gamma
 integer,          parameter  :: lu = 20
 integer                      :: ierr1

 write(*,"(a)") ' Writing '//trim(filename)//' with setup info'
 open(unit=lu,file=filename,status='replace',form='formatted')
 write(lu,"(a)") '# '//trim(tagline)
 write(lu,"(a)") '# input file for Phantom MHD linear wave test setup'

 write(lu,"(/,a)") '# MHD wave tests'
 call write_inopt(iselect,'iselect',' which wave test to run',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing iselect'

 call write_inopt(rotated,'rotated',' rotate wave vector?',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing rotated'

 write(lu,"(/,a)") '# resolution'
 call write_inopt(nx,'nx','resolution (number of particles in x) for -xleft < x < xshock',lu,ierr1)
 if (ierr1 /= 0) write(*,*) 'ERROR writing nx'

 write(lu,"(/,a)") '# Equation-of-state properties'
 call write_inopt(gamma,'gamma','adiabatic index',lu,ierr1)

 close(unit=lu)

end subroutine write_setupfile

!------------------------------------------
!+
!  Read setup parameters from input file
!+
!------------------------------------------
subroutine read_setupfile(filename,gamma,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db,read_inopt
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 real,             intent(out) :: gamma
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*, '(1x,2a)') 'setup_wave: reading setup options from ',trim(filename)

 nerr = 0
 call read_inopt(nx,'nx',db,min=8,errcount=nerr)
 call read_inopt(iselect,'iselect',db,min=0,errcount=nerr)
 call read_inopt(rotated,'rotated',db,errcount=nerr)
 call read_inopt(gamma,'gamma',db,min=1.,errcount=nerr)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_wave: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile

end module setup

