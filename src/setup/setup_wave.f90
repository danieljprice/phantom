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
!  DEPENDENCIES: boundary, dim, io, kernel, mpiutils, options, part,
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
 use part,         only:labeltype,set_particle_type,igas,idust,dustfrac
 use physcon,      only:pi
 use kernel,       only:radkern
 use dim,          only:maxvxyzu,use_dust,maxp
 use options,      only:use_dustfrac
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
!
! default options
!
 npartx  = 64
 ntypes  = 1
 rhozero = 1.
 massfac = 1.
 cs      = 1.
 ampl    = 1.d-4
 use_dustfrac = .false.
 if (id==master) then
    itype = 1
    print "(/,a,/)",'  >>> Setting up particles for linear wave test <<<'
    call prompt(' enter number of '//trim(labeltype(itype))//' particles in x ',npartx,8,int(maxp/144.))
    if (use_dust) then
       dust_method = 2
       dtg = 1.
       idrag = 2
       call prompt('Which dust method do you want? (1=one fluid,2=two fluid)',dust_method,1,2)
       if (dust_method == 1) then
          use_dustfrac = .true.
       else
          use_dustfrac = .false.
       endif
       if (use_dustfrac) K_code = 1000. ! for a more sensible better option
       call prompt('Enter dust to gas ratio',dtg,0.)
       call prompt('Enter constant drag coefficient',K_code,0.)
       if (use_dustfrac) then
          massfac = 1. + dtg
       else
          ntypes  = 2
       endif
    endif
 endif
 call bcast_mpi(npartx)
!
! boundaries
!
 xmini = -0.5
 xmaxi = 0.5
 length = xmaxi - xmini
 deltax = length/npartx
 ! try to give y boundary that is a multiple of 6 particle spacings in the low density part
 fac = 6.*(int((1.-epsilon(0.))*radkern/6.) + 1)
 deltay = fac*deltax*sqrt(0.75)
 deltaz = fac*deltax*sqrt(6.)/3.
 call set_boundary(xmin,xmax,-deltay,deltay,-deltaz,deltaz)
!
! general parameters
!
 time   = 0.
 kwave  = 2.*pi/length
 denom = length - ampl/kwave*(cos(kwave*length)-1.0)
!
! setup particles in the closepacked lattice
!
 if (maxvxyzu >= 4) then
    gamma = 5./3.
 else
    gamma  = 1.
 endif

 npart = 0
 npart_total = 0
 npartoftype(:) = 0

 overtypes: do itypes=1,ntypes
    select case (itypes)
    case(1)
       itype = igas
       if (id==master) call prompt('enter '//trim(labeltype(itype))//&
                            ' density (gives particle mass)',rhozero,0.)

       if (id==master) call prompt('enter sound speed in code units (sets polyk)',cs,0.)
       if (maxvxyzu < 4) then
          call bcast_mpi(cs)
          polyk = cs**2
          print*,' polyk = ',polyk
       else
          polyk = 0.
       endif
       if (id==master) call prompt('enter perturbation amplitude',ampl)
       call bcast_mpi(ampl)
       if (use_dustfrac) then
          call bcast_mpi(dtg)
       endif
    case(2)
       itype = idust
       rhozero = dtg*rhozero
    end select

    call bcast_mpi(rhozero)

    npart_previous = npart

    if (itype == igas) then
       call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax, &
                         hfact,npart,xyzh,nptot=npart_total,rhofunc=rhofunc)
       xmin_dust = xmin + dust_shift*deltax
       xmax_dust = xmax + dust_shift*deltax
       ymin_dust = ymin + dust_shift*deltax
       ymax_dust = ymax + dust_shift*deltax
       zmin_dust = zmin + dust_shift*deltax
       zmax_dust = zmax + dust_shift*deltax
    else
       call set_unifdis('closepacked',id,master,xmin_dust,xmax_dust,ymin_dust, &
                         ymax_dust,zmin_dust,zmax_dust,deltax, &
                         hfact,npart,xyzh,nptot=npart_total,rhofunc=rhofunc)
    endif

    !--set which type of particle it is
    do i=npart_previous+1,npart
       call set_particle_type(i,itype)

       vxyzu(1,i) = ampl*sin(kwave*(xyzh(1,i)-xmin))
       vxyzu(2:3,i) = 0.
!
!--perturb internal energy if not using a polytropic equation of state
!  (do this before density is perturbed)
!
       if (maxvxyzu >= 4) then
          if (gamma > 1.) then
             uuzero = cs**2/(gamma*(gamma-1.))
             przero = (gamma-1.)*rhozero*uuzero
          else
             uuzero = 3./2.*cs**2
             przero = cs**2*rhozero
          endif
          vxyzu(4,i) = uuzero + przero/rhozero*ampl*sin(kwave*(xyzh(1,i)-xmin))
       endif
!
!--one fluid dust: set dust fraction on gas particles
!
       if (use_dustfrac) then
          if (itype==igas) then
             call set_dustfrac(dtg,dustfrac(:,i))
          else
             dustfrac(:,i) = 0.
          endif
       endif
    enddo

    npartoftype(itype) = npart - npart_previous
    if (id==master) print*,' npart = ',npart,npart_total

    totmass = massfac*rhozero*dxbound*dybound*dzbound
    if (id==master) print*,' box volume = ',dxbound*dybound*dzbound,' rhozero = ',rhozero

    massoftype(itype) = totmass/npartoftype(itype)
    if (id==master) print*,' particle mass = ',massoftype(itype)

 enddo overtypes

contains

real function rhofunc(x)
 real, intent(in) :: x

 rhofunc = 1. + ampl*sin(kwave*(x - xmin))

end function rhofunc

end subroutine setpart

end module setup
