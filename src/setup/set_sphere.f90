!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: spherical
!
!  DESCRIPTION:
!   This module sets up spherical particle distributions
!   By default this is done by cropping and stretching cubes
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: random, stretchmap, unifdis
!+
!--------------------------------------------------------------------------
module spherical
 use unifdis, only:set_unifdis,mask_prototype,mask_true
 implicit none

 public  :: set_sphere

 real, parameter :: pi = 4.*atan(1.)

 integer, parameter :: &
   ierr_notinrange    = 1, &
   ierr_not_converged = 2, &
   ierr_unknown       = 3

 private

contains
!-----------------------------------------------------------------------
!+
!  This subroutine positions particles on a sphere - uniform
!  density by default, or if a function rhofunc(r) or a tabulated
!  array rhotab(:) is passed, with an arbitrary radial density profile rho(r)
!
!  The function is assumed to be a real function with a single argument:
!
!  real function rhofunc(r)
!   real, intent(in) :: r
!
!   rhofunc = 1./r**2
!
!  end function rhofunc
!
!  The table rhotab(:) is assumed to contain density values
!  on equally spaced radial bins between rmin and rmax
!+
!-----------------------------------------------------------------------
subroutine set_sphere(lattice,id,master,rmin,rmax,delta,hfact,np,xyzh, &
                      rhofunc,rhotab,rtab,xyz_origin,nptot,dir,exactN,np_requested,mask)
 use stretchmap, only:set_density_profile
 character(len=*), intent(in)    :: lattice
 integer,          intent(in)    :: id,master
 integer,          intent(inout) :: np
 real,             intent(in)    :: rmin,rmax,hfact
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(inout) :: delta
 real,             external,      optional :: rhofunc
 real,             intent(in),    optional :: rhotab(:), rtab(:)
 integer,          intent(in),    optional :: dir
 integer,          intent(in),    optional :: np_requested
 integer(kind=8),  intent(inout), optional :: nptot
 real,             intent(in),    optional :: xyz_origin(3)
 logical,          intent(in),    optional :: exactN
 procedure(mask_prototype), optional :: mask
 procedure(mask_prototype), pointer  :: my_mask
 integer,          parameter     :: maxits = 20
 real,             parameter     :: tol    = 1.e-9
 integer                         :: i,npin,icoord,ierr
 real                            :: xmin,xmax,ymin,ymax,zmin,zmax,vol_sphere
 logical                         :: use_sphereN

 xmin = -rmax
 xmax =  rmax
 ymin = -rmax
 ymax =  rmax
 zmin = -rmax
 zmax =  rmax
 npin =  np
 use_sphereN = .false.
 if (present(exactN) .and. present(np_requested)) then
    if ( exactN ) use_sphereN = .true.
 endif
 if (present(mask)) then
    my_mask => mask
 else
    my_mask => mask_true
 endif
 !
 !--Create a sphere of uniform density
 !
 if (lattice=='random' .and. present(np_requested)) then
    call set_sphere_mc(id,master,rmin,rmax,hfact,np_requested,np,xyzh,ierr,&
                       nptot,my_mask)
 elseif ( use_sphereN ) then
    vol_sphere = 4.0/3.0*pi*rmax**3
    call set_unifdis_sphereN(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,delta,&
                             hfact,np,np_requested,xyzh,rmax,vol_sphere,nptot,my_mask,ierr)
 else
    call set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax, &
                     zmin,zmax,delta,hfact,np,xyzh,.false.,&
                     rmin=rmin,rmax=rmax,nptot=nptot,verbose=.false.,centre=.true.)
 endif

 !
 ! allow stretch mapping to give arbitrary density profiles
 ! in any of the coordinate directions as specified by
 ! the input function rhofunc, or in a table rhotab
 !
 icoord = 1 ! default direction is radial
 if (present(dir)) then
    if (dir >= 1 .and. dir <= 3) icoord = dir
 endif
 if (present(rhofunc)) then
    call set_density_profile(np,xyzh,min=rmin,max=rmax,rhofunc=rhofunc,&
                             start=npin,geom=3,coord=icoord)
 elseif (present(rhotab) .and. present(rtab)) then
    call set_density_profile(np,xyzh,min=rmin,max=rmax,rhotab=rhotab,xtab=rtab,&
                             start=npin,geom=3,coord=icoord)
 elseif (present(rhotab)) then
    call set_density_profile(np,xyzh,min=rmin,max=rmax,rhotab=rhotab,start=npin,geom=3,coord=icoord)
 endif

 if (present(xyz_origin)) then
    if (id==master) write(*,"(1x,a,3(es10.3,1x))") 'shifting origin to ',xyz_origin(:)
    do i=npin+1,np
       ! shift positions and velocities to specified origin
       xyzh(1,i) = xyzh(1,i) + xyz_origin(1)
       xyzh(2,i) = xyzh(2,i) + xyz_origin(2)
       xyzh(3,i) = xyzh(3,i) + xyz_origin(3)
    enddo
 endif

end subroutine set_sphere

subroutine set_sphere_mc(id,master,rmin,rmax,hfact,np_requested,np,xyzh, &
                         ierr,nptot,mask)
!
! set up uniform spherical particle distribution using Monte Carlo particle
! placement. Particles are placed in pairs so the distribution
! is symmetric about the origin
!
 use random,     only:ran2
 use stretchmap, only:set_density_profile
 integer,          intent(in)    :: id,master,np_requested
 integer,          intent(inout) :: np   ! number of actual particles
 real,             intent(in)    :: rmin,rmax,hfact
 real,             intent(out)   :: xyzh(:,:)
 integer,          intent(out)   :: ierr
 integer(kind=8),  intent(inout), optional :: nptot
 procedure(mask_prototype) :: mask
 integer :: i,npin,iseed,maxp
 real    :: vol_sphere,rr,phi,theta,mr,dir(3)
 real    :: sintheta,costheta,sinphi,cosphi,psep
 integer(kind=8) :: iparttot

 npin = np
 iparttot = npin
 vol_sphere = 4./3.*pi*(rmax**3 - rmin**3)
 ! use mean particle spacing to set initial smoothing lengths
 psep = (vol_sphere/real(np_requested))**(1./3.)
 iseed = -1978
 maxp  = size(xyzh(1,:))
 ierr  = 1

 do i=npin+1,npin+np_requested,2
    !
    ! get random mass coordinate i.e. m(r)
    !
    mr = ran2(iseed)
    !
    ! invert to get mass coordinate from radial coordinate, i.e. r(m)
    !
    rr = mr**(1./3.)*rmax ! uniform density
    !
    ! get a random position on sphere
    !
    phi = 2.*pi*(ran2(iseed) - 0.5)
    costheta = 2.*ran2(iseed) - 1.
    theta    = acos(costheta)
    sintheta = sin(theta)
    sinphi   = sin(phi)
    cosphi   = cos(phi)
    dir  = (/sintheta*cosphi,sintheta*sinphi,costheta/)
    !
    ! add TWO particles, symmetric around the origin
    !
    iparttot = iparttot + 1
    if (mask(iparttot)) then
       np = np + 1
       if (np > maxp) then
          print*,' ERROR: np > maxp'
          return
       endif
       xyzh(1:3,np) = rr*dir
       xyzh(4,np)   = hfact*psep
    endif
    iparttot = iparttot + 1
    if (mask(iparttot)) then
       np = np + 1
       if (np > maxp) then
          print*,' ERROR: np > maxp'
          return
       endif
       xyzh(1:3,np) = -rr*dir
       xyzh(4,np)   = hfact*psep
    endif
 enddo
 ierr = 0
 if (present(nptot)) nptot = iparttot
 if (id==master) write(*,"(1x,a,i10,a)") 'placed ',np-npin,' particles in random-but-symmetric sphere'

end subroutine set_sphere_mc

!-----------------------------------------------------------------------
!+
!  If setting up a uniform sphere, this will iterate to get
!  approximately the desired number of particles
!  Note: Since this is only for a sphere, the call to this is the
!        same as if a sphere were being created by set_unifdis
!+
!-----------------------------------------------------------------------
subroutine set_unifdis_sphereN(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                    hfact,npart,nps_requested,xyzh,r_sphere,v_sphere,npart_total,mask,ierr)
 character(len=*), intent(in)    :: lattice
 integer,          intent(in)    :: id,master
 integer,          intent(inout) :: npart
 integer,          intent(in)    :: nps_requested
 real,             intent(in)    :: xmin,xmax,ymin,ymax,zmin,zmax,hfact,r_sphere,v_sphere
 real,             intent(out)   :: psep,xyzh(:,:)
 integer(kind=8),  intent(inout) :: npart_total
 integer,          intent(out)   :: ierr
 procedure(mask_prototype)       :: mask
 integer(kind=8)                 :: npart_local
 integer                         :: nps_lo,nps_hi,npr_lo,npr_hi,test_region,iter,nold(4)
 integer                         :: npin,npmax,npart0,nx,np,dn
 logical                         :: iterate_to_get_nps
 !
 !--Initialise values
 !
 ierr          = 0
 test_region   = 10
 npmax         = size(xyzh(1,:))
 nps_lo        = 0
 nps_hi        = npmax
 npr_lo        = 0
 npr_hi        = npmax
 np            = nps_requested
 dn            = nps_requested/20
 iter          = 0
 npin          = npart
 npart_local   = npart_total
 npart0        = 0
 nold          = 0
 iterate_to_get_nps = .true.

 print*, ' set_sphere: Iterating to form sphere with approx ',nps_requested,' particles'
 !
 !--Perform the iterations
 !
 do while (iterate_to_get_nps .and. iter < 100)
    iter = iter + 1
    nx   = int(np**(1./3.)) - 1           ! subtract 1 because of adjustment due to periodic BCs
    psep = (v_sphere)**(1./3.)/real(nx) ! particle separation in sphere
    if (lattice=='closepacked') psep = psep*sqrt(2.)**(1./3.)         ! adjust psep for close-packed lattice
    npart       = npin
    npart_total = npart_local
    call set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                    hfact,npart,xyzh,.false.,rmax=r_sphere,nptot=npart_total,verbose=.false.,mask=mask)
    if (nold(1)==np .and. nold(2)==npart0 .and. nold(3)==nps_lo .and. nold(4)==nps_hi) iterate_to_get_nps = .false.
    if (nps_lo > 0 .and. nps_hi < npmax) then
       nold(1) = np
       nold(2) = npart0
       nold(3) = nps_lo
       nold(4) = nps_hi
    endif
    npart0 = npart - npin
    if (npart0==nps_requested) then
       iterate_to_get_nps = .false.
    elseif (npart0 < nps_requested .and. nps_hi==npmax) then
       ! initialising for the case where npart0 is too small
       nps_lo = npart0
       npr_lo = np
       np     = np + dn
    elseif (npart0 > nps_requested .and. nps_lo==0) then
       ! initialising for the case where npart0 is too large
       nps_hi = npart0
       npr_hi = np
       np     = np - dn
    else
       if (nps_lo==0 .or. nps_hi==npmax) then
          ! finalise the boundaries, and begin testing upwards
          if (nps_lo == 0    ) then
             nps_lo = npart0
             npr_lo = np
          elseif (nps_hi == npmax) then
             nps_hi = npart0
             npr_hi = np
          endif
          np = npr_lo + (npr_hi - npr_lo)/3
          test_region = -1
       elseif (test_region == -1) then
          if (npart0 <= nps_lo .or. npart0 > nps_requested) then
             test_region = 1
             np     = npr_lo +  2*(npr_hi - npr_lo)/3
          else
             nps_lo = npart0
             npr_lo = np
             np     = npr_lo + (npr_hi - npr_lo)/3
          endif
       elseif (test_region == 1) then
          if (npart0 < nps_hi .and. npart0 > nps_requested ) then
             nps_hi = npart0
             npr_hi = np
             np     = npr_lo + 2*(npr_hi - npr_lo)/3
          endif
       else
          print "(a)",' ERROR: set_sphere: This option should not be possible'
          ierr = ierr_unknown
          return
       endif
    endif
 enddo
 if (iter >= 100 .or. nps_lo > nps_requested .or. nps_requested > nps_hi) then
    print "(a)",' ERROR: set_sphere: Failed to converge to the correct number of particles in the sphere'
    ierr = ierr_not_converged
 else
    if (nps_lo > nps_requested .or. nps_requested > nps_hi) then ! sanity check
       print "(a)",' ERROR: set_sphere: Did not converge for number of particles in the sphere'
       ierr = ierr_notinrange
       return
    endif
    ! always use more particles than requested
    npart = nps_hi
    if (nps_requested - nps_lo < nps_hi - nps_requested) then
       write(*,'(a,I8,a,F5.2,a)') " set_sphere: The closest number of particles to the requested number is " &
                    ,nps_lo,", which is ",float(nps_requested-nps_lo)/float(nps_requested)*100.0 &
                    ,"% than less requested."
       write(*,'(a)') " set_sphere: We will not use fewer than the requested number of particles."
    endif
    write(*,'(a,I8,a,F5.2,a)') " set_sphere: Using " &
              , nps_hi," particles, which is ",float(nps_hi - nps_requested)/float(nps_requested)*100.0 &
              ,"% more than requested."
  endif
  write(*,'(a,I10,a)') ' set_sphere: Iterations complete: added ',npart,' particles in sphere'

end subroutine set_unifdis_sphereN

end module spherical
