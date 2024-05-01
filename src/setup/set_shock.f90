!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setshock
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: stretchmap, unifdis
!
 implicit none

 public :: set_shock, adjust_shock_boundaries, fsmooth

 private

contains
!-----------------------------------------------------------------------
!+
!  utility routine to set up a 3D particle distribution to represent
!  a discontinuous density jump located at x=xshock, with:
!
!    density = rholeft  for x <= xshock
!    density = rhoright for x > xshock
!+
!-----------------------------------------------------------------------
subroutine set_shock(latticetype,id,master,itype,rholeft,rhoright,xmin,xmax,ymin,ymax,zmin,zmax,&
                     xshock,dxleft,hfact,smooth_fac,npart,xyzh,massoftype,iverbose,ierr,mask)
 use unifdis, only:set_unifdis,get_ny_nz_closepacked,is_closepacked,mask_prototype,mask_true
 use stretchmap, only:rho_func
 character(len=*), intent(in) :: latticetype
 integer, intent(in) :: id,master,itype,iverbose
 real,    intent(in) :: rholeft,rhoright,xshock,xmin,xmax,ymin,ymax,zmin,zmax
 real,    intent(in) :: dxleft,hfact,smooth_fac
 integer, intent(inout) :: npart
 real,    intent(out)   :: xyzh(:,:),massoftype(:)
 integer, intent(out)   :: ierr
 procedure(mask_prototype), optional :: mask
 procedure(mask_prototype), pointer  :: my_mask
 procedure(rho_func), pointer :: density_func
 integer :: npartold,ny,nz
 real :: totmass,volume,dxright
 real :: xminleft(3),xmaxleft(3),xminright(3),xmaxright(3)
 logical, parameter :: periodic=.true.
 !
 ! sanity check some of the input to avoid seg faults
 !
 ierr = 0
 if (itype <= 0 .or. itype > size(massoftype)) then
    ierr = 1
    return
 endif
 if (present(mask)) then
    my_mask => mask
 else
    my_mask => mask_true
 endif

 density_func => rhosmooth

 !
 ! set limits of the different domains
 !
 xminleft(:)  = (/xmin,ymin,zmin/)
 xmaxleft(:)  = (/xmax,ymax,zmax/)
 xminright(:) = (/xmin,ymin,zmin/)
 xmaxright(:) = (/xmax,ymax,zmax/)
 npartold = npart

 if (abs(rholeft-rhoright) > epsilon(0.)) then
    ! if there is a density jump, divide the x axis into two halves at xshock
    xmaxleft(1)  = xshock
    xminright(1) = xshock
    if (id==master) write(*,'(1x,3(a,es16.8))') 'shock: left half  ',xminleft(1), ' to ',xmaxleft(1), ' with dx_left  = ',dxleft
    dxright = dxleft*(rholeft/rhoright)**(1./3.) ! NB: dxright is corrected for closepacked

    ! set up a uniform lattice for left half
    if (smooth_fac > 0.) then
       call set_unifdis(latticetype,id,master,xminleft(1),xmaxleft(1),xminleft(2), &
            xmaxleft(2),xminleft(3),xmaxleft(3),dxleft,hfact,npart,xyzh,&
            periodic,rhofunc=density_func,mask=my_mask)
    else
       call set_unifdis(latticetype,id,master,xminleft(1),xmaxleft(1),xminleft(2), &
            xmaxleft(2),xminleft(3),xmaxleft(3),dxleft,hfact,npart,xyzh,&
            periodic,mask=my_mask)
    endif

    ! set particle mass
    volume            = product(xmaxleft-xminleft)
    totmass           = volume*rholeft
    massoftype(itype) = totmass/real(npart-npartold)
    if (id==master) write(*,'(1x,a,es16.8)') 'shock: particle mass = ',massoftype(itype)

    if (is_closepacked(latticetype)) then
       ! now adjust spacing on right hand side to get correct density given the particle mass
       volume  = product(xmaxright - xminright)
       totmass = volume*rhoright
       call get_ny_nz_closepacked(dxright,xminright(2),xmaxright(2),xminright(3),xmaxright(3),ny,nz)
       dxright = (xmaxright(1) - xminright(1))/((totmass/massoftype(itype))/(ny*nz))
    endif

    ! now set up box for right half
    if (id==master) write(*,'(1x,3(a,es16.8))') 'shock: right half ',xminright(1),' to ',xmaxright(1),' with dx_right = ',dxright

    if (smooth_fac > 0.) then
       call set_unifdis(latticetype,id,master,xminright(1),xmaxright(1), &
            xminright(2),xmaxright(2),xminright(3),xmaxright(3),dxright,hfact,&
            npart,xyzh,periodic,npy=ny,npz=nz,rhofunc=density_func,mask=my_mask) ! set right half
    else
       call set_unifdis(latticetype,id,master,xminright(1),xmaxright(1), &
            xminright(2),xmaxright(2),xminright(3),xmaxright(3),dxright,hfact,&
            npart,xyzh,periodic,npy=ny,npz=nz,mask=my_mask) ! set right half
    endif

 else  ! set all of volume if densities are equal
    write(*,'(3(a,es16.8))') 'shock: uniform density  ',xminleft(1), ' to ',xmaxright(1), ' with dx  = ',dxleft
    call set_unifdis(latticetype,id,master,xminleft(1),xmaxleft(1),xminleft(2), &
                     xmaxleft(2),xminleft(3),xmaxleft(3),dxleft,hfact,npart,xyzh,periodic,mask=my_mask)
    volume           = product(xmaxleft-xminleft)
    dxright          = dxleft
    totmass          = rholeft*volume
    massoftype(itype) = totmass/real(npart-npartold)
 endif

contains
!----------------------------------------------
!+
!  Callback function to setup smooth density
!+
!----------------------------------------------
real function rhosmooth(x)
 real, intent(in) :: x

 rhosmooth = fsmooth(x,xshock,min(dxleft,dxright),smooth_fac,rholeft,rhoright)

end function rhosmooth

end subroutine set_shock

!------------------------------------------------------------------------
!+
!  Function specifying smoothing of various quantities across the shock
!+
!------------------------------------------------------------------------
real function fsmooth(x,x0,psep,fac,fl,fr)
 real, intent(in) :: x,x0,psep,fac,fl,fr
 real :: delta,exx
 real, parameter :: dsmooth = 20.

 if (fac > 0.)  then
    delta = (x - x0)/(fac*psep)
 else
    delta = 2.*sign(1.0,x-x0)*dsmooth
 endif

 if (delta > dsmooth) then
    fsmooth = fr
 elseif (delta < -dsmooth) then
    fsmooth = fl
 else
    exx = exp(delta)
    fsmooth = (fl + fr*exx)/(1. + exx)
 endif

end function fsmooth

!-----------------------------------------------------------------------
!+
!  Adjust the shock boundaries to allow for inflow/outflow
!+
!-----------------------------------------------------------------------
subroutine adjust_shock_boundaries(dxleft,dxright,radkern,vxleft,vxright, &
                                   densleft,densright,tmax,ndim,&
                                   xmin,xmax,ymin,ymax,zmin,zmax,use_closepacked)
 real,    intent(in)    :: dxleft,radkern,vxleft,vxright,densleft,densright,tmax
 real,    intent(out)   :: dxright
 real,    intent(inout) :: xmin,xmax,ymin,ymax,zmin,zmax
 integer, intent(in)    :: ndim
 logical, intent(in)    :: use_closepacked
 real :: fac

 if (vxleft > tiny(vxleft)) then
    xmin = xmin - vxleft*tmax
 endif
 dxright = dxleft*(densleft/densright)**(1./ndim) ! NB: dxright here is only approximate
 if (vxright < -tiny(vxright)) then
    xmax = xmax - vxright*tmax
 endif
 ! try to give y boundary that is a multiple of 6 particle spacings in the low density part
 fac = -6.*(int(1.99*radkern/6.) + 1)*max(dxleft,dxright)
 if (use_closepacked) then
    ymin   = fac*sqrt(0.75)
    zmin   = fac*sqrt(6.)/3.
 else
    ymin   = fac
    zmin   = fac
 endif
 ymax  = -ymin
 zmax  = -zmin

end subroutine adjust_shock_boundaries

end module setshock
