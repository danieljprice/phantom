!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module geodesic
!
! geodesic
!
! :References: None
!
! :Owner: Fitz) Hu
!
! :Runtime parameters: None
!
! :Dependencies: cons2primsolver, eos, extern_gr, io, metric_tools, part,
!   timestep
!

implicit none

public :: integrate_geodesic

private

contains

subroutine integrate_geodesic(pmass,xyzh,vxyzu,dens,pr,gamma,temp,pxyzu,dist,time)
 use extern_gr,      only:get_grforce
 use metric_tools,   only:pack_metric,pack_metricderivs
 use eos,            only:ieos,equationofstate
 use cons2primsolver,only:conservative2primitive
 use io,             only:warning,fatal
 use timestep,       only:bignumber,xtol,ptol
 use part,           only:rhoh,ien_type
 real, intent(inout) :: xyzh(:),vxyzu(:),pxyzu(:)
 real, intent(inout) :: dens,pr,gamma,temp,pmass
 real, intent(in), optional :: dist,time
 real :: metrics(0:3,0:3,2),metricderivs(0:3,0:3,3),fext(3)
 real :: t,tend,v,dt,hdt
 integer, parameter :: itsmax = 50
 integer :: its,ierr
 real :: xyz(3),pxyz(3),eni,vxyz(1:3),uui,rho,spsoundi,pondensi
 real :: vxyz_star(3),xyz_prev(3),pprev(3),fstar(3)
 real :: pmom_err,x_err
 logical :: converged

 xyz       = xyzh(1:3)
 pxyz      = pxyzu(1:3)
 eni       = pxyzu(4)
 vxyz      = vxyzu(1:3)
 uui       = vxyzu(4)
 rho       = rhoh(xyzh(4),pmass)
 if (present(dist)) then
   v = sqrt(dot_product(vxyz,vxyz))
   tend = dist/v
 elseif (present(time)) then
   tend = time
 else
   tend = 0.
 endif

 call pack_metric(xyz,metrics(:,:,:))
 call pack_metricderivs(xyz,metricderivs(:,:,:))
 call get_grforce(xyzh(:),metrics(:,:,:),metricderivs(:,:,:),vxyz,dens,uui,pr,fext(1:3),dt)

 t = 0.
 dt = min(0.1,tend*0.1,dt)
 hdt = 0.5*dt
 
 do while (t <= tend)
    t    = t + dt
    pxyz = pxyz + hdt*fext


    ! Note: grforce needs derivatives of the metric,
    ! which do not change between pmom iterations
    its = 0
    converged = .false.
    pmom_iterations: do while (its <= itsmax .and. .not. converged)
       its   = its + 1
          
       pprev = pxyz
       call conservative2primitive(xyz,metrics(:,:,:),vxyz,dens,uui,pr,&
                                       temp,gamma,rho,pxyz,eni,ierr,ien_type)
       if (ierr > 0) call warning('cons2primsolver [in integrate_geodesic (a)]','did not converge')
       call get_grforce(xyzh(:),metrics(:,:,:),metricderivs(:,:,:),vxyz,dens,uui,pr,fstar)

       pxyz = pprev + hdt*(fstar - fext)
       pmom_err = maxval(abs(pxyz - pprev))
       if (pmom_err < ptol) converged = .true.
       fext = fstar
    enddo pmom_iterations
    if (its > itsmax ) call warning('integrate_geodesic',&
                              'max # of pmom iterations',var='pmom_err',val=pmom_err)

    call conservative2primitive(xyz,metrics(:,:,:),vxyz,dens,uui,pr,temp,&
                                    gamma,rho,pxyz,eni,ierr,ien_type)
    if (ierr > 0) call warning('cons2primsolver [in integrate_geodesic (b)]','did not converge')
    xyz = xyz + dt*vxyz
    call pack_metric(xyz,metrics(:,:,:))

    its        = 0
    converged  = .false.
    vxyz_star = vxyz

    xyz_iterations: do while (its <= itsmax .and. .not. converged)
       its         = its+1
       xyz_prev    = xyz
  
       call conservative2primitive(xyz,metrics(:,:,:),vxyz_star,dens,uui,&
                                       pr,temp,gamma,rho,pxyz,eni,ierr,ien_type)
       if (ierr > 0) call warning('cons2primsolver [in integrate_geodesic (c)]','did not converge')
       xyz  = xyz_prev + hdt*(vxyz_star - vxyz)
       x_err = maxval(abs(xyz-xyz_prev))
       if (x_err < xtol) converged = .true.
       vxyz = vxyz_star
       ! UPDATE METRIC HERE
       call pack_metric(xyz,metrics(:,:,:))
    enddo xyz_iterations
    call pack_metricderivs(xyz,metricderivs(:,:,:))
    if (its > itsmax ) call warning('integrate_geodesic','Reached max number of x iterations. x_err ',val=x_err)

    call equationofstate(ieos,pondensi,spsoundi,dens,xyzh(1),xyzh(2),xyzh(3),temp,uui)
    pr = pondensi*dens

    call get_grforce(xyzh(:),metrics(:,:,:),metricderivs(:,:,:),vxyzu(1:3),dens,vxyzu(4),pr,fext(1:3),dt)
    pxyzu(1:3) = pxyzu(1:3) + hdt*fext(1:3)
    call conservative2primitive(xyz,metrics(:,:,:),vxyz,dens,uui,&
                                    pr,temp,gamma,rho,pxyz,eni,ierr,ien_type)
    xyzh(1:3) = xyz(1:3)
    pxyzu(1:3) = pxyz(1:3)
    vxyzu(1:3) = vxyz(1:3)
    vxyzu(4) = uui
 enddo
end subroutine integrate_geodesic

end module geodesic
