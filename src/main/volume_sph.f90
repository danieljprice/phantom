!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module volume_sph
!
! Generalized volumes X_i=m_i/R_i, R_i=sum_j m_j W_ij(h_i), rho_i=R_i C_i,
! C_i=sum_j X_j W_ij(h_i). The h constraint uses a polynomial neighbour count.
! Pressure coefficients differentiate the complete discrete thermal energy, including
! the dependence of neighbouring X_j on R_j and h_j.
!
! :References: Garcia-Senz, Cabezon & Escartin (2022), A&A 659, A175
!   (volume estimator; the central variational forces here are not ISPH).
!
! :Owner: Francesco Radica
!
! :Runtime parameters: None
!
! :Dependencies: boundary, derivutils, dim, kdtree, kernel, neighkdtree, part
!
 use part, only:rho,particle_mass,gradh,dhdn,luminosity
 use dim, only:periodic,igradomega,track_lum
 use kernel, only:wkern,grkern,cnormk,radkern2
 use kdtree, only:inodeparts,inoderange
 use neighkdtree, only:active_leaves,nactive_leaves,get_neighbour_list,listneigh
 use derivutils, only:get_count_kernel,cnorm_count
 implicit none
 private
 real, allocatable, public :: rho0(:),pressure_a(:),pressure_b(:)
 real, allocatable, public :: pressure_count(:),hdot_count(:)
 real, allocatable :: drho0dt(:)
 public :: volume_density,volume_pressure,volume_work

contains

subroutine volume_density(npart,xyzh)
 integer, intent(in) :: npart
 real, intent(in) :: xyzh(:,:)
 integer :: ia,icell,ip,i,j,n,nneigh
 real :: cache(1,0),dr(3),r2,w,dwdh,c

 if (allocated(rho0)) then
    if (size(rho0)<npart) deallocate(rho0,pressure_a,pressure_b,drho0dt,pressure_count,hdot_count)
 endif
 if (.not.allocated(rho0)) allocate(rho0(npart),pressure_a(npart),pressure_b(npart),drho0dt(npart), &
                                   pressure_count(npart),hdot_count(npart))
 rho0(1:npart) = rho(1:npart)
 drho0dt(1:npart) = 0.

 !$omp parallel do default(none) schedule(runtime) &
 !$omp shared(nactive_leaves,active_leaves,inoderange,inodeparts,xyzh,rho0,rho,particle_mass) &
 !$omp private(ia,icell,ip,i,j,n,nneigh,cache,dr,r2,w,dwdh,c)
 do ia=1,nactive_leaves
    icell = active_leaves(ia)
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,cache,0)
    do ip=inoderange(1,icell),inoderange(2,icell)
       i = abs(inodeparts(ip))
       if (xyzh(4,i)<=0.) cycle
       c = 0.
       do n=1,nneigh
          j = listneigh(n)
          if (xyzh(4,j)<=0.) cycle
          dr = separation(xyzh(1:3,i),xyzh(1:3,j))
          r2 = dot_product(dr,dr)
          call volume_kernel(r2,xyzh(4,i),w,dwdh)
          c = c + particle_mass(j)/rho0(j)*w
       enddo
       rho(i) = rho0(i)*c
    enddo
 enddo
 !$omp end parallel do
end subroutine volume_density

subroutine volume_pressure(xyzh,vxyzu,pressure)
 real, intent(in) :: xyzh(:,:),vxyzu(:,:),pressure(:)
 integer :: ia,icell,ip,i,j,n,nneigh
 real :: cache(1,0),dr(3),r2,w,dwdh,wt,dwt,adjoint,ch,c,a,b,rh,rdot,ndot,hn,omega,q,g,gc,wc,projv

 !$omp parallel do default(none) schedule(runtime) &
 !$omp shared(nactive_leaves,active_leaves,inoderange,inodeparts,xyzh,pressure,rho0,rho,particle_mass) &
 !$omp shared(gradh,pressure_a,pressure_b,pressure_count,hdot_count,drho0dt,vxyzu) &
 !$omp private(ia,icell,ip,i,j,n,nneigh,cache,dr,r2,w,dwdh,wt,dwt,adjoint,ch,c,a,b,rh,rdot,ndot,hn,omega,q,g,gc,wc,projv)
 do ia=1,nactive_leaves
    icell = active_leaves(ia)
    ! The adjoint uses W_ji(h_j), so the union of both supports is required.
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,cache,0,getj=.true.)
    do ip=inoderange(1,icell),inoderange(2,icell)
       i = abs(inodeparts(ip))
       if (xyzh(4,i)<=0.) cycle
       adjoint = 0.
       ch = 0.
       rh = 0.
       rdot = 0.
       ndot = 0.
       do n=1,nneigh
          j = listneigh(n)
          if (xyzh(4,j)<=0.) cycle
          dr = separation(xyzh(1:3,i),xyzh(1:3,j))
          r2 = dot_product(dr,dr)
          call volume_kernel(r2,xyzh(4,i),w,dwdh)
          call volume_kernel(r2,xyzh(4,j),wt,dwt)
          ch = ch + particle_mass(j)/rho0(j)*dwdh
          rh = rh + particle_mass(j)*dwdh
          adjoint = adjoint + particle_mass(j)*pressure(j)*rho0(j)/rho(j)**2*wt
          if (r2>tiny(r2) .and. r2<radkern2*xyzh(4,i)**2) then
             q = sqrt(r2)/xyzh(4,i)
             projv = dot_product(dr,vxyzu(1:3,i)-vxyzu(1:3,j))/sqrt(r2)
             g = grkern(q*q,q)*cnormk/xyzh(4,i)**4
             call get_count_kernel(q,wc,gc)
             gc = gc*cnorm_count/xyzh(4,i)**4
             rdot = rdot + particle_mass(j)*g*projv
             ndot = ndot + gc*projv
          endif
       enddo
       c = rho(i)/rho0(i)
       a = pressure(i)/(rho(i)*rho0(i)) - adjoint/rho0(i)**2
       b = pressure(i)/(rho(i)*c)
       hn = dhdn(xyzh(4,i))
       omega = real(gradh(igradomega,i))
       pressure_a(i) = a
       pressure_b(i) = b
       pressure_count(i) = (a*rh+b*ch)*hn*omega
       hdot_count(i) = hn*omega*ndot
       drho0dt(i) = rdot + rh*hdot_count(i)
    enddo
 enddo
 !$omp end parallel do
end subroutine volume_pressure

subroutine volume_work(xyzh,vxyzu,pressure,fxyzu,do_work,thermal_fraction,dtlimit)
 real, intent(in) :: xyzh(:,:),vxyzu(:,:),pressure(:),thermal_fraction
 real, intent(inout) :: fxyzu(:,:),dtlimit
 logical, intent(in) :: do_work
 integer :: ia,icell,ip,i,j,n,nneigh
 real :: cache(1,0),dr(3),r2,w,dwdh,q,g,c,cdot,hdot,rhodot,limit,pdv

 limit = dtlimit
 !$omp parallel do default(none) schedule(runtime) reduction(min:limit) &
 !$omp shared(nactive_leaves,active_leaves,inoderange,inodeparts,xyzh,vxyzu,pressure,rho0,rho,particle_mass) &
 !$omp shared(drho0dt,hdot_count,fxyzu,do_work,thermal_fraction,track_lum,luminosity) &
 !$omp private(ia,icell,ip,i,j,n,nneigh,cache,dr,r2,w,dwdh,q,g,c,cdot,hdot,rhodot,pdv)
 do ia=1,nactive_leaves
    icell = active_leaves(ia)
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,cache,0)
    do ip=inoderange(1,icell),inoderange(2,icell)
       i = abs(inodeparts(ip))
       if (xyzh(4,i)<=0.) cycle
       if (do_work .or. track_lum) then
          cdot = 0.
          hdot = hdot_count(i)
          do n=1,nneigh
             j = listneigh(n)
             if (xyzh(4,j)<=0.) cycle
             dr = separation(xyzh(1:3,i),xyzh(1:3,j))
             r2 = dot_product(dr,dr)
             call volume_kernel(r2,xyzh(4,i),w,dwdh)
             g = 0.
             if (r2>tiny(r2) .and. r2<radkern2*xyzh(4,i)**2) then
                q = sqrt(r2)/xyzh(4,i)
                g = grkern(q*q,q)*cnormk/xyzh(4,i)**4 * &
                    dot_product(dr,vxyzu(1:3,i)-vxyzu(1:3,j))/sqrt(r2)
             endif
             cdot = cdot + particle_mass(j)/rho0(j)*(g + dwdh*hdot - w*drho0dt(j)/rho0(j))
          enddo
          c = rho(i)/rho0(i)
          rhodot = c*drho0dt(i) + rho0(i)*cdot
          pdv = pressure(i)/rho(i)**2*rhodot
          if (do_work) fxyzu(4,i) = fxyzu(4,i) + pdv
          if (track_lum) luminosity(i) = luminosity(i) + real(particle_mass(i)*pdv,kind=kind(luminosity))
       endif
       if (size(vxyzu,1)>=4) then
          if (vxyzu(4,i)>0. .and. fxyzu(4,i)<-tiny(1.)) &
             limit = min(limit,-thermal_fraction*vxyzu(4,i)/fxyzu(4,i))
       endif
    enddo
 enddo
 !$omp end parallel do
 dtlimit = limit
end subroutine volume_work

pure function separation(xi,xj) result(dr)
 use boundary, only:dxbound,dybound,dzbound
 real, intent(in) :: xi(3),xj(3)
 real :: dr(3),box(3)
 dr = xi-xj
 if (periodic) then
    box = [dxbound,dybound,dzbound]
    dr = dr - box*anint(dr/box)
 endif
end function separation

pure subroutine volume_kernel(r2,h,w,dwdh)
 real, intent(in) :: r2,h
 real, intent(out) :: w,dwdh
 real :: q2,q
 w = 0.
 dwdh = 0.
 q2 = r2/h**2
 if (q2>=radkern2) return
 q = sqrt(q2)
 w = cnormk/h**3*wkern(q2,q)
 dwdh = -(3.*w + cnormk/h**3*q*grkern(q2,q))/h
end subroutine volume_kernel

end module volume_sph
