!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: slab
!
!  DESCRIPTION:
!   This module sets up particles in a thin slab, i.e. 3D box but with
!   small aspect ratio in the z direction. Useful for performing 2D
!   test problems in 3D
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, unifdis
!+
!--------------------------------------------------------------------------
module slab
 implicit none

contains
!----------------------------------------------------------------
!+
!  setup particles in thin slab geometry, mainly useful
!  for performing 2D test problems in 3D
!+
!----------------------------------------------------------------
subroutine set_slab(id,master,nx,xmini,xmaxi,ymini,ymaxi,deltax,hfact,np,xyzh,lattice)
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound
 use unifdis,      only:set_unifdis
 use domain,       only:i_belong
 integer,          intent(in)    :: id,master,nx
 integer,          intent(inout) :: np
 real,             intent(in)    :: xmini,xmaxi,ymini,ymaxi,hfact
 real,             intent(out)   :: xyzh(:,:),deltax
 character(len=*), intent(in), optional    :: lattice
 real :: dz
 character(len=20) :: mylattice
!
! use close packed lattice by default
!
 if (present(lattice)) then
    mylattice = lattice
 else
    mylattice = 'closepacked'
 endif
!
! set z direction boundary to achieve nx x nx x 12 particles by default
! the position of the z boundary depends on the lattice choice
!
 select case(mylattice)
 case('closepacked')
    dz = 4.*sqrt(6.)/nx
 case default
    deltax = (xmaxi-xmini)/nx
    dz = 6.*deltax
 end select
 call set_boundary(xmini,xmaxi,ymini,ymaxi,-dz,dz)
 deltax = dxbound/nx
!
! set particle lattice
!
 call set_unifdis(mylattice,id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                  hfact,np,xyzh,.true.,mask=i_belong)

end subroutine set_slab

end module slab
