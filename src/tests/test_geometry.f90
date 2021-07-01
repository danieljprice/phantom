!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module testgeometry
!
! Unit tests of the geometry module
!
! :References:
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: geometry, io, testutils
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_geometry

 private

contains
!--------------------------------------------
!+
!  Unit tests of the geometry module
!+
!--------------------------------------------
subroutine test_geometry(ntests,npass)
 use geometry, only:maxcoordsys,coord_transform,labelcoordsys,vector_transform,&
                    get_coord_limits,small_number,igeom_planetwake
 integer, intent(inout) :: ntests,npass
 real    :: xin(3),xout(3),xtmp(3),tol
 real    :: vecin(3),vecout(3),vectmp(3),xmin(3),xmax(3),rad
 integer :: igeom,ndim,ierr,ndiff(2),itest

 if (id==master) write(*,"(/,a)") '--> TESTING GEOMETRY MODULE'
!
!--check that forward followed by inverse transform
!  returns original result
!
 ndim = 3
 xin   = (/5.,6.,7./)
 vecin = (/2.,3.,4./)

 do itest=1,2
    select case(itest)
    case(2)
       xin = (/0.,0.,0./)
       if (id==master) write(*,"(/,a)") '--> checking transforms for coords and vectors are reversible at r=0'
    case default
       if (id==master) write(*,"(/,a)") '--> checking transforms for coords and vectors are reversible'
    end select
    ndiff = 0
    do igeom=1,maxcoordsys
       ! forward
       call coord_transform(xin,ndim,1,xout,ndim,igeom,ierr)
       call vector_transform(xin,vecin,ndim,1,vecout,ndim,igeom,ierr)

       ! reverse
       call coord_transform(xout,ndim,igeom,xtmp,ndim,1,ierr)
       call vector_transform(xout,vecout,ndim,igeom,vectmp,ndim,1,ierr)

       tol = small_number
       if (igeom==igeom_planetwake) tol = 1e-13
       call checkval(3,xtmp,xin,tol,ndiff(1),trim(labelcoordsys(igeom)))
       call checkval(3,vectmp,vecin,tol,ndiff(2),trim(labelcoordsys(igeom)))
    enddo
    call update_test_scores(ntests,ndiff,npass)
 enddo
!
! check that routines that allegedly do the same thing actually do
!
 if (id==master) write(*,"(/,a)") '--> testing get_coord_limits matches coord_transform'
 ndiff = 0
 do igeom=1,maxcoordsys
    xmin = 0.
    xmax = 1.
    rad = 0.1
    call get_coord_limits(rad,xin,xout,xmin,xmax,igeom)
    call coord_transform(xin,ndim,1,xtmp,ndim,igeom,ierr)
    call checkval(3,xout,xtmp,tol,ndiff(1),trim(labelcoordsys(igeom)))
 enddo
 call update_test_scores(ntests,ndiff,npass)

 if (id==master) write(*,"(/,a)") '<-- GEOMETRY TEST COMPLETE'

end subroutine test_geometry

end module testgeometry
