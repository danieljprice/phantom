!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: sortutils
!
!  DESCRIPTION:
!  contains low level sorting utilities
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module sortutils
 implicit none
 public :: indexx,indexxfunc,r2func,r2func_origin,set_r2func_origin

 private
 real, private :: x0,y0,z0

contains

!----------------------------------------------------------------
!+
!  function returning radius squared given xyzh
!+
!----------------------------------------------------------------
real function r2func(xyzh)
 real, intent(in) :: xyzh(4)

 r2func = xyzh(1)**2 + xyzh(2)**2 + xyzh(3)**2

end function r2func

!----------------------------------------------------------------
!+
!  function returning radius squared given xyzh and origin pos
!+
!----------------------------------------------------------------
subroutine set_r2func_origin(xorigin,yorigin,zorigin)
 real, intent(in) :: xorigin,yorigin,zorigin

 x0 = xorigin
 y0 = yorigin
 z0 = zorigin

end subroutine set_r2func_origin

real function r2func_origin(xyzh)
 real, intent(in) :: xyzh(4)

 r2func_origin = (xyzh(1)-x0)**2 + (xyzh(2)-y0)**2 + (xyzh(3)-z0)**2

end function r2func_origin
!----------------------------------------------------------------
!+
!  standard sorting routine using Quicksort
!+
!----------------------------------------------------------------
subroutine indexx(n, arr, indx)
 integer, parameter :: m=7, nstack=500
 integer, intent(in)  :: n
 real,    intent(in)  :: arr(n)
 integer, intent(out) :: indx(n)

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer :: istack(nstack)
 real :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l < m) then
    do j = l + 1, ir
       indxt = indx(j)
       a = arr(indxt)
       do i = j - 1, 1, -1
          if (arr(indx(i)) <= a) goto 2
          indx(i + 1) = indx(i)
       enddo
       i = 0
2      indx(i + 1) = indxt
    enddo
    if (jstack==0) return
    ir = istack(jstack)
    l = istack(jstack - 1)
    jstack = jstack - 2
 else
    k = (l + ir)/2
    itemp = indx(k)
    indx(k) = indx(l + 1)
    indx(l + 1) = itemp
    if (arr(indx(l + 1)) > arr(indx(ir))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(ir)
       indx(ir) = itemp
    endif
    if (arr(indx(l)) > arr(indx(ir))) then
       itemp = indx(l)
       indx(l) = indx(ir)
       indx(ir) = itemp
    endif
    if (arr(indx(l + 1)) > arr(indx(l))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(l)
       indx(l) = itemp
    endif
    i = l + 1
    j = ir
    indxt = indx(l)
    a = arr(indxt)

3   continue
    i = i + 1
    if (arr(indx(i)) < a) goto 3
4   continue
    j = j - 1
    if (arr(indx(j)) > a) goto 4
    if (j < i) goto 5
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    goto 3

5   indx(l) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    if (jstack > nstack) then
       print*,'fatal error!!! stacksize exceeded in sort'
       print*,'need to set parameter nstack higher in subroutine indexx '
       stop
    endif
    if (ir - i + 1 >= j - l) then
       istack(jstack) = ir
       istack(jstack - 1) = i
       ir = j - 1
    else
       istack(jstack) = j - 1
       istack(jstack - 1) = l
       l = i
    endif
 endif

 goto 1
end subroutine indexx

!----------------------------------------------------------------
!+
!  customised low-memory sorting routine using Quicksort
!  sort key value on-the-fly by calling the function func
!  which can be any function of the particle positions
!+
!----------------------------------------------------------------
subroutine indexxfunc(n, func, xyzh, indx)
 integer, parameter :: m=7, nstack=500
 integer, intent(in)  :: n
 real, external :: func
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(out) :: indx(n)

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer :: istack(nstack)
 real :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l < m) then
    do j = l + 1, ir
       indxt = indx(j)
       a = func(xyzh(:,indxt))
       do i = j - 1, 1, -1
          if (func(xyzh(:,indx(i))) <= a) goto 2
          indx(i + 1) = indx(i)
       enddo
       i = 0
2      indx(i + 1) = indxt
    enddo
    if (jstack==0) return
    ir = istack(jstack)
    l = istack(jstack - 1)
    jstack = jstack - 2
 else
    k = (l + ir)/2
    itemp = indx(k)
    indx(k) = indx(l + 1)
    indx(l + 1) = itemp
    if (func(xyzh(:,indx(l+1))) > func(xyzh(:,indx(ir)))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(ir)
       indx(ir) = itemp
    endif
    if (func(xyzh(:,indx(l))) > func(xyzh(:,indx(ir)))) then
       itemp = indx(l)
       indx(l) = indx(ir)
       indx(ir) = itemp
    endif
    if (func(xyzh(:,indx(l+1))) > func(xyzh(:,indx(l)))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(l)
       indx(l) = itemp
    endif
    i = l + 1
    j = ir
    indxt = indx(l)
    a = func(xyzh(:,indxt))

3   continue
    i = i + 1
    if (func(xyzh(:,indx(i))) < a) goto 3
4   continue
    j = j - 1
    if (func(xyzh(:,indx(j))) > a) goto 4
    if (j < i) goto 5
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    goto 3

5   indx(l) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    if (jstack > nstack) then
       print*,'fatal error!!! stacksize exceeded in sort'
       print*,'need to set parameter nstack higher in subroutine indexxfunc '
       stop
    endif
    if (ir - i + 1 >= j - l) then
       istack(jstack) = ir
       istack(jstack - 1) = i
       ir = j - 1
    else
       istack(jstack) = j - 1
       istack(jstack - 1) = l
       l = i
    endif
 endif

 goto 1
end subroutine indexxfunc

end module sortutils
