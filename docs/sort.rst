Sorting particles during analysis
=================================

Phantom has some useful sorting routines if you need to get the
particles in a particular order. These are in the module
“utils_sort.f90” in src/main

Sorting by radius with origin set to 0,0,0
------------------------------------------

An example of how to sort by radius can be found in
src/setup/phantomsetup.F90:

::

       use dim,       only:maxp
       use sortutils, only:indexxfunc
       integer :: iorder(max)
       ...
       call indexxfunc(npart,r2func_origin,xyzh,iorder)

this returns the integer array “iorder” which can be used to access the
rank of each particle in radius. For example, to print the particles in
order of radius we would then use

::

       do i=1,npart
         j = iorder(i)
         r = sqrt(dot_product(xyzh(1:3,j),xyzh(1:3,j))
         print*,'r = ',r,' particle ',j
       enddo

Using an origin other than 0,0,0
--------------------------------

We can optionally set the origin using the “set_r2func_origin” function.
Our modified version would be:

::

       use sortutils, only:set_r2func_origin,indexxfunc
       integer :: iorder(max)
       real :: x0(3)
       ...
       x0(:) = 0.
       call set_r2func_origin(x0(1),x0(2),x0(3))
       call indexxfunc(npart,r2func_origin,xyzh,iorder)

Using the centre of mass as the origin
--------------------------------------

In phantomsetup we set the origin to be the centre of mass of the
particles, by calling the “set_centreofmass” subroutine:

::

       call get_centreofmass(x0,v0,npart,xyzh,vxyzu)

Sorting by some other quantity
------------------------------

To sort by a general quantity (e.g. a real array “var”), just use the
“indexx” routine from the same module:

::

      real :: var(maxp)
      call indexx(npart,var,iorder)
