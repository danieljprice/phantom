Phantom coding philosophy/developer guidelines
==============================================

What Phantom is not
~~~~~~~~~~~~~~~~~~~

Phantom is **not** a code for testing algorithms. It is a “take the best
and make it fast” production code for high resolution simulations in 3D.

What Phantom should be
~~~~~~~~~~~~~~~~~~~~~~

-  Few options
-  Fast
-  Low-memory

Developer guidelines
~~~~~~~~~~~~~~~~~~~~

-  Avoid repeated code. Do not cut-and-paste. Use subroutines and
   functions.
-  Group related functionality together. New physics = changing one
   file.
-  Actively delete old/obsolete code.
-  Add unit tests for important modules.
-  Automate testing as much as possible.
-  Please follow the :doc:`house style <styleguide>` in your coding

Much of what I try to follow comes from Joel on Software:

| `12 steps to better
  code <http://www.joelonsoftware.com/articles/fog0000000043.html>`__
| `Daily builds are your
  friend <http://www.joelonsoftware.com/articles/fog0000000023.html>`__

At some point I would like to follow these as much as possible:

`European Fortran 90 coding
standards <http://research.metoffice.gov.uk/research/nwp/numerical/fortran90/f90_standards.html>`__

Rambling justification for above
--------------------------------

`ndspmhd <http://users.monash.edu.au/~dprice/ndspmhd>`__, is the code I
use for testing algorithms. That code has many different options, all
specified at runtime. It has dynamically allocated arrays. It can be
compiled in 1, 2 or 3 dimensions. It’s very flexible. But it is bloated,
has everything I’ve ever tried in it, and would be hard to optimise or
parallelise. Phantom is meant to take the “best” of things I’ve tried in
ndspmhd and implement them fast and parallelised in 3D.

In the most important routines like densityforce, this means I will
sacrifice readability for speed. The initial design goal was to have a
very low memory footprint. That is memory should **only** be allocated
if it is strictly necessary for the physics. Hence most of the physics
is chosen at compile-time rather than runtime. Hence also Phantom uses
static rather than dynamically-allocated arrays.
