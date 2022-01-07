How to set up and run a common envelope binary simulation
=========================================================

Polytropic model
----------------

Use SETUP=star:

If not specified, use the default options.

::

   1. $PHANTOM_DIR/scripts/writemake.sh star > Makefile

A. setup the polytropic star

::

   2. make setup
   3. ./phantomsetup poly (option 2 Polytrope, number of particles = 1000, Relax star automatically = yes)
   4. vim  poly.setup, to change :
   write_rho_to_file =           T    ! write density profile to file

  5. ./phantomsetup poly.setup
  6. make
  7. vim poly.in (to change values like tmax=500.00, dtmax=50.00 and nfulldump = 1)
  8. ./phantom poly.in


B. include a sink particle inside the polytropic star

::

  9. make moddump
  10. ./phantommoddump poly_00002 core_00000.tmp 0.0 (option 4, mass core = 0.01, hsoft = 0.1). Replace poly_00002 by the best-suited outcome of step 8 and replace core_00000.tmp by any other name you may like.
  11. vim core.in (optional to change values like tmax or dtmax)
  12. ./phantom core.in

C. Setup the binary system

::

  13. ./phantommoddump core_00002 binary_00000.tmp 0.0 (option 1, companion mass = 0.4, orbital separation = 0.8, softening length for companion = 0.1)
  14. vim binary.in (optional, tmax=200.00, dtmax=0.100)
  15. ./phantom binary.in

D. Analysis

::

  16. make analysis
  17. ./phantomanalysis binary_000* (option 1, corotating frame? No) This step is made to obtain the orbital parameters of the binary system in the file separation_vs_time.ev
  18. splash separation_vs_time.ev


Real star model
---------------
