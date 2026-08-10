Injecting particles from another simulation
===========================================

Initial setup
-------------

To ensure particle masses and units are consistent between the donor and
receiver simulations, it is recommended to use ``phantommoddump`` with the
``inject_sim`` module to set up the new simulation.

::

   make SRCINJECT=inject_sim.f90
   make moddump SRCINJECT=inject_sim.f90
   ./phantommoddump YOUR_EXISTING_SIMULATION YOUR_NEW_SIMULATION TIME

``phantommoddump`` may prompt for parameters on the first run. If so, run
the same command again after answering the prompts.

At the end of these steps you will have an initial dump for the new
simulation and a corresponding ``.in`` file.

See also :doc:`moddump-recipes` for other moddump workflows.

Content of the .in file
-----------------------

Options controlling particle injection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   # options controlling particle injection
       start_dump =  'dump_00000'   ! dumpfile to start for injection
         r_inject =     5.000E+14   ! radius to inject tde outflow (in cm)
       final_dump =  'dump_02000'   ! stop injection after this dump

``start_dump``
   Dump file to begin injection. The code starts injection when the
   simulation time reaches the time of ``start_dump``. After each injected
   dump is used, the dump number in ``.in`` is incremented automatically when
   a full output dump is written.

If dumps live in another directory, use a **relative** path with quotation
marks, for example::

       start_dump =  'PATH/TO/YOUR/OTHER/DIR/dump_00000'

``r_inject``
   Injection radius (cm). For TDE outflow setups, particles that cross this
   radius from inside to outside in the donor simulation are injected into the
   new simulation.

``final_dump``
   Stop injection after this dump is reached. If ``start_dump`` includes a path,
   ``final_dump`` must use the same path prefix.
