How MPI works in Phantom
========================

Phantom is parallelised with openMP within each MPI task, and with MPI between tasks.
This page gives the overall design, and where to find each part. For how to run with MPI,
see :doc:`running phantom with MPI </user-guide/mpi>`.

Domain decomposition
--------------------
The particles are divided between tasks by *sections* of the geometric tree, i.e. the same
tree of midpoint splits along the longest axis that is built for gravity on one task. The
sections are refined until none carries more than a fraction of a task's share of the work
(particles weighted by how often they are active), and each task owns a contiguous run of
sections in depth-first order. Any number of tasks works.

- src/main/mpi_sections.F90: builds the sections, cuts them between tasks (``cut_sections``)
  and renews them on full steps if the work is out of balance by more than 10 per cent
- src/main/mpi_balance.F90: moves particles to the task that owns their section

The tree
--------
Each task builds its local tree from its own sections (``maketree`` in src/main/kdtree.F90).
The cells above the sections are then built identically on every task from the gathered
section roots, with exact sizes, so that the gravity tree, and hence the forces, are the
same as on one task to round-off.

All of this, and everything else that spans tasks, is in src/main/mpi_kdtree.F90:

- ``maketreeglobal``: assigns owners, moves particles, builds the local tree and the top of the
  global tree
- ``extend_local_tree(0,nactive)``: copies the nodes of each local tree near other tasks into the
  global tree, and appends the global tree to the local tree (the *combined tree*), so that each
  local leaf's branch runs up to the global root
- ``extend_local_tree(1)``: before force, imports the nodes of other tasks that the tree walks
  of this task could open (the *locally essential tree*), and the particles of other tasks'
  leaves of large reach
- ``sync_hmax_mpi``: updates the smoothing lengths in the global tree after density

Density and force
-----------------
Each task loops over its active cells as on one task. A cell whose neighbours are partly
on other tasks is sent to those tasks, which compute its interactions with their particles
and send back the partial sums; the owner adds them up. With self-gravity, force uses the
dual tree walk on the combined tree, and a cell is sent only to the tasks owning the leaves
the walk paired it with. Momentum and angular momentum are conserved exactly, as on one task.

- src/main/dens.F90, src/main/force.F90: the cell loops, and what is sent where
- src/main/mpi_derivs.F90: the cell exchange (non-blocking sends, receives polled by one
  thread at a time)
- src/main/mpi_dens.F90, src/main/mpi_force.F90: the cell data types and their packing

What the log shows
------------------
At each full dump the log has a table of the time each task spends in each phase
(``d_lcomp``, ``f_lcomp``: its own cells; ``d_rcomp``, ``f_rcomp``: cells of other tasks;
``d_wait``, ``f_wait``: waiting for other tasks; ``f_let``: the import before force), and the
cells and data exchanged. Load imbalance shows as large differences in step time minus waiting.

Testing
-------
The tests of the global and combined tree (kdtree), of the tree gravity (gravity) and of the
derivatives (derivs) also check the MPI code when run on several tasks, e.g.::

   make phantomtest MPI=yes
   mpirun -np 4 ./bin/phantomtest kdtree

Test with task counts that are not powers of two (e.g. 3) as well.
