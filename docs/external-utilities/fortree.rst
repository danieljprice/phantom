Fortree
=======

Fortree is a python app that makes visualisation trees of Fortran code.

Before you start
----------------

Make sure you have `Python 3 <https://www.python.org/downloads/>`__ and
`Graphviz <https://graphviz.gitlab.io/download/>`__ installed on your
machine.

1. Get Fortree
--------------

Available here: https://github.com/EstherTaillifet/fortree

2 . Edit the initialisation file
--------------------------------

Edit init.txt. Specify the name of the part of the code you want to
visualise and if it is a PROGRAM a MODULE of a ROUTINE. Specify the path
to the src of Phantom and the path of the file that defines the part of
the code you wan’t to display. This is the starting point for Fortree.
Choose the type of Fotree you want to make. CALL will make a call tree
of the code. You can choose to display all the calls or only the call of
routines defined in the code. You can display all the levels or a number
of level you may choose. DEF will make a definition tree and is thus
only available for PROGRAM or MODULE. DEP will make a module dependency
tree and is thus only available for MODULE.

::

   TREE_ROOT_NAME                  evolve                                                          # Name of the part of the code to fortree.
   TREE_ROOT_TYPE                  MODULE                                                          # What is it ? Values: "PROGRAM", "MODULE" or "ROUTINE".
   ROOT_FILE_PATH                  /Phantom/src/main/evolve.F90                        # Path to the root file.
   DIRECTORY_TO_PARSE              /Users/esther/Research/Code/Phantom/src                         # Directory to your fortran files. 

   FORTREE_TYPE                    CALL_TREE       # What type of tree do you need ? CALL_TREE, DEF_TREE or DEP_TREE (modules dependencies, for modules only). 
   If you chose CALL_TREE:
   SHOW_ONLY_DEF                   YES             # Do you whant to display only the calls of functions defined in your code ? "YES" or "NO" (= show all the calls).
   N_LEVELS                        3               # How many levels do you want to display ? Integer (min val = 2) or "ALL".

   OUTPUT_NAME                     evolve_call     # A name for your outputs ? (Default = fortree)

3 . Run Fortree
---------------

::

   ~$ Python3 fortree.py init.txt
   ===============================================================
   =========================  FORTREE  ===========================
   ===============================================================
   ---------------------------------------------------------------
   Directory =  Phantom/src
   Root path =  Phantom/src/main/evolve.F90
   Tree root type =  MODULE
   Tree root name =  evolve
   Render type =  CALL_TREE
   Show only def =  True
   N levels to render =  3
   Output name =  evolve_call
   ---------------------------------------------------------------
   ---------------------------------------------------------------
   Render done.
   Tree file: evolve_call.eps 
   ---------------------------------------------------------------
   ===============================================================
   Fortree took 7.529879808425903 seconds to run.
   ===============================================================

4 . Enjoy your tree
-------------------

The output file is in .eps. |evolve_call_3.png|

Please report if you find any bugs (esther.taillifet@univ-lyon1.fr).
Thanks a lot.

.. |evolve_call_3.png| image:: https://bitbucket.org/repo/MyzKMr/images/1678860101-evolve_call_3.png

