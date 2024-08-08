Splash
======

Splash is a free and open source visualisation tool for SPH data, developed closely alongside Phantom.

- Docs: https://splash-viz.readthedocs.io/
- Repo: https://github.com/danieljprice/splash

Examples
--------

::
    splash -r density dump_0*

Plots column density render of all snapshots from a simulation

::
    splash -r density dump_0* --movie

Make an mp4 movie of the above

::
    splash *.ev

Plot energy vs time files, including automatic recognition of column labels
