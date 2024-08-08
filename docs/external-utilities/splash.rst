Splash
======

Splash is a free and open source visualisation tool for SPH data, developed closely alongside Phantom.

- Docs: https://splash-viz.readthedocs.io/
- Repo: https://github.com/danieljprice/splash

Examples
--------

Plot column density render of all snapshots from a simulation:

    splash -r density dump_0*

Make an mp4 movie of the above, like so:

    splash -r density dump_0* --movie

which produces something like:

   https://zenodo.org/records/11438154/files/Priceetal24_figure1_logdensity_schwarzschild_4m_adiabatic.mp4

splash can also be used to plot the energy vs time files:

   splash *.ev

including automatic recognition of column labels
