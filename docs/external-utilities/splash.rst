Splash
======

Splash is a free and open source visualisation tool for SPH data, developed closely alongside Phantom.

- Docs: https://splash-viz.readthedocs.io/
- Repo: https://github.com/danieljprice/splash

Examples
--------

Plot column density render of all snapshots from a simulation::

    splash -r density dump_0*

Make an mp4 movie of the above, like so::

    splash -r density dump_0* --movie

which produces something like

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
       <iframe width="560" height="315" src="https://www.youtube.com/embed/bPurbeQNgvI" title="YouTube video player" frameborder="0" allow="accelerometer; encrypted-media; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
    </div>

splash can also be used to plot the energy vs time files::

   splash *.ev

including automatic recognition of column labels
