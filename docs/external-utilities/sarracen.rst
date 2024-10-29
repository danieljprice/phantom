Sarracen
========

Sarracen is a Python library for smoothed particle hydrodynamics (SPH)
analysis and visualization which can directly read the native
binary format of Phantom.

- Docs: https://sarracen.readthedocs.io/
- Repo: https://github.com/ttricco/sarracen

Examples
--------

To read phantom data with saracen, use::

    import sarracen

    sdf, sdf_sinks = sarracen.read_phantom('dumpfile')

This gives you a pandas dataframe containing all the data in the snapshot.
You can access individual arrays using their labels, e.g.::

    sdf['vmag'] = np.sqrt(sdf['vx']**2 + sdf['vy']**2 + sdf['vz']**2)

You can easily render images using SPH interpolation similar to the functionality
in SPLASH. For example, to plot column density from a 3D simulation::

    sdf.render('rho')

or to plot a slice in the x-y plane with the inferno colour map, use::

    sdf.render('rho',xsec=0.,cmap='inferno')

For other examples, see the `documentation <https://sarracen.readthedocs.io/>`__
