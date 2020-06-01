Phantom
=======

> The Phantom Smoothed Particle Hydrodynamics code

About
-----

Phantom is a 3D Smoothed Particle Hydrodynamics and Magnetohydrodynamics code for astrophysics. It was written and developed by Daniel Price with contributions from many others (see AUTHORS). It is designed to be a fast 3D SPH code with a low memory footprint, for production runs. It is not a code for testing algorithms (use NDSPMHD instead).

Status
------
![testkd](https://github.com/danieljprice/phantom/workflows/testkd/badge.svg) ![mpi](https://github.com/danieljprice/phantom/workflows/mpi/badge.svg) ![GR](https://github.com/danieljprice/phantom/workflows/GR/badge.svg)

Links
-----

- Project homepage: http://phantomsph.bitbucket.io/
- Code repository: https://github.com/danieljprice/phantom/
- Documentation: https://phantomsph.readthedocs.org/
- Code paper: http://adsabs.harvard.edu/abs/2018PASA...35...31P

Code structure
--------------

The Phantom source code is structured as follows:

| ---------------- | ----------------------------------------------------- |
| `build/Makefile` | main Makefile for compiling Phantom and all utilities |
| `src/main`       | source for main code                                  |
| `src/setup`      | source for Phantomsetup utility                       |
| `src/tests`      | source for unit tests and the Phantom testsuite       |
| `src/utils`      | source for optional utilities and analysis            |

Getting help
------------

If you need help, please try the following, in order:

1. Check the [documentation](https://phantomsph.readthedocs.org/).
2. File an issue, as a [bug report](https://github.com/danieljprice/phantom/issues/new) or [feature request](https://github.com/danieljprice/phantom/issues/new), using the issue tracker.

Slack
-----

We welcome general discussion about Phantom, smoothed particle hydrodynamics,
and astrophysics at the [Phantom Slack](https://phantomsph.slack.com/).

Citation
--------

Please cite [Price et al. 2018](http://adsabs.harvard.edu/abs/2018PASA...35...31P) when using phantom. Wherever possible, please try to also cite original references for the algorithms you are using. A partial list can be found in `docs/phantom.bib` file, or by reading the relevant sections of the phantom paper.

Licence
-------

See LICENCE file for usage and distribution conditions.

Copyright (c) 2007-2020 Daniel Price and contributors (see AUTHORS file).

Release notes
-------------

For CHANGES see the release notes: https://phantomsph.readthedocs.io/en/latest/releasenotes.html.
