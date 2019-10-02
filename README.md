Phantom
=======

> The Phantom Smoothed Particle Hydrodynamics code

About
-----

Phantom is a 3D Smoothed Particle Hydrodynamics and Magnetohydrodynamics code for astrophysics. It was written and developed by Daniel Price with contributions from many others (see AUTHORS). It is designed to be a fast 3D SPH code with a low memory footprint, for production runs. It is not a code for testing algorithms (use NDSPMHD instead).

Links
-----

- Project homepage: http://phantomsph.bitbucket.io/
- Code repository: https://bitbucket.org/danielprice/phantom/
- Documentation: https://phantomsph.readthedocs.org/

Code structure
--------------

The Phantom source code is structured as follows:

|                  |                                                       |
| ---------------- | ----------------------------------------------------- |
| `build/Makefile` | main Makefile for compiling Phantom and all utilities |
| `src/main`       | source for main code                                  |
| `src/setup`      | source for Phantomsetup utility                       |
| `src/tests`      | source for unit tests and the Phantom testsuite       |
| `src/utils`      | source for optional utilities and analysis            |

Citation
--------

See `docs/phantom.bib` file for relevant papers to cite when using Phantom.

Licence
-------

See LICENCE file for usage and distribution conditions.

Copyright (c) 2007-2019 Daniel Price and contributors (see AUTHORS file).

Release notes
-------------

For CHANGES see the release notes: https://phantomsph.readthedocs.io/en/latest/releasenotes.html.
