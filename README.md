Phantom
=======

> The Phantom Smoothed Particle Hydrodynamics code

About
-----

Phantom is a 3D Smoothed Particle Hydrodynamics and Magnetohydrodynamics code for astrophysics. It was written and developed by Daniel Price with contributions from many others (see [AUTHORS](https://github.com/danieljprice/phantom/blob/master/AUTHORS)). It is designed to be a fast, 3D SPH code for production runs. It is not a code for testing algorithms. Key features include:

- Compressible hydrodynamics and magnetohydrodynamics
- Advanced shock capturing using state-of-the-art methods
- Sink particles with 4th order symplectic integration of gravitational dynamics
- Self-gravity using a momentum-conserving fast multipole method and adaptive gravitational force softening
- Multi-species dust-gas mixtures (dust-as-mixture or dust-as-particles)
- Dust growth and porosity evolution
- General relativistic hydrodynamics in Kerr, Schwarzschild and Minkowski metrics
- Radiation hydrodynamics using either implicit flux-limited diffusion or live-coupled Monte Carlo Radiative transfer
- Adaptive particle refinement
- Wind injection, dust nucleation and radiative acceleration
- Automated mapping of 1D stellar evolution models into 3D stars using Relax-O-Matic^TM
- N-body regularisation for fast simulation of star clusters
- Easy and flexible setup of binary, triple and multiple stars with discs and planets in arbitrary configurations
- Direct coupling with Monte Carlo radiative transfer for easy comparison with observations


Status
------

[![build](https://github.com/danieljprice/phantom/actions/workflows/build.yml/badge.svg)](https://github.com/danieljprice/phantom/actions/workflows/build.yml)
[![test](https://github.com/danieljprice/phantom/actions/workflows/test.yml/badge.svg)](https://github.com/danieljprice/phantom/actions/workflows/test.yml)
[![mpi](https://github.com/danieljprice/phantom/actions/workflows/mpi.yml/badge.svg)](https://github.com/danieljprice/phantom/actions/workflows/mpi.yml)
[![mcfost](https://github.com/danieljprice/phantom/actions/workflows/mcfost.yml/badge.svg)](https://github.com/danieljprice/phantom/actions/workflows/mcfost.yml)
[![Documentation](https://readthedocs.org/projects/phantomsph/badge/?version=latest)](https://phantomsph.readthedocs.io/en/latest/?badge=latest)

Links
-----

- Project homepage: http://phantomsph.github.io/
- Code repository: https://github.com/danieljprice/phantom/
- Documentation: https://phantomsph.readthedocs.org/
- Code paper: http://adsabs.harvard.edu/abs/2018PASA...35...31P

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

Getting help
------------

If you need help, please try the following, in order:

1. Check the [documentation](https://phantomsph.readthedocs.org/).
2. If you encounter a bug, [file an issue](https://github.com/danieljprice/phantom/issues/new)
3. If you want to request a feature, [start a discussion](https://github.com/danieljprice/phantom/discussions/new)
4. If you need help on how to use phantom, [start a discussion](https://github.com/danieljprice/phantom/discussions/new)

We welcome general discussion about Phantom, Smoothed Particle Hydrodynamics,
and astrophysics at the [Phantom Slack](https://phantomsph.slack.com/). However, please use the github issues for support requests.
Please note this is a community-supported code and that we rely on user contributions and community engagement to deal with issues raised on github.

Contributing
------------
We welcome contributions, including (but not limited to):

1. Starring or [forking](https://github.com/danieljprice/phantom/fork) the repository to let us know you are using the code.
2. Providing help to other users on github discussions or on the Slack channel.
3. Code, via [pull request](https://github.com/danieljprice/phantom/pulls). Please read developer section of user guide for guidelines. Our wish list is the [issue tracker](https://github.com/danieljprice/phantom/issues/).
4. Documentation, also by [pull request](https://github.com/danieljprice/phantom/pulls). Docs can be edited in the docs/ directory of the main code.
5. Suggestions for features or bug reports, via the [issue tracker](https://github.com/danieljprice/phantom/issues/new). Please file bugs via github rather than by email.

Citation
--------

Please cite [Price et al. (2018)](http://adsabs.harvard.edu/abs/2018PASA...35...31P) when using Phantom. Wherever possible, please try to also cite original references for the algorithms you are using. A partial list can be found in `docs/phantom.bib` file, or by reading the relevant sections of the paper.

Other things
-------------

For CHANGES see the release notes: https://phantomsph.readthedocs.io/en/latest/releasenotes.html.
See LICENSE file for usage and distribution conditions.

Copyright (c) 2007-2025 Daniel Price and contributors (see AUTHORS file).

