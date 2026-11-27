
Running a simulation with stellar wind and dust formation
=========================================================

The wind and dust formation algorithms are described in `Siess et al. (2022)`, and algortihms for the radiation field in `Esseldeurs et al. (2023)`

If you find a bug, please send me an email at lionel.siess@ulb.be


Initial setup
-------------

Note that given the mass loss rate and particle mass, spheres of particles will be periodically injected, meaning that the number of particles in the simulation will keep increasing.

::

   $PHANTOM_DIR/scripts/writemake.sh wind > Makefile
   make; make setup
   ./phantomsetup wind

For an isothermal wind, use SETUP=isowind

::

   $PHANTOM_DIR/scripts/writemake.sh isowind > Makefile
   make; make setup
   ./phantomsetup wind


The first time it is run, ``phantomsetup`` asks for the number of stars, the type of wind, and then the properties of each star.
For each star you first choose a type, which sets the defaults for the questions that follow (press return to keep a default):

- an AGB star (1.2 Msun, Teff=3000 K, Reff=1 au, Mdot=1e-7 Msun/yr), the default for the primary star;
- a hot star with a trans-sonic wind (1.2 Msun, L=2e4 Lsun, Reff=2 au, Mdot=1e-5 Msun/yr, no wind velocity);
- a main-sequence star without a wind (1 Msun, L=1 Lsun, Reff=1 Rsun), the default for companions.

You are then asked for the mass, accretion radius, which two of L, Teff and Reff to give, and, if the star has a wind,
the wind velocity, temperature, acceleration parameter alpha and the equatorial velocity of the star.
For a binary you are asked for the semi-major axis, eccentricity and initial true anomaly (180, i.e. apastron, by default);
a hierarchical system gets default orbits that you edit in the .setup file. Finally you are asked how to set the resolution (see below).

At the end of these instructions, a wind.setup and wind.in file are created. Each file contains specific options that are described below.
Edit them and rerun ``./phantomsetup wind`` to get to the final setting.

Content of the .setup file
--------------------------

::

   # input file for wind setup routine
   # set any two of lum, Teff and Reff (third one is derived)

   # stars
                  nstar =           1    ! number of stars (1=single, 2=binary, >=3 hierarchical system)

   # primary star
           primary_mass =       1.600    ! primary star mass (Msun)
           primary_racc =       2.500    ! primary star accretion radius (au)
            primary_lum =       7250.    ! primary star luminosity (Lsun)
           primary_Reff =       2.517    ! primary star effective radius (au)
           primary_Rinj =           0    ! wind injection radius (e.g. 2.5 au, 0=Reff)
           primary_mdot =   3.000E-06    ! primary wind mass loss rate (in Msun/yr)
          primary_vwind =        16.5    ! primary wind velocity (in km/s)
      primary_wind_temp =       1500.    ! primary wind temperature (K, 0 = use Teff)
          primary_alpha =       0.000    ! primary alpha parameter
            primary_veq =       0.000    ! primary equatorial velocity (in km/s)

   # resolution
            iresolution =           0    ! particle mass from 0=mass,1=particles_per_shell,2=sonic point
      mass_of_particles =   1.000E-11    ! particle mass (Msun)

   # equation of state
             wind_gamma =       1.300    ! adiabatic index (initial if Krome chemistry used)

The .setup file contains the properties of each star and of its wind, the orbit if there is a companion, and the resolution.
Each star is a sink particle. Its properties, e.g. its luminosity, are used to calculate the radiation pressure.
Companions are added by increasing ``nstar``. Each extra star gets its own block of options (``secondary_*``, then ``star3_*``, ``star4_*``, ...), and star i is sink particle i.
For ``nstar = 2`` the orbit is given by the generic orbit options (``itype``, ``a``, ``e``, ``i``, ``O``, ``w``, ``f``), where the semi-major axis can be given with units (e.g. ``a = 143 au``).
For ``nstar >= 3`` the stars are placed in a hierarchical system set by the ``hierarchy`` string and the orbit of each level (e.g. ``1_a``, ``12_e``); the header of each orbit in the .setup file names its two components (e.g. ``# binary 12: secondary + star3``). Star i is the i-th label in the hierarchy string, so the order of the string sets the order of the stars (see :doc:`hierarchicalsystems` for the hierarchical labels).
All orbits start at apastron by default (true anomaly ``f = 180``).
Single-star and binary setup files from older versions, with ``icompanion_star = 0`` or ``1``, ``semi_major_axis``, ``eccentricity`` and ``true_anomaly``, are still read, with a warning.
To convert such a file to the new options, run ``phantomsetup`` once with the number of stars on the command line, e.g. ``./phantomsetup wind --nstar=2``.
Triples from older versions (``icompanion_star = 2``) are not converted automatically: keep the original file and convert a copy to ``nstar = 3`` with the hierarchical options.
Any of the stars can have a wind (set ``*_mdot > 0``), and each wind is computed from the 1D wind profile of its own star.

The wind of each star is set by:

- ``primary_mdot``: the mass loss rate. This is an input and is never changed by the code.
- ``primary_vwind``: the wind velocity at the injection radius. Set it to zero for a trans-sonic wind (``wind_type = 1``), where the code finds the velocity.
- ``primary_wind_temp``: the wind temperature at the injection radius. Zero means the stellar effective temperature is used.
- ``primary_Rinj``: the radius at which the wind is injected, with units (e.g. ``2.5 au``). Zero means
  the wind is injected at the stellar effective radius ``primary_Reff``. A larger value injects the wind above the photosphere; the stellar radius
  used for the dust temperature and radiation is still computed from the luminosity and effective temperature.

Note also that

.. math::

      \textrm{primary_lum} = 4\pi\times\textrm{primary_Reff}^2\times\sigma\times\textrm{primary_Teff}^4

so you provide 2 out of these 3 variables and the third, if omitted or set to zero, is computed from the formula above.
``phantomsetup`` asks which two to give, and writes only those two to the .setup file. To change which two are given,
replace one of them with the third (e.g. ``primary_lum`` with ``primary_Teff``).
If you provide all three, the luminosity and effective radius are kept, the effective temperature is recomputed,
and a warning is given if it differs from the value you provided.

Choosing the resolution
-----------------------

The particle mass sets the resolution. Particles are injected in spherical shells, and the distance between shells is kept close to the
distance between particles on a shell, so the particles are spaced evenly in all directions. The particle spacing at the injection radius is then

.. math::

      \Delta \approx \left(\frac{m}{\rho}\right)^{1/3}, \qquad \rho = \frac{\dot{M}}{4\pi R^2 v},

where :math:`m` is the particle mass, :math:`\dot{M}` the mass loss rate, :math:`v` the injection velocity and :math:`R` the injection radius.
A smaller particle mass resolves the density profile of the wind better. To change the resolution, e.g. for a convergence study, change the particle mass
by a factor of 8 (or 10), which changes the linear resolution :math:`\Delta` by a factor of 2. The number of particles per shell, :math:`N \approx 4\pi R^2/\Delta^2`,
follows from the particle mass:

.. math::

      N = \left(\frac{\sqrt{4\pi}\, s\, R\, \dot{M}}{m\, v}\right)^{2/3},

where :math:`s` is ``wind_shell_spacing``. Shells are injected every :math:`\Delta t = m N / \dot{M}`.
The only constraint is that the wind must be resolved azimuthally: a warning is given if there are fewer than 12 particles per shell.

The particle mass is set only by ``phantomsetup``, using ``iresolution``:

- ``iresolution = 0`` (recommended, the default): give the particle mass directly, ``mass_of_particles`` (in Msun). All stars with winds share this one particle mass.
- ``iresolution = 1``: give the number of particles per shell of the wind of the first star with a wind, ``particles_per_shell``, and the particle mass is computed from the formula above. This is mainly useful to reproduce old calculations (see below).
- ``iresolution = 2``: give the number of particle spacings between the injection radius and the sonic point, ``nzones_sonic_point``. This ensures the acceleration zone of a trans-sonic wind is resolved. It solves the 1D wind equations, sets ``wind_type = 1`` in the .in file, and is only available for a single star with a non-isothermal wind.

``phantomsetup`` prints the particle mass and, for each star with a wind, the number of particles per shell, so you can check the resolution before running.
Since the number of particles per shell must be an integer, the actual shell spacing differs slightly from ``wind_shell_spacing``; it is printed at the start of the run.
For example, for the star above, decreasing the particle mass by successive factors of 8:

================= ======================= ====================== =======================
particle mass     particles per shell     time between shells    particles per 1000 yr
================= ======================= ====================== =======================
1.9e-9 Msun       252                     0.16 yr                1.6e6
2.4e-10 Msun      1008                    0.08 yr                1.3e7
3.0e-11 Msun      4032                    0.04 yr                1.0e8
3.8e-12 Msun      16128                   0.02 yr                8.0e8
================= ======================= ====================== =======================

Migrating old input files
~~~~~~~~~~~~~~~~~~~~~~~~~

The runtime parameter ``iwind_resolution`` has been removed and is ignored with a warning, because its meaning changed over time.
Restarting from an existing dump file gives the same number of particles per shell and the same mass loss rate as before.
To set up an old calculation again, use the following in the .setup file:

=========================================================== ==========================================================================
Old setting                                                 New setting
=========================================================== ==========================================================================
``iwind_resolution = q`` (before January 2025)              ``iresolution = 1``, ``particles_per_shell = 40q(q-1)+12``,
                                                            i.e. 252 for q=3, 492 for q=4, 812 for q=5. Shells are now
                                                            Fibonacci spheres instead of icosahedra, so results agree
                                                            statistically but not bit for bit.
``iwind_resolution = N`` (from January 2025)                ``iresolution = 1``, ``particles_per_shell = N``. Runs between
                                                            March 2025 and this change usually used N+1 particles per
                                                            shell and a mass loss rate multiplied by ((N+1)/N)^(3/2);
                                                            the run log shows both (``Mdot [Msun/yr] : ... --> ...``).
                                                            Set both to those values to reproduce such runs exactly.
``iwind_resolution = 0``, ``wind_type = 0``                 ``iresolution = 0`` with the same ``mass_of_particles``.
``iwind_resolution = 0``, ``wind_type = 1``                 ``iresolution = 2``, ``nzones_sonic_point = 8``.
=========================================================== ==========================================================================

Content of the .in file
-----------------------

Options controlling particle injection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   # options for injecting/removing particles
              wind_type =           0    ! wind type 0=std, 1=transonic, 2=pulsation, 3=jets
           rfill_domain =       0.000    ! outer radius of the background density profile
     wind_shell_spacing =       1.000    ! desired ratio of sphere spacing to particle spacing
      iboundary_spheres =           5    ! number of boundary spheres (integer)
         outer_boundary =         30.    ! delete gas particles outside this radius (au)
                  rkill =      -1.000    ! deactivate particles outside this radius (<0 is off)

Here's a brief description of each of them (technical details can be found in `Siess et al. (2022)`)::

              wind_type =           0    ! wind type 0=std, 1=transonic, 2=pulsation, 3=jets

``wind_type = 0`` uses the wind velocity given in the .setup file.
``wind_type = 1`` lets the code find the trans-sonic solution. In this case you need a high wind temperature (coronal wind, > 10,000 K) so the pressure gradient can overcome the stellar gravity.
``wind_type = 2`` adds a pulsating inner boundary (``pulsation_period``, ``piston_velocity``) and ``wind_type = 3`` launches jets (``jet_edge_velocity``, ``jet_opening_angle``)::

           rfill_domain =       0.000    ! outer radius of the background density profile

fill the domain up to this radius (au) with the 1D wind solution at the start. This option can limit the effect of boundary conditions::

     wind_shell_spacing =       1.000    ! desired ratio of sphere spacing to particle spacing

the ratio between the distance between 2 consecutive shells and the distance between 2 particles on a shell.
Its value should be kept close to unity::

      iboundary_spheres =           5    ! number of boundary spheres (integer)

set the number of shells that serve as inner boundary condition for the wind::

         outer_boundary =         30.    ! delete gas particles outside this radius (au)

To limit the number of particles, delete from the memory the particles that go beyond ``outer_boundary`` (in astronomical units).
This option is slightly different from ``rkill`` where in this case the particles are declared dead and remain allocated.


Options controlling dust
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   # options controlling dust
          idust_opacity =           2    ! compute dust opacity (0=off,1=on (bowen), 2 (Gail))
              kappa_gas =   2.000E-04    ! constant gas opacity (cm²/g)
          wind_CO_ratio =       2.000    ! wind initial C/O ratio

::

          idust_opacity =           2    ! compute dust opacity (0=off,1=on (bowen), 2 (Gail))

set the type of dust formalism. Nucleation is only available with ``idust_opacity = 2``

::

              kappa_gas =   2.000E-04    ! constant gas opacity (cm²/g)

default gas opacity. Only activated if ``idust_opacity > 0``

::

          wind_CO_ratio =       2.000    ! wind initial C/O ratio

set the C/O ratio of the ejected wind material. For the moment only C-rich chemistry (C/O > 1) is implemented. Option only available with ``idust_opacity = 2``


Options controlling radiation pressure from sink particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   # options controling radiation pressure from sink particles
        isink_radiation =           3    ! sink radiation pressure method (0=off,1=alpha,2=dust,3=alpha+dust)
              alpha_rad =       1.000    ! fraction of the gravitational acceleration imparted to the gas
             iget_tdust =           1    ! dust temperature (0:Tdust=Tgas 1:T(r) 2:Flux dilution 3:Lucy 4:MCfost)
        iray_resolution =          -1    ! set the number of rays to 12*4**iray_resolution (deactivated if <0)
              tdust_exp =         0.5    ! exponent of the dust temperature profile

::

        isink_radiation =           3    ! sink radiation pressure method (0=off,1=alpha,2=dust,3=alpha+dust)

set how radiation pressure is accounted for. The star's effective gravity is given by

.. math::

              g_\mathrm{eff} = \frac{Gm}{r^2} \times (1-\alpha_\mathrm{rad}-\Gamma)

alpha is an ad-hoc parameter that allows the launching of the wind in case of a cool wind for example when dust is not accounted for.
Gamma is the Eddington factor that depends on the dust opacity. gamma is therefore <> 0 only when dust is activated (``idust_opacity > 0``)::

              alpha_rad =       1.000    ! fraction of the gravitational acceleration imparted to the gas

parameter entering in the above equation for the effective gravity::

             iget_tdust =           1    ! dust temperature (0:Tdust=Tgas 1:T(r) 2:Flux dilution 3:Lucy 4:MCfost)

defines how the dust temperature is calculated. By default one assumes Tdust = Tgas but other options are available as well.
Options 1-3 use analytical prescriptions, and option 4 uses full 3D RT using the MCfost code (under development!)::

        iray_resolution =          -1    ! set the number of rays to 12*4**iray_resolution (deactivated if <0)

If ``iget_tdust = 1-3``, the dust temperature profile is then given by an analytical prescription.
In these prescriptions (see `Esseldeurs et al. (2023)`), there is directional dependance, where the resolution of this directional dependance is set by iray_resolution.

::

              tdust_exp =         0.5    ! exponent of the dust temperature profile

If ``iget_tdust = 1``, the dust temperature profile is then given by

.. math::

              T_\mathrm{dust}(r) = T_\mathrm{star}*(R_\mathrm{star}/r)^\textrm{tdust_exp}

where T_star and R_star are the stellar (effective) temperature and radius as defined in the .setup file



**Have fun :)**
