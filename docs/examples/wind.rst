
Running a simulation with stellar wind and dust formation
=========================================================

The wind and dust formation algorithms are described in `Siess et al. (2022)`, and algortihms for the radiation field in `Esseldeurs et al. (2023)`

If you find a bug, please send me an email at lionel.siess@ulb.be


Initial setup
-------------

Note that given the mass loss rate and particle's mass, spheres of particles will be periodically injected, meaning that the number of particles in the simulation will keep increasing.

::

   $PHANTOM_DIR/scripts/writemake.sh wind > Makefile
   make; make setup
   ./phantomsetup wind

For an isothermal wind, use SETUP=isowind

::

   $PHANTOM_DIR/scripts/writemake.sh isowind > Makefile
   make; make setup
   ./phantomsetup wind


At the end of these instructions, a wind.setup and wind.in file are created. Each file contains specific options that are described below.
Note that you may need to run ``./phantomsetup wind`` a few times to get to the final setting. 

Content of the .setup file
--------------------------

::

   # input file for wind setup routine
           primary_mass =       1.500    ! primary star mass (Msun)
           primary_racc =       1.000    ! primary star accretion radius (au)
            primary_lum =   2.000E+04    ! primary star luminosity (Lsun)
           primary_Teff =       3000.    ! primary star effective temperature (K)
           primary_Reff =       1.000    ! primary star effective radius (au)
        icompanion_star =           0    ! set to 1 for a binary system
      mass_of_particles =   1.000E-11    ! mass resolution (Msun)
             wind_gamma =       1.666    ! adiabatic index (initial if Krome chemistry used)


The .setup file contains the stellar properties and sets the mass of the particle (see however  ``iwind_resolution``).
Each star is considered as a sink particles and its properties, e.g. its luminosity, will be used to calculate the radiation pressure.
Companions can be added using the icompanion_star parameter.

Note also that 

.. math::

      \textrm{primary_lum} = 4\pi\times\textrm{primary_Reff}^2\times\sigma\times\textrm{primary_Teff}^4 
      
so you only need to provide 2 out of these 3 variables. 

- If you set one of these variables to zero, it will be  recalculated according to the previous formula. 

- If you provide all the quantites, the radius will be recalculated

Content of the .in file
-----------------------

Options controlling particle injection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  # options controlling particle injection
             sonic_type =           0    ! find transonic solution (1=yes,0=no)
          wind_velocity =         15.    ! injection wind velocity (km/s, if sonic_type = 0)
     wind_inject_radius =       1.100    ! wind injection radius (au, if 0 take Rstar)
         wind_mass_rate =   1.000E-05    ! wind mass loss rate (Msun/yr)
       wind_temperature =       2500.    ! wind temperature at the injection point (K)
       iwind_resolution =          10    ! if<>0 set number of particles on the sphere, reset particle mass
           nfill_domain =           0    ! number of spheres used to set the background density profile
     wind_shell_spacing =       1.000    ! desired ratio of sphere spacing to particle spacing
      iboundary_spheres =           5    ! number of boundary spheres (integer)
         outer_boundary =         50.    ! delete gas particles outside this radius (au)

Here’s a brief description of each of them (remember that technical details can be found in `Siess et al. (2023)`::

             sonic_type =           0    ! find transonic solution (1=yes,0=no)

decide whether you set the initial wind velocity (``sonic_type = 0``) or if you let the code find the trans-sonic solution.
In this latter case, you need a high wind temperature (coronal wind) so the pressure gradient can overcome the stellar gravity::

          wind_velocity =         15.    ! injection wind velocity (km/s, if sonic_type = 0)

set the launching wind velocity (if sonic_type = 0)::

     wind_inject_radius =       1.100    ! wind injection radius (au, if 0 take Rstar)

set the distance from the star's center where the wind is launched. If set to zero, the stellar surface is assumed::

         wind_mass_rate =   1.000E-05    ! wind mass loss rate (Msun/yr)

set the mass loss rate::

       wind_temperature =       2500.    ! wind temperature at the injection point (K)

set the wind temperature. For trans-sonic solution, this value needs to be high (> 10,000 K)::

       iwind_resolution =          10    ! if<>0 set number of particles on the sphere, reset particle mass

set the number of particles to be launched and given the mass loss rate determines the particle's mass.
If set to zero, the particle mass defined in the .setup file is used and the code finds the corresponding number of particles to be launched::

           nfill_domain =           0    ! number of spheres used to set the background density profile

set a background density profile. This option can limit the effect of boundary conditions. The larger nfill_domain, the bigger the domain::

     wind_shell_spacing =       1.000    ! desired ratio of sphere spacing to particle spacing

set the resolution of the simulation.
This parameters gives the ratio between the distance of 2 particles on an ejected sphere and the distance between 2 consecutive spheres.
Its value should be kept close to unity that::

      iboundary_spheres =           5    ! number of boundary spheres (integer)

set the number of shells that serve as inner boundary condition for the wind::

         outer_boundary =         50.    ! delete gas particles outside this radius (au)

To limit the number of particles, delete from the memory the particles that go beyond ``outer_boundary`` (in astronomical unit).
This option is slightly different from ``rkill`` where in this case the particles are declared dead and remained allocated.


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

defines how the dust temperature is calculated. By default one assumes Tdust = Tgas but other options are availabe as well. 
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
