Setting up hierarchical (triple) systems embedded in accretion discs
====================================================================

**DISCLAIMER**: Although in the future it will be possible to set up a generic hierarchical systems with discs, up to now the only tested and supported configuration consists in an hierarchical triple system surrounded by a circumtriple disc.

In the *disc* and *dustydisc* environment it is possible to set up an hierarchical triple system setting *nsinks* = 3 in the .setup file.
The hierarchical triple system initial condition is built by means of the osculating elements of the outer and of the inner binary orbit. Firstly *phantomsetup* builds a binary with the parameters of the outer orbit. Then, one of the sinks of the outer binary is substituted with an inner binary of the same mass of the substituted sink. The center of mass of the inner binary follows the orbit of the substituted sink. The subroutine devoted to this is called *set_multiple* and can be called as much times as needed, in order to build a generic hierarchical system. 

The following additional parameters are needed in the .setup file, in order to define the triple orbit:

- *q2*, the mass ratio of the inner binary;
- *subst*, an integer that refers to the outer binary's sink to be substituted (see hierarchical index section);
- *binary2_a*, *binary2_e*, *binary2_i* , *binary2_O*, *binary2_w*, *binary2_f*, respectively the semi-major axis, eccentricity, inclination, PA of the ascending node and argument of periapsis of the inner binary orbit;
- *accr2a*, *accr2b* respectively the accretion radius of the primary and of the secondary of the inner binary.

The *subst* index can be specified as positive or negative. In the former case, the orientation of the orbital plane of the inner binary will refer to the orbital plane of the outer binary. In the latter case, the orientation of the inner binary's orbital plane will refer to the sky (as the outer binary orbital plane).

The output interval of Phantom is expressed in terms of the outer binary period. If you want the output interval as fraction of the inner binary's orbital period, you have to express *deltat* in the .setup file as negative.

Up to now, if *subst* is positive it is possible to set hierarchical triples with *binary_O* = 0 only.


Hierarchical index
------------------

The following is the nomenclature used in *setup_disc.f90* in order to refer to different hierarchical levels and stars in a general hierarchical system. In particular the *subst* parameter in the .setup file is the hierarchical index of the star that has to be substituted.

::
   
          __ 1 __               Single star
        /         \
      11           12           Binary 
      |           /  \
      11        121  122        Hierarchical triple
     /  \        |    |
   111  112     121  122        Hierarchical quadruple
           ...

The hierarchical indexes (1, 11, 12, 112, ...) are used to identify the hierarchical position of the star in the system. In the code they are used to specify the sink to be substituted with a binary when calling the *set_multiple* routine in order to build an hierarchical system.

HIERARCHY file
--------------
When *phantomsetup* builds up an hierarchical system it stores some information about the initial condition of each star in the system and about its hierarchical level. Each line of the HIERARCHY file refers to a different star (that could have been already substituted).

For each star are listed, in order:

- xyzmh_ptmass index, i.e. the index of the sink in the sinks array (if 0 the sink has been substituted)
- the hierarchical index of the star
- the mass of the star
- the mass of the companion of the star (either it is a single star or an hierarchical subsystem)
- the semi-major axis of the star's orbit
- the eccentricity of the star's orbit
- the orbital period of the star
- the inclination of the star's orbit
- the argument of pericenter of the star's orbit
- the ascending node longitude of the star's orbit
