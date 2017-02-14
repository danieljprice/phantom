-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Setting up live galactic discs
A. R. Pettitt
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

There are two supported setup routines, and support for creating gadget snapshots and phantom-readable ascii files.
Below are some notes on running and using these initial conditions.

arpic:
-Uses setup procedure based on Hernquist 1993 method. Designed for live disc+bulge system, in a static halo.
-Run makegal.sh to run galsetupBD.py and make a bulge-disc galaxy system in a static halo. This requires bulge.set to be present. This will create a set of asciifiles that phantomsetup reads in.
-You can easily replace the live halo with a static one, just make sure that phantom is set to call the correct bulge potential.
-If you want no bulge component then use the routine galsetupD.py, which also requires a static halo. You can use galsetupBD.py with a near zero bulge mass but it will take considerable time to set the zero-mass bulge particles.
-If you want no bulge and no halo (why?) then I advise using galsetupD.py and simply setting Mhalo to a very small value in the code itself, the setup should still be fine as no halo particles are set so no time is wasted initialising them.
-NOTE: the setup is in python so and very fast, will on the order of an hour to make a 1million particle system.
-Make sure that phantom is set to use the appropriate NFW potential (same mass and scale lengths as arpic), or else the rotation curve will not hold.
-asii2snap.f can be used to make a gadget snapshot form these asciifiles if you prefer (useful for testing).


magalie:
-Run the publicly available magalie routine (Boily et al. 2001, found in NEMO toolbox; Teuben 1995) for bulge-halo-disc systems. This provides m.dat and magalie.in data and parameter files. magalie2ascii.py then creates the same format ascii files for phantom.
-magalie2snap.f can be used to make a gadget snapshot form the m.dat and magalie.in files.


phantom:
-The arpic and magalie routines makes things in SI, which is converted to cgs then code units by phantomsetup.
-Use the galdisc setup to read in these asciifiles. It will ask a few questions (such as do you want gas included?). Make sure you say yes to using live stars or the code will set up gas itself and no read in any of these files.


gadget:
-Will simply read in the snapshots made by ascii2snap.f/magalie2snap.f and run with it. You should make sure the units are the same as used in the setup file, though they are set to gadget defaults as is, so should be fine.
-NOTE: don’t expect gadget and phantom to give the same answers. Gadget2 used fixed gravitational softening for the stars and dark matter, whereas phantom uses an adaptive method. Gadget3 also contains this method, but is not publicly available (though Gizmo is, which is based on Gadget3; Hopkins 2015). 
