Frequently Asked Questions
==========================

What is the lowest alpha viscosity for a disc that Phantom can simulate?
------------------------------------------------------------------------

The default way of including disc viscosity in SPH is to modify the
artificial viscosity to mimic a physical viscosity. Details are given in
`Lodato & Price
(2010) <http://adsabs.harvard.edu/abs/2010MNRAS.405.1212L>`__ or in the
Phantom paper. The Shakura-Sunyaev alpha is related to alpha_AV
according to

::

   alpha_SS = 1/10 alpha_AV <h>/H

where is the mean smoothing length in an annulus and H is the disc scale
height at that annulus.

1. Lowest meaningful value of alpha_AV is 0.1 (see e.g. Meru & Bate
   2010; lower than this and you will get spreading independent of the
   value of alpha_AV due to particle jitter).

2. For a disc simulated with 10^6 SPH particles we typically resolve the
   scale height by 2-5 resolution lengths, giving h/H ~ 0.2 – 0.5.

3. So therefore the lowest alpha_SS you can simulate at this resolution
   is 2 - 5 x 10^-3. To reach a lower alpha_SS, you have to lower /H,
   which means using more particles.

If your desire is simply to minimise the numerical viscosity in the
code, rather than represent a physical disc viscosity, you should simply
turn off the DISC_VISCOSITY flag. This way you are just using the shock
viscosity with the switch that minimises dissipation away from shocks.
Using the disc viscosity flag turns off the switches (see Lodato & Price
2010).

How does the code identify accreted particles?
----------------------------------------------

::

   use part, only:isdead_or_accreted
   if (.not.isdead_or_accreted(xyzh(4,i)) then
      print 'particle ',i,' is not accreted'
   endif

How can I remove particles from the simulation?
-----------------------------------------------

to kill particles (not the same as accreted), use

::

   call kill_particle(i)

to remove dead particles from the array

::

   call shuffle_part(npart)

How do I make a movie?
----------------------

Type /png instead of /xw at the splash command prompt. This generates a
series of images splash_0000.png splash_0001.png etc.

use the script provided to make a movie:

::

   ~/splash/scripts/movie.sh

which just executes ffmpeg as follows

::

   ffmpeg -i splash_%04d.png -r 10 -vb 50M -bt 100M -vcodec mpeg4 -vf setpts=4.*PTS movie.mp4
