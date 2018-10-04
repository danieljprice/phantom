#!/bin/sh

#Set the disc postions (gas+disc)
python galsetupBD.py bulge.set p d
#Set the bulge postions
python galsetupBD.py bulge.set p b
#Set the disc velocities (gas+disc)
python galsetupBD.py bulge.set v d
#Set the bulge velocities
python galsetupBD.py bulge.set v b
#finish up
python galsetupBD.py bulge.set f d
