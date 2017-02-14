#!/bin/sh

#Set the disc postions (gas+disc)
python galsetupD.py disc.set p
#Set the disc velocity
python galsetupD.py disc.set v
#finish up
python galsetupD.py disc.set f
