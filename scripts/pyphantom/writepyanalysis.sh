#!/bin/bash
#
# This script writes a wrapper module for python analysis module in your working directory
#
# Written by David Liptai 2018
#

FNAME=$PWD/libanalysis.py
PYPHANTOM=${0/writepyanalysis.sh}  # This gets you the directory for where the python analysis module is

echo "import sys"                                  >  $FNAME
echo "sys.path.insert(0, '"$PYPHANTOM"')"          >> $FNAME
echo "from phantomanalysis import PhantomAnalysis" >> $FNAME
