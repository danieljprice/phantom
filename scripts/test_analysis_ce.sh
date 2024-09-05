#!/bin/bash
#
# script to perform unit tests of the analysis_common_envelope module
# checks that answering the prompts a certain way produces expected output
#
# contributed by Miguel Gonzalez-Bolivar, Feb 2023
#

# delete everything from previous tests
rm -f *.ev *txt

# grab the data file from the server if it doesn't already exist
file=binary_01000
if [ ! -f $file ]; then
   curl -k https://zenodo.org/records/13163487/files/binary_01000 -o binary_01000; err=$?;
   if [ $err -gt 0 ]; then
      exit $err;
   fi
fi

# perform phantomanalysis tests
./phantomanalysis $file > /dev/null << SEP
1
no
SEP

./phantomanalysis $file > /dev/null << BOUND
2
no
2
1.667
0.6984
0.0142
BOUND

./phantomanalysis $file > /dev/null << ENERGIES
3
no
2
1.667
ENERGIES
