#!/bin/bash
#
#--given a list of .setup files
#  creates the directories, runs phantom setup
#  and prepares the run for queue submission
#
if [ $# -lt 1 ]; then
   echo 'Usage: '$0' setupfile(s)'
else
   scriptdir=`dirname $0`;

   for x in $@;
   do
       dir=${x/.setup/};
       infile=${x/.setup/.in};
       if [ ! -d $dir ]; then
          echo 'creating directory '$dir;
          mkdir $dir;
       fi
       if [ -e $x ]; then
          echo "moving setup file $x into $dir/";
          mv $x $dir;
       fi
       echo "entering directory $dir";
       cd $dir;
       echo 'writing '$dir'/Makefile';
       ${scriptdir}/writemake.sh > Makefile;
       if [ X$QSYS == X ]; then
          qfile='run.qscript';
       else
          qfile='run.'$QSYS;
       fi
       echo 'writing '$dir'/'$qfile;
       make qscript INFILE=${infile} > ${qfile};
       if [ -e ../phantom ]; then
          echo 'copying ../phantom into run directory';
          cp ../phantom .;
       else
          echo 'making phantom';
          make;
       fi
       if [ -e ../phantomsetup ]; then
          echo 'copying ../phantomsetup into run directory';
          cp ../phantomsetup .;
       else
          echo 'making phantomsetup';
          make setup;
       fi
       ./phantomsetup $x;
       cd -;
   done
fi
