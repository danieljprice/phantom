#!/bin/bash
#
#--given a list of input files (or equivalently, directory names)
#  creates the directories, writes a Makefile in each and compiles the code
#
if [ $# -lt 1 ]; then
   echo 'Usage: '$0' inputfile(s)'
else
   scriptdir=`dirname $0`;

   for x in $@;
   do
       dir=${x/.in/};
       if [ ! -d $dir ]; then
          echo 'creating directory '$dir;
          mkdir $dir;
       fi
       if [ -e $x ]; then
          echo "moving input file $x into $dir/";
          mv $x $dir;
       fi
       echo "entering directory $dir";
       cd $dir;
       echo 'writing '$dir'/Makefile';
       $scriptdir/writemake.sh > Makefile;
       if [ X$QSYS == X ]; then
          qfile='run.qscript';
       else
          qfile='run.'$QSYS;
       fi
       echo 'writing '$dir'/'$qfile;
       make qscript INFILE=$x > $qfile;
       if [ -e ../phantom ]; then
          echo 'copying ../phantom into run directory';
          cp ../phantom .;
       else
          echo 'making phantom';
          make;
       fi
       dumpfile=`cat "$x" | grep dumpfile | sed "s/dumpfile =//g" | sed "s/\\!.*//g" | sed "s/\ //g"`
       if [ -e "../$dumpfile" ]; then
          echo "copying ../$dumpfile into run directory";
          cp ../$dumpfile .;
       else
          echo "WARNING: $dumpfile could not be copied into run directory";
       fi 
       cd -;
   done
fi
