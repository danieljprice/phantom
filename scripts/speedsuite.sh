#!/bin/bash
phantomdir=~/phantom
#tests='sedov kozai warp turb128'
tests='sedov'
export IND_TIMESTEPS=yes
mpiscaling='0 1 2 3 4';
#mpiscaling='0';
ompscaling='0 1 2 3 4 5 6';
#ompscaling='0';
for test in $tests; do
 scalingfile="$test-scaling.txt";
 if [ -s $scalingfile ]; then
    echo "moving old $scalingfile";
    mv $scalingfile $scalingfile.bak;
 fi
 # write header
 echo "# nmpi nomp ntotal wall(s) cpu(s)" > $scalingfile;
 for nomp in $ompscaling; do
  if [ $nomp -gt 0 ]; then
     export OPENMP=yes;
     export OMP_NUM_THREADS=$nomp;
     nomps=$nomp;
     ompdir="-omp-$nomp";
  else
     export OPENMP=no;
     nomps=1;
     ompdir='';
  fi
  for nmpi in $mpiscaling; do
      if [ $nmpi -gt 0 ]; then
         mpidir="-mpi-$nmpi";
         mpicmd="mpiexec -np $nmpi";
         mpicmdsetup="mpiexec -np 1";
         export MPI=yes;
         nmpis=$nmpi;
      else
         mpidir='';
         mpicmd='';
         mpicmdsetup='';
         export MPI=no;
         nmpis=1;
      fi
      dir="$test$mpidir$ompdir";
      if [ ! -d $dir ]; then
         mkdir $dir;
      fi
      # check if it is a directory
      if [ -d $dir ]; then
         echo "--- $dir ---";
         cd $dir;
         case $test in
         sedov )
            setup=sedov;
            echo "50" > $test.setup;;
         kozai )
            setup=binarydisc;;
         warp )
            setup=warp;;
         turb128 )
            setup=turbdrive;;
         esac
         $phantomdir/scripts/writemake.sh $setup > Makefile;
         make setup >& $test.makesetup;
         make >& $test.makelog;
         if [ -e phantomsetup ]; then
            #
            #--run the simulation
            #
            $mpicmdsetup ./phantomsetup $test < $test.setup >& $test.setuplog;
            $mpicmd ./phantom "$test.in" >& $test.log;

            #
            #--validate the result
            #
            energy=`grep 'Etot' $test.log | tail -1`;
            echo $energy;
            #
            #--get timings
            #
            timing=`grep 'Since code start' $test.log | tail -1`;
            wall=`grep 'Since code start' $test.log | tail -1 | cut -d: -f 3`;
            cpu=`grep 'Since code start' $test.log | tail -1 | cut -d: -f 4`;
            wall=${wall/s cpu/};
            cpu=${cpu/s cpu\/wall/};
            echo $timing;
            echo $nmpi $nomp $(( nmpis*nomps )) $wall $cpu >> ../$scalingfile;
         else
            echo "error: could not build phantomsetup";
         fi
         cd ..;
      else
         echo "error: $test directory does not exist, skipping";
      fi
  done
 done
done
