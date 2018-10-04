#!/bin/bash
#
# Wrapper script to perform nightly build checks
# for the PHANTOM code on a host machine
# Contains machine-specific environment settings
#
dir=~/phantom-nightly
#
# the following two routines prevent the nightly checks
# from running until the previous nightly checks have
# finished
#
lockfile='LOCKFILE'
get_lock()
{
  if [ -e $lockfile ]; then
     echo "$lockfile exists";
     # check if it is more than 48 hours old, if so delete
     if test `find $lockfile -mmin +2880`; then
        echo "deleting $lockfile more than 48 hours old";
        rm $lockfile;
        return 0;
     else
        return 1;
     fi
  else
     date >> $lockfile
     return 0;
  fi
}
release_lock()
{
  rm $lockfile;
}
load_modules()
{
  source /etc/profile.d/modules.sh
  source /etc/profile.d/modulecmd.sh
  source ~/.bashrc
  module purge
  source ~/.modules_buildbot
  module list
}

check_ifort()
{
  cat << EOF > hello.f90
program hello
 print*,'hello world'
end program hello
EOF
  ifort -o hello hello.f90;
  return $?;
}

check_git()
{
  if type -p git; then
     echo "git exists";
     return 0;
  else
     echo "ERROR: git not found, cannot run nightly checks";
     return 1;
  fi
}

run_nightly_build()
{
  check_ifort;
  if [ $? -eq 0 ]; then
     echo "ifort OK, running nightly build checks";
     ./nightly.sh;
  else
     echo "ERROR with ifort on MSG, build checks not run";
  fi
}

run_bots()
{
  echo "running bots...";
  cd $dir/phantom-bots/scripts;
  ./bots.sh --commit >& $dir/bots.log;
  if [ $? -eq 0 ]; then
     echo "bots ran OK";
  else
     echo "error running bots";
  fi
}
run_stats()
{
  cd $dir/phantom-bots/scripts;
  ./stats.sh $dir >& $dir/stats.log;
  if [ $? -eq 0 ]; then
     echo "stats ran OK";
  else
     echo "error running stats";
  fi
  cd $dir;
  ./phantom-bots/scripts/get_buildbot_timings.sh 20*.html > buildbot_timing.txt;
  ./phantom-bots/scripts/parse-timings.pl 20*.html > testbot_timing.txt;
  ./phantom-bots/scripts/plotstats.sh >> $dir/stats.log;
  cd $dir
}
if [ -d $dir ]; then
   cd $dir;
else
   echo "ERROR: $dir does not exist"
   exit;
fi
get_lock;
if [ $? -eq 0 ]; then
   echo "got lock";
else
   echo "aborting, previous build not finished...";
   exit;
fi
load_modules;
check_git;
if [ $? -eq 0 ]; then
   run_bots;
   run_nightly_build;
   run_stats;
else
   echo "aborting...";
fi
release_lock;
