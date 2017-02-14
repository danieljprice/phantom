#!/bin/bash
dir=~/phantom-nightly
load_modules()
{
  source /etc/profile.d/modules.sh
  source ~/.profile
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
  cd $dir;
}
cd $dir;
load_modules;
check_git;
if [ $? -eq 0 ]; then
   run_bots;
   run_nightly_build;
else
   echo "aborting...";
fi
