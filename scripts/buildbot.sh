#!/bin/bash
#
# The Phantom build-bot
#
# Cycles through all defined setups in build/Makefile
# and checks that these compile. Results are collated into
# tables for the html pages.
#
# You can run this script yourself from the scripts directory as follows:
# cd phantom/scripts; ./buildbot.sh
#
# Outputs are put in the phantom/logs directory
#
# Written by Daniel Price, 2012-2017, daniel.price@monash.edu
#
if [ X$SYSTEM == X ]; then
   echo "Error: Need SYSTEM environment variable set to check PHANTOM build";
   echo "Usage: $0 [max idim to check] [url]";
   exit;
fi
if [ $# -gt 0 ]; then
   maxdim=$1;
   if (($maxdim > 0)) && (($maxdim < 2000000000)); then
      echo "Using maxdim = $maxdim";
   else
      echo "Usage: $0 [max idim to check] [url]";
      exit;
   fi
else
   maxdim=11000000;
fi
pwd=$PWD;
phantomdir="$pwd/../";
listofcomponents='main utils setup';
#listofcomponents='setup';
#
# change the line below to exclude things that depend on external libraries from the build
#
nolibs='MESAEOS=no'
if [ $# -gt 1 ]; then
   url=$2;
   echo "url = $url";
else
   url='';
fi
if [ ! -e $phantomdir/scripts/$0 ]; then
   echo "Error: This script needs to be run from the phantom/scripts directory";
   exit;
fi
#
# make the subdirectory "logs" if it does not exist
#
if [ ! -d $phantomdir/logs ]; then
   cd $phantomdir;
   mkdir logs;
   cd $pwd;
fi
htmlfile="$phantomdir/logs/build-status-$SYSTEM.html";
faillog="$phantomdir/logs/build-failures-$SYSTEM.txt";
faillogsetup="$phantomdir/logs/setup-failures-$SYSTEM.txt";
#
# delete old log files
#
if [ -e $htmlfile ]; then
   rm $htmlfile;
fi
if [ -e $faillog ]; then
   rm $faillog;
fi
if [ -e $faillogsetup ]; then
   rm $faillogsetup;
fi
#
# utility routine for printing results to html file
#
red="#FF0000";
amber="#FF6600";
green="#009900";
white="#FFFFFF";
pass='pass';
fail='fail';
warn='warn';
ncheck=0;
ntotal=0;
nfail=0;
nwarn=0;
print_result()
{
  text=$1;
  result=$2;
  ncheck=$(( ncheck + 1 ));
  ntotal=$(( ntotal + 1 ));
  case $result in
  $pass )
     colour=$green;;
  $warn )
     nwarn=$(( nwarn + 1 ));
     colour=$amber;;
  $fail )
     nfail=$(( nfail + 1 ));
     colour=$red;;
  * )
     colour=$white;;
  esac
  echo "<td bgcolor=\"$colour\">$text</td>" >> $htmlfile;
}
#
# procedure for checking phantomsetup utility works properly
#
check_phantomsetup ()
{
   myfail=0;
   setup=$1;
   if [ -e ./bin/phantomsetup ]; then
      print_result "exists" $pass;
   else
      print_result "FAIL: phantomsetup does not exist" $fail;
      myfail=$(( myfail + 1 ));
   fi
   dirname="test-phantomsetup";
   cd /tmp/;
   rm -rf $dirname;
   mkdir $dirname;
   cd /tmp/$dirname;
   cp $phantomdir/bin/phantomsetup .;
   myinput="\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
   prefix="myrun";
   echo -e "$myinput" > myinput.txt;
   sed '/-e/d' myinput.txt > mycleanin.txt
   ./phantomsetup $prefix < mycleanin.txt > /dev/null; err=$?;
   if [ $err -eq 0 ]; then
      print_result "runs" $pass;
   else
      print_result "FAIL: requires input other than 'Enter'" $fail;
      myfail=$(( myfail + 1 ));
   fi
   ./phantomsetup $prefix < myinput.txt > /dev/null;
   ./phantomsetup $prefix < myinput.txt > /dev/null;
   if [ -e "$prefix.setup" ]; then
      print_result "creates .setup file" $pass;
   else
      print_result "no .setup file" $warn;
   fi
   if [ -e "$prefix.in" ]; then
      print_result "creates .in file" $pass;
   else
      print_result "FAILED to create .in file after 3 attempts" $fail;
      myfail=$(( myfail + 1 ));
   fi
   if [ $myfail -gt 0 ]; then
      echo $setup >> $faillogsetup;
   fi
}
#
# get list of targets, components and setups to check
#
allsetups=`grep 'ifeq ($(SETUP)' $phantomdir/build/Makefile | grep -v skip | cut -d, -f 2 | cut -d')' -f 1`
for component in $listofcomponents; do
case $component in
 'setup')
   text="$component runs, creates .setup and .in files with no unspecified user input";
   listofsetups=$allsetups;
   listoftargets='setup';;
 'utils')
   text="$component build";
   listofsetups='test';
   listoftargets='utils';;
 *)
   text='build';
   listofsetups=$allsetups;
   listoftargets='phantom setup analysis moddump';;
esac
#
# write html header
#
echo "<h2>Checking Phantom $text, SYSTEM=$SYSTEM</h2>" >> $htmlfile;
echo "Build checked: "`date` >> $htmlfile;
echo "<table>" >> $htmlfile;
#
#--loop over each setup
#
for setup in $listofsetups; do
   cd $phantomdir;
   dims="`make --quiet SETUP=$setup getdims`";
   dims=${dims//[^0-9]/};  # number only
   if [ "X$dims" == "X" ]; then
      dims=$maxdim;
   fi
   if (($dims <= $maxdim)); then
      maxp='';
   else
      maxp="MAXP=$maxdim";
      echo "compiling $setup with $maxp";
   fi
   echo "<tr>" >> $htmlfile;
   for target in $listoftargets; do
      ntotal=$((ntotal + 1));
      err=0;
      newwarn=1;
      rm -f warnings.tmp;
      errorlog="./logs/make-$target-errors-$setup-$SYSTEM.txt";
      makeout="/dev/null";
      errorlogold="./logs/make-$target-errors-$setup-$SYSTEM-old.txt";
      if [ -e $errorlog ]; then
         mv $errorlog $errorlogold;
      fi
      colour=$white; # white (default)
      ncheck=$((ncheck + 1));
      printf "Checking $setup ($target)... ";
      if [ "$component"=="setup" ]; then
         rm -f $phantomdir/bin/phantomsetup;
         #make clean >& /dev/null;
      fi
      #case $component in
      #'utils')
      #  make cleanutils >& /dev/null;;
      #*)
      #  make clean >& /dev/null;;
      #esac
      make SETUP=$setup $nolibs $maxp $target 1> $makeout 2> $errorlog; err=$?;
      #--remove line numbers from error log files
      sed -e 's/90(.*)/90/g' -e 's/90:.*:/90/g' $errorlog | grep -v '/tmp' > $errorlog.tmp && mv $errorlog.tmp $errorlog;
      if [ $err -eq 0 ]; then
         echo "OK";
         colour=$green;
         text='OK';
      else
         echo "FAILED"; grep Error $errorlog;
         colour=$red;
         text='**FAILED**';
         nfail=$((nfail + 1));
         echo $setup >> $faillog;
      fi
      if [ -e $errorlogold ]; then
         diff --unchanged-line-format="" --old-line-format="" --new-line-format="%L" $errorlogold $errorlog | tail -20 > warnings.tmp
         if [ -s warnings.tmp ]; then
            newwarn=1;
         else
            newwarn=0;
         fi
      fi
      if [ $newwarn -eq 1 ] && [ $err -eq 0 ]; then
         colour=$amber;
         text='**NEW WARNINGS**';
      fi
      htext='';
      echo "<td bgcolor=\"$colour\">$setup ($target)</td>" >> $htmlfile;
      errors='';
      if [ -e $errorlog ] && [ "X$component" != "Xsetup" ]; then
         grep Error $errorlog > errors.tmp;
         while read line; do
            errors+="$line<br/>";
         done < errors.tmp;
      fi
      if [ -e errors.tmp ]; then
         rm errors.tmp;
      fi
      ref="`basename $errorlog`";
      if [ "X$url" == "X" ]; then
         href=$ref;  # link to local file
      else
         href="$url/$ref"; # link to online file
      fi
      if [ $err -gt 0 ]; then
         echo "<td><a href=\"$href\">error log</a></td><td>$errors</td>" >> $htmlfile;
      else
         if [ $err -lt 0 ]; then
            echo "<td></td><td></td>" >> $htmlfile;
         else
            echo "<td><a href=\"$href\">warnings</a></td><td></td>" >> $htmlfile;
         fi
      fi
      if [ $newwarn -eq 1 ] && [ "X$component" != "Xsetup" ]; then
         hwarn='';
         nwarn=$((nwarn + 1));
         if [ -e warnings.tmp ]; then
            while read line; do
                hwarn+="$line<br/>";
            done < warnings.tmp
         fi
         echo "<td>NEW WARNINGS:<br/>$hwarn</td>" >> $htmlfile;
      else
         echo "<td>$htext</td>" >> $htmlfile;
      fi
      if [ "X$target" == "Xsetup" ] && [ "X$component" == "Xsetup" ]; then
         check_phantomsetup $setup;
      fi
   done
   echo "</tr>" >> $htmlfile;
   cd $pwd;
done
echo "</table>" >> $htmlfile;
echo "<p>Checked $ncheck of $ntotal; <strong>$nfail failures</strong>, $nwarn new warnings</p>" >> $htmlfile;
done
