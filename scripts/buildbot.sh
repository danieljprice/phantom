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
# Written by Daniel Price, 2012-2015, daniel.price@monash.edu
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
#
# change the line below to exclude things that depend on external libraries from the build
#
nolibs='MESAEOS=no'
if [ $# -gt 1 ]; then
   url=$2;
else
   url='';
fi
if [ ! -e $phantomdir/scripts/$0 ]; then
   echo "Error: This script needs to be run from the phantom/scripts directory";
   exit;
fi
echo "url = $url";
htmlfile="$phantomdir/logs/build-status-$SYSTEM.html";
faillog="$phantomdir/logs/build-failures-$SYSTEM.txt";
if [ -e $htmlfile ]; then
   rm $htmlfile;
fi
if [ -e $faillog ]; then
   rm $faillog;
fi
listofcomponents='main utils';
for component in $listofcomponents; do
case $component in
 'utils')
   text=$component;
   listofsetups='test';
   listoftargets='utils';;
 *)
   text='';
   listofsetups=`grep 'ifeq ($(SETUP)' $phantomdir/build/Makefile | cut -d, -f 2 | cut -d')' -f 1 | sed '/tracers/d' | sed '/mcfost/d'`
   listoftargets='phantom setup analysis moddump';;
esac
echo "<h2>Checking Phantom $text build, SYSTEM=$SYSTEM</h2>" >> $htmlfile;
echo "Build checked: "`date` >> $htmlfile;
ncheck=0;
ntotal=0;
nfail=0;
nwarn=0;
#
# write html header
#
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
      colour="#FFFFFF"; # white (default)
      ncheck=$((ncheck + 1));
      printf "Checking $setup ($target)... ";
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
         colour="#009900";  # green
         text='OK';
      else
         echo "FAILED"; grep Error $errorlog;
         colour="#FF0000";  # red
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
         colour="#FF6600"; # amber
         text='**NEW WARNINGS**';
      fi
      htext='';
      echo "<td bgcolor=\"$colour\">$setup ($target)</td>" >> $htmlfile;
      errors='';
      if [ -e $errorlog ]; then
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
      if [ $newwarn -eq 1 ]; then
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
   done
   echo "</tr>" >> $htmlfile;
   cd $pwd;
done
echo "</table>" >> $htmlfile;
echo "<p>Checked $ncheck of $ntotal; <strong>$nfail failures</strong>, $nwarn new warnings</p>" >> $htmlfile;
done
