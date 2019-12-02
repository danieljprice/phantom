#!/bin/bash
#
# The Phantom test-bot
#
# This is a wrapper for the Phantom test suite ("make test")
# Runs all possible tests and collates the results into
# tables for the html  pages.
#
# Also checks for build failures in the test suite
#
# You can run this script yourself from the scripts directory as follows:
# cd phantom/scripts; ./testbot.sh
#
# Outputs are put in the phantom/logs directory
#
# Written by Daniel Price, 2013-2017, daniel.price@monash.edu
#
if [ X$SYSTEM == X ]; then
   echo "Error: Need SYSTEM environment variable set to run PHANTOM tests";
   echo "Usage: $0 [url]";
   exit;
fi
pwd=$PWD;
phantomdir="$pwd/../";
if [ $# -gt 0 ]; then
   url=$1;
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
htmlfile="$phantomdir/logs/test-status-$SYSTEM.html";
faillog="$phantomdir/logs/test-failures-$SYSTEM.txt";
if [ -e $htmlfile ]; then
   rm $htmlfile;
fi
if [ -e $faillog ]; then
   rm $faillog;
fi
white="#FFFFFF";
red="#FF0000";
amber="#FF6600";
green="#009900";
for debug in no yes; do # perform each test with and without debug flags
if [ $debug == "yes" ]; then
   tag='-debug';
   debugflag=' DEBUG=yes';
   echo $debugflag;
else
   tag='';
   debugflag='';
fi
echo "<h2>Checking Phantom testsuite, SYSTEM=$SYSTEM$debugflag</h2>" >> $htmlfile;
echo "Test suite checked: "`date` >> $htmlfile;
echo "<table>" >> $htmlfile;
ncheck=0;
nfail=0;
listofsetups="test testkd test2 testcyl testgrav testdust testnimhd testgrowth";
for setup in $listofsetups; do
    cd $phantomdir;
    errorlog="./logs/test-results-$setup-$SYSTEM$tag.txt";
    errorlogold="./logs/test-results-$setup-$SYSTEM$tag-old.txt";
    if [ -e $errorlog ]; then
       mv $errorlog $errorlogold;
    fi
    colour=$white;
    result='failed to run';
    timing='';
    printf "Checking $setup...";
    ncheck=$(( ncheck + 1 ));
    err=0;
    herr='';
    if [ -e errors.tmp ]; then
       rm errors.tmp;
    fi
    if [ $setup == "testgrav" ]; then
       arg="gravity ptmass";
    elif [ $setup == "testdust" ]; then
       arg="dust";
    elif [ $setup == "testnimhd" ]; then
       arg="nimhd";
    elif [ $setup == "testgrowth" ]; then
       arg="dustgrowth";
    else
       arg="";
    fi
    make SETUP=$setup $debugflag >& $errorlog; err=$?;
    if [ $err -eq 0 ]; then
       # build was OK, proceed to run test suite
       ./bin/phantom test $arg >& $errorlog;
       if [ -s $errorlog ]; then
          grep 'TEST SUITE PASSED' $errorlog >& /dev/null; pass=$?;
          passes=`grep 'PASSED:' $errorlog`;
          fails=`grep 'FAILED:' $errorlog`;
          errors='';
          grep 'PASSED' $errorlog >& /dev/null; runerr=$?;
          if [ $pass -eq 0 ]; then
             text=$passes;
             result="$passes"
             colour=$green;
             tcpu=`grep 'total cpu time' $errorlog | cut -d= -f 2`;
             tcpu=${tcpu/(/};
             twall=`grep 'total wall time' $errorlog | cut -d= -f 2`;
             twall=${twall/(/};
             timing="completed in $twall ($tcpu cpu)";
          elif [ $runerr -eq 0 ]; then
          #  "successfully" failed, i.e. completed but with test failures
             grep -B 5 ^' FAILED' $errorlog > errors.tmp;
             if [ -e errors.tmp ]; then
                while read line; do
                  herr+="$line<br/>";
                done < errors.tmp
             fi
             result="$fails $passes"
             text=$fails;
             colour=$red;
             nfail=$(( nfail + nerr ));
          else
             text="FAILED TO COMPLETE";
             result=$text;
             fails="1";
             colour=$amber;
          fi
          tmp=${fails/FAILED:/};
          nerr=${tmp/of*/};
          nfail=$(( nfail + nerr ));
          echo $text;
       else
          colour=$amber;
          text='***FAILED TO RUN***';
          nfail=$(( nfail + 1 ));
          echo $text;
       fi
    else
       result='failed to build';
       colour=$red;
       text='***FAILED TO BUILD***';
       nfail=$(( nfail + 1 ));
    fi
    if [ $nfail -gt 0 ]; then
       echo "$setup $result" >> $faillog;
    fi
    echo "<tr><td bgcolor=$colour>$setup</td>" >> $htmlfile;
    ref="`basename $errorlog`";
    if [ "X$url" == "X" ]; then
       href=$ref;  # link to local file
    else
       href="$url$ref"; # link to file on web server
    fi
    if [ -e errors.tmp ]; then
       echo "<td>$fails<br/>$passes</td><td>$herr</td><td></td><td><a href=\"$href\">show details</a></td></tr>" >> $htmlfile;
    else
       echo "<td>$result</td><td>$timing</td><td>$herr</td><td><a href=\"$href\">show details</a></td></tr>" >> $htmlfile;
    fi
    cd $pwd;
done
echo "</table>" >> $htmlfile;
echo "<p>Performed $ncheck test runs; <strong>$nfail failures</strong></p>" >> $htmlfile;
done
