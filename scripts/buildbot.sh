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
# Written by Daniel Price, 2012-2023, daniel.price@monash.edu
#
if [ X$SYSTEM == X ]; then
   echo "Error: Need SYSTEM environment variable set to check PHANTOM build";
   echo "Usage: $0 [max idim to check] [url]";
   exit;
fi

# Default arguments
maxdim=11000000;
url='';
batch=1;
nbatch=1;

while [[ "$1" == --* ]]; do
  case $1 in
    --maxdim)
      shift;
      maxdim=$1; # max idim to check
      ;;

    --url)
      shift;
      url=$1; # url for results
      ;;

    --parallel)
      shift;
      batch=$1; # the batch number being run
      shift;
      nbatch=$1; # total number of batches to divide work into
      # Example:
      # --parallel 1 10
      # will divide the tests into 10 batches, and "1" specifies that the first batch is being run
      ;;

    *)
      badflag=$1
      ;;
  esac
  shift
done

if [[ "$badflag" != "" ]]; then
   echo "ERROR: Unknown flag $badflag"
   exit
fi


if (($maxdim > 0)) && (($maxdim < 2000000000)); then
   echo "Using maxdim = $maxdim";
else
   echo "Usage: $0 [max idim to check] [url]";
   exit;
fi

echo "url = $url";

pwd=$PWD;
phantomdir="$pwd/../";
listofcomponents='main setup analysis utils';
#listofcomponents='setup'
#
# get list of targets, components and setups to check
#
allsetups=`grep 'ifeq ($(SETUP)' $phantomdir/build/Makefile_setups | grep -v skip | cut -d, -f 2 | cut -d')' -f 1`
#allsetups='star'
setuparr=($allsetups)
batchsize=$(( ${#setuparr[@]} / $nbatch + 1 ))
offset=$(( ($batch-1) * $batchsize ))
allsetups=${setuparr[@]:$offset:$batchsize}
echo "Batch ${batch}: ${allsetups}"
#
# change the line below to exclude things that depend on external libraries from the build
#
nolibs='MESAEOS=no'
#
# disallow compiler warnings?
#
nowarn='NOWARN=yes'
#
# check we are running the script from the scripts directory
#
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
htmlfile_all="$phantomdir/logs/build-status-$SYSTEM.html";
faillog="$phantomdir/logs/build-failures-$SYSTEM.txt";
faillogsetup="$phantomdir/logs/setup-failures-$SYSTEM.txt";
failloganalysis="$phantomdir/logs/analysis-failures-$SYSTEM.txt";
#
# delete old log files
#
if [ -e $htmlfile_all ]; then
   rm $htmlfile_all;
fi
if [ -e $faillog ]; then
   rm $faillog;
fi
if [ -e $faillogsetup ]; then
   rm $faillogsetup;
fi
if [ -e $failloganalysis ]; then
   rm $failloganalysis;
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
  echo $text;
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
      print_result "phantomsetup exists" $pass;
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
   #
   # run ./phantomsetup prefix, answering any questions with "Enter"
   #
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
   #
   # run phantomsetup up to 3 times to successfully create/rewrite the .setup file
   #
   infile="${prefix}.in"
   ./phantomsetup $prefix < myinput.txt > /dev/null;
   ./phantomsetup $prefix < myinput.txt > /dev/null;
   if [ -e "$prefix.setup" ]; then
      print_result "creates .setup file" $pass;
      #test_setupfile_options "$prefix" "$prefix.setup" $infile;
   else
      print_result "no .setup file" $warn;
   fi
   if [ -e "$infile" ]; then
      print_result "creates .in file" $pass;
      #
      # if creating the .in file succeeds, try to run phantom
      # on the .in file with nmax=0
      #
      if [ -e $phantomdir/bin/phantom ]; then
         print_result "phantom exists" $pass;
         cp $phantomdir/bin/phantom .;
         #
         # set nmax=0 in the .in file
         #
         sed 's/nmax = .*/nmax = 0/g' ${infile} > ${infile}.bak
         mv ${infile}.bak ${infile}
         #
         # run phantom on the .in file
         #
         ./phantom $infile > /dev/null; err=$?;
         if [ $err -eq 0 ]; then
            print_result "./phantom $infile runs ok" $pass;
            dumpfile="${prefix}_00000"
            if [ -s $dumpfile ]; then
               print_result "${dumpfile} successfully created" $pass;
            else
               print_result "FAIL: did not create ${dumpfile}" $fail;
               myfail=$(( myfail + 1 ));
            fi
         else
            print_result "FAIL: ./phantom $infile fails with error" $fail;
            myfail=$(( myfail + 1 ));
         fi
      else
         print_result "FAIL: phantom does not exist" $fail;
         myfail=$(( myfail + 1 ));
      fi
   else
      print_result "FAILED to create .in file after 3 attempts" $fail;
      myfail=$(( myfail + 1 ));
   fi
   if [ $myfail -gt 0 ]; then
      echo $setup >> $faillogsetup;
   fi
}
#
# check that all possible values of certain
# variables in the .setup file work
#
test_setupfile_options()
{
   myfail=0;
   setup=$1;
   setupfile=$2;
   infile=$3;
   range=''
   if [ "X$setup"=="Xstar" ]; then
      param='iprofile'
      range='1 2 3 4 5 6 7'
   fi
   for x in $range; do
       valstring="$param = $x"
       echo "checking $valstring"
       sed "s/$param.*=.*$/$valstring/" $setupfile > ${setupfile}.tmp
       cp ${setupfile}.tmp $setupfile
       rm $infile
       ./phantomsetup $setupfile < /dev/null > /dev/null;
       ./phantomsetup $setupfile < /dev/null;

       if [ -e $infile ]; then
          print_result "successful phantomsetup with $valstring" $pass;
       else
          print_result "FAIL: failed to create .in file with $valstring" $fail;
          myfail=$(( myfail + 1 ));
          echo $setup $valstring >> $faillogsetup;
       fi
   done
}
#
# unit tests for phantomanalysis utility
# (currently only exist for SETUP=star)
#
check_phantomanalysis ()
{
   myfail=0;
   setup=$1;
   if [ -e ./bin/phantomanalysis ]; then
      print_result "exists" $pass;
   else
      print_result "FAIL: phantomanalysis does not exist" $fail;
      myfail=$(( myfail + 1 ));
   fi
   dirname="test-phantomanalysis";
   cd /tmp/;
   if [ -d $dirname ]; then
      # do not wipe entire directory to avoid repeatedly downloading data files
      rm -f $dirname/phantomanalysis;
      rm -f $dirname/*.ev $dirname/*.txt;
   else
      mkdir $dirname;
   fi
   cd /tmp/$dirname;
   cp $phantomdir/bin/phantomanalysis .;
   if [ "X$setup" == "Xstar" ]; then
      echo "performing analysis unit tests for SETUP=$setup"
      $pwd/test_analysis_ce.sh; err=$?;
      python $pwd/test_analysis_ce.py; err=$?;
   else
      #echo "there are no analysis unit tests for SETUP=$setup"
      err=0;
   fi
   if [ $err -eq 0 ]; then
      print_result "runs and passes analysis tests" $pass;
   else
      print_result "FAIL: did not pass phantomanalysis tests" $fail;
      myfail=$(( myfail + 1 ));
   fi
   if [ $myfail -gt 0 ]; then
      echo $setup >> $failloganalysis;
   fi
}

listofsetups=$allsetups;
firstsetup=${listofsetups%% *};
lastsetup=${listofsetups##* };
#
#--loop over each setup
#
for setup in $listofsetups; do
   #
   #--loop over each component (phantom,phantomsetup,phantomanalysis,etc)
   #
   for component in $listofcomponents; do
   case $component in
    'analysis')
      text="$component runs, creates output files successfully";
      if [ "$setup" != "star" ]; then continue; fi # skip unless SETUP=star
      listoftargets='analysis';;
    'setup')
      text="$component runs, creates .setup and .in files with no unspecified user input";
      listoftargets='setup';;
    'utils')
      text="$component build";
      if [ "$setup" != "test" ]; then continue; fi # skip unless SETUP=test
      listoftargets='utils';;
    *)
      text='build';
      listoftargets='phantom setup analysis moddump phantomtest';;
   esac

   # output needed for github actions
   echo "::group:: make SETUP=${setup} checking phantom ${text}";

   htmlfile=${htmlfile_all/.html/-${component}.html}
   #
   # write html header
   #
   if [[ "$setup" == "$firstsetup" ]]; then
      if [ -e $htmlfile ]; then
         rm $htmlfile;
      fi
      #echo "opening $htmlfile for output"
      echo "<h2>Checking Phantom $text, SYSTEM=$SYSTEM</h2>" >> $htmlfile;
      echo "Build checked: "`date` >> $htmlfile;
      echo "<table>" >> $htmlfile;
   fi
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

      mydebug=''
      printf "Checking $setup ($target)... ";
      if [[ "$component" == "setup" ]]; then
         rm -f $phantomdir/bin/phantomsetup;
         mydebug='DEBUG=yes' # compile phantomsetup with DEBUG=yes for setup test
         #make clean >& /dev/null;
      fi
      if [[ "$setup" == "blob" ]]; then
         mynowarn='';
         echo "allowing warnings for SETUP=blob"
      else
         mynowarn=$nowarn;
      fi
      make SETUP=$setup $nolibs $mynowarn $maxp $target $mydebug 1> $makeout 2> $errorlog; err=$?;
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
         # also build phantom main binary
         echo "compiling phantom with SETUP=$setup"
         make SETUP=$setup $nolibs $mynowarn $maxp $mydebug 1>> $makeout 2>> $errorlog; err=$?;
         check_phantomsetup $setup;
      elif [ "X$target" == "Xanalysis" ] && [ "X$component" == "Xanalysis" ]; then
         check_phantomanalysis $setup;
      fi
   done
   echo "</tr>" >> $htmlfile;
   cd $pwd;
   if [ "$GITHUB_ACTIONS" == "true" ]; then echo "::endgroup::"; fi
#
# close html file on last setup
#
   if [[ "$setup" == "$lastsetup" ]]; then
      #echo "closing $htmlfile"
      echo "</table>" >> $htmlfile;
      echo "<p>Checked $ncheck of $ntotal; <strong>$nfail failures</strong>, $nwarn new warnings</p>" >> $htmlfile;
      cat $htmlfile >> $htmlfile_all;
   fi
done
done
if [ "$nfail" -gt 0 ] && [ "$RETURN_ERR" == "yes" ]; then exit 1; fi
