#!/bin/bash
#--------------------------------------------------------
#
#  Wrapper script to plot nightly code statistics using
#  google charts API.
#
#  Written by: Daniel Price Jan. 2018 (revised Feb 2019)
#
#--------------------------------------------------------
webdir=$PWD/web;
outdir=$webdir/nightly/stats/;
if [ -d $outdir ]; then
   dir=$outdir;
else
   dir='.';
fi
outfile="$dir/authorcount.js";
echo "writing to $outfile";
./make_google_chart.sh author_count.txt "Number of authors" \
                                        "Unique contributors in git" \
                                        "Number of authors" > $outfile
outfile="$dir/buildstatus.js";
echo "writing to $outfile";
./make_google_chart.sh build_status.txt "Nightly build status" \
                                        "0=buildtestfail, 1=buildfail 2=testfail 3=ok" \
                                        "Build score" > $outfile

outfile="$dir/testbottiming.js";
echo "writing to $outfile";
./make_google_chart.sh testbot_timing.txt "Testbot timings" \
                                       "Walltime used for nightly tests (s)" \
                                       "test(ifort)" "testkd(ifort)" "test(gfortran)" "testkd(gfortran)" "--max=200"> $outfile
outfile="$dir/buildbottiming.js";
echo "writing to $outfile";
./make_google_chart.sh buildbot_timing.txt "Buildbot timings" \
                                       "Time to compile all SETUP=blah options (minutes)" \
                                        "ifort" "gfortran" "--max=300" > $outfile
outfile="$dir/setupcount.js";
echo "writing to $outfile";
./make_google_chart.sh setup_count.txt "Number of unique SETUP= and SYSTEM= options" \
                                       "(fewer is better)" \
                                       "SETUP=" "SYSTEM=" > $outfile
outfile="$dir/codecount.js";
echo "writing to $outfile";
./make_google_chart.sh code_count.txt "Lines of code" \
                                      "Excluding comment and blank lines" \
                                      "main" "setup" "tests" "utils" > $outfile
outfile="$dir/ifdefcount.js";
echo "writing to $outfile";
./make_google_chart.sh ifdef_count.txt "Number of unique #ifdef options" \
                                       "(fewer is better)" \
                                       "number of #ifdef options" > $outfile
outfile="$dir/subcount.js";
echo "writing to $outfile";
./make_google_chart.sh sub_count.txt "Number of modules, subroutines and functions"\
                                     "more subroutines is good" \
                                     "modules" "subroutines" "functions" > $outfile;
