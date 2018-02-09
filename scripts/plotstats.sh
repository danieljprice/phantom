#!/bin/bash
#--------------------------------------------------------
#
#  Script to plot nightly code statistics using
#  google charts API. 
#
#  Written by: Daniel Price Jan. 2018
#
#--------------------------------------------------------
webdir=$PWD/web;
outdir=$webdir/nightly/stats/;
#
# generic routines to write code to print google charts
#
print_chart_header()
{
  echo "    google.charts.load('current', {'packages':['line']});";
  echo "    google.charts.setOnLoadCallback(drawChart);";
  echo;
  echo "    function drawChart() {";
  echo;
  echo "      var data = new google.visualization.DataTable();";
  echo "      data.addColumn('date', 'Date');";
}
print_column_header()
{
  echo "      data.addColumn('number', '$1');";
}
print_data_header()
{
  echo "      data.addRows([";
}
parse_datetag()
{
  tag=$1;
  year=${tag:0:4};
  mon=${tag:4:2};
  mon=${mon#0};
  mon=$(( mon - 1 ));
  printf -v month "%02d" $mon;
  day=${tag:6:2};
  echo "new Date($year,$month,$day)";
}
parse_dateiso()
{
  tag=$1;
  year=${tag:0:4};
  mon=${tag:5:2};
  mon=${mon#0};
  mon=$(( mon - 1));
  printf -v month "%02d" $mon;
  day=${tag:8:2};
  # time
  tag=$2;
  hour=${tag:0:2};
  min=${tag:3:2};
  sec=${tag:6:2};
  echo "new Date($year,$month,$day,$hour,$min,$sec)";
}
print_entry()
{
  nDate=$(parse_datetag $1);
  echo "        [$nDate, $2],";
}
print_twocol_entry()
{
  nDate=$(parse_datetag $1);
  echo "        [$nDate, $2, $3],";
}
print_last_entry()
{
  nDate=$(parse_datetag $1);
  echo "        [$nDate, $2]";
}
print_data_footer()
{
  echo "      ]);";
}
print_chart_footer()
{
  echo;
  echo "      var options = {";
  echo "       chart: {";
  echo "         title: '$1',";
  echo "         subtitle: '$2'";
  echo "        },";
  if [ "X$4" != "X" ]; then
     echo "   axes: { ";
     echo "       y: { ";
     echo "          all: {       ";
     echo "             range: {  ";
     echo "                max:$4,";
     echo "                min: 0 ";
     echo "             } ";
     echo "          }";
     echo "       }";
     echo "    },";
  fi
  echo "        width: 900,";
  echo "        height: 500";
  echo "      };"
  echo;
  echo "      var chart = new google.charts.Line(document.getElementById('$3'));";
  echo;
  echo "      chart.draw(data, google.charts.Line.convertOptions(options));";
  echo "    }";
}
#
#  routines specific to the graphs we want to make
#
graph_author_stats()
{
   datafile='author_count.txt';
   if [ -e $datafile ]; then
      print_chart_header;
      print_column_header "Number of authors";
      print_data_header;
      while read -r datetag val; do
         print_entry $datetag $val;
      done < $datafile
      print_data_footer;
      print_chart_footer "Number of authors"  "Unique contributors in git" "author_count";
   else
      echo "ERROR: could not find data file $datafile";
   fi
}
graph_build_status()
{
   datafile='build_status.txt';
   if [ -e $datafile ]; then
      print_chart_header;
      print_column_header "Build score";
      print_data_header;
      while read -r datetag val; do
         print_entry $datetag $val;
      done < $datafile
      print_data_footer;
      print_chart_footer "Nightly build status"  "0=buildtestfail, 1=buildfail 2=testfail 3=ok" "build_status";
   else
      echo "ERROR: could not find data file $datafile";
   fi
}
graph_timing_data()
{
   datafile='testbot_timing.txt'
   if [ -e $datafile ]; then
      print_chart_header;
      print_column_header "test(ifort)";
      print_column_header "testkd(ifort)";
      print_column_header "test(gfortran)";
      print_column_header "testkd(gfortran)";
      print_data_header
      cat $datafile;
      print_data_footer;
      print_chart_footer "Testbot timings" "Walltime used for nightly tests (s)" "testbot_timing" "200";
   fi
}
graph_buildbot_data()
{
   datafile='buildbot_timing.txt'
   if [ -e $datafile ]; then
      print_chart_header;
      print_column_header "ifort";
      print_column_header "gfortran";
      print_data_header
      cat $datafile;
      print_data_footer;
      print_chart_footer "Buildbot timings" "Time to compile all SETUP=blah options (minutes)" "buildbot_timing" "300";
   fi
}
graph_code_count()
{
   datafile='code_count.txt'
   if [ -e $datafile ]; then
      print_chart_header;
     # print_column_header "total"
      print_column_header "main"
      print_column_header "setup"
      print_column_header "tests"
      print_column_header "utils"
      print_data_header
      while read -r mydate mytime tz nmain nsetup ntests nutils; do
          #ntotal=$(( nmain + nsetup + ntests + nutils ));
          echo "    [$(parse_dateiso $mydate $mytime $tz),$nmain,$nsetup,$ntests,$nutils],";
      done < $datafile
      print_data_footer
      print_chart_footer "Lines of code" "Excluding comment and blank lines" "code_count"
   else
      echo "ERROR: could not find data file $datafile";
   fi
}
graph_ifdef_count()
{
   datafile='ifdef_count.txt'
   if [ -e $datafile ]; then
      print_chart_header;
      print_column_header "number of #ifdef options";
      print_data_header;
      while read -r datetag val; do
         print_entry $datetag $val;
      done < $datafile
      print_data_footer;
      print_chart_footer "Number of unique #ifdef options"  "(fewer is better)" "ifdef_count";
   else
      echo "ERROR: could not find data file $datafile";
   fi
}
graph_setup_count()
{
   datafile='setup_count.txt'
   if [ -e $datafile ]; then
      print_chart_header;
      print_column_header "SETUP=";
      print_column_header "SYSTEM=";
      print_data_header;
      while read -r datetag vala valb; do
         print_twocol_entry $datetag $vala $valb;
      done < $datafile
      print_data_footer;
      print_chart_footer "Number of unique SETUP= and SYSTEM= options"  "(fewer is better)" "setup_count";
   else
      echo "ERROR: could not find data file $datafile";
   fi
}
graph_sub_count()
{
   datafile='sub_count.txt'
   if [ -e $datafile ]; then
      print_chart_header;
     # print_column_header "total"
      print_column_header "modules"
      print_column_header "subroutines"
      print_column_header "functions"
      print_data_header
      while read -r mydate mytime tz nmod nsub nfunc; do
          echo "    [$(parse_dateiso $mydate $mytime $tz),$nmod,$nsub,$nfunc],";
      done < $datafile
      print_data_footer
      print_chart_footer "Number of modules, subroutines and functions" "more subroutines is good" "sub_count"
   else
      echo "ERROR: could not find data file $datafile";
   fi
}
#
# main script
#
if [ -d $outdir ]; then
   dir=$outdir;
else
   dir='.';
fi
outfile="$dir/authorcount.js";
echo "writing to $outfile";
graph_author_stats > $outfile;

outfile="$dir/buildstatus.js";
echo "writing to $outfile";
graph_build_status > $outfile;

outfile="$dir/testbottiming.js";
echo "writing to $outfile";
graph_timing_data > $outfile;

outfile="$dir/buildbottiming.js";
echo "writing to $outfile";
graph_buildbot_data > $outfile;

outfile="$dir/setupcount.js";
echo "writing to $outfile";
graph_setup_count > $outfile;

outfile="$dir/codecount.js";
echo "writing to $outfile";
graph_code_count > $outfile;

outfile="$dir/ifdefcount.js";
echo "writing to $outfile";
graph_ifdef_count > $outfile;

outfile="$dir/subcount.js";
echo "writing to $outfile";
graph_sub_count > $outfile;
