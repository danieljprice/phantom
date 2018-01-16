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
print_entry()
{
  nDate=$(parse_datetag $1);
  echo "        [$nDate, $2],";
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
if [ -d $outdir ]; then
   outfile="$outdir/authorcount.js";
   graph_author_stats > $outfile;
   echo "writing to $outfile";
   outfile="$outdir/buildstatus.js";
   graph_build_status > $outfile;
   outfile="$outdir/testbottiming.js";
   graph_timing_data > $outfile;
else
   graph_author_stats;
   echo; echo; echo;
   graph_build_status;
   echo; echo; echo;
   graph_timing_data;
fi
