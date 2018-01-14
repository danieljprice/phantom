#!/bin/bash
#--------------------------------------------------------
#
#  Script to plot nightly code statistics using
#  google charts API. 
#
#  Written by: Daniel Price Jan. 2018
#
#--------------------------------------------------------
pwd=$PWD;
phantomdir="$pwd/../";
if [ ! -e $phantomdir/scripts/$0 ]; then
   echo "Error: This script needs to be run from the phantom/scripts directory";
   exit;
fi
webdir=$dir/web;
#
# generic routines to write code to print google charts
#
print_chart_header()
{
  echo "  google.charts.load('current', {'packages':['line']});";
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
  echo "        width: 900,";
  echo "        height: 500";
  echo "      };"
  echo;
  echo "      var chart = new google.charts.Line(document.getElementById('$tag'));";
  echo;
  echo "      chart.draw(data, google.charts.Line.convertOptions(options));";
}
#
#  routines specific to the graphs we want to make
#
graph_author_stats()
{
   authordata='author_count.txt';
   if [ -e $authordata ]; then
      print_chart_header;
      print_column_header "Number of authors";
      print_data_header;
      while read -r datetag val; do
         print_entry $datetag $val;
      done < $authordata
      print_data_footer;
      print_chart_footer "Number of authors"  "Unique contributors in git" "author_count";
   fi
}
outdir=$webdir/nightly/stats/;
if [ -d $outdir ]; then
   graph_author_stats > $outdir/phantomstats.js;
else
   graph_author_stats;
fi
