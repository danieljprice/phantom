#!/bin/bash
#--------------------------------------------------------
#
#  Script to convert text data files to plots via
#  google charts API.
#
#--------------------------------------------------------
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
# figure out which date format the file uses for time entries
get_date_format()
{
  format=0;
  string=$1;
  if [[ $string == ????-??-?? ]]; then
     format=1;
  elif [[ $string == *"new"* ]]; then
     format=2;
  fi
  echo $format;
}
#
# wrapper routine
#
if [ $# -lt 3 ]; then
   echo "Usage: $0 datafile title subtitle label1 label2 label3..."
else
   file=$1;
   max=0;
   if [ -e $file ]; then
      title=$2;
      subtitle=$3;
      mytag=${file%.*}
      print_chart_header;
      # print column headers
      i=0;
      for arg in "$@"; do
          i=$(( i + 1 ));
          if [ $i -gt 3 ]; then
             if [[ $arg == '--max='* ]]; then
                max=${arg#--max=};
                #echo "MAX=$max";
             else
                print_column_header "$arg";
             fi
          fi
      done
      dateline=`tail -1 $file`
      mydateformat=`get_date_format $dateline`
      print_data_header;
      #echo "DATE FORMAT IS $mydateformat";
      case $mydateformat in
      2)
          cat $file;;
      1)
          while read -r mydate mytime tz val; do
             val=`echo ${val} | sed 's/ /,/g'`; # add commas
             echo "    [$(parse_dateiso $mydate $mytime $tz),$val],";
          done < $file;;
      *)
         while read -r datetag val; do
            val=`echo ${val} | sed 's/ /, /g'`; # add commas
            echo "        [$(parse_datetag $datetag), $val],";
            #print_entry $datetag $val;
         done < $file;;
      esac
      print_data_footer
      if [ $max -gt 0 ]; then
         print_chart_footer "$title"  "$subtitle" "$mytag" "$max";
      else
         print_chart_footer "$title"  "$subtitle" "$mytag";
      fi
   else
      echo "ERROR: could not find data file $file";
   fi
fi
