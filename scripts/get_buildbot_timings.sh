#
# figure out how long the buildbot took from the output .html files
#
convert_date()
{
   local UNAME=$(uname);
   if [[ "$UNAME" == "Darwin" ]]; then
      date -j -f "%a %b %d %T %Z %Y" "$1" "+%s";
   else
      date -d "$1" "+%s";
   fi
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
get_buildbot_timings()
{
  file=$1;
  if [ ${file:0:2} == "20" ]; then # check filename starts with year [20xx]
     dtag=${file/.html/};
     grep 'Build checked' $file | sed 's/.*checked: //' > dates.tmp
     t1=$(convert_date "`sed -n 1p dates.tmp`");
     t2=$(convert_date "`sed -n 2p dates.tmp`");
     t3=$(convert_date "`sed -n 3p dates.tmp`");
     t4=$(convert_date "`sed -n 4p dates.tmp`");
     timemsg=$(( (t2 - t1)/60 ));
     timegfc=$(( (t4 - t3)/60 ));
     #echo $timemsg $timegfc;
     if [ $timemsg -ge 0 ]; then
        echo "     [$(parse_datetag $dtag),$timemsg,$timegfc],";
     fi
     rm dates.tmp;
  fi
}
for x in "$@"; do
    get_buildbot_timings $x;
done
