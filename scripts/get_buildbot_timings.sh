#
# figure out how long the buildbot took from the output .html files
#
convert_date()
{
   date -j -f "%a %b %d %T %Z %Y" "$1" "+%s";
}
get_buildbot_timings()
{
  file=$1;
  if [ ${file:0:2} == "20" ]; then # check filename starts with year [20xx]
     dtag=${file/.html/};
     year=${dtag:0:4};
     mon=${dtag:4:2};
     day=${dtag:6:2};
     timediff=2;
     grep 'Build checked' $file | sed 's/.*checked: //' > dates.tmp
     t1=$(convert_date "`sed -n 1p dates.tmp`");
     t2=$(convert_date "`sed -n 2p dates.tmp`");
     t3=$(convert_date "`sed -n 3p dates.tmp`");
     t4=$(convert_date "`sed -n 4p dates.tmp`");
     timemsg=$(( (t2 - t1)/60 ));
     timegfc=$(( (t4 - t3)/60 ));
     #echo "GOT t3=$t3 t4=$t4 tgfc=$timegfc"
     #echo $timemsg $timegfc;
     if [ $timemsg -ge 0 ]; then
        echo "     [new Date($year,$(( mon - 1 )),$day),$timemsg,$timegfc],";
     fi
     rm dates.tmp;
  fi
}
for x in "$@"; do
    get_buildbot_timings $x;
done
