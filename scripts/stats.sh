#!/bin/bash
#
# The Phantom stats-bot
#
# Collects and prints statistics about code
# e.g. lines of code, number of authors
#
# Written by Daniel Price, 2018-, daniel.price@monash.edu
#
pwd=$PWD;
phantomdir="$pwd/../";
outdir=$1;
if [ ! -e $phantomdir/scripts/$0 ]; then
   echo "Error: This script needs to be run from the phantom/scripts directory";
   exit;
fi
datetag=`date "+%Y%m%d"`;
datetagiso=`date "+%Y-%m-%d %H:%M:%S %z"`;
#-------------
# subroutines
#-------------
get_author_count()
{
  nauth=`cd $phantomdir; git shortlog -s -n HEAD | wc -l`;
  echo $nauth;
}
count_code()
{
  dir=$1;
  ncode=`sed '/^!/d;/^\s*$/d;' $dir/*.*90 | wc -l`;
  #ntotal=`cat *.*90 | wc -l`;
  #ncomments=$(( ntotal - ncode ));
  echo "$ncode";
}
count_matches()
{
  n=`cd $phantomdir; grep "$1" src/*/*.*90 | cut -d':' -f 2 | wc -l`;
  echo "$n";
}
count_unique_matches()
{
  n=`cd $phantomdir; grep "$1" src/*/*.*90 | cut -d':' -f 2 | sort -u | wc -l`;
  echo "$n";
}
get_subroutine_count()
{
  nsub=$(count_matches 'end subroutine');
  nmod=$(count_matches 'end module');
  nfunc=$(count_matches 'end function');
  echo "$nmod $nsub $nfunc";
}
get_lines_of_code()
{
  str='';
  for dir in main setup tests utils; do
     str+="$(count_code $phantomdir/src/$dir) ";
  done
  echo "$str";
}
#
# subroutine below can be used to reconstruct code count
# from entire git history. Not used by default as
# we simply append to this file each night (but can
# be used to reconstruct file if accidentally deleted)
#
remake_stats_from_git_history()
{
  cd $phantomdir;
  git log --format="%ad %h" --date=iso > git_history;
  while read -r date mytime tz hash; do
     git checkout $hash
     ncode=$(get_lines_of_code);
     echo $date $mytime $tz $ncode >> code_count.txt;
     subcount="$(get_subroutine_count)";
     echo $date $mytime $tz $subcount >> sub_count.txt;
     # use simpler date format
     mydate=${date/-/};
     mydate=${mydate/-/};
     nauthors=$(get_author_count);
     echo $mydate $nauthors >> author_count.txt;
     nifdef="$(count_unique_matches '#ifdef')";
     echo $mydate $nifdef >> ifdef_count.txt;
  done < git_history
  rm git_history;
  git checkout master;
}
#
# reconstruct build status from the ok-20180101 git tags
#
get_build_status_from_git_tags()
{
  for x in `cd $phantomdir; git tag`; do
      dtag=${x##*-};
      status=${x%%-*};
      case $status in
      'buildtestfail')
         score=0;;
      'buildfail')
         score=1;;
      'testfail')
         score=2;;
      'ok')
         score=3;;
      *) 
         score=-1;;
      esac
      if [ $score -ge 0 ]; then
         echo $dtag $score;
      fi
  done
}
# uncomment the following to recreate stats from entire git history
# otherwise we just give instant stats
#remake_stats_from_git_history;
#exit;
#
# run these in turn
#
nauthors=$(get_author_count);
ncode="$(get_lines_of_code)";
nifdef="$(count_unique_matches '#ifdef')";
subcount="$(get_subroutine_count)";
echo "Lines of code: $ncode";
echo "Number of modules, subroutines, functions: $subcount";
echo "Number of #ifdef statements : $nifdef";
echo "Number of authors           : $nauthors";
if [ "X$outdir" != "X" ]; then
   echo "Writing to $outdir/author_count.txt";
   echo $datetag $nauthors >> $outdir/author_count.txt;
   echo "Writing to $outdir/code_count.txt";
   echo $datetagiso $ncode >> $outdir/code_count.txt;
   echo "Writing to $outdir/ifdef_count.txt";
   echo $datetag $nifdef >> $outdir/ifdef_count.txt;
   echo "Writing to $outdir/sub_count.txt";
   echo $datetagiso $subcount >> $outdir/sub_count.txt;
   echo "Writing build status to $outdir/build_status.txt";
   get_build_status_from_git_tags | sort > $outdir/build_status.txt;
fi
