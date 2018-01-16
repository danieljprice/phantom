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
  nauth=`cd $phantomdir; git shortlog -s -n | cut -f 2 | wc -l`;
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
count_code_in_git_history()
{
  cd $phantomdir;
  git log --format="%ad %h" --date=iso > git_history;
  while read -r date mytime tz hash; do
     git checkout $hash
     ncode=$(get_lines_of_code);
     echo $date $mytime $tz $ncode >> code_count.txt;
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
#
# run these in turn
#
nauthors=$(get_author_count);
ncode="$(get_lines_of_code)";
echo "Lines of code: $ncode";
echo "Number of authors: $nauthors";
if [ "X$outdir" != "X" ]; then
   echo "Writing to $outdir/author_count.txt";
   echo $datetag $nauthors >> $outdir/author_count.txt;
   echo "Writing to $outdir/code_count.txt";
   echo $datetagiso $ncode >> $outdir/code_count.txt;
   echo "Writing build status to $outdir/build_status.txt";
   get_build_status_from_git_tags | sort > $outdir/build_status.txt;
fi
