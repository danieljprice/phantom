#!/bin/bash
#
# script used to create phantom official release tarballs
# should have copy of phantom checked out in phantom/
# and copy of phantom website checked out in web-phantom/
#
vernum=1.0
datetag=`date "+%Y%m%d"`;
distfile="phantom-v${vernum}.tar.gz"
urlgitrepo='https://bitbucket.org/danielprice/phantom'
phantomdir="$PWD/phantom"
EDITOR='nedit'
if [ ! -e $phantomdir ]; then
   echo "Usage: phantom_release.sh [run from dir with phantom/ as a subdirectory]"
   echo "e.g.:"; echo; echo "$ cd"; echo "$ ~/phantom/scripts/phantom_release.sh"; echo;
   exit;
fi
webdir="$PWD/web-phantom"
if [ ! -e $webdir ]; then
   echo "no directory web-phantom; require copy of phantom website"
   checkout_website
fi
checkout_website () {
   echo "cloning website"
   git clone https://bitbucket.org/phantomsph/phantomsph.bitbucket.org web-phantom
}
make_distfiles()
{
   cd $phantomdir;
   echo "creating $distfile"
   git archive master --format=tar.gz --prefix="phantom/" > $distfile;
   git log --oneline > ChangeLog;
}
get_distfile_sha()
{
   shachecksum=`openssl sha256 $distfile | cut -d'=' -f 2`;
   shafile=$distfile.sha256.txt;
   echo "sha256: $shachecksum";
   echo $shachecksum > $shafile;
}
update_website()
{
   cd $webdir;
   echo "$distfile -> $webdir/releases"
   cp $phantomdir/$distfile releases/
   echo "$shafile -> $webdir/releases"
   cp $phantomdir/$shafile releases/
   cp $phantomdir/ChangeLog releases/
   today=`date +"%d/%m/%Y"`
   gitsha=`cd $phantomdir; git rev-parse --short HEAD`;
   echo "<tr><td>$today</td><td><a href=\"releases/$distfile\">$distfile</a></td><td>16Mb</td><td><a href=\"releases/$shafile\">sha256</a></td><td><a href=\"$urlgitrepo/changeset/$gitsha\">$gitsha</a></td><td><a href=\"releases/releasenotes.txt\">release notes</a></td></tr>" >> releases.html
   if [ "X$EDITOR" != "X" ]; then
      $EDITOR releases/releasenotes.txt &
      $EDITOR index.html releases.html &
   fi
}
push_website()
{   
   cd $webdir
   git pull
   cd $webdir/releases/
   $phantomdir/scripts/indexhtml.pl
   git add $distfile $shafile releasenotes.txt ChangeLog index.html
   git commit -m "[deploy] version $vernum"
   git push
}
tag_code()
{
   cd $phantomdir;
   gittag="v${vernum}-$datetag"
   git tag $gittag;
   git push --tags;
}
stage_release()
{ 
   make_distfiles
   get_distfile_sha
   update_website
}
deploy_release()
{
   push_website
   tag_code
}
stage_release
if [ "$1" == "deploy" ]; then
   echo "deploy release"
   deploy_release  # only do this step when you are sure
fi
