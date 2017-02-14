#!/bin/bash
#
#  Wrapper script that runs the nightly checks on the code.
#  Interfaces with the central git repository.
#
#  This script in turn runs the buildbot (buildbot.sh)
#  and the testbot (testbot.sh) and emails the results
#  to the relevant users.
#
#  Written by Daniel Price, daniel.price@monash.edu
#
# script settings
dir=$PWD
wikidir=$dir/wiki
codedir=$dir/phantom
wiki=$wikidir/nightly/Build.wiki
url="https://bitbucket.org/danielprice/phantom"
admin="daniel.price@monash.edu";
systems="msg gfortran";
mailfile="$dir/mail.tmp";
mailtmp="$dir/mailcontent.html";
incoming="$dir/incoming.txt";
tagfile="$dir/gittag.tmp";
gittag='';
summary=[];
changesets='';
mailto='';
names='';
sendmail=1;
# work out who to blame
get_incoming_changes ()
{
   echo "--- getting incoming changes ---";
   git fetch
   git log ..@{u} --format="commit:%h%nauthor: %aN <%aE>%nsummary:%s" > $incoming;
   # if no changes, do not run build checks
   if [ -s $incoming ]; then
      echo "got incoming changes, running nightly checks...";
   else
      echo "no changes found, nightly checks not run";
      exit;
   fi
   git log ..@{u} --format="%aN <%aE>" | sort -u > $dir/users.list;
   # extract list of people to send mail to
   mailto='';
   names='';
   while read line; do
      mailto+="$line,";
      names+="`echo $line | cut -d' ' -f 1`, ";
   done < $dir/users.list
   echo "mailto: $mailto";
   # get list of incoming changesets
   changesets=`git log ..@{u} --format="%h"`;
#   changesets=`grep commit $incoming | cut -d':' -f 2`;
   git log ..@{u} --format="%s" > $dir/summary.txt;
   i=0;
   while read line; do
      summary[$((i++))]=`echo $line`
   done < $dir/summary.txt
}
pull_changes ()
{
   # pull the changes
   echo "--- pulling changes ---";
   echo `date` >> $dir/gitpull.log;
   git pull >> $dir/gitpull.log; err=$?;
   if [ $err -ne 0 ]; then
      echo "error pulling changes, buildbot not run";
      exit;
   fi
}
run_buildbot ()
{
   # run the buildbot and testbot scripts
   echo "--- running buildbot ---";
   cd $codedir/scripts
   for sys in $systems; do
      export SYSTEM=$sys;
      echo "SYSTEM=$SYSTEM";
      ./testbot.sh "$url/wiki/nightly/";
      ./buildbot.sh 17000000 "$url/wiki/nightly/";
   done
}
pull_wiki ()
{
   # pull existing wiki files from server
   echo "--- pulling wiki changes ---";
   cd $wikidir
   git pull
}
write_wiki_tag_and_mail_files ()
{
   echo "--- writing wiki and mail files ---"
   # write header for nightly build results
   cd $dir;
   # initialise html file
   echo "<h2>Changesets tested since last build</h2>" > $mailtmp
   echo "<table>" >> $mailtmp
   # initialise wiki file
   cat << EOF > $wiki
=Phantom nightly build results=
Performed `date`

==Changesets tested since last build==
EOF
   i=0;
   for changeset in $changesets; do
       ref=$((i++));
       echo "<tr><td><a href=\"$url/changeset/$changeset\">$changeset</a></td><td>${summary[$ref]}</td></tr>" >> $mailtmp;
       echo "| [[$url/changeset/$changeset|$changeset]] | ${summary[$ref]}" >> $wiki;
   done
   echo >> $wiki;
   echo "</table>" >> $mailtmp;  # start blank mail file
#
#--write content of wiki files and determine if there were errors/failures
#
   cd $codedir/logs;
   warnings='';
   errors='';
   failtext="build";
   prev='';
   for sys in $systems; do
       export SYSTEM=$sys;
       cat test-status-$SYSTEM.wiki >> $wiki
       cat test-status-$SYSTEM.html >> $mailtmp
       cat build-status-$SYSTEM.wiki >> $wiki
       cat build-status-$SYSTEM.html >> $mailtmp
       files=`ls make-*errors*-$SYSTEM.txt make-*errors*-$SYSTEM-debug.txt test-results*-$SYSTEM.txt test-results*-$SYSTEM-debug.txt`
       for x in $files; do
           wikifile="$wikidir/nightly/${x/.txt/.wiki}";
           echo '{{{' > $wikifile;
           cat $x >> $wikifile;
           echo '}}}' >> $wikifile;
       done
       faillog=$codedir/logs/build-failures-$SYSTEM.txt;
       faillogtest=$codedir/logs/test-failures-$SYSTEM.txt;
       if [ -e $faillog ]; then
          fails='';
          for fail in `cat $faillog`; do
             fails+="$fail ";
          done
          if [ "X$fails" != "X" ]; then
             if [ "X$prev" != "X" ]; then
                errors+=" and";
             fi
             errors+=" build failure for SETUP=$fails with SYSTEM=$SYSTEM";
             prev="True";
          fi
       fi
       if [ -e $faillogtest ]; then
          if [ -e $faillog ]; then
             errors+=" and";
             failtext="build and test";
          else
             if [ "X$prev" != "X" ]; then
                errors+=" and";
             fi
             failtext="test";
          fi
          errors+=" testsuite failures";
          prev="True";
       fi
       grep -q 'NEW WARNINGS' $codedir/logs/build-status-$SYSTEM.wiki; gotwarn=$?;
       if [ $gotwarn -eq 0 ]; then
          if [ "X$prev" != "X" ]; then
             warnings+=" and";
          fi
          warnings+=" new compiler warnings with SYSTEM=$SYSTEM";
          prev="True":
       fi
   done
#
#--compose an email with a relevant summary of results
#  also compose a git tag to tag the revision in the repo with the results
#
   msg='';
   text='Please ';
   gotissues=0;
   sendmail=1;
   datetag=`date "+%Y%m%d"`;
   if [ "X$errors" != "X" ]; then
      gotissues=1;
      msg+="$errors";
      text+="<strong>immediately fix</strong> $failtext failures that were your fault";
      if [ "X$warnings" != "X" ]; then
         msg+="$warnings";
         text+=", and avoid introducing new warnings";
      fi
      gittag="${failtext/ and /}fail-$datetag"
   else
      if [ "X$warnings" != "X" ]; then
         gotissues=1;
         msg+="$warnings";
         text+="avoid introducing new warnings where possible";
      fi
      gittag="ok-$datetag"
   fi
   echo $mytag > $tagfile;
   if [ $gotissues -eq 0 ]; then
      preamble='Congratulations! Your changes to Phantom in the last 24 hours passed all tests.';
      text='Please give yourself a pat on the back';
      msg=' Congratulations!';
   else
      preamble="I detected ${msg/build failure/<strong>build failure</strong>} due to changes pushed in the last 24 hours.";
      preamble="${preamble/testsuite failures/<strong>testsuite failures</strong>}";
   fi
   if [ "X$mailto" == "X" ]; then
      echo "BLANK users list in mail, sending to admin only";
      mailto=$admin;
      names="phantom-admin,";
   fi
#
#--compose the email
#
   if [ $sendmail -eq 1 ]; then
      subject="[buildbot/phantom]:$msg";
      echo "sending mail to $mailto";
      echo "subject=$subject";
      echo "warnings=$warnings";
      cat << EOM > $mailfile
To: $mailto
Cc: $admin
From: Phantom buildbot <$admin>
Subject: $subject
Content-Type: text/html
<p>==This is an automated email==</p>

<p>Dear $names</p>

<p>$preamble Details below, or on the <a href="$url/wiki/nightly/Build">wiki page</a>.</p>

<p>$text.</p>

<p>This version of the code has been tagged as $gittag.</p>

<p>With warmest regards from the bowels of the Monash Campus Cluster,</p>
<em>The Phantom buildbot.</em>
<br/>
<br/>
<hr/>
EOM
      cat $mailtmp >> $mailfile;
   fi
}
tag_code_and_push_tags()
{
   cd $codedir
   git tag $gittag
   git push --tags
}
send_email ()
{
  echo "--- sending results email ---";
  if [ $sendmail -eq 1 ]; then
     cat $mailfile | /usr/sbin/sendmail -t;
  fi
}
commit_and_push_wiki ()
{
   echo "--- commit and push to wiki ---";
   # commit and push changes to wiki
   cd $wikidir
   git add nightly/*.wiki
   git status
   git commit -m "[buildbot]: results `date`"
   git push
}
# enter code directory
cd $codedir
get_incoming_changes
pull_changes
run_buildbot
pull_wiki
write_wiki_tag_and_mail_files
tag_code_and_push_tags
send_email
commit_and_push_wiki
cd $dir
echo "--- finished buildbot `date` ---";
