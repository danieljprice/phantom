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
dir=$PWD
datetag=`date "+%Y%m%d"`;
#---------------------
#  script settings
#---------------------
#webdir=$dir/web
webdir=$dir/phantomsph.bitbucket.org
codedir=$dir/phantom
benchdir=$dir/phantom-benchmarks
url="https://phantomsph.bitbucket.io/"
urlgitrepo="https://bitbucket.org/danielprice/phantom";
#url="http://users.monash.edu.au/~dprice/phantom";
#webserver="users.monash.edu.au:WWW/phantom/";
urllogs="https://users.monash.edu.au/~dprice/phantom/nightly/logs/"
sender="daniel.price@infra.monash.edu";
admin="daniel.price@monash.edu";
systems="msg gfortran";
mailfile="$dir/mail.tmp";
htmlfile="$dir/$datetag.html";
incoming="$dir/incoming.txt";
tagfile="$dir/gittag.tmp";
#---------------------
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
   # get list of incoming changesets
   changesets=`git log ..@{u} --format="%h"`;
#   changesets=`grep commit $incoming | cut -d':' -f 2`;
   git log ..@{u} --format="%s" > $dir/summary.txt;
}
extract_names_of_users ()
{
   # extract list of people to send mail to from users.list
   mailto='';
   names='';
   while read line; do
      mailto+="$line,";
      names+="`echo $line | cut -d' ' -f 1`, ";
   done < $dir/users.list
   echo "mailto: $mailto";
}
extract_changeset_list ()
{
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
   cd $benchdir; git pull >& gitpull.log
   cd $codedir;
}
run_buildbot ()
{
   # run the buildbot and testbot scripts
   echo "--- running buildbot ---";
   cd $codedir/scripts
   for sys in $systems; do
      export SYSTEM=$sys;
      echo "SYSTEM=$SYSTEM";
      export PHANTOM_DIR=$codedir; # so setup tests can find data files
      ./testbot.sh "$urllogs/nightly/logs/";
      ./buildbot.sh 17000000 "$urllogs/nightly/logs/";
   done
}
run_benchmarks ()
{
   # run performance suite
   echo "--- running benchmarks ---";
   cd $benchdir;
   export PHANTOM_DIR=$codedir;
   for sys in $systems; do
       export SYSTEM=$sys;
       echo "SYSTEM=$sys";
       ./run-benchmarks.sh;
   done
   ./plot-benchmarks.sh;
}
pull_wiki ()
{
   # pull existing web files from server
   echo "--- pulling website changes ---";
   cd $webdir
   git pull
}
write_htmlfile_gittag_and_mailfile ()
{
   echo "--- writing html and mail files ---"
   # write header for nightly build results
   cd $dir;
   # initialise html file
   echo "<h2>Changesets tested since last build</h2>" > $htmlfile
   echo "<table>" >> $htmlfile
   # initialise html file
   i=0;
   for changeset in $changesets; do
       ref=$((i++));
       echo "<tr><td><a href=\"$urlgitrepo/changeset/$changeset\">$changeset</a></td><td>${summary[$ref]}</td></tr>" >> $htmlfile;
   done
   echo "</table>" >> $htmlfile;  # start blank mail file
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
       cat test-status-$SYSTEM.html >> $htmlfile
       cat build-status-$SYSTEM.html >> $htmlfile
       cat $benchdir/opt-status-$SYSTEM.html >> $htmlfile
       files=`ls make-*errors*-$SYSTEM.txt make-*errors*-$SYSTEM-debug.txt test-results*-$SYSTEM.txt test-results*-$SYSTEM-debug.txt`
       weblogdir="$webdir/nightly/logs"
       if [ -d "$weblogdir" ]; then
          for x in $files; do
              webfile="$weblogdir/$x";
              cp $x $webfile;
          done
       else
          echo "ERROR: cannot copy logs to $weblogdir: directory does not exist";
       fi
       faillog=$codedir/logs/build-failures-$SYSTEM.txt;
       faillogtest=$codedir/logs/test-failures-$SYSTEM.txt;
       faillogsetup=$codedir/logs/setup-failures-$SYSTEM.txt;
       if [ -e $faillog ]; then
          fails='';
          for fail in `cat $faillog | sort -u`; do
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
       if [ -e $faillogsetup ]; then
          if [ -e $faillog ]; then
             errors+=" and";
             failtext+=" and setup";
          else
             if [ "X$prev" != "X" ]; then
                errors+=" and";
             fi
             failtext="setup";
          fi
          errors+=" setup failures";
          prev="True";
       fi
       grep -q 'NEW WARNINGS' $codedir/logs/build-status-$SYSTEM.html; gotwarn=$?;
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
   if [ "X$errors" != "X" ]; then
      gotissues=1;
      msg+="$errors";
      text+="<strong>immediately fix</strong> $failtext failures that were your fault";
      if [ "X$warnings" != "X" ]; then
         msg+="$warnings";
         text+=", and avoid introducing new warnings";
      fi
      failtext=${failtext/ and /};
      gittag="${failtext/ and /}fail-$datetag";
      gittag="${gittag/ and /}";
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
      preamble="Detected ${msg}";
      preamblehtml="I detected ${msg/build failure/<strong>build failure</strong>} due to changes pushed in the last 24 hours.";
      preamblehtml="${preamble/testsuite failures/<strong>testsuite failures</strong>}";
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

<p>$preamblehtml Details below, or on the <a href="$url/nightly/">web page</a>.</p>

<p>$text.</p>

<p>This version of the code has been tagged as $gittag.</p>

<p>With warmest regards from the bowels of the Monash Campus Cluster,</p>
<em>The Phantom buildbot.</em>
<br/>
<br/>
<hr/>
EOM
      cat $htmlfile >> $mailfile;
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
post_to_slack ()
{
  message=$1;
  webhookurl="https://hooks.slack.com/services/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"
  channel="#commits"
  username="buildbot"
  json="{\"channel\": \"$channel\", \"username\": \"$username\", \"text\": \"$message\", \"icon_emoji\": \":ghost:\"}"

  curl -X POST --data-urlencode "payload=$json" $webhookurl
}
commit_and_push_to_website ()
{
   echo "--- commit and push to web server / git repo ---";
   # commit and push changes to web server
   cp $htmlfile $webdir/nightly/build;
   cp $benchdir/performance.html ${webdir}/nightly/opt/index.html;
   cp $benchdir/*.js ${webdir}/nightly/opt/;
   cd $webdir/nightly/build;
   if [ -e $htmlfile ]; then
      cp $htmlfile index.html;
   fi
   cd $webdir;
   #rsync -avz nightly/ $webserver/nightly/;
   git pull
   git add nightly/*.html
   git add nightly/opt/*
   git add nightly/stats/*
   git status
   git commit -m "[nightly]: results `date`"
   git push
}
# enter code directory
cd $codedir
get_incoming_changes
extract_names_of_users
extract_changeset_list
pull_changes
run_buildbot
run_benchmarks
#pull_wiki
write_htmlfile_gittag_and_mailfile
#longmessage="${names/,/}: $preamble $text"
#echo "$longmessage"
message="status: <$url/nightly/build/$datetag.html|$gittag>"
post_to_slack "$message"
tag_code_and_push_tags
send_email
commit_and_push_to_website
cd $dir
echo "--- finished buildbot `date` ---";
