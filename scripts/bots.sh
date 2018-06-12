#!/bin/bash
#
# Script to automatically perform various maintenance
# tasks on Phantom source files and commit the 
# changes to the code
#
# Current bots implemented are:
#
# [tab-bot]: automatically removes tabs from source files
# [format-bot]: replaces F77-style .gt.,.ge. etc. with >, >=
# [header-bot]: updates the module headers in source files
#               (relies on external perl script)
# [author-bot]: updates AUTHORS file
# [indent-bot]: adjusts indentation (relies on external findent program)
#
# Written by Daniel Price 2014-
#
tmpdir="/tmp/";
pwd=$PWD;
phantomdir="$pwd/../";
if [ ! -s $phantomdir/scripts/$0 ]; then
   echo "Error: This script needs to be run from the phantom/scripts directory";
   exit;
fi
scriptdir="$phantomdir/scripts";
codedir="../";
if [ ! -d $codedir ]; then
   echo "Error running bots: $codedir does not exist";
   exit;
fi
headerfile="$scriptdir/HEADER-module";
programfile="$scriptdir/HEADER-program";
authorsfile="AUTHORS";
writegitfile="writegitinfo.f90"
if [ ! -s $headerfile ]; then
   echo "WARNING: cannot find $headerfile: header-bot cannot run";
fi
if [ ! -s $programfile ]; then
   echo "WARNING: cannot find $programfile: header-bot cannot run";
fi
if [ ! -s $codedir/$authorsfile ]; then
   echo "WARNING: cannot find $authorsfile: author-bot cannot run";
fi
docommit=0;
applychanges=0;
if [[ "$1" == "--apply" ]]; then
   applychanges=1;
fi
if [[ "$1" == "--commit" ]]; then
   docommit=1;
   applychanges=1;
fi
cd $codedir;
if [[ $docommit == 1 ]]; then
   git pull;
fi
allfiles='';
bots_to_run='tabs gt shout header whitespace authors endif indent';
#bots_to_run='shout';
for edittype in $bots_to_run; do
    filelist='';
    case $edittype in
    'authors' )
       dirlist=".";
       goback='-';
       listoffiles="$authorsfile";;
    * )
       dirlist="src/*";
       goback="../../";
       listoffiles="*.*90";;
    esac
    for dir in $dirlist; do
        if [ -d $dir ]; then
           cd $dir;
           #echo $srcdir;
           for file in $listoffiles; do
               out="$tmpdir/$file"
               case $edittype in
               'tabs' )
                 sed 's/	/        /g' $file > $out;;
               'whitespace' )
                 sed 's/ *$//' $file > $out;;
               'gt' )
                 sed -e 's/\.gt\./ \> /g' \
                     -e 's/\.GT\./ \> /g' \
                     -e 's/\.lt\./ \< /g' \
                     -e 's/\.LT\./ \< /g' \
                     -e 's/\.le\./ \<= /g' \
                     -e 's/\.LE\./ \<= /g' \
                     -e 's/\.ge\./ \>= /g' \
                     -e 's/\.GE\./ \>= /g' \
                     -e 's/\.eq\./==/g' \
                     -e 's/\.EQ\./==/g' \
                     -e 's/\.NE\./ \/= /g' \
                     -e 's/\.ne\./ \/= /g' $file > $out;;
               'shout' )
                 sed -e 's/SQRT(/sqrt(/g' \
                     -e 's/NINT(/nint(/g' \
                     -e 's/STOP/stop/g' \
                     -e 's/ATAN/atan/g' \
                     -e 's/ACOS(/acos(/g' \
                     -e 's/ASIN(/asin(/g' \
                     -e 's/COS(/cos(/g' \
                     -e 's/SIN(/sin(/g' \
                     -e 's/EXP(/exp(/g' \
                     -e 's/LOG(/log(/g' \
                     -e 's/READ(/read(/g' \
                     -e 's/ENDDO/enddo/g' \
                     -e 's/END DO/end do/g' \
                     -e 's/END PARALLEL/end parallel/g' \
                     -e 's/END SUBROUTINE/end subroutine/g' \
                     -e 's/END MODULE/end module/g' \
                     -e 's/CONTAINS/contains/g' \
                     -e 's/$OMP PARALLEL/$omp parallel/g' \
                     -e 's/$OMP DO SCHEDULE/$omp do schedule/g' \
                     -e 's/$OMP/$omp/g' \
                     -e 's/WRITE(/write(/g' \
                     -e 's/^USE /use /g' \
                     -e 's/^ USE / use /g' \
                     -e 's/^  USE / use /g' \
                     -e 's/^CALL /call /g' \
                     -e 's/CHARACTER(/character(/g' \
                     -e 's/SUBROUTINE /subroutine /g' \
                     -e 's/^MODULE /module /g' \
                     -e 's/RETURN/return/g' \
                     -e 's/DOT_PRODUCT/dot_product/g' \
                     -e 's/ DIMENSION(/ dimension(/g' \
                     -e 's/CONTINUE/continue/g' \
                     -e 's/ INTEGER/ integer/g' \
                     -e 's/ REAL/ real/g' $file > $out;;
               'endif' )
                 sed -e 's/end if/endif/g' \
                     -e 's/end do/enddo/g' $file > $out;;
               'header' )
                 $scriptdir/header.pl --headerfile=$headerfile --programfile=$programfile --replace $file > $out;;
               'authors' )
                  if [ "$file"=="$authorsfile" ]; then
                     echo '#-------------------------------------------------------#' > $out;
                     echo '# Contributors to Phantom, updated automatically using  #' >> $out;
                     echo '#                                                       #' >> $out;
                     echo '#  git shortlog -s -n -e | cut -f 2                     #' >> $out;
                     echo '#                                                       #' >> $out;
                     echo '# Edit .mailmap if your name or email are wrong         #' >> $out;
                     echo '#-------------------------------------------------------#' >> $out;
                     git shortlog -s -n -e HEAD | grep -v global | cut -f 2 >> $out;
                  else
                     cat $file > $out;
                  fi;;
               'indent' )
                  if type -p findent; then
                     findent -r1 -m1 -c3 -Rr -C- -k- -j1 < $file > $out;
                  fi;;
               esac
               if [ -s $out ]; then
                  if [[ `diff -q $out $file` ]]; then
                     filelist+=" $dir/$file";
                     if [[ $applychanges == 1 ]]; then
                        cp $out .;
                     else
                        echo "--- $file ---";
                        diff $out $file;
                     fi
                  fi
                  rm $out;
               else
                  echo "ERROR in bot $edittype - no output file $out";
               fi
           done
           cd $goback;
        fi
    done
    case $edittype in
    'tabs' )
      msg='[tab-bot] tabs removed';;
    'whitespace' )
      msg='[space-bot] whitespace at end of lines removed';;
    'gt' )
      msg='[format-bot] obsolete .gt. .lt. .ge. .le. .eq. .ne. replaced';;
    'shout' )
      msg='[format-bot] F77-style SHOUTING removed';;
    'endif' )
      msg='[format-bot] end if -> endif; end do -> enddo';;
    'header' )
      msg='[header-bot] updated file headers';;
    'authors' )
      msg='[author-bot] updated AUTHORS file';;
    'indent' )
      msg='[indent-bot] standardised indentation';;
    esac
    if [[ "X$filelist" != "X" ]]; then 
       if [[ $docommit == 0 ]]; then
          echo "$msg";
       fi
       echo "Modified files = $filelist";
    fi
    if [[ $docommit == 1 ]]; then
       git commit -m "$msg" $filelist;
    fi
    allfiles+=$filelist;
done
if [[ $docommit == 1 ]]; then
   git push;
else
   if [[ "X$allfiles" != "X" ]]; then 
      if [[ $applychanges == 0 ]]; then
         echo; echo "To apply changes use:";
         echo; echo "$0 --apply";
      fi
      echo; echo "To apply and commit changes use:";
      echo; echo "$0 --commit"; echo;
   else
      echo "No changes";
   fi
fi
