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
cd "${0%/*}"
tmpdir="/tmp/";
pwd=$PWD;
phantomdir="$pwd/../";
scriptdir="$phantomdir/scripts";
codedir="../";
if [ ! -d $codedir ]; then
   echo "Error running bots: $codedir does not exist";
   exit 100;
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
doindent=1;
gitstaged=0;
input_file='';
bot_names='';
while [[ "$1" == --* ]]; do
  case $1 in
    --apply)
      applychanges=1;
      ;;

    --commit)
      docommit=1;
      applychanges=1;
      ;;

    --no-indent)
      doindent=0;
      ;;

   --staged-files-only)
      gitstaged=1;
      ;;

   --file)
      shift
      input_file=$1;
      break;
      ;;

   --files)
      shift
      input_files="$*";
      break;
      ;;

    --only)
      shift
      bot_names=$1;
      ;;

    *)
      badflag=$1
      ;;
  esac
  shift
done

if [[ "$badflag" != "" ]]; then
   echo "ERROR: Unknown flag $badflag"
   exit
fi

if [[ $gitstaged == 1 && $docommit == 1 ]]; then
   echo "--staged-files-only and --commit cannot both be used because "
   echo "this will commit your git changes with the automated commit message"
   exit
fi

if [[ $doindent == 1 ]]; then
   if ! command -v findent > /dev/null; then
      echo "ERROR: findent not found, please install:                   ";
      echo "       https://www.ratrabbit.nl/ratrabbit/findent/index.html";
      echo "                                                            ";
      echo "       or disable indent-bot using --no-indent              ";
      exit
   fi
fi
cd $codedir;
if [[ $docommit == 1 ]]; then
   git pull;
fi
get_only_files_in_git()
{
   mylist='';
   for file in $1; do
       git ls-files --error-unmatch $file >& /dev/null; giterr=$?;
       if [ $giterr == 0 ]; then
          mylist+="$file ";
       fi
   done
   # skip case where no f90 files present, bash returns literal '*.*90'
   if [[ "$mylist" != '*.*90 ' ]]; then
      echo $mylist
   fi
}
allfiles='';
bots_to_run='tabs gt shout header whitespace authors endif';
if [[ $doindent == 1 ]]; then
   bots_to_run="${bots_to_run} indent";
fi
#bots_to_run='shout';

#
# if --only flag is given, override list of bots_to_run
#
if [[ "$bot_names" != "" ]]; then
   bots_to_run="${bot_names}"
   echo ">> Running only ${bots_to_run} bots via --only flag"
fi
#
# cycle over all the possible bots and run them...
#
modified=0
for edittype in $bots_to_run; do
    filelist='';
    case $edittype in
    'authors' )
       dirlist=".";
       goback='-';
       filenamepattern="$authorsfile";;
    * )
       dirlist="src/*";
       goback="../../";
       filenamepattern="*.*90";;
    esac
    for dir in $dirlist; do
        if [ -d $dir ]; then
           cd $dir;
           if [[ $gitstaged == 1 && "$filenamepattern" == "*.*90" ]]; then
             files=`git diff --name-only --cached --relative -- "./*.*90" | tr '\n' ' '`;
           else
             files=$filenamepattern
           fi
           if [[ "$input_file" != "" && "$edittype" != "authors" ]]; then
           # Only process the input file
             if [[ "$dir" == "$(dirname $input_file)" ]]; then
               myfiles=$(basename $input_file)
             else
               myfiles=""
             fi
           elif [[ "$input_files" != "" ]]; then
             # if using pre-commit to pass in a list of files, then there is no need to check
             # if these files are in git
             myfiles="$files"
           else
             myfiles=`get_only_files_in_git "$files"`
           fi
           for file in $myfiles; do
               if [[ "$input_files" != "" && "$input_files" != *"$dir/$file"* && "$edittype" != "authors" ]]; then
                 continue
               fi
               if [[ "$file" == "libphantom-evolve.F90" ]]; then # skip as causes indent-bot trouble
                  continue
               fi
               out="$tmpdir/$file"
#               echo "FILE=$file OUT=$out";
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
                     -e 's/REAL(/real(/g' \
                     -e 's/DBLE(/dble(/g' \
                     -e 's/ STOP/ stop/g' \
                     -e 's/ATAN/atan/g' \
                     -e 's/ACOS(/acos(/g' \
                     -e 's/ASIN(/asin(/g' \
                     -e 's/COS(/cos(/g' \
                     -e 's/SIN(/sin(/g' \
                     -e 's/EXP(/exp(/g' \
                     -e 's/LOG(/log(/g' \
                     -e 's/READ(/read(/g' \
                     -e 's/OPEN(/open(/g' \
                     -e 's/OPEN (/open(/g' \
                     -e 's/CLOSE(/close(/g' \
                     -e 's/CLOSE (/close(/g' \
                     -e 's/INDEX(/index(/g' \
                     -e 's/ANY(/any(/g' \
                     -e 's/, STATUS=/,status=/g' \
                     -e 's/,STATUS=/,status=/g' \
                     -e 's/, FORM=/,form=/g' \
                     -e 's/,FORM=/,form=/g' \
                     -e 's/TRIM(/trim(/g' \
                     -e 's/IF (/if (/g' \
                     -e 's/) THEN/) then/g' \
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
                     -e 's/KIND=/kind=/g' \
                     -e 's/ INTEGER/ integer/g' \
                     -e 's/ PARAMETER/ parameter/g' \
                     -e 's/INTENT(IN)/intent(in)/g' \
                     -e 's/INTENT(OUT)/intent(out)/g' \
                     -e 's/INTENT(INOUT)/intent(inout)/g' \
                     -e 's/intent(IN)/intent(in)/g' \
                     -e 's/intent(OUT)/intent(out)/g' \
                     -e 's/intent(INOUT)/intent(inout)/g' \
                     -e 's/INT(/int(/g' \
                     -e 's/MIN(/min(/g' \
                     -e 's/ABS(/abs(/g' \
                     -e 's/NOT(/not(/g' \
                     -e 's/IAND(/iand(/g' \
                     -e 's/IEOR(/ieor(/g' \
                     -e 's/MODULO(/modulo(/g' \
                     -e 's/HUGE(/huge(/g' \
                     -e 's/TINY(/tiny(/g' \
                     -e 's/SELECTED_REAL_KIND(/selected_real_kind(/g' \
                     -e 's/SELECTED_INT_KIND(/selected_int_kind(/g' \
                     -e 's/SELECTED_REAL_kind(/selected_real_kind(/g' \
                     -e 's/SELECTED_INT_kind(/selected_int_kind(/g' \
                     -e 's/KIND(/kind(/g' \
                     -e 's/ REAL/ real/g' $file > $out;;
               'endif' )
                 sed -e 's/end if/endif/g' \
                     -e 's/end do/enddo/g' \
                     -e 's/else if/elseif/g' \
                     -e 's/open (/open(/g' \
                     -e 's/, file = /,file=/g' \
                     -e 's/, file=/,file=/g' \
                     -e 's/, status = /,status=/g' \
                     -e 's/, status=/,status=/g' \
                     -e 's/, iostat = /,iostat=/g' \
                     -e 's/, iostat=/,iostat=/g' \
                     -e 's/, access = /,access=/g' \
                     -e 's/, access=/,access=/g' \
                     -e 's/, form = /,form=/g' \
                     -e 's/, form=/,form=/g' \
                     -e 's/, action = /,action=/g' \
                     -e 's/, action=/,action=/g' \
                     -e 's/, iomsg = /,iomsg=/g' \
                     -e 's/, iomsg=/,iomsg=/g' \
                     -e 's/(unit =/(unit=/g' \
                     -e 's/if(/if (/g' \
                     -e 's/)then/) then/g' $file > $out;;
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
                  if command -v findent > /dev/null; then
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
      msg='[format-bot] end if -> endif; end do -> enddo; if( -> if (';;
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
       modified=$((modified + 1))
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
exit $modified
