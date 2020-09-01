#!/usr/bin/env bash
# script to auto-generate .rst files for Fortran source
# with links to the autodoc output from Sphinx
subdir='main'
mydir=$PWD;
cd ../../../src/$subdir;
#files=`ls *.*90 | grep -v 'setup_' | grep -v 'libsetup'`
files=`ls *.f90 | grep -v 'inject_' | grep -v 'kernel_' | grep -v 'metric_' | grep -v 'bondi'`
for x in $files; do
    y=${x/.f90/}
    z=${y/.F90/}
    outfile=$mydir/"${z}.rst"
    echo > $outfile
    echo "$z" >> $outfile
    echo "=========================" >> $outfile
    echo -e "\n.. f:autosrcfile:: ../../../src/$subdir/$z.f90\n" >> $outfile
done;
cd $mydir
index="index.rst"
indext="${index}.tmp"
rm -i $index
echo "$subdir" > $indext
echo "===========" >> $indext
echo -e ".. toctree:: \n   :maxdepth: 2\n\n" >> $indext
for x in *.rst; do
    echo "   ${x/.rst/}" >> $indext
done
mv $indext $mydir/$index
