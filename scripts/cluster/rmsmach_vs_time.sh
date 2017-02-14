#!/bin/bash
#
# A bash script to read time and rms mach number from logfiles,
# and output to a file.
#
# Example use:
#		./rmsmach_vs_time.sh cluster01.log cluster02.log
#


#
#--List of the logfiles that need to be combined
#
#listARRAY=(cluster01.log cluster02.log cluster05.log cluster06.log)
listARRAY=("$@")
list=${listARRAY[*]}
l=${#listARRAY[*]} 							# length of the array
l=$(($l-1))

#
#--If no logfiles have been provided, complain.
#
if [ "$#" -lt "1" ]; then
	echo "Please provide logfiles as input:"
	echo "e.g. rmsmach_vs_time.sh cluster01.log cluster02.log"
	exit 1;
fi

#
#--Make temporary folder for temp files
#
#folder="tempfolderforrmsmach"
#mkdir $folder

#
#--Names for temporary files
#
#time="$folder/time.temp"
#rms="$folder/rmsmach.temp"
time="123time.temp123"
rms="123rmsmach.temp123"
out="rmsmach_vs_time.data"

#
#--Store the length of every log file (except the last one)
#
prev=0
for i in $(seq 0 $l); do
	len[$i]=$(($(grep -c 'TIME =' ${listARRAY[$i]})+$prev))
	prev=len[$i]
done

#
#--Find the times and rms mach numbers in all the log files
#
grep -h 'TIME =' $list | cut -d : -f1 | cut -c22-32 > $time
grep -h 'RMS Mach #=' $list | tr -d 'RMS Mach #=' > $rms

#
#--Delete the duplicate values of rms Mach no. that come from the top of every continued log file
#
for i in ${len[*]}; do
	sed -i ""$i"d" $rms
	#sed ""$i"d" $rms > ".some-temporary-file1234" && mv ".some-temporary-file1234" $rms
done

#
#--Merge time and rms Mach no. files together as columns (i.e. col1=time col2=rmsmach)
#
paste $time $rms > $out

#
#--Remove temporary files and folder
#
#rm -rf $folder
rm $time $rms

#
#--Open in splash
#
#splashprefix="rmsmach_vs_time_splash"
#printf "time\nRMS mach no." > "$splashprefix.columns"
#asplash $out -p $splashprefix
