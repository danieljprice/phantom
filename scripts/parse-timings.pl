#!/usr/bin/perl
#
# script to parse .html files from nightly build results
# in order to extract timing data
#
# Daniel Price, Jan 2018
#
#
# Usage
#
if ($#ARGV <= 0) {
   die "Usage: $0 20180101.html\n";
}
foreach my $file (@ARGV) {
   open my $fh, '<', $file or die "can't open $file\n";
   my ($year,$mon,$day) = ($file =~ m/(\d\d\d\d)(\d\d)(\d\d)\.html/);
   my %t;
   my @keys=qw( test-msg testkd-msg test-gfortran testkd-gfortran );
   while (<$fh>) {
      #print "$_\n";
      my $line="$_";
      if ($line =~ m/completed in/ && $line !~ m/debug/) {
	 #print "$line";
	 my ($name,$system) = ($line =~ m/test-results-(.*)-(.*)\.txt/);
	 my $min,$sec;
	 if ($line =~ m/completed in\s*(\d*\.?\d*)\s+min,\s+(\d*\.?\d*)\s+s/) {
            $min=$1;
	    $sec=$2;
	 } else {
            $min=0;
	    ($sec)=($line =~ m/completed in\s*(\d*\.?\d*)\s*s/)
	 }
	 my $time=$min*60 + $sec;
	 my $tag="$name-$system";
	 $t{$tag} = $time;
	 #print "$name($system) $time";
      }
   }
   $mon = $mon - 1;
   if ($t{'test-msg'} > 0 && $t{'test-msg'} < 350) {
      print "        [new Date($year,$mon,$day)";
      for my $key (@keys) {
          if ($t{$key}) {
	     print ",$t{$key}";
	  } else {
	     print ",500.0";
	  }
      }
      print "],\n";
   }
   #for my $key (keys % t) {
   #    print "$key is $t{$key}\n";
   #}
}

