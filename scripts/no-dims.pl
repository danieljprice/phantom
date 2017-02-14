#!/usr/bin/perl
#
# @(#) Perl script to remove dimension() statements in Fortran 90 code
# @(#) Also aligns the intent() :: statements so they line up neatly
#
use Text::Balanced qw<extract_bracketed>;
my @file = '';
while (<STDIN>) {
  my @block = ();
  my $gotblock = 0;
  #
  # match variable declaration block
  #
  while ( m/^(\s*)([real,integer,logical,character].*)\:\:\s*(.*)/ ) {
     my $pad = $1; # spaces
     my $argline = $2; # everything to left of  ::
     my $varline = $3; # everything to right of ::
     $argline =~ s/(\S)$/$1 /g; # ensure there is a space between variables and ::
     my $dims = '';
     my $newargline = '';
     my $newvarline = '';
     my $newline = '';
     my $comment = '';
     # only apply to intent/dimension lines and do not touch line
     # if it contains a line continuation
     if ( ($argline =~ m/intent|dimension/g) and not ($varline =~ m/\&|\=/g) ) {
        if ($argline =~ m/dimension(\(.*\))/) {
           my $tmp = $1;
           $dims = extract_bracketed( $tmp, '()' );
        }

        # remove dimension statement from argline string
        # (this leaves blank argument, see below)
        $argline =~ s/,dimension/, dimension/g; # make sure there is a space
        my $str  = quotemeta("dimension$dims");
        $argline =~ s/$str//g;
        
        # make sure there are spaces between variables
        $argline =~ s/(.*?),([a-z]+)/$1, $2/g;

        # get remaining arguments from splitting argline
        my @args = split(', ',$argline);
        my $count=0;
        foreach $arg (@args) {
           if ( !($arg =~ m/^\s*$/) ) { # skip blank arguments
              $arg =~ s/\s//g;  # remove all spaces
              if ($count==0) {
                 $newargline = "$arg";
              } else {
                 $newargline = "$newargline, $arg";
              }
              $count++;
           }
        }
        
        if ( $varline =~ m/(.*?)(\s*\!.*)/ ) {
           $varline = $1;
           $comment = $2;
        }
        my @vars = split(',',$varline);
        $count=0;
        foreach $var (@vars) {
           $var =~ s/\s*$//g; # strip spaces at end of variable
           if ($count==0) {
              $newvarline = "$var$dims";
           } else {
              if ($var =~ m/\w+/) {
                 # only add space if argument contains word character
                 $newvarline = "$newvarline,$var$dims";
              } else {
                 $newvarline = "$newvarline,$var$dims";
              }
           }
           $count++;
        }
        $newline = "$pad$newargline :: $newvarline$comment\n";
     } else {
        # make sure there are spaces between variables
        $argline =~ s/(.*?),([a-z]+)/$1, $2/g;

        # regurgitate line, with spacing adjusted
        $newline = "$pad$argline\:\: $varline\n";
     }
     push @block,$newline;
     $_ = <>;
     $gotblock = 1;
  }
  if ($gotblock==1) {
     # pad the block to same length
     my @newblock = ();
     my $padlen = 0;
     my $padleni = 0;
     my $padlenc = 0;
     # work out alignment for real/integer/etc. and intent statements
     foreach $line (@block) {
       if ( $line =~ m/(.*)(intent.*\))(.*)(\:\:.*)/ ) {
          my $len = length($1);
          my $leni = length($2);
          if ($len > $padlen) { $padlen = $len; }
          if ($leni > $padleni) { $padleni = $leni; }
          
          # get max length of "intent()," part (with comma)
          if ( $line =~ m/.*(intent.*?\),)/ ) {
             my $lenc = length($1);
             if ($lenc > $padlenc) { $padlenc = $lenc; }
          }
       }
     }
     # align all the real/integer/etc and intent() :: statements
     foreach $line (@block) {
       if ( $line =~ m/(.*)(intent.*\))(.*)(\:\:.*)/ ) {
          my $vari=$1;
          my $in=$2;
          my $ex=$3;
          my $rest=$4;
          # if string contains "intent(), blah..." pad after comma
          if ( $ex =~ m/^\,.*/ ) {
             $ex =~ s/^\,//;  # strip comma
             $newline = sprintf("%-*s%-*s%s%s\n",$padlen,$vari,$padlenc,"$in,",$ex,$rest);             
          } else {
          # otherwise pad as normal
             $newline = sprintf("%-*s%-*s%s%s\n",$padlen,$vari,$padleni,$in,$ex,$rest);
          }
       } else {
          $newline = $line;
       }
       push @newblock, $newline;
     }
     push @file, @newblock;
     $gotblock = 0;
  }
  push @file, $_;
}
print @file;
