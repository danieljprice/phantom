#!/usr/bin/perl
#
# @(#) Perl script to remove tabs
#
while (<STDIN>) {
  s/\t/        /;
  print "$_";
}
