#!/usr/bin/perl
# prints a random word from the dictionary
my @wordlist = `cat /usr/share/dict/words`;
my $randomno = int(rand($#wordlist));
print @wordlist[$randomno];
