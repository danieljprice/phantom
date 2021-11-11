#!/usr/bin/perl
# prints a random word from the dictionary
my $dict='/usr/share/dict/words';
my @wordlist;
if (-e $dict) {
   @wordlist = `cat $dict`;
} else {
   @wordlist = <DATA>;
}
my $randomno = int(rand($#wordlist));
print @wordlist[$randomno];

__DATA__
echidna
wombat
platypus
emu
kangaroo
possum
gecko
koala
dingo
numbat
quokka
quoll
croc
shark
fruitbat
flyingfox
snake
redback
bilby
bandicoot
bettong
dugong
thylacine
potoroo
dunnart
cockatoo
kookaburra
rosella
pademelon
sugarglider
