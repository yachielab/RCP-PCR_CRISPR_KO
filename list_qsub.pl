#!/usr/bin/perl -w 

use strict;

my @sh_files = @ARGV;

for my $sh (@sh_files){
    print "qsub -o /dev/null -e /dev/null $sh\n";
}
