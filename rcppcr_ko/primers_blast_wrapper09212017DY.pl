

#!/usr/bin/perl -w

use strict;
my $db =shift;
my $seq_dir = shift; ##### FULL PATH REQIRED!
$seq_dir =~ s/\/$//g;

my @inputs = @ARGV; ##### FULL PATH REQUIRED!

my $num=0;
my $q  =0;
for my $input (@inputs){
    my $name = 'na';
    if($input =~ /\/([^\/]+)\.fna/){
	$name = $1;
    }

    $q++ if(!($num % 30));

    my $savefile = $seq_dir . '/blast/sh.blast/primers_blast'.$q.'.sh';
    open SAVE,">>$savefile";

    my $command = 'blastn -task blastn -strand plus -db' .$db.  '-outfmt "10 qseqid sseqid qstart qend sstart send evalue bitscore btop" -gapopen 2 -gapextend 2 ';
    my $blastsave = $seq_dir . '/blast/out.blast/'.$name.'.blast';
    $command .= '-query '.$input.' -out '.$blastsave;
    print SAVE "$command\n";

    $num++;

    close SAVE;
}
