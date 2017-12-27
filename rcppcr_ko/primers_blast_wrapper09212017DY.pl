

#!/usr/bin/perl -w

use strict;

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

    my $savefile = $seq_dir . '/blast/sh.primers_blast/primers_blast'.$q.'.sh';
    open SAVE,">>$savefile";
    
    my $command = 'blastn -task blastn -strand plus -db /home/t14905dy/projects/RCP-PCR/KO_clone/Data/db/fasta/const-seq.fna -outfmt "10 qseqid sseqid qstart qend sstart send evalue bitscore btop" -gapopen 2 -gapextend 2 ';
    my $blastsave = $seq_dir . '/blast/out.primers_blast/'.$name.'.blast';
    $command .= '-query '.$input.' -out '.$blastsave;
    print SAVE "$command\n";
    
    $num++;

    close SAVE;
}
