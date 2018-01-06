#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $primers_fasta  = shift;
my $bar2num_file   = shift;
my $targets        = shift;


my $seq_dir = shift; ##### FULL PATH REQIRED!
$seq_dir =~ s/\/$//g;

my $R1_id = shift;
my $R2_id = shift;

my $seq_DB  = &get_seq($primers_fasta);
my $bar2num = &get_bar2num($bar2num_file);

my $primers_R1_blast = $seq_dir.'/blast/out.primers_blast/'. $R1_id.'.blast';
my $primers_R2_blast = $seq_dir.'/blast/out.primers_blast/'. $R2_id.'.blast';
my $fasta_R1         = $seq_dir.'/fragmented_fasta/'.  $R1_id.'.fna';
my $fasta_R2         = $seq_dir.'/fragmented_fasta/'.  $R2_id.'.fna';

my $seq_reads;
$seq_reads->{R1} = &get_seq($fasta_R1);
$seq_reads->{R2} = &get_seq($fasta_R2);

#print Dumper $seq_reads;


my @files = ($primers_R1_blast,$primers_R2_blast);

my $data;
for my $file (@files){

    my $read_direction = 'R1';
    if($file =~ /\_R2\_/){
	$read_direction = 'R2';
    }

    open FILE,$file;
    while(<FILE>){
	chomp;
	my @array = split /\,/,$_;
	if($#array>5){
	    my ($read,$template,$str_on_read,$end_on_read,$str,$end,$evalue,$btop) =
		@array[0,1,2,3,4,5,6,8];

	    if($str_on_read<$end_on_read && $str == 1 && $end == $seq_DB->{$template}->{length}){ ########## THINK ABOUT THIS
		my $element;
		$element->{str}    = $str_on_read;
		$element->{end}    = $end_on_read;
		$element->{evalue} = $evalue;
		$element->{btop}   = $btop;
		$element->{Target_s}=$str;
		$element->{Target_e}=$end;
		$data->{$read}->{$read_direction}->{$template} = $element;
	    }
	}
    }
    close FILE;
}

#print Dumper $data;


my $data2;
while(my($read,$value) = each %$data){
    my $R1     = $value->{R1};
    my $R2     = $value->{R2};
    my $seq;
    #my $KO_region = 0;
    $seq->{R1} = $seq_reads->{R1}->{$read}->{seq};
    $seq->{R2} = $seq_reads->{R2}->{$read}->{seq};
    #Dumper $R1;



    # PS1.0 and PS2.0 at the reasonable positions? (end at <40bp)
    # Row and Column priming sites from the reasonable positions? (start at <50bp)
    # Which category? DB-BC / DB-lox / AD-BC / AD-lox?
    # Specific category assignment?
    # Row, Column and Plate tags?

    my ($category,$P1_seq,$P2_seq,$R_seq,$C_seq,$locusR1,$locusR2,$orientation) = &assign_category($value,$R1,$R2,$seq,$targets);
    #print Dumper $KO_region;
    #print $category, $P1_seq, $P2_seq, $R_seq, $C_seq;

    my $P1_num = ($P1_seq && length($P1_seq) < 12)? &barcode_matching($bar2num,$P1_seq):0;
    my $P2_num = ($P2_seq && length($P2_seq) < 12)? &barcode_matching($bar2num,$P2_seq):0;
    my $R_num  = ($R_seq  && length($R_seq)  < 12)? &barcode_matching($bar2num,$R_seq) :0;
    my $C_num  = ($C_seq  && length($C_seq)  < 12)? &barcode_matching($bar2num,$C_seq) :0;


    if($category && $P1_num*$P2_num*$R_num*$C_num){

	$P1_num = sprintf "%.1f", $P1_num/10;
	$P1_num =~ s/\.//g;
	$P2_num = sprintf "%.1f", $P2_num/10;
	$P2_num =~ s/\.//g;
	$R_num  = sprintf "%.1f", $R_num/10;
	$R_num  =~ s/\.//g;
	$C_num  = sprintf "%.1f", $C_num/10;
	$C_num  =~ s/\.//g;

       	my ($P_TAG,$R_TAG,$C_TAG) = ("P$P1_num\-P$P2_num","R$R_num","C$C_num");
	#print "$category\t$P_TAG\t$R_TAG\t$C_TAG\n";

	if ( $orientation > 0){
	$data2->{$category}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$read}->{"Plus"}=$locusR1;
	#$data2->{$category}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$read}->{"Minus"}=$locusR2;
	}

	if ( $orientation < 0){
	    #$data2->{$category}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$read}->{"Minus"}=$locusR1;
	    $data2->{$category}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$read}->{"Plus"}=$locusR2;
	}
    }}


print Dumper $data2;


sub btop_matching(){
    my $mut1 = shift;
    my $mut2  = shift;
    my $len      = shift;
    my $s1   = shift;
    my $s2    = shift;
    my @mut_arr;
    #print $main_mut;
    #print $sub_mut;
    my $b1 = &get_btop($mut1);
    my $b2 = &get_btop($mut2);
    print Dumper $b1;
    print Dumper $b2;
    #print $len;
    @mut_arr = &mut_match($len,$b1,$b2)
}


sub mut_match(){
    my $t_size = shift;
    my $m1     = shift;
    my $m2     = shift;
    #print $t_size;
    my @del = ( 0 ) x$t_size;
    print @del;




}




sub get_btop(){
    my $btop = shift;
    my $btop_num    = $btop;
    my $btop_string = $btop;

    $btop_num    =~ s/[^\d]+/\|/g;
    $btop_string =~ s/\d+/\|/g;


    my @array1 = split /\|/,$btop_num;
    my @array2 = split /\|/,$btop_string;


    my @array;

    my $n1 = 0;
    my $n2 = 0;

    while(1){
        last if(!$array1[$n1] && !$array2[$n2]);

        if($array1[$n1]){
            push @array,$array1[$n1];
        }
        if($array2[$n2]){
            push @array,$array2[$n2];
        }else{
            $n2++;
            if($array2[$n2]){
                push @array,$array2[$n2];
            }
        }
$n1++;
$n2++;
    }

    my $btop_data;

    my $s_i = 0;
    for my $element (@array){
        if($element =~ /^\d+$/){
            $s_i += $element;
}else{
            for(my $j=0;$j<length($element);$j+=2){
                my $unit = substr $element,$j,2;
                my $q_stat = substr $unit,0,1;
my $s_stat = substr $unit,1,1;
my $stat   = 'N_A';
                if($q_stat eq '-'){
                    $s_i++;
                    $btop_data->{del}->{$s_i} = $s_stat;
                }elsif($s_stat eq '-'){
                    #my $s_i_next = $s_i+1;
                    $btop_data->{ins}->{$s_i} .= "$q_stat";
                }else{
                    $s_i++;
                    $btop_data->{snp}->{$s_i} = "$s_stat$q_stat";
                }
            }
        }
    }
    $btop_data->{ssize} = $s_i;
    return $btop_data;
}


sub barcode_matching(){
    my $bar2num = shift;
    my $seq     = shift;

    my @bar_count;
    for my $key (keys %$bar2num){

        my $count;
        for(my $k=0; $k<length($key)-1; $k++){
            my $k_nuc = substr $key,$k,2;
            for(my $i=0; $i<length($seq)-1; $i++){
                my $i_nuc = substr $seq,$i,2;

                my $rel = $k-$i;
                if($k_nuc eq $i_nuc){
                    $count->{$rel} ||= 0;
                    $count->{$rel}  += 1;
                }else{
                    $count->{$rel} ||= 0;
                    $count->{$rel}  += 0;
                }
            }
        }

        my @array = sort{$b <=> $a} values %$count;
        my $max_count = shift @array;
        my $element;
        $element->{bar} = $key;
        $element->{max} = $max_count;
        push @bar_count, $element;
    }
    my $val;

    @bar_count = sort{$b->{max} <=> $a->{max}} @bar_count;
    if($bar_count[0]->{max} == $bar_count[1]->{max}){
        $val = 0;
    }elsif($bar_count[0]->{max}>5){ ########## THINK ABOUT THIS THRESHOLD
        my $called_bar = $bar_count[0]->{bar};
        $val = $bar2num->{$called_bar};
    }else{
        $val = 0;
    }

    return $val;
}














############################################################
# functions
###########################################################


sub assign_category(){

    my $value = shift;
    my $R1    = shift;
    my $R2    = shift;
    my $seq   = shift;
    my $targets = shift;


    #print Dumper $R1;
    #print Dumper $R2;



    # PS1.0 and PS2.0 at the reasonable positions? (end at <40bp)
    my $goto1 = 0;


    my $P_SEQ1 = 0;
    my $P_SEQ2 = 0;

    if($R1->{'PS1.0-primer'} && $R2->{'PS2.0-primer'}){
	if($R1->{'PS1.0-primer'}->{end} < 40 && $R2->{'PS2.0-primer'}->{end} < 40){
	    $P_SEQ1 = substr $seq->{R1},$R1->{'PS1.0-primer'}->{str}-10,9;
	    $P_SEQ2 = substr $seq->{R2},$R2->{'PS2.0-primer'}->{str}-10,9;
	    $goto1 = 1;
	}
    }

    # Row and Column priming sites from the reasonable positions? (start at <50bp)
    # Which Target loci?
    my @Target ;
    my $R_seq;
    my $C_seq;
    my $locusR1;
    my $locusR2;
    my $LEN;
    my $orientation;

    if($goto1){
	open FILE,$targets;
	while(<FILE>){
	    chomp;
	    my @values = split(/,/, $_);
	    my $_ = $values[0];
	    my $target_len = length($values[1]);
	    #print "$_"."$target_len\n";


	    if(exists($R1->{$_."_Frd"}) && exists($R2->{$_."_Rvs"}) && exists($R1->{"DBU1-primer"}) && exists($R2->{"DBD2-primer"})){
		if($R1->{'PS1.0-primer'}->{end} < $R1->{$_."_Frd"}->{str} && $R1->{$_."_Frd"}->{str} < 70 &&
		   $R2->{'PS2.0-primer'}->{end} < $R2->{$_."_Rvs"}->{str} && $R2->{$_."_Rvs"}->{str} < 70){

		    $R_seq->{$_} = substr $seq->{R1},$R1->{'PS1.0-primer'}->{end},$R1->{"DBU1-primer"}->{str}-$R1->{'PS1.0-primer'}->{end}+1-2;
		    $C_seq->{$_} = substr $seq->{R2},$R2->{'PS2.0-primer'}->{end},$R2->{"DBD2-primer"}->{str}-$R2->{'PS2.0-primer'}->{end}+1-2;
		    if(exists($R1->{$_."_Target"})){ #&& exists($R2->{"c".$_."_Target"})){
			#print Dumper $R1->{$_."_Target"};
			#print Dumper $R2->{"c".$_."_Target"};
			if( $R1->{$_."_Target"}->{Target_s} ==1 && $R1->{$_."_Target"}->{Target_e} == $target_len ){
			    #$R2->{"c".$_."_Target"}->{Target_s} ==1 && $R2->{"c".$_."_Target"}->{Target_e} == $target_len ){
			push @Target, $_;
			$locusR1 = $R1->{$_."_Target"}->{btop};
			#$locusR2 = $R2->{"c".$_."_Target"}->{btop};
			$orientation = 1;

			}}
		    if(exists($R2->{$_."_Target"}) ){#&& exists($R1->{"c".$_."_Target"})){
			if( $R2->{$_."_Target"}->{Target_s} ==1 && $R2->{$_."_Target"}->{Target_e}==$target_len){ #&&
                            #$R1->{"c".$_."_Target"}->{Target_s} ==1 && $R1->{"c".$_."_Target"}->{Target_e}==$target_len ){
			#print Dumper $R1->{"c".$_."_Target"};
			#print Dumper $R2->{$_."_Target"};
			$locusR2 = $R1->{$_."_Target"}->{btop};
                        #$locusR1 = $R2->{"c".$_."_Target"}->{btop};
			$orientation = -1;
			push @Target, $_;

			}}
		    #push @Target, $_;
		}}
	}


    }
    #print @Target;
    # Specific category assignment?
    my $category = 0;
    my $R_SEQ    = 0;
    my $C_SEQ    = 0;

    if( scalar @Target  == 1){
	$category = $Target[0];
	#print $category."\n";
	$R_SEQ = $R_seq->{$category};
	$C_SEQ = $C_seq->{$category};


    }






    return ($category,$P_SEQ1,$P_SEQ2,$R_SEQ,$C_SEQ,$locusR1,$locusR2,$orientation);
}







sub complement(){
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATGCatgc/TACGtacg/;
    return $seq;
}

sub get_bar2num(){
    my $file = shift;
    my $bar2num;
    open FILE,$file;
    while(<FILE>){
	chomp;
	my @array = split /\,/,$_;
	if($#array){
	    $bar2num->{$array[1]} = $array[0];
	}
    }
    close FILE;

    return $bar2num;
}


sub get_tag(){
    my $file = shift;
    my $tag_strct;
    #print $file;
    open FILE,$file;
    while(<FILE>){
	chomp;
	my @array = split /\,/,$_ ;
	#print @array;
	if($#array == 3){
	    $tag_strct->{$array[1]}->{$array[2]}->{$array[3]} = $array[0];
	}
    }
    close FILE;

    return $tag_strct;
}

sub get_seq(){
    my @files = @_;

    my $seq;
    for my $file (@files){

	my $tag = 'N_A';
	open FILE, $file;
	while(<FILE>){
	    chomp;
	    if($_ =~ /\>(.+)$/){
		$tag = $1;
	    }else{
		$seq->{$tag}->{seq} .= $_;
		$seq->{$tag}->{seq} =~ s/[^a-zA-Z]//g;
	    }
	}
	close FILE;

	while(my($tag,$value) = each %$seq){
	    $seq->{$tag}->{length} = length $seq->{$tag}->{seq};
	}
    }

    return $seq;
}
