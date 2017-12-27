#!/usr/bin/perl -w

use Data::Dumper;

my @files = @ARGV;
#print @files;
my $big_data;
for my $file (@files){
    my $data = do $file;
    while(my($key1,$value1) = each %$data){
	while(my($key2,$value2) = each %$value1){
	    while(my($key3,$value3) = each %$value2){
		while(my($key4,$value4) = each %$value3){
		    while(my($key5,$value5) = each %$value4){
			while(my($key6,$value6) = each %$value5){
			    #print $key1, $key2, $key3, $key4, $key5, $key6;

			    $big_data->{$key1}->{$key2}->{$key3}->{$key4}->{$key5}->{$key6} ||=0;
			    $big_data->{$key1}->{$key2}->{$key3}->{$key4}->{$key5}->{$key6}+= $value6;

	       
				    
			    }
			}
		    }
		
	    }
	}
    }
}


print Dumper $big_data;
