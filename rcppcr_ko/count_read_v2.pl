#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $data = do shift;


my $count;

while(my($category,$value1) = each %$data){
    while(my($P_TAG,$value2) = each %$value1){
	while(my($R_TAG,$value3) = each %$value2){
	    while(my($C_TAG,$value4) = each %$value3){
		while(my($read,$value5) = each %$value4){
		    while(my($strand,$value6) = each %$value5){
			my $btop = $value6;
			$count->{$category}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$strand}->{$btop} ||= 0;
			$count->{$category}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$strand}->{$btop} ++;
		    }




		}
	    }
	}
    }
}
print Dumper $count;
