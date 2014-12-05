#!/usr/bin/perl
use strict;
use warnings;


open(FH, $ARGV[0]);
while(<FH>){
	chomp;
	my @data = split("\t", $_);
	if ($. == 1){
	}
	else{
		$data[0] =~ s/\//,/g;
		my @head = split(",", $data[0]);	
		if($data[10] =~ /EXON/){
			print "$head[0]\t$head[1]\t$head[3]\t$head[4]\t$data[13]\t$data[14]\n";
		}
	}
}
close FH;
