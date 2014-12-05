#!/usr/bin/perl
use strict;
use warnings;

open(FH, "$ARGV[0]");
my %LIST;
while(<FH>){
	chomp;
	$LIST{$_} = "$_";
}
close FH;
open(FH,$ARGV[1]);
while(<FH>){
	chomp;
	my @d = split("\t",$_);
	if($d[6] =~ /Gene/){
		print "$_\t$ARGV[2]\n";
	}
	else {
		if(exists$LIST{$d[6]}){
			print "$_\tyes\n";
		}
		else{
			print "$_\t-\n";
		}
	}
}
close FH;
