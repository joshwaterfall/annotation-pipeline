#!/usr/bin/perl
use strict;
use warnings;
my $annovar =$ARGV[0];
my $sift = $ARGV[1];
my $pph = $ARGV[2];

my %SIFT;
open(FH, $sift);
while(<FH>){
        chomp;
        my @d =split("\t",$_);
#1       888659  T       C       TOLERATED       1
#1       1158631 A       G       N/A     N/A
#1       1249187 G       A       N/A     N/A
#1       1430985 C       T       DAMAGING        0.02
#1       1647871 T       C       N/A     N/A
	$SIFT{"chr$d[0]\t$d[1]\t$d[1]\t$d[2]\t$d[3]"} = "$d[4]\t$d[5]";
}
close FH;

my %PPH;
open(FH, $pph);
while(<FH>){
        chomp;
        my @d =split("\t",$_);
#chr1    13330627        G       A       benign  neutral 0.001
#chr1    13037860        A       G       benign  neutral 0
#chr1    13037864        T       A       benign  neutral 0
#chr1    15707775        A       G       benign  neutral 0.136
#chr1    15707775        A       G       unknown none
	if($_ !~ /unknown/){
		$PPH{"$d[0]\t$d[1]\t$d[1]\t$d[2]\t$d[3]"} = "$d[4]\t$d[5]\t$d[6]";
	}
}
close FH;

open(FH, $annovar);
while(<FH>){
	chomp;
	$_ =~ s/\t\t/\t0\t/g;
	$_ =~ s/\t\t/\t0\t/g;
	$_ =~ s/\t\t/\t0\t/g;
	$_ =~ s/\t\t/\t0\t/g;
	my @d =split("\t",$_);
	my $key ="$d[0]\t$d[1]\t$d[2]\t$d[3]\t$d[4]";
	if($key =~ /^Chr\t/i){
		print "$_\tSIFT Prediction\tSIFT Score\tPPH2 Prediction\tPPH2 Class\tPPH2 Probability\n";
	}
	else{
		print "$_\t";
		if(exists$SIFT{$key}){
			print "$SIFT{$key}\t";
		}
		else{
			print "-\t-\t";
		}
		if(exists$PPH{$key}){
			print "$PPH{$key}\n";
		}
		else{
			print "-\t-\t-\n";
		}
	}
}
close FH
