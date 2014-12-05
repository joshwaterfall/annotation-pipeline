#!/usr/bin/perl
use strict;
use warnings;

if(!$ARGV[0] and !$ARGV[1] and !$ARGV[2]){
	print STDERR "Usage:\n";
	print STDERR "\t$0 tNB008.annotations tNB008 5\n";
	die;
}



open(ANN_FH, "$ARGV[0]");
my %ANNOVAR;
while(<ANN_FH>){
        chomp;
        my @local = split ("\t", $_);
	if($ARGV[2] ==5){
		my $key = join "\t", @local[0..4];
		my $end = @local - 1 ;
		my $value = join "\t", @local[5..$end];
		$ANNOVAR{"$key"} = $value;
	}
	elsif($ARGV[2] ==4){
		my $key = "$local[0]\t$local[1]\t$local[3]\t$local[4]";
                my $end = @local;
		my $value = join "\t", @local[5..$end];
                $ANNOVAR{"$key"} = $value;
	}
}
close ANN_FH;

open (ORI,"$ARGV[1]");

while (<ORI>){
        chomp;
	my @temp = split("\t", $_);
	my $val;
	my $vcf;	
	my $end = @temp - 1 ;
	if ($ARGV[2] ==5){
		$val = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
		$vcf = join "\t", @temp[5..$end];
	}
	elsif($ARGV[2] ==4){
		$val = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]";
		$vcf = join "\t", @temp[4..$end];
        }
	print "$val\t$ANNOVAR{$val}\t$vcf\n";
}
