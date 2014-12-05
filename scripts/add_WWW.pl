#!/usr/bin/perl
use strict;
use warnings;
my $link1 = "http://www.cbioportal.org/public-portal/cross_cancer.do?cancer_study_id=all&data_priority=0&case_ids=&gene_set_choice=user-defined-list&gene_list=";
my $link2 = "%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit";
my $link3 = "http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=";
my $link4 = "#histo";
open(FH, $ARGV[0]);
while(<FH>){
	chomp;
	my @D = split("\t", $_);
	if($D[6] =~ /Gene/){
		print "$_\tcBIO Link\tCOSMIC Link\n";
	}
	else{
		my $A = $link1.$D[6].$link2;
		my $B = $link3.$D[6].$link4;
		print "$_\t$A \t$B \n";
	}
}
