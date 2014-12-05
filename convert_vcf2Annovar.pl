#!/usr/local/bin/perl
use warnings;
use strict;
#
# 12-04-2014
# This was written for converting the multiple sample vcf file into flat file containing the 
#   genotype, read counts and vaf for every sample
# Its parsing single alleleic event perfectly but there are some issues with multi-allelic calls.   
#


open (VAR, $ARGV[0]);
while (<VAR>) {
	chomp;
	if (m/^##/) {
		next;
	}
	my @field=split(/\t/,$_);
	@field >=8 or die "Error: invalid record found $_\n";
	my ($chr, $start, $ID, $ref, $alt, $quality_score, $filter, $info, $format, @sample) = @field;
	my ($end) = ($start);
	if ($chr =~ /^#CHR/i) {		#format specification line
		print "Chr\tStart\tEnd\tRef\tAlt\tQUAL\tFILTER\tINFO";
		for my $i (9 .. @field-1){
			print "\t$field[$i].Genotype\tTotalCoverage\tRefCoverage\tVarCoverage\tVariant Allele Freq";
		}
		print "\n";
		next;
	}
	if(length($ref) eq 1 and length($alt) eq 1 and $format =~ /GT:AD/){ # 1 base SNP
		print "$chr\t$start\t$start\t$ref\t$alt\t$quality_score\t$filter\t$info";
		for my $i (9 .. @field-1){
			my ($gt,$tot,$ref,$alt,$vaf) = processSNP($field[$i]);
			print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
		}
		print "\n";
	}
	elsif(length($ref) eq 1 and $alt =~ /,/ and length($alt) eq 3 and $format =~ /GT:AD/){ # BiAllele SNP
		my @allAllele = split(",", $alt);
		for my $all (@allAllele){
			print "$chr\t$start\t$start\t$ref\t$all\t$quality_score\t$filter\t$info";
			for my $i (9 .. @field-1){
				if ($field[$i] =~ /0\/0/ or $field[$i] =~ /0\/1/ or $field[$i] =~ /1\/1/ or $field[$i] =~ /\.\/\./){
					my ($gt,$tot,$ref,$alt,$vaf) = processSNP($field[$i]);
					print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
				}
				else{
					my ($gt,$tot,$ref,$alt,$vaf) = processMultiAllele($field[$i]); # Need Fixing
					print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
				}
			}
			print "\n";
		}
	}
	elsif(length($ref) >= length($alt) and $alt !~ /,/ and $format =~ /GT:AD/){
		my $head = substr ($ref, 0, length ($alt));
		my ($newstart, $newend);
		my ($newref, $newalt);
		if ($head eq $alt) {
			($newstart, $newend) = ($start+length ($head), $start + length ($ref)-1);
			($newref, $newalt) = (substr($ref, length($alt)), '-');
		}
		else {
			($newstart, $newend) = ($start, $start+length($alt)-1);
			($newref, $newalt) = ($ref, $alt);
		}
		print "$chr\t$newstart\t$newend\t$newref\t$newalt\t$quality_score\t$filter\t$info";
		for my $i (9 .. @field-1){
			my ($gt,$tot,$ref,$alt,$vaf) = processSNP($field[$i]);
			print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
		}
		print "\n";
	}
	elsif (length ($ref) < length ($alt) and $alt !~ /,/ and $format =~ /GT:AD/) {         #insertion or block substitution
		my $head = substr ($alt, 0, length ($ref));
		my ($newstart, $newend);
		my ($newref, $newalt);
		if ($head eq $ref) {
			($newstart, $newend) = ($start+length($ref)-1, $start+length($ref)-1);
			($newref, $newalt) = ('-', substr($alt, length($ref)));
		}
		else {
			($newstart, $newend) = ($start, $start+length($ref)-1);
			($newref, $newalt) = ($ref, $alt);
		}
		print "$chr\t$newstart\t$newend\t$newref\t$newalt\t$quality_score\t$filter\t$info";
		for my $i (9 .. @field-1){
			my ($gt,$tot,$ref,$alt,$vaf) = processSNP($field[$i]);
			print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
		}
		print "\n";
	}
	elsif(length($ref) >= length($alt) and $alt =~ /,/ and $format =~ /GT:AD/){
		my @allAllele = split(",", $alt);
		for my $all (@allAllele){
			my $head = substr ($ref, 0, length ($all));
			my ($newstart, $newend);
			my ($newref, $newalt);
			if ($head eq $all) {
				($newstart, $newend) = ($start+length ($head), $start + length ($ref)-1);
				($newref, $newalt) = (substr($ref, length($all)), '-');
			}
			else {
				($newstart, $newend) = ($start, $start+length($all)-1);
				($newref, $newalt) = ($ref, $all);
			}
			print "$chr\t$newstart\t$newend\t$newref\t$newalt\t$quality_score\t$filter\t$info";
			for my $i (9 .. @field-1){
				if ($field[$i] =~ /0\/0/ or $field[$i] =~ /0\/1/ or $field[$i] =~ /1\/1/ or $field[$i] =~ /\.\/\./){
					my ($gt,$tot,$ref,$alt,$vaf) = processSNP($field[$i]);
					print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
				}
				else{
					my ($gt,$tot,$ref,$alt,$vaf) = processMultiAllele($field[$i]);
					print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
				}
			}
			print "\n";
		}
	}
	elsif (length ($ref) < length ($alt) and $alt =~ /,/ and $format =~ /GT:AD/) {
		my @allAllele = split(",", $alt);
		for my $all (@allAllele){
			my $head = substr ($all, 0, length ($ref));
			my ($newstart, $newend);
			my ($newref, $newalt);
			if ($head eq $ref) {
				($newstart, $newend) = ($start+length($ref)-1, $start+length($ref)-1);
				($newref, $newalt) = ('-', substr($all, length($ref)));
			}
			else {
				($newstart, $newend) = ($start, $start+length($ref)-1);
				($newref, $newalt) = ($ref, $all);
			}
			print "$chr\t$newstart\t$newend\t$newref\t$newalt\t$quality_score\t$filter\t$info";
			for my $i (9 .. @field-1){
				if ($field[$i] =~ /0\/0/ or $field[$i] =~ /0\/1/ or $field[$i] =~ /1\/1/ or $field[$i] =~ /\.\/\./){
					my ($gt,$tot,$ref,$alt,$vaf) = processSNP($field[$i]);
					print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
				}
				else{
					my ($gt,$tot,$ref,$alt,$vaf) = processMultiAllele($field[$i]);
					print "\t$gt\t$tot\t$ref\t$alt\t$vaf";
				}
			}
			print "\n";
		}
	}
	else{
#		die || print "$_\n";
	}
}
sub processMultiAllele{
	my ($sampleInfo) = @_;
	my($gt,$tot,$ref,$alt,$vaf);
	my @temp = split(":", $sampleInfo);
	$gt = getGenotype($temp[0]);
	if($temp[1] =~ /\./ and length($temp[1]) eq 1){
		$gt = getGenotype($temp[0]);
		return ($gt, "0", "0", "0", "0");
	}
	elsif($temp[0] =~ /0\/2/){
		my @as = split(",", $temp[1]);
		$tot = $as[0] + $as[2];
		$ref = $as[0];
		$alt = $as[2];
		$vaf = getVAF($alt,$tot);
#		$vaf = "**";
	}
	elsif($temp[0] =~ /1\/2/){
		my @as = split(",", $temp[1]);
		$tot = $as[1] + $as[2];
		$ref = $as[1];
		$alt = $as[2];
		$vaf = getVAF($alt,$tot);
#		$vaf = "**";
	}
	elsif($temp[0] =~ /2\/2/){
                my @as = split(",", $temp[1]);
                $tot = $as[2];
                $ref = "0";
                $alt = $as[2];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /0\/3/){
                my @as = split(",", $temp[1]);
                $tot = $as[0] + $as[3];
                $ref = $as[0];
                $alt = $as[3];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /1\/3/){
                my @as = split(",", $temp[1]);
                $tot = $as[1] + $as[3];
                $ref = $as[1];
                $alt = $as[3];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /2\/3/){
                my @as = split(",", $temp[1]);
                $tot = $as[2] + $as[3];
                $ref = $as[2];
                $alt = $as[3];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /3\/3/){
		my @as = split(",", $temp[1]);
		$tot = $as[3];
		$ref = $as[0];
		$alt = $as[3];
		$vaf = getVAF($alt,$tot);
#		$vaf = "**";
	}
	elsif($temp[0] =~ /3\/4/){
                my @as = split(",", $temp[1]);
                $tot = $as[3] +$as[4];
                $ref = $as[3];
                $alt = $as[4];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /0\/4/){
                my @as = split(",", $temp[1]);
                $tot = $as[0] +$as[4];
                $ref = $as[0];
                $alt = $as[4];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /0\/5/){
                my @as = split(",", $temp[1]);
                $tot = $as[0] +$as[5];
                $ref = $as[0];
                $alt = $as[5];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /0\/6/){
                my @as = split(",", $temp[1]);
                $tot = $as[0] +$as[6];
                $ref = $as[0];
                $alt = $as[6];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /1\/4/){
                my @as = split(",", $temp[1]);
                $tot = $as[1] +$as[4];
                $ref = $as[1];
                $alt = $as[4];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /4\/5/){
                my @as = split(",", $temp[1]);
                $tot = $as[4] +$as[5];
                $ref = $as[4];
                $alt = $as[5];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /2\/6/){
                my @as = split(",", $temp[1]);
                $tot = $as[2] +$as[6];
                $ref = $as[2];
                $alt = $as[6];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /4\/6/){
                my @as = split(",", $temp[1]);
                $tot = $as[4] +$as[6];
                $ref = $as[4];
                $alt = $as[6];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /5\/5/){
                my @as = split(",", $temp[1]);
                $tot = $as[0] +$as[5];
                $ref = $as[0];
                $alt = $as[5];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /1\/6/){
                my @as = split(",", $temp[1]);
                $tot = $as[1] +$as[6];
                $ref = $as[1];
                $alt = $as[6];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /3\/6/){
                my @as = split(",", $temp[1]);
                $tot = $as[3] +$as[6];
                $ref = $as[3];
                $alt = $as[6];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /6\/6/){
                my @as = split(",", $temp[1]);
                $tot = $as[0] +$as[6];
                $ref = $as[0];
                $alt = $as[6];
 		$vaf = getVAF($alt,$tot);
#               $vaf = "**";
        }
	elsif($temp[0] =~ /5\/6/){
                my @as = split(",", $temp[1]);
                $tot = $as[5] +$as[6];
                $ref = $as[5];
                $alt = $as[6];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /4\/4/){
                my @as = split(",", $temp[1]);
                $tot = $as[0] +$as[4];
                $ref = $as[0];
                $alt = $as[4];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /2\/4/){
                my @as = split(",", $temp[1]);
                $tot = $as[2] +$as[4];
                $ref = $as[2];
                $alt = $as[4];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /2\/5/){
                my @as = split(",", $temp[1]);
                $tot = $as[2] +$as[5];
                $ref = $as[2];
                $alt = $as[5];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /3\/5/){
                my @as = split(",", $temp[1]);
                $tot = $as[3] +$as[5];
                $ref = $as[3];
                $alt = $as[5];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	elsif($temp[0] =~ /1\/5/){
                my @as = split(",", $temp[1]);
                $tot = $as[1] +$as[5];
                $ref = $as[1];
                $alt = $as[5];
		$vaf = getVAF($alt,$tot);
#                $vaf = "**";
        }
	else{
		print "\n\n$sampleInfo\n\n";
		die;
	}
	return ($gt, $tot, $ref, $alt, $vaf);
}
sub getVAF{
	my($alt,$tot) =@_;
	my $vaf;
	if($alt >0){
                $vaf = sprintf("%.2f",($alt/$tot));
                return ($vaf);
        }
        else{
                return ("0");
        }
	
}
sub processSNP{
	my ($sampleInfo) = @_;
	my($gt,$tot,$ref,$alt,$vaf);
	if($sampleInfo =~ /\.\/\./){
		return ("NotCovered", "0", "0", "0", "0");
	}
	else{
		my @temp = split(":", $sampleInfo);
		if($temp[1] =~ /\./ and length($temp[1]) eq 1){
			$gt = getGenotype($temp[0]);
			return ($gt, "0", "0", "0", "0");
		}
		$gt = getGenotype($temp[0]);
		$tot = getTotal($temp[1]);
		($ref,$alt,$vaf) = getAlleleCount($temp[1]);
		return ($gt, $tot, $ref, $alt, $vaf);
	}
}
sub getAlleleCount{
	my ($in) = @_;
	my @temp = split(",", $in);
	my $total =0;
	for (@temp){
		$total += $_;
	}
	if($temp[1] >0){
		my $vaf = sprintf("%.2f",($temp[1]/$total));
		return ($temp[0], $temp[1], $vaf);
	}
	else{
		return ($temp[0], $temp[1], "0");
	}
}

sub getTotal{
	my ($in) = @_;
	my @temp = split(",", $in);
	my $total =0;
	for (@temp){
		$total += $_;
	}
	return "$total";
}


sub getGenotype{
	my ($in) = @_;
	if($in =~ /\.\/\./){
		return "NotCovered";
	}
	elsif($in =~ /0\/1/){
		return "Het";
	}
	elsif($in =~ /1\/1/ or $in =~ /2\/2/){
		return "Homo";
	}
	elsif($in =~ /0\/0/){
		return "Ref";
	}
	elsif($in =~ /0\/2/ or $in =~ /1\/2/ or $in =~ /0\/3/ or $in =~ /1\/3/ or $in =~ /2\/3/ or $in =~ /3\/3/ or $in =~ /3\/4/ or $in =~ /0\/5/ or $in =~ /1\/5/ or $in =~ /2\/5/ or $in =~ /3\/5/ or $in =~ /4\/5/ or $in =~ /5\/5/ or $in =~ /0\/6/ or $in =~ /0\/4/ or $in =~ /4\/6/ or $in =~ /4\/5/ or $in =~ /1\/4/){
		return "*Het*";
	}
	else{
		return "*Het*";
	}
}
