#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: 505 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2012-05-17 10:54:12 -0700 (Thu, 17 May 2012) $';

our ($verbose, $help, $man);
our ($variantfile);
our ($outfile, $format, $includeinfo, $snpqual, $snppvalue, $coverage, $maxcoverage, $chr, $chrmt, $altcov, $allelicfrac, $fraction, $species, 
	$filterword, $confraction, $allallele, $withzyg);

our %iupac = (R=>'AG', Y=>'CT', S=>'CG', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-'); ### <<< FOR 5500SOLiD LifeScope ( S=>'GC' is replaced by S=>'CG')
our %iupacrev = reverse %iupac; ### <<< FOR 5500SOLiD LifeScope

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'outfile=s'=>\$outfile, 'format=s'=>\$format, 'includeinfo'=>\$includeinfo,
	'snpqual=f'=>\$snpqual, 'snppvalue=f'=>\$snppvalue, 'coverage=i'=>\$coverage, 'maxcoverage=i'=>\$maxcoverage, 'chr=s'=>\$chr, 'chrmt=s'=>\$chrmt, 
	'fraction=f'=>\$fraction, 'altcov=i'=>\$altcov, 'allelicfrac'=>\$allelicfrac,
	'species'=>\$species, 'filter=s'=>\$filterword, 'confraction=f'=>\$confraction, 'allallele!'=>\$allallele, 'withzyg'=>\$withzyg) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($variantfile) = @ARGV;

$chrmt ||= 'M';

if (not $format) {
	$format = 'pileup';
	print STDERR "NOTICE: the default --format argument is set as 'pileup'\n";
}

if (defined $outfile) {
	open (STDOUT, ">$outfile") or die "Error: cannot write to output file $outfile: $!\n";
}

defined $snpqual and $format eq 'pileup' || $format eq 'vcf4' || pod2usage ("Error in argument: the --snpqual is supported only for the 'pileup' or 'vcf4' format");
if (not defined $snpqual and $format eq 'pileup') {
	$snpqual = 20;
	print STDERR "NOTICE: the default --snpqual argument for pileup format is set as 20\n";
}

if (not defined $snppvalue) {
	$snppvalue = 1;		#by default, do not use any of the P-value cutoff in filtering out GFF3-SOLID files (this is differnt from handling pileup files)
}

if (not defined $coverage) {
	$coverage = 0;
}

if (defined $fraction) {
	$format eq 'pileup' or $format eq 'vcf4' or pod2usage ("Error in argument: the '--fraction' argument is supported for the pileup or vcf4 format only");
	$format eq 'vcf4' and print STDERR "NOTICE: the --fraction argument works ONLY on indels for vcf4 format\n";
	$fraction >= 0 and $fraction <=1 or pod2suage ("Error in argument: the --fraction argument must be between 0 and 1 inclusive");
} else {
	$fraction = 0;
}

if (defined $confraction) {
	$format eq 'vcf4' and print STDERR "NOTICE: the --confraction argument works ONLY on indels for vcf4 format\n";
	$confraction >= 0 and $fraction <=1 or pod2suage ("Error in argument: the --confraction argument must be between 0 and 1 inclusive");
} else {
	$confraction = 0;
}

if (defined $altcov) {
	$format eq 'pileup' or pod2usage ("Error in argument: the '--altcov' argument is supported for the '--format pileup' only");
	$altcov < $coverage or pod2suage ("Error in argument: the --altcov argument must be less than --coverage");
	$altcov > 0 or pod2suage ("Error in argument: the --altcov argument must be a positive integer");
}

if ($allallele) {
	$format eq 'vcf4' or pod2usage ("Error in argument: the '--allallele' argument is only supported for the '--format vcf4'");
}

if ($format eq 'pileup') {
	convertPileup ($variantfile);
} elsif ($format eq 'cg') {
	convertCG ($variantfile);
} elsif ($format eq 'soap') {
	print STDERR "WARNING: the support for '--format soap' is not well developed yet and may contain bugs for indel analysis.\n";
	convertSOAP ($variantfile);
} elsif ($format eq 'maq') {
	print STDERR "WARNING: the support for '--format maq' is not well developed yet and may contain bugs.\n";
	convertMAQSNP ($variantfile);
} elsif ($format eq 'casava') {
	if (not defined $chr) {
		pod2usage ("Error in argument: please supply --chr argument for the '--format casava'");
	}
	convertCASAVA ($variantfile, $chr);
} elsif ($format eq 'vcf4') {
	convertVCF4 ($variantfile);
} elsif ($format eq 'annovar') {
	convertANNOVAR ($variantfile);
} else {
	pod2usage ("Error in argument: the --format $format is not currently supported. Please contact ANNOVAR developer for adding the support");
}


sub convertPileup {
	my ($variantfile) = @_;
	my ($countline, $countvar, $counthom, $counthet, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0/;
	
	if ($variantfile eq 'stdin') {
		*VAR = *STDIN;
	} else {
		open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	}
	print STDERR "NOTICE: Column 6-9 in output are heterozygosity status, SNP quality, total reads, reads with mutation\n";

	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		my $hom = 'hom';
		my @field = split (/\t/, $_);
		@field >= 10 or die "Error: invalid record found in pileupfile $variantfile (at least 10 fields expected): <$_>\n";
		my ($chr, $pos, $wt, $call, @other) = @field;
		my ($cons_qual, $snp_quality, $readcount, $readallele) = @field[4,5,7,8];
		$chr =~ s/^chr//;
		$wt = uc $wt;					#wt may or may not in upper case, it depends on the input FASTA file
		$call = uc $call;				#usually call is in upper case
		$readallele = uc $readallele;			#lower case indicate the opposite strand
		
		$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
		
		$snp_quality >= $snpqual or next;		#quality of the variant call
		$readcount >= $coverage or next;		#read count of the variant call
		$maxcoverage and $readcount <= $maxcoverage || next;	#maximum coverage of the variant call
		
		if ($wt eq '*') {				#indel
			#example:
			#1       970271  *       +C/+C   39      106     44      5       +C      *       1       4       0       0       0
			#1       1548977 *       */+CCG  29      29      42      3       *       +CCG    2       1       0       0       0
			#1       1674810 *       */+C    24      24      42      6       *       +C      5       1       0       0       0
			#1       968466  *       -CT/-CT 53      339     55      5       -CT     *       5       0       0       0       0
			#1       1093600 *       -GAAA/* 29      29      53      3       -GAAA   *       1       2       0       0       0
			#1       1110101 *       */-A    41      41      17      6       *       -A      5       1       0       0       0
			#1       1215395 *       */-TC   26      26      32      4       *       -TC     3       1       0       0       0
			my @obs = split (/\//, $call);		#example: '+AG/+AG' as homozygotes, '*/-TA' or '*/+T' as heterozygotes
			@obs == 2 or die "Error: pileup record contains more than two alternative alleles: <$_>\n";
			my ($end, $ref, $alt);
			my ($indelreadcount);			#number of reads supporting the indel
			
			
			if ($obs[0] eq $obs[1]) {
				#something weird may occur in SamTools output: 22      15231121        *       */*     360     0       32      158     *       +GA     156     2       0       0       0
				$obs[0] eq '*' and next;	
	
				#for deletions, SAMTOOLS represent deletion using a location before the first deleted base in the reference sequence coordinate system
				#for example, a deletion in Samtools output is "1       109266688       *       */-CTT  1429    1429    58      43      *       -CTT    24      19      0       0       0"
				#the correct ANNOVAR input (for rs35029887) should be "1       109266689       109266691       CTT     -       het     1429"
				#insertions are fine without change; for example, check rs5745466 in Genome Browser; samtools report "1       76119508        *       +AT/+AT"
				#for this insertion, ANNOVAR input file (for rs5745466) becomes "1       76119508        76119508        -       AT      hom     1601"

				if ($obs[0] =~ m/^\-/) {
					$pos++;			#add 1 to position in deletion
				}
				
				$indelreadcount = calculateIndelReadCount ($obs[0], \@field);
				$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $indelreadcount >= $altcov || next;
				
				if ($chr eq $chrmt or $allelicfrac) {
					$hom = sprintf ("%.3f", $indelreadcount/$readcount);
				}
				($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[0]);
				print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
				$counthom++;
			} else {
				$hom = 'het';
				if ($obs[0] =~ m/^[\-\+]/) {
					$obs[0] =~ m/^\-/ and $pos++;
					($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[0]);
					$indelreadcount = calculateIndelReadCount ($obs[0], \@field);
					$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
					defined $altcov and $indelreadcount >= $altcov || next;
					
					if ($chr eq $chrmt or $allelicfrac) {
						$hom = sprintf ("%.3f", $indelreadcount/$readcount);
					}
					print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
					$counthet++;
				}
				if ($obs[1] =~ m/^[\-\+]/) {
					$obs[1] =~ m/^\-/ and $pos++;
					($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[1]);
					$indelreadcount = calculateIndelReadCount ($obs[1], \@field);
					$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
					defined $altcov and $indelreadcount >= $altcov || next;
					
					if ($chr eq $chrmt or $allelicfrac) {
						$hom = sprintf ("%.3f", $indelreadcount/$readcount);
					}
					print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
					$counthet++;
				}
			}
			$countindel++;
		} else {
			#1       798494  G       A       36      36      58      3       AAA     bbb
			#1       798677  T       K       33      33      52      26      ,$.,,G.GG,.,......,..G,,...     b`bbBaJIbFbZWaTNQbb_VZcbbb
			#1       856182  G       A       42      42      50      5       AAAAA   B\bbb
			#1       861034  A       M       48      48      49      14      ,$,.,..,cc.c.,c bbBbb`]BFbHbBB
			#1       864289  T       K       22      22      56      6       .g,,g,  BH^_BB
			
			$wt eq $call and next;			#this is not a SNP
			my $obs = $iupac{$call} or die "Error: invalid best call ($call) in <$_>\n";
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] ne $obs[1]) {
				$hom = 'het';
			}
				
			
			if ($obs[0] eq $wt) {			#obs[0] is guaranteed to be an alternative allele
				@obs = @obs[1,0];
			}
			if ($wt eq 'A' and $obs[0] eq 'G' or $wt eq 'G' and $obs[0] eq 'A' or $wt eq 'C' and $obs[0] eq 'T' or $wt eq 'T' and $obs[0] eq 'C') {
				unless ($wt ne $obs[0] and $wt ne $obs[1] and $obs[0] ne $obs[1]) {
					$countti++;
				}
				
			} else {
				unless ($wt ne $obs[0] and $wt ne $obs[1] and $obs[0] ne $obs[1]) {
					$counttv++;
				}
			}
			
			my $mutallelecount;
			
			if ($obs[1] eq $wt) {			#het SNP
				if ($chr eq $chrmt or $allelicfrac) {
					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
				}
				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
				
				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
				$counthet++;
			} elsif ($obs[1] ne $obs[0]) {		#het SNP but both differ from reference allele
				if ($chr eq $chrmt or $allelicfrac) {
					$hom = calculateAllelicFraction ($obs[1], $field[8], $readcount);
				}
				$mutallelecount = calculateMutAlleleCount ($obs[1], $readallele);
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
				
				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[1], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
				if ($chr eq $chrmt) {
					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
				}
				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
				
				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
				$counthet++;
				$counthet++;
			} else {				#homo SNP
				if ($chr eq $chrmt or $allelicfrac) {
					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
				}
				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
				
				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
				$counthom++;
			}
			$countsnp++;
		}
		$countvar++;
	}
	my $triallelic = $countsnp-$countti-$counttv;
	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
	print STDERR "NOTICE: Among ${\($counthet+$counthom)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
}

sub calculateIndelReadCount {
	my ($obs, $field) = @_;
	#make sure to use upper case in the comparison, for example:
	#chr10   130133  *       */-ca   189     189     59      31      *       -ca     27      4       0       0       0
	if ($obs eq uc $field->[8]) {
		return $field->[10];
	} elsif ($obs eq uc $field->[9]) {
		return $field->[11];
	} else {
		die "Error: invalid record in pileup file (indel counts cannot be inferred): <$obs> vs <@$field>\n";
	}
}

sub calculateMutAlleleCount {
	my ($allele, $string) = @_;	#they should be already both in upper case
	$string =~ s/\^.//g;		#^ is followed by mapping quality
	$string =~ s/\$//g;
	$string =~ s/[+-]1[^\d]//g;	#1 followed by a non-number
	$string =~ s/[+-]2..//g;
	$string =~ s/[+-]3...//g;
	$string =~ s/[+-]4....//g;
	$string =~ s/[+-]5.....//g;
	$string =~ s/[+-]6......//g;
	$string =~ s/[+-]7.......//g;
	$string =~ s/[+-]8........//g;
	$string =~ s/[+-]9.........//g;
	$string =~ s/[+-]10..........//g;
	
	#make sure to use upper case letter
	my @string = split (//, uc $string);
	my $count = 0;
	for my $i (0 .. @string-1) {
		$allele eq $string[$i] and $count++;
	}
	return $count;
}

sub calculateAllelicFraction {
	my ($obs, $readbase, $readcount) = @_;
	my @readbase = split (//, $readbase);
	my $count=0;
	for my $i (0 .. @readbase-1) {
		uc $obs eq uc $readbase[$i] and $count++;
	}
	my $hom = $count/$readcount;
	length ($hom) > 5 and $hom > 0.001 and $hom = sprintf ("%.3f", $hom);
	return $hom;
}

sub recalculateEndRefObs {		#recalculate end position, reference allele and observed allele
	my ($end, $ref, $obs) = @_;
	if ($obs =~ m/^\-(\w+)/) {	#deletion
		$end += (length ($1)-1);
		$ref = $1;
		$obs = '-';
	} elsif ($obs =~ m/^\+(\w+)/) {	#insertion
		$ref = '-';
		$obs = $1;
	} else {
		die "Error: cannot handle $end, $ref, $obs\n";
	}
	return ($end, $ref, $obs);
}

sub convertCG { # Written By Rajesh Patidar rajbtpatidar@gmail.com
	my ($variantfile) = @_;
	
	unless (open (VAR, $variantfile)){
		die "Error: cannot read from variant file $variantfile: $!\n";
		print STDERR "NOTICE: Converting variants from $variantfile\n";
	}
#	>locus	ploidy	chromosome	begin	end
#	zygosity	varType	reference	allele1Seq	allele2Seq allele1VarScoreVAF	
#	allele2VarScoreVAF	allele1VarScoreEAF	allele2VarScoreEAF	allele1VarQuality	allele2VarQuality
#	allele1HapLink	allele2HapLink	allele1XRef	allele2XRef evidenceIntervalId	
#	allele1ReadCount	allele2ReadCount	referenceAlleleReadCount	totalReadCount	allele1Gene
#	allele2Gene	pfam	miRBaseId	repeatMasker segDupOverlap
#	relativeCoverageDiploid	calledPloidy	relativeCoverageNondiploid	calledLevel	relativeCoverageSomaticNondiploid
#	somaticCalledLevel	bestLAF	lowLAF	highLAF	allele1ReadCount-N1
#	allele2ReadCount-N1	referenceAlleleReadCount-N1	totalReadCount-N1	locusDiffClassification	somaticCategory
#	somaticRank	somaticScore	somaticQuality
#728     2       chr1    38906   38907   hom     snp     C       T       T       110     172     110     172     VQHIGH  VQHIGH                  dbsnp.108:rs3874156;dbsnp.131:rs75829199;dbsnp.86:rs806726      dbsnp.108:rs3874156;dbsnp.131:rs75829199;dbsnp.86:rs806726      355     4       4       5       9       645520:NR_026818.1:FAM138A:TSS-UPSTREAM:UNKNOWN-INC     645520:NR_026818.1:FAM138A:TSS-UPSTREAM:UNKNOWN-INC                     MLT1E1A-int:ERVL-MaLR:38.5      6       1.14    N       1.06    1.000   0.77    1.029   0.21    0.18    0.26    8       8       3       11      alt-identical;alt-identical
	my $foundheader=0;
	my $countline = 0;
	my $prechr =0;
	my $prestart=0;
	my $preend=0;
	my $preref=0;
	my $preobs=0;
	my $prevartype=0;
	
	while (<VAR>) {
		chomp $_;
		if (m/^>locus/) {
			$foundheader++;
		}
		if (not $foundheader) {
			$countline > 50 and die "Error: invalid CG-var file format for $variantfile (>locus record is not found within the first 50 lines)\n";
			next;
		}
		if($_ =~ /^>locus/){
			print "Chrt\tStart\tEnd\tRef\tObs\t$_\n";
		}
		my $info = $_;
		my @line = split (/\t/, $_);
		my ($locus, $ploidy, $chr, $start, $end, $zygocity, $vartype, $ref, $obs) = split (/\t/, $_);
		$vartype eq 'ins' or $start++;		
		#CG use zero-start, half-open coordinate. Insertion does not need to be processed 
		#example, "16476   2       2       chr1    751820  751820  ins             T       49              dbsnp:rs59038458"
		$obs eq '' and $obs = '-';
		$ref eq '' and $ref = '-';
		if ($vartype =~ m/^snp|ins|del|delins|sub$/) {		#in new versions of the files, "sub" is used instead of "delins".
			if ($chr eq $prechr and $start eq $prestart and $end eq $preend and $obs eq $preobs) {		#homozygous mutation
				print "$chr\t$start\t$end\t$ref\t$obs\t$info\n";
				($prechr, $prestart, $preend, $prevartype, $preref, $preobs) = qw/0 0 0 0 0 0/;
			} else {
				if ($prestart and $preend) {
					print "$prechr\t$prestart\t$preend\t$preref\t$preobs\t$info\n";
				}
				($prechr, $prestart, $preend, $prevartype, $preref, $preobs) = ($chr, $start, $end, $vartype, $ref, $obs);
			}
		}
	}
}


sub convertSOAP {
	my ($variantfile) = @_;
	my ($countline, $countvar, @other);
	
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		
		my @field = split (/\t/, $_);
		if (@field == 18) {		#snp file
			my ($chr, $pos, $wt, $call, @other) = @field;
			$chr =~ s/^chr//;
	
			$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
	
			my $obs = $iupac{$call} or die "Error: invalid best call in <$_>\n";
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] eq $wt and $obs[1] eq $wt) {
				die "Error: reference alleles are identical to observed alleles: <$_>\n";
			} elsif ($obs[0] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] ne $obs[0]) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
				$countvar++;
			} else {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
			}
		} elsif (@field == 6) {		#indel file
			my ($chr, $pos, $strand, $indellen, $call, $homo) = @field;
			$homo eq 'Homo' and $homo = 'hom';
			$homo eq 'Hete' and $homo = 'het';
			$chr =~ s/^chr//;
	
			$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
	
			if ($indellen =~ m/^\+(\d+)$/) {		#insertion
				length ($call) == $1 or die "Error: invalid record found in SOAPindel file: <$_>\n";
				print join("\t", $chr, $pos, $pos, '-', $call, $homo), "\n";
			} elsif ($indellen =~ m/^\-(\d+)$/) {		#deletion
				length ($call) == $1 or die "Error: invalid record found in SOAPindel file: <$_>\n";
				print join("\t", $chr, $pos, $pos+$1-1, $call, '-', $homo), "\n";
			} else {
				die "Error: invalid record found in SOAPindel file: <$_>\n";
			}
		} else {
			die "Error: invalid record found in $variantfile (18 or 6 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
		}
		$countvar++;
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
}


sub convertANNOVAR {
	my ($variantfile) = @_;
	my ($countline, $countvar, $countsnp, $invalid) = (0, 0, 0, 0);
	my ($countti, $counttv) = (0, 0);
	
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	while (<VAR>) {
		$countline++;
		s/[\r\n]+$//; 
		my @field = split (/\s+/, $_);
		@field >= 5 or die "Error: invalid record found in annovar input file (at least 5 tab or space delimited fields expected): <$_>\n";
		my ($chr, $start, $end, $ref, $obs) = @field;
		if ($ref =~ m/^[ACGT]$/ and $obs =~ m/^[ACGT]$/) {
			if ($ref eq 'A' and $obs eq 'G' or $ref eq 'G' and $obs eq 'A' or $ref eq 'C' and $obs eq 'T' or $ref eq 'T' and $obs eq 'C') {
				$countti++;
			} else {
				$counttv++;
			}
			$countsnp++;
		}
		
		($ref, $obs) = (uc $ref, uc $obs);
		$chr =~ s/^chr//;
		if ($chr =~ m/[^\w\.]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {		#chr name could contain . (example: GL000212.1)
			$invalid++;
		} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
			or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
			or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
			or $start =~ m/[^\d]/ 			#start is not a number
			or $end =~ m/[^\d]/ 			#end is not a number
			or $start > $end			#start is more than end
			or $ref ne '0' and $end-$start+1 != length ($ref) 	#length mismatch with ref
			or $ref eq '-' and $start != $end	#length mismatch for insertion
			) {
			$invalid++;
		}
		print "$_\n";
		$countvar++;
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
	$invalid and print STDERR "WARNING: $invalid input lines have invalid formats\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions\n";
}

sub convertMAQSNP {
	my ($variantfile) = @_;
	my ($countline, $countvar, @other);
	
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		
		my @field = split (/\t/, $_);
		my @other = ();
		if (@field == 12) {					#SNPs
			my ($chr, $pos, $wt, $call, @other) = @field;
			$chr =~ s/^chr//;
	
			$includeinfo and @other = @field;			#unless -includeinfo is set, the other will not be printed
	
			my $obs = $iupac{$call} or die "Error: invalid best call in <$_>\n";
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] eq $wt and $obs[1] eq $wt) {
				die "Error: reference alleles are identical to observed alleles: <$_>\n";
			} elsif ($obs[0] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] ne $obs[0]) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
				$countvar++;
			} else {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
			}
			$countvar++;
		} elsif (@field == 13) {				#indels; the deletion start site do not need changes; the duplication start site need additional investigation by ANNOVAR developers
			my ($chr, $pos, $type, $numread, $call, @other) = @field;
			$chr =~ s/^chr//;
	
			$includeinfo and @other = @field;			#unless -includeinfo is set, the other will not be printed
	
			my @obs = split (/:/, $call);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] =~ m/^\-(\d+)/) {		#deletion
				my $len = $1;
				print $chr, "\t", $pos, "\t", $pos+$len-1, "\t", $obs[1], "\t", '-', "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[0] =~ m/^(\d+)/) {	#insertion
				my $len = $1;
				print $chr, "\t", $pos-1, "\t", $pos-1, "\t", '-', "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";	#2011jul12: changed pos to pos-1 for insertions
			}
			$countvar++;
		} else {
			die "Error: invalid record found in $variantfile (12 or 13 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
		}
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
}

sub convertCASAVA {
	my ($variantfile, $chr) = @_;
	my ($countline, $countvar, @other);
	
	my ($intype);
	my ($pos_index, $call_index, $reference_index, $type_index, $score_index, $total_index, $used_index);
	my ($ref_indel_index, $quality_index, $maxgtype_index, $bp1_reads_index, $ref_reads_index, $indel_reads_index, $other_reads_index);
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		my @field;

		if (m/^#/) {
			s/^#//;
			if (s/^\$\sCOLUMNS\s//) {
				@field = split (/\s+/, $_);
			} else {
				@field = split (/\t/, $_);
			}
			if (m/\bposition\b/ or m/\bpos\b/) {
				for my $i (0 .. @field-1) {
					if ($field[$i] eq 'position' or $field[$i] eq 'pos') {
						$pos_index = $i;
					} elsif ($field[$i] eq 'max_gt|poly_site') {		#this has priority over max_gt per se
						$intype = 'snp';
						print STDERR "NOTICE: Automatically detected input type as $intype\n";
						$call_index = $i;
					} elsif ($field[$i] eq 'modified_call' or $field[$i] eq 'max_gt') {
						$intype = 'snp';
						$call_index = $i;
					} elsif ($field[$i] eq 'reference' or $field[$i] eq 'ref') {
						$reference_index = $i;
					} elsif ($field[$i] eq 'type') {
						$type_index = $i;
					} elsif ($field[$i] eq 'score') {
						$score_index = $i;
					} elsif ($field[$i] eq 'total') {
						$total_index = $i;
					} elsif ($field[$i] eq 'used') {
						$used_index = $i;
					} elsif ($field[$i] eq 'ref/indel') {
						$intype = 'indel';
						print STDERR "NOTICE: Automatically detected input type as $intype\n";
						$ref_indel_index = $i;
					} elsif ($field[$i] eq 'Q(max_gt|poly_site)') {	#this has priority over Q(max)
						$quality_index = $i;
					} elsif ($field[$i] eq 'Q(max_gt)') {	#this has priority over Q(max)
						$quality_index = $i;
					} elsif ($field[$i] eq 'Q(indel)') {
						$quality_index = $i;
					} elsif ($field[$i] eq 'max_gtype') {
						$maxgtype_index = $i;
					} elsif ($field[$i] eq 'bp1_reads') {
						$bp1_reads_index = $i;
					} elsif ($field[$i] eq 'ref_reads') {
						$ref_reads_index = $i;
					} elsif ($field[$i] eq 'indel_reads') {
						$indel_reads_index = $i;
					} elsif ($field[$i] eq 'other_reads') {
						$other_reads_index = $i;
					}
				}
			}
			next;
		}
		
		##$ COLUMNS seq_name pos bcalls_used bcalls_filt ref Q(snp) max_gt Q(max_gt) max_gt|poly_site Q(max_gt|poly_site) A_used C_used G_used T_used
		#chr21.fa	9411785	1	0	G	11	GT	3	GT	3	0	0	0	1
		#chr21.fa	9414658	1	0	T	10	CT	3	CT	3	0	1	0	0
		#chr21.fa	9415181	2	0	C	52	TT	5	TT	5	0	0	0	2
		#chr21.fa	9415317	2	0	C	6	CT	6	CT	34	0	1	0	1

		
		$intype or die "Error: unable to recognize the correct type of the input file (make sure that header line is present in $variantfile)\n";
		@field = split (/\t/, $_);
		
		if ($intype eq 'snp') {					#SNPs
			defined $pos_index and defined $reference_index and defined $call_index or die "Error: unalbe to find the position, reference and modified_call column header in $variantfile\n";
			my ($pos, $wt, $obs) = @field[$pos_index, $reference_index, $call_index];
			my (@other);
			defined $pos and defined $wt and defined $obs or die;
			$includeinfo and @other = @field;
			
			length ($obs) == 1 and $obs .= $obs;
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed allele $obs should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] eq $wt and $obs[1] eq $wt) {
				die "Error: reference alleles are identical to observed alleles: <$_>\n";
			} elsif ($obs[0] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] ne $obs[0]) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
				$countvar++;
			} else {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
			}
			$countvar++;
		} elsif ($intype eq 'indel') {				#indels
			defined $pos_index and defined $ref_indel_index and defined $maxgtype_index or die "Error: unable to find the pos, ref_indel and max_gtype column header in $variantfile\n";
			my ($pos, $call, $hom, @other) = @field[$pos_index, $ref_indel_index, $maxgtype_index];
			$includeinfo and @other = @field;

			#hg19 coordinate below; insertion needs position adjustment!!! deletion is fine
			#948847  1I      CCTCAGGCTT      -/A     ATAATAGGGC      969     hom     47      het     22      0       16      6       A       1       2
			#978604  2D      CACTGAGCCC      CT/--   GTGTCCTTCC      251     hom     20      het     8       0       4       4       CT      1       0
			#1276974 4I      CCTCATGCAG      ----/ACAC       ACACATGCAC      838     hom     39      het     18      0       14      4       AC      2       4
			#1289368 2D      AGCCCGGGAC      TG/--   GGAGCCGCGC      1376    hom     83      het     33      0       25      9       TG      1       0
			#185137455     11I10M2I        TATGTGTCCT      -----------TTTTTTATTT--/AAATGATAGACTTTTTTTTTTAA ATTTCAGAAA      1126    het     988     hom    45       20      24      7       N/A     0       0
			#1276931 2D41M4I CACACACATG      CACACACACGCACACACGTGCAATGTGAAAACACCTCATGCAG----/--CACACACGCACACACGTGCAATGTGAAAACACCTCATGCAGACAC ACACATGCAC      548     hom     16      het     8       0       11      11      N/A     0       0
			
			my @obs = split (/\//, $call);
			@obs == 2 or die "Error: observed indel allele $call should correspond to two alleles: <$_>\n";
			if ($obs[0] =~ m/^\-+$/) {		#insertion
				my $len = length ($obs[0]);
				print $chr, "\t", $pos-1, "\t", $pos-1, "\t", '-', "\t", $obs[1], "\t", $hom, "\t", join ("\t", @other), "\n";
			} elsif ($obs[1] =~ m/^\-+$/) {		#deletion
				my $len = length ($obs[0]);
				print $chr, "\t", $pos, "\t", $pos+$len-1, "\t", $obs[0], "\t", '-', "\t", $hom, "\t", join ("\t", @other), "\n";
			} elsif (length ($obs[0]) eq length ($obs[1])) {	#block substitution
				$obs[0] =~ s/\-//g;
				$obs[1] =~ s/\-//g;
				print $chr, "\t", $pos, "\t", $pos+length($obs[0])-1, "\t", $obs[0], "\t", $obs[1], "\t", $hom, "\t", join ("\t", @other), "\n";
			} else {
				die "Error: invalid record found in indel line: <$_>\n";
			}
			$countvar++;
		} else {
			die "Error: invalid record found in $variantfile (11 or 15 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
		}
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
}

sub convertVCF4 {
	my ($variantfile) = @_;
	
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0 0/;
	
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		$countline++;
		
		if (m/^##fileformat=VCFv(\d+\.)/) {
			$1<4 and print STDERR "ERROR: Input file is not in VCF version 4 format but is $_" and exit;
		}
		if (m/^##UnifiedGenotyper/) {
			$source_program = 'gatksnp';
			print STDERR "NOTICE: Detected that the VCF4 file is generated by GATK UnifiedGenotyper\n";
			$includeinfo or print STDERR "NOTICE: column 6-10 represent heterozygosity status, quality score, read depth, RMS mapping quality, quality by depth\n";
			$fraction and print STDERR "WARNING: the --fraction argument will be ignored for GATK SNP calls!!!\n";
			$confraction and print STDERR "WARNING: the --confraction argument will be ignored for GATK SNP calls!!!\n";
		}
		if (m/^##IndelGenotyper/) {
			$source_program = 'gatkindel';
			print STDERR "NOTICE: Detected that the VCF4 file is generated by GATK IndelGenotyper\n";
			$includeinfo or print STDERR "NOTICE: column 6-10 represent heterozygosity status, quality score, read depth, read count supporting indel call, RMS mapping quality\n";
		}
		
		if (not m/^#/ and not $source_program) {	#finished reading header line but did not detect the source program
			$includeinfo or print STDERR "NOTICE: for SNPs, column 6 and beyond MAY BE heterozygosity status, quality score, read depth, RMS mapping quality, quality by depth, if these information can be recognized automatically\n";
			$includeinfo or print STDERR "NOTICE: for indels, column 6 and beyond MAY BE heterozygosity status, quality score, read depth, read count supporting indel call, RMS mapping quality, if these information can be recognized automatically\n";
			$source_program = 'unknown';
		}
		
		m/^##/ and next;		#skip comment lines
		s/[\r\n]+$//;		#delete trailing new lines
		my @Rajesh = split("\t",$_);
		my $lastelement = $#Rajesh;
		if($_ =~ /#CHROM/){
			print "Chr\tStart\tEnd\tRef\tAlt\t";
			my $otherinfo = join "\t", @Rajesh[5..$lastelement];
			print "$otherinfo\n";
		}
		else{
		my $otherinfo = join "\t", @Rajesh[5..$lastelement];	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
	
		#format description: http://www.1000genomes.org/wiki/Analysis/vcf4.0
		#standard VCF4 should have 8 columns, but some software may produce more columns (for example, for genotype calls). The first 8 columns should follow the specification
		
		#example of VCF4 generated by GATK SNP caller
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
		#1       55      .       T       G       34.82   .       DP=2;Dels=0.00;HRun=0;HaplotypeScore=0.00;MQ=14.16;MQ0=0;QD=17.41;SB=-10.00     GT:DP:GL:GQ     0/1:1:-6.66,-0.30,-0.00:1.76
		#1       2646    .       G       A       40.91   .       DP=4;Dels=0.00;HRun=0;HaplotypeScore=0.00;MQ=7.50;MQ0=3;QD=10.23;SB=-10.00      GT:DP:GL:GQ     0/1:1:-7.27,-0.30,-0.00:1.76
		
		#example of VCF4 generated by GATK indel caller
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
		#1       2525324 .       G       GC      .       PASS    AC=5,5;DP=12;MM=4.8,3.7142856;MQ=29.0,42.285713;NQSBQ=33.0,46.463768;NQSMM=0.24,0.20289855;SC=0,5,1,6  GT       0/1
		#1       3553372 .       GC      G       .       PASS    AC=6,6;DP=6;MM=0.8333333,0.0;MQ=60.0,0.0;NQSBQ=63.533333,0.0;NQSMM=0.0,0.0;SC=0,6,0,0   GT      1/0
		#1       6093011 .       CG      C       .       PASS    AC=31,31;DP=32;MM=0.7096774,2.0;MQ=59.64516,60.0;NQSBQ=64.192184,39.666668;NQSMM=0.0,0.11111111;SC=23,8,0,1     GT      1/0
		
		#example of VCF4 generated by 1000G
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		#1       533     .       G       C       .       PASS    AA=.;AC=6;AN=120;DP=423
		#1       41342   .       T       A       .       PASS    AA=.;AC=29;AN=120;DP=188
		#1       41791   .       G       A       .       PASS    AA=.;AC=5;AN=120;DP=192
		#1       44449   .       T       C       .       PASS    AA=C;AC=2;AN=120;DP=166
		#1       44539   rs2462492       C       T       .       PASS    AA=T;AC=2;AN=120;DP=131    
		
		#example of VCF4 generated by 1000G
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		#1       1000153 .       TCACA   T       100     PASS    AF=0.115095;HP=1;NF=16;NR=13;NS=52;CA=0;DP=615
		#1       1000906 .       CA      C       48      PASS    AF=0.0772696;HP=1;NF=2;NR=9;NS=51;CA=0;DP=281
		#1       1000950 rs60561655;-/G  CG      C       100     PASS    AF=0.447771;HP=5;DB;NF=10;NR=20;NS=50;CA=M;DP=291
		#1       1010786 rs36095298;-/G,mills,venter     A       AG      100     PASS    AF=0.774334;HP=1;DB;NF=21;NR=27;NS=51;CA=0;DP=306
		#1       1026158 .       T       TGGGGG  100     PASS    AF=0.115637;HP=1;NF=5;NR=2;NS=52;CA=0;DP=591
                
                #example of VCF4 generated by SamTools mpileup (Note that GT was not the first field in the FORMAT string)
                ##fileformat=VCFv4.0
		##samtoolsVersion=0.1.16 (r963:234)
		##fileformat=VCFv4.0
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  1247MFL0003.NOVO.srt.bam
		#chr1    14574   .       A       G       3.54    .       "DP=3;AF1=0.4999;CI95=0.5,0.5;DP4=1,0,2,0;MQ=21;PV4=1,1,1,1"    PL:GT:GQ        "31,0,34:0/1:32"
		#chr1    14930   .       A       G       37      .       "DP=19;AF1=0.5;CI95=0.5,0.5;DP4=7,5,5,1;MQ=25;PV4=0.6,6.3e-05,1,0.23"   PL:GT:GQ        "67,0,103:0/1:70"
		#chr1    16495   .       G       C       28      .       "DP=4;AF1=0.5;CI95=0.5,0.5;DP4=0,0,4,0;MQ=32"   PL:GT:GQ        "70,12,70:0/1:58"
		#chr1    59040   .       T       C       4.77    .       "DP=4;AF1=0.4999;CI95=0.5,0.5;DP4=0,2,2,0;MQ=22;PV4=0.33,0.21,1,1"      PL:GT:GQ        "33,0,39:0/1:35"
		#chr1    69270   .       A       G       46      .       "DP=20;AF1=0.5;CI95=0.5,0.5;DP4=2,0,18,0;MQ=24;PV4=1,1,1,0.28"  PL:GT:GQ        "94,18,100:0/1:78"
		#chr1    69511   .       A       G       24      .       "DP=5;AF1=0.5;CI95=0.5,0.5;DP4=1,0,2,1;MQ=25;PV4=1,0.46,1,0.039"        PL:GT:GQ        "54,0,57:0/1:55"
                
		#reserved VCF4 sub-fields in the INFO field
		#    * AA ancestral allele
		#    * AC allele count in genotypes, for each ALT allele, in the same order as listed
		#    * AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
		#    * AN total number of alleles in called genotypes
		#    * BQ RMS base quality at this position
		#    * CIGAR cigar string describing how to align an alternate allele to the reference allele
		#    * DB dbSNP membership
		#    * DP combined depth across samples, e.g. DP=154
		#    * END end position of the variant described in this record (esp. for CNVs)
		#    * H2 membership in hapmap2
		#    * MQ RMS mapping quality, e.g. MQ=52
		#    * MQ0 Number of MAPQ == 0 reads covering this record
		#    * NS Number of samples with data
		#    * SB strand bias at this position
		#    * SOMATIC indicates that the record is a somatic mutation, for cancer genomics
		#    * VALIDATED validated by follow-up experiment


		#SAMtools/BCFtools specific information
		#SAMtools/BCFtools may write the following tags in the INFO field in VCF/BCF.
		#Tag	Description
		#I16	16 integers:
		#1	#reference Q13 bases on the forward strand 	2	#reference Q13 bases on the reverse strand
		#3	#non-ref Q13 bases on the forward strand 	4	#non-ref Q13 bases on the reverse strand
		#5	sum of reference base qualities 	6	sum of squares of reference base qualities
		#7	sum of non-ref base qualities 	8	sum of squares of non-ref base qualities
		#9	sum of ref mapping qualities 	10	sum of squares of ref mapping qualities
		#11	sum of non-ref mapping qualities 	12	sum of squares of non-ref mapping qualities
		#13	sum of tail distance for ref bases 	14	sum of squares of tail distance for ref bases
		#15	sum of tail distance for non-ref bases 	16	sum of squares of tail distance for non-ref
		#INDEL	Indicating the variant is an INDEL.
		#DP	The number of reads covering or bridging POS.
		#DP4	Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling. Sum can be smaller than DP because low-quality bases are not counted.
		#PV4	P-values for 1) strand bias (exact test); 2) baseQ bias (t-test); 3) mapQ bias (t); 4) tail distance bias (t)
		#FQ	Consensus quality. If positive, FQ equals the phred-scaled probability of there being two or more different alleles. If negative, FQ equals the minus phred-scaled probability of all chromosomes being identical. Notably, given one sample, FQ is positive at hets and negative at homs.
		#AF1	EM estimate of the site allele frequency of the strongest non-reference allele.
		#CI95	Equal-tail (Bayesian) credible interval of the site allele frequency at the 95% level.
		#PC2	Phred-scaled probability of the alternate allele frequency of group1 samples being larger (,smaller) than of group2 samples.
		#PCHI2	Posterior weighted chi^2 P-value between group1 and group2 samples. This P-value is conservative.
		#QCHI2	Phred-scaled PCHI2
		#RP	Number of permutations yeilding a smaller PCHI2

		#example of triallelic variants generated by mpileup/bcftools
		#1       156706559       .       A       C,G     114     .	DP=20;AF1=1;CI95=1,1;DP4=0,0,1,19;MQ=60;FQ=-63  GT:PL:GQ	1/2:237,126,90,162,0,138:99
		#6       31129642        .       A       G,C     76      .	DP=31;AF1=1;CI95=1,1;DP4=0,0,28,3;MQ=60;FQ=-75  GT:PL:GQ	1/2:255,194,146,164,0,119:99
		#1       11297762        .       T       C,A     98      .	DP=19;AF1=1;CI95=1,1;DP4=0,0,17,1;MQ=60;FQ=-78  GT:PL:GQ	1/1:131,51,0,120,28,117:99
		
		my @field=split(/\t/,$_);
		@field >=8 or die "Error: invalid record found in VCF4 file (at least 8 tab-delimited fields expected): <$_>\n";
		my ($chr, $start, $ID, $ref_allele, $mut_allele, $quality_score, $filter, $info, $format, $sample) = @field;
		my ($end);
		my ($mut_allele2, $zygosity);
		
		if ($filterword) {		#ensure that the filter field contains the filterword
			$filter =~ m/\b$filterword\b/i or next;
		}
		
		$info =~ s/^"//; $info =~ s/"$//;
		
		#sometimes the alleles are not in the same case
		#chr1    1869771 1869774 actc    aCTctc          43.5    13      INDEL;DP=13;AF1=0.5;CI95=0.5,0.5;DP4=0,4,4,0;MQ=37;PV4=0.029,0.45,1,0.46
		$ref_allele = uc $ref_allele;
		$mut_allele = uc $mut_allele;
		
		#if ($ID eq '.' || $ID =~ /^rs/) {		#per MISHIMA, Hiroyuki suggestion (vcf4's third column (ID column) are not always ".")
		#	$end = $start;				#this block is commented out on 2011feb19
		#}
		
		if ($mut_allele eq '.') {			#no variant call was made at this position
			next;
		}
		
		if ($mut_allele =~ m/([^,]+),([\w,]+)/) {	#there could be more than two alternative alleles
			$mut_allele = $1;
			$mut_allele2 = $2;
		}
		
		if(length($ref_allele)==1 && length($mut_allele)==1) {  	### output snv
			if ($ref_allele =~ m/[^ACGTacgt]/ ) {
				print STDERR "WARNING: invalid allele record found in VCF4 file (ACGT expected): <$ref_allele> and <$mut_allele> in line <$_>\n";
				$ref_allele = 0;
			}
			if ( $mut_allele =~ m/[^ACGTacgt]/) {
				print STDERR "WARNING: invalid allele record found in VCF4 file (ACGT expected): <$ref_allele> and <$mut_allele> in line <$_>\n";
				$mut_allele = 0;
			}
				
			my ($unfiltered_read_depth) = $info =~ /DP=(\d+)/;
			my ($MappingQuality) = $info =~ /MQ=([^;]+)/; 
			my ($QualityByDepth) = $info =~ /QD=([^;]+)/;		
			

			
			if ($coverage) {
				defined $unfiltered_read_depth and $unfiltered_read_depth >= $coverage || next;
				if ($maxcoverage) {
					defined $unfiltered_read_depth and $unfiltered_read_depth <= $maxcoverage || next;
				}
			}
			
			if ($snpqual) {
				defined $QualityByDepth and $QualityByDepth >= $snpqual || next;		#the QD was used here as quality score
			}			
			
			if (defined $format) {
				my @format = split (/:/, $format);
				undef $gtpos;
				for my $i (0 .. @format-1) {
					if ($format[$i] eq 'GT') {
						$gtpos = $i;
						last;
					}
				}
				if (defined $sample and defined $gtpos) {
					my @sample = split (/:/, $sample);
					if ($sample[$gtpos] =~ m#^0/1# or $sample[$gtpos] =~ m#^1/0#) {
						$zygosity = 'het';
						$counthet++;
					} elsif ($sample[$gtpos] =~ m#^1/1#) {
						$zygosity = 'hom';
						$counthom++;
					} else {
						$zygosity = 'unknown';
						$countunknown++;
					}
				} else {		#sometimes the input VCF file does not contain the GT field!!!
					$zygosity = 'unknown';
					$countunknown++;
				}
			} else {
				$zygosity = 'unknown';
				$countunknown++;
			}

			#the subject is called as homozygous for the first alternative allele (genotype 1/1. i.e. C/C), but since there was one read containing A, samtools still keep both alleles in the VCF file (but gives a very low probabilities for it).
			#1       11297762        .       T       C,A     98      . DP=19;AF1=1;CI95=1,1;DP4=0,0,17,1;MQ=60;FQ=-78  GT:PL:GQ 1/1:131,51,0,120,28,117:99			
			if ($mut_allele2 and $zygosity eq 'hom') {
				$mut_allele2 = '';
			}

			if (not $mut_allele2) {
				if ($ref_allele eq 'A' and $mut_allele eq 'G' or $ref_allele eq 'G' and $mut_allele eq 'A' or $ref_allele eq 'C' and $mut_allele eq 'T' or $ref_allele eq 'T' and $mut_allele eq 'C') {
					$countti++;
					
				} else {
					$counttv++;
				}
			}
			
			#print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', $includeinfo ? "\t$otherinfo" : '', "\n";	#commented Sep 2011
			if ($includeinfo) {
				print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, $withzyg?"\t$zygosity":"", "\t", $otherinfo, "\n";
			} else {
				print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', "\n";
			}
			
			if ($allallele) {
				if ($mut_allele2) {
					my @mut_allele2 = split (/,/, $mut_allele2);
					for my $i (0 .. @mut_allele2-1) {
						if ($includeinfo) {
							print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele2[$i], $withzyg?"\t$zygosity":"", "\t", $otherinfo, "\n";
						} else {
							print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele2[$i], "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', "\n";
						}
					}
				}
			}
			
			$countsnp++;
		} elsif (length($ref_allele) > 1 || length($mut_allele) > 1) {  ### output indel
			my ($indel_read_depth1, $indel_read_depth2) = $info =~ /\bAC=([^,;]+),([^,;]+)/;		#number of reads supporting consensus indel, any indel
			my ($unfiltered_read_depth) = $info =~ /\bDP=(\d+)/;
			my ($MappingQuality) = $info =~ /\bMQ=([^;]+)/;
			my ($QualityByDepth) = $info =~ /\bQD=([^;]+)/;	
					
			if ($coverage) {
				defined $unfiltered_read_depth and $unfiltered_read_depth >= $coverage || next;
				if ($maxcoverage) {
					defined $unfiltered_read_depth and $unfiltered_read_depth <= $maxcoverage || next;
				}
			}
			
			if ($snpqual) {
				defined $QualityByDepth and $QualityByDepth >= $snpqual || next;		#the QD was used here as quality score
			}
			
			if (defined $indel_read_depth1 and $unfiltered_read_depth) {
				$indel_read_depth1/$unfiltered_read_depth >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				if ($indel_read_depth2) {
					$indel_read_depth1/$indel_read_depth2 >= $confraction or next;
				}
			}
			

			#example VCF4 records below:
			#20      2       .       TCG     T       .       PASS    DP=100
			#Chr1    5473    .       AT      ATT     23.5    .       INDEL;DP=16;AF1=0.5;CI95=0.5,0.5;DP4=4,2,3,1;MQ=42;PV4=1,0.41,0.042,0.24
			#Chr1    6498    .       ATTTT   ATTTTT  53.5    .       INDEL;DP=9;AF1=1;CI95=1,1;DP4=0,0,5,3;MQ=28
			
			if(length($ref_allele) > length ($mut_allele)) { 		# deletion or block substitution
				my $head = substr($ref_allele, 0, length ($mut_allele));
				if ($head eq $mut_allele) {
					print $chr,"\t";
					print $start+length($head),"\t";
					print $start+length($ref_allele)-1,"\t";
					
					$ref_allele = substr ($ref_allele, length ($mut_allele));
					print $ref_allele,"\t";
					print "-";
				} else {
					print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele;
				}
			} elsif(length($mut_allele) >= length ($ref_allele)) { 		# insertion or block substitution
				my $head = substr ($mut_allele, 0, length ($ref_allele));
				if ($head eq $ref_allele) {
					print $chr,"\t";	
					print $start+length($ref_allele)-1,"\t";
					print $start+length($ref_allele)-1,"\t";
					
					$mut_allele = substr ($mut_allele, length ($ref_allele));
					print "-\t";
					print $mut_allele;
				} else {
					print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele;
				}
			}
			

			if (defined $format) {
				my @format = split (/:/, $format);
				undef $gtpos;
				for my $i (0 .. @format-1) {
					if ($format[$i] eq 'GT') {
						$gtpos = $i;
						last;
					}
				}
				if (defined $sample and defined $gtpos) {
					my @sample = split (/:/, $sample);
					if ($sample[$gtpos] =~ m#^0/1# or $sample[$gtpos] =~ m#^1/0#) {
						$zygosity = 'het';
						$counthet++;
					} elsif ($sample[$gtpos] =~ m#^1/1#) {
						$zygosity = 'hom';
						$counthom++;
					} else {
						$zygosity = 'unknown';
						$countunknown++;
					}
				}
			} else {
				$zygosity = 'unknown';
				$countunknown++;
			}

			#commented out on 2011May21
			#if (defined $sample) {
			#	if ($sample =~ m#^0/1# or $sample =~ m#^1/0#) {
			#		print "\thet";
			#		$counthet++;
			#	} elsif ($sample =~ m#^1/1#) {
			#		print "\thom";
			#		$counthom++;
			#	} # BEGIN ARQ
			#	elsif ($sample =~ m#^./.#) {
			#		print "\tunknown";
			#		$countunknown++;
			#	} # END ARQ
			#}
			
			
			if ($includeinfo) {
				 print $withzyg?"\t$zygosity":"", "\t", $otherinfo;
			} else {
				print "\t$zygosity";
				defined $quality_score and print "\t", $quality_score;
				defined $unfiltered_read_depth and print "\t", $unfiltered_read_depth;
				
				#defined $indel_read_depth1 and print "\t", $indel_read_depth1;		#commented out Nov 2011
				defined $MappingQuality and print "\t", $MappingQuality;
				defined $QualityByDepth and print "\t", $QualityByDepth;		#added in Nov 2011 to be consistent with SNP output
				#$includeinfo and print "\t", $otherinfo;	#commented Sep 2011
			}
			print "\n";
			$countindel++;


			#do the same thing again, exactly like above, except that we work on second mutation;
			#in the future, consider rewrite this paragraph to make the code more elegant	
			if ($allallele and $mut_allele2) {
				my @mut_allele2 = split (/,/, $mut_allele2);
				for my $mut_allele2 (@mut_allele2) {
					if(length($ref_allele) > length ($mut_allele2)) { 		# deletion or block substitution
						my $head = substr($ref_allele, 0, length ($mut_allele2));
						if ($head eq $mut_allele2) {
							print $chr,"\t";
							print $start+length($head),"\t";
							print $start+length($ref_allele)-1,"\t";
							
							$ref_allele = substr ($ref_allele, length ($mut_allele2));
							print $ref_allele,"\t";
							print "-";
						} else {
							print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
						}
					} elsif(length($mut_allele2) > length ($ref_allele)) { 		# insertion or block substitution
						my $head = substr ($mut_allele2, 0, length ($ref_allele));
						if ($head eq $ref_allele) {
							print $chr,"\t";	
							print $start+length($ref_allele)-1,"\t";
							print $start+length($ref_allele)-1,"\t";
							
							$mut_allele = substr ($mut_allele2, length ($ref_allele));
							print "-\t";
							print $mut_allele2;
						} else {
							print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
						}
					} else {		#identical length of alleles
						print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
					}

					if (defined $sample) {
						if ($sample =~ m#^0/1# or $sample =~ m#^1/0#) {
							$zygosity = "het";
							$counthet++;
						} elsif ($sample =~ m#^1/1#) {
							$zygosity =  "hom";
							$counthom++;
						} # BEGIN ARQ
						elsif ($sample =~ m#^./.#) {
							$zygosity = "unknown";
							$countunknown++;
						} # END ARQ
					}
										
					if ($includeinfo) {
						print $withzyg?"\t$zygosity":"", "\t", $otherinfo;
					} else {
						print "\t", $zygosity;
						print "\t", $quality_score;
						defined $unfiltered_read_depth and print "\t", $unfiltered_read_depth;
						
						defined $indel_read_depth1 and print "\t", $indel_read_depth1;
						defined $MappingQuality and print "\t", $MappingQuality;
						#$includeinfo and print "\t", $otherinfo;
					}
					print "\n";

				}
			}
		} # CLOSE BY RAJESH IF FOR CHROMOSOME 
		}
		$countvar++;
	}
	my $triallelic = $countsnp-$countti-$counttv;
	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
	print STDERR "NOTICE: Among ${\($counthet+$counthom+$countunknown)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
}


=head1 SYNOPSIS

 convert2annovar.pl [arguments] <variantfile>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --format <string>		input format (default: pileup)
            --outfile <file>		output file name (default: STDOUT)
            --snpqual <float>		quality score threshold in pileup file (default: 20)
            --snppvalue <float>		SNP P-value threshold in GFF3-SOLiD file (default: 1)
            --coverage <int>		read coverage threshold in pileup file (default: 0)
            --maxcoverage <int>		maximum coverage threshold (default: none)
            --includeinfo		include supporting information in output
            --chr <string>		specify the chromosome (for CASAVA format)
            --chrmt <string>		chr identifier for mitochondria (default: M)
            --altcov <int>		alternative allele coverage threshold (for pileup format)
            --allelicfrac		print out allelic fraction rather than het/hom status (for pileup format)
            --fraction <float>		minimum allelic fraction to claim a mutation (for pileup/vcf4_indel format)
            --filter <string>		output variants with this filter (case insensitive, for vcf4 format)
            --confraction <float>	minimum consensus indel / all indel fraction (for vcf4 format)
            --allallele			print all alleles when multiple calls are present (for vcf4 format)
            --withzyg			print zygosity when -includeinfo is used (for vcf4 format)

 Function: convert variant call file generated from various software programs 
 into ANNOVAR input format
 
 Example: convert2annovar.pl -format pileup -outfile variant.query variant.pileup
          convert2annovar.pl -format cg -outfile variant.query variant.cg
          convert2annovar.pl -format soap variant.snp > variant.avinput
          convert2annovar.pl -format maq variant.snp > variant.avinput
          convert2annovar.pl -format casava -chr 1 variant.snp > variant.avinput
          convert2annovar.pl -format vcf4 variantfile > variant.avinput
          convert2annovar.pl -format vcf4 -filter pass variantfile > variant.avinput

 Version: $LastChangedDate: 2012-05-17 10:54:12 -0700 (Thu, 17 May 2012) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--format>

the format of the input files.

=item B<--outfile>

specify the output file name. By default, output is written to STDOUT.

=item B<--snpqual>

quality score threshold in the pileup file, such that variant calls with lower 
quality scores will not be printed out in the output file. When VCF4 file is 
used, this argument works on the Quality-by-Depth measure, rather than the raw 
quality measure.

=item B<--coverage>

read coverage threshold in the pileup file, such that variants calls generated 
with lower coverage will not be printed in the output file.

=item B<--includeinfo>

specify that the output should contain additional information in the input line. 
By default, only the chr, start, end, reference allele, observed allele and 
homozygosity status are included in output files.

=item B<--chr>

specify the chromosome for CASAVA format

=item B<--chrmt>

specify the name of mitochondria chromosome (default is MT)

=item B<--altcov>

the minimum coverage of the alternative (mutated) allele to be printed out in 
output

=item B<--fraction>

specify the minimum fraction of alternative allele, to print out the mutation. 
For example, a site has 10 reads, 3 supports alternative allele. A -fraction of 
0.4 will not allow the mutation to be printed out.

=item B<--species>

specify the species from which the sequencing data is obtained. For the GFF3-
SOLiD format, when species is human, the chromosome 23, 24 and 25 will be 
converted to X, Y and M, respectively.

=item B<--filter>

for VCF4 file, only print out variant calls with this filter annotated. For 
example, if using GATK VariantFiltration walker, you will see PASS, 
GATKStandard, HARD_TO_VALIDATE, etc in the filter field. Using 'pass' as a 
filter is recommended in this case.

=item B<--confraction>

consesus indel fraction, calculated as reads supporting consensus indels divided 
by reads supporting any indels

=item B<--allallele>

print all alleles for mutations at a locus, rather than the first allele, if the 
input VCF4 file contains multiple alternative alleles for a mutation. By 
default, this option is off. When it is on, two lines will be printed out in the 
output, and both will have the same quality scores as VCF4 does not provide 
separate quality scores for individual alleles.

=back

=head1 DESCRIPTION

This program is used to convert variant call file generated from various 
software programs into ANNOVAR input format. Currently, the program can handle 
Samtools genotype-calling pileup format, Solid GFF format, Complete Genomics 
variant format, SOAP format. These formats are described below.

=over 8

=item * B<pileup format>

The pileup format can be produced by the Samtools genotyping calling subroutine. 
Note that the phrase pileup format can be used in several instances, and here I 
am only referring to the pileup files that contains the actual genotype calls. 

Using SamTools, given an alignment file in BAM format, a pileup file with 
genotype calls can be produced by the command below:

	samtools pileup -vcf ref.fa aln.bam> raw.pileup
	samtools.pl varFilter raw.pileup > final.pileup

ANNOVAR will automatically filter the pileup file so that only SNPs reaching a 
quality threshold are printed out (default is 20, use --snpqual argument to 
change this). Most likely, users may want to also apply a coverage threshold, 
such that SNPs calls from only a few reads are not considered. This can be 
achieved using the -coverage argument (default value is 0).

An example of pileup files for SNPs is shown below:

	chr1 556674 G G 54 0 60 16 a,.....,...,.... (B%A+%7B;0;%=B<:
	chr1 556675 C C 55 0 60 16 ,,..A..,...,.... CB%%5%,A/+,%....
	chr1 556676 C C 59 0 60 16 g,.....,...,.... .B%%.%.?.=/%...1
	chr1 556677 G G 75 0 60 16 ,$,.....,...,.... .B%%9%5A6?)%;?:<
	chr1 556678 G K 60 60 60 24 ,$.....,...,....^~t^~t^~t^~t^~t^~t^~t^~t^~t B%%B%<A;AA%??<=??;BA%B89
	chr1 556679 C C 61 0 60 23 .....a...a....,,,,,,,,, %%1%&?*:2%*&)(89/1A@B@@
	chr1 556680 G K 88 93 60 23 ..A..,..A,....ttttttttt %%)%7B:B0%55:7=>>A@B?B;
	chr1 556681 C C 102 0 60 25 .$....,...,....,,,,,,,,,^~,^~. %%3%.B*4.%.34.6./B=?@@>5.
	chr1 556682 A A 70 0 60 24 ...C,...,....,,,,,,,,,,. %:%(B:A4%7A?;A><<999=<<
	chr1 556683 G G 99 0 60 24 ....,...,....,,,,,,,,,,. %A%3B@%?%C?AB@BB/./-1A7?

The columns are chromosome, 1-based coordinate, reference base, consensus base, 
consensus quality, SNP quality, maximum mapping quality of the reads covering 
the sites, the number of reads covering the site, read bases and base qualities.

An example of pileup files for indels is shown below:

	seq2  156 *  +AG/+AG  71  252  99  11  +AG  *  3  8  0

ANNOVAR automatically recognizes both SNPs and indels in pileup file, and process them correctly.

An example of the Complete Genomics format is shown below:

	#BUILD  1.5.0.5
	#GENERATED_AT   2009-Nov-03 19:52:21.722927
	#GENERATED_BY   dbsnptool
	#TYPE   VAR-ANNOTATION
	#VAR_ANN_SET    /Proj/Pipeline/Production_Data/REF/HUMAN-F_06-REF/dbSNP.csv
	#VAR_ANN_TYPE   dbSNP
	#VERSION        0.3
	
	>locus  ploidy  haplotype       chromosome      begin   end     varType reference       alleleSeq       totalScore      hapLink xRef
	1       2       all     chr1    0       959     no-call =       ?                       
	2       2       all     chr1    959     972     =       =       =                       
	3       2       all     chr1    972     1001    no-call =       ?                       
	4       2       all     chr1    1001    1008    =       =       =                       
	5       2       all     chr1    1008    1114    no-call =       ?                       
	6       2       all     chr1    1114    1125    =       =       =                       
	7       2       all     chr1    1125    1191    no-call =       ?                       
	8       2       all     chr1    1191    1225    =       =       =                       
	9       2       all     chr1    1225    1258    no-call =       ?                       
	10      2       all     chr1    1258    1267    =       =       =                       
	12      2       all     chr1    1267    1275    no-call =       ?                       
	13      2       all     chr1    1275    1316    =       =       =                       
	14      2       all     chr1    1316    1346    no-call =       ?                       
	15      2       all     chr1    1346    1367    =       =       =                       
	16      2       all     chr1    1367    1374    no-call =       ?                       
	17      2       all     chr1    1374    1388    =       =       =                       
	18      2       all     chr1    1388    1431    no-call =       ?                       
	19      2       all     chr1    1431    1447    =       =       =                       
	20      2       all     chr1    1447    1454    no-call =       ?                       

The following information is provided in documentation from Complete Genomics, that describes the var-ASM format.

	1. locus. Identifier of a particular genomic locus
	2. ploidy. The ploidy of the reference genome at the locus (= 2 for autosomes, 2 for pseudoautosomal regions on the sex chromosomes, 1 for males on the non-pseudoautosomal parts of the sex chromosomes, 1 for mitochondrion, '?' if varType is 'no-ref' or 'PAR-called-in-X'). The reported ploidy is fully determined by gender, chromosome and location, and is not inferred from the sequence data.
	3. haplotype. Identifier for each haplotype at the variation locus. For diploid genomes, 1 or 2. Shorthand of 'all' is allowed where the varType field is one of 'ref', 'no-call', 'no-ref', or 'PAR-called-in-X'. Haplotype numbering does not imply phasing; haplotype 1 in locus 1 is not necessarily in phase with haplotype 1 in locus 2. See hapLink, below, for phasing information.
	4. chromosome. Chromosome name in text: 'chr1','chr2', ... ,'chr22','chrX','chrY'. The mitochondrion is represented as 'chrM'. The pseudoautosomal regions within the sex chromosomes X and Y are reported at their coordinates on chromosome X.
	5. begin. Reference coordinate specifying the start of the variation (not the locus) using the half-open zero-based coordinate system. See section 'Sequence Coordinate System' for more information.
	6. end. Reference coordinate specifying the end of the variation (not the locus) using the half-open zero-based coordinate system. See section 'Sequence Coordinate System' for more information.
	7. varType. Type of variation, currently one of:
		snp: single-nucleotide polymorphism
		ins: insertion
		del: deletion
		sub: Substitution of one or more reference bases with the bases in the allele column
		'ref' : no variation; the sequence is identical to the reference sequence on the indicated haplotype
		no-call-rc: 'no-call reference consistent 'one or more bases are ambiguous, but the allele is potentially consistent with the reference
		no-call-ri: 'no-call reference inconsistent' one or more bases are ambiguous, but the allele is definitely inconsistent with the reference
		no-call: an allele is completely indeterminate in length and composition, i.e. alleleSeq = '?'
		no-ref: the reference sequence is unspecified at this locus.
		PAR-called-in-X: this locus overlaps one of the pseudoautosomal regions on the sex chromosomes. The called sequence is reported as diploid sequence on Chromosome X; on chromosome Y the sequence is reported as varType = 'PAR-called-in-X'.
	8. reference. The reference sequence for the locus of variation. Empty when varType is ins. A value of '=' indicates that the user must consult the reference for the sequence; this shorthand is only used in regions where no haplotype deviates from the reference sequence.
	9. alleleSeq. The observed sequence at the locus of variation. Empty when varType is del. '?' isused to indicate 0 or more unknown bases within the sequence; 'N' is used to indicate exactly one unknown base within the sequence.'=' is used as shorthand to indicate identity to the reference sequence for non-variant sequence, i.e. when varType is 'ref'.
	10. totalScore. A score corresponding to a single variation and haplotype, representing the confidence in the call.
	11. hapLink. Identifier that links a haplotype at one locus to haplotypes at other loci. Currently only populated for very proximate variations that were assembled together. Two calls that share a hapLink identifier are expected to be on the same haplotype,
	12. xRef. Field containing external variation identifiers, currently only populated for variations corroborated directly by dbSNP. Format: dbsnp:[rsID], with multiple entries separated by the semicolon (;).

In older versions of the format specification, the sub keyword used to be insdel 
keyword. ANNOVAR takes care of this.

=item * B<SOAPsnp format>

An example of the SOAP SNP caller format is shown below:

	chr8  35782  A  R  1  A  27  1  2  G  26  1  2  5   0.500000  2.00000  1  5   
	chr8  35787  G  R  0  G  25  4  6  A  17  2  4  10  0.266667  1.60000  0  5   

The following information is provided in documentation from BGI who developed 
SOAP suite. It differs slightly from the description at the SOAPsnp website, and 
presumably the website is outdated.

	Format description:(left to right)
	1. Chromosome name
	2. Position of locus
	3. Nucleotide at corresponding locus of reference sequence
	4. Genotype of sequencing sample
	5. Quality value
	6. nucleotide with the highest probability(first nucleotide)
	7. Quality value of the nucleotide with the highest probability
	8. Number of supported reads that can only be aligned to this locus 
	9. Number of all supported reads that can be aligned to this locus
	10. Nucleotide with higher probability 
	11. Quality value of nucleotide with higher probability 
	12. Number of supported reads that can only be aligned to this locus 
	13. Number of all supported reads that can be aligned to this locus 
	14. Total number of reads that can be aligned to this locus 
	15. Order and quality value
	16. Estimated copy number for this locus 
	17. Presence of this locus in the dbSNP database. 1 refers to presence and 0 refers to inexistence
	18. The distance between this locus and another closest SNP

=item * B<SOAPindel format>

The current version of ANNOVAR handles SoapSNP and SoapIndel automatically via a 
single argument '--format soap'. An example of SOAP indel caller format is shown 
below:

	chr11   44061282        -       +2      CT      Hete
	chr11   45901572        +       +1      C       Hete
	chr11   48242562        *       -3      TTC     Homo
	chr11   57228723        *       +4      CTTT    Homo
	chr11   57228734        *       +4      CTTT    Homo
	chr11   57555685        *       -1      C       Hete
	chr11   61482191        -       +3      TCC     Hete
	chr11   64608031        *       -1      T       Homo
	chr11   64654936        *       +1      C       Homo
	chr11   71188303        +       -1      T       Hete
	chr11   75741034        +       +1      T       Hete
	chr11   76632438        *       +1      A       Hete
	chr11   89578266        *       -2      AG      Homo
	chr11   104383261       *       +1      T       Hete
	chr11   124125940       +       +4      CCCC    Hete
	chr12   7760052 *       +1      T       Homo
	chr12   8266049 *       +3      ACG     Homo

I do not see a documentation describing this format yet as of September 2010.

=item B<--SOAPsv format>

An example is given below:

	Chr2 Deletion 42894 43832 43167 43555 388 0-0-0 FR 41

An explanation of the structural variation format is given below:

	Format description (from left to right)
	1. Chromosome name
	2. Type of structure variation
	3. Minimal value of start position in cluster
	4. Maximal value of end position in cluster
	5. Estimated start position of this structure variation
	6. Estimated end position of this structure variation
	7. Length of SV
	8. Breakpoint of SV (only for insertion)
	9. Unusual matching mode (F refers to align with forward sequence, R refers
	to align with reverse
	sequence)
	10. number of paired-end read which support this structure variation

=item * B<MAQ format>

MAQ can perform alignment and generate genotype calls, including SNP calls and 
indel calls. The format is described below:

For indel header: The output is TAB delimited with each line consisting of chromosome, start 
position, type of the indel, number of reads across the indel, size of the indel 
and inserted/deleted nucleotides (separated by colon), number of indels on the 
reverse strand, number of indels on the forward strand, 5' sequence ahead of the 
indel, 3' sequence following the indel, number of reads aligned without indels 
and three additional columns for filters.

An example is below:

	chr10   110583  -       2       -2:AG   0       1       GCGAGACTCAGTATCAAAAAAAAAAAAAAAAA        AGAAAGAAAGAAAAAGAAAAAAATAGAAAGAA        1       @2,     @72,   @0,
	chr10   120134  -       8       -2:CA   0       1       CTCTTGCCCGCTCACACATGTACACACACGCG        CACACACACACACACACATCAGCTACCTACCT        7       @65,62,61,61,45,22,7,   @9,12,13,13,29,52,67,   @0,0,0,0,0,0,0,
	chr10   129630  -       1       -1:T    1       0       ATGTTGTGACTCTTAATGGATAAGTTCAGTCA        TTTTTTTTTAGCTTTTAACCGGACAAAAAAAG        0       @       @      @
	chr10   150209  -       1       4:TTCC  1       0       GCATATAGGGATGGGCACTTTACCTTTCTTTT        TTCCTTCCTTCCTTCCTTCCCTTTCCTTTCCT        0       @       @      @
	chr10   150244  -       2       -4:TTCT 0       1       CTTCCTTCCTTCCTTCCCTTTCCTTTCCTTTC        TTCTTTCTTTCTTTCTTTCTTTTTTTTTTTTT        0       @       @      @
	chr10   159622  -       1       3:AGG   0       1       GAAGGAGGAAGGACGGAAGGAGGAAGGAAGGA        AGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGA        0       @       @      @
	chr10   206372  -       2       2:GT    1       0       ATAATAGTAACTGTGTATTTGATTATGTGTGC        GTGTGTGTGTGTGTGTGTGTGTGTGCGTGCTT        1       @37,    @37,   @8,
	chr10   245751  -       11      -1:C    0       1       CTCATAAATACAAGTCATAATGAAAGAAATTA        CCACCATTTTCTTATTTTCATTCATTTTTAGT        10      @69,64,53,41,30,25,22,14,5,4,   @5,10,21,33,44,49,52,60,69,70,  @0,0,0,0,0,0,0,0,0,0,
	chr10   253066  -       1       2:TT    0       1       TATTGATGAGGGTGGATTATACTTTAGAACAC        TATTCAAACAGTTCTTCCACATATCTCCCTTT        0       @       @      @
	chr10   253455  -       2       -3:AAA  1       0       GTTGCACTCCAGCCTGGCGAGATTCTGTCTCC        AAAAAAAAAAAAAAAAATTGTTGTGAAATACA        1       @55,    @19,   @4,

For snp output file: Each line consists of chromosome, position, reference base, 
consensus base, Phred-like consensus quality, read depth, the average number of 
hits of reads covering this position, the highest mapping quality of the reads 
covering the position, the minimum consensus quality in the 3bp flanking regions 
at each side of the site (6bp in total), the second best call, log likelihood 
ratio of the second best and the third best call, and the third best call.

An example is below:

	chr10   83603   C       T       28      12      2.81    63      34      Y       26      C
	chr10   83945   G       R       59      61      4.75    63      62      A       47      G
	chr10   83978   G       R       47      40      3.31    63      62      A       21      G
	chr10   84026   G       R       89      22      2.44    63      62      G       49      A
	chr10   84545   C       T       54      9       1.69    63      30      N       135     N
	chr10   85074   G       A       42      5       1.19    63      38      N       108     N
	chr10   85226   A       T       42      5       1.00    63      42      N       107     N
	chr10   85229   C       T       42      5       1.00    63      42      N       112     N
	chr10   87518   A       G       39      4       3.25    63      38      N       9       N
	chr10   116402  T       C       39      4       1.00    63      38      N       76      N


=item * B<CASAVA format>

An example of Illumina CASAVA format is given below:

	#position       A       C       G       T       modified_call   total   used    score           reference       type
	14930   3       0       8       0       GA      11      11      29.10:11.10             A       SNP_het2
	14933   4       0       7       0       GA      11      11      23.50:13.70             G       SNP_het1
	14976   3       0       8       0       GA      11      11      24.09:9.10              G       SNP_het1
	15118   2       1       4       0       GA      8       7       10.84:6.30              A       SNP_het2

An example of the indels is given below:

	# ** CASAVA depth-filtered indel calls **
	#$ CMDLINE /illumina/pipeline/install/CASAVA_v1.7.0/libexec/CASAVA-1.7.0/filterIndelCalls.pl--meanReadDepth=2.60395068970547 --indelsCovCutoff=-1 --chrom=chr1.fa /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0000.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0001.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0002.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0003.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0004.txt
	#$ CHROMOSOME chr1.fa
	#$ MAX_DEPTH undefined
	#
	#$ COLUMNS pos CIGAR ref_upstream ref/indel ref_downstream Q(indel) max_gtype Q(max_gtype) max2_gtype bp1_reads ref_reads indel_reads other_reads repeat_unit ref_repeat_count indel_repeat_count
	948847  1I      CCTCAGGCTT      -/A     ATAATAGGGC      969     hom     47      het     22      0       16      6       A       1       2
	978604  2D      CACTGAGCCC      CT/--   GTGTCCTTCC      251     hom     20      het     8       0       4       4       CT      1       0
	1276974 4I      CCTCATGCAG      ----/ACAC       ACACATGCAC      838     hom     39      het     18      0       14      4       AC      2       4
	1289368 2D      AGCCCGGGAC      TG/--   GGAGCCGCGC      1376    hom     83      het     33      0       25      9       TG      1       0

=item * B<VCF4 format>

VCF4 can be used to describe both population-level variation information, or for 
reads derived from a single individual.

One example of the indel format for one individual is given below:

	##fileformat=VCFv4.0
	##IGv2_bam_file_used=MIAPACA2.alnReAln.bam
	##INFO=<ID=AC,Number=2,Type=Integer,Description="# of reads supporting consensus indel/any indel at the site">
	##INFO=<ID=DP,Number=1,Type=Integer,Description="total coverage at the site">
	##INFO=<ID=MM,Number=2,Type=Float,Description="average # of mismatches per consensus indel-supporting read/per reference-supporting read">
	##INFO=<ID=MQ,Number=2,Type=Float,Description="average mapping quality of consensus indel-supporting reads/reference-supporting reads">
	##INFO=<ID=NQSBQ,Number=2,Type=Float,Description="Within NQS window: average quality of bases from consensus indel-supporting reads/from reference-supporting reads">
	##INFO=<ID=NQSMM,Number=2,Type=Float,Description="Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads">
	##INFO=<ID=SC,Number=4,Type=Integer,Description="strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
	##IndelGenotyperV2=""
	##reference=hg18.fa
	##source=IndelGenotyperV2
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Miapaca_trimmed_sorted.bam      
	chr1    439     .       AC      A       .       PASS    AC=5,5;DP=7;MM=7.0,3.0;MQ=23.4,1.0;NQSBQ=23.98,25.5;NQSMM=0.04,0.0;SC=2,3,0,2   GT      1/0
	chr1    714048  .       T       TCAAC   .       PASS    AC=3,3;DP=9;MM=3.0,7.1666665;MQ=1.0,10.833333;NQSBQ=23.266666,21.932203;NQSMM=0.0,0.15254237;SC=3,0,3,3 GT      0/1
	chr1    714049  .       G       GC      .       PASS    AC=3,3;DP=9;MM=3.0,7.1666665;MQ=1.0,10.833333;NQSBQ=23.233334,21.83051;NQSMM=0.0,0.15254237;SC=3,0,3,3  GT      0/1
	chr1    813675  .       A       AATAG   .       PASS    AC=5,5;DP=8;MM=0.4,1.0;MQ=5.0,67.0;NQSBQ=25.74,25.166666;NQSMM=0.0,0.033333335;SC=4,1,1,2       GT      0/1
	chr1    813687  .       AGAGAGAGAGAAG   A       .       PASS    AC=5,5;DP=8;MM=0.4,1.0;MQ=5.0,67.0;NQSBQ=24.54,25.2;NQSMM=0.02,0.06666667;SC=4,1,1,2    GT      1/0


=back

The code was written by Dr. Kai Wang and modified by Dr. Germán Gastón Leparc. 
Various users have provided sample input files for many SNP callin software, for 
the development of conversion subroutines. We thank these users for their 
continued support to improve the functionality of the script.

For questions or comments, please contact kai@openbioinformatics.org.

=cut
