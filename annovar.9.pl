#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#Version Log
# 12/05/2014
# Started pusshing to github
#
# 12/04/2014
# Cleaned up the code. Removed most of the hardcoded paths.
# Archieved old scripts 	
# 
# 04/10/2014
# Added new TCGA annotation
#
# 08/27/2012
# Modified some of the components
# refer to README for information on input file type
# 
#
my $host =`hostname`;
chomp $host;
if($host =~ /helix/){
	print STDERR "\n\n\n\tThis pipeline can only be ran on BIOWULF\n\n\n";
	exit;
}
my $file;
my $cs = 5;
my $node= "norm";
my $code_dir = "/data/khanlab/apps/annovar";
my $data_dir = "/data/khanlab/ref/annovar/humandb";
GetOptions(
		'infile=s'	=>\$file,
		'cs=s'		=>\$cs,
		'q=s'		=>\$node,
		help		=>sub {pod2usage(1);},
	  )or pod2usage(2);

=head1 SYNOPSIS

annovar.pl[options]
Usage:
	--infile	give the file containig output 
	--cs		coordinate system
	--q		ccr if you want to use CCR nodes	
This script is run annovar http://www.openbioinformatics.org/annovar/ and add the result columns to the output file 
Usage:
	/data/khanlab/apps/annovar/annovar.pl -infile /data/khanlab/apps/annovar/working/testData -cs 5
	
-infile is a file in which the first 5 columns are following 
Chr	Start	End	Ref	Alt # User commnets(other columns)
The program will output header to the annovar result if the headings of the first columns are same as above.

-cs cooridnate system. there are 2 ways your input can be handled:
	Chr	Position	Ref	Alt # User comments(other columns)

More detailed documentation can be found in README file in the /data/khanlab/apps/annovar/


Contact Rajesh Patidar rajbtpatidar@gmail.com rajesh.patidar@nih.gov if you are having trouble using it:

=cut

my $DIR ='';
my $FILE='';
if (!$file){
	print STDERR "Give the file -infile <file containing snv list in residue or space based coordinate system>\n";
	die;	
}
if ($cs eq '5'){
	print STDERR "\n\nconsidered that input have Chr,Start,End,Ref,Alt columns \n";
	print STDERR "Your final output will be $file.annotated.txt\n";
}
elsif($cs eq '4'){
	print STDERR "4 = Residue based coordinate system i.e. Chr,Pos,Ref,Alt columns \n";
	print STDERR "Your final output will be $file.annotated.txt\n";
}
else{
	print STDERR "The default coordinate system is SPACE based -cs 5 but you can also give RESIDUE based\n";
	print STDERR "options 4 or 5 allowed only\n";
	die;
}

my $dir = `pwd`;
chomp $dir;


if(-e "$dir/$file"){
	print STDERR "\nYour input file is $dir/$file\n";
	$DIR = $dir;
	$FILE =$file;
}
elsif(-e $file){
	print STDERR "\nYour input file is $file\n";
	if ($file =~ /(.*)\/(.*)/){
		$DIR = $1;
		$FILE =$2;
	}
}
else{
	print STDERR "File $file is not present is your current directory\n";
	print STDERR "Terminating...\n";
	die;
}
#########################
## Make Input Files    ##
#########################
`cd $DIR; $code_dir/scripts/MakeAnnotationInputs.pl $FILE $cs`;
##########################
####     ANNOVAR    ######
##########################
unless (open (OUT, ">$DIR/annovar_$FILE.sh")){
	print "I tried writing in $DIR but failed.. \n\nPlease check the permissions  !!\n";
	die;
}
print OUT "#!/bin/sh\n";
print OUT "#\n";
print OUT "# This script will be submitted on cluster \n";
print OUT "#\n";
print OUT "#PBS -N annovar\n";
print OUT "#\n";
print OUT "cd $DIR\n";
print OUT "$code_dir/scripts/multi_db.pl --input $FILE.anno --output $FILE.annovar --autothread --dbtype gene --dbtype band --dbtype snp138 --dbtype 1000g2014oct_all --dbtype 1000g2014oct_eur --dbtype 1000g2014oct_afr --dbtype 1000g2014oct_amr --dbtype 1000g2014oct_eas --dbtype 1000g2014oct_sas --dbtype esp6500_all --dbtype esp6500_ea --dbtype esp6500_aa --dbtype cg69 --dbtype nci60 --dbtype cadd --dbtype clinvar_20140702 --dbtype cosmic70\n";
print OUT "rm -rf $FILE.anno annovar_$FILE.sh\n";
close OUT;
###########################
###	SIFT4.0.3b 	###
###########################
unless (open (OUT, ">$DIR/SIFT_$FILE.sh")){
	print "I tried writing in $DIR but failed.. \n\nPlease check the permissions  !!\n";
	die;
}
print OUT "#!/bin/sh\n";
print OUT "#\n";
print OUT "# This script will be submitted on cluster \n";
print OUT "#\n";
print OUT "#PBS -N ansift\n";
print OUT "#\n";
print OUT "cd $DIR\n";
print OUT "clearscratch\n";
print OUT "module load SIFT\n";
print OUT "SIFT_exome_nssnvs.pl -i $FILE.sift -d /fdb/sift/Human_db_37 -o \$SIFT_SCRATCHDIR -z $DIR/$FILE.sift_predictions.tsv\n";
print OUT "$code_dir/scripts/ParseSIFT.pl $DIR/$FILE.sift_predictions.tsv >$DIR/$FILE.sift.out\n";
print OUT "rm -rf $DIR/$FILE.sift_predictions.tsv $FILE.sift SIFT_$FILE.sh\n";
close OUT;
###########################
### 	PolyPhen2	###
###########################

unless (open (OUT, ">$DIR/PPH_$FILE.sh")){
	print "I tried writing in $DIR but failed.. \n\nPlease check the permissions  !!\n";
	die;
}
my $user = `whoami`;
chomp $user;
print OUT "#!/bin/sh\n";
print OUT "#\n";
print OUT "# This script will be submitted on cluster \n";
print OUT "#\n";
print OUT "#PBS -N anpph\n";
print OUT "#\n";
print OUT "cd $DIR\n";
print OUT "/usr/local/polyphen-2.2.2r405/bin/mapsnps.pl -g hg19 -U -y $FILE.pph2 $FILE.pph\n"; # changed to qsub by David Hoover
print OUT "/usr/local/polyphen-2.2.2r405/bin/pph_swarm.pl --block -d /spin1/scratch/$user/tmp -o $FILE.pph2.out $FILE.pph2 --queue $node\n";
print OUT "$code_dir/scripts/ParsePPH2.pl $FILE.pph2.out >$FILE.polyphen2.out\n";
print OUT "rm -rf $FILE.pph2 anpph.e* anpph.o* swb*.o* swb*.e* pph_roundup.e* pph_roundup.o* $FILE.pph2.out PPH_$FILE.sh\n";

close OUT;
###########################
### Join All the results ##
###########################
unless (open (OUT, ">$DIR/join_$FILE.sh")){
	print "I tried writing in $DIR but failed.. \n\nPlease check the permissions  !!\n";
	die;
}
$code_dir = "$code_dir/scripts";
#/data/khanlab/apps/annovar/scripts/addAnnotation.pl
print OUT "#!/bin/sh\n";
print OUT "#\n";
print OUT "# This script will be submitted on cluster \n";
print OUT "#\n";
print OUT "#PBS -N anjoin\n";
print OUT "#\n";
print OUT "cd $DIR\n";
print OUT "$code_dir/join.anno.sift.pph2.pl $FILE.annovar $FILE.sift.out $FILE.polyphen2.out >$FILE.a.c.s.p\n";
print OUT "$code_dir/addAnnotation.pl  $data_dir/hg19_hgmd.2014.3.txt $FILE.a.c.s.p >$FILE.a.c.s.p.h\n";
print OUT "$code_dir/GeneAnnotation.pl $data_dir/hg19_CancerGeneCensus.08.25.14.txt $FILE.a.c.s.p.h   CancerGeneCensus.08.25.14 >$FILE.a.c.s.p.h.c\n";
print OUT "$code_dir/GeneAnnotation.pl $data_dir/reportableGenes                    $FILE.a.c.s.p.h.c ACMG_reportableGenes      >$FILE.a.c.s.p.h.c.A\n";
print OUT "$code_dir/addAnnotation.pl $data_dir/hg19_MATCH_v1_08_2014.txt       $FILE.a.c.s.p.h.c.A >$FILE.a.c.s.p.h.c.A.m\n";
print OUT "$code_dir/addAnnotation.pl $data_dir/hg19_MATCH_v2_09_2014.txt       $FILE.a.c.s.p.h.c.A.m >$FILE.a.c.s.p.h.c.A.m.m\n";
print OUT "$code_dir/addAnnotation.pl $data_dir/hg19_MCG.08.25.14.txt           $FILE.a.c.s.p.h.c.A.m.m >$FILE.a.c.s.p.h.c.A.m.m.m\n";
print OUT "$code_dir/addAnnotation.pl $data_dir/PediatricGenome.11.24.14.txt    $FILE.a.c.s.p.h.c.A.m.m.m > $FILE.a.c.s.p.h.c.A.m.m.m.t\n";
print OUT "$code_dir/addAnnotation.pl $data_dir/hg19_Clinseq_951.txt            $FILE.a.c.s.p.h.c.A.m.m.m.t >$FILE.a.c.s.p.h.c.A.m.m.m.t.c\n";
print OUT "$code_dir/addAnnotation.pl $data_dir/hg19_exac02.txt                 $FILE.a.c.s.p.h.c.A.m.m.m.t.c >$FILE.a.c.s.p.h.c.A.m.m.m.t.c.e\n";
print OUT "$code_dir/addAnnotation.pl $data_dir/hg19_exac02.1.txt               $FILE.a.c.s.p.h.c.A.m.m.m.t.c.e >$FILE.a.c.s.p.h.c.A.m.m.m.t.c.e.e\n";
print OUT "$code_dir/add_WWW.pl                                                 $FILE.a.c.s.p.h.c.A.m.m.m.t.c.e.e >$FILE.a.c.s.p.h.c.A.m.m.m.t.c.e.e.l\n";
print OUT "$code_dir/joinAllAnnotations.pl $FILE.a.c.s.p.h.c.A.m.m.m.t.c.e.e.l $FILE $cs >$FILE.annotated.txt\n";
print OUT " rm -rf $FILE.a.c.* $FILE.annovar $FILE.sift.out $FILE.polyphen2.out $FILE.pph join_$FILE.sh anpph* ansift* annovar.e* annovar.o* \n";
print OUT "chmod 760 $FILE.annotated.txt\n";
print OUT "chgrp khanlab $FILE.annotated.txt\n";
print OUT "echo -e \"Annotation Pipeline finished on File $FILE\\n Please collect result from $DIR\\n\\nRegards, \\nOncogenomics Section\\nNCI\" |mutt  -s \"Annotation Pipeline Status\"  `whoami`\@mail.nih.gov\n";
close OUT;
if($node =~ /norm/){
	my $annovar =`qsub -l nodes=1:g24:c24 $DIR/annovar_$FILE.sh`;
	my $SIFT    =`qsub -l nodes=1:g24:c24 $DIR/SIFT_$FILE.sh`;
	my $PPH     =`qsub -l nodes=1:g24:c24 $DIR/PPH_$FILE.sh`;
	chomp $annovar;
	chomp $SIFT;
	chomp $PPH;
	my $join    =`qsub -l nodes=1:g24:c16 -W depend=afterany:$annovar:$SIFT:$PPH $DIR/join_$FILE.sh`;
	print "$join\n";
}
elsif($node =~ /ccr/){
	my $annovar =`qsub -q ccr -l nodes=1 $DIR/annovar_$FILE.sh`;
	my $SIFT    =`qsub -q ccr -l nodes=1 $DIR/SIFT_$FILE.sh`;
	my $PPH     =`qsub -q ccr -l nodes=1 $DIR/PPH_$FILE.sh`;
	chomp $annovar;
	chomp $SIFT;
	chomp $PPH;
	my $join    =`qsub -q ccr -l nodes=1 -W depend=afterany:$annovar:$SIFT:$PPH $DIR/join_$FILE.sh`;
	print "$join\n";

}
else{
	print "Node $node is not allowed\n";
	die;
}
