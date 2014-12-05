#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
Getopt::Long::Configure("bundling"); # allows option bundling
use Proc::Background;
use File::Spec;
use File::Basename;
require File::Spec::Unix;
use File::Temp qw(tempfile);
use strict;

$|=1;


my %OPT;
my %PAR;
my %ANNOVAR_OPTIONS;
my %ANNOVAR_FLAGS;
$PAR{exe}="/usr/local/apps/ANNOVAR/2014-07-14/annotate_variation.pl --buildver hg19";
$PAR{defaultData} = "/data/khanlab/ref/annovar/humandb";
set_global_hash();
get_options();

# Add the user-defined ANNOVAR options and flags
foreach my $i (sort keys %ANNOVAR_OPTIONS) { $PAR{exe} .= " --$i $ANNOVAR_OPTIONS{$i}" if defined($ANNOVAR_OPTIONS{$i}); }
foreach my $i (sort keys %ANNOVAR_FLAGS) { $PAR{exe} .= " --$i" if defined($ANNOVAR_FLAGS{$i}); }

# Run all of this in $OPT{scratch}
chdir $OPT{scratch};

populate_hg18snpA();
populate_hg19snpA();
populate_hg19cadd();

my $HR = set_global_hash();

my $FINAL; # summation of all annovar runs
# add each variant as a blank entry
open FILE, "<$OPT{input}";
while (<FILE>) {
  chomp;
  my @t = split /\s+/,$_;
  $FINAL->{$t[0]}{$t[1]}{$t[2]}{$t[3]}{$t[4]}{0} = 1;
}
close FILE;

my @cmds;

foreach my $dbtype (sort keys %{$HR}) {
# see if its already done
  my $done = 1;
  foreach my $outfile (@{$HR->{$dbtype}{outfiles}}) {
    if (! -e $outfile) { 
      print STDERR "$outfile is NOT DONE\n" unless $OPT{debug};
      undef $done; 
    }
    else {
      print STDERR "$outfile is done\n" unless $OPT{debug};
    }
  }
  if (! $done) {
    push @cmds,$HR->{$dbtype}{cmd};
    print $HR->{$dbtype}{cmd}."\n" if $OPT{debug};
  }
}

if ($OPT{showcmds}) {
  foreach my $c (@cmds) {
    print "-"x80,"\n$c\n";
  }
  print "-"x80,"\n";
}

exit if $OPT{debug};

pseudo_swarm($OPT{threads},@cmds); # run the commands

# extract the columns
foreach my $dbtype (@{$PAR{dbtypes_ordered}}) { # maintain order
  foreach my $outfile (@{$HR->{$dbtype}{outfiles}}) {
    if ($ANNOVAR_FLAGS{otherinfo} && defined $HR->{$dbtype}{out}{$outfile}{otherinfolabels}) {
      foreach my $otherinfolabel (@{$HR->{$dbtype}{out}{$outfile}{otherinfolabels}}) {
        push @{$PAR{heads}},$otherinfolabel if ((not defined $PAR{heads}) || (!(grep /^$otherinfolabel$/,@{$PAR{heads}})));
      }
      otherinfo_extract_cols($outfile,$HR->{$dbtype}{out}{$outfile});
    }
    else {
      foreach my $label (@{$HR->{$dbtype}{out}{$outfile}{labels}}) {
        push @{$PAR{heads}},$label if ((not defined $PAR{heads}) || (!(grep /^$label$/,@{$PAR{heads}})));
      }
      extract_cols($outfile,$HR->{$dbtype}{out}{$outfile});
    } 
  }
}

# create header row
open OUTFILE, ">$OPT{output}" or warn "Can't open $OPT{output}\n";

printf OUTFILE ("%s\t%s\n","Chr\tStart\tEnd\tRef\tAlt",(join "\t",(@{$PAR{heads}})));

# print out final spreadsheet
foreach my $chr (sort {$a <=> $b} keys %{$FINAL}) {
  foreach my $start (sort {$a <=> $b} keys %{$FINAL->{$chr}}) {
    foreach my $end (sort {$a <=> $b} keys %{$FINAL->{$chr}{$start}}) {
      foreach my $ref (sort keys %{$FINAL->{$chr}{$start}{$end}}) {
        foreach my $alt (sort keys %{$FINAL->{$chr}{$start}{$end}{$ref}}) {
          my $hr = $FINAL->{$chr}{$start}{$end}{$ref}{$alt};
          print OUTFILE "$chr\t$start\t$end\t$ref\t$alt";
          foreach my $head (@{$PAR{heads}}) {
            print OUTFILE "\t".$hr->{$head};
          }
          print OUTFILE "\n";
        }
      }
    }
  }
}
close OUTFILE;

#==============================================================================
sub otherinfo_extract_cols
{
  my ($file,$hr) = @_;

  my @labels = @{$hr->{labels}};
  my @otherinfolabels = @{$hr->{otherinfolabels}};
  my @values;
  foreach my $label (@labels) {
    push @values,$hr->{$label};
  }

  open FILE, "<$file" or warn "Can't open $file\n";
  while (<FILE>) {
    chomp;
    my @tmp = split /\t/,$_;
    foreach my $i (0 .. $#labels) {
      if (defined $tmp[$values[$i]]) {
        my @tmp2 = split /,/,$tmp[$values[$i]];
        foreach my $j (0 .. $#otherinfolabels) {
          $FINAL->{$tmp[$hr->{Chr}]}{$tmp[$hr->{Start}]}{$tmp[$hr->{End}]}{$tmp[$hr->{Ref}]}{$tmp[$hr->{Alt}]}{$otherinfolabels[$j]} = $tmp2[$j];
        }
      }
    }
  }
  close FILE;
}
#==============================================================================
sub extract_cols
{
  my ($file,$hr) = @_;

  my @labels = @{$hr->{labels}};
  my @values;
  foreach my $label (@labels) {
    push @values,$hr->{$label};
  }

  open FILE, "<$file" or warn "Can't open $file\n";
  while (<FILE>) {
    chomp;
    my @tmp = split /\t/,$_;
    foreach my $i (0 .. $#labels) {
      if (defined $tmp[$values[$i]]) {
        $FINAL->{$tmp[$hr->{Chr}]}{$tmp[$hr->{Start}]}{$tmp[$hr->{End}]}{$tmp[$hr->{Ref}]}{$tmp[$hr->{Alt}]}{$labels[$i]} = $tmp[$values[$i]]; 
      }
    }
  }
  close FILE;
}
#==============================================================================
sub get_options {
  $OPT{data} = $PAR{defaultData} unless defined $OPT{data};
  &GetOptions(
    "h"  => \$OPT{help},
    "help"  => \$OPT{help},
    "d"=> \$OPT{debug},
    "debug"=> \$OPT{debug},
    "verbose"=> \$OPT{verbose},
    "data=s"=> \$OPT{data},
    "input=s"=> \$OPT{input},
    "output=s"=> \$OPT{output},
    "scratch=s"=> \$OPT{scratch},
    "threads=i"=> \$OPT{threads},
    "generic=s@"=> \$OPT{generic},
    "dbtype=s@"=> \$OPT{dbtypes},
    "autothread"=> \$OPT{autothread},
    "showcmds"=> \$OPT{showcmds},

# ANNOVAR_OPTIONS

    "batchsize=i"=>\$ANNOVAR_OPTIONS{batchsize}, # batch size for processing variants per batch (default: 5m)
    "genomebinsize=i"=>\$ANNOVAR_OPTIONS{genomebinsize}, # bin size to speed up search (default: 100k for -geneanno, 10k for -regionanno)
    "expandbin=i"=>\$ANNOVAR_OPTIONS{expandbin}, # check nearby bin to find neighboring genes (default: 2m/genomebinsize)
    "neargene=i"=>\$ANNOVAR_OPTIONS{neargene}, # distance threshold to define upstream/downstream of a gene
    "exonicsplicing"=>\$ANNOVAR_FLAGS{exonicsplicing}, # report exonic variants near exon/intron boundary as 'exonic;splicing' variants
    "score_threshold=f"=>\$ANNOVAR_OPTIONS{score_threshold}, # minimum score of DB regions to use in annotation
    "reverse"=>\$ANNOVAR_FLAGS{reverse}, # reverse directionality to compare to score_threshold
    "normscore_threshold=f"=>\$ANNOVAR_OPTIONS{normscore_threshold}, # minimum normalized score of DB regions to use in annotation
    "rawscore"=>\$ANNOVAR_FLAGS{rawscore}, # output includes the raw score (not normalized score) in UCSC Browser Track
    "minqueryfrac=i"=>\$ANNOVAR_OPTIONS{minqueryfrac}, # minimum percentage of query overlap to define match to DB (default: 0)
    "splicing_threshold=i"=>\$ANNOVAR_OPTIONS{splicing_threshold}, # distance between splicing variants and exon/intron boundary (default: 2)
    "indel_splicing_threshold=i"=>\$ANNOVAR_OPTIONS{indel_splicing_threshold}, # if set, use this value for allowed indel size for splicing variants (default: --splicing_threshold)
    "maf_threshold=f"=>\$ANNOVAR_OPTIONS{maf_threshold}, # filter 1000G variants with MAF above this threshold (default: 0)
    "sift_threshold=f"=>\$ANNOVAR_OPTIONS{sift_threshold}, # SIFT threshold for deleterious prediction for -dbtype avsift (default: 0.05)
    "precedence=s"=>\$ANNOVAR_OPTIONS{precedence}, # comma-delimited to specify precedence of variant function (default: exonic>intronic...)
    "indexfilter_threshold=f"=>\$ANNOVAR_OPTIONS{indexfilter_threshold}, # controls whether filter-based annotation use index if this fraction of bins need to be scanned (default: 0.9)

  ) || print_options("ERROR");

  $OPT{scratch} = "/scratch" if (not defined $OPT{scratch}); 
  $ANNOVAR_FLAGS{otherinfo} = 1; # this is the default

  print_options() if $OPT{help};

  $OPT{threads} = 1 if (not defined $OPT{threads});
  die "No input file given!\n" if !$OPT{input};
  die "No output file given!\n" if !$OPT{output};

# find absolute paths for files
  $OPT{scratch} = File::Spec->rel2abs($OPT{scratch});
  $OPT{input} = File::Spec->rel2abs($OPT{input});
  $OPT{output} = File::Spec->rel2abs($OPT{output});

  die "Can't open $OPT{input}\n" if (!-r $OPT{input});
  die "Output file $OPT{output} already exists!\n" if (-e $OPT{output});
  die "scratch dir $OPT{scratch} unusable!\n" if (!-d $OPT{scratch});
  die "scratch dir $OPT{scratch} unusable!\n" if (!-w $OPT{scratch});
  die "scratch dir $OPT{scratch} unusable!\n" if (!-r $OPT{scratch});
  die "scratch dir $OPT{scratch} unusable!\n" if (!-x $OPT{scratch});

  if ((defined $OPT{threads}) && ($OPT{threads} < 1)) { print_options("threads must be > 0!"); }
  if ((defined $OPT{threads}) && ($OPT{threads} > 16)) {
    warn "threads must be <= 16\n";
    $OPT{threads} = 16;
  }
  if (defined $OPT{autothread}) {
    $OPT{threads} = howManyProc();
  }
  if (!-d $OPT{data}) { die "No such dir $OPT{data}!\n"; }
  $OPT{data} = File::Spec->rel2abs($OPT{data});

  my @tmp;
  if (defined $OPT{generic}) {
    foreach my $i (@{$OPT{generic}}) {
      push @tmp,File::Spec->rel2abs($i);
    }
    @{$OPT{generic}}=@tmp;
  }
}
#==============================================================================
sub print_options {
  my $mess = shift;
  print "\n".$mess."\n\n" if $mess;
  print <<EOF1;
usage: $0 --input file --output file [ options ]
EOF1

  exit if $mess;

  print <<EOF2;

Run annovar against multiple dbs.  Currently the build is hg19 and these
dbtypes are used:
EOF2

  print "\n  --geneanno:\n\n";

  foreach my $dbname (sort keys %{$PAR{dbnames}}) {
    if ($PAR{dbnames}{$dbname}{type} eq "geneanno") {
      print "    $dbname\n";
    }
  }

  print "\n  --filter:\n\n";
  
  foreach my $dbname (sort keys %{$PAR{dbnames}}) {
    if ($PAR{dbnames}{$dbname}{type} eq "filter") {
      print "    $dbname\n";
    }
  }

  print "\n  --regionanno:\n\n";

  foreach my $dbname (sort keys %{$PAR{dbnames}}) {
    if ($PAR{dbnames}{$dbname}{type} eq "regionanno") {
      print "    $dbname\n";
    }
  }

  print <<EOF3;

Annovar is run against all these dbtypes, and a tab-delimited file containing
all the results is written.
    
options:

  --input      input file for annovar
  --output     output file
  --threads    run with more than one thread (max = 16)
  --autothread run with as many threads as are available on the node
  --dbtype     pick and choose dbtypes (can be more than one)
                 not selecting a dbtype defaults to all dbtypes
                 (see above for list)

  --generic    path to a single filter database file, must end in '.txt'

  --data       use a different data dir
                 (default = $OPT{data})
  --scratch    scratch dir where work is done
                 (default = $OPT{scratch})

  --showcmds   display the annotate_variation.pl commands
  -h, --help   show this menu
  -d, --debug  run in debug mode

Additionally, the following annotate_variation.pl options are also available:

  --batchsize
  --genomebinsize
  --expandbin
  --neargene
  --exonicsplicing
  --score_threshold
  --reverse
  --normscore_threshold
  --rawscore
  --minqueryfrac
  --splicing_threshold
  --indel_splicing_threshold
  --maf_threshold
  --sift_threshold
  --precedence
  --indexfilter_threshold

The annotate_variation.pl option --otherinfo is used by default.

Last updated: Apr 29, 2014 (David Hoover)

EOF3
  exit;
}
#==============================================================================
sub pseudo_swarm
#
# This subroutine accepts two arguments: a list of system commands, and a
# maximum number to run at any given time.  It submits each in the list to
# the system until the maximum number of jobs is reached.  When one of the
# previously submitted jobs finishes, it submits the next, unless there are
# no more jobs to submit.  When all jobs are completed, it returns.
{
  my $n = shift;

# populate initial pids

  my %jobs;
  while ($n) { $jobs{$n} = Proc::Background->new(shift) if (@_); $n--; }

# if there are any pids to go keep an eye on each pid

  while (%jobs) {
    foreach my $n (sort keys %jobs) {

# if one pid dies, replace it

      if (!$jobs{$n}->alive) {
        delete $jobs{$n};
        $jobs{$n} = Proc::Background->new(shift) if (@_);
      }
    }
    sleep 2;  # you are getting sleepy
  }
  sleep 2;    # you are getting very sleepy

  return;
}
#==============================================================================
sub set_global_hash
{

  my $outtag = "1";
  my (undef,$tmptag) = tempfile("XXXXXXXX",OPEN=>0); # create a random tag
# --geneanno types
  $PAR{dbnames}{"gene"} = {("type" => "geneanno")};
  $PAR{dbtypes}{gene} = {(
    "cmd" => "$PAR{exe} --geneanno --outfile ${tmptag}_$outtag --dbtype gene $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.variant_function" => {(
        "gene.Region" => 0,
        "Gene" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["gene.Region","Gene"],
      )},
      "${tmptag}_$outtag.exonic_variant_function" => {
        "gene.VariantType" => 1,
        "gene.AminoAcidChange" => 2,
        "Chr" => 3, "Start" => 4, "End" => 5, "Ref" => 6, "Alt" => 7,
        "labels" => ["gene.VariantType","gene.AminoAcidChange"],
      },
    )},
    "outfiles" => ["${tmptag}_$outtag.variant_function","${tmptag}_$outtag.exonic_variant_function"],
  )};

  $outtag++;
  $PAR{dbnames}{"ensGene"} = {("type" => "geneanno")};
  $PAR{dbtypes}{ensGene} = {(
    "cmd" => "$PAR{exe} --geneanno --outfile ${tmptag}_$outtag --dbtype ensGene $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.variant_function" => {(
        "ensGene.Region" => 0,
        "ensGene" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["ensGene.Region","ensGene"],
      )},
      "${tmptag}_$outtag.exonic_variant_function" => {
        "ensGene.VariantType" => 1,
        "ensGene.AminoAcidChange" => 2,
        "Chr" => 3, "Start" => 4, "End" => 5, "Ref" => 6, "Alt" => 7,
        "labels" => ["ensGene.VariantType","ensGene.AminoAcidChange"],
      },
    )},
    "outfiles" => ["${tmptag}_$outtag.variant_function","${tmptag}_$outtag.exonic_variant_function"],
  )};

  $outtag++;
  $PAR{dbnames}{"knownGene"} = {("type" => "geneanno")};
  $PAR{dbtypes}{knownGene} = {(
    "cmd" => "$PAR{exe} --geneanno --outfile ${tmptag}_$outtag --dbtype knownGene $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.variant_function" => {(
        "knownGene.Region" => 0,
        "knownGene" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["knownGene.Region","knownGene"],
      )},
      "${tmptag}_$outtag.exonic_variant_function" => {
        "knownGene.VariantType" => 1,
        "knownGene.AminoAcidChange" => 2,
        "Chr" => 3, "Start" => 4, "End" => 5, "Ref" => 6, "Alt" => 7,
        "labels" => ["knownGene.VariantType","knownGene.AminoAcidChange"],
      },
    )},
    "outfiles" => ["${tmptag}_$outtag.variant_function","${tmptag}_$outtag.exonic_variant_function"],
  )};

# --filter types

# 1000 genome stuff
  $outtag++;
  $PAR{dbnames}{"1000g2014oct_all"} = {("type" => "filter")};
  $PAR{dbtypes}{"1000g2014oct_all"} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype 1000g2014oct_all $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_ALL.sites.2014_10_dropped" => {(
        "1000g08.2014" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["1000g08.2014"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_ALL.sites.2014_10_dropped"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"1000g2014oct_afr"} = {("type" => "filter")};
  $PAR{dbtypes}{"1000g2014oct_afr"} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype 1000g2014oct_afr $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_AFR.sites.2014_10_dropped" => {(
        "1000g08.2014_afr" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["1000g08.2014_afr"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_AFR.sites.2014_10_dropped"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"1000g2014oct_amr"} = {("type" => "filter")};
  $PAR{dbtypes}{"1000g2014oct_amr"} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype 1000g2014oct_amr $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_AMR.sites.2014_10_dropped" => {(
        "1000g08.2014_amr" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["1000g08.2014_amr"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_AMR.sites.2014_10_dropped"],
  )};
 # HERE 
  $outtag++;
  $PAR{dbnames}{"1000g2014oct_eas"} = {("type" => "filter")};
  $PAR{dbtypes}{"1000g2014oct_eas"} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype 1000g2014oct_eas $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_EAS.sites.2014_10_dropped" => {(
        "1000g08.2014_eas" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["1000g08.2014_eas"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_EAS.sites.2014_10_dropped"],
  )};

  $outtag++;
  $PAR{dbnames}{"1000g2014oct_sas"} = {("type" => "filter")};
  $PAR{dbtypes}{"1000g2014oct_sas"} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype 1000g2014oct_sas $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_SAS.sites.2014_10_dropped" => {(
        "1000g08.2014_sas" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["1000g08.2014_sas"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_SAS.sites.2014_10_dropped"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"1000g2014oct_eur"} = {("type" => "filter")};
  $PAR{dbtypes}{"1000g2014oct_eur"} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype 1000g2014oct_eur $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_EUR.sites.2014_10_dropped" => {(
        "1000g08.2014_eur" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["1000g08.2014_eur"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_EUR.sites.2014_10_dropped"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"esp6500_all"} = {("type" => "filter")};
  $PAR{dbtypes}{esp6500_all} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype esp6500_all $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_esp6500_all_dropped" => {(
        "NHLBI6500_ALL" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["NHLBI6500_ALL"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_esp6500_all_dropped"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"esp6500_aa"} = {("type" => "filter")};
  $PAR{dbtypes}{esp6500_aa} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype esp6500_aa $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_esp6500_aa_dropped" => {(
        "NHLBI6500_AA" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["NHLBI6500_AA"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_esp6500_aa_dropped"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"esp6500_ea"} = {("type" => "filter")};
  $PAR{dbtypes}{esp6500_ea} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype esp6500_ea $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_esp6500_ea_dropped" => {(
        "NHLBI6500_EA" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["NHLBI6500_EA"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_esp6500_ea_dropped"],
  )};

# snp132 breakup 
  $PAR{dbnames}{"snp132"} = {("type" => "filter")};
  foreach my $n (@{$PAR{hg19snpA}}) {
    $outtag++;
    $PAR{dbtypes}{"snp132${n}"} = {(
      "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype snp132${n} $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
      "out" => {(
        "${tmptag}_$outtag.hg19_snp132${n}_dropped" => {(
          "dbsnp132" => 1,
          "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
          "labels" => ["dbsnp132"],
        )},
      )},
      "outfiles" => ["${tmptag}_$outtag.hg19_snp132${n}_dropped"],
    )};
  } 
# snp135 breakup 
  $PAR{dbnames}{"snp135"} = {("type" => "filter")};
  foreach my $n (@{$PAR{hg19snpA}}) {
    $outtag++;
    $PAR{dbtypes}{"snp135${n}"} = {(
      "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype snp135${n} $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
      "out" => {(
        "${tmptag}_$outtag.hg19_snp135${n}_dropped" => {(
          "dbsnp135" => 1,
          "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
          "labels" => ["dbsnp135"],
        )},
      )},
      "outfiles" => ["${tmptag}_$outtag.hg19_snp135${n}_dropped"],
    )};
  } 
# snp137 breakup 
  $PAR{dbnames}{"snp137"} = {("type" => "filter")};
  foreach my $n (@{$PAR{hg19snpA}}) {
    $outtag++;
    $PAR{dbtypes}{"snp137${n}"} = {(
      "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype snp137${n} $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
      "out" => {(
        "${tmptag}_$outtag.hg19_snp137${n}_dropped" => {(
          "dbsnp137" => 1, "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
          "labels" => ["dbsnp137"],
        )},
      )},
      "outfiles" => ["${tmptag}_$outtag.hg19_snp137${n}_dropped"],
    )};
  } 
# snp138 breakup 
  $PAR{dbnames}{"snp138"} = {("type" => "filter")};
  foreach my $n (@{$PAR{hg19snpA}}) {
    $outtag++;
    $PAR{dbtypes}{"snp138${n}"} = {(
      "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype snp138${n} $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
      "out" => {(
        "${tmptag}_$outtag.hg19_snp138${n}_dropped" => {(
          "dbsnp138" => 1, "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
          "labels" => ["dbsnp138"],
        )},
      )},
      "outfiles" => ["${tmptag}_$outtag.hg19_snp138${n}_dropped"],
    )};
  } 

  $PAR{dbnames}{"cadd"} = {("type" => "filter")};
  foreach my $n (@{$PAR{hg19cadd}}) {
    $outtag++;
    $PAR{dbtypes}{"cadd${n}"} = {(
      "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype cadd${n} $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
      "out" => {(
        "${tmptag}_$outtag.hg19_cadd${n}_dropped" => {(
          "caddRAW" => 1,
          "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
          "labels" => ["caddRAW"],
          "otherinfolabels" => ["caddRAW","caddPHRED"],
        )},
      )},
      "outfiles" => ["${tmptag}_$outtag.hg19_cadd${n}_dropped"],
    )};
  } 

#LJB23 (LJB version 2.3) non-synonymous variants annotation
# ljb23_all
  $outtag++;
  $PAR{dbnames}{"ljb23_all"} = {("type" => "filter")};
  $PAR{dbtypes}{ljb23_all} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype ljb23_all $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_ljb23_all_dropped" => {(
        "LJB23_all" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["LJB23_all"],
        "otherinfolabels" => [ "LJB23_SIFT_score","LJB23_SIFT_score_converted","LJB23_SIFT_pred",
          "LJB23_Polyphen2_HDIV_score","LJB23_Polyphen2_HDIV_pred",
          "LJB23_Polyphen2_HVAR_score","LJB23_Polyphen2_HVAR_pred",
          "LJB23_LRT_score","LJB23_LRT_score_converted","LJB23_LRT_pred",
          "LJB23_MutationTaster_score","LJB23_MutationTaster_score_converted","LJB23_MutationTaster_pred",
          "LJB23_MutationAssessor_score","LJB23_MutationAssessor_score_converted","LJB23_MutationAssessor_pred",
          "LJB23_FATHMM_score","LJB23_FATHMM_score_converted","LJB23_FATHMM_pred",
          "LJB23_RadialSVM_score","LJB23_RadialSVM_score_converted","LJB23_RadialSVM_pred",
          "LJB23_LR_score","LJB23_LR_pred",
          "LJB23_GERP++",
          "LJB23_PhyloP",
          "LJB23_SiPhy",],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_ljb23_all_dropped"],
  )};

# avsift
  $outtag++;
  $PAR{dbnames}{"avsift"} = {("type" => "filter")};
  $PAR{dbtypes}{avsift} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype avsift $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_avsift_dropped" => {(
        "avsift" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["avsift"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_avsift_dropped"],
  )};

  $outtag++;
  $PAR{dbnames}{"cosmic70"} = {("type" => "filter")};
  $PAR{dbtypes}{cosmic70} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype cosmic70 $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_cosmic70_dropped" => {(
        "cosmic70" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["cosmic70"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_cosmic70_dropped"],
  )};

  $outtag++;
  $PAR{dbnames}{"clinvar_20140702"} = {("type" => "filter")};
  $PAR{dbtypes}{clinvar_20140702} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype clinvar_20140702 $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_clinvar_20140702_dropped" => {(
        "clinvar_20140702" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["clinvar_20140702"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_clinvar_20140702_dropped"],
  )};

  $outtag++;
  $PAR{dbnames}{"nci60"} = {("type" => "filter")};
  $PAR{dbtypes}{nci60} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype nci60 $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_nci60_dropped" => {(
        "nci60" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["nci60"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_nci60_dropped"],
  )};

  $outtag++;
  $PAR{dbnames}{"cg69"} = {("type" => "filter")};
  $PAR{dbtypes}{cg69} = {(
    "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype cg69 $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_cg69_dropped" => {(
        "cg69" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["cg69"],
      )},
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_cg69_dropped"],
  )};

# --regionanno types
  $outtag++;
  $PAR{dbnames}{"band"} = {("type" => "regionanno")};
  $PAR{dbtypes}{band} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype band $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_cytoBand" => {(
        "CytoBand" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["CytoBand"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_cytoBand"],
  )};

  $outtag++;
  $PAR{dbnames}{"dgv"} = {("type" => "regionanno")};
  $PAR{dbtypes}{dgv} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype dgv $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_dgv" => {(
        "dbGenomicVariants" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["dbGenomicVariants"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_dgv"],
  )};

  $outtag++;
  $PAR{dbnames}{"segdup"} = {("type" => "regionanno")};
  $PAR{dbtypes}{segdup} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype segdup $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_genomicSuperDups" => {(
        "genomicSuperDuplicate" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["genomicSuperDuplicate"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_genomicSuperDups"],
  )};

  $outtag++;
  $PAR{dbnames}{"gwascatalog"} = {("type" => "regionanno")};
  $PAR{dbtypes}{gwascatalog} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype gwascatalog $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_gwasCatalog" => {(
        "GWAS" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["GWAS"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_gwasCatalog"],
  )};

  $outtag++;
  $PAR{dbnames}{"mce46way"} = {("type" => "regionanno")};
  $PAR{dbtypes}{mce46way} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype mce46way $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_phastConsElements46way" => {(
        "MostConservedEle_46_Vertebrates" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["MostConservedEle_46_Vertebrates"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_phastConsElements46way"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"tfbs"} = {("type" => "regionanno")};
  $PAR{dbtypes}{tfbs} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype tfbs $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_tfbsConsSites" => {(
        "TFBS" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["TFBS"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_tfbsConsSites"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"mirna"} = {("type" => "regionanno")};
  $PAR{dbtypes}{mirna} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype mirna $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_wgRna" => {(
        "miRNA" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["miRNA"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_wgRna"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"evofold"} = {("type" => "regionanno")};
  $PAR{dbtypes}{evofold} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype evofold $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_evofold" => {(
        "EvoFold" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["EvoFold"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_evofold"],
  )};
  
  $outtag++;
  $PAR{dbnames}{"mirnatarget"} = {("type" => "regionanno")};
  $PAR{dbtypes}{mirnatarget} = {(
    "cmd" => "$PAR{exe} --regionanno --outfile ${tmptag}_$outtag --dbtype mirnatarget $OPT{input} $OPT{data} &> /dev/null ; cat ${tmptag}_$outtag.log",
    "out" => {(
      "${tmptag}_$outtag.hg19_targetScanS" => {(
        "miRNA_TargetSite" => 1,
        "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
        "labels" => ["miRNA_TargetSite"],
      )}
    )},
    "outfiles" => ["${tmptag}_$outtag.hg19_targetScanS"],
  )};
  
  if (defined $OPT{generic}) {
    foreach my $path (@{$OPT{generic}}) {
      my $dir = dirname($path);
      my $file = basename($path);
      my $name;
      if ($file=~/(.*)\.txt$/) {
        $name = $1 
      }
      if ($dir && $file && $name) {
        $outtag++;
        $PAR{dbtypes}{$name} = {(
          "type" => "filter",
          "cmd" => "$PAR{exe} --filter --outfile ${tmptag}_$outtag --dbtype generic --genericdbfile $file $OPT{input} $dir &> /dev/null ; cat ${tmptag}_$outtag.log",
          "out" => {(
            "${tmptag}_$outtag.hg19_generic_dropped" => {(
              "$name" => 1,
              "Chr" => 2, "Start" => 3, "End" => 4, "Ref" => 5, "Alt" => 6,
              "labels" => ["$name"],
            )},
          )},
          "outfiles" => ["${tmptag}_$outtag.hg19_generic_dropped"],
        )};
        push @{$OPT{dbtypes}},$name;
      }
    }
  }

  my $hr;

# default ordering
  @{$PAR{dbtypes_ordered}} = sort keys %{$PAR{dbtypes}};
  @{$PAR{default_dbtypes}} = @{$PAR{dbtypes_ordered}};


# maintain ordering given by user
  if (defined $OPT{dbtypes}) {
    undef $PAR{dbtypes_ordered};
    foreach my $i (@{$OPT{dbtypes}}) {
      if ($i =~ /^snp132$/) { # snp132 breakup
        foreach my $n (@{$PAR{hg19snpA}}) {
          push @{$PAR{dbtypes_ordered}},"snp132${n}";
          $hr->{"snp132${n}"} = $PAR{dbtypes}{"snp132${n}"};
        }
      }
      elsif ($i =~ /^snp135$/) { # snp135 breakup
        foreach my $n (@{$PAR{hg19snpA}}) {
          push @{$PAR{dbtypes_ordered}},"snp135${n}";
          $hr->{"snp135${n}"} = $PAR{dbtypes}{"snp135${n}"};
        }
      }
      elsif ($i =~ /^snp137$/) { # snp137 breakup
        foreach my $n (@{$PAR{hg19snpA}}) {
          push @{$PAR{dbtypes_ordered}},"snp137${n}";
          $hr->{"snp137${n}"} = $PAR{dbtypes}{"snp137${n}"};
        }
      }
      elsif ($i =~ /^snp138$/) { # snp138 breakup
        foreach my $n (@{$PAR{hg19snpA}}) {
          push @{$PAR{dbtypes_ordered}},"snp138${n}";
          $hr->{"snp138${n}"} = $PAR{dbtypes}{"snp138${n}"};
        }
      }
      elsif ($i =~ /^cadd$/) { # cadd breakup
        foreach my $n (@{$PAR{hg19cadd}}) {
          push @{$PAR{dbtypes_ordered}},"cadd${n}";
          $hr->{"cadd${n}"} = $PAR{dbtypes}{"cadd${n}"};
        }
      }
      elsif (grep /^$i$/,@{$PAR{default_dbtypes}}) {
        push @{$PAR{dbtypes_ordered}},$i;
        $hr->{$i} = $PAR{dbtypes}{$i};
      }
    }   

# Make sure that the dbtypes chosen actually exist!

    my $bomb;
    foreach my $i (@{$OPT{dbtypes}}) {
      if (!(grep /^$i/,@{$PAR{default_dbtypes}})) {
        print "dbtype $i isn't available for this script\n";
        $bomb=1;
      }
    }
    exit if $bomb;
  }
  else {
    $hr = $PAR{dbtypes};
  }

  return $hr;
}
#==============================================================================
sub howManyProc
{
  open FILE, "</proc/cpuinfo";
  my $h;
  while (<FILE>) {
    $h->{phys}->{$1}=1 if (m/^(physical id.*)/);
    $h->{core}->{$1}=1 if (m/^(core id.*)/);
    $h->{proc}->{$1}=1 if (m/^(processor.*)/);
  }
  close FILE;
  my $cpu = (scalar(keys %{$h->{phys}})) * (scalar(keys %{$h->{core}}));
  $cpu = (scalar(keys %{$h->{proc}})) if !$cpu;
  return $cpu;
}
#==============================================================================
# hg19 cadd breadkup hash
sub populate_hg19cadd
{
  @{$PAR{hg19cadd}} = qw(
_1 _2 _3 _4 _5 _6 _7 _8 _9 _10 _11 _12 _13 _14 _15 _16 _17 _18 _19 _20 _21 _22 _M _X _Y
_GL000191.1 _GL000192.1 _GL000193.1 _GL000194.1 _GL000195.1 _GL000196.1 _GL000197.1 _GL000198.1 _GL000199.1
_GL000200.1 _GL000201.1 _GL000202.1 _GL000203.1 _GL000204.1 _GL000205.1 _GL000206.1 _GL000207.1 _GL000208.1
_GL000209.1 _GL000210.1 _GL000211.1 _GL000212.1 _GL000213.1 _GL000214.1 _GL000215.1 _GL000216.1 _GL000217.1
_GL000218.1 _GL000219.1 _GL000220.1 _GL000221.1 _GL000222.1 _GL000223.1 _GL000224.1 _GL000225.1 _GL000226.1
_GL000227.1 _GL000228.1 _GL000229.1 _GL000230.1 _GL000231.1 _GL000232.1 _GL000233.1 _GL000234.1 _GL000235.1
_GL000236.1 _GL000237.1 _GL000238.1 _GL000239.1 _GL000240.1 _GL000241.1 _GL000242.1 _GL000243.1 _GL000244.1
_GL000245.1 _GL000246.1 _GL000247.1 _GL000248.1 _GL000249.1
);
}
#==============================================================================
# hg19 snp breadkup hash
sub populate_hg19snpA
{
  @{$PAR{hg19snpA}} = qw(
_chr1
_chr1_gl000191_random
_chr1_gl000192_random
_chr2
_chr3
_chr4_ctg9_hap1
_chr4_gl000193_random
_chr4_gl000194_random
_chr4
_chr5
_chr6
_chr6_apd_hap1
_chr6_cox_hap2
_chr6_dbb_hap3
_chr6_mann_hap4
_chr6_mcf_hap5
_chr6_qbl_hap6
_chr6_ssto_hap7
_chr7
_chr7_gl000195_random
_chr8
_chr8_gl000196_random
_chr8_gl000197_random
_chr9
_chr9_gl000198_random
_chr9_gl000199_random
_chr9_gl000200_random
_chr9_gl000201_random
_chr10
_chr11
_chr11_gl000202_random
_chr12
_chr13
_chr14
_chr15
_chr16
_chr17
_chr17_ctg5_hap1
_chr17_gl000203_random
_chr17_gl000204_random
_chr17_gl000205_random
_chr17_gl000206_random
_chr18
_chr18_gl000207_random
_chr19
_chr19_gl000208_random
_chr19_gl000209_random
_chr20
_chr21
_chr21_gl000210_random
_chr22
_chrM
_chrUn_gl000211
_chrUn_gl000212
_chrUn_gl000213
_chrUn_gl000214
_chrUn_gl000215
_chrUn_gl000216
_chrUn_gl000217
_chrUn_gl000218
_chrUn_gl000219
_chrUn_gl000220
_chrUn_gl000221
_chrUn_gl000222
_chrUn_gl000223
_chrUn_gl000224
_chrUn_gl000225
_chrUn_gl000226
_chrUn_gl000227
_chrUn_gl000228
_chrUn_gl000229
_chrUn_gl000230
_chrUn_gl000231
_chrUn_gl000232
_chrUn_gl000233
_chrUn_gl000234
_chrUn_gl000235
_chrUn_gl000236
_chrUn_gl000237
_chrUn_gl000238
_chrUn_gl000239
_chrUn_gl000240
_chrUn_gl000241
_chrUn_gl000242
_chrUn_gl000243
_chrUn_gl000244
_chrUn_gl000245
_chrUn_gl000246
_chrUn_gl000247
_chrUn_gl000248
_chrUn_gl000249
_chrX
_chrY
);
}
#==============================================================================
# hg18 snp breakup hash
sub populate_hg18snpA
{
  @{$PAR{hg18snpA}} = qw(
_chr1
_chr2
_chr3
_chr4
_chr5
_chr6
_chr7
_chr8
_chr9
_chr10
_chr11
_chr12
_chr13
_chr14
_chr15
_chr16
_chr17
_chr18
_chr19
_chr20
_chr21
_chr22
_chrX
_chrY
_chrM
);
}
#==============================================================================
