#!/usr/bin/perl 


use strict;
use File::Basename qw(basename);
use FindBin qw($Bin);
use Getopt::Long;

my ($Outdir,$Name,$Help);
GetOptions(
	"outdir:s"=>\$Outdir,
	"name:s"=>\$Name,
	"help"=>\$Help,
);
$Outdir ||= "./";
$Outdir=~s/\/$//;
$Name ||= "para";

die "Usage:\n\t--outdir\n\t--name\n\t--help\nperl $0 <genome.fa> <gene.gff> [options]\n" if (@ARGV < 2);

my $file_genome=shift;
my $file_gff=shift;
my $genome=basename($file_genome);
my $gff=basename($file_gff);

my $pro_dir="./script/snap-2013-11-29/snap";

###########
my $ch_form="/PUBLIC/software/DENOVO/bio/annotation/pipeline_v2.0/scripts/training_SNAP_train_chform.pl";


`perl $ch_form --outdir $Outdir $file_genome $file_gff`;
###########
`./script/snap-2013-11-29/fathom $Outdir/$gff.ann $Outdir/$genome.dna -gene-stats 1>>$Outdir/train.log 2>>$Outdir/train.erro`;
`./script/snap-2013-11-29/fathom $Outdir/$gff.ann $Outdir/$genome.dna -validate 1>>$Outdir/train.log 2>>$Outdir/train.erro`;
`./script/snap-2013-11-29/fathom $Outdir/$gff.ann $Outdir/$genome.dna -categorize 1000 1>>$Outdir/train.log 2>>$Outdir/train.erro`;
`./script/snap-2013-11-29/fathom uni.ann uni.dna -export 1000 -plus 1>>$Outdir/train.log 2>>$Outdir/train.erro`;
`mkdir $Outdir/$Name` unless (-e "$Outdir/$Name");

chomp(my $curr_dir=`pwd`);
chdir "$Outdir/$Name";
`./script/snap-2013-11-29/forge $curr_dir/export.ann $curr_dir/export.dna 1>>$Outdir/train.log 2>>$Outdir/train.err`;
chdir "$curr_dir";
`./script/snap-2013-11-29/hmm-assembler.pl $file_genome $Outdir/$Name > $Outdir/$Name.hmm 2>>$Outdir/train.err`;

