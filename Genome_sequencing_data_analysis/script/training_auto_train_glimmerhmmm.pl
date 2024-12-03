#!/usr/bin/perl 

use strict;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Getopt::Long;
use Data::Dumper;

die "perl $0 <gff> <genome.fa> [options]\n\t-d name	name = name of training directory\n\t-outdir\tgive a dir to print the output
	-i\ti1,i2,...,in
		isochores to be considered (e.g. if two isochores are desired between
		0-40% GC content and 40-100% then the option should be: -i 0,40,100;
		default is -i 0,40,100 )\n" if(@ARGV < 2);

my ($Dname,$Outdir,$i);
GetOptions(
	"d:s"=>\$Dname,
	"outdir:s"=>\$Outdir,
	"i:s"=>\$i,
);
$Outdir ||= "./";
$Outdir=~s/\/$//;
$i ||= "0,40,100";
$Dname ||= "GHMM_train";
$Dname="$Outdir/$Dname";
my $para="-d $Dname -i $i";

my $gff=shift;
my $genome=shift;
my $gffname=basename($gff);
my $genomename=basename($genome);


`perl ./script/training_glimmerhmm_chform_train.pl $gff > $Outdir/$gffname.cds`;
`perl ./script/trainGlimmerHMM $genome $Outdir/$gffname.cds $para`;
