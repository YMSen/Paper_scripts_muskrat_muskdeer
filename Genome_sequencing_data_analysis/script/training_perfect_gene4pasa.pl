#!/usr/bin/perl

=head1 Name
	perfect_gene.pl

=head1 Command-line Option
	--sco	set the score of genes.dafault frac=95
	--start	set the max numbers of the amino acid extending the start.default star=10
	--stop	set the max numbers of the amino acid extending the end.default end=10
	--outdir
	--help
	
=head1 Usage
	perl perfect_gene.pl [--sco] [--start] [--stop] [--outdir] <genome_fa> <gff> [<gff> ...] 

=head1 Example
	perl perfect_gene.pl --sco 99 --start 5 --stop 5 seqence.fa genewise1.gff genewise2.gff genewise3.gff

=cut



use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);

my ($Start,$Stop,$Outdir,$Help);
GetOptions(
	"start:s"=>\$Start,
	"stop:s"=>\$Stop,
	"outdir:s"=>\$Outdir,
	"help"=>\$Help
);
$Start ||= 10;
$Stop ||= 10;
$Outdir ||= "./";
$Outdir =~ s/\/$//;

die `pod2text $0` if ($Help || @ARGV < 2);

my ($file_genome,@file_gff)=@ARGV;
my $file_genome_name=basename($file_genome);

my $mult_exon		="./script/training_mult-exon_err.pl";
my $fishInWinter	="./script/fishInWinter.pl";
my $add_start_stop	="./script/training_add_start_stop.pl";

for(my $i=0;$i < @file_gff;$i++){
	my $file_gff_sub=$file_gff[$i];
	my $file=basename($file_gff_sub);
	`perl $mult_exon $Outdir/$file $file_genome > $Outdir/$file.stop`;
	`perl $fishInWinter --bf gff --ff gff --except $Outdir/$file.stop $Outdir/$file_gff_sub > $Outdir/$file.stop.un`;
	`perl $add_start_stop $file_genome $Outdir/$file.stop.un $Start $Stop > $Outdir/$file.stop.un.$Start\_$Stop.gff`;
}

`cat $Outdir/*$Start\_$Stop.gff > $Outdir/$file_genome_name.gff`;
