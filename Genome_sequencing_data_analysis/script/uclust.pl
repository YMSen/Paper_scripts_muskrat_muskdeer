#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

use File::Basename qw(basename dirname);

my ($select_spe) = "RepeatModeler;Piler;RepeatScout;LTR_FINDER";
GetOptions(
    "type:s" => \$select_spe
    );
my @denovo = split/;/,$select_spe;
my $file = shift;
my $outdir = dirname($file);
$outdir =~ s{/+$}{};
my $clu = "$outdir/new.uc.result";
my $id = "$outdir/del.list";
`./usearch -cluster_fast $file -uc $clu -id 0.8 -query_cov 0.8 -target_cov 0.8 -minqt 0.8`;

open IA,"$clu" || die "can't open IA,$!\n";
open OUT,">$id" || die "can't open OUT,$!\n";
while (<IA>){
	my @arr = split;
	next unless ($arr[0]=~/H/);
	if($arr[8]=~/^RepeatModeler|^Piler|^RepeatScout|^LTR_FINDER/ && $arr[9]!~/^RepeatModeler|^Piler|^RepeatScout|^LTR_FINDER/){
		print OUT "$arr[8]\n";
	}
	if($arr[9]=~/^RepeatModeler|^Piler|^RepeatScout|^LTR_FINDER/ && $arr[8]!~/^RepeatModeler|^Piler|^RepeatScout|^LTR_FINDER/){   
		print OUT "$arr[9]\n";
	}
}
close IA;
close OUT;

open IB, "$file" || die "Can't open IB : $!\n";
	$/ = ">";<IB>;
	while (my $seq = <IB>) {
		chomp $seq;
		open IC,"$id"||  die "Can't open IC : $!\n";
		$/ = "\n";
		my $flag=0;
		while(my $ID=<IC>){
			chomp($ID);
				#print "$ID\n";
			if ($seq=~/^$ID/){
				#print "match $ID\n";
				$flag=1;
			}
		}
		if ($flag == 0){
			print ">$seq\n";
			}
		
		$/ = ">";
	}
	$/ = "\n";
	close IB;
	close IC;	
		
	
