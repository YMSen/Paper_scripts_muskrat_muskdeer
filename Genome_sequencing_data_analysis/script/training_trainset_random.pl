#!usr/bin/perl -w
use strict;
die "Usage:<gff> <set_number>\n" unless @ARGV==2;
my $gff=shift;
my $lim=shift;

my %gene;
open IN,$gff;
  while(<IN>){
    chomp;
    my @b=split /\t/;
    my $name = $1 if ($b[8] =~ /=(\S+);/);
    $gene{$name} .= "$_\n";
  }
close IN;

my $num;
open OUT,">$gff.random.gff";
foreach my $k (keys %gene){
  $num++;
  if($num > $lim){
    close OUT;
	exit;
  }
  print OUT $gene{$k};
}
