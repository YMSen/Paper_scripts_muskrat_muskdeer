#! usr/bin/perl -w
use strict;

die "<head><file>\n" unless @ARGV==2;
open(IN,$ARGV[0])||die;
my %scaffold;
while(<IN>)
{
	chomp;
	my @cut=split /\s+/;
	$scaffold{$cut[0]}=1; 
#	print $cut[0]."\n";
}
close IN;
open(IN,$ARGV[1])||die;
while(<IN>)
{
	chomp;
	my @cut=split /\s+/;
	if(!exists $scaffold{$cut[0]})
	{
		print  $_."\n";
	}
}

close IN;



