#!/usr/bin/perl

my $file = shift;
my $output;
open OUT, ">$file.table" || die "fail";
open IN, $file || die "fail";
while (<IN>) {
        chomp;
	$output = "";
        my @t = split /\t/;
        my $len = $t[3]-$t[2]+1;
        $output .= "$t[0]\t$t[5]\t$t[4]\t$t[7]\t$t[8]\t$len\n";
	print OUT $output;
}
close IN;
#open OUT, ">$file.table" || die "fail";
#print OUT $output;
close OUT;

