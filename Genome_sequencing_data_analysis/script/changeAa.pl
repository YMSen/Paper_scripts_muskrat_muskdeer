#!/usr/bin/perl -w
use strict;
my $in=shift;
open IN,$in;
$/=">";$/=<IN>;$/="\n";
while(<IN>){
	chomp;
	(my $id=$_)=~s/\s+.*$//;
	$/=">";
        my $seq=<IN>;
        chomp $seq;
        $seq=~ tr/[a-z]/[A-Z]/;
        $/="\n";
	print ">".$_."\n".$seq;
}
