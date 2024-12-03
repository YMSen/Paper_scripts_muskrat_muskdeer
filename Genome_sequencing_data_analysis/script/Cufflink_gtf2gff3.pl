#!/usr/bin/perl -w

use strict;
my $gtf_change=shift;
open IN,"$gtf_change";
my $id;
my $target;
while(<IN>){
	chomp;
	if ($_ =~ /^\S+/){
        my @arry=split /\t/,$_;
	    if($arry[8]=~ /transcript_id "([^"]+)"(.+)/){
		$id =$1;
		$target =$2;
	    }
	    print "$arry[0]\tCUFF\tcDNA_match\t$arry[3]\t$arry[4]\t$arry[5]\t$arry[6]\t$arry[7]\tID=$id; Target=$id$target\n" if($arry[2]=~/exon/);
     }else{
         next;
     }
}
