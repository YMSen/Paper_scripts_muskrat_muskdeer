#!/usr/bin/perl -w
use strict;

die "use <list><gff>" unless @ARGV==2;

open(IN,"$ARGV[0]")||die;
my %idList;
while(<IN>){
       chomp;
       my $id=(split /\t/)[0];
       $idList{$id}=1;
}
close IN;

open(IN,"$ARGV[1]")||die;
while(<IN>){
       next if(/^\s+/);
       chomp;
       my $id;
       if(/mRNA.*ID=([^;]+)/){
               $id=$1;
       }elsif(/CDS.*Parent=([^;]+)/){
               $id=$1;
       }else{
		next;}
       if(exists $idList{$id}){
                 print $_."\n";
       }
}
close IN;
