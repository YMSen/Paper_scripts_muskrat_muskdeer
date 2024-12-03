#!/usr/bin/perl -w
use strict;

die "use <gff3>" unless @ARGV==1;

my $gff3=shift;
my %mRNA;

open(IN,$gff3)||die;
while(<IN>){
       chomp;
       next if(/^$/);
       next if(/^#/);
       next if(/^\s+/);
       my @cut=split /\t/;
       my $id;
       if($cut[2] eq 'CDS'){
                  ($id)=$cut[8]=~/(Parent=.*)$/;
                  $id=~s/Parent/ID/;
                  if(!exists $mRNA{$id}{s} || $mRNA{$id}{s}>$cut[3]){$mRNA{$id}{s}=$cut[3];}
                  if(!exists $mRNA{$id}{l} || $mRNA{$id}{l}<$cut[4]){$mRNA{$id}{l}=$cut[4];}
       }
}
close IN;

open(IN,$gff3)||die;
while(<IN>){
       chomp;
       next if(/^$/);
       next if(/^#/);
       next if(/^\s+/);
       my @cut=split /\t/;
       my $id;
       my $flag=0;
       if($cut[2] eq 'mRNA'){
                  $flag=1;
                  ($id)=$cut[8]=~/(ID=[^;]+)/;
                  $cut[3]=$mRNA{$id}{s};
                  $cut[4]=$mRNA{$id}{l}; 
                  
       }elsif($cut[2] eq 'CDS'){
                  $flag=1;
                  ($id)=$cut[8]=~/(Parent=.*)$/;
       }
       if($flag==1){
   	    my $line=join("\t",@cut[0..7]);
	    print "$line\t$id;\n";
       }
}

close IN;
