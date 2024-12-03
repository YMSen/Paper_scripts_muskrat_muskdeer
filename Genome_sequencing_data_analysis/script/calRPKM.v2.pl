#!/usr/bin/perl
use strict;
use warnings;

##----------------------------------------------------------------------------
## Gene length and reads number information is needed
##
## rpkm =    ( No. of reads in the gene) *1000,000,000
##       -------------------------------------------------
##        (No. of all reads in genes)*(the gene length)
##----------------------------------------------------------------------------

die "Usage:perl $0 <readcount> <gtf> <out>" unless @ARGV==3;

open RC, "$ARGV[0]" or die $!;
open GTF, "$ARGV[1]" or die $!;
open OUT, ">$ARGV[2]" or die $!;

my %rc;
my %trans;
my $total = 0;
while(<RC>){
    chomp;
    last if(/no_feature/);
    my @tmp = split /\s+/, $_;
    ${$rc{$tmp[0]}} = [$tmp[1], 1e20, 0];
    $total += $tmp[1];
}
die "Error:readcount is zero!" if($total == 0);
while(<GTF>){
    chomp;
    my @tmp = split /\t/, $_;
    next if ($tmp[2] ne "exon");
    my $name = $1 if($tmp[8] =~ /gene_id "(.*?)";/);
    my $transName = $1 if($tmp[8] =~ /transcript_id "(.*?)";/);
    ${$trans{$name}{$transName}}=0;
    if(exists $rc{$name}){
#        ${$rc{$name}}->[1] = $tmp[3] if ($tmp[3]<${$rc{$name}}->[1]);
#        ${$rc{$name}}->[2] = $tmp[4] if ($tmp[4]>${$rc{$name}}->[2]);
         ${$rc{$name}}->[2] += $tmp[4] - $tmp[3] +1; #Only exon is used
    }
}

my %transNum;
foreach my $transkey(keys %trans){
    $transNum{$transkey} = 0;
    foreach my $transkeys2(keys %{$trans{$transkey}}){
        $transNum{$transkey} += 1;
    }
}


foreach my $key(keys %rc){
#   my $length = ${$rc{$key}}->[2] - ${$rc{$key}}->[1];
    $transNum{$key} = 1 unless(exists $transNum{$key}); 
    my $length = int(${$rc{$key}}->[2] / $transNum{$key} + 0.5);
    next if($length <= 0);
    my $rpkm = (${$rc{$key}}->[0] * 1e9)/($total * $length);
    print OUT "$key\t${$rc{$key}}->[0]\t$rpkm\n";
    
}
close RC;
close GTF;
close OUT;
