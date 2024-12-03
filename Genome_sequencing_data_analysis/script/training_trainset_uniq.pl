#!/usr/bin/perl -w
use strict;

die "Usage:<set_gff>\n" unless @ARGV==1;
my $set_gff=shift;

my (%uniq,%ids,%idss,%hash,%final);

open IN,$set_gff;
  while(<IN>){
    chomp;
    my @a = split /\t/,$_;
    my $only=$1 if ($a[8] =~ /=(.+?)[\-;]/);
    my $id = $1 if ($a[8] =~ /=(\S+?);/);
    if($a[2] eq 'CDS'){
      $uniq{$only}{$id}++;
    }
  }
close IN;

foreach my $k(sort keys %uniq){
  foreach my $kk(keys %{$uniq{$k}}){
    my $v=$uniq{$k}{$kk};
    if(!$ids{$k}){
       $ids{$k}=$v;
       $idss{$k}=$kk;
    }elsif($ids{$k}){
       if($v>$ids{$k}){
          $ids{$k}=$v;
          $idss{$k}=$kk;
       }else{
          next;
       }
    }    
  }
}

foreach my $ks(keys %idss){
  my $vs=$idss{$ks};
  $final{$vs}++;
}

open IN1,$set_gff;
  while(<IN1>){
    chomp;
    my @a = split /\t/,$_;
    my $id = $1 if ($a[8] =~ /=(\S+?);/);
    if($final{$id}){
      $hash{$id}{out} .= "$_\n";
    }
  }
close IN1;

open OUT,">$set_gff.uniq.gff";
foreach my $k (sort keys %hash){
  print OUT $hash{$k}{out};
}
close OUT;

